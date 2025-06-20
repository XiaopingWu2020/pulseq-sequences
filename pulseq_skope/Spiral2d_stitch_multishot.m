%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

close all; clc; clearvars;
addpath(genpath(fullfile('..','pulseq')));
addpath(genpath(fullfile('..','utils'))); 
addpath(genpath(fullfile('.','src'))); % for helper functions
addpath(genpath(fullfile('.','plot'))); 

%% Scanner specs
seq_name     = "Spiral2d_7T";          
seq_type     = "2D Spiral Stitch MultiShot";
scanner_type = "Siemens dotplus 10.5T SC72CD"; %"Siemens dotplus 10.5T HG"; %"Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);

% Gmax and Smax for conventional gradients.
% NOTE: the Gmax and Smax should be set in Set_scanner.m
% Let's use the gradSafetyMargin param to be conservative. 
% sk.sys.maxGrad = 60*1e-3 * sk.sys.gamma;  % [Hz]
% sk.sys.maxSlew = 170     * sk.sys.gamma;  % [Hz]

% sequence parameters
seq_params = struct('trigChannel', 'osc0'); % 'osc0','osc1','ext1'
seq_params.fov                 = 150e-3; % [m], 230e-3; 
seq_params.N                   = 500;    % Define FOV and resolution
seq_params.accelerationFactor  = 1;     % acceleration factor
seq_params.alpha               = 10;     % [degree] flip angle
seq_params.alpha_fatsat        = 110;    % [degree] flip angle for fat saturation
seq_params.thickness           = 2e-3;   % [m], slice
seq_params.Nslices             = 1;      
seq_params.sliceGap            = 1;      % slice gap in fraction of slice thickness. 
seq_params.TE                  = 5e-3;   % [s]
seq_params.TR                  = 500e-3; % [s]
seq_params.nRepeats            = 1;      % 
seq_params.maxAdcSegmentLength = 1000;   % number of ADC samples per segment preferrably <=1000 and has to be <=8192.

% spiral design
seq_params.isTimeOptimal            = true; % true for time optimal variable density spiral design. false for archimedean spiral design
% seq_params.isRotationallyVariant    = true; % only in effect when seq_params.isTimeOptimal= true;

% Let's use nSpiralInterleaves as number of shots so that the
% accelerationFactor specified will give the true acceleration factor. 
seq_params.nSpiralInterleaves       = 30;      % number of Interleaves in trajectory design

% For spiral readout design, Gmax and Smax were scaled by gradSafetyMargin due to hardware limitations.
seq_params.gradSafetyMargin         = 0.75;   %0.94

% skope relevant
seq_params.probeType                = 'F';    % 'H' for proton, 'F' for fluorine
seq_params.probeRadius              = 0.4e-3; % radius of the field probe in m. 
seq_params.signalCutoff             = 0.41;   % signal cutoff level.
seq_params.nSegments2measure        = 0;      %1;% number of gradient segments to measure and stitch. When set to 0, number of segments will be determined automatically.
seq_params.stitchMode               = 'concurrent'; 
% 1. concurrent: No shot is repeated. In each shot, a given gradient segment is measured. 
%   This scheme is compatible with concurrent field monitoring. 
% 2. sequential: Each shot is repeated to measure all gradient segments in a row
%   before moving to the next shot. Perhaps should not be used in practice. 
% 3. interleaved: Entire sequence is repeated to measure only one gradient segment at a time throughout the sequence. 
% Both sequential and interleaved schemes are applicable to any sequence, but at the cost of requiring a calibration longer than
% the native scan.

seq_params.interSessionDelay        = 6;      % [s], only used for interleaved mode. 
seq_params.nInterleaves             = 1;      % number of excitations per dynamic.
seq_params.skopeMinTR               = 110e-3; % [s]
seq_params.gradFreeDelay            = 200e-6; % [s]
seq_params.nPrescans                = 0;      % for sync between skope and scanner
seq_params.useSingleAdcSegment4Sync = false;  % 
% seq_params.doFastestPrescan         = false;  % found that the sync TR was also restricted by skope interleave TR at least on our systems true: use minimum TR; false: use a TR greater than skope interleaveTR.

%%% sequence will calculate the following parameters:
%%% seq_params.skopeInterleaveTR, seq_params.skopeAcqDuration, seq_params.trigger2AdcTime 

sk.Prepare(seq_params);
sk.SkopeReport();

fn= ['xw_sp2d-',num2str(1e3*sk.seq_params.resolution,2),'mm-r',num2str(sk.seq_params.accelerationFactor)];

grad0 = sk.seq_params.grad0;
save([fn,'.mat'],'grad0')

if sk.seq_params.nSegments2measure> 1
    fn= [fn,'-',num2str(sk.seq_params.nSegments2measure),'segs'];
end

if sk.seq_params.nPrescans>0
    if sk.seq_params.useSingleAdcSegment4Sync
        fn= [fn,'-1adc4sync'];    
    end
else
    fn= [fn,'-noSync'];
end

switch sk.seq_params.stitchMode
    case 'interleaved'
        fn= [fn,'-itlv'];
    case 'sequential'
        fn= [fn,'-seq'];
    otherwise
end

sk.Write(fn);

sk.Write('external');

% Debug
sk.Check_timing();

nTRs=100;
sk.Plot(nTRs);

if 0
    sk.TestReport();
end
%% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = sk.seq.calculateKspacePP( ...
    'blockRange', [1, sk.seq_params.end_block(1)]);
% plot k-spaces
figure; plot(t_ktraj, ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
fig = plot_kspace(ktraj, ktraj_adc);

%% PNS calc
warning('OFF', 'mr:restoreShape');
[pns_ok, pns_n, pns_c, tpns] = sk.seq.calcPNS('MP_GPA_K2259_2000V_650A_SC72CD_EGA.asc'); % TERRA-XJ

if (pns_ok)
    fprintf('PNS check passed successfully\n');
else
    fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
end

%% plot the stitching pattern
[fig] = plot_stitch_multishot(sk.seq_params.grad0, sk.seq_params.triggerDelays, ...
    sk.seq_params.probeType, sk.seq_params.probeRadius, sk.seq_params.signalCutoff, ...
    sk.seq_params.nSegments2measure);
% print(fig, '-dpng', '-loose', '-r300', '-image', '2DSpiral_MultiShot_Stitch_Variable.png');