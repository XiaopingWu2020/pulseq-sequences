%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

close all; clc; clearvars;
addpath(genpath(fullfile('..','pulseq')));
addpath(genpath(fullfile('..','utils'))); 
addpath(genpath(fullfile('.','src'))); % for helper functions
addpath(genpath(fullfile('.','plot'))); 

%% Scanner specs
seq_name     = "Epi2d_7T"; %"Epi2d_7T";          
seq_type     = "2D EPI Stitch";
scanner_type = "Siemens Terra 7T SC72CD"; %"Siemens dotplus 10.5T SC72CD"; %"Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);
% sequence parameters
seq_params = struct('trigChannel', 'osc0'); % 'osc0','osc1','ext1'
seq_params.fov                      = 200e-3; % [m]
seq_params.N                        = 400; %153; %92; %115;%230; %64; % Define FOV and resolution

seq_params.accelerationFactor       = 6;      % acceleration factor
seq_params.partialFourier           = 1;      % partial Fourier factor: 1: full sampling 0: start with ky=0

seq_params.alpha                    = 25;     % [degree] flip angle
seq_params.alpha_fatsat             = 110;    % [degree] flip angle for fat saturation
seq_params.thickness                = 2e-3;   % slice thickness
seq_params.Nslices                  = 1;      
seq_params.sliceGap                 = 1;      %10; % slice gap in fraction of slice thickness. 
seq_params.TE                       = 22e-3;  % s
seq_params.TR                       = 500e-3; %150e-3; % s
seq_params.nRepeats                 = 1;
seq_params.nNavigators              = 0;
seq_params.readoutTime              = 13e-4;   %10e-4; % 1/bandwidthPerPixel

% skope relevant
seq_params.probeType                = 'H';    % 'H' for proton, 'F' for fluorine
seq_params.probeRadius              = 0.4e-3; % radius of the field probe in m. 
seq_params.signalCutoff             = 0.35;   % signal cutoff level.
seq_params.nSegments2measure        = 0;      %1;% number of gradient segments to measure and stitch. When set to 0, number of segments will be determined automatically.
seq_params.isHighRes                = true;   % works when nSegments2measure is set to 0. true for "variable-segment" , false for "constant-segment"
seq_params.stitchMode               = 'concurrent';  % only concurrent has been implemented by now.
% 1. concurrent: No shot is repeated. In each shot, a given gradient segment is measured. 
%   This scheme is compatible with concurrent field monitoring. 
% 2. sequential: Each shot is repeated to measure all gradient segments in a row
%   before moving to the next shot. Perhaps should not be used in practice. 
% 3. interleaved: Entire sequence is repeated to measure only one gradient segment at a time throughout the sequence. 
% Both sequential and interleaved schemes are applicable to any sequence, but at the cost of requiring a calibration longer than
% the native scan.

seq_params.interSessionDelay        = 6;      % [s], only used for interleaved mode. 
seq_params.nInterleaves             = 1;      % number of excitations per dynamic.
seq_params.skopeMinTR               = 110e-3; % s
seq_params.gradFreeDelay            = 200e-6; % s
seq_params.nPrescans                = 5;      % for sync between skope and scanner
seq_params.useSingleAdcSegment4Sync = false;  %true; % 
% seq_params.doFastestPrescan         = false;  % true: use minimum TR; false: use a TR greater than skope interleaveTR.

%%% sequence will calculate the following parameters:
%%% seq_params.skopeInterleaveTR, seq_params.skopeAcqDuration, seq_params.trigger2AdcTime 

sk.Prepare(seq_params);
sk.SkopeReport();

fn= ['xw_ep2d-',num2str(1e3*sk.seq_params.resolution,2),'mm-r',num2str(sk.seq_params.accelerationFactor)];

grad0 = sk.seq_params.grad0;
save([fn,'.mat'],'grad0')

if sk.seq_params.nNavigators>0
    sk.Write('xw_ep2d-navi');
else
    if sk.seq_params.nPrescans>0
        if sk.seq_params.useSingleAdcSegment4Sync
            sk.Write([fn,'-1adc4sync'])
        else
            if sk.seq_params.partialFourier<1
                sk.Write(fn);
            else
                sk.Write([fn,'-noPF']);
            end
        end
    else
        sk.Write('xw_ep2d-lean')
    end

end

% Debug
sk.Check_timing();

nTRs=5;
sk.Plot(nTRs);

if 0
sk.TestReport();
end
%% PNS calc
% warning('OFF', 'mr:restoreShape');
% [pns_ok, pns_n, pns_c, tpns] = sk.seq.calcPNS('MP_GPA_K2259_2000V_650A_SC72CD_EGA.asc'); % TERRA-XJ
% 
% if (pns_ok)
%     fprintf('PNS check passed successfully\n');
% else
%     fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
% end
%% k-space trajectory calculation
% [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = sk.seq.calculateKspacePP();
% % plot k-spaces
% % figure; plot(t_ktraj, ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
% fig = plot_kspace(ktraj, ktraj_adc);

%% plot the stitching pattern
[fig] = plot_stitch(sk.seq_params.grad0, sk.seq_params.triggerDelays, ...
    sk.seq_params.probeType, sk.seq_params.probeRadius, sk.seq_params.signalCutoff, ...
    sk.seq_params.nSegments2measure);
% print(fig, '-dpng', '-loose', '-r300', '-image', [fn,'.png');