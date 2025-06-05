%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

close all; clc; clearvars;
root = '.';
addpath(genpath(fullfile(root,'pulseq')));
addpath(genpath(fullfile(root,'pulseq_skope'))); 
addpath(genpath(fullfile(root,'utils'))); 


% trigChannel: 'ext1' used for 7T
% for dynamic field monitoring [H], set RF to ~0 
fa  = 1e-100;     % flip angle for excitation
faf = 1e-100;     % flip angle for fat saturation

% % for imaging 
% fa  = 90; 
% faf = 110;

maxgrad = 51; % empirical value, for PNS consideration
gradSafetyMargin = maxgrad / 70;      

%% Scanner specs
seq_name     = "Spiral2d_7T";          
seq_type     = "2D Spiral Stitch";
scanner_type = "Siemens Terra 7T SC72CD"; %"Siemens dotplus 10.5T HG"; %"Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);
sk.sys.maxGrad = 60*1e-3 * sk.sys.gamma;  % [Hz], restriction for other gradients, except for readout 
sk.sys.maxSlew = 170     * sk.sys.gamma;  % [Hz]
% sequence parameters
seq_params = struct('trigChannel', 'ext1'); % 'osc0','osc1','ext1'
seq_params.fov                 = 200*1e-3; % [m], 230e-3; 
seq_params.N                   = 400;      % 375;%300;%188;%153; %230; %384; %288;%230; %256; %220; % Define FOV and resolution
seq_params.accelerationFactor  = 4;      % 30;%16;%5;%4; % acceleration factor
seq_params.alpha               = fa;     % flip angle
seq_params.alpha_fatsat        = faf;    % flip angle for fat saturation
seq_params.thickness           = 2e-3;   % [m], slice
seq_params.Nslices             = 1;      
seq_params.sliceGap            = 1;      % 10; % slice gap in fraction of slice thickness. 
seq_params.TE                  = 5e-3;   % [s]
seq_params.TR                  = 500e-3; % [s], 200e-3; 200e-3; 
seq_params.nRepeats            = 1;      % 5; 3; 100;
seq_params.maxAdcSegmentLength = 1000;   % 2200; %1000; % number of ADC samples per segment preferrably <=1000 and has to be <=8192.


% spiral design
seq_params.isTimeOptimal            = true; % true for time optimal variable density spiral design. false for archimedean spiral design
% seq_params.isRotationallyVariant    = true; % only in effect when seq_params.isTimeOptimal= true;
seq_params.nSpiralInterleaves       = 1;
seq_params.gradSafetyMargin         = gradSafetyMargin;   % grad safety margin for spiral

% skope relevant
seq_params.probeType                = 'H';    % 'H' for proton, 'F' for fluorine
seq_params.probeRadius              = 0.4e-3; % radius of the field probe in m. 
seq_params.signalCutoff             = 0.70;   % signal cutoff level.
seq_params.nSegments2measure        = 0;      %1;% number of gradient segments to measure and stitch. When set to 0, number of segments will be determined automatically.
seq_params.stitchMode               = 'concurrent'; % 'concurrent: No shot is repeated. In each shot, a given gradient segment is measured. This scheme is compatible with concurrent field monitoring. 
% sequential: Each shot is repeated to measure all gradient segments in a row
% before moving to the next shot. Perhaps should not be used in practice. 
% interleaved: Entire sequence is repeated to measure only one gradient segment at a time throughout the sequence. 
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
%%
res = round(1e3*sk.seq_params.resolution, 2);
a = fix(res);
b = (res - a)*100;
if mod(b, 10) == 0
    b = b/10;
end

fn= ['7T_',num2str(a),'p',num2str(b),'_',num2str(seq_params.N),'_r',num2str(sk.seq_params.accelerationFactor)];
fn= [fn,'_seg',num2str(sk.seq_params.nSegments2measure)];
fn = [fn, '_max', num2str(maxgrad)];
if fa > 1 
    fn= [fn, '_fa', num2str(seq_params.alpha)];
end 
%% PNS calc
warning('OFF', 'mr:restoreShape');
[pns_ok, pns_n, pns_c, tpns] = sk.seq.calcPNS('MP_GPA_K2259_2000V_650A_SC72CD_EGA.asc'); % TERRA-XJ

if (pns_ok)
    fprintf('PNS check passed successfully\n');
else
    fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
end
nTRs=1; sk.Plot(nTRs);
%% 
sk.Write(fn);

sk.Check_timing();

% sk.TestReport();



