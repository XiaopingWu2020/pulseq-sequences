%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

close all; clc; clear all;
addpath(genpath('../pulseq/')); % Pulseq path

%% Scanner specs
seq_name = "Epi2d_7T"; %"Epi2d_7T";          
seq_type = "2D EPI";
scanner_type = "Siemens dotplus 10.5T HG"; %"Siemens dotplus 10.5T SC72CD"; %"Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);
% sequence parameters
seq_params = struct('trigChannel', 'osc0'); % 'osc0','osc1','ext1'
seq_params.fov=230e-3; 
seq_params.N= 288; %153; %92; %115;%230; %64; % Define FOV and resolution

seq_params.accelerationFactor= 6; % acceleration factor
seq_params.partialFourier= 0.5;       % partial Fourier factor: 1: full sampling 0: start with ky=0

seq_params.alpha=25;                  % flip angle
seq_params.thickness=2e-3;            % slice
seq_params.Nslices=1;      
seq_params.sliceGap= 1;%10; % slice gap in fraction of slice thickness. 
seq_params.TE = 22e-3; % s
seq_params.TR = 400e-3; %150e-3; % s
seq_params.nRepeats = 2;
seq_params.nNavigators = 0;
seq_params.readoutTime = 5e-4; %10e-4; % 1/bandwidthPerPixel
% skope relevant
seq_params.nInterleaves= 1; % number of excitations per dynamic.
seq_params.skopeMinTR= 110e-3; % s
seq_params.gradFreeDelay= 200e-6; % s
seq_params.nPrescans= 20; % for sync between skope and scanner
seq_params.useSingleAdcSegment4Sync= false; %true; % 
%seq_params.doFastestPrescan= false; % true: use minimum TR; false: use a TR greater than skope interleaveTR.

%%% sequence will calculate the following parameters:
%%% seq_params.skopeInterleaveTR, seq_params.skopeAcqDuration, seq_params.trigger2AdcTime 

sk.Prepare(seq_params);
sk.SkopeReport();

if sk.seq_params.nNavigators>0
    sk.Write('xw_ep2d-navi');
else

    fn= ['xw_ep2d-',num2str(1e3*sk.seq_params.resolution,2),'mm-r',num2str(sk.seq_params.accelerationFactor)];
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
sk.Write('external');

% Debug
sk.Check_timing();

nTRs=50;
sk.Plot(nTRs);

if 0
sk.TestReport();
end