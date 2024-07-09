%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

close all; clc; clearvars;
addpath(genpath(fullfile('..','pulseq')));
% addpath('~/matlab/skope_TCP_client/methods')
% addpath('~/matlab/skope_MatlabAqSysDataImport/bin')

%% Scanner specs
seq_name = "Epi2d_moco"; %"Epi2d_7T";          
seq_type = "2D EPI MoCo"; % used to specify the correct prep function
scanner_type = "Siemens Terra 7T SC72CD"; %"Siemens dotplus 10.5T SC72CD"; %"Siemens dotplus 10.5T HG"; %"Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);
% sequence parameters
seq_params = struct('trigChannel', 'osc0'); % 'osc0','osc1','ext1'

seq_params.fov=220e-3; %152e-3; %210e-3; % 180 is a good estimate of the head size in LR direction. 
seq_params.N= 110;%220; %608; %380;%275;%115;%230; %64; % Define FOV and resolution

seq_params.accelerationFactor= 3; % acceleration factor
seq_params.partialFourier= 0.5; %1; %0.75;%0.25; %0.5 (pf6/8);       % partial Fourier factor: 1: full sampling 0: start with ky=0

%seq_params.alpha=10;%25;                  % flip angle
seq_params.thickness= 2e-3; %120e-3;            % slice
seq_params.Nslices=1;      
seq_params.sliceGap= 1;%10; % slice gap in fraction of slice thickness. 
seq_params.TE = 1e-3*16;%[16 20]; %22e-3; % s
seq_params.TR = 1; %400e-3; %150e-3; % s
seq_params.nRepeats = 10; %180;%2;%1500;%3000; %200;
seq_params.nNavigators = 0;
seq_params.readoutTime = 0.70e-3; %0.52e-3;%0.38e-3; % min with HG for 0.8 mm iso % 1/bandwidthPerPixel

% motion detection
seq_params.trigChannelMoco= 'ext1'; %'osc0'; %'ext1';
seq_params.gradToneDuration= 1e-3; %0.6e-3; %1e-3; % duration of tones in second. 
seq_params.toneFrequencies= 1e3*[6, 7, 8]; %[2, 4, 6]; %1e3*[6, 7, 8]; % nominal tone frequencies [fx, fy, fz] in hz
seq_params.toneAmplitude= 0.01; %4; %3.71; % constant max tone amplitude for gx, gy, and gz in mT/m.

seq_params.fieldControlDuration= 1e-3; % duration (after position encoding gradient tones) during which background field is monitored. 
% by default, insert gradient tones after image readout.
% seq_params.maximizeTimeEfficiency= false; %true; % this flag only has effect for multiple echo acquisition. 
% % false: insert gradient tones between image readout and spoilers. true;
% % insert gradient tones superimposed gradient spoilers after image readout.
%
% === field tracking scheme ===
%'T': gradient tones; 'TF': gradient tones and background field; 'GTF': gradient tones, background field and dynamic field.
%'F': background field; 'GT': dynamic field and tones.
%'G': dynamic field.
seq_params.fieldTrackingScheme= "G"; %"GTF"; %'GTF'; 
% ========================

% skope relevant
seq_params.nInterleaves= 1; % number of excitations per dynamic.
seq_params.skopeMinTR= 110e-3; % s
seq_params.gradFreeDelay= 200e-6; % s
seq_params.nPrescans= 10;%20; % for sync between skope and scanner
seq_params.useSingleAdcSegment4Sync= false; %true; % 
%seq_params.doFastestPrescan= false; % true: use minimum TR; false: use a TR greater than skope interleaveTR.

%%% sequence will calculate the following parameters:
%%% seq_params.skopeInterleaveTR, seq_params.skopeAcqDuration, seq_params.trigger2AdcTime 

%%% tcp control
% 10.5t
%host = '134.84.19.164'; %'skope-105T'. if empty, no measurement automation for field monitoring.
%host_moco = '134.84.19.37'; %'skope-terra'. if empty, no measurement automation for motion detection.

% terra
% host = '134.84.19.51'; %'skope-terra'. if empty, no measurement automation for field monitoring.
%
sk.Prepare(seq_params);
sk.SkopeReport();

fn= ['xw_ep2d_moco-',num2str(1e3*sk.seq_params.resolution,2),'mm-r',num2str(sk.seq_params.accelerationFactor)];    
if sk.seq_params.nNavigators>0
    fn= [fn,'-navi'];
else

    if sk.seq_params.nPrescans>0
        if sk.seq_params.useSingleAdcSegment4Sync
            fn= [fn,'-1adc4sync'];
        else
            if sk.seq_params.partialFourier==1
                fn= [fn,'-noPF'];
            end
        end
    else
        fn= [fn,'-lean'];
    end

end
sk.Write(fn);

sk.Write('external');

% Debug
sk.Check_timing();

nTRs=1;%50;
sk.Plot(nTRs);

if 0
    sk.TestReport();
end

% %% tcp control
% % field monitoring
% if exist('host','var')&&~isempty(host)
%     disp('                          ')
%     input(['If you like to do realtime field tracking for motion calibration, please go set it up on a different client. \n' ...
%         'Press Enter when you are ready to proceed...'],'s');
%     
%     input(['If you like to do field control, please go set it up on a different client. \n' ...
%         'Press Enter when you are ready to proceed...'],'s');
% 
%     tcpAutomateScan([fn,'.seq'], host, false);
% end
% 
% % motion tracking
% if exist('host_moco','var')&& ~isempty(host_moco)
%     disp('                          ')
%     input(['If you like to do motion tracking, please go set it up on a different client. \n' ...
%         'Press Enter when you are ready to proceed...'],'s');
% 
%     tcpAutomateScan([fn,'.seq'], host_moco, true);
% end
