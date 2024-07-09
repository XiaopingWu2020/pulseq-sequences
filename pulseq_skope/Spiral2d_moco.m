%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

close all; clc; clearvars;
addpath(genpath(fullfile('..','pulseq')));
addpath(genpath(fullfile('.','src'))); % for helper functions
% addpath('~/matlab/skope_TCP_client/methods')
% addpath('~/matlab/skope_MatlabAqSysDataImport/bin')

%% Scanner specs
seq_name = "Spiral2d_moco";          
seq_type = "2D Spiral MoCo";
scanner_type = "Siemens Terra 7T SC72CD"; %"Siemens dotplus 10.5T SC72CD"; %"Siemens dotplus 10.5T HG LP"; %"Siemens dotplus 10.5T HG"; %"Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);
% sequence parameters
seq_params = struct('trigChannel', 'osc0'); % 'osc0','osc1','ext1'

seq_params.fov= 220e-3;%152e-3;%210e-3; %256e-3; %230e-3; 220 is about the 95th percentile of the distance between tip of nose and back of head.  150 is a good estimate of the brain height.
seq_params.N= 220;%275;%220; %210; %110; %840; %525; %380; %608; %263;%275; %315;%575;%288;%115;%128;%153; %230; %384; %288;%230; %256; %220; % Define FOV and resolution
seq_params.accelerationFactor= 5; %6; %5; % acceleration factor

% spiral design
seq_params.isTimeOptimal= true; % true for time optimal variable density spiral design. false for archimedean spiral design
%seq_params.isRotationallyVariant= true; % only in effect when seq_params.isTimeOptimal= true;
seq_params.nSpiralInterleaves= 1;
seq_params.gradSafetyMargin= 0.75;%0.85; % 0.9 1;%0.96;% 

%seq_params.alpha=10;%25;                  % flip angle
seq_params.thickness= 2e-3; %120e-3;            % slice
seq_params.Nslices=1; %3;      
seq_params.sliceGap= 5;%10; %1; %10; % slice gap in fraction of slice thickness. 
seq_params.TE = 1e-3.*16;%[16 20]; %22e-3; % s
seq_params.TR = 1;%1;%400e-3; % s % 0.5 sec tr was found too short for realtime motion tracking, leading to timeout while sending data error. however, using 1 sec tr for ~10 min scan also resulted in timeout while sending data error.
seq_params.nRepeats = 5;%1500;%600;%1500; %3000;%200;
seq_params.maxAdcSegmentLength= 1000; %2200; %1000; % number of ADC samples per segment preferrably <=1000 and has to be <=8192.

% motion detection
seq_params.trigChannelMoco= 'ext1';
seq_params.gradToneLength= 1e-3; %0.6e-3; %1e-3; % length of tones in second. 
seq_params.toneFrequencies= 1e3*[6 7 8]; %[2, 4, 6]; %1e3*[6, 7, 8]; % nominal tone frequencies [fx, fy, fz] in hz
seq_params.toneAmplitude= 0.01; %4; %3.71; % constant max tone amplitude for gx, gy, and gz in mT/m.
seq_params.image2gradToneDelay= 4e-3; %4e-3; % delay between image readout gradients and the position encoding gradient tones. 
% for single echo, insert gradient tones between excitation and image readout. 
% for multiple echoes, insert gradient tones after image readout.
seq_params.maximizeTimeEfficiency= false; %true; % this flag only has effect for multiple echo acquisition. 
% false: insert gradient tones between image readout and spoilers. 
% true; insert gradient tones superimposed gradient spoilers after image readout.
%
% === field tracking scheme ===
%'T': gradient tones; 'TF': gradient tones and background field; 'TFG': gradient tones, background field and dynamic field.
%'F': background field; 'FG': background and dynamic field.
%'G': dynamic field.
seq_params.fieldTrackingScheme= 'G'; %'TFG'; %'TFG'; 
% ========================

% skope relevant
seq_params.nInterleaves= 1; % number of excitations per dynamic.
seq_params.skopeMinTR= 110e-3; % s
seq_params.gradFreeDelay= 200e-6; % s
seq_params.nPrescans= 5; %10;%20; % for sync between skope and scanner
seq_params.useSingleAdcSegment4Sync= false; %false; % 
%seq_params.doFastestPrescan= false; % found that the sync TR was also restricted by skope interleave TR at least on our systems true: use minimum TR; false: use a TR greater than skope interleaveTR.

%%% sequence will calculate the following parameters:
%%% seq_params.skopeInterleaveTR, seq_params.skopeAcqDuration, seq_params.trigger2AdcTime 

%%% tcp control
% 105t
%host = '134.84.19.164'; %'skope-105T'. if empty, no measurement automation for field monitoring.
%host_moco = '134.84.19.37'; %'skope-terra'. if empty, no measurement automation for motion detection.

% terra
% host = '134.84.19.51'; %'skope-terra'. if empty, no measurement automation for field monitoring.
%
sk.Prepare(seq_params);
sk.SkopeReport();

fn= ['xw_sp2d_moco-',num2str(1e3*sk.seq_params.resolution,2),'mm-r',num2str(sk.seq_params.accelerationFactor),...
    '-te',num2str(round(1e3*sk.seq_params.TE)),'-tr',num2str(sk.seq_params.TR),'-sli',num2str(sk.seq_params.Nslices)];
if sk.seq_params.nPrescans>0
    if sk.seq_params.useSingleAdcSegment4Sync
        fn= [fn,'-1adc4sync'];    
    end
else
    fn= [fn,'-noSync'];
end
sk.Write(fn);

sk.Write('external');

% Debug
sk.Check_timing();

nTRs=2;
sk.Plot(nTRs);

if 0
    sk.TestReport();
end

% %% tcp control
% % field monitoring
% if exist('host','var')&& ~isempty(host)
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

