%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

close all; clc; clearvars;
addpath(genpath('../pulseq/')); % Pulseq path

%% Scanner specs
seq_name = "Gre2d_7T";          
seq_type = "2D GRE";
scanner_type = "Siemens dotplus 10.5T SC72CD"; %"Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);
% sequence parameters
seq_params = struct('trigChannel', 'osc0'); % 'osc0','osc1','ext1'
seq_params.fov= 220e-3; %210e-3; %230e-3; 
seq_params.N= 220; %116; % Define FOV and resolution
seq_params.partialFourier= 1; %6/8;
seq_params.reordering=  'ascending'; %'ascending'; % 'descending'
%seq_params.alpha= 20; %25;                  % flip angle
seq_params.thickness= 2e-3; %2e-3;            % slice
seq_params.sliceGap= 5;%10; % slice gap in fraction of slice thickness. 
seq_params.Nslices=1;    
seq_params.nRepeats = 1;
seq_params.TE =  1e-3*(1:2); %[2.71 3.39]*1e-3; %[5 5+3*1.3] * 1e-3; %1e-3*[1 2]; % s
seq_params.TR = 200e-3; %100e-3; % s
seq_params.readoutTime= 0.5e-3; %3.2e-3;              % ADC duration
seq_params.gradSafetyMargin= 0.8; %1;%0.96;% 

% skope relevant
seq_params.nInterleaves= 1;
seq_params.skopeMinTR= 110e-3; % s
seq_params.gradFreeDelay= 200e-6; % s
seq_params.nPrescans= 10; % for sync between skope and scanner
%%% sequence will calculate the following parameters:
%%% seq_params.skopeInterleaveTR, seq_params.skopeAcqDuration, seq_params.trigger2AdcTime 

sk.Prepare(seq_params);
sk.SkopeReport();

%fn= ['xw_gre2d-',num2str(1e3*sk.seq_params.resolution,2),'mm-r',num2str(sk.seq_params.accelerationFactor)];
fn= ['xw_gre2d-',num2str(1e3*sk.seq_params.resolution,2),'mm-',sk.seq_params.reordering,'-sli',num2str(sk.seq_params.Nslices)];
       
if sk.seq_params.nPrescans>0
    sk.Write(fn);
else
    sk.Write([fn,'-noSync'])
end
sk.Write('external');

% Debug
sk.Check_timing();
nTRs=50;
sk.Plot(nTRs);

if 0
sk.TestReport();
end

%%% tcp control
%host = '134.84.19.164'; %'skope-105T'. if empty, no measurement automation for field monitoring.

%% tcp control
% field monitoring
if exist('host','var')&& ~isempty(host)
    disp('                          ')
    input(['If you like to do realtime field tracking for motion calibration, please go set it up on a different client. \n' ...
        'Press Enter when you are ready to proceed...'],'s');
    
    input(['If you like to do field control, please go set it up on a different client. \n' ...
        'Press Enter when you are ready to proceed...'],'s');

    tcpAutomateScan([fn,'.seq'], host, false);
end

% motion tracking
if exist('host_moco','var')&& ~isempty(host_moco)
    disp('                          ')
    input(['If you like to do motion tracking, please go set it up on a different client. \n' ...
        'Press Enter when you are ready to proceed...'],'s');

    tcpAutomateScan([fn,'.seq'], host_moco, true);
end


