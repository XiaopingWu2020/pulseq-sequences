%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

clearvars;
addpath(genpath(fullfile('..','pulseq')));
%addpath('~/matlab/skope_TCP_client/methods')

%% Scanner specs
seq_name = "GIRF_calib";          
seq_type = "GIRF Calibration";
scanner_type = "Siemens Terra 7T SC72CD"; % "Siemens dotplus 10.5T HG LP"; % "Siemens dotplus 10.5T SC72CD"; Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);
seq_params = struct('trigChannel', 'osc0'); % 'osc0','osc1','ext1'
seq_params.TR = 2;%4; % [s]    
seq_params.nRepeats = 25;%50; % number of repeats. 
seq_params.useArbitraryGrad= false; % false to use mr.makeTrapezoid and true to use mr.makeArbitraryGrad

seq_params.gradSafetyMargin= 0.9;
%%% tcp control
% host = '134.84.19.164'; %'skope-105T'. if empty, no measurement automation for field monitoring.
%host_moco = '134.84.19.37'; %'skope-terra'. if empty, no measurement automation for motion detection.

%%
sk.Prepare(seq_params);

fn= 'girf_calib';
if sk.seq_params.useArbitraryGrad
    fn= [fn, '_arbGrad'];
end
sk.Write(fn);
sk.Write('external');

% Debug
sk.Check_timing();
sk.Plot(36);

% %% tcp control
% % field monitoring
% if exist('host','var')&& ~isempty(host)
%     tcpAutomateScan([fn,'.seq'], host, false);
% end


