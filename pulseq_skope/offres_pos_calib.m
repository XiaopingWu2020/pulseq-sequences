%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

close all; clc; clear all;
addpath(genpath('../pulseq/')); % Pulseq path

%% Scanner specs
seq_name = "local_offres_calib";          
seq_type = "Local Off-res Calibration";
scanner_type = "Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);
seq_params = struct('trigChannel', 'osc0'); % 'osc0','osc1','ext1'
seq_params.TR= 5; % s
seq_params.nRepeats= 10;%20; %1;

sk.Prepare(seq_params);

fn= 'offres_calib';
if sk.seq_params.nRepeats>1
    fn= [fn,'-x',num2str(sk.seq_params.nRepeats)];
end

sk.Write(fn);

% Debug
sk.Check_timing();
sk.Plot(10);

