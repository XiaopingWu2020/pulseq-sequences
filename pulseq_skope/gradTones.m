%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

close all; clc; clear all;
addpath(genpath('../pulseq/')); % Pulseq path

%% Scanner specs
seq_name = "gradTones";          
seq_type = "Gradient Tones";
scanner_type = "Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);
seq_params = struct('trigChannel', 'osc0'); % 'osc0','osc1','ext1'
seq_params.TR = 200e-3; % [s], set to > 110 ms which is the skope minTR.     
seq_params.N_rep = 1500;%3000;%20; % number of repeats. 
seq_params.gradToneLength= 1e-3; % length of tones in second. 1 ms for 10.5T motion probes. 

sk.Prepare(seq_params);
sk.Write('gradTones-x1500');

% Debug
sk.Check_timing();
sk.Plot(2);

