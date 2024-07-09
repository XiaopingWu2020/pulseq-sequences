%% Probe Environment Adjustment
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

close all; clc; clear all;
addpath(genpath('../pulseq/')); % Pulseq path

%% Scanner specs
seq_name = "local_ec_calib";          
seq_type = "Local Eddy-current Calibration";
scanner_type = "Siemens Terra 7T SC72CD";

sk = Skope(seq_name, seq_type, scanner_type);
seq_params = struct('N_blips_axis', 32, 'trigChannel','osc0'); % number of blips per excitation and axis
sk.Prepare(seq_params);
sk.Write('ec_calib');

% Debug
sk.Check_timing();
sk.Plot(0.2);

