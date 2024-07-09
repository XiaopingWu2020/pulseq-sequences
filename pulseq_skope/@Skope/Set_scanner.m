%% Set scanner parameters
% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function [sys, scanner] = Set_scanner(~, scanner_type)
    scanner = [];
    switch scanner_type
        case "Siemens Terra 7T SC72CD"
            scanner.type = scanner_type;
            
            scanner.maxGrad = 70;
            scanner.maxGrad_unit = 'mT/m';

            scanner.maxSlew = 200;
            scanner.maxSlew_unit = 'T/m/s';

            scanner.B0 = 7;
            scanner.effProbeLifetime= 36e-3; % it was found with cranberry that skope acq durations longer than this value would result in noisy field measurement especially toward the tail of the kspace trajectory. 
            scanner.T1GM= 1.78; % estimated gray matter T1 in sec

        case "Siemens dotplus 10.5T SC72CD"
            scanner.type = scanner_type;
            
            scanner.maxGrad = 70;
            scanner.maxGrad_unit = 'mT/m';

            scanner.maxSlew = 200;
            scanner.maxSlew_unit = 'T/m/s';

            scanner.B0 = 10.5;
            scanner.effProbeLifetime= 30e-3; %43e-3; % found for our 16 clip on field clips.
            scanner.T1GM= 2.1; % estimated gray matter T1 in sec

        case "Siemens dotplus 10.5T HG" % assuming more powerful GPA is in use
            scanner.type = scanner_type;

            scanner.maxGrad = 175;
            scanner.maxGrad_unit = 'mT/m';

            scanner.maxSlew = 900;
            scanner.maxSlew_unit = 'T/m/s';

            scanner.B0 = 10.5;
            scanner.effProbeLifetime= 30e-3; %43e-3; % found for our 16 clip on field clips.
            scanner.T1GM= 2.1; % estimated gray matter T1 in sec

        case "Siemens dotplus 10.5T HG LP" %assuming current body gradient GPA is in use
            scanner.type = scanner_type;

            scanner.maxGrad = 117;
            scanner.maxGrad_unit = 'mT/m';

            scanner.maxSlew = 800; 
            scanner.maxSlew_unit = 'T/m/s';

            scanner.B0 = 10.5;
            scanner.effProbeLifetime= 30e-3; %43e-3; % found for our 16 clip on field clips.
            scanner.T1GM= 2.1; % estimated gray matter T1 in sec

        case "Canon Galan 3T"
            % TODO

        case "Philips XX 3T"
            % TODO

        case "Philips XX 7T"
            % TODO
    end
    
    sys = mr.opts(  'MaxGrad',  scanner.maxGrad,  'GradUnit', scanner.maxGrad_unit, ...
                    'MaxSlew',  scanner.maxSlew,  'SlewUnit', scanner.maxSlew_unit, ...
                   'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6, ...
                   'B0', scanner.B0, 'effProbeLifetime', scanner.effProbeLifetime, ...
                   'T1GM', scanner.T1GM);
end