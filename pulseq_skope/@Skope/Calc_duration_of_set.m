% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function t = Calc_duration_of_set(~, set)
% Calculate the total duration of all event objects contained in a set
    
    t = 0;
    
    for i = set
        t = t + mr.calcDuration(i);
    end      
end 