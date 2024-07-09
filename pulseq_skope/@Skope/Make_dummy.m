% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function dummy = Make_dummy(~)
% Prepare a gradient with 0 amplitude
% to create a [SHAPES] paragraph required by PulSeq IDEA
dummy_waveform_ = zeros(1, 2);
dummy = mr.makeArbitraryGrad('x', dummy_waveform_);      
end 