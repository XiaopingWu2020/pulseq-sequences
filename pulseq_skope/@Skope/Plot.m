% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function Plot(this, nTRs)
% Default PulSeq plot

if isfield(this.seq_params,'TR')
    this.seq.plot('timeRange', [0 nTRs]*this.seq_params.TR, 'TimeDisp', 'ms', 'label', 'lin');
else
    this.seq.plot('timeRange', [0 nTRs]*this.seq_params.total_time, 'TimeDisp', 'ms', 'label', 'lin');
end

end