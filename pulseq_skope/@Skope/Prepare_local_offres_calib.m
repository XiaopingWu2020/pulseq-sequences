% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function Prepare_local_offres_calib(this)
% Prepare local off-resonance and position calibration

    %% Set parameters
    this.seq_params.T_flattop = 1100e-3; % [s]
    grad_amp_ = 2.5; % [mT/m]
    this.seq_params.grad_amp = mr.convert(grad_amp_, 'mT/m', 'Hz/m', 'gamma', this.sys.gamma);
    this.seq_params.T_trig_delay = 990e-3; % trigger delay [s]
    this.seq_params.T_inter = 100e-3; % delay between event-blocks [s]

    nRepeats= this.seq_params.nRepeats;
    TR= this.seq_params.TR;

    %% Prepare event objects and combine to eventblocks
    area_flattop = this.seq_params.grad_amp * this.seq_params.T_flattop; 

    mr_trig = mr.makeDigitalOutputPulse(this.seq_params.trigChannel,'duration', this.sys.gradRasterTime, 'delay', this.seq_params.T_trig_delay);
    mr_inter = mr.makeDelay(this.seq_params.T_inter);
      
    mr_noG = mr.makeTrapezoid('x', 'FlatTime', this.seq_params.T_flattop, 'FlatArea', 0); % could be replaced by delay
    mr_Gx = mr.makeTrapezoid('x', 'FlatTime', this.seq_params.T_flattop, 'FlatArea', (-1)*area_flattop);
    mr_Gy = mr.makeTrapezoid('y', 'FlatTime', this.seq_params.T_flattop, 'FlatArea', area_flattop);
    mr_Gz = mr.makeTrapezoid('z', 'FlatTime', this.seq_params.T_flattop, 'FlatArea', area_flattop);

    TRfill= TR- (mr.calcDuration(mr_noG,mr_trig)+ mr_inter.delay+ ...
        mr.calcDuration(mr_Gx,mr_trig)+ mr_inter.delay+ ...
        mr.calcDuration(mr_Gy,mr_trig)+ mr_inter.delay+ ...
        mr.calcDuration(mr_Gz,mr_trig)+ mr_inter.delay);
    TRfill= this.Rasterize(TRfill, "grad");
    mr_TRfill= mr.makeDelay(TRfill);

    this.seq_params.nDynamics= 4* nRepeats;
    this.seq_params.skopeAcqDuration= 50e-3;
    this.seq_params.skopeInterleaveTR= 200e-3;
    this.seq_params.nInterleaves= 1;
    this.seq_params.nPrescans= 0;

    T_tot = 0; % counter: total time
    for idx=1: nRepeats
        % no-grad
        this.seq.addBlock(mr_trig, mr_noG);
        this.seq.addBlock(mr_inter);
        T_tot = T_tot + mr.calcDuration(mr_noG,mr_trig) + mr_inter.delay;


        this.seq.addBlock(mr_trig, mr_Gx);
        this.seq.addBlock(mr_inter);
        T_tot = T_tot + mr.calcDuration(mr_Gx,mr_trig) + mr_inter.delay;

        this.seq.addBlock(mr_trig, mr_Gy);
        this.seq.addBlock(mr_inter);
        T_tot = T_tot + mr.calcDuration(mr_Gy,mr_trig) + mr_inter.delay;

        this.seq.addBlock(mr_trig, mr_Gz);
        this.seq.addBlock(mr_inter);
        T_tot = T_tot + mr.calcDuration(mr_Gz,mr_trig) + mr_inter.delay;

        this.seq.addBlock(mr_TRfill);
        T_tot = T_tot + mr_TRfill.delay;

    end

    this.seq_params.TR= TR;
    assert(nRepeats * this.seq_params.TR - T_tot < 1e-12);

    % Required by PulSeq IDEA
    this.seq.addBlock(this.Make_dummy());
end 