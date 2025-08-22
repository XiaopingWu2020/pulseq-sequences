function Prepare_girf_calib(this)
% Prepare gradient impulse response function (GIRF) calibration
% S. Vannesjo, et al. MRM 2013

    %% Set parameters
    this.seq_params.axis = ['x', 'y', 'z'];   
    this.seq_params.delay_trig2blip = 5e-3; % [s] delay between trigger and blip
    nRepeats= this.seq_params.nRepeats;
    
    %this.seq_params.relax = 0.9; % relax system limits to provide 180 T/m/s slew rate. 
    gradSafetyMargin= this.seq_params.gradSafetyMargin;

    % Compliance to gradient raster time
    % for time2peak = 1e-6*(50:10:160); % [s]
    this.seq_params.N_blips_axis= 12;
    dur_base_ = 100e-6; % [s]; blip duration (factor must be even)
    dur_incr_ = 20e-6; % [s]; increment to adjust blip duration (factor must be even)
    this.seq_params.dur_base = this.Rasterize(dur_base_, "grad");
    this.seq_params.dur_incr = this.Rasterize(dur_incr_, "grad");
    assert(~mod((this.seq_params.dur_base / 2), this.sys.gradRasterTime));
    assert(~mod((this.seq_params.dur_incr / 2), this.sys.gradRasterTime));

    if (mod(this.seq_params.N_blips_axis, 2) ~= 0)
        error("Choose even number of blips per axis 'N_blips_axis'!"); % balancing condition
    end

    this.seq_params.nDynamics= length(this.seq_params.axis)* this.seq_params.N_blips_axis* nRepeats;
    this.seq_params.skopeAcqDuration= 50e-3;
    this.seq_params.skopeInterleaveTR= this.seq_params.TR- 100e-3;
    this.seq_params.nInterleaves= 1;
    this.seq_params.nPrescans= 0;

    %% Prepare event objects
    % Trigger
    mr_trig = mr.makeDigitalOutputPulse(this.seq_params.trigChannel,'duration', this.sys.gradRasterTime);
    
    % Delay
    mr_delay_trig2blip = mr.makeDelay(this.seq_params.delay_trig2blip);

    % Blips
    mr_blips = {};
    sign= 1;
    for ax = this.seq_params.axis
        for idx= 0:(this.seq_params.N_blips_axis-1)

            dur = this.seq_params.dur_base + idx * this.seq_params.dur_incr;

            % use relaxed slew-rate
            mr_blip_ = this.Make_blip(ax, this.scanner.maxSlew * gradSafetyMargin   , this.scanner.maxSlew_unit, dur, sign);
            assert(dur - mr.calcDuration(mr_blip_) < 1e-15); % sanity check

            if this.seq_params.useArbitraryGrad
                blip_= this.createTrapGrad(mr_blip_);
                mr_blip_ = mr.makeArbitraryGrad(ax,blip_,this.sys); % so that we know exactly what input waveform is applied.
            end
            
            mr_blips = [mr_blips(:)', {mr_blip_}];
        end
    end
   
    if this.seq_params.useArbitraryGrad
        save girf_calib_arbGrad mr_blips
    else
        save girf_calib mr_blips
    end
   
    %% Combines event objects to eventblocks
    T_tot = 0; % counter: total time

    for i= 1:nRepeats
        for mr_blip_=mr_blips

            % trigger block
            this.seq.addBlock(mr_trig);
            T_tot = T_tot + mr_trig.duration;

            % delay
            this.seq.addBlock(mr_delay_trig2blip);
            T_tot = T_tot + mr_delay_trig2blip.delay;

            this.seq.addBlock(mr_blip_);
            T_tot = T_tot + mr.calcDuration(mr_blip_);

            % delay
            iDelay= this.seq_params.TR - ...
                (mr_trig.duration+ mr_delay_trig2blip.delay + mr.calcDuration(mr_blip_));
            iDelay= this.Rasterize(iDelay, "grad");
            mr_delay = mr.makeDelay(iDelay);
            this.seq.addBlock(mr_delay);
            T_tot = T_tot + mr_delay.delay;
        end
    end

    assert(nRepeats * this.seq_params.TR - T_tot < 1e-15);

    % Required by PulSeq IDEA
    this.seq.addBlock(this.Make_dummy());
end 