function Prepare_grad_tones(this)
% Prepare gradient tones for motion tracking. 
% M. Haeberlin, et al. MRM 2015

    %% Set parameters
    % this.seq_params.axis = ['x', 'y', 'z'];   
    TR= this.seq_params.TR; % sec
    gradFreeDelay= 200e-6; %this.seq_params.gradFreeDelay; % [s] 

    % tone specs
    tw= this.seq_params.gradToneLength; %1e-3; %4.8e-3; % length of tones in second. 
    tw= this.Rasterize(tw, "grad");
    dt= this.sys.gradRasterTime;

    f0= 1e3*[6, 7, 8]; % nominal tone frequencies [fx, fy, fz] in hz
    A0= 4;%3.71; % constant max tone amplitude for gx, gy, and gz in mT/m. 
    A = mr.convert(A0, 'mT/m', 'Hz/m', 'gamma', this.sys.gamma);

    f= round(tw* f0)./tw; % frequency adjustment to ensure harmonics across the duration of tones. 
    t_vec= dt* (0:round(tw/dt)); 

    % grad tones
    gx= A* sin(2*pi*f(1)* t_vec);
    gy= A* sin(2*pi*f(2)* t_vec);
    gz= A* sin(2*pi*f(3)* t_vec);

    gt.gx= gx;
    gt.gy= gy;
    gt.gz= gz;
    gt.t_vec= t_vec;    
    gt.unit= 'Hz/m';
    gt.gamma= this.sys.gamma;
    gt.A0= A0;
    gt.f0= f0;
    save gradtones_naive gt;
    disp('-> gradient tones waveforms saved in gradtones.mat...')
    
    % todo: check gradient tones to ensure they do not violate slew rate
    % limits

    %
    mr_gx = mr.makeArbitraryGrad('x',(-1)* gx,this.sys);
    mr_gy = mr.makeArbitraryGrad('y',gy,this.sys);
    mr_gz = mr.makeArbitraryGrad('z',gz,this.sys);

    this.seq_params.nDynamics= this.seq_params.N_rep;
    this.seq_params.skopeAcqDuration= gradFreeDelay+ mr.calcDuration(mr_gx,mr_gy,mr_gz)+ 1e-3;
    this.seq_params.skopeInterleaveTR= TR- 100e-6;
    this.seq_params.nInterleaves= 1;
    this.seq_params.nPrescans= 0;

    %% Prepare event objects
    % Trigger
    mr_trig = mr.makeDigitalOutputPulse(this.seq_params.trigChannel,'duration', this.sys.gradRasterTime);
    mr_gradFreeDelay = mr.makeDelay(gradFreeDelay);

    % TRfill
    TRfill= TR- mr_gradFreeDelay.delay- mr.calcDuration(mr_gx, mr_gy, mr_gz);
    TRfill= this.Rasterize(TRfill, "grad");
    mr_TRfill= mr.makeDelay(TRfill);
    
    %% Combines event objects to eventblocks
    T_tot = 0; % counter: total time

    for i= 1:this.seq_params.N_rep
            % trigger block
            this.seq.addBlock(mr_trig, mr_gradFreeDelay);
            T_tot = T_tot + mr_gradFreeDelay.delay;

            % gradient tones
            this.seq.addBlock(mr_gx, mr_gy, mr_gz);
            T_tot = T_tot + mr.calcDuration(mr_gx, mr_gy, mr_gz);

            % TRfill
            this.seq.addBlock(mr_TRfill);
            T_tot = T_tot + mr_TRfill.delay;
    end

    this.seq_params.TR= TR;
    assert(this.seq_params.N_rep * this.seq_params.TR - T_tot < 1e-12);

    % Required by PulSeq IDEA
    this.seq.addBlock(this.Make_dummy());
end 