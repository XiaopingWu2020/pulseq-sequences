function Prepare_2dGre(this)
% Prepare 2D GRE sequence
% this is based on writeGradientEcho_label.m from pulseq.
%
this.seq_params.resolution= this.seq_params.fov/ this.seq_params.N;
fov = this.seq_params.fov;
Nx= this.seq_params.N; Ny= Nx;
%alpha= this.seq_params.alpha;
thickness= this.seq_params.thickness;
Nslices= this.seq_params.Nslices;
TEs= this.seq_params.TE; TE= TEs(1);
TR= this.seq_params.TR;
gradFreeDelay= this.seq_params.gradFreeDelay;
sliceGap= this.seq_params.sliceGap;
partialFourier= this.seq_params.partialFourier;  
if partialFourier< 0.5
    error('-> Partial Fourier Factor cannot be less than 0.5 !!!')
else
    partFourierFactor= 2*(partialFourier-0.5);
end
reordering= this.seq_params.reordering;
Nreps = this.seq_params.nRepeats;

gradSafetyMargin= this.seq_params.gradSafetyMargin;
this.sys.maxGrad= gradSafetyMargin.* this.sys.maxGrad;
this.sys.maxSlew= gradSafetyMargin.* this.sys.maxSlew;

rfSpoilingInc=117;              % RF spoiling increment
readoutTime= this.seq_params.readoutTime; %1e-3; %3.2e-3;              % ADC duration
%
ro_os=2; %1;                   % oversampling factor (in contrast to the product sequence we don't really need it)
% Trigger
mr_trig = mr.makeDigitalOutputPulse(this.seq_params.trigChannel,'duration', this.sys.gradRasterTime);
mr_gradFreeDelay = mr.makeDelay(gradFreeDelay);

% Create alpha-degree slice selection pulse and gradient
alpha= 1;
[mr_rf, mr_gz] = mr.makeSincPulse(alpha*pi/180,'Duration',3e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'system',this.sys);

% Define other gradients and ADC events
Ny_pre=floor(partFourierFactor*Ny/2-1); % PE steps prior to ky=0, excluding the central line
Ny_post=round(Ny/2+1); % PE lines after the k-space center including the central line
Ny_meas=Ny_pre+Ny_post;
deltak=1/fov;
mr_gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',readoutTime,'system',this.sys);

% actual_area=mr_gx.area-mr_gx.amplitude/mr_gx.riseTime*blip_dur/2*blip_dur/2/2-mr_gx.amplitude/mr_gx.fallTime*blip_dur/2*blip_dur/2/2;
% mr_gx.amplitude= mr_gx.amplitude/actual_area*kWidth;
% mr_gx.area = mr_gx.amplitude*(mr_gx.flatTime + mr_gx.riseTime/2 + mr_gx.fallTime/2);
% mr_gx.flatArea = mr_gx.amplitude*mr_gx.flatTime;

% calculate mr_adc
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be qite different from Nx and
% readoutTime/Nx, respectively.
adcDwellNyquist=deltak/mr_gx.amplitude/ro_os;
% round-down dwell time to 100 ns
adcDwell=floor(adcDwellNyquist*1e7)*1e-7;
adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of mr_adc samples need to be divisible by 4
% MZ: no idea, whether ceil,round or floor is better for the adcSamples...
mr_adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',mr_gx.riseTime);
% realign the mr_adc with respect to the gradient
time_to_center=mr_adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
mr_adc.delay=round((mr_gx.riseTime+mr_gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
% this rounding actually makes the sampling points on odd and even readouts
% to appear misalligned. However, on the real hardware this misalignment is
% much stronger anyways due to the grdient delays

% negate x gradient to be consistent with IDEA
if mr_gx.channel=='x'
    mr_gx.amplitude= (-1)* mr_gx.amplitude;
    mr_gx.area= (-1)* mr_gx.area;
    mr_gx.flatArea= (-1)* mr_gx.flatArea;
end

mr_gxPre = mr.makeTrapezoid('x','Area',-mr_gx.area/2,'Duration',1e-3,'system',this.sys);
mr_gzReph = mr.makeTrapezoid('z','Area',-mr_gz.area/2,'Duration',1e-3,'system',this.sys);
phaseAreas0 = ((0:Ny-1)-Ny/2)*deltak; % phase area should be Kmax for clin=0 and -Kmax for clin=Ny... strange
gyPre = mr.makeTrapezoid('y','Area',max(abs(phaseAreas0(:))),'Duration',mr.calcDuration(mr_gxPre),'system',this.sys);

% assuming upper kspace is omitted if there is partial Fourier
switch reordering
    case 'ascending'
        sampleIdx= 0:(Ny_meas-1);
    case 'descending'
        sampleIdx= [0 (Ny_meas-1):-1:1];
    otherwise
        error('-> Specified reordering not supported...')
end
phaseAreas= phaseAreas0(sampleIdx+ 1);

% gradient spoiling
mr_gxSpoil=mr.makeTrapezoid('x','Area', 2*Nx*deltak,'system',this.sys);
mr_gzSpoil=mr.makeTrapezoid('z','Area', 4/thickness,'system',this.sys);

% Calculate timing
minTE= mr_gz.fallTime + mr_gz.flatTime/2 + mr.calcDuration(mr_gzReph) ...
    + gradFreeDelay + mr.calcDuration(mr_gxPre)+ mr.calcDuration(mr_gx)/2;
if TE< minTE
    TE= minTE+ this.sys.gradRasterTime;
    disp('-> TE too short. Minimum allowed TE is used.')
end
delayTE = this.Rasterize(TE- minTE, "grad");
%
nContrasts= length(TEs);
if nContrasts> 1
    dTE= TEs(2)- TEs(1);
else
    dTE= 0;
end
if dTE< mr.calcDuration(mr_gx)
    disp('-> dTE too short. Increased to minimum required dTE...')
    dTE= mr.calcDuration(mr_gx)+ this.sys.gradRasterTime;
end
delayDTE= this.Rasterize(dTE- mr.calcDuration(mr_gx),'grad');

minTR1= this.seq_params.skopeMinTR+ 100e-6;
minExcitTR= mr.calcDuration(mr_gz)/2+ TE+ ( (nContrasts-1) + 0.5 )*mr.calcDuration(mr_gx) ...
    + (nContrasts-1)*delayDTE ...
    + mr.calcDuration(mr_gxSpoil,gyPre,mr_gzSpoil);
minTR2= minExcitTR* Nslices;
minTR= max(minTR1, minTR2);
if TR< minTR
    TR= minTR;
    disp('-> TR too short. Minimum allowed TR is used.')
end
delayTR= this.Rasterize(TR/Nslices- minExcitTR, "grad");
if delayTR==0
    delayTR= this.sys.gradRasterTime;
end
excitTR= minExcitTR+ delayTR;
TR= excitTR* Nslices;
this.seq_params.skopeInterleaveTR= TR- 100e-6;
if excitTR >= minTR1
    this.seq_params.skopeInterleaveTR= excitTR- 100e-6;
end

this.seq_params.TR= TR;
TEs= TE+ dTE*(0: (nContrasts-1));
this.seq_params.TE= TEs;
this.seq_params.nContrasts= nContrasts;
this.seq_params.echoSpacing= dTE;

%%% calculate Ernst angle per TR and T1
T1= this.sys.T1GM;
TR= this.seq_params.TR;
alpha= acosd(exp(-TR./ T1));
mr_rf.signal= alpha.* mr_rf.signal;
this.seq_params.alpha= alpha;

%
trig2AdcTime= gradFreeDelay+ mr.calcDuration(mr_gxPre)+ mr_adc.delay;
skopeAcqDuration= trig2AdcTime+ mr.calcDuration(mr_adc) + (nContrasts-1)*mr.calcDuration(mr_gx) ...
    + (nContrasts-1)*delayDTE + 1e-3;
this.seq_params.trigger2AdcTime= trig2AdcTime;
this.seq_params.skopeAcqDuration= skopeAcqDuration;
this.seq_params.nDynamics= ceil(Nreps* Nslices* Ny_meas / (this.seq_params.nInterleaves));
%minTR_sync= skopeAcqDuration+ 0.5e-3;

this.seq_params.adcDwellTime= mr_adc.dwell; %


rf_phase=0;
rf_inc=0;

% all LABELS / counters and flags are automatically initialized to 0 in the beginning, no need to define initial 0's
% so we will just increment LIN after the ADC event (e.g. during the spoiler)

%this.seq.addBlock(mr.makeLabel('SET','REV', 1)); % left-right swap fix (needed for 1.4.0)

%this.seq.addBlock(mr.makeDelay(1)); % older scanners like Trio may need this
% dummy delay to keep up with timing

% Prescans for sync
if this.seq_params.nPrescans> 0   
    this.seq.addBlock(mr.makeLabel('SET','DUM', 1));
    disp('-> Using native ADC for prescans...')
        syncMinTR= trig2AdcTime + nContrasts*(mr.calcDuration(mr_gx)+ delayDTE);
        syncDelay= this.seq_params.skopeInterleaveTR- syncMinTR+ 100e-6 - 1.5e-3;
        if this.seq_params.doFastestPrescan
            syncDelay= 0;
        end
    
        for i=1: this.seq_params.nPrescans
            this.seq.addBlock(mr.makeLabel('SET','ECO', 0));
            % insert trigger
            this.seq.addBlock(mr_trig, mr_gradFreeDelay);
            %
            this.seq.addBlock(mr_gxPre,gyPre);

            for iTE=1:nContrasts
                %this.seq.addBlock(mr.makeLabel('SET','REV', mr_gx.amplitude<0));
                this.seq.addBlock(mr.makeLabel('SET','REV', mr_gx.amplitude>0));
                this.seq.addBlock(mr_gx,mr_adc);
                this.seq.addBlock(mr.makeDelay(delayDTE));
                mr_gx.amplitude= -mr_gx.amplitude;
                this.seq.addBlock(mr.makeLabel('INC','ECO',1));
            end

            this.seq.addBlock(mr.makeDelay(1.5e-3));

            if syncDelay>0
                this.seq.addBlock(mr.makeDelay(this.Rasterize(syncDelay,'grad')));
            end
        end
        %
        this.seq.addBlock(mr.makeDelay(3)) % 3 sec delay allowing skope to switch back into normal field monitoring mode.
        this.seq.addBlock(mr.makeLabel('SET','DUM', 0));
end

for r=1:Nreps
    %this.seq.addBlock(mr_trigPhy, mr.makeLabel('SET','SLC', 0));
    % loop over phase encodes

    for i=1:Ny_meas
        % loop over slices and define sequence blocks
        this.seq.addBlock(mr.makeLabel('SET','LIN', sampleIdx(i)));
        this.seq.addBlock(mr.makeDelay(delayTR),mr.makeLabel('SET','SLC', 0));

        for s=1:Nslices
            this.seq.addBlock(mr.makeLabel('SET','ECO', 0));

            mr_rf.freqOffset=mr_gz.amplitude*(1+sliceGap)*thickness*(s-1-(Nslices-1)/2);
            mr_rf.phaseOffset=-2*pi*mr_rf.freqOffset*mr.calcRfCenter(mr_rf);% compensate for the slice-offset induced phase

            mr_rf.phaseOffset=rf_phase/180*pi;
            mr_adc.phaseOffset=rf_phase/180*pi;
            rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase=mod(rf_phase+rf_inc, 360.0);
            %
            this.seq.addBlock(mr_rf,mr_gz);
            this.seq.addBlock(mr_gzReph);
            this.seq.addBlock(mr.makeDelay(delayTE));
            % insert trigger
            this.seq.addBlock(mr_trig, mr_gradFreeDelay);
            %
            gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',mr.calcDuration(mr_gxPre),'system',this.sys);
            this.seq.addBlock(mr_gxPre,gyPre);

            for iTE=1:nContrasts
                %this.seq.addBlock(mr.makeLabel('SET','REV', mr_gx.amplitude<0));
                this.seq.addBlock(mr.makeLabel('SET','REV', mr_gx.amplitude>0));
                this.seq.addBlock(mr_gx,mr_adc);
                this.seq.addBlock(mr.makeDelay(delayDTE));
                mr_gx.amplitude= -mr_gx.amplitude;
                this.seq.addBlock(mr.makeLabel('INC','ECO',1));
            end
            if rem(nContrasts,2)~=0
                mr_gx.amplitude = -mr_gx.amplitude;
            end

            gyPre.amplitude=-gyPre.amplitude;
            this.seq.addBlock(mr_gxSpoil,gyPre,mr_gzSpoil);

            this.seq.addBlock(mr.makeDelay(delayTR),mr.makeLabel('INC','SLC', 1));

        end

    end
    this.seq.addBlock(mr.makeLabel('INC','REP', 1)); %
end
% Required by PulSeq IDEA
% this.seq.addBlock(this.Make_dummy());

%
% %% plot sequence and k-space diagrams
%
% seq.plot('timeRange', [0 32]*TR, 'TimeDisp', 'ms', 'label', 'lin');
%
% % k-space trajectory calculation
% [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%
% % plot k-spaces
% figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
% hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
% title('k-space components as functions of time');
% figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
% axis('equal'); % enforce aspect ratio for the correct trajectory display
% hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
% title('2D k-space');
%

end