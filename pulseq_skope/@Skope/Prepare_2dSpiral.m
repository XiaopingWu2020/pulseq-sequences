function Prepare_2dSpiral(this)
% Prepare 2D spiral sequence
% this is based on writeSpiral.m from pulseq.
%

this.seq_params.isSpiral = true;
this.seq_params.resolution= this.seq_params.fov/ this.seq_params.N;
fov = this.seq_params.fov/ this.seq_params.accelerationFactor;
N= this.seq_params.N/ this.seq_params.accelerationFactor;
alpha= this.seq_params.alpha;
thickness= this.seq_params.thickness;
Nslices= this.seq_params.Nslices;
TE= this.seq_params.TE;
gradFreeDelay= this.seq_params.gradFreeDelay;
sliceGap= this.seq_params.sliceGap;
safetyMargin= this.seq_params.gradSafetyMargin;

Nreps = this.seq_params.nRepeats;
adcSamplesPerSegment= this.seq_params.maxAdcSegmentLength;

Oversampling= 2; % by looking at the periphery of the spiral I would say it needs to be at least 2
phi= 0;%pi/2; % orientation of the readout e.g. for interleaving

% Create fat-sat pulse
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6* this.sys.B0* this.sys.gamma;
mr_rfFatSat = mr.makeGaussPulse(110*pi/180,'system',this.sys,'Duration',8e-3,'dwell',10e-6,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
mr_rfFatSat.phaseOffset=-2*pi*mr_rfFatSat.freqOffset*mr.calcRfCenter(mr_rfFatSat); % compensate for the frequency-offset induced phase
mr_gzFatSat = mr.makeTrapezoid('z',this.sys,'delay',mr.calcDuration(mr_rfFatSat),'Area',0.1/1e-4); % spoil up to 0.1mm

% Create slice selection pulse and gradient
[mr_rf, mr_gz, mr_gzReph] = mr.makeSincPulse(alpha*pi/180,'system',this.sys,'Duration',2e-3,...
    'SliceThickness',thickness,'apodization',0.42,'timeBwProduct',4,'use','excitation');

% define the trigger to play out
%mr_trigPhy=mr.makeTrigger('physio1','duration', 2000e-6); % duration after
mr_trig = mr.makeDigitalOutputPulse(this.seq_params.trigChannel,'duration', this.sys.gradRasterTime);
mr_gradFreeDelay = mr.makeDelay(gradFreeDelay);

% Define gradients and ADC events
deltak=1/fov;
kRadius = round(N/2);
kSamples=round(2*pi*kRadius)*Oversampling;
adcSamplesDesired=kRadius*kSamples;
spiral_grad_shape= this.designArchimedean(deltak, kRadius, kSamples, safetyMargin);

% calculate ADC
% round-down dwell time to 100 ns
adcTime = this.sys.gradRasterTime*size(spiral_grad_shape,2);
% actually it is trickier than that: the (Siemens) interpreter sequence
% per default will try to split the trajectory into segments <=1000 samples
% and every of these segments will have to have duration aligned to the
% gradient raster time
adcSegments=round(adcSamplesDesired/adcSamplesPerSegment);
adcSamples=adcSegments*adcSamplesPerSegment;
adcDwell=round(adcTime/adcSamples/100e-9)*100e-9; % on Siemens adcDwell needs to be aligned to 100ns (if my memory serves me right)
adcSegmentDuration=adcSamplesPerSegment*adcDwell; % with the 1000 samples above and the 100ns alignment we automatically fullfill the segment alignment requirement
% update segment count
adcSegments=floor(adcTime/adcSegmentDuration);
if adcSegments> 128
    adcSamplesPerSegment= ceil(adcSegments* adcSamplesPerSegment/ 128/ 100)*100;
    
    adcSegments=round(adcSamplesDesired/adcSamplesPerSegment);
    adcSamples=adcSegments*adcSamplesPerSegment;
    adcDwell=round(adcTime/adcSamples/100e-9)*100e-9; % on Siemens adcDwell needs to be aligned to 100ns (if my memory serves me right)
    adcSegmentDuration=adcSamplesPerSegment*adcDwell; % with the 1000 samples above and the 100ns alignment we automatically fullfill the segment alignment requirement
    % update segment count
    adcSegments=floor(adcTime/adcSegmentDuration);

    this.seq_params.maxAdcSegmentLength = adcSamplesPerSegment; %
    fprintf('NOTE: ADC segment length too short, increased to %d \n', this.seq_params.maxAdcSegmentLength);
end
if mod(adcSegmentDuration, this.sys.gradRasterTime)>eps
    error('ADC segmentation model results in incorrect segment duration');
end

adcSamples=adcSegments*adcSamplesPerSegment;
mr_adc = mr.makeAdc(adcSamples,'Dwell',adcDwell);%

this.seq_params.adcSegmentDuration= adcSegmentDuration; 
this.seq_params.adcLength= adcSamples;
this.seq_params.nAdcSegments= adcSegments; %
fprintf('Number of ADC segments: %d \n', this.seq_params.nAdcSegments);
this.seq_params.adcDwellTime= adcDwell; %
fprintf('ADC dwell time (us): %f \n', 1e6*this.seq_params.adcDwellTime);
this.seq_params.effectiveAdcTime= adcSamples* adcDwell/ adcTime;
fprintf('Effective ADC time (ADC time/gradient duration): %f \n', this.seq_params.effectiveAdcTime);


% extend spiral_grad_shape by repeating the last sample
% this is needed to accomodate for the ADC tuning delay
spiral_grad_shape = [spiral_grad_shape spiral_grad_shape(:,end)];

grad0.shape= spiral_grad_shape;
grad0.dt= this.sys.gradRasterTime;
grad0.unit= 'Hz/m';
grad0.gamma= this.sys.gamma;
fn= ['xw_sp2d-',num2str(1e3*this.seq_params.resolution,2),'mm-r',num2str(this.seq_params.accelerationFactor)];
save([fn,'.mat'],'grad0')

% readout grad
mr_gx = mr.makeArbitraryGrad('x',(-1)* spiral_grad_shape(1,:),this.sys);
mr_gy = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:),this.sys);

% spoilers
mr_gzSpoil=mr.makeTrapezoid('z',this.sys,'Area',N*4/fov);
mr_gxSpoil=mr.makeExtendedTrapezoid('x','times',[0 mr.calcDuration(mr_gzSpoil)],'amplitudes',(-1)*[spiral_grad_shape(1,end),0]); %todo: make a really good spoiler
mr_gySpoil=mr.makeExtendedTrapezoid('y','times',[0 mr.calcDuration(mr_gzSpoil)],'amplitudes',[spiral_grad_shape(2,end),0]); %todo: make a really good spoiler

% because of the ADC alignment requirements the sampling window possibly
% extends past the end of the trajectory (these points will have to be
% discarded in the reconstruction, which is no problem). However, the
% ramp-down parts and the Z-spoiler now have to be added to the readout
% block otherwise there will be a gap inbetween
% mr_gzSpoil.delay=mr.calcDuration(mr_gx);
% mr_gxSpoil.delay=mr_gzSpoil.delay;
% mr_gySpoil.delay=mr_gzSpoil.delay;
% mr_gxCombined=mr.addGradients({mr_gx,mr_gxSpoil}, this.sys);
% mr_gyCombined=mr.addGradients({mr_gy,mr_gySpoil}, this.sys);
%mr_gzCombined=mr.addGradients([gzReph,mr_gzSpoil], this.sys);


% Calculate timing
minTE= mr_gz.fallTime + mr_gz.flatTime/2 + mr.calcDuration(mr_gzReph)+ gradFreeDelay;
if TE< minTE
    TE= minTE;
    disp('-> TE too short. Minimum allowed TE is used.')
end
delayTE = this.Rasterize(TE- minTE, "grad");
if delayTE==0
    delayTE= this.sys.gradRasterTime;
end
TE= minTE+ delayTE;
this.seq_params.TE= TE;

minExcitTR= mr.calcDuration(mr_rfFatSat)+ mr.calcDuration(mr_gzFatSat) ...
    + mr.calcDuration(mr_gz)/2+ TE+ mr.calcDuration(mr_gx,mr_gy,mr_adc) ...
    + mr.calcDuration(mr_gzSpoil,mr_gxSpoil,mr_gySpoil);
delayTR= this.calculateTiming(minExcitTR);
%
trig2AdcTime= gradFreeDelay;
skopeAcqDuration= trig2AdcTime + mr.calcDuration(mr_gx,mr_gy,mr_adc) + 1e-3;
this.seq_params.trigger2AdcTime= trig2AdcTime;
this.seq_params.skopeAcqDuration= skopeAcqDuration;
this.seq_params.nDynamics= ceil(Nreps* Nslices / (this.seq_params.nInterleaves));

if this.seq_params.skopeAcqDuration> this.sys.effProbeLifetime
    fprintf('-> NOTE: skope acq duration is longer than effective probe lifetime of %d ms...\n', 1e3* this.sys.effProbeLifetime);
    disp('-> Please consider adjusting imaging parameters (eg, decrease resolution or increase undersampling rate')
    disp('-> to avoid noisy field measurement toward the tail of the kspace trajectory...')
end
    
% Prescans for sync
if this.seq_params.nPrescans>0
    
    disp('=> Adding sync scans...')
    this.seq.addBlock(mr.makeLabel('SET','DUM', 1)); % set dummy scans

    if this.seq_params.useSingleAdcSegment4Sync
        disp('-> Using a single ADC segment for prescans...')
        mr_adc4sync = mr.makeAdc(adcSamplesPerSegment,'Dwell',adcDwell);
        nGradSamples= round(adcSegmentDuration/ this.sys.gradRasterTime)+ 1;
        mr_gx4sync = mr.makeArbitraryGrad('x',(-1)*spiral_grad_shape(1,1:nGradSamples),this.sys);
        mr_gy4sync = mr.makeArbitraryGrad('y',spiral_grad_shape(2,1:nGradSamples),this.sys);

        syncMinTR= trig2AdcTime + mr.calcDuration(mr_gx4sync,mr_gy4sync,mr_adc4sync);
        syncDelay= this.seq_params.skopeInterleaveTR- syncMinTR+ 100e-6 - 1.5e-3;
        if this.seq_params.doFastestPrescan
            syncDelay= 0;
        end

        for i=1: this.seq_params.nPrescans
            % insert trigger
            this.seq.addBlock(mr_trig, mr_gradFreeDelay);
            %
            this.seq.addBlock(mr_gx4sync,mr_gy4sync,mr_adc4sync);

            this.seq.addBlock(mr.makeDelay(1.5e-3));

            if syncDelay>0
                this.seq.addBlock(mr.makeDelay(this.Rasterize(syncDelay,'grad')));
            end
        end
    else % use native ADC
        disp('-> Using native ADC for prescans...')
        syncMinTR= trig2AdcTime + mr.calcDuration(mr_gx,mr_gy,mr_adc);
        syncDelay= this.seq_params.skopeInterleaveTR- syncMinTR+ 100e-6 - 1.5e-3;
        if this.seq_params.doFastestPrescan
            syncDelay= 0;
        end
        for i=1: this.seq_params.nPrescans
            % insert trigger
            this.seq.addBlock(mr_trig, mr_gradFreeDelay);
            %
            this.seq.addBlock(mr_gx,mr_gy,mr_adc);
            this.seq.addBlock(mr.rotate('z',phi,mr_gxSpoil,mr_gySpoil));

            this.seq.addBlock(mr.makeDelay(1.5e-3));

            if syncDelay>0
                this.seq.addBlock(mr.makeDelay(this.Rasterize(syncDelay,'grad')));
            end
        end
    end
    %
    this.seq.addBlock(mr.makeDelay(3)) % 3 sec delay allowing skope to switch back into normal field monitoring mode.
    this.seq.addBlock(mr.makeLabel('SET','DUM', 0));
end

% Define sequence blocks
for r=1:Nreps
    %this.seq.addBlock(mr_trigPhy, mr.makeLabel('SET','SLC', 0));
    this.seq.addBlock(mr.makeLabel('SET','SLC', 0));
    for s=1:Nslices
        this.seq.addBlock(mr_rfFatSat,mr_gzFatSat);
        mr_rf.freqOffset=mr_gz.amplitude*(1+sliceGap)*thickness*(s-1-(Nslices-1)/2);
        mr_rf.phaseOffset=-2*pi*mr_rf.freqOffset*mr.calcRfCenter(mr_rf); % compensate for the slice-offset induced phase
        this.seq.addBlock(mr_rf,mr_gz);
        this.seq.addBlock(mr_gzReph);

        this.seq.addBlock(mr.makeDelay(delayTE));

        % insert trigger
        this.seq.addBlock(mr_trig, mr_gradFreeDelay);

        this.seq.addBlock(mr.rotate('z',phi,mr_gx,mr_gy,mr_adc));
        this.seq.addBlock(mr.rotate('z',phi,mr_gxSpoil,mr_gySpoil,mr_gzSpoil));

        this.seq.addBlock(mr.makeDelay(delayTR)); 
        this.seq.addBlock(mr.makeLabel('INC','SLC', 1));
    end
    this.seq.addBlock(mr.makeLabel('INC','REP', 1)); %
end

disp('-> Finished...')

end


%
% function gos= designArchimedean(this,deltak,kRadius,kSamples)
% % define k-space parameters
% %readoutTime = 4.2e-4;
% 
% % calculate a raw Archimedian spiral trajectory
% clear ka;
% ka(kRadius*kSamples+1)=1i; % init as complex
% for c=0:kRadius*kSamples
%     r=deltak*c/kSamples;
%     a=mod(c,kSamples)*2*pi/kSamples;
%     ka(c+1)=r*exp(1i*a);
% end
% ka=[real(ka); imag(ka)];
% % calculate gradients and slew rates
% [ga, sa]=mr.traj2grad(ka);
% 
% % limit analysis
% safety_magrin= 0.94; % we need that  otherwise we just about violate the slew rate due to the rounding errors
% dt_gcomp=abs(ga)/(this.sys.maxGrad*safety_magrin)*this.sys.gradRasterTime;
% dt_gabs=abs(ga(1,:)+1i*ga(2,:))/(this.sys.maxGrad*safety_magrin)*this.sys.gradRasterTime;
% dt_scomp=sqrt(abs(sa)/(this.sys.maxSlew*safety_magrin))*this.sys.gradRasterTime;
% dt_sabs=sqrt(abs(sa(1,:)+1i*sa(2,:))/(this.sys.maxSlew*safety_magrin))*this.sys.gradRasterTime;
% 
% figure;plot([dt_gabs; max(dt_gcomp); dt_sabs; max(dt_scomp)]');title('time stepping defined by gradient and slew-rate');
% 
% dt_smooth=max([dt_gabs;dt_sabs]);
% dt_rough=max([dt_gcomp;dt_scomp]);
% 
% % apply the lower limit not to lose the trajectory detail
% nppr= 4; %8;
% dt_min=nppr*this.sys.gradRasterTime/kSamples; % we want at least 4 points per revolution
% dt_smooth0=dt_smooth;
% dt_rough0=dt_rough;
% dt_smooth(dt_smooth<dt_min)=dt_min;
% dt_rough(dt_rough<dt_min)=dt_min;
% 
% figure;plot([dt_smooth0; dt_smooth; dt_rough0; dt_rough]');title('combined time stepping');
% 
% t_smooth=[0 cumsum(dt_smooth,2)];
% t_rough=[0 cumsum(dt_rough,2)];
% 
% kopt_smooth=interp1(t_smooth, ka', (0:floor(t_smooth(end)/this.sys.gradRasterTime))*this.sys.gradRasterTime)';
% kopt_rough=interp1(t_rough, ka', (0:floor(t_rough(end)/this.sys.gradRasterTime))*this.sys.gradRasterTime)';
% 
% % analyze what we've got
% fprintf('duration orig %d us\n', round(1e6*this.sys.gradRasterTime*length(ka)));
% fprintf('duration smooth %d us\n', round(1e6*this.sys.gradRasterTime*length(kopt_smooth)));
% fprintf('duration rough %d us\n', round(1e6*this.sys.gradRasterTime*length(kopt_rough)));
% 
% [gos, sos]=mr.traj2grad(kopt_smooth);
% [gor, sor]=mr.traj2grad(kopt_rough);
% 
% figure;plot([gos;abs(gos(1,:)+1i*gos(2,:))]');title('gradient with smooth (abs) constraint')
% figure;plot([gor;abs(gor(1,:)+1i*gor(2,:))]');title('gradient with rough (component) constraint')
% 
% figure;plot([sos;abs(sos(1,:)+1i*sos(2,:))]');title('slew rate with smooth (abs) constraint')
% figure;plot([sor;abs(sor(1,:)+1i*sor(2,:))]');title('slew rate with rough (component) constraint')
% figure;plot(kopt_smooth(1,:),kopt_smooth(2,:));title('k trajectory')
% end
