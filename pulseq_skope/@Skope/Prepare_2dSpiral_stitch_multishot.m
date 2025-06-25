function Prepare_2dSpiral_stitch_multishot(this)
% Prepare 2D spiral sequence
% this is based on writeSpiral.m from pulseq.
%

this.seq_params.isSpiral = true;
this.seq_params.resolution= this.seq_params.fov/ this.seq_params.N;

isHighRes= false; 
if this.seq_params.resolution<= 0.5e-3  
    isHighRes= true;
end

% multi-shot spiral
nShots         = this.seq_params.nSpiralInterleaves;

% setting for sequence parameters
fov           = this.seq_params.fov / this.seq_params.accelerationFactor;
N             = this.seq_params.N / this.seq_params.accelerationFactor;
alpha         = this.seq_params.alpha;
alpha_fatsat  = this.seq_params.alpha_fatsat;
thickness     = this.seq_params.thickness;
Nslices       = this.seq_params.Nslices;
TE            = this.seq_params.TE;
gradFreeDelay = this.seq_params.gradFreeDelay;
sliceGap      = this.seq_params.sliceGap;
nsegs2measure = this.seq_params.nSegments2measure;
safetyMargin  = this.seq_params.gradSafetyMargin;
Nreps         = this.seq_params.nRepeats;
adcSamplesPerSegment = this.seq_params.maxAdcSegmentLength;

% setting for dephsing signal model
probeType     = this.seq_params.probeType;
probeRadius   = this.seq_params.probeRadius;
signalCutoff  = this.seq_params.signalCutoff;



phi= 0;%pi/2; % orientation of the readout e.g. for interleaving

% Create fat-sat pulse
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6* this.sys.B0* this.sys.gamma;
mr_rfFatSat = mr.makeGaussPulse(alpha_fatsat*pi/180,'system',this.sys,'Duration',8e-3,'dwell',10e-6,...
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

if ~this.seq_params.isTimeOptimal % design Archimedean spiral
    Oversampling= 2; % by looking at the periphery of the spiral I would say it needs to be at least 2
    kRadius = round(N/2);
    kSamples=round(2*pi*kRadius)*Oversampling;
    adcSamplesDesired=kRadius*kSamples;
    deltak=1/fov;
    spiral_grad_shape= this.designArchimedean(deltak, kRadius, kSamples, safetyMargin);
else% design time optimal spiral
    %[k,g,s,time] = designVariableDensitySpiral(this, Nitlv, isRotationallyInvariant, res, fov, radius)
    Oversampling= 10;
    %nShots= this.seq_params.nSpiralInterleaves; 
    res= this.seq_params.resolution;
    FOV= [fov fov]; radius= [0 1];
    isRotationallyVariant= false; %this.seq_params.isRotationallyVariant;
    [~,vdspiral]= this.designVariableDensitySpiral(nShots,isRotationallyVariant,res,FOV,radius,safetyMargin);
    spiral_grad_shape= vdspiral.';
    adcSamplesDesired= Oversampling.* size(spiral_grad_shape,2);
end


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
[~, nGrad] = size(spiral_grad_shape);

t_grad0    = linspace(0.5, nGrad-0.5, nGrad)*this.seq.gradRasterTime;
t_grad     = (0:nGrad) * this.seq.gradRasterTime;
t_adc      = linspace(0.5, mr_adc.numSamples-0.5, mr_adc.numSamples)*mr_adc.dwell;

grad_shape = interp1(t_grad0, spiral_grad_shape', t_grad, 'makima', 0);
k = (grad_shape(2:end, :) + grad_shape(1:end-1, :)) .* (this.sys.gradRasterTime / 2);
k = cumsum(k, 1);
k = [zeros(1,3); k];
k_adc = interp1(t_grad, k, t_adc, 'makima', 0);

spiral_grad_shape_all = zeros(nShots, size(spiral_grad_shape, 1), size(spiral_grad_shape, 2));
for ishot = 1:nShots
    theta = (ishot-1)/nShots*2*pi;
    R = [cos(theta), -sin(theta), 0;
         sin(theta),  cos(theta), 0;
                  0,           0, 0];

    spiral_grad_shape_all(ishot, :, :) = R * spiral_grad_shape;
end

grad0.shape     = spiral_grad_shape;
grad0.shape_all = spiral_grad_shape_all;
grad0.dt        = this.sys.gradRasterTime;
grad0.unit      = 'Hz/m';
grad0.gamma     = this.sys.gamma;
this.seq_params.grad0 = grad0;
% fn= ['xw_sp2d-',num2str(1e3*this.seq_params.resolution,2),'mm-r',num2str(this.seq_params.accelerationFactor)];
% save([fn,'.mat'],'grad0')

% readout grad
mr_gx = mr.makeArbitraryGrad('x', (-1)* spiral_grad_shape(1,:), 'system', this.sys, 'first', 0, 'last', 0);
mr_gy = mr.makeArbitraryGrad('y',       spiral_grad_shape(2,:), 'system', this.sys, 'first', 0, 'last', 0);

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
disp('=====================')
if ~nsegs2measure        % determine number of segments automatically
    if isHighRes % mostly dephasing dominated.
        triggerDelays= calcTriggerDelays(grad0, probeType, probeRadius, signalCutoff);
        segDuration= max(diff(triggerDelays));
        nsegs2measure= length(triggerDelays);
        disp('-> Auto mode: High resolution scenario: ')
    else % low resolution, mostly T2* decay dominated
        nsegs2measure= ceil(mr.calcDuration(mr_gx,mr_gy,mr_adc)./ (0.9*this.sys.effProbeLifetime) );
        segDuration= this.Rasterize(mr.calcDuration(mr_gx,mr_gy,mr_adc)./ nsegs2measure,"grad");
        triggerDelays= (0: nsegs2measure-1).* segDuration;
        disp('-> Auto mode: Long readout scenario: ')
    end

else % manually specify segment number
    segDuration= this.Rasterize(mr.calcDuration(mr_gx,mr_gy,mr_adc)./ nsegs2measure,"grad");
    triggerDelays= (0: nsegs2measure-1).* segDuration;
    disp('-> Manual mode: ')
end
fprintf('-> Number of gradient segments: %d \n', nsegs2measure);
this.seq_params.nSegments2measure= nsegs2measure;
this.seq_params.triggerDelays= triggerDelays;
this.seq_params.readoutGradientDuration= mr.calcDuration(mr_gx,mr_gy,mr_adc); 

trig2AdcTime= gradFreeDelay;
skopeAcqDuration= trig2AdcTime + segDuration + 1e-3;
this.seq_params.trigger2AdcTime= trig2AdcTime;
this.seq_params.skopeAcqDuration= skopeAcqDuration;

this.seq_params.nDynamics= ceil(Nreps * Nslices * nShots * nsegs2measure / (this.seq_params.nInterleaves));

this.seq_params.segmentDuration= segDuration;
fprintf('-> Max gradient segment duration (ms): %2.2f \n', 1e3.*segDuration);
disp('=====================')

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

this.seq_params.end_block = [];
switch this.seq_params.stitchMode

    case 'concurrent'
        % define a segment indices array indicating which gradient segment to measure at a
        % given TR
        segIndexArray= zeros(1,Nreps* Nslices);
        if nsegs2measure>1

            if Nreps*Nslices/nsegs2measure< 1
                Nreps= ceil(nsegs2measure./Nslices);    this.seq_params.nRepeats= Nreps;
                disp('-> NOTE: number of repetitions is increased to map out the entire readout...')
            end

            segIndexArray= repmat(0:(nsegs2measure-1), 1, ceil(Nreps*Nslices/nsegs2measure));
        end

        % this.seq_params.nDynamics= ceil(Nreps* Nslices / (this.seq_params.nInterleaves));

        % Define sequence blocks
        
        for r=1:Nreps
            %this.seq.addBlock(mr_trigPhy, mr.makeLabel('SET','SLC', 0));

            this.seq.addBlock(mr.makeLabel('SET', 'LIN', 0));
            for ishot = 1:nShots
                mr_gx.waveform = -spiral_grad_shape_all(ishot, 1, :);
                mr_gy.waveform =  spiral_grad_shape_all(ishot, 2, :);

                this.seq.addBlock(mr.makeLabel('SET','SLC', 0));
                
                counter = 1;
                for s=1:Nslices
                   
                    this.seq.addBlock(mr_rfFatSat,mr_gzFatSat);
                    mr_rf.freqOffset=mr_gz.amplitude*(1+sliceGap)*thickness*(s-1-(Nslices-1)/2);
                    mr_rf.phaseOffset=-2*pi*mr_rf.freqOffset*mr.calcRfCenter(mr_rf); % compensate for the slice-offset induced phase
                    this.seq.addBlock(mr_rf,mr_gz);
                    this.seq.addBlock(mr_gzReph);

                    this.seq.addBlock(mr.makeDelay(delayTE));

                    if segIndexArray(counter)== 0
                        this.seq.addBlock(mr_trig, mr_gradFreeDelay)
                        this.seq.addBlock(mr.rotate('z',phi,mr_gx,mr_gy,mr_adc));
                    else
                        this.seq.addBlock(mr_gradFreeDelay);

                        mr_trig.delay= triggerDelays(segIndexArray(counter)+ 1)- mr_gradFreeDelay.delay;
                        % disp(triggerDelays(segIndexArray(counter)+ 1))
                        this.seq.addBlock(mr.rotate('z',phi,mr_gx,mr_gy,mr_adc, mr_trig));
                    end

                    this.seq.addBlock(mr.rotate('z',phi,mr_gxSpoil,mr_gySpoil,mr_gzSpoil));

                    this.seq.addBlock(mr.makeDelay(delayTR));

                    this.seq.addBlock(mr.makeLabel('INC','SLC', 1));

                    counter= counter+ 1;
                end

                this.seq.addBlock(mr.makeLabel('INC', 'LIN', 1));

                if ishot == nShots
                    this.seq_params.end_block = [this.seq_params.end_block, size(this.seq.blockDurations, 2)];
                end
     
            end

            this.seq.addBlock(mr.makeLabel('INC','REP', 1)); %

        end

    case 'sequential'

    % Define sequence blocks
    for r=1:Nreps
        %this.seq.addBlock(mr_trigPhy, mr.makeLabel('SET','SLC', 0));
        this.seq.addBlock(mr.makeLabel('SET','SLC', 0));

        for s=1:Nslices
            this.seq.addBlock(mr.makeLabel('SET', 'LIN', 0));
            for ishot = 1:nShots
                mr_gx.waveform = -spiral_grad_shape_all(ishot, 1, :);
                mr_gy.waveform =  spiral_grad_shape_all(ishot, 2, :);
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
    
                counter= 1;
                while counter< nsegs2measure
                    this.seq.addBlock(mr_rfFatSat,mr_gzFatSat);
                    %         mr_rf.freqOffset=mr_gz.amplitude*(1+sliceGap)*thickness*(s-1-(Nslices-1)/2);
                    %         mr_rf.phaseOffset=-2*pi*mr_rf.freqOffset*mr.calcRfCenter(mr_rf); % compensate for the slice-offset induced phase
                    this.seq.addBlock(mr_rf,mr_gz);
                    this.seq.addBlock(mr_gzReph);
    
                    this.seq.addBlock(mr.makeDelay(delayTE));
    
                    % insert trigger
                    this.seq.addBlock(mr_gradFreeDelay);
    
                    mr_trig.delay= segDuration.* counter- mr_gradFreeDelay.delay;
                    this.seq.addBlock(mr.rotate('z',phi,mr_gx,mr_gy,mr_adc, mr_trig));
                    this.seq.addBlock(mr.rotate('z',phi,mr_gxSpoil,mr_gySpoil,mr_gzSpoil));
    
                    this.seq.addBlock(mr.makeDelay(delayTR));
    
                    counter= counter+ 1;
                end
                this.seq.addBlock(mr.makeLabel('INC', 'LIN', 1));
            end
            this.seq.addBlock(mr.makeLabel('INC','SLC', 1));
        end

        this.seq.addBlock(mr.makeLabel('INC','REP', 1)); %
    end

    case 'interleaved'
    % first segment

    % sequence blocks

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

    % other segments
    counter= 1;
    while counter< nsegs2measure
        this.seq.addBlock(mr.makeDelay(this.seq_params.interSessionDelay));

        % sequence blocks

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
                this.seq.addBlock(mr_gradFreeDelay);

                mr_trig.delay= segDuration.* counter- mr_gradFreeDelay.delay;
                this.seq.addBlock(mr.rotate('z',phi,mr_gx,mr_gy,mr_adc, mr_trig));
                this.seq.addBlock(mr.rotate('z',phi,mr_gxSpoil,mr_gySpoil,mr_gzSpoil));

                this.seq.addBlock(mr.makeDelay(delayTR));

                this.seq.addBlock(mr.makeLabel('INC','SLC', 1));
            end

            this.seq.addBlock(mr.makeLabel('INC','REP', 1)); %
        end

        counter= counter+ 1;

    end

    otherwise
        disp('-> Acquisition scheme specified is not supported... returned...')

end

disp('-> Finished...')

end

%%%%=====================
function triggerDelays= calcTriggerDelays(grad0, probeType, probeRadius, signalCutoff)
%%% smartly segment the gradient for ultrahigh resolution
% r= 0.4e-3; % radius of the field probe in m. Our motion probe is 0.8 mm in radius per Cameron's email..
% signalCutoff= 0.35; % signal cutoff level.

% signal decay due to dephasing inside a voxel follows sinc(a) where a is the max phase angle relative to the middle of the voxel.
x=0:0.005:1;
maxPhase= x(find(sinc(x)>signalCutoff, 1, 'last' ))* pi; % in rads.
maxK= maxPhase./ probeRadius;

% grad= hz2tesla(grad0.shape(1:2,:));
grad = grad0.shape(1:2,:) / grad0.gamma;
dt= grad0.dt;
ktraj= calc_ktraj_from_grad(grad, probeType, false, dt);
ktrajAbs= sqrt(ktraj(1,:).^2+ ktraj(2,:).^2);

idxCutoff= 1;
idx2measure= [];
while ~isempty(idxCutoff)

    idxCutoff= find(ktrajAbs> maxK,1,'first');

    if ~isempty(idxCutoff)
        idx2measure= [idx2measure idxCutoff];
        grad= grad(:, idxCutoff:end);
        ktraj= calc_ktraj_from_grad(grad, probeType, false, dt);
        ktrajAbs= sqrt(ktraj(1,:).^2+ ktraj(2,:).^2);
    end
end

triggerDelays= [0 cumsum(idx2measure-1).* dt]; 

end


