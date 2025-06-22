function Prepare_2dEpi_stitch(this)
% Prepare 2D EPI Stitch sequence
% this is based on writeEpi_label.m and writeEpiRS.m from pulseq.


% setting for dephsing signal model
probeType     = this.seq_params.probeType;
probeRadius   = this.seq_params.probeRadius;
signalCutoff  = this.seq_params.signalCutoff;

nsegs2measure = this.seq_params.nSegments2measure;
isHighRes     = this.seq_params.isHighRes;

this.seq_params.resolution = this.seq_params.fov / this.seq_params.N;
fov_y                      = this.seq_params.fov / this.seq_params.accelerationFactor;
fov_x                      = this.seq_params.fov;
Nx                         = this.seq_params.N; 
Ny                         = Nx / this.seq_params.accelerationFactor;
alpha                      = this.seq_params.alpha;
alpha_fatsat               = this.seq_params.alpha_fatsat;
thickness                  = this.seq_params.thickness;
Nslices                    = this.seq_params.Nslices;
TE                         = this.seq_params.TE;
gradFreeDelay              = this.seq_params.gradFreeDelay;
sliceGap                   = this.seq_params.sliceGap;

Nreps                      = this.seq_params.nRepeats;
Navigator                  = this.seq_params.nNavigators;

pe_enable                  = 1;  % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
ro_os                      = 2; %1; % oversampling factor (in contrast to the product sequence we don't really need it)
readoutTime                = this.seq_params.readoutTime;     % this controls the readout bandwidth
partFourierFactor          = this.seq_params.partialFourier;  % partial Fourier factor: 1: full sampling 0: start with ky=0

% Create fat-sat pulse
sat_ppm     = -3.45;
sat_freq    = sat_ppm*1e-6* this.sys.B0* this.sys.gamma;
mr_rfFatSat = mr.makeGaussPulse(alpha_fatsat*pi/180,'system',this.sys,'Duration',8e-3,'dwell',10e-6,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
mr_rfFatSat.phaseOffset=-2*pi*mr_rfFatSat.freqOffset*mr.calcRfCenter(mr_rfFatSat); % compensate for the frequency-offset induced phase
mr_gzFatSat = mr.makeTrapezoid('z',this.sys,'delay',mr.calcDuration(mr_rfFatSat),'Area',0.1/1e-4, 'duration', 1.5e-3); % spoil up to 0.1mm

% Create slice selection pulse and gradient
[mr_rf, mr_gz, mr_gzReph] = mr.makeSincPulse(alpha*pi/180,'system',this.sys,'Duration',3e-3,...
    'SliceThickness',thickness,'apodization',0.42,'timeBwProduct',4,'use','excitation');

% define the trigger to play out
mr_trig = mr.makeDigitalOutputPulse(this.seq_params.trigChannel,'duration', this.sys.gradRasterTime);
mr_gradFreeDelay = mr.makeDelay(gradFreeDelay);

% Define other gradients and mr_adc events
deltak_y = 1/fov_y;
deltak_x = 1/fov_x;
kWidth   = Nx*deltak_x;

% Phase blip in shortest possible time
% we round-up the duration to 2x the gradient raster time
blip_dur = ceil(2*sqrt(deltak_y/this.sys.maxSlew)/(2*this.sys.gradRasterTime))*(2*this.sys.gradRasterTime); 
% the split code below fails if this really makes a trpezoid instead of a triangle...
% we use negative blips to save one k-space line on our way towards the k-space center
mr_gy    = mr.makeTrapezoid('y', this.sys, 'Area', -deltak_y, 'Duration', blip_dur);

% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end each equal to a half of blip_dur
% the area between the blips should be defined by kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slewrate and then scale down the amlitude to fix the area
extra_area      = blip_dur / 2 * blip_dur / 2 * this.sys.maxSlew; % check unit!;
mr_gx           = mr.makeTrapezoid('x', this.sys, 'Area', kWidth+extra_area, 'duration', readoutTime+blip_dur);
actual_area     = mr_gx.area - mr_gx.amplitude/mr_gx.riseTime*blip_dur/2*blip_dur/2/2 ...
    - mr_gx.amplitude/mr_gx.fallTime*blip_dur/2*blip_dur/2/2;
mr_gx.amplitude = mr_gx.amplitude / actual_area * kWidth;
mr_gx.area      = mr_gx.amplitude * (mr_gx.flatTime + mr_gx.riseTime/2 + mr_gx.fallTime/2);
mr_gx.flatArea  = mr_gx.amplitude * mr_gx.flatTime;

% calculate mr_adc
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be qite different from Nx and
% readoutTime/Nx, respectively.
adcDwellNyquist = deltak_x / mr_gx.amplitude / ro_os;

adcDwell        = floor(adcDwellNyquist/this.sys.adcRasterTime)*this.sys.adcRasterTime; % round-down dwell time to 100 ns
adcSamples      = floor(readoutTime/adcDwell/4)*4; % on Siemens the number of mr_adc samples need to be divisible by 4
% MZ: no idea, whether ceil, round or floor is better for the adcSamples...
mr_adc          = mr.makeAdc(adcSamples, 'Dwell', adcDwell, 'Delay',blip_dur/2);
% realign the mr_adc with respect to the gradient
time_to_center  = mr_adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
mr_adc.delay    = round((mr_gx.riseTime + mr_gx.flatTime/2 - time_to_center)/1e-6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
% this rounding actually makes the sampling points on odd and even readouts
% to appear misalligned. However, on the real hardware this misalignment is
% much stronger anyways due to the grdient delays

% FOV positioning requires alignment to grad. raster... -> TODO

% split the blip into two halves and produce a combined synthetic gradient
gy_parts                       = mr.splitGradientAt(mr_gy, blip_dur/2, this.sys);
[mr_gyBlipup, mr_gyBlipdown,~] = mr.align('right', gy_parts(1), 'left', gy_parts(2), mr_gx);
mr_gyBlipdownup                = mr.addGradients({mr_gyBlipdown, mr_gyBlipup}, this.sys);

% pe_enable support
mr_gyBlipup.waveform     = mr_gyBlipup.waveform     * pe_enable;
mr_gyBlipdown.waveform   = mr_gyBlipdown.waveform   * pe_enable;
mr_gyBlipdownup.waveform = mr_gyBlipdownup.waveform * pe_enable;

% phase encoding and partial Fourier
Ny_pre  = round(partFourierFactor*Ny/2 - 1); % PE steps prior to ky=0, excluding the central line
Ny_post = round(Ny/2 + 1); % PE lines after the k-space center including the central line
Ny_meas = Ny_pre + Ny_post;

% negate x gradient to be consistent with IDEA
if mr_gx.channel == 'x'
    mr_gx.amplitude = (-1) * mr_gx.amplitude;
    mr_gx.area      = (-1) * mr_gx.area;
    mr_gx.flatArea  = (-1) * mr_gx.flatArea;
end

% Pre-phasing gradients
mr_gxPre = mr.makeTrapezoid('x', this.sys, 'Area', -mr_gx.area/2  );
mr_gyPre = mr.makeTrapezoid('y', this.sys, 'Area', Ny_pre*deltak_y);
% prolong the duration slightly to decrease PNS
mr_gxPre = mr.makeTrapezoid('x', this.sys, 'Area', -mr_gx.area/2  , 'duration', mr.calcDuration(mr_gxPre)*1.2);
mr_gyPre = mr.makeTrapezoid('y', this.sys, 'Area', Ny_pre*deltak_y, 'duration', mr.calcDuration(mr_gyPre)*1.2);

% if Navigator==0
%     [mr_gxPre,mr_gyPre]=mr.align('right',mr_gxPre,mr_gyPre);
% end

% relax the PE prepahser to reduce stimulation
mr_gyPre           = mr.makeTrapezoid('y', this.sys, 'Area', mr_gyPre.area, 'Duration', mr.calcDuration(mr_gxPre,mr_gyPre));
mr_gyPre.amplitude = mr_gyPre.amplitude*pe_enable;

mr_gzSpoil = mr.makeTrapezoid('z', this.sys, 'Area', deltak_x*Nx*4);

% Calculate timing
if Navigator>0
    minTE= mr_gz.fallTime + mr_gz.flatTime/2 + mr.calcDuration(mr_gzReph) ...
        + gradFreeDelay + mr.calcDuration(mr_gxPre)+ Navigator*mr.calcDuration(mr_gx) ...
        + mr.calcDuration(mr_gyPre) + (Ny_pre+0.5)*mr.calcDuration(mr_gx);
else % no Navi
    minTE= mr_gz.fallTime + mr_gz.flatTime/2 + mr.calcDuration(mr_gzReph) ...
        + gradFreeDelay + mr.calcDuration(mr_gxPre,mr_gyPre)+ (Ny_pre+0.5)*mr.calcDuration(mr_gx);
end
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
%
%minTR1= this.seq_params.skopeMinTR+ 100e-6;
minExcitTR = mr.calcDuration(mr_rfFatSat)+ mr.calcDuration(mr_gzFatSat) ...
    + mr.calcDuration(mr_gz)/2+ TE+ (Ny_post-0.5)*mr.calcDuration(mr_gx) ...
    + mr.calcDuration(mr_gzSpoil);
delayTR    = this.calculateTiming(minExcitTR);

%
if Navigator>0
    trig2AdcTime = gradFreeDelay + mr.calcDuration(mr_gxPre) + Navigator*mr.calcDuration(mr_gx)...
        + mr.calcDuration(mr_gyPre) + mr_adc.delay;
else % no Navi
    trig2AdcTime = gradFreeDelay + mr.calcDuration(mr_gxPre,mr_gyPre) + mr_adc.delay;
end
skopeAcqDuration                 = (trig2AdcTime- mr_adc.delay) + Ny_meas*mr.calcDuration(mr_gx) + 1e-3;
this.seq_params.trigger2AdcTime  = trig2AdcTime;
this.seq_params.skopeAcqDuration = skopeAcqDuration;
this.seq_params.nDynamics        = ceil(Nreps* Nslices / (this.seq_params.nInterleaves));

if this.seq_params.skopeAcqDuration> this.sys.effProbeLifetime
    fprintf('-> NOTE: skope acq duration is longer than effective probe lifetime of %d ms...\n', 1e3* this.sys.effProbeLifetime);
    disp('-> Please consider adjusting imaging parameters (eg, decrease resolution or increase undersampling rate')
    disp('-> to avoid noisy field measurement toward the tail of the kspace trajectory...')
end

if Navigator>0
    trig2AdcTimePhaseCorrection= gradFreeDelay+ mr.calcDuration(mr_gxPre)+ mr_adc.delay;
    this.seq_params.trigger2AdcTimePhaseCorrection= trig2AdcTimePhaseCorrection;
end

echoSpacing                     = mr.calcDuration(mr_gx);
this.seq_params.echoSpacing     = echoSpacing;
this.seq_params.echoTrainLength = Ny_meas;


%% go through a whole readout to get the gradient shape
seq = mr.Sequence(this.sys);

if Navigator>0
    seq.addBlock(mr_gxPre);
    for n=1:Navigator
        seq.addBlock(mr_gx, mr_adc);
        mr_gx.amplitude = -mr_gx.amplitude;   % Reverse polarity of read gradient
    end
    seq.addBlock(mr_gyPre);% lin/nav/avg reset
else % no Navigator
    seq.addBlock(mr_gxPre, mr_gyPre);% lin/nav/avg reset
end

for i=1:Ny_meas
    if i==1
        seq.addBlock(mr_gx, mr_gyBlipup    , mr_adc); % Read the first line of k-space with a single half-blip at the end
    elseif i==Ny_meas
        seq.addBlock(mr_gx, mr_gyBlipdown  , mr_adc); % Read the last line of k-space with a single half-blip at the beginning
    else
        seq.addBlock(mr_gx, mr_gyBlipdownup, mr_adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
    end
    mr_gx.amplitude = -mr_gx.amplitude;   % Reverse polarity of read gradient
end
[gw,~,~,t_adc] = seq.waveforms_and_times();
% figure; title('gradient wave form, in T/m');
% plot(gw{1}(1,:),gw{1}(2,:),gw{2}(1,:),gw{2}(2,:),gw{3}(1,:),gw{3}(2,:)); % plot the entire gradient shape
% find beginning and end times and resample GWs to a regular sampling raster
tf = [];
tl = [];
for i = 1:3
    if size(gw{i}, 2) >0
        tf(end+1) = gw{i}(1,   1);
        tl(end+1) = gw{i}(1, end);
    end
end
nt_min = floor(min(tf) / seq.gradRasterTime + eps); 
nt_max = ceil( max(tl) / seq.gradRasterTime - eps);
% shift raster positions to the centers of the raster periods
nt_min = nt_min + 0.5;
nt_max = nt_max - 0.5;
if nt_min < 0.5
    nt_min = 0.5;
end
t_axis     = (nt_min:nt_max) * seq.gradRasterTime;
grad_shape = zeros(length(t_axis), 3);
for i = 1:3
    if size(gw{i},2) > 0
        grad_shape(:,i) = interp1(gw{i}(1,:), gw{i}(2,:), t_axis, 'linear', 0);
    end
end

grad0.shape  = grad_shape';
grad0.dt     = this.sys.gradRasterTime;
grad0.unit   = 'Hz/m';
grad0.gamma  = this.sys.gamma;
grad0.t_grad = t_axis;
grad0.t_adc  = t_adc;

this.seq_params.grad0 = grad0;
% fn= ['xw_sp2d-',num2str(1e3*this.seq_params.resolution,2),'mm-r',num2str(this.seq_params.accelerationFactor)];
% save([fn,'.mat'],'grad0')


%% 
disp('=====================')

grad_duration = length(grad0.t_grad) * grad0.dt;
if ~nsegs2measure        % determine number of segments automatically
    if isHighRes % mostly dephasing dominated.
        triggerDelays = calcTriggerDelays(grad0, probeType, probeRadius, signalCutoff);
        segDuration   = max(diff(triggerDelays));
        nsegs2measure = length(triggerDelays);
        disp('-> Auto mode: High resolution scenario: ')
    else % low resolution, mostly T2* decay dominated
        nsegs2measure = ceil(grad_duration./ (0.9*this.sys.effProbeLifetime) );
        segDuration   = this.Rasterize(grad_duration./ nsegs2measure,"grad");
        triggerDelays = (0: nsegs2measure-1).* segDuration;
        disp('-> Auto mode: Long readout scenario: ')
    end
else % manually specify segment number
    segDuration   = this.Rasterize(grad_duration./ nsegs2measure,"grad");
    triggerDelays = (0: nsegs2measure-1).* segDuration;
    disp('-> Manual mode: ')
end
fprintf('-> Number of gradient segments: %d \n', nsegs2measure);
this.seq_params.nSegments2measure       = nsegs2measure;
this.seq_params.triggerDelays           = triggerDelays;
this.seq_params.readoutGradientDuration = grad_duration;  % the whole EPI readout

trig2AdcTime     = gradFreeDelay;
skopeAcqDuration = trig2AdcTime + segDuration + 1e-3;
this.seq_params.trigger2AdcTime  = trig2AdcTime;
this.seq_params.skopeAcqDuration = skopeAcqDuration;

this.seq_params.nDynamics= ceil(Nreps* Nslices* nsegs2measure / (this.seq_params.nInterleaves));

this.seq_params.segmentDuration= segDuration;
fprintf('-> Max gradient segment duration (ms): %2.2f \n', 1e3.*segDuration);
disp('=====================')

if this.seq_params.skopeAcqDuration> this.sys.effProbeLifetime
    fprintf('-> NOTE: skope acq duration is longer than effective probe lifetime of %d ms...\n', 1e3* this.sys.effProbeLifetime);
    disp('-> Please consider adjusting imaging parameters (eg, decrease resolution or increase undersampling rate')
    disp('-> to avoid noisy field measurement toward the tail of the kspace trajectory...')
end
   
%% Prescans for sync
if this.seq_params.nPrescans>0

    disp('=> Adding sync scans...')
    this.seq.addBlock(mr.makeLabel('SET','DUM', 1)); % set dummy scans

    if this.seq_params.useSingleAdcSegment4Sync
        disp('-> Using a single ADC segment for prescans...')

        syncMinTR = trig2AdcTime + mr.calcDuration(mr_gx);
        syncDelay = this.seq_params.skopeInterleaveTR - syncMinTR + 100e-6 - 1.5e-3;
        if this.seq_params.doFastestPrescan
            syncDelay= 0;
        end

        for idx=1: this.seq_params.nPrescans
            % insert trigger
            this.seq.addBlock(mr_trig, mr_gradFreeDelay);
            %
            if Navigator>0
                this.seq.addBlock(mr_gxPre);
                for n=1:Navigator
                    this.seq.addBlock(mr_gx);
                    mr_gx.amplitude = -mr_gx.amplitude;   % Reverse polarity of read gradient
                end
                this.seq.addBlock(mr_gyPre);%
            else % no Navigator
                this.seq.addBlock(mr_gxPre, mr_gyPre);%
            end

            this.seq.addBlock(mr_gx, mr_adc);

            this.seq.addBlock(mr.makeDelay(1.5e-3));

            if syncDelay>0
                this.seq.addBlock(mr.makeDelay(this.Rasterize(syncDelay,'grad')));
            end
        end

    else % use native ADC
        disp('-> Using native ADC for prescans...')
        syncMinTR= trig2AdcTime + Ny_meas*mr.calcDuration(mr_gx);
        syncDelay= this.seq_params.skopeInterleaveTR- syncMinTR+ 100e-6 - 1.5e-3;
        if this.seq_params.doFastestPrescan
            syncDelay= 0;
        end

        for idx=1: this.seq_params.nPrescans
            % insert trigger
            this.seq.addBlock(mr_trig, mr_gradFreeDelay);
            %
            if Navigator>0
                this.seq.addBlock(mr_gxPre);
                for n=1:Navigator
                    this.seq.addBlock(mr_gx);
                    mr_gx.amplitude = -mr_gx.amplitude;   % Reverse polarity of read gradient
                end
                this.seq.addBlock(mr_gyPre);%
            else % no Navigator
                this.seq.addBlock(mr_gxPre, mr_gyPre);%
            end

            for i=1:Ny_meas

                if i==1
                    this.seq.addBlock(mr_gx, mr_gyBlipup,     mr_adc); % Read the first line of k-space with a single half-blip at the end
                elseif i==Ny_meas
                    this.seq.addBlock(mr_gx, mr_gyBlipdown,   mr_adc); % Read the last line of k-space with a single half-blip at the beginning
                else
                    this.seq.addBlock(mr_gx, mr_gyBlipdownup, mr_adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                end
                mr_gx.amplitude = -mr_gx.amplitude;   % Reverse polarity of read gradient
            end

            if rem(Navigator+Ny_meas,2)~=0
                mr_gx.amplitude = -mr_gx.amplitude;
            end

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

%% Define sequence blocks
switch this.seq_params.stitchMode
    case 'concurrent'
        % define a segment indices array indicating which gradient segment to measure at a given TR
        segIndexArray= zeros(1,Nreps* Nslices);
        if nsegs2measure>1
            if Nreps*Nslices/nsegs2measure< 1
                Nreps= ceil(nsegs2measure./Nslices);    this.seq_params.nRepeats= Nreps;
                disp('-> NOTE: number of repetitions is increased to map out the entire readout...')
            end
            segIndexArray= repmat(0:(nsegs2measure-1), 1, ceil(Nreps*Nslices/nsegs2measure));
        end

        this.seq_params.nDynamics= ceil(Nreps* Nslices / (this.seq_params.nInterleaves));

        counter = 1;
        for r=1:Nreps
            this.seq.addBlock(mr.makeLabel('SET','SLC', 0));
            for s=1:Nslices
                this.seq.addBlock(mr_rfFatSat,mr_gzFatSat);
                mr_rf.freqOffset=mr_gz.amplitude*(1+sliceGap)*thickness*(s-1-(Nslices-1)/2);
                mr_rf.phaseOffset=-2*pi*mr_rf.freqOffset*mr.calcRfCenter(mr_rf); % compensate for the slice-offset induced phase
                this.seq.addBlock(mr_rf,mr_gz);
                this.seq.addBlock(mr_gzReph, ...
                    mr.makeLabel('SET','NAV',1),...
                    mr.makeLabel('SET','LIN', round(Ny/2)));
        
                this.seq.addBlock(mr.makeDelay(delayTE));


                TriggerDelay = triggerDelays(segIndexArray(counter) + 1);
                % Reset duration of current TR
                TimeDelay = 0 ; % [s], the trigger start time is [mr_gradFreeDelay.delay] seconds ahead of the start of the readout gradient

                if TriggerDelay - TimeDelay <= mr_gradFreeDelay.delay - mr_trig.duration
                    mr_trig.delay = TriggerDelay - TimeDelay;
                    this.seq.addBlock(mr_trig, mr_gradFreeDelay);
                    % fprintf('TR %d, Trigger delay: %.3f ms\n', r, TriggerDelay*1e3);
                else
                    this.seq.addBlock(mr_gradFreeDelay);
                end

                TimeDelay = TimeDelay + this.seq.blockDurations(end);

        
                if Navigator>0
                    if TriggerDelay > TimeDelay && TriggerDelay - TimeDelay <= mr.calcDuration(mr_gxPre) - mr_trig.duration 
                        mr_trig.delay = TriggerDelay - TimeDelay;
                        this.seq.addBlock(mr_trig, mr_gxPre);
                        % fprintf('TR %d, Trigger delay: %.3f ms\n', r, TriggerDelay*1e3);
                    else
                        this.seq.addBlock(mr_gxPre);
                    end
                    TimeDelay = TimeDelay + this.seq.blockDurations(end);

                    for n=1:Navigator
                        if TriggerDelay > TimeDelay && TriggerDelay - TimeDelay <= mr.calcDuration(mr_gx,mr_adc) - mr_trig.duration 
                            mr_trig.delay = TriggerDelay - TimeDelay;
                            this.seq.addBlock(mr_trig, mr_gx,mr_adc, ...
                                mr.makeLabel('SET','REV', mr_gx.amplitude<0), mr.makeLabel('SET','SEG', mr_gx.amplitude<0), mr.makeLabel('SET','AVG',n==Navigator));
                            % fprintf('TR %d, Trigger delay: %.3f ms\n', r, TriggerDelay*1e3);
                        else
                            this.seq.addBlock(mr_gx,mr_adc, ...
                                mr.makeLabel('SET','REV', mr_gx.amplitude<0), mr.makeLabel('SET','SEG', mr_gx.amplitude<0), mr.makeLabel('SET','AVG',n==Navigator));
                        end

                        mr_gx.amplitude = -mr_gx.amplitude;   % Reverse polarity of read gradient
                        TimeDelay = TimeDelay + this.seq.blockDurations(end);
                    end
                    if TriggerDelay > TimeDelay && TriggerDelay - TimeDelay <= mr.calcDuration(mr_gyPre) - mr_trig.duration 
                        mr_trig.delay = TriggerDelay - TimeDelay;
                        this.seq.addBlock(mr_trig, mr_gyPre, ...
                            mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','NAV', 0), mr.makeLabel('SET','AVG', 0) );% lin/nav/avg reset
                        % fprintf('TR %d, Trigger delay: %.3f ms\n', r, TriggerDelay*1e3);
                    else
                        this.seq.addBlock(mr_gyPre, ...
                            mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','NAV', 0), mr.makeLabel('SET','AVG', 0) );% lin/nav/avg reset
                    end
                    TimeDelay = TimeDelay + this.seq.blockDurations(end);
                else % no Navigator
                    if TriggerDelay > TimeDelay && TriggerDelay - TimeDelay <= mr.calcDuration(mr_gxPre, mr_gyPre) - mr_trig.duration 
                        mr_trig.delay = TriggerDelay - TimeDelay;
                        this.seq.addBlock(mr_trig, mr_gxPre, mr_gyPre, ...
                            mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','NAV', 0), mr.makeLabel('SET','AVG', 0) );% lin/nav/avg reset
                        % fprintf('TR %d, Trigger delay: %.3f ms\n', r, TriggerDelay*1e3);
                    else
                        this.seq.addBlock(mr_gxPre, mr_gyPre, ...
                            mr.makeLabel('SET','LIN', 0), mr.makeLabel('SET','NAV', 0), mr.makeLabel('SET','AVG', 0) );% lin/nav/avg reset
                    end
                    TimeDelay = TimeDelay + this.seq.blockDurations(end);
                end
        
                for i=1:Ny_meas
                    this.seq.addBlock(mr.makeLabel('SET','REV', mr_gx.amplitude<0), ...
                        mr.makeLabel('SET','SEG', mr_gx.amplitude<0));

                    if TriggerDelay > TimeDelay && TriggerDelay - TimeDelay <= mr.calcDuration(mr_gx,mr_gyBlipup,mr_adc) - mr_trig.duration 
                        mr_trig.delay = TriggerDelay - TimeDelay;
                        if i==1
                            this.seq.addBlock(mr_trig, mr_gx, mr_gyBlipup    , mr_adc); % Read the first line of k-space with a single half-blip at the end
                        elseif i==Ny_meas
                            this.seq.addBlock(mr_trig, mr_gx, mr_gyBlipdown  , mr_adc); % Read the last line of k-space with a single half-blip at the beginning
                        else
                            this.seq.addBlock(mr_trig, mr_gx, mr_gyBlipdownup, mr_adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        end
                        % fprintf('TR %d, Trigger delay: %.3f ms\n', r, TriggerDelay*1e3);
                    else
                        if i==1
                            this.seq.addBlock(mr_gx, mr_gyBlipup    , mr_adc); % Read the first line of k-space with a single half-blip at the end
                        elseif i==Ny_meas
                            this.seq.addBlock(mr_gx, mr_gyBlipdown  , mr_adc); % Read the last line of k-space with a single half-blip at the beginning
                        else
                            this.seq.addBlock(mr_gx, mr_gyBlipdownup, mr_adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        end
                    end
                    
                    TimeDelay = TimeDelay + this.seq.blockDurations(end);

                    this.seq.addBlock(mr.makeLabel('INC','LIN', 1));  %
                    
                    mr_gx.amplitude = -mr_gx.amplitude;   % Reverse polarity of read gradient
        
                end
                this.seq.addBlock(mr_gzSpoil);
                if rem(Navigator+Ny_meas,2)~=0
                    mr_gx.amplitude = -mr_gx.amplitude;
                end
                this.seq.addBlock(mr.makeDelay(delayTR),mr.makeLabel('INC','SLC', 1)); %
                counter= counter + 1;
            end
            this.seq.addBlock(mr.makeLabel('INC','REP', 1)); %
        end
    otherwise
        disp('-> Acquisition scheme specified is not supported... returned...')
end
disp('-> Finished...')

end


function triggerDelays = calcTriggerDelays(grad0, probeType, probeRadius, signalCutoff)
    %%% smartly segment the gradient for ultrahigh resolution
    % signal decay due to dephasing inside a voxel follows sinc(a) where a is the max phase angle relative to the middle of the voxel.
    x        = 0:0.005:1;
    maxPhase = x(find(sinc(x)>signalCutoff, 1, 'last' ))* pi; % in rads.
    maxK     = maxPhase./ probeRadius;
    
    grad     = grad0.shape(1:2,:) / grad0.gamma;
    dt       = grad0.dt;
    ktraj    = calc_ktraj_from_grad(grad, probeType, false, dt);
    ktrajAbs = sqrt(ktraj(1,:).^2 + ktraj(2,:).^2);
    
    idxCutoff   = 1;
    idx2measure = [];
    while ~isempty(idxCutoff)
        idxCutoff = find(ktrajAbs> maxK,1,'first');
        if ~isempty(idxCutoff)
            idx2measure = [idx2measure idxCutoff];
            grad        = grad(:, idxCutoff:end);
            ktraj       = calc_ktraj_from_grad(grad, probeType, false, dt);
            ktrajAbs    = sqrt(ktraj(1,:).^2+ ktraj(2,:).^2);
        end
    end
    triggerDelays= [0 cumsum(idx2measure-1).* dt]; 
end