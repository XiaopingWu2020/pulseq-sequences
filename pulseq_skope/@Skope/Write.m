% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function Write(this, varargin)
% Write to file
   
   if (nargin == 2)
       name = varargin{1};
   else
       name = this.seq_name; % default name
   end
   
   this.seq.setDefinition('Name', name);
   this.seq.setDefinition('skope_trigChannel', this.seq_params.trigChannel);
   
   if this.isImagingSeq
       fov = this.seq_params.fov; N= this.seq_params.N;
       thickness = this.seq_params.thickness;
       Nslices = this.seq_params.Nslices;
       sliceGap= this.seq_params.sliceGap;
       flipangle= this.seq_params.alpha;

       this.seq.setDefinition('FOV', [fov fov (1+sliceGap)*thickness*Nslices]);
       this.seq.setDefinition('matrixSize', [N N Nslices]);
       this.seq.setDefinition('nSlices', Nslices);
       this.seq.setDefinition('Slice_Thickness', thickness);
       this.seq.setDefinition('resolution', fov/N);
       this.seq.setDefinition('flipAngle', flipangle);
   end

   if isfield(this.seq_params,'nRepeats')
       this.seq.setDefinition('nRepetitions', this.seq_params.nRepeats);
   end

   % skope
   if isfield(this.seq_params,'delay_trig2blip')
       this.seq.setDefinition('skope_girfDelay', this.seq_params.delay_trig2blip);
   end
   if isfield(this.seq_params,'trigChannelMoco')
       this.seq.setDefinition('skope_trigChannel_moco', this.seq_params.trigChannelMoco);
   end
   if isfield(this.seq_params,'skopeMinTR')
       this.seq.setDefinition('skope_minimalTR', this.seq_params.skopeMinTR);
   end
   if isfield(this.seq_params,'skopeInterleaveTR')
       this.seq.setDefinition('skope_interleaveTR', this.seq_params.skopeInterleaveTR);
   end
   if isfield(this.seq_params,'skopeAcqDuration')
       this.seq.setDefinition('skope_acqDuration', this.seq_params.skopeAcqDuration);
   end
   if isfield(this.seq_params,'skopeAcqDurationMoco')
       this.seq.setDefinition('skope_acqDuration_moco', this.seq_params.skopeAcqDurationMoco);
   end
   if isfield(this.seq_params,'gradFreeDelay')
       this.seq.setDefinition('skope_gradFreeDelay', this.seq_params.gradFreeDelay);
   end
   if isfield(this.seq_params,'trigger2AdcTimePhaseCorrection')
       this.seq.setDefinition('skope_trig2AdcTimePhaseCorrection', this.seq_params.trigger2AdcTimePhaseCorrection);
   end
   if isfield(this.seq_params,'trigger2AdcTime')
       this.seq.setDefinition('skope_trig2AdcTime', this.seq_params.trigger2AdcTime);
   end
   if isfield(this.seq_params,'trigger2fieldCtrlTime')
       this.seq.setDefinition('skope_trig2fieldControlTime', this.seq_params.trigger2fieldCtrlTime);
   end
   if isfield(this.seq_params,'fieldControlDuration')
       this.seq.setDefinition('skope_fieldControlDuration', this.seq_params.fieldControlDuration);
   end
   if isfield(this.seq_params,'trigger2toneTime')
       this.seq.setDefinition('skope_trig2gradToneTime', this.seq_params.trigger2toneTime);
   end
   if isfield(this.seq_params,'gradToneDuration')
       this.seq.setDefinition('skope_gradToneDuration', this.seq_params.gradToneDuration);
   end
   if isfield(this.seq_params,'nDynamics')
       this.seq.setDefinition('skope_nDynamics', this.seq_params.nDynamics);
   end
   if isfield(this.seq_params,'nPrescans')
       this.seq.setDefinition('skope_nPrescans', this.seq_params.nPrescans);
   end
   if isfield(this.seq_params,'nInterleaves')
       this.seq.setDefinition('skope_nInterleaves', this.seq_params.nInterleaves);
   end
   if isfield(this.seq_params,'nSegments2measure')
       this.seq.setDefinition('skope_nSegments2measure', this.seq_params.nSegments2measure);
   end
   if isfield(this.seq_params,'segmentDuration')
       this.seq.setDefinition('skope_maxSegmentDuration', this.seq_params.segmentDuration);
   end
   if isfield(this.seq_params,'stitchMode')
       this.seq.setDefinition('skope_stitchMode', this.seq_params.stitchMode);
   end
   if isfield(this.seq_params,'triggerDelays')
       this.seq.setDefinition('skope_triggerDelays', this.seq_params.triggerDelays);
   end
   
   if isfield(this.seq_params,'echoSpacing')
       this.seq.setDefinition('echoSpacing', this.seq_params.echoSpacing);
   end
   if isfield(this.seq_params,'prePhaserDuration')
       this.seq.setDefinition('prePhaserDuration', this.seq_params.prePhaserDuration);
   end
   if isfield(this.seq_params,'TR')
       this.seq.setDefinition('TR', this.seq_params.TR);
   end
   if isfield(this.seq_params,'TE')
       this.seq.setDefinition('TE', this.seq_params.TE);
   end
   if isfield(this.seq_params,'isSpiral')
       this.seq.setDefinition('isSpiral', this.seq_params.isSpiral);
   end
   if isfield(this.seq_params,'echoTrainLength')
       this.seq.setDefinition('echoTrainLength', this.seq_params.echoTrainLength);
   end
   if isfield(this.seq_params,'nContrasts')
       this.seq.setDefinition('nContrasts', this.seq_params.nContrasts);
   end
   if isfield(this.seq_params,'nSpokes')
       this.seq.setDefinition('nSpokes', this.seq_params.nSpokes);
   end
   if isfield(this.seq_params,'readoutGradientDuration')
       this.seq.setDefinition('readoutGradientDuration', this.seq_params.readoutGradientDuration);
   end

   % ADC
   if isfield(this.seq_params,'adcLength')
       this.seq.setDefinition('adc_length', this.seq_params.adcLength);
   end
   if isfield(this.seq_params,'maxAdcSegmentLength')
       this.seq.setDefinition('adc_segmentLength', this.seq_params.maxAdcSegmentLength);
   end
   if isfield(this.seq_params,'nAdcSegments')
       this.seq.setDefinition('adc_nSegments', this.seq_params.nAdcSegments);
   end
   if isfield(this.seq_params,'adcDwellTime')
       this.seq.setDefinition('adc_DwellTime', this.seq_params.adcDwellTime);
   end
   if isfield(this.seq_params,'adcSegmentDuration')
       this.seq.setDefinition('adc_segmentDuration', this.seq_params.adcSegmentDuration);
   end
   if isfield(this.seq_params,'effectiveAdcTime')
       this.seq.setDefinition('adc_effectiveTime', this.seq_params.effectiveAdcTime);
   end

   if isfield(this.seq_params,'useSingleAdcSegment4Sync')
       this.seq.setDefinition('useSingleAdcSegment4Prescan', this.seq_params.useSingleAdcSegment4Sync);
   end
   if isfield(this.seq_params,'doFastestPrescan')
       this.seq.setDefinition('doFastestPrescan', this.seq_params.doFastestPrescan);
   end
   if isfield(this.seq_params,'accelerationFactor')
       this.seq.setDefinition('accelerationFactor', this.seq_params.accelerationFactor);
   end
   if isfield(this.seq_params,'partialFourier')
       this.seq.setDefinition('partialFourier', this.seq_params.partialFourier);
   end
   if isfield(this.seq_params,'reordering')
       this.seq.setDefinition('reordering', this.seq_params.reordering);
   end

   this.seq.write(join([name, '.seq'],''));

end 