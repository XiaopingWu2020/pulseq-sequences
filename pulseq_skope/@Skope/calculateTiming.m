function delayTR= calculateTiming(this, minExcitTR)
% calculate TR related timing
skopeInterleaveTRDiffFromTR= 10e-3; %100e-6;
TR= this.seq_params.TR;
Nslices= this.seq_params.Nslices;
%
minTR1= this.seq_params.skopeMinTR+ 100e-6;
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
this.seq_params.skopeInterleaveTR= TR- skopeInterleaveTRDiffFromTR;

nSkipsMinimum= ceil(minTR1./ excitTR)- 1;
nSkips= this.seq_params.nInterleaves- 1;
if nSkips>= nSkipsMinimum
    disp('-> Using the number of interleaves specified...')
    this.seq_params.skopeInterleaveTR= (nSkips+ 1)* excitTR- skopeInterleaveTRDiffFromTR;
else
    disp('-> Number of interleaves specified too small...')
    disp('-> Using the minimum required number of interleaves...')
    this.seq_params.skopeInterleaveTR= (nSkipsMinimum+ 1)* excitTR- skopeInterleaveTRDiffFromTR;
    this.seq_params.nInterleaves= nSkipsMinimum+ 1;
end

this.seq_params.TR= TR;
    
end 