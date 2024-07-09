% Author: Sebastian Rosenzweig <sebastian.rosenzweig@skope.ch>
% (c) 2021 Skope Magnetic Resonance Technologies AG

function Prepare(this, varargin)
% Prepare function wrapper

   if nargin > 1
       this.seq_params = varargin{1};
   end
   this.seq_params.doFastestPrescan= false;
   
   switch(this.seq_type)
       case("Local Eddy-current Calibration")
           this.Prepare_local_ec_calib();
    
       case("Local Off-res Calibration")
           this.Prepare_local_offres_calib();

       case("GIRF Calibration")
           this.Prepare_girf_calib();

       case("Gradient Tones")
           this.Prepare_grad_tones();

       case ("2D GRE")
           this.Prepare_2dGre();
           this.isImagingSeq = true;

       case ("2D EPI")
           this.Prepare_2dEpi();
           this.isImagingSeq = true;

       case ("2D EPI MoCo")
           this.Prepare_2dEpi_moco();
           this.isImagingSeq = true;

       case ("2D Spiral")
           this.Prepare_2dSpiral();
           this.isImagingSeq = true;

       case ("2D Spiral MoCo")
           this.Prepare_2dSpiral_moco();
           this.isImagingSeq = true;

       case ("2D Spiral Stitch")
           this.Prepare_2dSpiral_stitch();
           this.isImagingSeq = true;

       otherwise
           error("Unknown seq_type");  
   end
   
   this.seq_params.total_time = sum(this.seq.blockDurations);
end 