# pulseq-sequences
This is a matlab toolbox for creating MR pulse sequences in the pulseq framework that are compatible with field monitoring and motion tracking using NMR probes. 
The toolbox is built based on pulseq sequences for field monitoring provided by Skope MRT (see below for more detail) and pulseq (see the `pulseq` folder for more detail). 
Skope compatible sequences are all placed under the `pulseq_skope` folder. 
The time optimal spiral gradient design is devised using an implementation generously shared by Michael Lustig at https://people.eecs.berkeley.edu/~mlustig/Software.html. 

If you use the code in our toolbox, please consider citing the following paper:
- Jinyuan Zhang, Zihao Zhang*, Zhentao Zuo, Rong Xue, Yan Zhuo, Cameron Cushing, Alexander Bratch, Edward Auerbach, Andrea Grant, Jing An, Kamil Ugurbil, Xiaoping Wu. "Data stitching for dynamic field monitoring with NMR probes", MRM 2025. In press.



### Copyright & License Notice
This software is copyrighted by the Regents of the University of Minnesota. It can be freely used for educational and research purposes by non-profit institutions and US government agencies only. 
Other organizations are allowed to use this software only for evaluation purposes, and any further uses will require prior approval. The software may not be sold or redistributed without prior approval. 
One may make copies of the software for their use provided that the copies, are not sold or distributed, are used under the same terms and conditions. 
As unestablished research software, this code is provided on an "as is'' basis without warranty of any kind, either expressed or implied. 
The downloading, or executing any part of this software constitutes an implicit agreement to these terms. These terms and conditions are subject to change at any time without prior notice.


# Pulseq sequences for field-monitoring

Includes sequences:

 - Off-resonance and position calibration
 - Local eddy current calibration

Please visit https://github.com/pulseq/pulseq for further information about Pulseq.

## Getting started

### Clone this repository  

    git clone https://github.com/SkopeMagneticResonanceTechnologies/Pulseq-Sequences.git

Add Pulseq as a submodule

    cd .\Pulseq-LocalEddyCurrentCalibration\
    git submodule add https://github.com/pulseq/pulseq
    
### Creating the Pulseq sequence files in MATLAB  

 - Open MATLAB 
 - Run offres_pos_calib.m to create the sequence file for the off-resonance and position calibration
 - Run local_ec_calib.m to create the sequence file for the local eddy current calibration

Note that the polarity of the x-axis has been reversed to be played out correctly on Siemens scanners.
