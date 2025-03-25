function [k,g,s,time] = designVariableDensitySpiral(this, Nitlv, isRotationallyVariant, res, fov, radius, safetyMargin)
% design time optimal variable density spiral using Lustig's method 
% described in Lustig et al, TMI 2008
% 
%	Input:
%		Nitlv	-	Number of interleves
%       isRotationallyVariant  -   Indicates type of solution. false for rotationally
%       invariant solution.
%		res	-	resolution (in m)
%		fov	- 	vector of fov (in m)
%		radius	-	vector of radius corresponding to the fov
%
%
% created by Xiaoping Wu, 3/13/2023
% modified by Jinyuan Zhang, 2025.03.25, Gmax and Smax for the readout gradient and other gradients could be different.

res= 1e3.* res; % mm
fov= 1e2.* fov; % cm

% maxGrad= mr.convert(this.sys.maxGrad,'Hz/m','mT/m','gamma',this.sys.gamma);
% maxSlew= 1e-3.*mr.convert(this.sys.maxSlew,'Hz/m','mT/m','gamma',this.sys.gamma);
maxGrad= this.scanner.maxGrad;
maxSlew= this.scanner.maxSlew;

Gmax= 1e-3.* tpm2gpcm(safetyMargin.* maxGrad); % mT/m -> G/cm
Smax= 1e-3.* tpm2gpcm(safetyMargin.* maxSlew); % T/m/s -> G/cm/ms
T= 1e3* this.sys.gradRasterTime; % ms
ds= 0.0002; % step size for integration. smaller for better numerical accuracy.

%%% basically a copy of vdSpiralDesign() here.
%[k,g,s,time,Ck] = vdSpiralDesign(Nitlv, rv, res, fov, radius, Gmax, Smax, T, ds, 'pchip');

kmax = 5/res; % 1/cm

% if length(radius)<2
% 	error(' radius must be at least length=2');
% end

dr = 1/500/max(fov/Nitlv);
r = 0:dr:kmax;   kmax = max(r);

fov =  interp1(radius*kmax,fov,r,'pchip');
dtheta = 2*pi*dr.*fov/Nitlv;
theta = cumsum(dtheta);

C = r.*exp(1i*theta);
C = [real(C(:)), imag(C(:)), C(:)*0];

[~,time,g,s,k] = minTimeGradient(C,isRotationallyVariant, 0, 0, Gmax, Smax,T,ds,0);

g= 1e3.* gpcm2tpm(g); % mT/m
s= 1e3.* gpcm2tpm(s); % T/m/s
k= 200.* pi.* k; % 1/cm -> rad/m
%================
% plot
t_vec= (1:size(g,1)).* T;
L = t_vec(end);
kmax1= 200.* pi.* kmax;

figure, subplot(2,2,1), plot(k(:,1), k(:,2)); title('k-space'); axis([-kmax1 kmax1 -kmax1 kmax1]);
xlabel('rad/m'); ylabel('rad/m')
subplot(2,2,2), plot(t_vec,g(:,1)); axis([0,L,-maxGrad,maxGrad]); title('gradient waveforms')
hold on, plot(t_vec,g(:,2), 'r'); xlabel('ms');ylabel('mT/m'); hold off
legend('gx', 'gy', 'Location', 'Northwest');
subplot(2,2,3), plot(t_vec, (g(:,1).^2 + g(:,2).^2).^0.5); axis([0 L 0 sqrt(2).*maxGrad]);
hold on, plot(t_vec, maxGrad.*ones(size(t_vec)), 'k--'); hold off
title('gradient magnitude');xlabel('ms');ylabel('mT/m');
legend('|g|', 'max g')
subplot(2,2,4), plot(t_vec(1:end-1), (s(:,1).^2 + s(:,2).^2).^0.5); 
hold on, plot(t_vec(1:end-1), maxSlew.*ones(size(t_vec(1:end-1))), 'k--'); hold off
title('slew-rate magnitude'); xlabel('ms');ylabel('T/m/s'); axis([0 L 0 sqrt(2).*maxSlew]);
legend('|s|', 'max slew')

g= mr.convert(g,'mT/m','Hz/m','gamma',this.sys.gamma);
end