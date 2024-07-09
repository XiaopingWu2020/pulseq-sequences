function gos= designArchimedean(this,deltak,kRadius,kSamples,safety_margin)
% define k-space parameters
%readoutTime = 4.2e-4;

% calculate a raw Archimedian spiral trajectory
% if nargin< 5
%     safety_margin= 0.94;
% end
clear ka;
ka(kRadius*kSamples+1)=1i; % init as complex
for c=0:kRadius*kSamples
    r=deltak*c/kSamples;
    a=mod(c,kSamples)*2*pi/kSamples;
    ka(c+1)=r*exp(1i*a);
end
ka=[real(ka); imag(ka)];
% calculate gradients and slew rates
[ga, sa]=mr.traj2grad(ka);

% limit analysis
%safety_magrin= 0.25; %0.94; % we need that  otherwise we just about violate the slew rate due to the rounding errors
dt_gcomp=abs(ga)/(this.sys.maxGrad*safety_margin)*this.sys.gradRasterTime;
dt_gabs=abs(ga(1,:)+1i*ga(2,:))/(this.sys.maxGrad*safety_margin)*this.sys.gradRasterTime;
dt_scomp=sqrt(abs(sa)/(this.sys.maxSlew*safety_margin))*this.sys.gradRasterTime;
dt_sabs=sqrt(abs(sa(1,:)+1i*sa(2,:))/(this.sys.maxSlew*safety_margin))*this.sys.gradRasterTime;

figure;plot([dt_gabs; max(dt_gcomp); dt_sabs; max(dt_scomp)]');title('time stepping defined by gradient and slew-rate');

dt_smooth=max([dt_gabs;dt_sabs]);
dt_rough=max([dt_gcomp;dt_scomp]);

% apply the lower limit not to lose the trajectory detail
nppr= 4; %8;
dt_min=nppr*this.sys.gradRasterTime/kSamples; % we want at least 4 points per revolution
dt_smooth0=dt_smooth;
dt_rough0=dt_rough;
dt_smooth(dt_smooth<dt_min)=dt_min;
dt_rough(dt_rough<dt_min)=dt_min;

figure;plot([dt_smooth0; dt_smooth; dt_rough0; dt_rough]');title('combined time stepping');

t_smooth=[0 cumsum(dt_smooth,2)];
t_rough=[0 cumsum(dt_rough,2)];

kopt_smooth=interp1(t_smooth, ka', (0:floor(t_smooth(end)/this.sys.gradRasterTime))*this.sys.gradRasterTime)';
kopt_rough=interp1(t_rough, ka', (0:floor(t_rough(end)/this.sys.gradRasterTime))*this.sys.gradRasterTime)';

% analyze what we've got
fprintf('duration orig %d us\n', round(1e6*this.sys.gradRasterTime*length(ka)));
fprintf('duration smooth %d us\n', round(1e6*this.sys.gradRasterTime*length(kopt_smooth)));
fprintf('duration rough %d us\n', round(1e6*this.sys.gradRasterTime*length(kopt_rough)));

[gos, sos]=mr.traj2grad(kopt_smooth);
[gor, sor]=mr.traj2grad(kopt_rough);

figure;plot([gos;abs(gos(1,:)+1i*gos(2,:))]');title('gradient with smooth (abs) constraint')
figure;plot([gor;abs(gor(1,:)+1i*gor(2,:))]');title('gradient with rough (component) constraint')

figure;plot([sos;abs(sos(1,:)+1i*sos(2,:))]');title('slew rate with smooth (abs) constraint')
figure;plot([sor;abs(sor(1,:)+1i*sor(2,:))]');title('slew rate with rough (component) constraint')
figure;plot(kopt_smooth(1,:),kopt_smooth(2,:));title('k trajectory')
end