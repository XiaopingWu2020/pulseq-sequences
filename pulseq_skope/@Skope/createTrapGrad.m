function grad= createTrapGrad(this,basegrad,sf)
% createTrapGrad create a trapezoidal gradient given pulseq trapezoidal gradient structure.
if nargin< 3
    sf = 1;
end
dt= this.sys.gradRasterTime;
amp= sf.* basegrad.amplitude;
rise= basegrad.riseTime./ sf;
flat= basegrad.flatTime./ sf;
fall= basegrad.fallTime./ sf;
%delay= igrad.delay;

% rise
ntrise= round(rise./dt);
ampinc= amp./ ntrise;
grise= ampinc:ampinc:amp;
% flat
ntflat= round(flat./dt);
gflat= amp*ones(1,ntflat);
% fall
ntfall= round(fall./dt);
ampinc= amp./ ntfall;
gfall= amp:-ampinc:ampinc;
% delay
%ntdelay= round(delay./dt);
if amp~= 0
    grad= [0 grise gflat gfall 0]; %
    %grad= circshift(grad,ntdelay,2);
else
    grad= zeros(1,ntrise+ ntflat+ ntfall+ 2);
end
% t_vec0= (0:length(grad0)-1).* dt;
% 
% nt= round(length(grad0).* dt./ tDwell);
% t_vec= (0:nt-1).* tDwell;
% grad= interp1(t_vec0, grad0, t_vec);
% grad(isnan(grad))= 0;
end