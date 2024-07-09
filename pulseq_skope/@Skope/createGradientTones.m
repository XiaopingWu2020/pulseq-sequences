function [gx,gy,gz,trigDelay]= createGradientTones(this, A0,f0,tw,basegrad)
%createGradientTones create gradient tones for position encoding in motion
%detection with nmr markers. 
%   Detailed explanation goes here
%

if nargin< 5
    basegrad=[];
end
trigDelay= 0;
dt= this.sys.gradRasterTime;
gamma= this.sys.gamma;
A = mr.convert(A0, 'mT/m', 'Hz/m', 'gamma', gamma); % todo: use axis specific amplitudes

f= round(tw* f0)./tw; % frequency adjustment to ensure harmonics across the duration of tones.
t_vec= dt* (0:round(tw/dt));

% grad tones
gx= A* sin(2*pi*f(1)* t_vec);
gy= A* sin(2*pi*f(2)* t_vec);
gz= A* sin(2*pi*f(3)* t_vec);

gt.gx= gx;
gt.gy= gy;
gt.gz= gz;
gt.t_vec= t_vec;
gt.unit= 'Hz/m';
gt.gamma= gamma;
gt.A0= A0;
gt.f0= f0;
gt.A= A;
gt.f= f;

if ~isempty(basegrad)
    sf= 0.3; % scaling down to leave room for superimposing gradient tones at various specs. 
    bgradx= this.createTrapGrad(basegrad.gx,sf);
    bgrady= this.createTrapGrad(basegrad.gy,sf);
    bgradz= this.createTrapGrad(basegrad.gz,sf);

    % scale so that the max combined gradient would remain unchanged.
    if max(abs(bgradx))~= 0
        bgradx= abs((max(bgradx)- max(gx))./ max(bgradx)).* bgradx;
    end
    if max(abs(bgrady))~= 0
        bgrady= abs((max(bgrady)- max(gy))./ max(bgrady)).* bgrady;
    end
    if max(abs(bgradz))~= 0
        bgradz= abs((max(bgradz)- max(gz))./ max(bgradz)).* bgradz;
    end

    if length(bgradx)> length(gx)
        npts2shift= round(0.5*(length(bgradx)- length(gx)));
        trigDelay= npts2shift* dt;

        gx(length(bgradx))= 0;
        gy(length(bgrady))= 0;
        gz(length(bgradz))= 0;
        gx= circshift(gx, npts2shift, 2);
        gy= circshift(gy, npts2shift, 2);
        gz= circshift(gz, npts2shift, 2);
    else
        bgradx(length(gx))= 0;
        bgrady(length(gy))= 0;
        bgradz(length(gz))= 0;
    end

    gx= gx+ bgradx;
    gy= gy+ bgrady;
    gz= gz+ bgradz;

    gt.gx_spoil= gx;
    gt.gy_spoil= gy;
    gt.gz_spoil= gz; 
end

gt.trigDelay= trigDelay;
save gradtones gt;
disp('-> gradient tones waveforms saved in gradtones.mat...')

end
