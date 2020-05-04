
%Width=1/2; %1/width must be even greater than 2. 
t=linspace(0,1,100000);
w0=2*pi*100;
%tp=mod(t,pi/w0/J0);
%Ft=(heaviside(w0/Width*J0*tp-(pi/Width/2-pi))-heaviside(w0/Width*J0*tp-(pi/Width/2+pi))).*sin(w0/Width*J0*tp);
dw=5384.88986487;
wrf=6283.18530718;
%fm=J0*w0/Width;

fm=753.982236862;

Phint=dw/fm*(1-cos(fm*t));
Fopt=cos(Phint+wrf*t+pi/2);
Fcos=cos(wrf*t);
figure
%plot(t,Fopt)
% hold on
% plot(t,Fcos)

