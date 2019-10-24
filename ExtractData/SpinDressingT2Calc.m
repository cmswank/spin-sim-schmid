function [T2,w]=SpinDressingT2Calc(B1,w,G1x,G0x,T,L)
%Gives T2 from Spin Dressing perturbation theory given various parameters. 
% B1 is the Dressing field strength. in G 
% w is the Dressing field frequency.  in rad/sec
%Gx is the Gradient in G/cm
%T is the temperature in Kelvin
%L is the length of the volume in the z direction in cm. 

if nargin <5 
    L=40;
end
if nargin < 4
    T=0.450;
end
if nargin <2
    w=6000;
end

if nargin < 1
    B1=0.38750592; %Gauss
end

if nargin < 3
    Gx=1E-6; %(1E-6 ppB0/cm)
end


D=1.6*T^(-7);%Diffusion coefficient for he3 in superfluid HE-II %he3 at room T 1.5*1000./P';
%m=5.007*10^(-24)*2.4;  %effective he-3 mass cgs
%kb=1.3807*10^(-16);    %boltzman constant cgs
%b=m/kb./T;
%tc=D.*b;
w1=20378.9*B1;
dw1=20378.9*G1x;
w0=20378.9*0.03;
dw0=20378.9*G0x;

%2nd order perturbation theory. 
R2=((w0.*dw1./w.*besselj(1,w1./w)).^2+dw0.^2.*besselj(0,w1./w).^2-2*besselj(0,w1./w).*besselj(1,w1./w).*w0./w.*dw1.*dw0).*(L.^4./120./D);

T2=1./R2;

%old estimation. 
% R2simp=20378.9.^2.*besselj(0,w1./w).^2*w0.^2/w.^2.*Gx.^2.*L.^4./120./D;
% T2simp=1./R2simp;
%T1=1./(20378.9^2./2.*Gx.^2.*2.*tc./b./w0.^2./(1+w0.^2.*tc.^2));

