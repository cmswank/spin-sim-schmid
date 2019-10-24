function [Sxtheorysd,fwsd,SxxMcgregor]= ConditionalDensityTest(wdress,B1,gamma,B0,Temp,Grad,Lx)

thpts=200; %resolution
%Temp=.300;%[0.1:0.1:0.6];%TT(1);%.3;%ones(100,1)*TT;%linspace(0.1,.001,.);%0.3;    %Temperature of he-3 superfluid he-4 solution (low X phonon dominated);

%Lx=0.4;%7.6;%.2;%10.2;%10.2; !!Changed to SI!!!     %length in cm of dimension to be correlated. 

%Grad=2E-6;%1E-4/100;%1E-7;%1E-5/100;  %(T/m) first number is G/cm, divide by 100 to make it T/m
Efield=7.5E7; %V/m

%gamma=2.037947093e8;
%This part is cgs...
%
m=5.007*10^(-24)*2.2;  %effective he-3 mass cgs
kb=1.3807*10^(-16);    %boltzman constant cgs
b=m/kb./Temp;         %1/speed^2
D3=1.6*Temp^(-7);%ones(100,1)*D33;% 1.6*Temp^(-7)/2;%     %measured diffusion coefficient
tauc=D3.*b;           %conversion to exact theory collision time
lam=1./tauc;            %conversion to rate
%%%comes out as a rate so its OK to use cgs here.
%%%BUT!!!!!!! We use b later so we MUST convert
b=b*100^2; %Convert speed to correct units!

%%%%%%%%%%Neutrons or Helium-3 Read Below

%b=1/9; %ucn spectrum approximation. !!!!!!!!!!Comment out for 3HE  !!!!!!!!!!!!!!!
%lam=0; %ucn approximation. !!!!!!!!!!Comment out for 3HE  !!!!!!!!!!!!!!!

%%%%Neutrons or Helium-3 END.

 w=ones(100,1)*linspace(0,800,thpts);%20378.9*ones(100,1)*B00';%.03;%ones(100,1)*linspace(0,20378.9*.045,thpts); %angular speed
 fw=w(1,:)/pi/2;                              %frequency (for plotting!)
 lx=(1:2:2*size(w,1))'*ones(1,size(w,2));        %fourier transform coefficient indicie (only positive indicies)
 q=pi*lx/Lx;                             %spacial wave number
% ghe3=sqrt(pi.*b/2).w*w0f/q.*Faddeeva_w(1i.*sqrt(b/2).*(lam+1i.*w)./q);  %the G function  
%  phe3=ghe3./(1-lam.*ghe3);                         %transformed probability density 
%  Sxtheory=sum(8.*Lx^2./pi^4./lx.^4.5*phe3,1); %!!!!!!this is 8 because only positive lx!!!    %transformed correlation function.  

%w0f=besselj(0,1.27);
%calculating factor for w0f (Bessel J0 when wdress & gamma*B1 >> w0) 
%this is in cgs for some reason. (doesn't matter its a factor)
%w=[3000 6000 1800 2100 2400 4200 10000 18000 30000 ];
% (in tesla) Bz = [1.9102418e-5, 3.8750592e-5, 1.10631580e-5, 1.31017478e-5, 1.51160267e-5, 2.69935194e-5, 6.477662320e-5, 1.167322440e-4, 1.946176552e-4]
%B1=1.10631580e-5*10^4;  %1.9102418e-5*10^4;
%wdress=1800;
gammaG=gamma/10^4;
B0=B0*10^4;
B1=B1*10^4;
[gammap]=SpinDressingHmatrix(B0,B1,wdress,gammaG,0,0,0);
w0f=gammap/gammaG; % w0 'fix' or factor, it fixes the energy levels....

disp(w0f);
%disp(besselj(0,gammaG*B1/wdress));
%w0f=besselj(0,gammaG*B1/wdress);
%spin dressing Relaxation (the cos(w*tau) gives the mean with w= [w-w0, w+w0] 
fwsd=fw*w0f;%besselj(0,1.11);
%w1=gamma*B1;
wp=wdress+w0f*w;
wm=-wdress+w0f*w;
wpp=wdress-w0f*w;
wmm=-wdress-w0f*w;
wp=w;
ghe3p=sqrt(pi.*b/2)./q.*Faddeeva_w(1i.*sqrt(b/2).*(lam+1i.*wp)./q);  %the G function
ghe3pp=sqrt(pi.*b/2)./q.*Faddeeva_w(1i.*sqrt(b/2).*(lam+1i.*wpp)./q);
phe3p=ghe3p./(1-lam.*ghe3p);    
phe3pp=ghe3pp./(1-lam.*ghe3pp);     %transformed probability density 
ghe3m=sqrt(pi.*b/2)./q.*Faddeeva_w(1i.*sqrt(b/2).*(lam+1i.*wm)./q);  %the G function  
ghe3mm=sqrt(pi.*b/2)./q.*Faddeeva_w(1i.*sqrt(b/2).*(lam+1i.*wmm)./q);  %the G function  
phe3m=ghe3m./(1-lam.*ghe3m);                         %transformed probability density 
phe3mm=ghe3mm./(1-lam.*ghe3mm);                         %transformed probability density 
phe3=2*pi*phe3p;%2*pi*(phe3p+phe3m);%don't use the other ones they are zero!  %transform of cos(wt) 
Sxtheorysd=sum(8.*Lx^2./pi^4./lx.^4.*phe3,1);

w0f=gammap/gammaG;
R2=real(Sxtheorysd);
%disp(Grad);
%T2 (ONLY USE 7.6 For this!!!!! (assuming other gradients are fine...)
T2=1./R2; %why 2? some reason? 
D3=D3*0.01^2; %diffusion in m^2/s
%diffusion theory T2 for non B0 Spin Dressing rf relaxation. 
%does it even have a j0?
SxxMcgregor=Lx^4/60/D3.*(1./(1+wp.^2.*Lx.^4./D3.^2./pi.^4));%+1./(1+wm.^2.*Lx.^4./D3.^2./pi.^4));
R2McGregorRF=Lx^4/60/D3.*(1./(1+wm.^2.*Lx.^4./D3.^2./pi.^4));%+1./(1+wp.^2.*Lx.^4./D3.^2./pi.^4));
T2McgregorRf=1./R2McGregorRF;
