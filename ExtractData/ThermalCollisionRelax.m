function [T2,fwsd,T2McgregorRf,lam]= ThermalCollisionRelax(wdress,B1,gamma,B0,Temp,Grad,Lx)

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
%  Sxtheory=sum(8.*Lx^2./pi^4./lx.^4.*phe3,1); %!!!!!!this is 8 because only positive lx!!!    %transformed correlation function.  

%w0f=besselj(0,1.27);
%calculating factor for w0f (Bessel J0 when wdress & gamma*B1 >> w0) 
%this is in cgs for some reason. (doesn't matter its a factor)
%w=[3000 6000 1800 2100 2400 4200 10000 18000 30000 ];
% (in tesla) Bz = [1.9102418e-5, 3.8750592e-5, 1.10631580e-5, 1.31017478e-5, 1.51160267e-5, 2.69935194e-5, 6.477662320e-5, 1.167322440e-4, 1.946176552e-4]
%B1=1.10631580e-5*10^4;  %1.9102418e-5*10^4;
%wdress=1800;%  [T21p,fwsdp,T21mcp]= ThermalCollisionRelaxJ0(w(n),B1(n),gamma,B0,T,G1(n),L);
%  T21(n)=T21p(157);
%  disp(T21(n));
%  T21mc(n)=T21mcp(157);

gammaG=gamma/10^4;
B0=B0*10^4;
B1=B1*10^4;
[gammap]=SpinDressingHmatrix(B0,B1,wdress,gammaG,0,0,0);
w0f=gammap/gammaG; % w0 'fix' or factor, it fixes the energy levels....
disp(w0f);
%spin dressing Relaxation (the cos(w*tau) gives the mean with w= [w-w0, w+w0] 
fwsd=fw*w0f;%besselj(0,1.11);
%w1=gamma*B1;
wp=w0f*gamma*B0-wdress;
wm=w0f*gamma*B0+wdress;
ghe3p=sqrt(pi.*b/2)./q.*Faddeeva_w(1i.*sqrt(b/2).*(lam+1i.*wp)./q);  %the G function  
phe3p=ghe3p./(1-lam.*ghe3p);                         %transformed probability density 
ghe3m=sqrt(pi.*b/2)./q.*Faddeeva_w(1i.*sqrt(b/2).*(lam+1i.*wm)./q);  %the G function  
phe3m=ghe3m./(1-lam.*ghe3m);                         %transformed probability density 
phe3=1/2*(phe3p+phe3m); %+phe3m;
Sxtheorysd=sum(8.*Lx^2./pi^4./lx.^4.*phe3,1);


%Sxtheory(1)=Lx^2/12; %we know this.... no need for numerical error.  

%spectrum of x vx correlation function, units are Length^2 (m^2) 
%dw=2*pi*fw.*imag(Sxtheory)+1/12*Lx^2;


%E-B frequency shift in Hz!!!! 
%multiply by 2 when comparing to phase shift from Monte Carlo,  
% found from +E compared to -E. becasue delta E  = 2 E . because it is in terms of the Fourier Transform -inf to +inf of e^(-1iwt). 
%dfm=-2*gamma^2*Efield/3E8^2*Grad*(2*pi*fw.*imag(Sxtheory)+1/12*Lx^2)/2/pi;

%spin dressing E-B frequency shift in Hz!!!!
%beginning factor of 2 because +E compared to -E
%dfmsd=-2*w0f*gamma^2*Efield/3E8^2*Grad*(2*pi*fwsd.*imag(Sxtheorysd)+1/12*Lx^2)/2/pi;


%E-B relaxation in Hz!!!!!
%Rfm=2*gamma^2*Efield/3E8^2*Grad*(2*pi*fw.*real(Sxtheory));

%B squared phase in Hz!!!!!
%dwB2=gamma^2/2*Grad^2*imag(Sxtheory)/2/pi; %in Hzcd cd

%B Squared Spin Dressing Relaxation 
%T1
R2=gamma^2*Grad^2.*(w0f.*real(Sxtheorysd));
%disp(Grad);
%T2 (ONLY USE 7.6 For this!!!!! (assuming other gradients are fine...)
T2=1./R2; %why 2? some reason? 
%T2McGregorRF=1/(Lx^4*gamma^2*(Grad)^2/120/(D3/100^2));
D3=D3*0.01^2; %diffusion in m^2/s
R2McGregorRF=w0f.^2.*Lx.^4.*gamma.^2.*Grad.^2./2/120./D3.*1./(1+wdress.^2.*Lx.^4./D3.^2./pi.^4);
T2McgregorRf=1./R2McGregorRF;
%E^2 Relaxation
%R2vE=gamma^2*Efield/3E8^4*fw.^2.*real(Sxtheory);


%taud=pi^2*Lx/D3;
%dth=1/2/pi/w(1,end);
%Sxtheory(1)=Lx^2/12;%real(Sxtheory(2)+(Sxtheory(2)-Sxtheory(3)));
%Rvvtheory=real(ifft(real((w(1,:).^2).*(Sxtheory))));
%Rxxth=real(ifft(real(Sxtheory)));
%Rvv2=diff(diff(Rxxth));

%Riccardo's code definitions. (slightly different)
% #define NUCLEARMAGNETON 3.1524512326e-14 // MeV / T
% #define EV 1.602176487e-19 // J
% #define HBAR 1.054571628e-34 // J s
% 
% #define NEUTRON_MAGNETICMOMENT = -1.9130427 // relative to nuclear magneton
% #define NEUTRON_GAMMA -1.832471850e8 // rad / s / T
% #define HE3_GAMMA -2.037947093e8 // rad / s / T
% thpts=300; %resolution    %Temperature of he-3 superfluid he-4 solution (low X phonon dominated);
% Lx=2.54;40;10.2;      %length in cm of dimension to be correlated. 
% m=5.007*10^(-24)*2.4;  %effective he-3 mass cgs
% kb=1.3807*10^(-16);    %boltzman constant cgs
% b=m/kb/Temp;         %1/speed^2
% D3=1.6*Temp^(-7);      %measured diffusion coefficient
% tauc=D3*b;           %conversion to exact theory collision time
% lam=1/tauc;            %conversion to rate
% w=ones(100,1)*linspace(0,100*2*pi,thpts); %angular speed
% fw=w(1,:)/pi/2;                              %frequency (for plotting!)
% lx=(1:2:200)'*ones(1,size(w,2));        %fourier transform coefficient indicie (only positive indicies)
% q=pi*lx/Lx;                             %spacial wave number
% ghe3=sqrt(pi.*b/2)./q.*faddeeva(1i.*sqrt(b/2).*(lam+1i.*w)./q,50);  %the G function  
% phe3=ghe3./(1-lam.*ghe3);                         %transformed probability density 
% Sxtheory=sum(8.*Lx^2./pi^4./lx.^4.*phe3,1);     %transformed correlation function.  
% %dw2=2*pi*fw.*imag(Sxtheory)+1/12*Lx^2;
% %hold on;
% %plot(fw,dw/dw(1),'b','LineWidth',2)
% %plot(fw,dw2/dw(1),'r','LineWidth',2)
% %plot(fw,(dw2+dw)/dw(1),'Color',[.6 0 .6],'LineWidth',2)
% 
% 
%plot(TT,gam*B0*imag(Sxtheory)+Lx^2/12); %this is 0 to infinity (no 1/2)
