thpts=200; %resolution
Temp=.45;%[0.1:0.1:0.6];%TT(1);%.3;%ones(100,1)*TT;%linspace(0.1,.001,.);%0.3;    %Temperature of he-3 superfluid he-4 solution (low X phonon dominated);

Lx=.102;%7.6;%.2;%10.2;%10.2; !!Changed to SI!!!     %length in cm of dimension to be correlated. 

Grad=1E-6/100;%1E-7;%1E-5/100;  %first number is G/cm, divide by 100 to make it T/m
Efield=7.5E8; %V/m

%This part is cgs...
m=5.007*10^(-24)*2.4;  %effective he-3 mass cgs
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
w0=20378.9*0.030*0.67; %about the effective dressed spin freq
w=ones(100,1)*linspace(80,120,thpts);%20378.9*ones(100,1)*B00';%.03;%ones(100,1)*linspace(0,20378.9*.045,thpts); %angular speed
w1=w0+w;%20378.9*ones(100,1)*B00';%.03;%ones(100,1)*linspace(0,20378.9*.045,thpts); %angular speed
w2=w0-w;%20378.9*ones(100,1)*B00';%.03;%ones(100,1)*linspace(0,20378.9*.045,thpts); %angular speed
fw=w(1,:)/pi/2;                              %frequency (for plotting!)
lx=(1:2:2*size(w,1))'*ones(1,size(w,2));        %fourier transform coefficient indicie (only positive indicies)
q=pi*lx/Lx;                             %spacial wave number
ghe31=sqrt(pi.*b/2)./q.*faddeeva(1i.*sqrt(b/2).*(lam+1i.*w1)./q,1000);  %the G function  
ghe32=sqrt(pi.*b/2)./q.*faddeeva(1i.*sqrt(b/2).*(lam+1i.*w2)./q,1000);  %the G function  
phe31=ghe31./(1-lam.*ghe31);
phe32=ghe32./(1-lam.*ghe32); %transformed probability density 
%for Gx*vxE we have this wierd bessel function, 
%for Gx^2 we don't. 
J0=1;%-4*1i*besselj(0,20378.9./w);
Sxtheory1=sum(J0.*8.*Lx^2./pi^4./lx.^4.*(phe31),1); %!!!!!!this is 8 because only positive lx!!!    %transformed correlation function.  
Sxtheory2=sum(J0.*8.*Lx^2./pi^4./lx.^4.*(phe32),1);
%Sxtheory(1)=Lx^2/12; %we know this.... no need for numerical error.  
Sxtheory=Sxtheory2;%+Sxtheory2;
%correlation function behaviour, no units
dw=2*pi*((w0+1i*fw).*imag(Sxtheory1)+(w0-1i*fw).*imag(Sxtheory2));%+1/12/12^2*Lx^2;

%E-B frequency shift in Hz!!!!
dfm=203789000^2*Efield/3E8^2*Grad*(2*pi*fw.*imag(Sxtheory))/2/pi;%+1/12*Lx^2)/2/pi;

%E-B relaxation in Hz!!!!!
Rfm=2*203789000^2*Efield/3E8^2*Grad*(2*pi*fw.*real(Sxtheory));

%B squared phase in Hz!!!!!
dwB2=203789000^2/2*Grad^2*imag(Sxtheory)/2/pi; %in Hz

%B Squared Relaxation
%T1
%R2=203789000^2/2*Grad^2*real(Sxtheory);
R2=203789000^2/2*Grad^2*(Sxtheory);
%T2 (ONLY USE 7.6 For this!!!!! (assuming other gradients are fine...)
T2=2/R2(1);
T2McGregor=1/(Lx^4*203789000^2*(Grad)^2/120/(D3/100^2));


%E^2 Relaxation
R2vE=203789000^2*Efield/3E8^4*fw.^2.*real(Sxtheory);


taud=pi^2*Lx/D3;
dth=1/2/pi/w(1,end);
Sxtheory(1)=Lx^2/12;%real(Sxtheory(2)+(Sxtheory(2)-Sxtheory(3)));
Rvvtheory=real(ifft(real((w(1,:).^2).*(Sxtheory))));
Rxxth=real(ifft(real(Sxtheory)));
Rvv2=diff(diff(Rxxth));

% 
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
