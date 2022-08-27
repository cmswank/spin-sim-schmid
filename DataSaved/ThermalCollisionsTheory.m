thpts=2000; %resolution
Temp=.3;%[0.2:0.1:0.6];%TT(1);%.3;%ones(100,1)*TT;%linspace(0.1,.001,.);%0.3;    %Temperature of he-3 superfluid he-4 solution (low X phonon dominated);

Lx=.4;%7.6;%.2;%10.2;%10.2; !!Changed to SI!!!     %length in cm of dimension to be correlated. 
gyro=203789000;

Grad=1E-5/100;%1E-7;%1E-5/100;  %first number is G/cm, divide by 100 to make it T/m
Efield=7.5E8; %V/m

%This part is cgs...
m=5.007*10^(-24)*2.4;  %effective he-3 mass cgs
kb=1.3807*10^(-16);    %boltzman constant cgs
b=m/kb./Temp;         %1/speed^2
D3=1.6*Temp^(-7);%ones(100,1)*D33;% 1.6*Temp^(-7)/2;%     %measured diffusion coefficient
tauc=D3.*b;           %conversion to exact theory collision time
lam=1./tauc;          %conversion to rate
%%%comes out as a rate so its OK to use cgs here.
%%%BUT!!!!!!! We use b later so we MUST convert
b=b*100^2; %Convert speed to correct units!

%%%%%%%%%%Neutrons or Helium-3 Read Below

%gyro=gyro/1.112;  % !!!!!!!!!!Comment out for 3HE  !!!!!!!!!!!!!!!
%b=1/5.2; %ucn spectrum approximation. !!!!!!!!!!Comment out for 3HE  !!!!!!!!!!!!!!!
%lam=1; %ucn approximation. !!!!!!!!!!Comment out for 3HE  !!!!!!!!!!!!!!!
 %lam = 1 is about diffuse bounce probability ~ 5%

%%%%Neutrons or Helium-3 END.
 gyroVE=gyro;
% if Lx >.3
%     gyro=gyro*besselj(0,gyro*4E-5/2/pi/1000);
% end
% 
% if Lx <.3
%     gyroVE=gyro*besselj(0,gyro*4E-5/2/pi/1000);
% end

w=ones(100,1)*linspace(0,2*pi*2000,thpts);%20378.9*ones(100,1)*B00';%.03;%ones(100,1)*linspace(0,20378.9*.045,thpts); %angular speed
fw=w(1,:)/pi/2;                              %frequency (for plotting!)
lx=(1:2:2*size(w,1))'*ones(1,size(w,2));        %fourier transform coefficient indicie (only positive indicies)
q=pi*lx/Lx;                             %spacial wave number
ghe3=sqrt(pi.*b/2)./q.*faddeeva(1i.*sqrt(b/2).*(lam+1i.*w)./q,10000);  %the G function  
phe3=ghe3./(1-lam.*ghe3);                         %transformed probability density 
Sxtheory=sum(8.*Lx^2./pi^4./lx.^4.*phe3,1); %!!!!!!this is 8 because only positive lx!!!    %transformed correlation function.  
%Sxtheory(1)=Lx^2/12; %we know this.... no need for numerical error.  

%goofball gyro logic to get the correct bessel function factor;

%correlation function behaviour, no units
dw=2*pi*fw.*imag(Sxtheory)+1/12*Lx^2;
vvr=-2^2*pi^2*fw.^2.*real(Sxtheory);
vv=-2^2*pi^2*fw.^2.*imag(Sxtheory);
%E-B frequency shift ("geometric phase") in Hz!!!! 
dfm=gyro.^2.*Efield/3E8^2*Grad*(2*pi*fw.*imag(Sxtheory)+1/12*Lx^2)/2/pi;

%E^2 frequency shift in Hz!!!
dE2=gyroVE^2/2*Efield^2/3E8^4.*(2^2*pi^2*fw.^2.*imag(Sxtheory)+2*pi*fw/12*Lx^2)/2/pi;

%dE2p=gyroVE^2/2*Efield^2/3E8^4.*(2^2*pi^2*fw.^2.*1i.*imag(Sxtheory)+2*pi*fw.*1i/12*Lx^2);




Svi=(2^2*pi^2*fw.^2.*imag(Sxtheory)+2*pi*fw/12*Lx^2);
%E-B relaxation in Hz!!!!!
Rfm=2*gyro^2*Efield/3E8^2*Grad*(2*pi*fw.*real(Sxtheory));

%B squared phase in Hz!!!!!
dB2=gyro^2/2*Grad^2*imag(Sxtheory)/2/pi; %in Hz

%B Squared Relaxation in Hz!!
%T1
R2=gyro^2*Grad^2*real(Sxtheory);
%T2 (ONLY USE 7.6 For this!!!!! (assuming other gradients are fine...)
T2=2/R2(1);
T2McGregor=1/(Lx^4*gyro^2*(Grad)^2/120/(D3/100^2));


%E^2 Relaxation
R2vE=gyroVE^2*Efield^2/3E8^4*2.^2*pi.^2*fw.^2.*real(Sxtheory);


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
