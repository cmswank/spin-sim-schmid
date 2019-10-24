%compare to theory, you asshole!
gamma=20378.9;
T=.450;
vmax=500;
phi0=pi/2;%phi0=0.8;
Gy=0.40597*9E-6; %for comparison,  %1 GY G/m = 100 GY G/cm
B1=0.40597;

Ly=40;
wrf=2*pi*1000;
w0=20378.9*0.03;
D3=1.6*T^(-7);
tt=linspace(0,10,1000);

%EDM signal
dn=1E-28;  %neutron EDM [e cm]. 
hbar=6.582122E-16; %hbar [e V s]
gammaE=2*dn/hbar;  %gyroelectric ratio. [rad cm/V/s] 
E0=7.5E4; %V/cm
dE0=2*E0; % E_+ - E_-
wdn=dE0*gammaE; %[rad/s]

J0=besselj(0,gamma*B1/wrf);
SignalEDMJ0=(J0*wdn.*tt);
SignalEDM=(wdn.*tt);
J0x=besselj(0,gamma*B1/wrf-gamma/1.112*B1/wrf);
J1x=besselj(1,gamma*B1/wrf-gamma/1.112*B1/wrf);
J1x3=besselj(1,gamma*B1/wrf);
J1xn=besselj(1,gamma*B1/1.112/wrf);
%THermalization theory correlation functions:
%%%%painful but why not. am I right?!?!
k=1.380658E-16;  %boltz' in cgs
m=3.0160293*2.25*1.66054E-24; %effective mass of helium-3 in SFHe-II.  
tc=m/k/T*D3; %collision time
wn3=w0*J0*(1+1/1.112); %sum of the dressed frequencies. 
wp=wn3*tc; %reduced frequency from thermal collision theory. 
tb=Ly*sqrt(m/k/T);
tcb=tc/tb;
%delta from thermal collision theory:
Delta=sqrt(2)*tcb/sqrt(wp)/(1+wp^2)*(sin(sqrt(wp)*(1+wp)/sqrt(2)/tcb)*(1+2*wp-wp^2)...
                                    +sinh(sqrt(wp)*(1-wp)/sqrt(2)/tcb)*(1-2*wp-wp^2))...
                                    /(cos(sqrt(wp)*(1+wp)/sqrt(2)/tcb)+cosh(sqrt(wp)*(1-wp)/sqrt(2)/tcb));

Svvinf=D3/(1+wp^2);

Svv=Svvinf*(1-Delta);

Sxxinf=tc*Ly^2*tcb^2/wp^2/(1+wp^2);

Sxx=Sxxinf*(1-Delta);

Sxx0=Ly^4/120/D3-1/12*Ly^2*tc;  %thermal collision bounded domain correction

Svv0=D3;

theoryNoise=zeros(length(tt),4);

%neutron Sxx(w=0). 
%Sneutron=gamma^2*w0^2/1.112^4*J1xn^2*6*Ly^3*1.0369/pi^4/vmax;

theoryNoise(:,1)=sqrt(1/2*(1+J0x)^2*(gamma^2*w0^2*(J1x3)^2*Ly^4/120/D3)*Gy^2/wrf^2*tt*sin(phi0)^2 +J1x^2*gamma^2*Gy^2/wrf^2*tt*D3*cos(phi0)^2);

%theoryNoise(:,1)=sqrt(1/4*(1+J0x)^2*(gamma^2*w0^2*(J1x3-J1xn)^2*Ly^4/120/D3)*Gy^2/wrf^2*t*sin(phi0)^2 +J1x^2*gamma^2*Gy^2/wrf^2*t*D3*cos(phi0)^2);

phi0=0;
theoryNoise(:,2)=sqrt(1/2*J1x^2*gamma^2*Gy^2/wrf^2*tt.*Svvinf.*cos(phi0)^2);%+1/2*(1+J0x)^2*J1x*gamma^2*Gy^2/wrf^2*t*Ly^4/120/D3*sin(phi0)^2);
%theoryNoise(:,2)=sqrt(1/2*J1x^2*gamma^2*Gy^2/wrf^2*t.*D3*cos(phi0)^2+1/2*(1+J0x)^2*J1x*gamma^2*Gy^2/wrf^2*t*Ly^4/120/D3*sin(phi0)^2);

%phi=0 noise from w3+wn noise. 
theoryNoise(:,3)=sqrt((1-J0x)^2*(gamma^2*w0^2*(J1x3)^2*Ly^4/120/D3)*Gy^2/wrf^2*tt);
theoryNoise(:,4)=sqrt((1-J0x)^2*(gamma^2*w0^2*(J1x3)^2*Sxx)*Gy^2/wrf^2*tt);
theoryNoiseRun=theoryNoise./sqrt(3E5);
