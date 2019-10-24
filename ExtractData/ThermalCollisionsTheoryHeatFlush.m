thpts=300;    %resolution
Temp=.450;    %Temperature of he-3 superfluid he-4 solution (low X phonon dominated);
Lx=200;       %10.2;      %length in cm of dimension to be correlated. 
vc=10;  %convection current cm/s (in a heatflush). 
m=5.007*10^(-24)*2.4;  %effective he-3 mass cgs
kb=1.3807*10^(-16);    %boltzman constant cgs
b=m/kb/Temp;         %1/speed^2
D3=1.6*Temp^(-7)/2;      %measured diffusion coefficient
tauc=D3*b;           %conversion to exact theory collision time
% lam=1/tauc;            %conversion to rate

w=ones(100,1)*linspace(0,150*2*pi,thpts); %angular speed
fw=w(1,:)/pi/2;                              %frequency (for plotting!)
lx=(1:2:200)'*ones(1,size(w,2));        %fourier transform coefficient indicie (only positive indicies)
lxvc=(-99:2:99)'*ones(1,size(w,2));
q=pi*lx/Lx;                             %spacial wave number
qvc=pi*lxvc/Lx;

% ghe3=sqrt(pi.*b/2)./q.*faddeeva(1i.*sqrt(b/2).*(lam+1i.*w)./q,50);  %the G function  
% phe3=ghe3./(1-lam.*ghe3);                         %transformed probability density 
% Sxtheory=sum(8.*Lx^2./pi^4./lx.^4.*phe3,1);     %transformed correlation function.  
% dw=2*pi*fw.*imag(Sxtheory)+1/12*Lx^2;
%a=.08; %axion length scale
% %thpts=300; %resolution    %Temperature of he-3 superfluid he-4 solution (low X phonon dominated);
% %Lx=40;%10.2;      %length in cm of dimension to be correlated. 
% m=5.007*10^(-24)*2.4;  %effective he-3 mass cgs
% kb=1.3807*10^(-16);    %boltzman constant cgs
% b=m/kb/Temp;         %1/speed^2
% D3=1.6*Temp^(-7);      %measured diffusion coefficient
% tauc=D3*b;           %conversion to exact theory collision time
%we can include a convection velocity in the damping term. 

lamp=1/tauc;            %conversion to rate
lam=1/tauc;       %rate with a convection current. 
lamvc=1/tauc+1i.*vc.*qvc;
mfp=1/sqrt(b).*tauc;

%w=ones(1000,1)*linspace(0,10000*2*pi,thpts); %angular speed
%fw=w(1,:)/pi/2;                              %frequency (for plotting!)
%lx=(1:1:1000)'*ones(1,size(w,2));        %fourier transform coefficient indicie (only positive indicies)
%q=pi*lx/Lx;   %spacial wave number

ghe3=sqrt(pi.*b/2)./q.*faddeeva(1i.*sqrt(b/2).*(lam+1i.*w)./q,50);  %the G function  
phe3=ghe3./(1-lamp.*ghe3);  %transformed probability density 
ghe3vctop=sqrt(pi.*b/2)./abs(qvc).*faddeeva(1i.*sqrt(b/2).*(lam+1i.*w)./abs(qvc),50);
ghe3vc=sqrt(pi.*b/2)./abs(qvc).*faddeeva(1i.*sqrt(b/2).*(lamvc+1i.*w)./abs(qvc),50);  %the G function  
phe3vc=ghe3vctop./(1-lamp.*ghe3vc);                         %transformed probability density 
Sxtheory=sum(8.*Lx^2./pi^4./lx.^4.*phe3,1);     %transformed correlation function.  
%Sxtheoryvc=sum(8.*Lx^2./pi^4./lx.^4.*phe3vc,1);     %transformed correlation function.  
Sxtheoryvc=sum(-Lx^2./pi^2./lxvc.^2.*exp(1i.*qvc./2*Lx).*phe3vc,1);     %transformed correlation function.  



%relaxation

gamma=20378.9;
Gx=3E-5; %1E-3 B0/cm at 100 Hz. 
Rx=gamma^2*Gx^2*Sxtheory;
Rxvc=gamma^2*Gx^2/2*real(Sxtheoryvc);


%SxHF=2*sqrt(1/b)
%Axion Field 1D
% Ax=(a.*exp(-Lx.*(1./a +1i.*q)).*(1+ exp(1i.*Lx.*q)).*(-1i-a.*q+exp(Lx./a).*(1i-a.*q)+exp(1i.*Lx.*q).*(-1i+a.*q)+exp(Lx.*(1./a + 1i.*q)).*(1i+a.*q)))./(2.*Lx.*(1+1i.*a.*q).*(1i+a.*q));
% Sax=sum(2*Ax.*Ax.*phe3,1);
% dw2=2*pi*fw.*imag(Sxtheory)+1/12*Lx^2;
% hold on;
% plot(fw,real(Sax)*real(Sxtheory(2))/real(Sax(2)),'r','LineWidth',2)
%plot(fw,real(Sxtheory),'b','LineWidth',2)
%plot(fw,dw/dw(1),'b','LineWidth',2)
%plot(fw,dw2/dw(1),'r','LineWidth',2)
%plot(fw,(dw2+dw)/dw(1),'Color',[.6 0 .6],'LineWidth',2)



