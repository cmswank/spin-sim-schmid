%from parameter file
width=8;
mu  = 5;
beta = width*pi/mu;
expo=(1+mu*1i);
B0 = 3e-6;
w_mean = B0*(1.83247185e8 + 2.037894730e8)/2;
gamma_n=1.83247172E8; % rad/s/T
gamma_3=2.037894585E8;
w_n = B0*gamma_n;
w_He = B0*gamma_3;
rscale=1.053845;
Bscale=1.9745e-3*3.8750592e-05;
t=linspace(-1.4,1.4,1E5);
zeropad=1E7;
%B1pulse=Bscale.*((rscale.*exp(1i.*w_n.*t).*cosh(beta.*t).^expo+1/rscale.*exp(1i.*w_He.*t).*cosh(beta.*t).^expo));
B1pulse=real(Bscale.*((rscale.*exp(1i.*w_n.*t).*sech(beta.*t).^expo+1/rscale.*exp(1i.*w_He.*t).*sech(beta.*t).^expo)));
S1pulse=fft(B1pulse,zeropad);
S1pulse=abs(S1pulse)./abs(S1pulse(27160));
%f=1/t(end):1/t(end):length(t)/t(end);
deltat=t(end)-t(end-1);
f=1/deltat/zeropad:1/deltat/zeropad:1/deltat; %zeropadded f

%%%Tailored pulses
%say delta f is 2 Hz, so delta t of pulse is 0.5 seconds.....  0.5/2.8;
%dressing factor decreases by ~6 its not clear how much the phase noise
%increases 
% f_3 = 97.3
% f_n=  87.5

%%%%% OPIMIZED PULSE.
%%%%% UNDRESSED PULSE. 
%%%%% Bcode = 0.0153505815 (code input*dress mag = B1 mag)
      %B1mag= 0.0153505815*3.8750592e-5
      %B1mag= 0.00594844120669248 G or in Tesla B1= 5.948441206692480e-07 T
      %b1alpha = 1.066615610224131  ratio of spectrum height to n vs He3
      %(ie b1alpha= .959099139*gamma_3/gamma_n;)
      b1alpha=.959099139*gamma_3/gamma_n; 
%%%Nominal pulse.
%STpulse=[zeros(1,32),b1alpha*1./(1+exp(-3*linspace(-5,5,10))),b1alpha*ones(1,4),...
%        ones(1,4),1./(1+exp(-3*linspace(5,-5,10))),zeros(1,32)];
%plength=0.5;

%%%%robust dressed spin tip pulse. 
STpulse=[zeros(1,32),1./(1+exp(-3*linspace(-5,5,10))),ones(1,4),...
        ones(1,4),1./(1+exp(-3*linspace(5,-5,10))),zeros(1,32)];  
plength=0.78;


zeropad=5E4;   
BTpulse=real(fftshift(ifft(STpulse,zeropad)));    


fp=1/plength:1/plength:length(STpulse)/plength;

dt=plength/zeropad:plength/zeropad:plength;
dt=dt';
ft=1/dt(1)/zeropad:1/dt(1)/zeropad:length(dt)/dt(1)/zeropad;

STactual=abs(fft(BTpulse));

%Normalize pulse for working with stuff. 

%BTpulse=BTpulse'./max(BTpulse);

%Crazy idea for robust pulse, (there is a zero frequency component...)
BTpulse=ones(size(BTpulse));

fileID = fopen('/data1/cmswank/spin-sim-xliu/BField/B1Pulse0.dat','w');
fwrite(fileID,BTpulse,'double');
%wow shit is like really really good...







