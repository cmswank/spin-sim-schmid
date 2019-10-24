datmaker=@(noisep,dt,lengtht,tend,ns)(abs(fft(noisep))*dt/sqrt(tend));

S2N=7; %Desired signal to noise. 
T2=400;



disp('|');
disp('|');

noiselen=1.5E6;
t=linspace(0,3*T2,noiselen);
dt=t(2);
% noise must scale with sample time width(dt). 
noisescale=sqrt(dt);

noisep=noisescale*random('norm',0,1/sqrt(dt),noiselen,1);
f=0:1/t(end):1/t(2);  %corresponding frequency resolution
dat=datmaker(noisep,dt,length(t),t(end),noisescale);  %%units of 1/sqrt(Hz)

%convert signal2noise to pure signal at desired ratio. 
Signal_converstion=S2N*rms(noisep); 
B0=.030; %30 mG
gamma=20378.9;
numint=15; %number of time intervals
maxslope=0.001*B0;
%minslope=0.001*B0;
slope=-maxslope+rand(numint,1)*2*maxslope;% slope in each time interval
w0=zeros(length(t),1);
signal=zeros(length(t),1);
FID=zeros(length(t),1);
for n=1:length(t)
    if n==1
        w0(n)=gamma*(B0+slope(1+floor((n-1)./length(t).*numint)).*t(n));%gamma*B0=97 Hz
      integ=0;
    else
        integ=integ+t(2)*trapz(w0(n-1:n));
        w0(n)=w0(n-1)+gamma*slope(1+floor((n-1)./length(t).*numint)).*dt;
    end
%     if mod(n,1000)==0
%     disp(n);
%     end
    signal(n)=cos(integ);
    %what we will measure
end

FID=Signal_converstion*signal.*exp(-t/T2);

time=length(t)/numint/10;%number of time points in each fft interval
freq=(0:1:time-1)./dt/time;%frequency domain
interval=1:time:1500000;%intial time of each fft interval
fpeak=zeros(1,length(interval));
for i=1:length(interval)
    sint=abs(fft(FID(interval(i):interval(i)+time-1)))*dt/sqrt(t(end));
    sint(time/2:end)=0;
    plot(freq,sint)
%     if i==1%length(int)
%         figure
%         plot(t(int(i):int(i)+time-1),FID(int(i):int(i)+time-1));
%     end
    hold on
    [value1,max1]=max(sint);
    
     fpeak(i)=max1;%(max1+max2)/2;
     fpeak(i)=(fpeak(i)-1)./dt/time;
    plot(freq(max1),value1,'r*')
end
figure
plot(t(interval),fpeak,'.');
hold on
% plot(t(interval(1)),fpeak(1),'.','color','green');
% n=2:length(interval);
% plot(t(interval(n)),fpeak(n-1)+(fpeak(n)-fpeak(n-1))/2,'.','color','blue');
%plot(t(int(n)),fpeak(n-1)+(fpeak(n)-fpeak(n-1))/2,'.','color','blue');
% plot(t(int(n)),fpeak(1+floor(int(n)./time))+(fpeak(int(n))-fpeak(1+floor(int(n)./time)))/2,'.');
plot(t(interval),w0(interval)/2/pi,'.','color','red');
%legend('first pt','recovered freq','signal freq')