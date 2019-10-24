%Noise Scaling Example

%this function generates randomly distributed noise according to a normal
%distribution with mean 0 and std 1.
%it can be used to generte noise for a signal, 

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

%+noisep;
%dat=datmaker(signal,noisep',dt,length(t),t(end),noisescale);  %%units of 1/sqrt(Hz)


%%NOTES ON S2N: Main principle is that the rms of the noise 
%is the integral of the transformed noise squared.
%
% the rms of the noise increases with the sqrt of the sample time width(dt). 
% the rms of the noise does not change with total observation time(t(end)). 

%the transformed noise deviation does not change with sample time width(dt). 
%the transformed noise deviation increases with the sqrt of the total observation
%time(t(end)). 


%Fid Decay decay            100 Hz
fvar=.0033;
tvar=2*pi*fvar;

B0=.030; %30 mG



Bvar=0.0001*B0;

gamma=20378.9;

f0=gamma*B0/2/pi;

Bvariation=Bvar.*t;%*sin(5*t/t(end))/50;%cos(t.*tvar);

w0=gamma*(B0+Bvariation);%gamma*B0=97 Hz

signal=cos(t.*w0);

FID=Signal_converstion*signal;%.*exp(-t/T2);%what we will measure


Vout=FID'+noisep;

% figure(22) 
% hold off
% plot(t,FID)
% hold on
% plot(t,noisep)
% figure(23)
% plot(t,noisep+FID')
% 
% figure(24)
% hold off
% 
% plot(f(1:200000),dat(1:200000));
% sfid=real(fft(FID))*dt/sqrt(t(end));
% hold on
% plot(f(1:200000),sfid(1:200000));
% hold off
% figure
% svout=real(fft(Vout))*dt/sqrt(t(end));
% plot(f(1:200000),svout(1:200000));

time=12500;%number of time points in each fft interval
freq=(0:1:time-1)./dt/time;%frequency domain
int=1:time:1500000;%intial time of each fft interval
fpeak=zeros(1,length(int));
for i=1:length(int)
    sint=abs(fft(Vout(int(i):int(i)+time-1)))*dt/sqrt(t(end));
    plot(freq,sint)
    if i==length(int)
        plot(t(int(i):int(i)+time-1),Vout(int(i):int(i)+time-1));
    end
    hold on
    [m,fpeak(i)]=max(sint);
    fpeak(i)=(fpeak(i)-1)./dt/time;
end
figure
plot(t(int),(fpeak-f0)/2+f0,'.');
hold on
plot(t(int),w0(int)/2/pi,'.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%find peaks of Vout
% [pks,locs]=findpeaks(FID);
% f1=zeros(1,length(locs)-1);%frequency we recovered
% for i = 1:length(locs)-1
%    f1(i)=1./((locs(i+1)-locs(i)).*dt);
% end
% figure
% plot(locs(1:length(locs)-1).*dt,f1,'.')
% figure
% plot(locs(1:length(locs)-1).*dt,f1)
% 
% [pks,locs]=findpeaks(Vout);
% k=1;
% f1=zeros(1,length(locs)-1);%frequency we recovered
% while k<length(locs)
%    while k<length(locs) && pks(k)<0
%        locs(k)=[];
%        pks(k)=[];
%        f1(k)=[];
%    end
%    if k>=length(locs)
%        break;
%    elseif locs(k+1)-locs(k)<=5 
%        if locs(k+1)>locs(k)
%            locs(k)=[];
%            pks(k)=[];
%        else
%            locs(k+1)=[];
%            pks(k+1)=[];
%        end
%        f1(k)=[];
%    elseif locs(k+1)-locs(k)>=17 %this interval cannot be used since it's not one cycle
%        f1(k)=[];
%        locs(k)=[];
%        pks(k)=[];
%    else
%        f1(k)=1./((locs(k+1)-locs(k)).*dt);
%        k=k+1;
%    end
% end
% 
% figure
% plot(t(locs(:)),w0(locs(:)));
% hold on
% plot(locs(1:length(locs)-1).*dt,f1)
% hold off
% figure
% plot(t,Vout)
% hold on
% plot(locs.*dt,Vout(locs))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% disp(['a ',num2str(dt),' s sample time for a total time of ', num2str(t(end)), 's']);
% disp(['the rms of the time sequence noise is ', num2str(rms(noisep))])
% disp(['the std of the transformed noise is ',num2str(std(dat))]);
% disp(['sqrt of the integral of transformed noise squared is  ',num2str(sqrt(sum(dat.^2)*f(2)))])
% 
% disp('|');
% disp('|');
% 
% noiselen=1E6;
% t=linspace(0,1,noiselen);
% dt=t(2);
% noisep=noisescale*random('norm',0,1/sqrt(dt),noiselen,1);
%   %this is the time samples. 
% f=0:1/t(end):1/t(2);  %corresponding frequency resolution
% %dat=real(2*fft(noisep)*dt)/sqrt(t(end));  %%units of 1/sqrt(Hz)
% dat=datmaker(noisep,dt,length(t),t(end),noisescale);  %%units of 1/sqrt(Hz)
% 
% disp(['a ',num2str(dt),' s sample time, for a total time of ', num2str(t(end)), 's']);
% disp(['the rms of the time sequence noise is ', num2str(rms(noisep))])
% disp(['the std of the transformed noise is ',num2str(std(dat))]);
% disp(['sqrt of the integral of transformed noise squared is  ',num2str(sqrt(sum(dat.^2)*f(2)))])
% 
% disp('|');
% disp('|');
% 
% noiselen=1E5;
% t=linspace(0,1,noiselen);
% dt=t(2);
% noisep=noisescale*random('norm',0,1/sqrt(dt),noiselen,1);
%   %this is the time samples. 
% f=0:1/t(end):1/t(2);  %corresponding frequency resolution
% dat=datmaker(noisep,dt,length(t),t(end),noisescale);  %%units of 1/sqrt(Hz) %%units of 1/sqrt(Hz)
% 
% disp(['a  ',num2str(dt),' s sample time, for a total time of ', num2str(t(end)), 's']);
% disp(['the rms of the time sequence noise is ', num2str(rms(noisep))])
% disp(['the std of the transformed noise is ',num2str(std(dat))]);
% disp(['sqrt of the integral of transformed noise squared is  ',num2str(sqrt(sum(dat.^2)*f(2)))])
% disp('|');
% disp('|');
