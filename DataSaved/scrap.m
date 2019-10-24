% %Signal1sl=zeros(size(Signal1long));
% for i = 1:size(Signallong1,2)
%      Signal9sl(:,i)=filter(a,b,Signallong9(:,i));
%     Signal7sl(:,i)=filter(a,b,Signallong7(:,i));
% end
% SG2llong=Signal7sl-Signal9sl;
% rmsSG2l=rms(SG2llong,2);
% stdSG2l=std(SG2llong');



clear Signal9s Signal7s SignalGlong;

%Signal1sl=zeros(size(Signal1long));
% for i = 1:size(Signallong1,2)
%      Signal9s(:,i)=smooth(Signallong9(:,i),50);
%      Signal7s(:,i)=smooth(Signallong7(:,i),50);
%      %Signal7sl(:,size(Signallong1,2)-i+1)=smooth(Signallong7(:,i),50);
%      SignalGlong(:,i)=smooth(Signallong7(:,i)-Signallong9(:,i),50)./2;
% end

%Signal1sl=zeros(size(Signal1long));
for i = 1:size(Signal4,2)
     Signal8s(:,i)=smooth(Signal8(:,i),50);
     Signal9s(:,i)=smooth(Signal9(:,i),50);
     %Signal7sl(:,size(Signallong1,2)-i+1)=smooth(Signallong7(:,i),50);
     SignalGlong(:,i)=smooth(Signal8(:,i)-Signal9(:,i),50)./2;
end

%transform to expected nEDM statistics Signal for 1 run. (radiation background free....) 
SG2slong=(Signal8s-Signal9s)./2; %
rmsSG2s=rms(SG2slong,2)./sqrt(3E5);
rmsGlong=rms(SignalGlong,2)./sqrt(3E5);
SG2slong=SG2slong.*sqrt(size(Signal4,2))./sqrt(3E5);
%stdSG2s=std(SG2slong');


% Signal1s=zeros(size(Signal1));
% for i = 1:size(Signal1,2)
%     Signal9s(:,i)=smooth(Signal9(:,i),50);
%     Signal7s(:,i)=smooth(Signal7(:,i),50);
% end
% SG2slong=Signal7s-Signal9s;
% rmsSG2=rms(SG2slong,2);
% stdSG2=std(SG2slong');
% SG2s=mean(SG2slong,2);

%random('norm',0,1,noiselen,1);
% t=linspace(0,10,1000);
% testlin=zeros(1000);
% for i = 1:1000
%     testlin(:,i)=random('norm',0,1,length(t),1)'+0.1*t;
%     testlin(:,i)=smooth(testlin(:,i),20);
% end
% 
% testerlin=rms(testlin,2);