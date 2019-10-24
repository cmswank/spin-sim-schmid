xyxyac=zeros(2*Bin-1,Event-1);
xyac=zeros(2*Bin-1,Event-1);
xxac=zeros(2*Bin-1,Event-1);
xxxxac=zeros(2*Bin-1,Event-1);
yyac=zeros(2*Bin-1,Event-1);
for i = 1:Event-1
   xxac(:,i)=xcorr(x(:,i),'biased');
   yyac(:,i)=xcorr(y(:,i),'biased');
   xxxxac(:,i)=xcorr(x(:,i).*x(:,i),'biased');
   xyac(:,i)=xcorr(x(:,i),y(:,i),'biased');
   xyxyac(:,i)=xcov(x(:,i).*y(:,i));%,'biased');
   if mod(i,200)==0
       disp([num2str(round(i/(Event-1)*100)), ' %'])
   end
end
 Xxxxac=mean(xxxxac,2);
 Xxac=mean(xxac,2);
 Xyxyac=mean(xyxyac,2);
 Xyac=mean(xyac,2);
 Yyac=mean(yyac,2);
 XYtest=xcov(Xxac((end+1)/2:end),Yyac((end+1)/2:end),'biased');
 disp('100 %');
 figure
hold on
%plot([-flip(t(1:end-1)./2);t(1:end)./2],max(Xyxyac)/max(Xxac).*Xxac)
%plot([-flip(t(1:end-1));t(1:end)],Yyac)
plot([-flip(t(1:end-1));t(1:end)],Xyxyac)

%plot(real(fft(Xyxyac)))
%plot(real(fft(max(Xyxyac)/max(XYtest).*XYtest)))

plot([-flip(t(1:end-1));t(1:end)],max(Xyxyac)/max(2*Xxac-XYtest).*(Xxac+Yyac-XYtest));
%plot([-flip(t(1:end-1));t(1:end)],Xxxxac)
%plot([-flip(t(1:end-1));t(1:end)],Xyac)
%plot([-flip(t(1:end-1)/2);t(1:end)/2],Xxac.^2)
%legend('XYactual', 'XX double correlation')
%legend('XX', 'YY', 'XYXY', 'XX-YY (error)', 'XX XX', 'XY','<XX>^2')
% 
% 
% 
% Bin = 5001;
% Event = 5000;
% runnum=1;
% 
% datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(runnum),'.dat');
% 
% % import data
%     fileID = fopen(datafile);
%     A = fread(fileID, 'double');
%     fclose(fileID);
% 
%     B = reshape(A, 10, Bin, Event);
% % the following are matrices (time (Bin+1), particles (Event-1))
%     sx1 = squeeze(B(1,:,2:Event));
%     sy1 = squeeze(B(2,:,2:Event));
%     sz1 = squeeze(B(3,:,2:Event));
%     x =  squeeze(B(4,:,2:Event));
%     y =  squeeze(B(5,:,2:Event));
%     z = squeeze(B(6,:,2:Event));
%     vx =  squeeze(B(7,:,2:Event));
%     vy =  squeeze(B(8,:,2:Event));
%     vz = squeeze(B(9,:,2:Event));
%     tlarge = squeeze(B(10,:,2:Event)); 
%     t = squeeze(tlarge(:,1));
%     
%     
  