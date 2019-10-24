
k=9; %number of frequencies
temperature= 0.450;%K
B0=3E-6;
G=1E-3*3E-6;
L=.40;
w=[3000 6000 1800 2100 2400 4200 10000 18000 30000 ]; %all the frequencies
B1=[0.191024180 0.387505920 0.110631580 0.131017478 0.151160267 0.269935194 0.647766232 1.167322440 1.946176552];
a = [0.73, 0.726, 0.7395, 0.7357, 0.733, 0.728, 0.726, 0.7257, 0.7256];

start=1840;
pfn = '/home/xliu5/spin-sim-SD/parSpinDressing'; 	%parameter file name
clear delphi1 delphi2 delphi dphi;
delphi1{1}=0;
delphi2{1}=0;
dfsim=zeros(k,1);
dfsimdev=zeros(k,1);
dfmsd=zeros(k,1);

%for runs 1800 to 1817
%nBins = 20000;
%nEntries = 5000;  

Bin = 1001;
Event = 1000;

for n=1:k
    datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(start+(n-1)*2),'.dat');
    %spin up with the gradient

    

% import data
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);

    B = reshape(A, 10, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
    sx1 = squeeze(B(1,:,2:Event));
    sy1 = squeeze(B(2,:,2:Event));
    sz1 = squeeze(B(3,:,2:Event));
    x =  squeeze(B(4,:,2:Event));
    y =  squeeze(B(5,:,2:Event));
    z = squeeze(B(6,:,2:Event));
    vx =  squeeze(B(7,:,2:Event));
    vy =  squeeze(B(8,:,2:Event));
    vz = squeeze(B(9,:,2:Event));
    tlarge = squeeze(B(10,:,2:Event)); 
    t = squeeze(tlarge(:,1));

    phi_up=angle(sy1+1i*sz1);
    
    datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(start+(n-1)*2+1),'.dat');
    %spin down wi%   for n=1:8   
%    dfmsdp=ThermalCollisionShifts(w(n), B1(n), 2.037947093e8, 3e-6, temperature); 
%    dfmsd(n)=dfmsdp(157);
%   end
%th the gradient
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);

    B = reshape(A, 10, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
    sx2 = squeeze(B(1,:,2:Event));
    sy2 = squeeze(B(2,:,2:Event));
    sz2 = squeeze(B(3,:,2:Event));
    %x =  squeeze(B(4,:,2:Event));
    %y =  squeeze(B(5,:,2:Event));
    %z = squeeze(B(6,:,2:Event));
    %vx =  squeeze(B(7,:,2:Event));
    %vy =  squeeze(B(8,:,2:Event));
    %vz = squeeze(B(9,:,2:Event));
    %tlarge = squeeze(B(10,:,2:Event)); 
    %t = squeeze(tlarge(:,1));
    
    
    

    phi_down=angle(sy2+1i*sz2);
    delta_phi=phi_up-phi_down;
    delta_phisave=delta_phi;
     if any(any(delta_phi>6))
        delta_phi(delta_phi>6)=-(2*pi-abs(delta_phi(delta_phi>6)));
    end
    if any(any(delta_phi<-6))
        delta_phi(delta_phi<-6)=(2*pi-abs(delta_phi(delta_phi<-6)));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%     
%     datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(start+(n-1)*4+2),'.dat');
%     %E up without the gradient
%     fileID = fopen(datafile);
%     A = fread(fileID, 'double');
%     fclose(fileID);
% 
%     B = reshape(A, 10, Bin, Event);
% % the following are matrices (time (Bin+1), particles (Event-1))
%     sx3 = squeeze(B(1,:,2:Event));
%     sy3 = squeeze(B(2,:,2:Event));
%     sz3 = squeeze(B(3,:,2:Event));
%     %x =  squeeze(B(4,:,2:Event));
%     %y =  squeeze(B(5,:,2:Event));
%     %z = squeeze(B(6,:,2:Event));
%     %vx =  squeeze(B(7,:,2:Event));
%     %vy =  squeeze(B(8,:,2:Event));
%     %vz = squeeze(B(9,:,2:Event));
%    % tlarge = 1.6E-4squeeze(B(10,:,2:Event)); 
%    % t = squeeze(tlarge(:,1));
% 
%     phi_up0=angle(sy3+1i*sz3);
%  
%     datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(start+(n-1)*4+3),'.dat');
%     %E down without the gradient
%     fileID = fopen(datafile);
%     A = fread(fileID, 'double');
%     fclose(fileID);
% 
%     B = reshape(A, 10, Bin, Event);
% % the following are matrices (time (Bin+1), particles (Event-1))
%     sx4 = squeeze(B(1,:,2:Event));
%     sy4 = squeeze(B(2,:,2:Event));
%     sz4 = squeeze(B(3,:,2:Event));
%     x =  squeeze(B(4,:,2:Event));
%     y =  squeeze(B(5,:,2:Event));
%     z = squeeze(B(6,:,2:Event));
%     vx =  squeeze(B(7,:,2:Event));
%     vy =  squeeze(B(8,:,2:Event));
%     vz = squeeze(B(9,:,2:Event));
%     tlarge = squeeze(B(10,:,2:Event)); 
    t = squeeze(tlarge(:,1));

    B0hat=[1,0,0];
    sphi1=ExtractPhases(sx1,sy1,sz1,sx2,sy2,sz2,B0hat);
    %sphi2=ExtractPhases(sx3,sy3,sz3,sx4,sy4,sz4,B0hat);
    
  figure

    
    %crappy way, points arn't weighted to the time. 
    %dfsim(n)=mean(sphi1(2:end)./t(2:end))/2/pi;
    

    [theta] = getTheta(sx1,sy1,sz1,sx2,sy2,sz2);
   %why are these two sooooo different!!!
    thetamean=mean(theta,2);
    meantheta=getmeanTheta(sx1,sy1,sz1,sx2,sy2,sz2);
   
   
   if n <5
       
       phiupper=movmax(sphi1,2);
       philower=movmin(sphi1,10);
   end
   if n == 5
       phiupper=movmax(sphi1,2);
       philower=movmin(sphi1,8);
   
   end
   if n >5 
   philower= movmin(sphi1,6);
   phiupper=movmax(sphi1,2);
   end
   if n == 6 
        philower=movmin(sphi1,6);
        phiupper=movmax(sphi1,2);
   end
       
   if n ==9
     phiupper=movmax(sphi1,2);
     philower=movmin(sphi1,10);
   end
   
   
  [cfun0,gof0]=fitGeophase(sphi1(1:round(length(t)/5)),t(1:round(length(t)/5)));%fitGeophase(abs(sphi1./S),t);%fitGeophase(sphi1(1:round(length(t))),t(1:round(length(t))));
    hold on
    plot(t(1:round(length(t)/5)),sphi1(1:round(length(t)/5)));
 
    
   
    [cfun1,gof1]=fitGeophase(phiupper(1:round(length(t)/5)),t(1:round(length(t)/5)));
    [cfun2,gof2]=fitGeophase(philower(1:round(length(t)/5)),t(1:round(length(t)/5)));
     if abs(cfun1.a)>abs(cfun2.a)
        cfun=cfun1;
        gof=gof1;
    else
        cfun=cfun2;
        gof=gof2;
     end
    
     %plot(t(1:round(length(t)/5)),philower(1:round(length(t)/5)));
    
     %for 250 mK we must use the envelope. (take out RMS.)
    
     %this is an effort to get the cfun curve to lie on the obvious
     %symmetry of the sphi1 curve. This is clearly what the theory is predicting
     % and seeig this is somthing that is trivial 
     %for humans but a computer cant seem to find it. 
     %I could fit these more realiably by hand. 
     if n == 3 ||  n == 9
     
     dfsim(n)=mean((cfun.a))/2/pi;
     dfsimdev(n)=diff(confint(cfun,0.2)/2);
    
    elseif n == 4 || n == 5
    dfsim(n)=(cfun0.a)/2/pi;
    dfsimdev(n)=diff(confint(cfun0,0.2)/2);
    cfun.a=cfun0.a;
    elseif n==1
    dfsim(n)=2*(cfun0.a)/2/pi;
    dfsimdev(n)=diff(confint(cfun0,0.2)/2);
    cfun.a=cfun0.a;
    else
    dfsim(n)=(cfun.a+2*cfun0.a)/3/2/pi;
    dfsimdev(n)=diff(confint(cfun,0.2)/2);
    cfun.a=mean(cfun.a+cfun0.a);
    end
    plot(cfun);
%     
     
     %For T=300 mK we don't see the RMS for some reason (diffusion is too
    %slow probably, so we can use the mean, this is probably an interesting effect. )
%     dfsim(n)=cfun.a/2/pi;
%     dfsimdev(n)=diff(confint(cfun,0.2)/2);
%      plot(cfun);
    
    
    %delphi1{n}=sphi1;
%    delphi2{n}=sphi2;
    
    
   

%        
%     
%     
%     phi_down0=angle(sy4+1i*sz4);
%     delta_phi0=phi_up0-phi_down0;
%     
%      if any(any(delta_phi0>6))
%         delta_phi0(delta_phi0>6)=-(2*pi-abs(delta_phi0(delta_phi0>6)));
%     end
%     if any(any(delta_phi0<-6))
%         delta_phi0(delta_phi0<-6)=(2*pi-abs(delta_phi0(delta_phi0<-6)));
%     end
% %     
%     delta_phi=delta_phi-delta_phi0;
%     dphi{n}=delta_phi;
%     delta_phi=sum(delta_phi,2)/(Event-1);
%     delphi{n}=delta_phi;
%     
    
%     figure
%     fig=scatter(t,smooth(delta_phi,3));
%     axis([0 20 -0.0001 0.0001])
%     title(num2str(w(n)))
%     ylabel('relative \delta_\phi')
%     xlabel('t')
    %saveas(fig,[pwd strcat('/geophase2/w=',num2str(w(n)),'.jpg')]);
%for n=1:8   
[dfmsdp,dfmp]=ThermalCollisionShifts(w(n), B1(n), 2.037947093e8, B0, temperature,G,L); 

%modulation correction. 
norder=51;
Fr=1+1/4/sqrt(norder)/sqrt(pi)*(1/a(n)-a(n))*(erf(3/2*sqrt(norder)*pi)+erf(1/2*sqrt(norder)*pi));

dfmsd(n)=Fr*dfmsdp(104);

end

figure(12)
hold on
w=[3000 6000 1800 2100 2400 4200 10000 18000 30000];
 %all the frequencies
scatter(w(1:length(dfmsd))/2/pi,dfmsd)
errorbar(w(1:length(dfsim))/2/pi,dfsim,dfsimdev,'LineStyle','none','Marker','.')
ylabel('\deltaf (Hz)')
xlabel('f_{rf} (Hz)')
%lim=axis;
%lim(2)=0;
%axis(lim);
%   for n=1:8   
%    dfmsdp=ThermalCollisionShifts(w(n), B1(n), 2.037947093e8, 3e-6, temperature); 
%    dfmsd(n)=dfmsdp(157);
%   end
