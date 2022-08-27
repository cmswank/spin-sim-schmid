%run 1 is SD noise bandstop 40-100 1E-3 sigma...
%run 2000 is unknown, can't seem to extract it. 
%run 4000 is SD noise  perfect bandstop 50-70
%run 6000 doesn't exist for some reason, must have deleted it on accident. 
%run 8000 is ellip filter 7th order bp 20-100 (0.1 dB ripple and -60 dB, this made it off critically dressed. but good results. other than the bounce. )
%run 10000 is 0.1 sec for better reslution on bounce 20-100 perfectbp
%run 12000 is ramped noise (0.1 s no noise, 0.1 s ramp,)
%run 14000 is short time ramped noise 0.1 s total, 0.01 s no noise, 0.01 s ramp. 
%run 16000 is short time, no ramp, high pass filter. 
%run 18000 is shor time, highpass, bandstop noise around B1 freq. +- 100 Hz. 
%run 20000 is short time high pass at 800  Hz,
%run 22000 is long time (10 s) perfect bandstop 20-100 Hz.  
%run 24000 is long time (10 s) perfect highpass 600 Hz.  
%run 26000 is short time perfect highpass 100 Hz. 
%run 28000 same as 26000 but with starting along y. (B0 is Bx, and B1 is Bz). 
%run 30000 same as 26 but with starting between z and y, 1/sqrt(2) for each
%run 32000 no filter, 1 second
%run 34000 200 Hz highpass, 1 second but 10001 time steps (10*1001). 
%run 36000 200 Hz highpass, 1800-2200 bandstop 1001 steps. 
%run 38000 Bandpass from 950-1050. 1 s 1001 bins
%run 40000 No Noise with 950-1050 Bandpass on B1.  
%run 42000 bandpass around B1 Noise 10001
%run 44000 ???? 
%run 46000 100-1025 bandpass. 10 seconds 10001 bins
%run 48000 200-1025 bandpass. 10 seconds 10001 bins
%run 50000 40-1025 bandstop. 

%now we changed run 8200 is power supply noise robust modulated dressing.
%
startnum=1082;
Bin = 1001;

Event = 2000;
runs=1;
Sxn=zeros(runs,Bin);
%Sxmean=zeros(Bin,1);
Syn=zeros(runs,Bin);
Szn=zeros(runs,Bin);
Sx3=zeros(runs,Bin);
%Sxmean=zeros(Bin,1);
Sy3=zeros(runs,Bin);
Sz3=zeros(runs,Bin);
%phaseData=zeros(runs,Bin);
for i = 1:runs
  
    runnum=startnum+(i-1);

%datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(runnum),'.dat');

% import data
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);

    B = reshape(A, 10, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
    sx1 = squeeze(B(1,:,1:Event));
    sy1 = squeeze(B(2,:,1:Event));
    sz1 = squeeze(B(3,:,1:Event));
    x =  squeeze(B(4,:,1:Event));
    y =  squeeze(B(5,:,1:Event));
    z = squeeze(B(6,:,1:Event));
    vx =  squeeze(B(7,:,1:Event));
    vy =  squeeze(B(8,:,1:Event));
    vz = squeeze(B(9,:,1:Event));
    tlarge = squeeze(B(10,:,1:Event)); 
    
   if Event>1
       t = squeeze(tlarge(:,1));
    Sx(i,:)=mean(sx1,2);
    Sy(i,:)=mean(sy1,2);
    Sz(i,:)=mean(sz1,2);
   else
       t=tlarge;
       Sxn(i,:)=sx1;
       Syn(i,:)=sy1;
       Szn(i,:)=sz1;
   end
   %%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%% Helium
   %%%%%%%%%%%%%%%%%%
   runnum=startnum+runs+i-1;  %different from above 2*(i-1). 

datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(runnum),'.dat');

% import data
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);

    B = reshape(A, 10, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
    sx1 = squeeze(B(1,:,1:Event));
    sy1 = squeeze(B(2,:,1:Event));
    sz1 = squeeze(B(3,:,1:Event));
    x =  squeeze(B(4,:,1:Event));
    y =  squeeze(B(5,:,1:Event));
    z = squeeze(B(6,:,1:Event));
    vx =  squeeze(B(7,:,1:Event));
    vy =  squeeze(B(8,:,1:Event));
    vz = squeeze(B(9,:,1:Event));
    tlarge = squeeze(B(10,:,1:Event)); 
    
   if Event>1
       t = squeeze(tlarge(:,1));
    Sx3(i,:)=mean(sx1,2);
    Sy3(i,:)=mean(sy1,2);
    Sz3(i,:)=mean(sz1,2);
   else
       t=tlarge;
       Sx3(i,:)=sx1;
       Sy3(i,:)=sy1;
       Sz3(i,:)=sz1;
   end
   
   
   
   
   %phaseData(i,:)=mean(acos(sx1.*mean(sx1,2)+sy1.*mean(sy1,2).*sz1.*mean(sz1,2)),2);
   %phaseNoise=std(acos(sx1.*mean(sx1,2)+sy1.*mean(sy1,2).*sz1.*mean(sz1,2))');
    %eval(['Sx',num2str(runnum),'=mean(sx1,2);']);
    %eval(['Sy',num2str(runnum),'=mean(sy1,2);']);
    %eval(['Sz',num2str(runnum),'=mean(sz1,2);']);
  
   %for ii=1:length(Sx)
    %   Sxmean(ii)=mean(Sx(ii:end));
   %end
    
    
    %eval(['Signal',num2str(runnum),'=sqrt(Sx',num2str(runnum),'.^2+Sy',num2str(runnum),'.^2+Sz',num2str(runnum),'.^2);']);
end
phaseData=acos(Sxn.*Sx3+Syn.*Sy3+Szn.*Sz3);

phaseNoise=std(phaseData);
B0hat=[1 0 0 ];
phiData=ExtractPhasesLong(Sxn,Syn,Szn,Sx3,Sy3,Sz3,B0hat);
phiNoise=std(phiData);
%phaseNoise=std(acos(sqrt(Sx.*mean(Sx)+Sy.*mean(Sy)+Sz.*mean(Sz))));
%plot(t,phaseNoise)
%plot(Sxmean)