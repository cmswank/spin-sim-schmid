
Bin = 1001;
Event = 1000;
%%%%%%%neutron runs are 0,1,4,5

%neutron run
for runnum=[0,2,4,6,8]
%runnum=1;
%helium run
runnumhe=runnum+1;

datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(runnum),'.dat');

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
  
    datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(runnumhe),'.dat');
  
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);

    B = reshape(A, 10, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
    sx2 = squeeze(B(1,:,2:Event));
    sy2 = squeeze(B(2,:,2:Event));
    sz2 = squeeze(B(3,:,2:Event));
    x =  squeeze(B(4,:,2:Event));
    y =  squeeze(B(5,:,2:Event));
    z = squeeze(B(6,:,2:Event));
    vx =  squeeze(B(7,:,2:Event));
    vy =  squeeze(B(8,:,2:Event));
    vz = squeeze(B(9,:,2:Event));
    tlarge = squeeze(B(10,:,2:Event)); 
    t = squeeze(tlarge(:,1));
  
    
    
    sndots3=sndots3long(sx1,sy1,sz1,sx2,sy2,sz2);
    %sndots3=sndots3Mean(sx1,sy1,sz1,sx2,sy2,sz2);
    %sndots3Eu=sndots3Mean(sx3,sy3,sz3,sx4,sy4,sz4);
    %signalEu=1-sndots3Eu;
    %signalEd=1-sndots3Ed;
%     Sx1=mean(sx1,2);
%     Sx2=mean(sx2,2);
%     Sy1=mean(sy1,2);
%     Sy2=mean(sy2,2);
%     Sz1=mean(sz1,2);
%     Sz2=mean(sz2,2);
%     S1=sqrt(Sx1.^2+Sy1.^2+Sz1.^2);
%     S2=sqrt(Sx2.^2+Sy2.^2+Sz2.^2);
%     costheta=sndots3./(S1.*S2);
    %signalEd=1-exp(t.*(1/T2+1/T2n)).*sndots3Ed;
    eval(['Signal',num2str(runnum),'=1-sndots3;']);
    
end 
  
   