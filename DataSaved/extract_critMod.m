%figure;
hold on;

for i=13
Bin = 20001;
Event = 1;
 runnum=i;

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
    %t = squeeze(tlarge(:,1));
   
    eval(['Sangle',num2str(runnum),'=acos(sx1.*mean(sx1,2)+sy1.*mean(sy1,2)+sz1.*mean(sz1,2));']);
    %eval(['plot(t,rms(Sangle0,2)-rms(Sangle',num2str(runnum),',2));']);
    %drawnow
    %histogram(Sangle(10,:),30)
end