

Bin = 10001;
Event = 10;
 startnum=3015;
 runs=1;
 Sx=zeros(runs,Bin);
Sy=zeros(runs,Bin);
Sz=zeros(runs,Bin);

for i = 1:runs
    runnum=startnum+i-1;

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
    t = squeeze(tlarge(:,1));
    Sx(i,:)=mean(sx1,2);
    Sy(i,:)=mean(sy1,2);
    Sz(i,:)=mean(sz1,2);
    %eval(['Sx',num2str(runnum),'=mean(sx1,2);']);
    %eval(['Sy',num2str(runnum),'=mean(sy1,2);']);
    %eval(['Sz',num2str(runnum),'=mean(sz1,2);']);
   disp(i);
    %eval(['Signal',num2str(runnum),'=sqrt(Sx',num2str(runnum),'.^2+Sy',num2str(runnum),'.^2+Sz',num2str(runnum),'.^2);']);
end