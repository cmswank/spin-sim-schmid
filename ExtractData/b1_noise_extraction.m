% dphi=zeros(size(sphi1));
% for i = 1:length(sphi1)
%    dphi(i)=mean(diff(sphi1(1:i)./t(2)));
% end
% 


k=100; %number of seeds
B0=3E-6;
%%%old gamma=2.037947093e8;
gamma=2.037894585E8;
w0 =gamma*B0;
T=.450;
L=0.4;
w=[3000 6000 1800 2100 10000 18000 30000 2400 4200]; %all the frequencies
B1 = [1.9102418e-5, 3.8750592e-5, 1.10631580e-5, 1.31017478e-5, 6.477662320e-5, 1.167322440e-4, 1.946176552e-4, 1.51160267e-5, 2.69935194e-5];
G0=3E-3*B0;
%G1=G0*besselj(0,gamma.*B1./w)./besselj(1,gamma.*B1./w).*w./w0;
G1=G0*1.509.*w./w0;
start=0;
pfn = '/home/xliu5/spin-sim-SD/parSpinDressing'; 	%parameter file name
clear delphi1 delphi2 delphi dphi;
delphi1{1}=0;
delphi2{1}=0;
dfsim=zeros(k,1);
dfsimdev=zeros(k,1);
Sf1=zeros(k,1);
Sf2=zeros(k,1);
Sf3=zeros(k,1);
T21=zeros(k,1);
T22=zeros(k,1);
T23=zeros(k,1);
Bin = 5001;
Event = 1;

sxh3=zeros(Bin,k);
syh3=zeros(Bin,k);
szh3=zeros(Bin,k);
sxn=zeros(Bin,k);
syn=zeros(Bin,k);
szn=zeros(Bin,k);
xh3=zeros(Bin,k);
yh3=zeros(Bin,k);
zh3=zeros(Bin,k);
xn=zeros(Bin,k);
yn=zeros(Bin,k);
zn=zeros(Bin,k);



for n=1:k
    datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(start+n-1),'.dat');
    %spin up with the gradient
    

% import data
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);

    B = reshape(A, 10, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
    sx1 = squeeze(B(1,:,1:Event));
    sy1 = squeeze(B(2,:,1:Event));
    sz1 = squeeze(B(3,:,1:Event));
    x1 =  squeeze(B(4,:,1:Event));
    y1 =  squeeze(B(5,:,1:Event));
    z1 = squeeze(B(6,:,1:Event));
    %vx =  squeeze(B(7,:,1:Event));
    %vy =  squeeze(B(8,:,1:Event));
    %vz = squeeze(B(9,:,1:Event));
    tlarge = squeeze(B(10,:,1:Event)); 
    %t = squeeze(tlarge(:,1));
    t=tlarge;
    
    sxh3(:,n)=sx1';
    syh3(:,n)=sy1';
    szh3(:,n)=sz1';
    xh3(:,n)=x1';
    yh3(:,n)=y1';
    zh3(:,n)=z1';
    %phi_up=angle(sy1+1i*sz1);
    
   % disp('hello');
    datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(start+n-1+k),'.dat');
    %spin down with the gradient
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);

    B = reshape(A, 10, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
    sx2 = squeeze(B(1,:,1:Event));
    sy2 = squeeze(B(2,:,1:Event));
    sz2 = squeeze(B(3,:,1:Event));
   x =  squeeze(B(4,:,1:Event));
   y =  squeeze(B(5,:,1:Event));
   z = squeeze(B(6,:,1:Event));
    %vx =  squeeze(B(7,:,2:Event));
    %vy =  squeeze(B(8,:,2:Event));
    %vz = squeeze(B(9,:,2:Event));
    %tlarge = squeeze(B(10,:,2:Event)); 
    %t = squeeze(tlarge(:,1));
    sxn(:,n)=sx2';
    syn(:,n)=sy2';
    szn(:,n)=sz2';
    xn(:,n)=x';
    yn(:,n)=y';
    zn(:,n)=z';
    
    
end
