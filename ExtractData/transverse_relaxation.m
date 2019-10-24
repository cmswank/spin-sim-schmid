% dphi=zeros(size(sphi1));
% for i = 1:length(sphi1)
%    dphi(i)=mean(diff(sphi1(1:i)./t(2)));
% end
% 


k=1; %number of frequencies
B0=3E-6;
gamma=2.037947093e8;
w0 =gamma*B0;
T=.450;
L=0.4;
w=[3000 6000 1800 2100 10000 18000 30000 2400 4200]; %all the frequencies
B1 = [1.9102418e-5, 3.8750592e-5, 1.10631580e-5, 1.31017478e-5, 6.477662320e-5, 1.167322440e-4, 1.946176552e-4, 1.51160267e-5, 2.69935194e-5];
G0=3E-3*B0;
%G1=G0*besselj(0,gamma.*B1./w)./besselj(1,gamma.*B1./w).*w./w0;
G1=G0*1.509.*w./w0;
start=10;
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
for n=1:k
    datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(start),'.dat');
    %spin up with the gradient
    Bin = 5001;
    Event = 1;
    

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
    vx =  squeeze(B(7,:,1:Event));
    vy =  squeeze(B(8,:,1:Event));
    vz = squeeze(B(9,:,1:Event));
    tlarge = squeeze(B(10,:,1:Event)); 
    %t = squeeze(tlarge(:,1));
    t=tlarge;

    phi_up=angle(sy1+1i*sz1);
    
    disp('hello');
    datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(start+(n-1)*3+1),'.dat');
    %spin down with the gradient
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
    
    datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(start+(n-1)*3+2),'.dat');
    %E up without the gradient
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);

    B = reshape(A, 10, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
    sx3 = squeeze(B(1,:,2:Event));
    sy3 = squeeze(B(2,:,2:Event));
    sz3 = squeeze(B(3,:,2:Event));
    %x =  squeeze(B(4,:,2:Event));
    %y =  squeeze(B(5,:,2:Event));
    %z = squeeze(B(6,:,2:Event));
    %vx =  squeeze(B(7,:,2:Event));
    %vy =  squeeze(B(8,:,2:Event));
    %vz = squeeze(B(9,:,2:Event));
   % tlarge = 1.6E-4squeeze(B(10,:,2:Event)); 
   % t = squeeze(tlarge(:,1));

    %phi_up0=angle(sy3+1i*sz3);
 

%    figure(12)
% hold on
% scatter(w(1:length(Sf1))/2/pi,Sf1);
% scatter(w(1:length(Sf2))/2/pi,Sf2);
% scatter(w(1:length(Sf3))/2/pi,Sf3);
% scatter(w(1:length(T21))/2/pi,T21);
% scatter(w(1:length(T22))/2/pi,T22);
% scatter(w(1:length(T23))/2/pi,T23); t = squeeze(tlarge(:,1));



S1=abs(sum(sy1+1i*sz1,2))/(Event-1);
S2=abs(sum(sy2+1i*sz2,2))/(Event-1);
S3=abs(sum(sy3+1i*sz3,2))/(Event-1);

% Sfun1=cfunT2(t,S1,90);
Sfun2=cfunT2(t,S2,30);
Sfun3=cfunT2(t,S3,30);


if n == 7    
%Sfun1=cfunT2(t,S1,180);
Sfun2=cfunT2(t,S2,60);
Sfun3=cfunT2(t,S3,60);
end
%Newest spin dressed T2 fitter. should be the best!
if n ==6 
 [Ts21a,T21upper,T21lower] = sparseT2fit(t,S1,Event,1E4,1E9,0.96);
 [Ts21b,T21upper,T21lower] = sparseT2fit(t,S1,Event,1E4,1E9,0.94);
 Ts21=1/2*(Ts21a+Ts21b);
elseif n==7
 [Ts21a,T21upper,T21lower] = sparseT2fit(t,S1,Event,1E4,1E9,0.55);
 [Ts21b,T21upper,T21lower] = sparseT2fit(t,S1,Event,1E4,1E9,0.6);
 Ts21=1/2*(Ts21a+Ts21b);
else
    [Ts21,T21upper,T21lower] = sparseT2fit(t,S1,Event,1E4,1E9,0.5);   
end
disp([Ts21,T21upper,T21lower]);
[Ts22,T22upper,T22lower] = sparseT2fit(t,S2,Event,1,1E3,0.5);
[Ts23,T23upper,T23lower] = sparseT2fit(t,S3,Event,1,1E3,0.5);
%figure(13)
%hold on 
Sf1(n)=Ts21;%Sfun1.a;
Sf2(n)=Sfun2.a;
Sf3(n)=Sfun3.a;
DSf1(n)=(T21upper-T21lower)/2;
DSf2(n)=(T22upper-T22lower)/2;
DSf3(n)=(T23upper-T23lower)/2;



figure
hold on
plot(t,S1);
plot(t,S2);
plot(t,S3);
plot(t,exp(-t/Ts21),'r');
plot(t,exp(-t/Ts22),'r');
plot(t,exp(-t/Ts23),'r');
 
 [T21p,fwsdp,T21mcp,lam] = ThermalCollisionRelaxSDAxis(w(n),B1(n),gamma,B0,T,G1(n),L);
 T21(n)=T21p(157);
 
end
%this function is in cgs... doh 


%T21=SpinDressingT2Calc(B1*10^4,w,G1*100,G0*100,T,L*100);
T22=SpinDressingT2Calc(B1*10^4,w,-G1*100,G0*100,T,L*100);
T23=SpinDressingT2Calc(B1*10^4,w,0,G0*100,T,L*100);
figure(37)
hold on
errorbar(w(1:length(Sf1))/2/pi,Sf1,Sf1/sqrt(Event),'LineStyle','none','Marker','.');
errorbar(w(1:length(Sf2))/2/pi,Sf2,Sf2/sqrt(Event),'LineStyle','none','Marker','.');
errorbar(w(1:length(Sf3))/2/pi,Sf3,Sf3/sqrt(Event),'LineStyle','none','Marker','.');
scatter(w(1:length(T21))/2/pi,T21);
scatter(w(1:length(T22))/2/pi,T22);
scatter(w(1:length(T23))/2/pi,T23);

 
% 
% figure
% plot(tpeak,S73peak)
% hold on
% plot(t,S73)

%    sdelta_phi=smooth(delta_phi,80);
%     figure
%     fig=scatter(t,sdelta_phi);
%     axis([0 20 -0.00001 0.00001])
%     title(num2str(w(n)))
%     ylabel('relative \delta_\phi')
%     xlabel('t')