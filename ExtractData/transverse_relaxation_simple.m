% dphi=zeros(size(sphi1));
% for i = 1:length(sphi1)
%    dphi(i)=mean(diff(sphi1(1:i)./t(2)));
% end
% 


k=9; %number of frequencies
B0=3E-6;
gamma=2.037947093e8;
w0=gamma*B0;
T=.450;
L=0.4;
w=[3000 6000 1800 2100 10000 18000 30000 2400 4200]; %all the frequencies
B1 = [1.9102418e-5, 3.8750592e-5, 1.10631580e-5, 1.31017478e-5, 6.477662320e-5, 1.167322440e-4, 1.946176552e-4, 1.51160267e-5, 2.69935194e-5];
G0=1E-1*B0;
G1=G0*besselj(0,gamma.*B1./w)./besselj(1,gamma.*B1./w).*w./w0;





%G1=G0*1.509.*w./w0;
start=2200;
pfn = '/home/xliu5/spin-sim-SD/parSpinDressing'; 	%parameter file name
clear delphi1 delphi2 delphi dphi;
delphi1{1}=0;
delphi2{1}=0;
dfsim=zeros(k,1);
dfsimdev=zeros(k,1);
Sf1=zeros(k,1);
DSf1=zeros(k,1);
Sf2=zeros(k,1);
Sf3=zeros(k,1);
T21=zeros(k,1);
T22=zeros(k,1);
T23=zeros(k,1);
T21mc=zeros(k,1);
Speak=zeros(k,1);
for n=1:k%7
    datafile = strcat('/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_',num2str(start+(n-1)),'.dat');
    %spin up with the gradient
    Bin = 1001;
    Event = 5000;
    

% import data
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);

    B = reshape(A, 10, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
    sx1 = squeeze(B(1,:,2:Event));
    sy1 = squeeze(B(2,:,2:Event));
    sz1 = squeeze(B(3,:,2:Event));
    x1 =  squeeze(B(4,:,2:Event));
    y1 =  squeeze(B(5,:,2:Event));% selected based on what you read.
    z1 = squeeze(B(6,:,2:Event));
    vx =  squeeze(B(7,:,2:Event));
    vy =  squeeze(B(8,:,2:Event));
    vz = squeeze(B(9,:,2:Event));
    tlarge = squeeze(B(10,:,2:Event)); 
    t = squeeze(tlarge(:,1));

    phi_up=angle(sy1+1i*sz1);


S1=abs(sum(sy1+1i*sz1,2))/(Event-1);
%S2=abs(sum(sy2+1i*sz2,2))/(Event-1);
%S3=abs(sum(sy3+1i*sz3,2))/(Event-1);
% id2=1;%round(length(t)/2);
% [Sfun1,Speak]=cfunT2(t(id2:end),S1(id2:end),70);
%Sfun2=cfunT2(t,S2,30);
%Sfun3=cfunT2(t,S3,30);

%Newest spin dressed T2 fitter. should be the best!
[T2,T2upper,T2lower] = sparseT2fit(t,S1,Event,1E1,1E5,0.5);


%if n <= 7 || n >=5   
%id2=round(length(t)/3);

%Sfun1=cfunT2(t(id2:end),S1(id2:end),120);
%Sfun3=cfunT2(t,S3,60);
%end

%figure(13)
%hold on 
Sf1(n)=T2;
%Sf2(n)=Sfun2.a;
%Sf3(n)=Sfun3.a;




 figure
 hold on
 axis([0 max(t) 0.995 1.0001])
 plot(t,S1,'b');
 plot(t,exp(-t/T2),'r');
 
 % % %%%plot(t,S2);
% % %%%plot(t,S3);
% %plot(Sfun1,'g');
% %plot(t,Speak,'r');
%  %%%%plot(Sfun2);
%%% plot(Sfun3);

%propagation of error through fit
DSf1(n)=(T2upper-T2lower)/2;
 
%total-ish error  
 %DSf1(n)=sqrt(DSf1(n).^2+Event);
 
 
%  %this is for along the b0 direction
%  [T21p,fwsdp,T21mcp,lam]= ThermalCollisionRelax(w(n),B1(n),gamma,B0,T,G1(n),L);
%  T21(n)=T21p(157);
%  T21mc(n)=T21mcp;
          %%%%%%%%% B1=B1(1).*ones(size(B1)); %this is for run 2260-2268
          %%%%%%%%% ABSOLUTEUTLY DESTROYS OTHER RUNS!!! DON"T USE!!!!
%  %%%%this is for not along b0 direction
 [T21p,fwsdp,T21mcp,lam] = ThermalCollisionRelaxJ0(w(n),B1(n),gamma,B0,T,G1(n),L);
 T21(n)=T21p(157);
 %disp(T21(n));
 T21mc(n)=T21mcp(1,157);
end



%this function is in cgs... doh 
%T21=SpinDressingT2Calc(B1*10^4,w,G1*100,0,T,L*100);
%T22=SpinDressingT2Calc(B1*10^4,w,-G1*100,G0*100,T,L*100);
%T23=SpinDressingT2Calc(B1*10^4,w,0,G0*100,T,L*100);
figure(54)
hold on
%scatter(w(1:length(Sf1))/2/pi,Sf1);
%scatter(w(1:length(Sf2))/2/pi,Sf2);
%scatter(w(1:length(Sf3))/2/pi,Sf3);
%
%scatter(w(1:length(T22))/2/pi,T22);
%scatter(w(1:length(T23))/2/pi,T23);

% %don't plot high frequencies (long runs not complete)
%  wp=w([1:5]);%,8,9]);
%  T21=T21([1:5]);%,8,9]);
%  T21mc=T21mc([1:5]);%,8,9]);
%  Sf1=Sf1([1:5]);%,8,9]);
%  DSf1=DSf1([1:5]);%,8,9]);


% %don't plot high frequencies
%  wp=w([1:5,8,9]);
%  T21=T21([1:5,8,9]);
%  T21mc=T21mc([1:5,8,9]);
%  Sf1=Sf1([1:5,8,9]);
%  DSf1=DSf1([1:5,8,9]);
%scatter(wp(1:length(T21))/2/pi,T21);
%scatter(wp(1:length(T21mc))/2/pi,T21mc,'Marker','*');
%errorbar(wp(1:length(Sf1))/2/pi,Sf1,Sf1./DSf1,'LineStyle','none','Marker','.');
 



scatter(w(1:length(T21))/2/pi,T21);
%scatter(w(1:length(T21mc))/2/pi,T21mc,'Marker','*');
errorbar(w(1:length(Sf1))/2/pi,Sf1,DSf1,'LineStyle','none','Marker','.');



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