trials=1000;
wrf=6283.18530718;
 t=0:2E-5:1;
 tsh=0:2E-4:1;
 [A,B,C,D] = butter(10,[40/2500 120/2500],'stop');
 %[A,B,C,D] = butter(10,[980/2500 1020/2500],'bandpass');
 sos=ss2sos(A,B,C,D);
 
 for i = 1:trials 
 noise=random('normal',0,1E-3,5001,1);
 fnoise=sosfilt(sos,noise);
 
 fsnoise=interp1(tsh,fnoise,t);
 B1=cos(wrf*t-pi/2)+fsnoise;
 fileID = fopen(['/data1/cmswank/spin-sim-xliu/BField/B1Pulse',num2str(i),'.dat'],'w');
fwrite(fileID,B1,'double');
fclose(fileID);
end
%GOOD LUCK!!!







