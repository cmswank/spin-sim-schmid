trials=1000;
wrf=6283.18530718;
 t=0:2E-5:1;
 tsh=0:2E-4:1;
 
[b,a] = ellip(7,0.1,60,[20*2*tsh(2) 100*2*tsh(2)],'stop');
 
 %[A,B,C,D] = butter(10,[20/2500 100/2500],'stop');
 %[A,B,C,D] = butter(10,[980/2500 1020/2500],'bandpass');
 %sos=ss2sos(A,B,C,D);
 
 %new way of producing noise free signal using fourier decomposition. 
 bndstp=[20,500];

 
 
 for i = 1:trials
    
     %noise std. 
     noise_std=1E-4; %ppB1
     %noise cutoff frequency. 
     noise_fstop=5000; %Hz
     %
    
  noise=random('normal',0,noise_std,5001,1);
 %noisespec=real(fft(noise));
 %noisespec=noisespec.*[ones(bndstp(1)-1,1);zeros(bndstp(2)-bndstp(1),1);ones(length(noisespec)-2*bndstp(2)+1,1);zeros(bndstp(2)-bndstp(1),1);ones(bndstp(1),1)];
 fnoise=noise;
 %fnoise=sosfilt(sos,noise);
 %fnoise=ifft(noisespec);
 
 fsnoise=pchip(tsh,fnoise,t)';
 %fsnoise=fnoise;
 
 nramp=[zeros(500,1);linspace(0,1,500)';ones(length(fsnoise(1001:end)),1)];
 
 
 B1=sin(wrf*t'); % cos(wrf*t'-pi/2); %aka sine 
 
 bstop=[ones(bndstp(1)*1/t(end)-1,1);zeros((bndstp(2)-bndstp(1))/t(end),1);ones(length(B1)-2*bndstp(2)/t(end)+1,1);zeros((bndstp(2)-bndstp(1))/t(end),1);ones(bndstp(1)/t(end),1)];
 
 hpass=[zeros(bndstp(1)*1/t(end)-1,1);zeros((bndstp(2)-bndstp(1))/t(end),1);ones(length(B1)-2*bndstp(2)/t(end)+1,1);zeros((bndstp(2)-bndstp(1))/t(end),1);zeros(bndstp(1)/t(end),1)];
 
 b1stop=[ones(900,1);zeros(200,1);ones(length(B1)-2*1100,1);zeros(200,1);ones(900,1)];
 
 
 %%%%Filter
 ffinal=hpass;%hpass;%.*B1stop
 
 bfinal= hpass;%hpass;
 
 B1spec=fft(B1).*bfinal;%hpass;%bstop;
 
 fnspec=fft(fsnoise).*ffinal;%hpass;%bstop;
  
 fsnoise=real(ifft(fnspec));
 
 %fsnoise=fsnoise.*nramp;
 
 %fsnoise=filter(b,a,fsnoise);
 
 %fsnoise=real(conv(fsnoise,ifft(bstop)));
 
 %fspec=real(fft(fsnoise));
 
 B1=real(ifft(B1spec));
 
 
 B1=B1+fsnoise;
 
 fileID = fopen(['/data1/cmswank/spin-sim-xliu/BField/B1Pulse',num2str(i),'.dat'],'w');
 
 fwrite(fileID,B1,'double');
 
 fclose(fileID);

 end
%GOOD LUCK!!!







