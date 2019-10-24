Rxxd=0;
Ryyd=0;
Rvyd=0;
Rzzd=0;
Rvzd=0;
Rzvy=0;
Rzvz=0;
Rxy=0;
Rzy=0;
Rzz=0;
wtbar=waitbar(0,'I`ll finish when blue gets to the end.');
yy=y;
yyp=yy;
xx=x;
vxx=vx;8e-7
vyy=vy;
zz=z;
vzz=vz;
mmm=size(x,2);
nnn=size(x,1);
for i = 1:mmm
    %yyp(:,mmm)=yy(:,i)-mean(yy(:,i));
   % yvpp=xcorr(yyp(:,i),vvyy(:,i),'biased');
   % zvpp=xcorr(zz(:,i),vvzz(:,i),'biased');
   %vypp=dt*xcorr(vvyy(:,i),yy(:,i));
    ypp=xcorr(yyp(:,i),'biased');
    xypp=xcorr(xx(:,i),yyp(:,i),'biased');
    zypp=xcorr(zz(:,i),yyp(:,i),'biased');
    zpp=xcorr(zz(:,i),'biased');
    xpp=xcorr(xx(:,i),'biased');
    zvypp=xcorr(zz(:,i),vyy(:,i),'biased');
    zvzpp=xcorr(zz(:,i),vzz(:,i),'biased');
    
   % Rvyg=Rvyg+2*yvpp(nnn:2*nnn-1)/mmm;%+vypp(nnn:2*nnn-1)/mmm;
    Ryyd=Ryyd+ypp(nnn:2*nnn-1)/mmm;
  %  Rvzg=Rvzg+2*zvpp(nnn:2*nnn-1)/mmm;
    Rxxd=Rxxd+xpp(nnn:2*nnn-1)/mmm;
    Rxy=Rxy+xypp(nnn:2*nnn-1)/mmm;
    Rzy=Rzy+zypp(nnn:2*nnn-1)/mmm;
    Rzz=Rzz+zpp(nnn:2*nnn-1)/mmm;
    Rzvy=Rzvy+zvypp(nnn:2*nnn-1)/mmm;
    Rzvz=Rzvz+zvzpp(1:2*nnn-1)/mmm;
    waitbar(i/mmm);
    
end
close(wtbar);
 
pad=1E7;
%Szvz=-fft(diff(Rzz),pad);
%Szvy=-fft(diff(Rzy),pad);
 Szvy=2*fft(Rzvy,pad)/pad;
 Szvz=2*fft(Rzvz,pad)/pad;
 f=1/t(2)/pad:1/t(2)/pad:1/t(2);
Sxvx=fft(diff(Rxxd),pad);

%    Syyd=2*fft(Ryyd)/length(Ryyd);
%  Svyd=-fft(diff(Ryyd));
%  
%  Sxx=2*fft(Rxxd,1E6)/1E6;
%  Sxy=2*fft(Rxy)/length(Rxy);
%  Svzd=-fft(diff(Ryyd));
 Sxvy=-fft(diff(Rxy),pad);
% 
