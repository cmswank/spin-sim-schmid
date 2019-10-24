%reset correlations
zz=z;
xx=x;
yy=y;
Ryyg=0;
Rxxg=0;
Rzzg=0;
Rvvx=0;
Rvvy=0;
Rvvz=0;

%%% Part 3 of 3. Calculate Correlation Functions
wtbar=waitbar(0,'I`ll finish when red gets to the end. 3 of 3.');
%yyp=yy;
mmm=size(xx,2);
nnn=size(xx,1);
dt=t(2);

for i = 1:mmm
    %yyp(:,mmm)=yy(:,i)-mean(yy(:,i));
   % yvpp=xcorr(yyp(:,i),vvyy(:,i),'biased');
   % zvpp=xcorr(zz(:,i),vvzz(:,i),'biased');
   %vypp=dt*xcorr(vvyy(:,i),yy(:,i));
   if ~isnan(sum(vvxx(:,i)))
   ypp=xcorr(yy(:,i),'biased');
     zpp=xcorr(zz(:,i),'biased');
    xpp=xcorr(xx(:,i),'biased');
    %vxpp=xcorr(vvxx(:,i),'biased');
    %vypp=xcorr(vvyy(:,i),'biased');
    %vzpp=xcorr(vvzz(:,i),'biased');
   
    % Rvyg=Rvyg+2*yvpp(nnn:2*nnn-1)/mmm;%+vypp(nnn:2*nnn-1)/mmm;
    Ryyg=Ryyg+ypp(nnn:2*nnn-1)/mmm;
  %  Rvzg=Rvzg+2*zvpp(nnn:2*nnn-1)/mmm;
    Rxxg=Rxxg+xpp(nnn:2*nnn-1)/mmm;
    Rzzg=Rzzg+zpp(nnn:2*nnn-1)/mmm;
    %Rvvy=Rvvy+vypp(nnn:2*nnn-1)/mmm;
    %Rvvz=Rvvz+vzpp(nnn:2*nnn-1)/mmm;
    %Rvvx=Rvvx+vxpp(nnn:2*nnn-1)/mmm;
    waitbar(i/mmm);
   end 
end
close(wtbar);
 Sxxg=2*fft(Rxxg)/length(Rxxg); %normalized position spectrum (cm^2)
 Svxg=-fft(diff(Rxxg)); %normalized velocity position (cm^2/s) 

 Syyg=2*fft(Ryyg)/length(Ryyg);
 Svyg=-fft(diff(Ryyg));   %(diff needs normalization 1/dt opposite as fft, so they cancel)
 Szzg=2*fft(Rzzg)/length(Rzzg);
 Svzg=-fft(diff(Rzzg));
