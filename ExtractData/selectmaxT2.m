function [sphip,tp,cfun,Speak] = selectmaxT2(t,sx1,sy1,sz1,sphi1)
%CFUNT2 Summary of this function goes here
%   Detailed explanation goes here
S73= abs(mean(sy1+1i*sz1,2));
theta=acos(mean(sx1));
sphip=sphi1;
 tp=t;
 [S73peakppp,bbb]=findpeaks(S73);
 tp=tp(bbb);
 sphip=sphip(bbb);
 
 [S73peakpp,bb]=findpeaks(S73peakppp);
 tp=tp(bb);
 sphip=sphip(bb);
% [S73peakp,b]=findpeaks(S73peakpp);
 Speak=S73peakpp;
 %sphip=sphip(b);
 %tp=tp(b);


 
 cfun=fitGeophase(sphip,tp);
end

