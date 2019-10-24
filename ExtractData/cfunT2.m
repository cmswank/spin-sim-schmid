function [cfun,Speak] = cfunT2(t,S73,n)
%CFUNT2 Summary of this function goes here
%   Detailed explanation goes here
%S73= abs(mean(sy1+1i*sz1,2));
%sphip=sphi1;
 tp=t;
 Speak=movmax(S73,n);
%  [S73peakppp,bbb]=findpeaks(S73);
%  tp=tp(bbb);
%  %sphip=sphip(bbb);
%  
%  [S73peakpp,bb]=findpeaks(S73peakppp);
%  tp=tp(bb);
%  %sphip=sphip(bb);
%  %[S73peakp,b]=findpeaks(S73peakpp);
%  Speak=S73peakpp;
% % tp=tp(b);


 
 cfun=fitT2(tp,Speak);
end

