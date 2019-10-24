function [phi,theta]=ExtractPhases(sx1,sy1,sz1,sx2,sy2,sz2,B0hat)
%get the phase difference w.r.t. th B0 direction and/or total angle between 2
%spins (or 2 groups of spins). 
%sx1(t,n) is an array, where t is time index and n is paticle number, respective for all spins. 
%B0hat = [B0x/B0mag,B0y/B0mag,B0z/B0mag]; B0x,B0y,B0z are a scalers and B0mag is the magnitude of B0. 

%relative phase two sets of spins given a B0 direction. 
B0p=B0hat;
B0sq=dot(B0p,B0p);
phi=zeros(size(sx1,1),1);
theta=zeros(size(sx1,1),1);
for i = 1:size(sx1,1)
spin1=[sx1(i,:)',sy1(i,:)',sz1(i,:)'];
spin2=[sx2(i,:)',sy2(i,:)',sz2(i,:)'];
B0=ones(size(spin1,1),1)*B0p;
s1xs2=cross(spin1,spin2);
s1db0=sum(spin1.*B0,2); %s1 dot B0 product
s2db0=sum(spin2.*B0,2); %s2 dot B0 product
B0dSxS=sum(B0.*s1xs2,2); %etc
snphi=B0dSxS./(sqrt(B0sq-s1db0.^2).*sqrt(B0sq-s2db0.^2));
phi(i)=mean(asin(snphi));
theta(i)=mean(acos(sum(spin1.*spin2,2)));
end



