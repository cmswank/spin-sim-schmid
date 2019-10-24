



[T2,fwsd,T2McgregorRf,lam,I2]= ThermalCollisionRelaxJ0(3000,1.9102418e-5,2.037947093e8,3E-6,0.45,G1(1),0.4);



% dphi=zeros(size(sphi1));
% for i = 1:length(sphi1)
%    dphi(i)=mean(diff(sphi1(1:i)./t(2)));
% end
% 
% % 
% % S73=abs(Ssub)/999;
% % 
% %  tp=t;
% %  [S73peakppp,bbb]=findpeaks(S73);
% %  tp=tp(bbb);
% %  [S73peakpp,bb]=findpeaks(S73peakppp);
% %  tp=tp(bb);
% %  [S73peakp,b]=findpeaks(S73peakpp);
% %  Speak=S73peakp;
 tp=tp(b);
 
 %cfun=fitT2(tp,Speak);
 %plot(cfun);
 
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