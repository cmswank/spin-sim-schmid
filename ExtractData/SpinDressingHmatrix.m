function [gammap]=SpinDressingHmatrix(B0,Bdress,wdress,gamma,E0,dn,dwf)
%This functions calculates the effective gyromagnetic ratio in an strong RF
%field with angular frequency wdress,
%with a uniform B0 and E field as a perturbation. wdress is an angular frequency.  
%B0 and Bdress are in Gauss, gamma is the gyromagnetic ratio,
%dn is EDM in ecm, E0 is an E field in V/cm. 


Bd=Bdress;
wd=wdress;
hbar=6.582122E-16; %hbar in eVs
%gamma=20378.9;     %gryomagnetic ratio rad/G
gammaE=2*dn/hbar;  %gryoelectric ratio rad/(V/cm)  
%wd=6000;
%B0=0.03;
%E0=-75E3;  %kV/cm
%Bd=.37505920;        %Best number for 6000 according to simulations!
x=gamma*Bd/wd;
y=dwf/wd+gammaE*E0/wd+gamma*B0/wd;


Hsdd1=20+y/2:-1:1+y/2;
Hsdd2=20-y/2:-1:1-y/2;

Hsdd=sort(-[Hsdd1,Hsdd2]);
Hsdd=-Hsdd;
Hsdp=diag(Hsdd);
Hsd=Hsdp;

for i = 1:size(Hsdp,1)
    for j = 1:size(Hsdp,2)
    
        if mod(i+j-5,4)==0 && abs(i-j)<=4
            Hsd(i,j)=x/4;
        end
    
    
    end
end

[~,HSD]=eig(Hsd);

delE=diff(diag(HSD));
%plot(delE);
gammap=delE(size(HSD,1)/2-1)*wd/B0;
%factordown=gammadown/gamma;
% disp(['y = ',num2str(y),', x = ',num2str(x)]);
% disp(['10th order perturbation dressing factor = ',num2str(factordown),', J0(x) = ',num2str(besselj(0,x))]);
% disp(['This is a ', num2str(100*(1-factordown/besselj(0,x))),'% difference']);
% 
% 
% %%%This is when the electric field included in the Hamiltonian
% hbar=6.582122E-16;
% gamma=20378.9*1E4;
% gammaE=2*1E-26/hbar;  %            E metters. we are ignoring hbar it looks like. 
% wd=6000;
% B0=0.03/1E4;
% E0=75E3;
% Bd=.37505920/1E4;        %Best number for 6000 according to simulations!
% x=gamma*Bd/wd;
% y=gammaE*E0/wd+gamma*B0/wd;
% 
% 
% Hsdd1=20+y/2:-1:1+y/2;
% Hsdd2=20-y/2:-1:1-y/2;
% 
% Hsdd=sort(-[Hsdd1,Hsdd2]);
% Hsdd=-Hsdd;
% Hsdp=diag(Hsdd);
% Hsd=Hsdp;
% for i = 1:size(Hsdp,1)
%     for j = 1:size(Hsdp,2)
%     
%         if mod(i+j-5,4)==0 && abs(i-j)<=4
%             Hsd(i,j)=x/4;
%         end
%     
%     
%     end
% end
% 
% [Evec,HSD]=eig(Hsd);In an assignment  A(:) = B, the number of elements in A and B must be the same.

%Error in geometric_phase (line 63)

% 
% delE2=diff(diag(HSD));
% %plot(delE2);
% gammaup=delE2(size(HSD,1)/2-1)*wd/B0;
% factorup=gammaup/gamma;
% % disp(['y = ',num2str(y),', x = ',num2str(x)]);
% % disp(['10th order perturbation dressing factor = ',num2str(factorup),', J0(x) = ',num2str(besselj(0,x))]);
% % disp(['This is a ', num2str(100*(1-factorup/besselj(0,x))),'% difference']);
% 
% 
% 
% 
% hbar=6.582122E-16; %in eVs
% gamma=20378.9*1E4;
% gammaE=2*1E-26/hbar;  %gyroelectric ratio... (for field in V/cm...)    assuming 1E-26 e cm        
% wd=6000;
% B0=0.03/1E4;
% E0=0;%-75E3; %75 kV/cm...
% Bd=.37505920/1E4;        %Best number for 6000 according to simulations!
% x=gamma*Bd/wd;
% y=gammaE*E0/wd+gamma*B0/wd;
% 
% 
% Hsdd1=20+y/2:-1:1+y/2;
% Hsdd2=20-y/2:-1:1-y/2;
% 
% Hsdd=sort(-[Hsdd1,Hsdd2]);
% Hsdd=-Hsdd;
% Hsdp=diag(Hsdd);
% Hsd=Hsdp;
% for i = 1:size(Hsdp,1)
%     for j = 1:size(Hsdp,2)
%     
%         if mod(i+j-5,4)==0 && abs(i-j)<=4
%             Hsd(i,j)=x/4;
%         end
%     
%     
%     end
% end
% 
% [Evec,HSD]=eig(Hsd);
% 
% delE3=diff(diag(HSD));
% %plot(delE3);
% gamma0=delE3(size(HSD,1)/2-1)*wd/B0;
% factor0=gamma0/gamma;
% 
% E0=75E3;
% 
% 
% disp(['Delta S Eup vs Edown ', num2str((gammaup - gammadown)*B0),' Number based off of J0 equaton ', num2str(2*gammaE*besselj(0,x)*E0),', Traditional QM correction', num2str(2*gammaE*factor0*E0)]);
% % %wd
% % 100000
% % 30000
% % 18000
% % 10000
% % 6000
% % 4200
% % 3000
% % 2400
% % 2100
% % 1800 
% 
% % %B best simulated
% % 6.488345280
% % 1.946176552
% % 1.167322440
% % 0.647766232
% % 0.387505920
% % 0.269935194
% % 0.191024180
% % 0.151160267
% % 0.131017478
% % 0.110631580 
% % %Percent error
% % -0.003
% % -0.020
% % -0.053
% % -0.168
% % -0.466
% % -0.957
% % -1.901
% % -3.012
% % -4.001
% % -5.570
% % 