

%make sure these agree
B0 = 30e-7;         
w_rf = 100000;
B_rf = 6.488345e-4;

%[~,peak_n] = max(abs(fft(pangle_n,nBins*100)));
%w_eff = (peak_n *2*pi/(runTime*100))/2;

w_n_t = -gamma_3*besselj(0,gamma_3*B_rf/w_rf)*B0(1,1);   % theoretical w_eff value for n
w_He_t = -gamma_n*besselj(0,gamma_n*B_rf/w_rf)*B0(1,1);   % theoretical w_eff value for He
w_mean = (w_n_t + w_He_t)/2;


%deltat = 0.1/(2*fm);
dphi = 0.8;             % ang sep of He and n


deltat1 = 0.05;
deltat2 = 0.05;
%fm = 1/(deltat*20);


%graphical representation
w_ = 0:100:(2*w_rf);
f = (gamma_n * besselj(0,gamma_n*B_rf./w_) - gamma_3 * besselj(0,gamma_3*B_rf./w_) )*deltat1*B0;
plot(w_,f)
hold on
plot([w_rf w_rf], [-1 1])
plot([0 w_rf*2],[2*dphi 2*dphi])
plot([0 w_rf*2],[-2*dphi -2*dphi])
figure;

%actually solving it
syms w
f_l = (gamma_n * besselj(0,gamma_n*B_rf/w) - gamma_3 * besselj(0,gamma_3*B_rf/w) )*B0*deltat1 - 2*dphi;
f_u = (gamma_n * besselj(0,gamma_n*B_rf/w) - gamma_3 * besselj(0,gamma_3*B_rf/w) )*B0*deltat2 + 2*dphi;

w1 = vpasolve(f_l,w, [w_rf/2 w_rf]);
dw1 = w1 - w_rf;
w2 = vpasolve(f_u,w, [w_rf 2*w_rf]);
dw2 = w2 - w_rf;


frac_mod = 0.1; % fraction of time spent modulating
fm = 1/(floor(0.5*((deltat1+deltat2)/frac_mod - deltat1 - deltat2)/(2*pi/w_mean))*2*(2*pi/w_mean)+deltat1+deltat2);

deltat1_m = 0.05:0.00005:0.053;

par_m = zeros(length(deltat1_m),5);
i = 1;
for deltat = deltat1_m
    f_1 = (gamma_n*besselj(0,gamma_n*B_rf/w) - gamma_3*besselj(0,gamma_3*B_rf/w) )*B0*deltat + 2*dphi;
    f_2 = (gamma_n*besselj(0,gamma_n*B_rf/w) - gamma_3*besselj(0,gamma_3*B_rf/w) )*B0*deltat - 2*dphi;
    w1 = vpasolve(f_1,w, [w_rf/3 w_rf]);
    w2 = vpasolve(f_2,w, [w_rf w_rf*3]);
    % If can't find a solution for a w1 or w2 will set it and phase to -1
    if isempty(w1) == 1
        w1 = -1;
        phase1 = -1;
    else
        phase1 = mod(-gamma_n*besselj(0,gamma_n*B_rf/w1)*B0*deltat,2*pi);
    end
    if isempty(w2) == 1
        w2 = -1;
        phase2 = -1;
    else
        phase2 = mod(-gamma_n*besselj(0,gamma_n*B_rf/w2)*B0*deltat,2*pi);
    end
    par_m(i,:) = [deltat w1 phase1 w2 phase2];
    i = i + 1;
end


plot(par_m(:,1),par_m(:,3));
hold on
plot(par_m(:,1),par_m(:,5));
legend('lower w','higher w')
plot([0.03 0.07], [2*pi-dphi 2*pi-dphi])
plot([0.03 0.07], [dphi dphi])

%deltat1 = 0.05247, 1.473251511e5, 
%deltat2 = 0.05065, 7.6377659e4

