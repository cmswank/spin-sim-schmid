
[~,peak_n] = max(abs(fft(s_n(3,:),nBins*100)));
w_n_act = -peak_n *2*pi/(runTime*100)
w_n_theory = gamma_n*besselj(0,gamma_n*B_rf/w_rf)*B0(1,1)

[~,peak_He] = max(abs(fft(s_He(3,:),nBins*100)));
w_He_act = -peak_He *2*pi/(runTime*100)
w_He_theory = gamma_3 * besselj(0,gamma_3*B_rf/w_rf)*B0(1,1)

adj = w_He_act - w_n_act;

%{
% adjusting w_rf
w_rf__ = (w_rf-50):0.01:(w_rf+50);
w_diff = gamma_3*besselj(0,gamma_3*B_rf./w_rf__)*B0(1,1) - ...
         gamma_n*besselj(0,gamma_n*B_rf./w_rf__)*B0(1,1);
plot(w_rf__,w_diff)
hold on
plot([w_rf-50 w_rf+50],[adj adj])
xlabel('w_{rf}')
ylabel('w_{He3} - w_n')
figure;

syms w_rf_
f_d = gamma_3*besselj(0,gamma_3*B_rf/w_rf_)*B0(1,1) - ...
      gamma_n*besselj(0,gamma_n*B_rf/w_rf_)*B0(1,1) + adj;
w_rf_adj = vpasolve(f_d,w_rf_,[w_rf-50 w_rf+50]);

w_He_adj = gamma_3*besselj(0,gamma_3*B_rf/w_rf_adj)*B0(1,1);
w_n_adj = gamma_n*besselj(0,gamma_n*B_rf/w_rf_adj)*B0(1,1);
%}


% adjusting B_rf
B_rf__ = (B_rf-(0.5e-4)):1e-8:(B_rf+(0.5e-4));
w_diff = gamma_3*besselj(0,gamma_3*B_rf__./w_rf)*B0(1,1) - ...
         gamma_n*besselj(0,gamma_n*B_rf__./w_rf)*B0(1,1);

plot(B_rf__,w_diff)
hold on
xlabel('B_{rf}')
ylabel('w_{He3} - w_{n}')
plot([B_rf-(0.5e-4) B_rf+(0.5e-4)],[adj adj])

syms B_rf_
f_d = gamma_3*besselj(0,gamma_3*B_rf_/w_rf)*B0(1,1) - ...
      gamma_n*besselj(0,gamma_n*B_rf_/w_rf)*B0(1,1);% + adj;
B_rf_adj = vpasolve(f_d,B_rf_,[B_rf-(2e-4) B_rf+(2e-4)]);

w_n_adj = gamma_n*besselj(0,gamma_n*B_rf_adj/w_rf)*B0(1,1)
w_He_adj = gamma_3*besselj(0,gamma_3*B_rf_adj/w_rf)*B0(1,1)


B_rf_adj
(B_rf_adj - B_rf)*0.5 + B_rf
