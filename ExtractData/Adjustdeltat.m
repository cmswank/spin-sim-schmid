[~,peak_n] = max(abs(fft(s_n(3,:),nBins*100)));
w_n_act = -peak_n *2*pi/(runTime*100)
w_n_theory = gamma_n*besselj(0,gamma_n*B_rf/w_rf)*B0(1,1);

[~,peak_He] = max(abs(fft(s_He(3,:),nBins*100)));
w_He_act = -peak_He *2*pi/(runTime*100)
w_He_theory = gamma_3 * besselj(0,gamma_3*B_rf/w_rf)*B0(1,1);

angle_dw1 = 0.0593;
angle_dw2 = 0;

