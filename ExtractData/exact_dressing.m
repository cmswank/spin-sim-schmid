
w_rf = 10000;
B_rf_r = 1.189*w_rf/gamma_n;    %B_rf rough solution

syms B_rf_
f_d = gamma_3*besselj(0,gamma_3*B_rf_/w_rf)*B0(1,1) - ...
      gamma_n*besselj(0,gamma_n*B_rf_/w_rf)*B0(1,1);
B_rf_e = vpasolve(f_d,B_rf_,[B_rf_r-(2e-4) B_rf_r+(2e-4)]) %B_rf exact sol.