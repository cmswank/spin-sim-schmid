
%this script requires that an unmodulated signal be preloaded

%make sure these agree
B0 = 30e-7;         
w_rf = 100000;
B_rf = 6.488345e-4;

%global gamma_n;
%global gamma_3;
%global change;
%global B0;
%global B_rf;

w_n_t = -gamma_3*besselj(0,gamma_3*B_rf/w_rf)*B0(1,1);   % theoretical w_eff value for n
w_He_t = -gamma_n*besselj(0,gamma_n*B_rf/w_rf)*B0(1,1);   % theoretical w_eff value for He
w_mean = (w_n_t + w_He_t)/2;

dphi = 0.8; % ang sep of He and n

change = -2*dphi; % this is +ve for n overtake and -ve for He overtake
                  % currently we start with n ahead therefore dw1 is He
                  % overtake

fun = @FindDress_func;
x0 = [1.6*w_rf, 0.05];
x = fsolve(fun,x0);
dw1 = x(1) - w_rf;
deltat1 = x(2);

change = 2*dphi;

x0 = [0.7*w_rf, 0.05];
x = fsolve(fun,x0);
dw2 = x(1) - w_rf;
deltat2 = x(2);


frac_mod = 0.1; % fraction of time spent modulating
fm = 1/(floor(0.5*((deltat1+deltat2)/frac_mod - deltat1 - deltat2)/(2*pi/w_mean))*2*(2*pi/w_mean)+deltat1+deltat2);

%w_n_dw1 = B0*gamma_n*besselj(0,gamma_n*B_rf/w1)
%w_3_dw1 = B0*gamma_3*besselj(0,gamma_3*B_rf/w1)
%w_n_dw2 = B0*gamma_n*besselj(0,gamma_n*B_rf/w2)
%w_3_dw2 = B0*gamma_3*besselj(0,gamma_3*B_rf/w2)

