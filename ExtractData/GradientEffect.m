
fm = 1;
w_rf = 100000;
B_rf = 6.4885411e-4;
wrf_amp = 16589.183;
n_pow = 17;
scale1 = 0.00835;
scale2 = 0;
phi = pi/2;

L_z = 0.4;
L_x = 0.038;

Bxdx = 0;
Bzdz = 0;

B_rf_plus = B_rf + L_z/2 * Bzdz;
B_rf_minus = B_rf - L_z/2 * Bzdz;
B0_plus = B0(1,1) + L_x/2 * Bxdx;
B0_minus = B0(1,1) - L_x/2 * Bxdx;

B_rf_m = [ B_rf B_rf_plus B_rf_minus ];
B0_m = [ B0(1,1) B0_plus B0_minus ];

step = 0.0001;
t = 0:step:1/fm;
wfm = 2*pi*fm;

angles = zeros(length(t),3);
phase_n = 0;
phase_He = 0;
j = 1;
for B_rf_i = B_rf_m
    for B0_i = B0_m
        i = 1;
        for t_i = t
            rad = (t_i+step/2)*wfm + phi;
            w_rf_t = w_rf + wrf_amp*ModulationFunc(rad, n_pow, scale1 );
            phase_n = phase_n + gamma_n*besselj(0,gamma_n*B_rf_i/w_rf_t)*B0_i*step;
            phase_He = phase_He + gamma_3*besselj(0,gamma_3*B_rf_i/w_rf_t)*B0_i*step;
            pd_m = mod(abs(phase_n - phase_He),2*pi);
            angles(i,j) = pd_m;
            i = i + 1;
        end
        j = j + 1;
    end
end

form = ModulationFunc(t(:)*wfm+phi, n_pow, scale1 );
plot(t,angles(:,1));
hold on
for k = 2:(j-1)
    plot(t,angles(:,k));
end
plot([0 1/fm], [0.8 0.8])
%plot(t,form)
xlabel('t /s')
ylabel('angle /rad')
%legend('wrf','wrf+','wrf-');



