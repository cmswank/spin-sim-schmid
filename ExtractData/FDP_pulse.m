

fm = 1;
w_rf = 100000;
B_rf = 1.189*w_rf/gamma_n;
deltat1 = 0.5*0.1/fm;


deltat2 = deltat1;

syms w
f_l = (gamma_n * besselj(0,gamma_n*B_rf/w) - gamma_3 * besselj(0,gamma_3*B_rf/w) )*B0(1,1)*deltat1 - 2*dphi;
f_u = (gamma_n * besselj(0,gamma_n*B_rf/w) - gamma_3 * besselj(0,gamma_3*B_rf/w) )*B0(1,1)*deltat2 + 2*dphi;
w1 = vpasolve(f_l,w, [w_rf 3*w_rf]);
dw1 = w1 - w_rf;
w2 = vpasolve(f_u,w, [w_rf/2 w_rf]);
dw2 = w2 - w_rf;

%{
step = 0.001;
t = 0:step:1/fm;
wfm = 2*pi*fm;

wrf_amp_r = 10000:10:18000;
pars_m = zeros(length(wrf_amp_r),5);

phase_diff_min = 100000000;
i=1;
for wrf_amp_i = wrf_amp_r
    phase_n = 0;
    phase_He = 0;
    pd_m_max = 0;
    pd_m_min = 0;
    for t_i = t
        w_rf_t = w_rf + wrf_amp_i*(cos((t_i+step/2)*wfm)^n_pow+scale1*cos((t_i+step/2)*wfm/2)+scale2*sin(2*(t_i+step/2)*wfm));
        phase_n = phase_n + gamma_n*besselj(0,gamma_n*B_rf/w_rf_t)*B0(1)*step;
        phase_He = phase_He + gamma_3*besselj(0,gamma_3*B_rf/w_rf_t)*B0(1)*step;
        pd_m = phase_n - phase_He;
        if pd_m > pd_m_max
            pd_m_max = pd_m;
        elseif pd_m < pd_m_min
            pd_m_min = pd_m;
        end 
    end
    phase_diff = abs(phase_n - phase_He);
    if phase_diff < phase_diff_min
        phase_diff_min = phase_diff;
        i_min = i;
    end
    pars_m(i,:) = [wrf_amp_i phase_n phase_He pd_m_max pd_m_min];
    i = i +1;
end

disp(pars_m(i_min,1))
disp(pars_m(i_min,4))
disp(pars_m(i_min,5))

plot(pars_m(:,1),pars_m(:,2));
hold on
plot(pars_m(:,1),pars_m(:,3));
legend('n','He')
hold off

disp('make sure mod function is correct')

%}

step = 0.001;
t = 0:step:1/fm;
angle_pred = zeros(length(t));
phase_n = 0;
phase_He = 0;
halfT = 0.5/fm;
i = 1;
for t_i = t
    halves = floor(t_i*2.0*fm);
    mod_t = t_i - halfT * halves;
    w_ = w_rf;
    if (mod(halves,2) == 0 && mod_t < deltat1)
        w_ = w_rf + dw1;
    elseif (mod(halves,2) == 1 && mod_t < deltat2)
        w_ = w_rf + dw2;
    end
    phase_n = phase_n + gamma_n*besselj(0,gamma_n*B_rf/w_)*B0(1)*step;
    phase_He = phase_He + gamma_3*besselj(0,gamma_3*B_rf/w_)*B0(1)*step;
    pd_m = mod(abs(phase_n - phase_He),2*pi);
    angle_pred(i) = pd_m;
    i = i + 1;
end

plot(t,angle_pred);
hold on
plot([0 1/fm], [0.8 0.8])
xlabel('t /s')
ylabel('angle /rad')

