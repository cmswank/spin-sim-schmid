

fm = 1;
w_rf = 10000;
B_rf = 6.48854119e-5;

n_pow = 51;
scale1 = 0.726;
phi = pi/2;

offset = 0;

step = 0.0001;
t = 0:step:1/fm;
wfm = 2*pi*fm;

wrf_amp_r = 4240:0.5:4250;
pars_m = zeros(length(wrf_amp_r),5);

phase_diff_min = 100000000;
i=1;
i_min = 1;
for wrf_amp_i = wrf_amp_r
    phase_n = 0;
    phase_He = 0;
    pd_m_max = 0;
    pd_m_min = 0;
    for t_i = t
        rad = (t_i+step/2)*wfm + phi;
        w_rf_t = w_rf + wrf_amp_i*ModulationFunc(rad, n_pow, scale1 );
        phase_n = phase_n + gamma_n*besselj(0,gamma_n*B_rf/w_rf_t)*B0(1)*step;
        phase_He = phase_He + gamma_3*besselj(0,gamma_3*B_rf/w_rf_t)*B0(1)*step;
        pd_m = (phase_n - phase_He);
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
figure;



wrf_amp = pars_m(i_min,1);
angle_pred = zeros(length(t),1);
phase_n = 0;
phase_He = 0;
i = 1;
for t_i = t
    rad = (t_i+step/2)*wfm + phi;
    w_rf_t = w_rf + wrf_amp*ModulationFunc(rad, n_pow, scale1 );
    phase_n = phase_n + gamma_n*besselj(0,gamma_n*B_rf/w_rf_t)*B0(1)*step;
    phase_He = phase_He + gamma_3*besselj(0,gamma_3*B_rf/w_rf_t)*B0(1)*step;
    angle_pred(i,1) = abs(phase_n - phase_He);
    i = i + 1;
end

plot(t,angle_pred);
hold on
plot([0 1/fm], [0.8 0.8])
xlabel('t /s')
ylabel('angle /rad')


