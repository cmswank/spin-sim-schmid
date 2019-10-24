
fm = 1;
w_rf = 100000;
B_rf = 6.4885411e-4;
wrf_amp =  16591;

n_pow = 17;
scale1 = 0.00835;

step = 0.0001;
t = 0:step:1/fm;
wfm = 2*pi*fm;

angle_pred = zeros(length(t));
phase_n = 0;
phase_He = 0;
pd_m_max = 0;
pd_m_min = 0;
i = 1;
for t_i = t
    w_rf_t = w_rf + wrf_amp*ModulationFunc(t_i+step/2, wfm, n_pow, scale1 );
    phase_n = phase_n + gamma_n*besselj(0,gamma_n*B_rf/w_rf_t)*B0(1)*step;
    phase_He = phase_He + gamma_3*besselj(0,gamma_3*B_rf/w_rf_t)*B0(1)*step;
    pd_m = phase_n - phase_He;
    if pd_m > pd_m_max
        pd_m_max = pd_m;
    elseif pd_m < pd_m_min
        pd_m_min = pd_m;
    end 
    i = i + 1;
    angle_pred(i) = abs(pd_m);
end

disp(pd_m_max)
disp(pd_m_min)

%{
plot(t,angle_pred);
hold on
plot([0 1/fm], [0.8 0.8])
xlabel('t /s')
ylabel('angle /rad')
%}

