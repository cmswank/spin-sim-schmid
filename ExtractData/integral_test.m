n = 51;
a = 0.732;
phi_mod = 0;
fm = 1;
wfm= 2*pi*fm;

step = 0.001;
t= 0:step:1;
out = zeros(length(t),2);
i=1;
integral = 0;
for t_i = t
    N = floor((wfm*t_i+phi_mod)/(2*pi));
    x = mod(wfm*t_i+phi_mod,2*pi);
    phase =  ModIntegral(x,a,n) - ModIntegral(0,a,n) ;%+N*ModIntegral(2*pi,a,n)- ModIntegral(phi_mod,a,n);
    integral = integral + 0.5*step*(ModulationFunc(t_i*wfm,n,a)+ModulationFunc((t_i+step)*wfm,n,a));
    out(i,1) = phase;
    out(i,2) = integral;
    i = i + 1;
end

plot(out(:,1))
hold on
plot(out(:,2))
legend('analytic','sum')




