function F = FindDress_func( x )
% This function is used in FindDressingParsAdv.m
% x(1) is w_rf
% x(2) is deltat
%{
global gamma_n;
global gamma_3;
global change;
global B0;
global B_rf;
%}
B0 = 30e-7;         
B_rf = 6.488345e-4;
change = -0.8*2;
gamma_n = -1.83247185e8;
gamma_3 = -2.037894730e8;

F(1) = (gamma_n*besselj(0,gamma_n*B_rf/x(1)) - gamma_3*besselj(0,gamma_3*B_rf/x(1)) )*B0*x(2) - change;
F(2) = mod(gamma_n*besselj(0,gamma_n*B_rf/x(1))*B0*x(2),2*pi)-change/2;
end

