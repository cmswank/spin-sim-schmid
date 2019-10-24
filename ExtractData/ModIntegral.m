function [ integral ] = ModIntegral( x, a, n )

%integral = -sqrt(pi/n)/(2*a)*( a*a*(erf(sqrt(n)*(x-3*pi/2))+ erf(3*pi*sqrt(n)/2) ) + erf(0.5*sqrt(n)*(pi-2*x)) - erf(pi*sqrt(n)/2));
integral = -0.5*sqrt(pi/n)*(1/a*erf(0.5*sqrt(n)*(pi-2*x))+a*erf(sqrt(n)*(x-3*pi/2)));

end

