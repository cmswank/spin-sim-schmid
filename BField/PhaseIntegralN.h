#ifndef PHASEINTEGRAL_H
#define PHASEINTEGRAL_H

#define PI 3.141592653589793238462643383

double phaseModIntegral(double x, double a, double n){
	using namespace std;
	// this function return the integral from 0 to x of the function 1/a*exp(-n*(x-pi/2)^2)-a*exp(-n*(x-3*pi/2^2))
	return -sqrt(PI/n)/(2*a*2*PI)*( a*a*(erf(sqrt(n)*(x-3*PI/2))+ erf(3*PI*sqrt(n)/2) ) + erf(0.5*sqrt(n)*(PI-2*x)) - erf(PI*sqrt(n)/2));
}

double phaseCriticalMod(double x, double a, double n){
	using namespace std;
	// this function return the integral from 0 to x of the function 1/a*exp(-n*(x-pi/2)^2)-a*exp(-n*(x-3*pi/2^2))
	return -sqrt(PI/n)/(2*a*2*PI)*( a*a*(erf(sqrt(n)*(x-3*PI/2))+ erf(3*PI*sqrt(n)/2) ) + erf(0.5*sqrt(n)*(PI-2*x)) - erf(PI*sqrt(n)/2));
}



#endif