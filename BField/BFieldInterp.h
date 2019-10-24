#ifndef __BFieldInterp_H__
#define __BFieldInterp_H__
#include "BField/BField.h"
//#include <iostream>
//#include <assert.h>
//#include <cstdlib>
//#include <fstream>
//#include <sstream>
//#include <stdio.h>
//#include <boost/thread.hpp>
//#include <boost/math/special_functions.hpp>
//#include <boost/random.hpp>
//#include <boost/random/uniform_real_distribution.hpp>


class BFieldInterp : public BFactor  {
private:
	

	double cubicInterpolate (double p[4], double x);
	double nCubicInterpolate (int n, double* p, double coordinates[]);
	//Switching to ROOT random number generator, because it must work with Riccardo's code!
	//double bicubicInterpolate (double p[4][4], double x, double y);
	//double tricubicInterpolate (double p[4][4][4], double x, double y, double z);
	//boost random engine. 
	/*typedef boost::mt19937 RNG;    // Mersenne Twister
 	typedef boost::random::uniform_real_distribution<double> DIST_flat;
	typedef boost::variate_generator<RNG, DIST_flat> RAND_flat;

	typedef boost::normal_distribution<double> DIST_velocity;   // Normal Distribution
	typedef boost::variate_generator<RNG,DIST_velocity> RAND_velocity;
	RNG rng;
	DIST_flat dist_flat;
	RAND_flat rand_flat;
	DIST_velocity dist_vel;
	RAND_velocity rand_vel;
    */
    //also root seems easier...
	//TRandom3* rootRand;
public:
	//this defines the data in flattened C type array, (last column moves fastest) and the step size for all arrays. 
	//void whitenoiseGen(double nmag);
    
	BFieldInterp(bool fromfile, int direction, double* star,double* ste, int N[3], double *datainput);
    int Num[3];
	double testtest;
	double *interpdata;
    double *step;//=new double[3];
	double *start;//=new double[3];
	int dsize;
	

	//this is the function you want to interp. must be from a rectangular grid. 
	//must be an interpolation (no checking done to increase speed, it will extrapolate)
	//this is a tricubic spline interpolation. 
	//if its too slow I guess it could be made into a linear with ease. 
	double interp1D(double* pos);
	double interp3D(double* pos);
	double getFactor(BFieldVars &vars);




};

#endif