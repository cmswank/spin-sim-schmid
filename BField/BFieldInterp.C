//g++ BFieldInterp.C -o InterpTest
#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
//#include <stdio.h>
#include <string>
//#include <boost/thread.hpp>
//#include <boost/math/special_functions.hpp>
//#include <boost/random.hpp>
//#include <boost/random/uniform_real_distribution.hpp>
//#include "TRandom3.h"
#include "BField/BFieldInterp.h"


//implimentation of getFactor (required as this is inherits Riccardo's BFactor class)
double BFieldInterp::getFactor(BFieldVars &vars) {
   //std::cout<<vars.pos[0]<<"  "<<vars.pos[1]<<"  "<<vars.pos[2]<<"  "<<this->interp3D(vars.pos)<<"\n";
	
	return this->interp3D(vars.pos);



}



/// either 1D for time or 3D for position. if you need both use two. 
BFieldInterp::BFieldInterp(bool fromfile,int direction, double *star, double *ste, int N[3], double *datainput):step(ste),start(star)
						  /*:rng(flat_seed),dist_flat(0.0f,1.0f),rand_flat(rng,dist_flat),\
							dist_vel(0.0f,6.0f),rand_vel(rng,dist_vel)*/{
	
		Num[0]=N[0];
		Num[1]=N[1];
		Num[2]=N[2];
		dsize=N[0]*N[1]*N[2];
		int interpdim=0;
		
		for(int i = 0; i<3; i++)
		{
			if(N[i]>1)
				interpdim=interpdim+1;
		}
		
		
		interpdata =(double*) malloc(dsize*sizeof *interpdata);
		step =(double*) malloc(interpdim*sizeof *step);
		start=(double*) malloc(interpdim*sizeof *start);
		
		for(int i = 0; i<interpdim;i++)
		{
			step[i]=ste[i];
			start[i]=star[i];
		}
		
		double* interptempdata;

		//put binary file here. Or build your own way to input data. I use matlab to turn comsol file into binary rectangular grid. 
		if(fromfile){
		if (direction==0){
		std::ifstream B1file ("/data1/cmswank/spin-sim-xliu/BField/Bdatax.dat", std::ios::binary | std::ios::in);
		char* buff=new char [dsize*sizeof(double)];
		B1file.seekg(0,std::ios::beg);
		B1file.read(buff, dsize*sizeof(double));
		interptempdata =(double*)buff;	
		B1file.close();
			}
		
		if (direction==1){
		std::ifstream B1file ("/data1/cmswank/spin-sim-xliu/BField/Bdatay.dat", std::ios::binary | std::ios::in);
		char* buff=new char [dsize*sizeof(double)];
		B1file.seekg(0,std::ios::beg);
		B1file.read(buff, dsize*sizeof(double));
		interptempdata =(double*)buff;
		B1file.close();
			}
		
		if (direction==2){
		std::ifstream B1file ("/data1/cmswank/spin-sim-xliu/BField/Bdataz.dat", std::ios::binary | std::ios::in);
		char* buff=new char [dsize*sizeof(double)];
		B1file.seekg(0,std::ios::beg);
		B1file.read(buff, dsize*sizeof(double));
		interptempdata =(double*)buff;
		B1file.close();
			}
		
			


		//char buff[sizeof(double)];
		 //this will read data1 in array, you need to make it to lon
		//int tempb1size= B1file.tellg();
		//std::cout<<tempb1size<<"\n";



		for(int i = 0; i<dsize; i++){
		//must transfer data, not just use pointer or it will be lost for whatever dumb c++ reason. 
		interpdata[i]=interptempdata[i];
		

		//if(i<100) std::cout<<interpdata[i]<<"    ";
		
		//if(i==25000) std::cout<<"\n"<<interpdata[i]<<"    ";
		}

		
	}
	else
	{   ////build from data not in a file. input into the class on construction.  

		for(int i = 0; i<dsize; i++){
		//must transfer data, not just use pointer or it will be lost for whatever dumb c++ reason. 
		interpdata[i]=datainput[i];
		

		//if(i<100) std::cout<<interpdata[i]<<"    ";
		
		//if(i==25000) std::cout<<"\n"<<interpdata[i]<<"    ";
		}
	}


	//int poopscoup;
	//std::cin>>poopscoup;

	

}








double BFieldInterp::cubicInterpolate (double p[4], double x) {
	
	//std::cout<<"p0 "<<p[0]<<" \n"<<"p1 "<<p[1]<<" \n"<<"p2 "<<p[2]<<" \n"<<"p3 "<<p[3]<<"\n";
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}


///can be used to construct magntude of field for pulse in time.  
//
double BFieldInterp::interp1D(double *pos)
{
	double tempx=((pos[0]-start[0])/step[0]);

	/*std::cout<<"interp input pos "<<pos[0]<<"\n";
	std::cout<<"interp start value "<<start[0]<<"\n";
	std::cout<<"interp step value "<<step[0]<<"\n";
	std::cout<<"temp x value "<<tempx<<"\n";*/
	double tempp[4];
	int indx=(int)(tempx);
	
	int offx=0;
	if(indx>1){
		offx=indx-1;
		if(indx>Num[0]-3) offx=(Num[0]-4);	
		tempx=tempx-(double)offx;
	}

	//std::cout<<"index "<<offx<<"\n";
	tempx=tempx-1.0;

    double tempos[1]={tempx};

	for(int iii = 0; iii<4; iii++){
		//std::cout<<"interp data "<<interpdata[offx+iii]<<"    \n";
		tempp[iii]=interpdata[offx+iii];
		//std::cout<<"tempdata "<<tempp[iii]<<"    \n";
	}

	//return 0.0;
	return nCubicInterpolate(1,(double*) tempp,tempos); // tricubicInterpolate(tempp,tempx,tempy,tempz);

}


///////////////////////////3D

///position interpolation
double BFieldInterp::interp3D(double *pos)
{
	//std::cout<<step[0]<<" "<<step[1]<<"  "<<step[2]<<"\n";
	//std::cout<<start[0]<<" "<<start[1]<<"  "<<start[2]<<"\n";
	//std::cout<<Num[0]<<" "<<Num[1]<<"  "<<Num[2]<<"\n";
	double tempx=((pos[0]-start[0])/step[0]);
	double tempy=((pos[1]-start[1])/step[1]);
	double tempz=((pos[2]-start[2])/step[2]);
	
	double tempp[4][4][4];
	int indx=(int)(tempx);
	int indy=(int)(tempy);
	int indz=(int)(tempz);
	int offx=0,offy=0,offz=0;
	if(indx>1){
		offx=indx-1;
		if(indx>Num[0]-3) offx=(Num[0]-4);	
		tempx=tempx-(double)offx;
	}

	if(indy>1){
		offy=indy-1;
		if(indy>Num[1]-3) offy=(Num[1]-4);	
		tempy=tempy-(double)offy;
	
	}
	
	if(indz>1){
		offz=indz-1;
		if(indz>Num[2]-3) offz=(Num[2]-4);	
		tempz=tempz-(double)offz;
	}

	tempx=tempx-1.0;
	tempy=tempy-1.0;
	tempz=tempz-1.0;
    double tempos[3]={tempx,tempy,tempz};

	for(int i = 0; i<4; i++){for(int ii = 0; ii<4; ii++){for(int iii = 0; iii<4; iii++){
		//    x   y   z
		tempp[i][ii][iii]=interpdata[(i+offx)*Num[1]*Num[2]+(ii+offy)*Num[2]+(offz+iii)];

	}}}

	
	return nCubicInterpolate(3,(double*) tempp,tempos);// tricubicInterpolate(tempp,tempx,tempy,tempz);

}


//a different way. 

/*
double BFieldInterp::tricubicInterpolate (double p[4][4][4], double x, double y, double z) {
	double arr[4];
	arr[0] = bicubicInterpolate(p[0], y, z);
	arr[1] = bicubicInterpolate(p[1], y, z);
	arr[2] = bicubicInterpolate(p[2], y, z);
	arr[3] = bicubicInterpolate(p[3], y, z);
	return cubicInterpolate(arr, x);
}
*/

double BFieldInterp::nCubicInterpolate (int n, double* p, double coordinates[]) {
	assert(n > 0);
	if (n == 1) {
		return cubicInterpolate(p, *coordinates);
	}
	else {
		double arr[4];
		int skip = 1 << (n - 1) * 2;
		arr[0] = nCubicInterpolate(n - 1, p, coordinates + 1);
		arr[1] = nCubicInterpolate(n - 1, p + skip, coordinates + 1);
		arr[2] = nCubicInterpolate(n - 1, p + 2*skip, coordinates + 1);
		arr[3] = nCubicInterpolate(n - 1, p + 3*skip, coordinates + 1);
		return cubicInterpolate(arr, *coordinates);
	}
}

/*void BFieldInterp::whitenoiseGen(double nmag){
	
	//DIST_velocity dist_vel(0.0f,nmag);
	//rand_vel.distribution()=dist_vel;
	
	//double tempdat[dsize];
	//std::cout<<"creating noise array of size "<<dsize<<"\n";



	for(int i = 0; i<dsize; i++){
		
		interpdata[i]=rootRand->Gaus(0.0f,nmag);
		
	}	
	//interpdata=tempdat;
}
*/


/////Testing the class.

/*
int main () {
	//3D test. 
	////fill data structure that will be interpolated.
	double dat[1000];
	int iiii=0;
	for(int i = 0; i<10; i++){
		for(int ii = 0; ii<10; ii++){
			for(int iii = 0; iii<10; iii++){
			
				dat[iiii]=(double)(i+ii+iii);
				iiii++;
			}
		}
	}


	double step[3];  ///step size (spacing) between points for each dimension, 
	for(int i = 0; i<3; i++) step[i]=1.0;

	double start[3];  //starting point for each dimension. 
	for(int i = 0; i<3; i++) start[i]=0.0;
	
	int N[3]={10,10,10}; //number of points for each dimension.
	///                                  
	
	BFieldInterp* Bxi=new BFieldInterp(false,start,step,N,dat);
	
	double position[3]={5.5,4.4,8.2};
	std::cout<<"3D interpolation gives "<<Bxi->interp3D(position)<<"\n";
	
       
//       1D test.

	
		////fill data structure that will be interpolated.
	double dat1d[1000];
	
	for(int i = 0; i<1000; i++){

		dat1d[i]=(double)i;
	}



	double step1d[1]={0.01};
	//for(int i = 0; i<3; i++) step[i]=1.0;

	double start1d[1]={0.0};
	//for(int i = 0; i<3; i++) start[i]=0.0;
	

	/// must be 3 (for 1D,2D and 3D if not using other dimensions put in 1.)
	int N1d[3]={1000,1,1};

	BFieldInterp* Bxi1d=new BFieldInterp(false,start1d,step1d,N1d,dat1d);
	//std::cout<<"this will seg fault now \n";
	//Bxi->whitenoiseGen(1E-2);
	
	for(int i; i<7;i++)
	{
		double position[1]={(double)i/7.0};
	std::cout<<"1D interpolation gives "<<Bxi1d->interp1D(position)<<"\n";
	}	
	
}
*/