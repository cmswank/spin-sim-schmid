//g++ cmsInterpnoiseGen.cpp -I/data1/cmswank/BoostCodeSwank/boost_1_64_0 -o interpNoiseTest -lboost_thread -lboost_system `root-config --cflags --glibs`
#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdio.h>
//#include <boost/thread.hpp>
//#include <boost/math/special_functions.hpp>
//#include <boost/random.hpp>
//#include <boost/random/uniform_real_distribution.hpp>
#include "TRandom3.h"
#include "cmsInterpnoiseGen.h"


//example of ROOT random (is this all?!?!?!?, looks easier than boost. 
//gRandom->Circle(*position, *(position+1), r);


cmsInterp::cmsInterp(double star[1],double ste[1], int N[3], int flat_seed,int freq_cut,int hpass)/*:rng(flat_seed),dist_flat(0.0f,1.0f),rand_flat(rng,dist_flat),\
							dist_vel(0.0f,6.0f),rand_vel(rng,dist_vel)*/{
		Num[0]=N[0];
		Num[1]=N[1];
		Num[2]=N[2];
		dsize=N[0]*N[1]*N[2];
		interpdata =(double*) malloc(dsize*sizeof *interpdata);
		start[0]=star[0];
		step[0]=ste[0];
		rnd_seed=flat_seed;
		this->freq_cut=freq_cut;
		this->hpass=hpass;
		this->rootRand= new TRandom3(flat_seed);
}




double cmsInterp::cubicInterpolate (double p[4], double x) {
	
	//std::cout<<"p0 "<<p[0]<<" \n"<<"p1 "<<p[1]<<" \n"<<"p2 "<<p[2]<<" \n"<<"p3 "<<p[3]<<"\n";
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}
/*
double cmsInterp::bicubicInterpolate (double p[4][4], double x, double y) {
	double arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);
}
*/


double cmsInterp::interp1D(double *pos)
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

/*
double cmsInterp::interp3D(double *pos)
{
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
	
		tempp[i][ii][iii]=interpdata[(i+offx)*Num[1]*Num[2]+(ii+offy)*Num[2]+(offz+iii)];

	}}}

	
	return nCubicInterpolate(3,(double*) tempp,tempos);// tricubicInterpolate(tempp,tempx,tempy,tempz);

}

/*
double cmsInterp::tricubicInterpolate (double p[4][4][4], double x, double y, double z) {
	double arr[4];
	arr[0] = bicubicInterpolate(p[0], y, z);
	arr[1] = bicubicInterpolate(p[1], y, z);
	arr[2] = bicubicInterpolate(p[2], y, z);
	arr[3] = bicubicInterpolate(p[3], y, z);
	return cubicInterpolate(arr, x);
}
*/

double cmsInterp::nCubicInterpolate (int n, double* p, double coordinates[]) {
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

void cmsInterp::whitenoiseGen(double nmag){
	
	//DIST_velocity dist_vel(0.0f,nmag);
	//rand_vel.distribution()=dist_vel;
	
	//double tempdat[dsize];
	//std::cout<<"creating noise array of size "<<dsize<<"\n";
	//double tempstep=std::sqrt(step[0]);
	//double tempinterpdata[dsize];

	///using biquadratic filter to simulate analog filtering around the resonance of he-3 and neutorns. 
	///this is using the definition found in MATLAB ('sosfilt' function). 
	//double b0[10]=[0.851562784340117, 1, 1, 1, 1, 1, 1, 1, 1, 1];
	//double b1[10]=[-1.70236224840832,-1.99876038823638,-1.99897732511815,-1.99810493616397,-1.99844246601447,-1.99722738899032,-1.99711162293135,-1.99745942235603,-1.99775196317439,-1.99811022445441];
    //double b2[10]=[0.852406387524562,1.00062149657846,1.00089415908201,1.00000000000000,1.00037893782208,0.999145240070884,0.998999999501539,0.999321513660361,0.999688789508662,0.999961510851104];
	//a0 are all one in the is formalism. 
	//double a1[10]=[-1.94058412426200,-1.94411442994908,-1.94702414507003,-1.95522176524410,-1.96194610847760,-1.96759196042180,-1.97819531998354,-1.98723515502367,-1.98271163565115,-1.99541171151643];
	//double a2[10]=[0.943918717791235,0.946409900199799,0.951366609744894,0.956711356790862,0.967080205168559,0.968622134370657,0.978993452772803,0.987918616653233,0.988308545436206,0.996047957191888]
	double sp[1]={1.0};
	//for(int i = 0; i<3; i++) step[i]=1.0;

	double sr[1]={0.0};
	//for(int i = 0; i<3; i++) start[i]=0.0;
	

	/// must be 3 (for 1D,2D and 3D if not using other dimensions put in 1.)
	int Numb[3]={2*freq_cut,1,1};

	cmsInterp* wnShort=new cmsInterp(sr,sp,Numb,rnd_seed,freq_cut,hpass);

	for(int i = 0; i<2*freq_cut;i++){
		wnShort->interpdata[i]=rootRand->Gaus(0.0f,nmag);
	}

	//std::cout<<"this will seg fault now \n";
	
	
	double dt=step[0];

	double timep[1];
	
	Double_t *randnoise=new Double_t[dsize];

	for(int i = 0; i<dsize; i++){

		timep[0]= ((double)i)*2.*((double)freq_cut)/((double)dsize);
		///normalized to SNR=1, (or simply multiply by 1/SNR to get the correct noise amplitude)   
		randnoise[i]=wnShort->interp1D(timep);
		
	}	
	
	//std::cout<<" In_noise "<<randnoise[0]<<", ";
	Int_t n_size = dsize;
	Double_t *re_noise = new Double_t[n_size];
	Double_t *im_noise = new Double_t[n_size];

	
	TVirtualFFT *fft = TVirtualFFT::FFT(1, &n_size, "R2C ES K");
	fft->SetPoints(randnoise);
	fft->Transform();
	fft->GetPointsComplex(re_noise,im_noise);
	//std::cout<<"re_noise "<<re_noise[82]<<", im_noise "<<im_noise[82]<<"\n";
	

	//fftw_complex* Snoise = (fftw_complex*) fftw_malloc(sizeof(double)*dsize);
	//fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(double)*dsize);
	//fftw_plan pFFT = fftw_plan_dft_r2c_1d(dsize, interpdata, Snoise, FFTW_ESTIMATE);
	//fftw_execute(pFFT);

	//perfect highpass filter...
	//hpass=[zeros(bndstp(1)*t(end)-1,1);zeros((bndstp(2)-bndstp(1))*t(end),1);ones(length(t)-2*bndstp(2)*t(end)+1,1);zeros((bndstp(2)-bndstp(1))*t(end),1);zeros(bndstp(1)*t(end),1)];
	for(int i =0; i<dsize; i++){
		if((int)i/dt/dsize<hpass || (int)(dsize-i)/dt/dsize<hpass){
			re_noise[i]=0.;
			im_noise[i]=0.;
			//std::cout<<i<<"  ";
			
		}
		
	}
	
	TVirtualFFT *ifft = TVirtualFFT::FFT(1, &n_size, "C2R ES K");
	
	ifft->SetPointsComplex(re_noise,im_noise);
	ifft->Transform();
	ifft->GetPoints(randnoise);
	
	//fftw_plan pIFFT = fftw_plan_dft_1d(dsize, out, Snoise, FFTW_BACKWARD, FFTW_ESTIMATE);
	//fftw_execute(pIFFT);

	for(int i =0; i<dsize; i++){
		interpdata[i]=randnoise[i]/((double)dsize);
	}

	//std::cout<<" Out_noise "<<interpdata[0]<<"\n";
	delete fft;
	delete ifft;
	//fftw_destroy_plan(pFFT);
	//fftw_destroy_plan(pIFFT);
	return;

}
/*
int main () {
	
	////fill data structure that will be interpolated.
	/*double dat[1000];
	int iiii=0;
	for(int i = 0; i<10; i++){
		for(int ii = 0; ii<10; ii++){
			for(int iii = 0; iii<10; iii++){
			
				dat[iiii]=(double)(i+ii+iii);
				iiii++;
			}
		}
	}


	double step[3];
	for(int i = 0; i<3; i++) step[i]=1.0;

	double start[3];
	for(int i = 0; i<3; i++) start[i]=0.0;
	
	/*int N[3]={10,10,10};
	cmsInterp* Bxi=new cmsInterp(dat,start,step,N);
	double position[3]={5.5,4.4,8.2};
	std::cout<<"interpolation gives "<<Bxi->interp3D(position)<<"\n";
	*/
	
		////fill data structure that will be interpolated.
	/*double dat[1000];
	int iiii=0;
	for(int i = 0; i<1000; i++){

		dat[i]=(double)i;
	}




	double step[1]={0.01};
	//for(int i = 0; i<3; i++) step[i]=1.0;

	double start[1]={0.0};
	//for(int i = 0; i<3; i++) start[i]=0.0;
	

	/// must be 3 (for 1D,2D and 3D if not using other dimensions put in 1.)
	int N[3]={1000,1,1};

	cmsInterp* Bxi=new cmsInterp(start,step,N,123);
	//std::cout<<"this will seg fault now \n";
	Bxi->whitenoiseGen(1E-2);
	
	for(int i; i<100;i++)
	{
		double position[1]={(double)i/102};
	std::cout<<"interpolation gives "<<Bxi->interp1D(position)<<"\n";
	}	

}*/
