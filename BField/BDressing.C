#include "BDressing.h"
#include <math.h>
//#include "DEIntegrator.h"
#include "PhaseIntegralN.h"
#include <complex>

#define PI 3.141592653589793238462643383


/**
 * Implementation of getFactor
 */
double BDressingFactor::getFactor(BFieldVars &vars) {
	///starting to adsfget spagetti! too many cooks...
  //Initial pulse stuff was Ezra webb. C swank did the dressing, and noise

  if (pulse == 1) {
    double time = vars.t;
    if (time < T_p) {
      return Bscale*cos(time*w_p+phi_p);
    }
    else if (time < (T_p + T_pause)) {
      return 0.0;
    }
    else {
      return cos(get_phase(time - (T_p + T_pause) ));
    }
  }
  else if (pulse == 2) {
    double time = vars.t;
    if (time < (T_p-T_crop)) {
      double mu  = 5;
      double beta = width*PI/mu;
      double B0 = 3e-6;
      //double w_mean = B0*(1.83247185e8 + 2.037894730e8)/2;
      double w_n = B0*1.83247185e8;
      double w_He = B0*2.037894730e8;
      std::complex<double> expo(-1,-mu);
      std::complex<double> i(0,1);
      double t = time - T_p/2.0;
      return Bscale*std::real(rscale*std::exp(i*w_n*t)*std::pow(std::cosh(beta*t),expo)+1/rscale*std::exp(i*w_He*t)*std::pow(std::cosh(beta*t),expo));
      //return Bscale*std::real(std::exp(i*w_mean*t)*std::pow(std::cosh(beta*t),expo));
    }
    else if (time < (T_p - T_crop + T_pause)) {
      return 0.0;
    }
    else {
      return cos(get_phase(time - (T_p - T_crop +T_pause)));
    }
  }
  //noise!!!!! right now its white (gaussian) but can be defined in the cmsInterpnoiseGen.C random distribution
  else if (pulse == 3) {
    double time = vars.t;
    double interptime[1]={time};
    if (time < (T_p-T_crop)) {
      double mu  = 5;
      double beta = width*PI/mu;
      double B0 = 3e-6;
      double w_mean = B0*(1.83247185e8 + 2.037894730e8)/2;
      //double w_n = B0*1.83247185e8;
      //double w_He = B0*2.037894730e8;
      std::complex<double> expo(-1,-mu);
      std::complex<double> i(0,1);
      double t = time - T_p/2.0;

     // double tempdoublein;
     // std::cout<<"Please enter time double \n";
     // std::cin>>tempdoublein;
      //double testtest[1]={(double)tempdoublein};
      //std::cout<<"interp class start "<<interpnoise->start[0]<<"\n";
      //std::cout<<"interp class step "<<interpnoise->step[0]<<"\n";
      //std::cout<<" \n"<<this->interpnoise->interpdata[10]<<"  \n";
      //std::cout<<" \n"<<this->interpnoise->interpdata[8000]<<"  \n";
      //std::cout<<" \n"<<this->interpnoise->interp1D(interptime)<<"  \n";
      
      //return Bscale*(std::real(rscale*std::exp(i*w_n*t)*std::pow(std::cosh(beta*t),expo)+1/rscale*std::exp(i*w_He*t)*std::pow(std::cosh(beta*t),expo)+interpnoise->interp1D(interptime)));
      return Bscale*std::real((std::exp(i*w_mean*t)*std::pow(std::cosh(beta*t),expo))+interpnoise->interp1D(interptime));
    }
    else if (time < (T_p - T_crop + T_pause)) {
      return 0.0;
    }
    else {
      return cos(get_phase(time - (T_p - T_crop +T_pause)))+interpnoise->interp1D(interptime);
    }
  }

    else if (pulse == 4) {
    double time = vars.t;
    double interptime[1]={time};
    if (time < (T_p-T_crop)) {
      /*double mu  = 5;
      double beta = width*PI/mu;
      double B0 = 3e-6;
      double w_mean = B0*(1.83247185e8 + 2.037894730e8)/2;
      double w_n = B0*1.83247185e8;
      double w_He = B0*2.037894730e8;
      std::complex<double> expo(-1,-mu);
      std::complex<double> i(0,1);
      double t = time - T_p/2.0;
      */

     // double tempdoublein;
     // std::cout<<"Please enter time double \n";
     // std::cin>>tempdoublein;
      //double testtest[1]={(double)tempdoublein};
      //std::cout<<"interp class start "<<interpnoise->start[0]<<"\n";
      //std::cout<<"interp class step "<<interpnoise->step[0]<<"\n";
      //std::cout<<" \n"<<this->interpnoise->interpdata[10]<<"  \n";
      //std::cout<<" "<<this->b1Pulse->interpdata[10]<<"  ";
      //std::cout<<" "<<this->b1Pulse->interp1D(interptime)<<" ";
      
      //return Bscale*(std::real(rscale*std::exp(i*w_n*t)*std::pow(std::cosh(beta*t),expo)+1/rscale*std::exp(i*w_He*t)*std::pow(std::cosh(beta*t),expo)+interpnoise->interp1D(interptime)));
      return Bscale*(b1Pulse->interp1D(interptime)+interpnoise->interp1D(interptime));
      }
    else if (time < (T_p - T_crop + T_pause)) {
      return 0.0;
      }
    else {
      return cos(get_phase(time - (T_p - T_crop +T_pause)));
      }
    }

    else if (pulse == 5) {
    //Double Pulses! Attempting to use robust dressing to lock phases together in B1 pulse. 
      ///putting it so after the pulse it stops dressing. 
    double time = vars.t;
    double interptime[1]={time};
      /*double mu  = 5;
      double beta = width*PI/mu;
      double B0 = 3e-6;
      double w_mean = B0*(1.83247185e8 + 2.037894730e8)/2;
      double w_n = B0*1.83247185e8;
      double w_He = B0*2.037894730e8;
      std::complex<double> expo(-1,-mu);
      std::complex<double> i(0,1);
      double t = time - T_p/2.0;
      */

     // double tempdoublein;
     // std::cout<<"Please enter time double \n";
     // std::cin>>tempdoublein;
      //double testtest[1]={(double)tempdoublein};
      //std::cout<<"interp class start "<<interpnoise->start[0]<<"\n";
      //std::cout<<"interp class step "<<interpnoise->step[0]<<"\n";
      //std::cout<<" \n"<<this->interpnoise->interpdata[10]<<"  \n";
      //std::cout<<" "<<this->b1Pulse->interpdata[10]<<"  ";
      //std::cout<<" "<<this->b1Pulse->interp1D(interptime)<<" ";
      
      //return Bscale*(std::real(rscale*std::exp(i*w_n*t)*std::pow(std::cosh(beta*t),expo)+1/rscale*std::exp(i*w_He*t)*std::pow(std::cosh(beta*t),expo)+interpnoise->interp1D(interptime)));
      if(time<T_p)
        return Bscale*(b1Pulse->interp1D(interptime))+cos(get_phase(vars.t))+Bscale*interpnoise->interp1D(interptime);//Bscale*(b1Pulse->interp1D(interptime))+interpnoise->interp1D(interptime));
      else
        return 0.0;//+interpnoise->interp1D(interptime));
      }

      else if (pulse == 6) {
    //Double Pulses! Attempting to use robust dressing to lock phases together in B1 pulse. 
      ///putting it so after the pulse it stops dressing. 
    double time = vars.t;
    double interptime[1]={time};
      /*double mu  = 5;
      double beta = width*PI/mu;
      double B0 = 3e-6;
      double w_mean = B0*(1.83247185e8 + 2.037894730e8)/2;
      double w_n = B0*1.83247185e8;
      double w_He = B0*2.037894730e8;
      std::complex<double> expo(-1,-mu);
      std::complex<double> i(0,1);
      double t = time - T_p/2.0;
      */

     // double tempdoublein;
     // std::cout<<"Please enter time double \n";
     // std::cin>>tempdoublein;
      //double testtest[1]={(double)tempdoublein};
      //std::cout<<"interp class start "<<interpnoise->start[0]<<"\n";
      //std::cout<<"interp class step "<<interpnoise->step[0]<<"\n";
      //std::cout<<" \n"<<this->interpnoise->interpdata[10]<<"  \n";
      //std::cout<<" "<<this->b1Pulse->interpdata[10]<<"  ";
      //std::cout<<" "<<this->b1Pulse->interp1D(interptime)<<" ";
      
      //return Bscale*(std::real(rscale*std::exp(i*w_n*t)*std::pow(std::cosh(beta*t),expo)+1/rscale*std::exp(i*w_He*t)*std::pow(std::cosh(beta*t),expo)+interpnoise->interp1D(interptime)));
      if(time<T_p){
        //std::cout<<b1Pulse->interp1D(interptime)<<"   ";
        return b1Pulse->interp1D(interptime);//+cos(get_phase(vars.t))+Bscale*interpnoise->interp1D(interptime);//Bscale*(b1Pulse->interp1D(interptime))+interpnoise->interp1D(interptime));
      }
      else
        return 0.0;//+interpnoise->interp1D(interptime));
      }


    else
      return cos(get_phase(vars.t));
}

/**
 * Implementation of get_phase
 */
double BDressingFactor::get_phase(double time) {
	 
	return time*this->wrf + this->phi;
}

/**
 * Implementation of get_phase
 */
double BDressingCosModFactor::get_phase(double time) {
   double phase = wrf*time + phi;
    double wfm = 2*PI*fm;
  	phase += wrf_amp/wfm * sin(wfm * time);
	
	return phase;
}


/**
 * Implementation of get_phase
 */

/*
double wfm;
double scale1_;

class ModulationFunction
{
public:
  double operator()(double t) const
  {
    return pow(cos(wfm*t)+scale1_*cos(2*wfm*t),101);
  }
};


double BDressingFuncModFactor::get_phase(double time) {
  double phase = phase_t;
  double step = (time - time_p);
  phase += wrf*step;
  wfm = 2*PI*fm;
  //scale1_ =  scale1;
  //ModulationFunction modfunc;
  //double integral = DEIntegrator<ModulationFunction>::Integrate(modfunc, time_p, time, 1);
  double slices = 5;
  double sum_of_slices = 0;
  for (int i=1;i<5;++i){
    sum_of_slices += pow(cos(wfm*(time_p+i*step/slices))+scale1*cos(2*wfm*(time_p+i*step/slices)),101);
  }
  double integral = 0.5*step*(1/slices)*(pow(cos(wfm*time)+scale1*cos(2*wfm*time),101)+pow(cos(wfm*time_p)+scale1*cos(2*wfm*time_p),101)+2*sum_of_slices);
  phase += wrf_amp * integral;
  this->time_p=time;
  this->phase_t=phase;
  return phase;
}
*/
//Ezra Webb function for modulation, Functions used here are in PHaseIntegralN.h
  double BDressingFuncModFactor::get_phase(double time) {
  double phase = wrf*time + phi;
  double wfm = 2*PI*fm;
  double N = floor((wfm*time+phi_mod)/(2*PI));
  double x = fmod(wfm*time+phi_mod,2*PI);
  phase += wrf_amp*(N*phaseModIntegral(2*PI,scale1,scale2) + phaseModIntegral(x,scale1,scale2) - phaseModIntegral(phi_mod,scale1,scale2));
  return phase;
}

/*

double BDressingFuncModFactor::get_phase(double time) {
  double phase = wrf*time;
  double wfm = 2*PI*fm;
  double xend = time*wfm;
  phase += wrf_amp * phaseIntegral(xend, scale1);
  // phase = wrf_amp * integral of periodic w function
  // eg. for w = wrf + wrf_amp*cos(wfm*time)
  // phase = wrf_amp * sin(wfm * time)/wfm
  return phase;
}
*/

/**
 * Implementation of getFactor
 *
double BDressingCosModFactor::getFactor(BFieldVars &vars) {
	
	return cos(this->get_phase(vars.t));
}*/





/**
 * Implementation of getFactor
 */
double BDressingCosBModFactor::getFactor(BFieldVars &vars) {
  double wfm = 2*PI*fm;
  double scale = 1. + amp * cos(wfm * vars.t);
  return scale * BDressingFactor::getFactor(vars);
}

/**
 * Implementation of getFactor
 */
double BDressingPulsedBModFactor::getFactor(BFieldVars &vars) {
  int halves = floor(vars.t * 2 * fm);
  double mod_t = vars.t - 0.5/fm * halves;
  double scale;
  if (mod_t < deltat)
    scale = (halves % 2 == 0) ? 1.+scale1 : 1.+scale2;
  else
    scale = 1.;
  return scale * BDressingFactor::getFactor(vars);
}

/**
 * Implementation of get_phase
 */
double BDressingPulsedFreqModFactor::get_phase(double time) {
  ///THis has been changed to Critical or Robust Modulation!
  ////Critical modulation is the idea I had where you take 
  ////advantage of the amount of dressing as a function of
  //// the angle of the spins because the field is only along one 
  /// dimension. 
  ////This means that we might be able to create an semi-infinite T2 
  /// by slowing the spins in the leading part of the "fan"
  /// and speeding up the slower part. 
  //dw1 is amplitude of critical mod, fm is frequeny of critical modulation.
  //dw2 is phase of critrical modulation. phi is phase of dressing?
  //wrf is the usual omega rf. 
  
  return dw1/fm*(1.0-cos(fm*time+dw2))+wrf*time+phi; 


//Old stuff (pulsed freq mod was shown not to converge correctly, please use funtion pulsed freq mod. 
/*  int halves = floor(time*2.0*fm);
  double halfT = 0.5/fm;
  double phase = phi + wrf*halves*halfT;
  if (halves % 2 == 0) {
    phase += halves/2 * (dw1*deltat1 + dw2*deltat2) ;
  }
  else {
    phase += ( (halves-1)/2 *(dw1*deltat1 + dw2*deltat2)+ dw1* deltat1 ) ; 
  }
  double mod_t = time - halfT * halves;
  phase += mod_t * wrf;
  if (halves % 2 == 0)
    phase += dw1 * ((mod_t < deltat1) ? mod_t : deltat1);
  else
    phase += dw2 * ((mod_t < deltat2) ? mod_t : deltat2);
  return phase;
*/
}

/**
 * Implementation of getFactor
 *
double BDressingPulsedFreqModFactor::getFactor(BFieldVars &vars) {
  return cos(this->get_phase(vars.t));
}*/

double BDressingPulsedFreqModFactor2::get_phase(double time) {
  double dw1 = wrf_amp/scale;
  double dw2 = -wrf_amp*scale;
  double phase = floor(time*fm)*(wrf/fm+dw1*deltat1+dw2*deltat2)+phi;
  double t = fmod(time,1/fm);
  if (t <= 0.5*deltat1)
    phase += t*(wrf+dw1);
  else if (t<= (1/(2*fm) -deltat2/2))
    phase += deltat1*0.5*dw1 + t*wrf;
  else if (t<= (1/(2*fm) +deltat2/2))
    phase += deltat1*0.5*dw1 + t*wrf + (t-(1/(2*fm) - deltat2/2))*dw2;
  else if (t <= (1/fm - deltat1/2))
    phase += deltat1*0.5*dw1 + t*wrf + deltat2*dw2;
  else
    phase += deltat1*0.5*dw1 + t*wrf + deltat2*dw2 + (t - (1/fm - 0.5*deltat1))*dw1; 
  return phase;
}
