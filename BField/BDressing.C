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
              {
                //std::cout<<"hello"<<std::endl;
              //std::cout<<"Bscale:  "<<Bscale<<" B1 pulse interp: "<<b1Pulse->interp1D(interptime)<<" noise "<<Bscale*interpnoise->interp1D(interptime)<<std::endl;
              //std::cout<<Bscale*(b1Pulse->interp1D(interptime))+cos(get_phase(vars.t))+Bscale*interpnoise->interp1D(interptime)<<std::endl;//Bscale*(b1Pulse->interp1D(interptime))+interpnoise->interp1D(interptime));
        return cos(get_phase(vars.t))+Bscale*(b1Pulse->interp1D(interptime));//+interpnoise->interp1D(interptime));//Bscale*(b1Pulse->interp1D(interptime))+interpnoise->interp1D(interptime));
        }
      else{
        //std::cout<<"goodbye"<<std::endl;
        return 0.0;//+interpnoise->interp1D(interptime));
        
        }
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
      //input params. 
      //[0]=pulse 6 here, [1]=width and or filenum, [2]=Bscale, [3]=T_p, [4]= T_crop; [5]=T_pause, [6]= tempb1interpNum. [7]= rng seed, [8]= white noise std. 
      //[9]=tempinterpstart (noise) [10]= tempinterpstep (noise), [11]= tempinterpNum, [12]=b1interp step, [13]= b1interp stop, 
      if(time<T_p){
         //std::cout<<b1Pulse->interp1D(interptime)<<std::endl;
         //std::cout<<time<<std::endl;
         //std::cout<<Bscale*interpnoise->interp1D(interptime)<<std::endl;
        return cos(get_phase(vars.t))+Bscale*interpnoise->interp1D(interptime);//Bscale*(b1Pulse->interp1D(interptime))+interpnoise->interp1D(interptime));
      }
      else
        return 0.0;//+interpnoise->interp1D(interptime));
      }

     /* else if (pulse == 7) {
    //Double Pulses! Attempting to use robust dressing to lock phases together in B1 pulse. 
      ///putting it so after the pulse it stops dressing.

    double time = vars.t;
    double interptime[1]={time};

        return cos(get_phase(vars.t))+interpnoise->interp1D(interptime);//+cos(get_phase(vars.t))+Bscale*interpnoise->interp1D(interptime);//Bscale*(b1Pulse->interp1D(interptime))+interpnoise->interp1D(interptime));
    
        //return 0.0;//+interpnoise->interp1D(interptime));
      }*/
      else if (pulse == 7) {
    // 7 is now the modlation that goes negative that can cancel brown noise. 
        //Thanks!
   //use as Bscale->Wamp(b) [2],// use as T-P->a [3],//use as b T_crop->n [4],
    //double time = vars.t;
    //double x = 2*PI*width*(vars.t-1./width*floor(width*vars.t));
    //double interptime[1]={time};
    //double a =;
    //double b =;
    //double n =;
        double time = vars.t;
        double interptime[1]={time};
    

     //double beta=1.681; 

     
     //3845.309408, 6283.18530718

        ///I don't know about this. I think I will omit noise for now. 
        //BModCRFunction(x,a,b,n)
        //std::cout<<BModCRFunction(x,T_p,Bscale,T_crop)<<"\n";
        //double fun=BModCRFunction(x,T_p,Bscale,T_crop);
        //return cos(width*vars.t*0.5-PI*0.25)>0? fun*cos(get_phase(vars.t))+interpnoise->interp1D(interptime):-fun*cos(get_phase(vars.t))+interpnoise->interp1D(interptime);//+cos(get_phase(vars.t))+Bscale*interpnoise->interp1D(interptime);//Bscale*(b1Pulse->interp1D(interptime))+interpnoise->interp1D(interptime));
        //std::cout<<Bscale*interpnoise->interp1D(interptime)<<std::endl;
        //std::cout<<cos(get_phase(vars.t))<<std::endl;
      return cos(get_phase(vars.t))+Bscale*interpnoise->interp1D(interptime);
        //return 0.0;//+interpnoise->interp1D(interptime));
      }

    else
      //std::cout<<"factor1 "<<vars.t<<"  "<<cos(get_phase(vars.t))<<std::flush<<std::endl;
    return cos(get_phase(vars.t));
}

/**
 * Implementation of get_phase
 */
double BDressingFactor::get_phase(double time) {
	 
	return time*this->wrf + this->phi;
}



// 
//** this is so BDressingCosModFactor can be Floquet friendly. 

int BDressingCosModFactor::findFM() {
// easier to set wrf_amp than fm. no solving needed for wrf_amp!
//double phie = (c * wfm * std::log(std::exp(0.2e1 * k * PI / wfm) + std::exp(k * t1)) - c * wfm * std::log(std::exp(0.2e1 * k * PI / wfm) + std::exp(k * t2)) + std::log(0.1e1 + std::exp(k * t2)) * c * wfm - std::log(0.1e1 + std::exp(k * t1)) * c * wfm + 0.2e1 * PI * b * k) / wfm / k;
double b=wrf; 
double k = 100.; 
this->wrf_amp =  k * (nmod * 0.2e1 * 0.3141592654e1 * fm - b) / (log(exp(k) + exp(k * fm * t1)) - log(exp(k) + exp(k * fm * t2)) - log(0.1e1 + exp(k * fm * t1)) + log(0.1e1 + exp(k * fm * t2)));
std::cout<<"frequency modulation amplitude "<<wrf_amp<<std::endl;
std::cout<<"base frequency "<<b<<std::endl;
return 1;
}

/**
 * Implementation of get_phase
 */// This is now constant power modulation
double BDressingCosModFactor::get_phase(double time) {
   //double phase = wrf*time + phi;
    double b=wrf; 
    double k = 100.;   
    //this->wrf_amp =  k * (6 * 0.2e1 * 0.3141592654e1 * fm - b) / (log(exp(k) + exp(k * fm * t1)) - log(exp(k) + exp(k * fm * t2)) - log(0.1e1 + exp(k * fm * t1)) + log(0.1e1 + exp(k * fm * t2)));
    //std::cout<<b<<"  "<<fm<<"  "<<wrf_amp<<std::endl;
    double t=time;
    //std::cout<<"b "<<b<<std::endl;
    double wfm = 2.*PI*fm;
    double tfloor;
    double tmod=2.*PI/wfm*std::modf(t/(2.*PI/wfm),&tfloor);
    double c=wrf_amp;
    //std::cout<<"why nan? "<<std::exp(fm*t1)<<std::endl;
    //2.*PI*1600.; //amplitude of oscillation
    //double b=wrf;//2.*PI*800.; //base dressing oscillation (minimum)
    //double k = 100*fm;
  	double phitm = (b * k * fm * tmod + c * log(exp(tmod * k * fm) + exp(k * fm * t1)) - c * log(exp(tmod * k * fm) + exp(k * fm * t2)) - c * log(0.1e1 + exp(k * fm * t1)) + c * log(0.1e1 + exp(k * fm * t2))) / k / fm;
    double phie =    (c * log(exp(k) + exp(k * fm * t1)) - c * log(exp(k) + exp(k * fm * t2)) - c * log(0.1e1 + exp(k * fm * t1)) + c * log(0.1e1 + exp(k * fm * t2)) + b * k) / k / fm;
    double phi=  tfloor*phie+phitm;

    //phase += wrf_amp/wfm * sin(wfm * time);
	
	return phi;
}

/// Consant power modulations
double BDressingCosModFactor::getFactor(BFieldVars &vars) {
  double t=vars.t;
  //double wfm = 2.*PI*fm;
  double tfloor;
  double tmod=1./fm*std::modf(t*fm,&tfloor);
  double c=wrf_amp;//2.*PI*1600.; //amplitude of oscillation
  double b=wrf;//2.*PI*800.; //base dressing oscillation (minimum)
  double k = 100.;
  double omega= b+c/(1.+std::exp(-k*fm*(tmod-t1)))-c/(1.+std::exp(-k*fm*(tmod-t2)));
  //std::cout<<"c "<<c<<" b "<<b<<" omega "<<omega<<std::endl;
  double scale = std::sqrt(b/omega);
  double Bfact=BDressingFactor::getFactor(vars);
  //std::cout<<"Scale "<<scale<<" cos(phase) "<<Bfact<<" total "<<scale*Bfact<<std::flush<<std::endl;
  return scale * Bfact;
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
  double scale = 1. + amp *0.5*(1+ cos(wfm * vars.t));
  return scale * BDressingFactor::getFactor(vars);
}

/**
 * Implementation of getFactor SD type 4? I'm going to turn this into robust dressing as well. 
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
  //old robust dressing doesn't really work 
  ///this is the prostitute phase fuction, I use it for whatever crazy idea I have had most recently. 
  //beta approx 1.681
  return dw1/fm*(1.0-cos(fm*time+dw2))+wrf*time+phi; 
  ///return wrf*time+phi;

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
