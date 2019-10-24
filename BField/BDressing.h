#ifndef __BDRESSING_H__
#define __BDRESSING_H__
#include "BList.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "cmsInterpnoiseGen.h"
#include "cmsB1PulseInterp.h"
#include "BFieldInterp.h"
class BDressingFactor : public BFactor
{
 private:
 //typedef boost::mt19937 RNG;    // Mersenne Twister
 //typedef boost::normal_distribution<float> DIST_norm;   // Normal Distribution
 //typedef boost::variate_generator<RNG,DIST_norm> NOISE;    // Noise generator. 

     //RNG rng;
     //DIST_norm dist_noise;
     //NOISE noise;
      //: rng(82),dist_noise(0.0f,1.0f),noise(rng,dist_noise) 
        TRandom3 *r;
 public:
      
 BDressingFactor() : BFactor() {this->r = new TRandom3(0);wrf=0;phi=0;};
  double wrf;
  double phi;
  double wrf_amp;
  double fm;
  double amp;
  double dw1;
  double dw2;
  double scale1;
  double scale2;
  double deltat;
  double deltat1;
  double deltat2;
  double phi_mod;
  double scale;
  double noise;
  double w_p;
  double Bscale;
  double T_p;
  double phi_p;
  double width;
  double T_crop;
  double T_pause;
  double rscale;
  int pulse;
  
  cmsInterp* interpnoise;
  cmsB1PulseInterp* b1Pulse;

  virtual double getFactor(BFieldVars &vars);
 //protected:
  virtual double get_phase(double time);
};

class BDressingCosModFactor : public BDressingFactor
{
 public:
 BDressingCosModFactor() : BDressingFactor() {};

  protected:
  virtual double get_phase(double time);
  //virtual double getFactor(BFieldVars &vars);
};

class BDressingCosBModFactor : public BDressingFactor
{
 public:
 BDressingCosBModFactor() : BDressingFactor() {};
  

  virtual double getFactor(BFieldVars &vars);
};

class BDressingPulsedBModFactor : public BDressingFactor
{
 public:
 BDressingPulsedBModFactor() : BDressingFactor() {};
   
  virtual double getFactor(BFieldVars &vars);
  //protected:
};

class BDressingPulsedFreqModFactor : public BDressingFactor
{
 public:
 BDressingPulsedFreqModFactor() : BDressingFactor() {};

  virtual double get_phase(double time);
  //virtual double getFactor(BFieldVars &vars);  
	//protected:
};

class BDressingFuncModFactor : public BDressingFactor
{
 public:
 BDressingFuncModFactor() : BDressingFactor() {};

  protected:
  virtual double get_phase(double time);
  //virtual double getFactor(BFieldVars &vars);
};

class BDressingPulsedFreqModFactor2 : public BDressingFactor
{
 public:
 BDressingPulsedFreqModFactor2() : BDressingFactor() {};

  virtual double get_phase(double time);
  //virtual double getFactor(BFieldVars &vars);  
  //protected:
};

#endif
