#ifndef __NEUTRON_H
#define __NEUTRON_H

#include "TROOT.h"
#include "TF1.h"

#include "Boundary.h"
#include "ParticleField.h"
#include "BField/BField.h"
#include "Scattering.h"
#include "BField/BFieldInterp.h"

#define NUCLEARMAGNETON 3.1524512326e-14 // MeV / T
#define EV 1.602176487e-19 // J
#define HBAR 1.054571628e-34 // J s

#define NEUTRON_MAGNETICMOMENT = -1.9130427 // relative to nuclear magneton
#define NEUTRON_GAMMA -1.83247172e8 // rad / s / T    OLD: 1.832471850e8 
#define HE3_GAMMA -2.037894585e8  // rad / s / T   (2.037947093e8) old stuff Ezra's numbers come from this.  // rad / s / T

class ParticleStatus {
 public:
  ParticleStatus();
  Int_t valid; // if not valid will not undo

  Double_t position[3]; // 3 component position
  Double_t velocity[3]; // 3-component velocity
  Double_t spin[3]; // 3-component spin (should always be normalized)

  Double_t angle; // angle of precession about main B0 field (rad)
  // (not reliable)
  Int_t bounces; // number of bounces from wall
  Int_t propagations; // number of propagations
  Int_t backsteps; // number of rewinds due to numerical error

  Double_t lifetime; // partial lifetime (seconds)
  Int_t lifetimeSeconds; // number of seconds
};

class LastPropagation {
 public:
  Int_t valid;

  Double_t time;
  Int_t bounce;
  Double_t angle;
  Double_t error;
  Double_t nextTstep;
};

enum AdvanceOption {
  ADVANCE_NO_STOP,
  ADVANCE_STOP_ON_BOUNCE
};

class Neutron {

 public:
  Double_t MIN_T_STEP;

 public:
  Int_t no_normalize;
  // If true, the spin is not normalized at each time step 

 protected:

  ParticleStatus current;

  // B0dir and B0mag created when a field is specified
  Double_t B0dir[3]; // direction of main B field

  Boundary *boundary;

  Scattering *scattering;

  

  void defaultConstructor();

  Double_t gamma; // stores gamma

  
  /*
   * execute simulation with interaction of edm with E field
   */

 public:
  Double_t edm; // edm in e*m
  
  ParticleField *field;
  // previous record kept to enable undo 
  ParticleStatus previous;

  // record of last propagation
  LastPropagation last;

  // what numerical method to use:
  // 0 - Exact for constant B field (DEFAULT)
  // 1 - 4th order RungeKutta
  // 2 - Approximate Exact solutions with error estimation
  // 5 - 5th order RungeKutta with 4th order error estimation
  Int_t kutta;

  Double_t maxError;
  Double_t maxAngle;

  Int_t nScat; // number of scattering (0 for neturons)

  Neutron();
  Neutron(Double_t *position, Double_t *velocity);
  Neutron(Double_t *position, Double_t *velocity, Double_t *spin);
  virtual ~Neutron() {};

  Double_t *GetPosition(Int_t prev = 0);
  Double_t *GetVelocity(Int_t prev = 0);
  Double_t *GetSpin(Int_t prev = 0);

  Int_t GetBounces();
  
  Int_t GetProp();
  Int_t GetBackSteps() {return current.backsteps;}

  Double_t GetLifetime();

  Double_t GetGamma();
  Double_t GetGammaE();

  Double_t GetAngle();

  Int_t GetLastValid();
  Double_t GetLastTime();
  Int_t GetLastBounce();
  Double_t GetLastAngle();
  Double_t GetLastError();

  void SetBoundary(Boundary *b);
  void SetField(ParticleField *f);
  void SetScatter(Scattering *scatter = 0);

  void SetGamma(Double_t gamma) {this->gamma = gamma;}
  void SetEDM(Double_t edm) {this->edm = edm;}

  // Set to be used only for initialization override
  void SetPosition(Double_t *pos);
  void SetVelocity(Double_t *vel);
  void SetSpin(Double_t *spin);
  
  Int_t UndoPropagate();
  void SavePropagate();

  Double_t AddTime(Double_t time);

  virtual Double_t Advance(Double_t time, AdvanceOption bounceOption = ADVANCE_NO_STOP);

  Double_t Propagate(Double_t &step);
  void DoPropagate(Double_t time);

  Double_t PrecessLoop(Double_t time);
  Double_t Precess(Double_t time);

  Double_t PrecessExact(Double_t time, Double_t *B0, Double_t *newSpin);

  Double_t RungeKutta(Double_t time, Double_t *newSpin);
  Double_t RungeKutta5(Double_t time, Double_t *newSpin);

  Double_t PrecessStep(Double_t time, Double_t *newSpin);

  Double_t PrecessApprox(Double_t time, Double_t *newSpin);
  Double_t PrecessApprox(Double_t time, Double_t *newSpin, 
			 Double_t *B1, Double_t *B2);
  Double_t PrecessApprox(Double_t time, 
			 Double_t *spin, Double_t *newSpin, 
			 Double_t *B1, Double_t *B2);

  Double_t CalculateAngle(Double_t *oldSpin, Double_t *newSpin);

};

class Helium3 : public Neutron {
public:
  Helium3();
  Helium3(Double_t *position, Double_t *velocity);
  Helium3(Double_t *position, Double_t *velocity, Double_t *spin);
  virtual Double_t Advance(Double_t time, AdvanceOption bounceOption = ADVANCE_NO_STOP);

  TF1 *vDistribution;
};

#endif
