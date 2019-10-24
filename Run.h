#ifndef _RUN_H_
#define _RUN_H_

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TF1.h"
//#include "cmsInterpnoiseGen.h"
//#include "BField/cmsB1PulseInterp.h"

#include "ParticleField.h"
#include "BField/BFieldInterp.h"
#include "Boundary.h"
#include "Scattering.h"
// the class used for reading parameter files
#include "Neutron.h"
#include "RunParameters.h"

const int RUN_VERSION[] = {0, 1, 6};

// Globals for keyboard interruption
// TODO: un-globalize
extern int RUN_NEEDS_SAVE;
extern int RUN_EXIT_ASAP;

class RunData {
 public:
  static const int store_n_bounces = 201;

  RunData() : id(0), n(0) {};
  RunData(Int_t n);
  virtual ~RunData();

  void Reset();
  void Initialize(Int_t n);
  void BranchTTree(TTree *tree);
  void SetTTreeBranches(TTree *tree);

  Int_t id;
  Int_t n;
  Int_t *bounces;
  Int_t *scat;
  Int_t *steps;
  Int_t *backsteps;
  Double_t *time;
  Double_t *posx;
  Double_t *posy;
  Double_t *posz;
  Double_t *velx;
  Double_t *vely;
  Double_t *velz;
  Double_t *speed;
  Double_t *spinx;
  Double_t *spiny;
  Double_t *spinz;

  Double_t *fieldx;
  Double_t *fieldy;
  Double_t *fieldz;

  Double_t bx[store_n_bounces];
  Double_t by[store_n_bounces];
  Double_t bz[store_n_bounces];
  Int_t bkind[store_n_bounces];

};

class ParticleFactory {
 public:
  ParticleFactory();

  virtual Neutron *newParticle();
  virtual Double_t *newPosition();
  virtual Double_t *newVelocity();
  virtual Double_t *newSpin();

  Double_t temp_pos[3];
  Double_t temp_vel[3];
  Double_t temp_spin[3];

  RunParameters parameters;
  ParticleField *field;
  Boundary *boundary;
  Scattering *scattering;
  TF1 *vDistribution[4];
};

class Run {
  int version_check;

 public:

  RunData neutronData;
  RunParameters parameters;

  ParticleFactory *factory;

  // simulation name
  Int_t runID;
  char name[200];
  char filename[200];

  Int_t nCurrent;  // the current particle ID

  //TRandom *gRandom; // override the standard gRandom routine

  TF1 *vDistribution[4];
  ParticleField *field;
  Boundary *boundary;
  Scattering *scattering;



  // Mean free path information
  Double_t mfp; // mean free path in meters
  Double_t tau_C; // time between collisions in seconds

  TTree *dataTree;

  TFile *file;
  TRandom *randomPars;


  TTree *createTree(const char* name);

 public:

  Run(Int_t runID);
  virtual ~Run();

  void InitializeRun();
  void InitializeTFile();
  void WriteTFile();

  void SaveNeutronRecord(Neutron *n1, int i);
  Int_t runNeutron();
  Int_t runAll();

  Int_t loadRun();
  int CheckVersion();

  static int FileExists(const char *filename);
  static void parse_argument_line(int argc, char **argv, int *options, char **filename);

 private:
  void useBoltzmannDistribution(Double_t temperature);
};

void signal_hanlde_int(int signal);
void set_signal_handles();

#endif
