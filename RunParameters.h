#ifndef __RUN_PARAMETERS_H__
#define __RUN_PARAMETERS_H__

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TFormula.h"

#include "Reading.h"

class RunParameters {
 public:
  RunParameters();
  virtual ~RunParameters();
  void BranchTTree(TTree *tree);
  void SetTTreeBranches(TTree *tree);

  void PrintParameters();
  void LoadParameters();
  void SaveParameters();
  void SaveParVector(Int_t id, Double_t x=0, Double_t y=0, Double_t z=0, const Char_t *description = 0);
  void SaveParVector(Int_t id, Double_t *vals, const Char_t *description = 0);
  int ParseParameterFile(const char *filename);
  int VerifyParameters();

  int check_assignment_length(Reading *param, int args);
  int assignment(Reading *param, Int_t args, Double_t *var);
  int assignment(Reading *param, Int_t args, Int_t *var);
  int assignment(Reading *param, Int_t args, UInt_t *var);
  int assignment(std::string *expression, Int_t args, Char_t *var);
  int field_assignment (const char *field_component, const char *expression);


  // class members
  // --------------------------------------------------
  TTree *tree;
  // storage for the TTree
  Double_t store_par1[3];
  Int_t store_par_id;
  Char_t store_par_info[1024];

  // Parameters stored in TTree:
  Int_t version[3];
  Int_t simType; /** simulation type
		     -1  use simGamma
		     0   neutron
		     1   He3 atom in superfluid He4
		 **/
  Double_t simGamma;
  Double_t simEDM;
  UInt_t seed;
  Int_t numerical_method;
  Double_t max_error;
  Double_t max_angle;
  Double_t totalTime;
  Int_t offsetTime;

  Int_t nBins;
  Int_t randBins;
  /** 0 - no (default)
      1 - yes (`nBins` random bins)
  **/
  // Velocity information
  Double_t speed;
  Double_t temperature;

  // Geometry
  Int_t geometry;
  Double_t radius;

  // Wall diffusivity:
  Double_t diffusion;

  // Gravity (m/s^2)
  Double_t gravity;

  // Box specific
  Double_t boxLow[3];
  Double_t boxHigh[3];

  //Spin Dressing - C. Swank addition
 Double_t SDparam[8];
 Double_t Pulseparam[14];
 Double_t Interpparam[10];
 std::string InterpDir;
 
 Double_t noise;
  // Spin settings
  Double_t spinSettings[3]; /** 0 - Random Spin
			        1 - T1 setup (spins aligned with B0)
			        2 - T2 setup (spins perp to B0)
		            **/

  // Uniform fields
  Double_t B0[3];
  Double_t E0[3];
  Double_t *uniform_gradient;
  Double_t *simple_gradient;
  Double_t *simple_quad_gradient;
  // End of parameters stored in TTree

  Int_t nEntries;
  TFormula *field_formula[3];
  Double_t vdistrexpo;
  int no_normalize;
};



#endif
