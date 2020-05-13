#ifndef __BFIELDROOT_H__
#define __BFIELDROOT_H__

//#include "ROOT.h"
#include "TF1.h"
#include "TFormula.h"
#include "BDressing.h"
#include "BField.h"
#include "BFieldInterp.h"

class BFieldROOT {}; // helps rootcint

/**
 * BFactor based on provided TF1 function
 * - the variable x in the TF1 function interpreted as time
 */
class BTF1_Factor : BFactor
{
 public:
  BTF1_Factor(TF1 *fun = 0) : BFactor(), function(fun) {};
  virtual double getFactor(BFieldVars &vars);
  TF1 *function;
 protected:
  
};

/**
 * BFactor based on provided TFormula function
 * - the variables x,y,z,t are used in the TFormula
 */
class BTFormula_Factor : public BFactor {
 public:
  BTFormula_Factor(TFormula *fun = 0) : BFactor(), formula4(fun) {};
  virtual double getFactor(BFieldVars &vars);
  TFormula *formula4;
 protected:
  
};

/**
 * BFactor based on provided TFormula function
 * - the variables x,y,z,t are used in the TFormula
 */
class BTSDFormula_Factor : public BFactor {
 public:
  BTSDFormula_Factor(TFormula *fun = 0, BDressingFactor *Bsdf=0) : BFactor(), formula4(fun), BSDF(Bsdf) {};
  virtual double getFactor(BFieldVars &vars);
  TFormula *formula4;
  BDressingFactor *BSDF;  
 protected:

  
};







#endif
