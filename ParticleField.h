#ifndef __PARTICLE_FIELD_H__
#define __PARTICLE_FIELD_H__

#include "TROOT.h"

#include "TFormula.h"
#include "BField/cmsB1PulseInterp.h"
#include "BField/BField.h"
#include "BField/BDressing.h"
#include "BField/BFieldGradient.h"
#include "BField/BList.h"
#include "BField/BFieldInterp.h"

/**
 * Implementation of BFieldAdder to deal with
 * particle fields. Provides:
 * - a uniform field B0
 * - motional vxE field from E0
 * - uniform gradient
 * - TFormula based fields
 */
class ParticleField : public BFieldAdder {

  Double_t C = 2.99792458e8;

 protected:
  BFieldConst *uniformField;
  BFieldConst *Efield;
  BFieldMotional *Bv;
  
  BFieldGradient *uniformGrad;

 public:
  ParticleField(Double_t *B0, Double_t *E = 0);
  BDressingFactor *BSDF;	
  BFieldInterp *BFIx;
  BFieldInterp *BFIy;
  BFieldInterp *BFIz;
  void SetUniformGrad(Double_t *grad);
  void SetSimpleGrad(Double_t *grad);
  void SetSimpleQuadGrad(Double_t *grad);


  void add_field_formula(char component, TFormula *formula);
  void add_field_formulaSD(char component, TFormula *formula, BDressingFactor *BSDF);	
  void add_field_interp(char component,BFieldInterp *BFI);
  Double_t *GetE();

};


#endif
