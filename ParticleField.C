#include "ParticleField.h"

#include "BField/BFieldROOT.h"

#include "Vector.h"

#include <math.h>
#include <iostream>

/**
 * Constructor.
 * - append the values of B0 and motional Bv field to
 *   vector in BFieldAdder
 */
ParticleField::ParticleField(Double_t *B0, Double_t *E) :
  BFieldAdder(), uniformField(0), Efield(0), Bv(0), uniformGrad(0)
{
  uniformField = new BFieldConst(B0);
  append(uniformField);
  double *notarealdouble =(double*) malloc(20000*sizeof *notarealdouble);
  Efield = new BFieldConst(E);
  Bv = new BFieldMotional(Efield);
  append(Bv);
}

/**
 * Add uniform gradient
 * @arg grad: 9-array (3x3 matrix)
 */
void ParticleField::SetUniformGrad(Double_t *grad) {
  uniformGrad = new BFieldGradient(grad);
  append(uniformGrad);
}

void ParticleField::SetSimpleGrad(Double_t *grad) {
  BFieldSimpleGradient *temp = new BFieldSimpleGradient(grad);
  append(temp);
}

void ParticleField::SetSimpleQuadGrad(Double_t *grad) {
  BFieldSimpleQuadraticGradient *temp = new BFieldSimpleQuadraticGradient(grad);
  append(temp);
}


void ParticleField::add_field_interp(char component, BFieldInterp *BFI) {
  Double_t dir[] = {0, 0, 0};
  switch (component) {
  case 'x':
    dir[0] = 1;
    break;
  case 'y':
    dir[1] = 1;
    break;
  case 'z':
    dir[2] = 1;
    break;
  }
  append(new BFieldScaler(new BFieldConst(dir), BFI));
}


/**
 * Add a single component TFormula field (e.g. x,y,z) C. Swank additon: These should now include spin dressing factor. 
 * That is to say ONLY the TFormula fields are multiplied by the spin dressing factor. -C. Swank
 * This function is only called if the parameter SpinDressing exists in the parameter file (SpinDressing[0]!=0);	 
 */
void ParticleField::add_field_formulaSD(char component, TFormula *formula, BDressingFactor *BSDF) {
  Double_t dir[] = {0, 0, 0};
  switch (component) {
  case 'x':
    dir[0] = 1;
    break;
  case 'y':
    dir[1] = 1;
    break;
  case 'z':
    dir[2] = 1;
    break;
  }
  append(new BFieldScaler(new BFieldConst(dir),
			  new BTSDFormula_Factor(formula,BSDF)));
}


/**
 * Add a single component TFormula field (e.g. x,y,z)
 */

void ParticleField::add_field_formula(char component, TFormula *formula) {
  Double_t dir[] = {0, 0, 0};
  switch (component) {
  case 'x':
    dir[0] = 1;
    break;
  case 'y':
    dir[1] = 1;
    break;
  case 'z':
    dir[2] = 1;
    break;
  }
  append(new BFieldScaler(new BFieldConst(dir),
			  new BTFormula_Factor(formula)));
}

Double_t *ParticleField::GetE() {
  static Double_t ret[3];
  BFieldVars temp_vars;
  Efield->getField(ret, temp_vars);
  return ret;
}
