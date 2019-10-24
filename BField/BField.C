#include "BField.h"

/**
 * Constructor
 */
BFieldVars::BFieldVars(double t, double *pos, double *vel) : t(t) {
  for (int i=0; i < 3; i++) {
    this->pos[i] = pos[i];
    this->vel[i] = vel[i];
  }
}

/**
 * Constructor
 */
BFieldConst::BFieldConst(double *field) {
  if (field) {
    for (int i=0; i<3; i++)
      this->field_const[i] = field[i];
  }
  else {
    for (int i=0; i<3; i++)
      this->field_const[i] = 0.;
  }
}

/**
 * redefinition of getField
 */
int BFieldConst::getField(double *out, BFieldVars &vars) {
  for (int i=0; i<3; i++)
    out[i] = this->field_const[i];
  return 1;
}

/**
 * dummy
 */
double BFactor::getFactor(BFieldVars &vars) {
  return this->factor;
}
