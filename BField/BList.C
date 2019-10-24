#include "BList.h"

/**
 * Constructor
 */
BFieldAdder::BFieldAdder(BField *A, BField *B) : BField() {
  if (A) {
    fields.push_back(A);
    if (B) {
      fields.push_back(B);
    }
  }
}

/**
 * Destructor deletes all fields stored in vector
 */
BFieldAdder::~BFieldAdder() {
  for (unsigned int i=0; i < fields.size(); i++) {
    if (fields[i]) {
      delete fields[i];
      fields[i] = 0;
    }
  }
}

/**
 * Add a BField to the vector of fields
 */
int BFieldAdder::append(BField *field) {
  this->fields.push_back(field);
  return 1;
}

/**
 * Implementation of getField()
 * - sums all the fields in the vector
 */
int BFieldAdder::getField(double *out, BFieldVars &vars) {
  double sum[] = {0., 0., 0.};
  double temp[3];
  int n = this->fields.size();
  for (int i=0; i<n; i++) {
    fields[i]->getField(temp, vars);
    for (int j=0; j<3; j++)
      sum[j] += temp[j];
  }
  for (int i=0; i<3; i++)
    out[i] = sum[i];
  return 1;
}

/**
 * Destructor implementation, deletes children
 */
BFieldScaler::~BFieldScaler() {
  if (base_field)
    delete base_field;
  if (factor)
    delete factor;
}

/**
 * Implementation of getField()
 * - scales the field by the factor
 */
int BFieldScaler::getField(double *out, BFieldVars &vars) {
  double temp[3];
  base_field->getField(temp, vars);
  double f = factor->getFactor(vars);
  out[0] = f * temp[0];
  out[1] = f * temp[1];
  out[2] = f * temp[2];
  return 1;
}

/**
 * Set factor
 */
void BFieldScaler::SetFactor(BFactor *f) {
  factor = f;
}
