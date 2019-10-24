#ifndef __BFIELD_GRADIENT_H__
#define __BFIELD_GRADIENT_H__

#include "BField.h"

/**
 * Field gradient
 * - specified by a 3x3 matrix, flattened to a 9 array
 * - Bgrad_{i} = G_{ij} * x_{j}
 */
class BFieldGradient : public BField {
 public:
  BFieldGradient(double *uniform_gradient);
  virtual int getField(double *out, BFieldVars &vars);
 protected:
  double uniform_gradient[9];
};

/**
 * Simple field gradient
 * - specified by 3 values dBx/dx, dBy/dy, dBz/dz
 */
class BFieldSimpleGradient : public BField {
 public:
  BFieldSimpleGradient(double dBxdx, double dBydy, double dBzdz);
  BFieldSimpleGradient(double *grad_vector);

  virtual int getField(double *out, BFieldVars &vars);
 protected:
  double grad[3];
};

/**
 * Simple quadratic field gradient
 * - specified by 3 values d2Bx/dx2, d2By/dy2, d2Bz/dz2
 */
class BFieldSimpleQuadraticGradient : public BFieldSimpleGradient {
 public:
  BFieldSimpleQuadraticGradient(double a, double b, double c) : BFieldSimpleGradient(a, b, c) {}
  BFieldSimpleQuadraticGradient(double *grad_vector) : BFieldSimpleGradient(grad_vector) {}

  virtual int getField(double *out, BFieldVars &vars);
 protected:
  static int sign_matrix[9];
};

/**
 * Motional field
 * - (E x v)/c^2 motional field where E_field is provided in V/m
 * - velocity vector
 */
class BFieldMotional : public BField {
 public:
  BFieldMotional() : E_field(0) {};
  BFieldMotional(BField *E_field) : E_field(E_field) {};
  virtual ~BFieldMotional();

  virtual int getField(double *out, BFieldVars &vars);
 protected:
  BField *E_field;
};

#endif
