#include "BFieldGradient.h"

#include "Vector.h"

/**
 * Constructor
 * @arg grad : a flattened 3x3 matrix (e.g. double[9])
 */
BFieldGradient::BFieldGradient(double *grad) {
  for (int i=0; i<9; i++) {
    uniform_gradient[i] = grad[i];
  }
}

/**
 * implementation of getField()
 * - uses position and gradient matrix to compute gradient field
 */
int BFieldGradient::getField(double *out, BFieldVars &vars) {
  out[0] = Vector::DotProduct(uniform_gradient, vars.pos);
  out[1] = Vector::DotProduct(uniform_gradient + 3, vars.pos);
  out[2] = Vector::DotProduct(uniform_gradient + 6, vars.pos);
  return 1;
}

/**
 * Constructor
 */
BFieldSimpleGradient::BFieldSimpleGradient(double dBxdx, double dBydy, double dBzdz) : BField() {
  grad[0] = dBxdx;
  grad[1] = dBydy;
  grad[2] = dBzdz;
}

/**
 * Constructor
 */
BFieldSimpleGradient::BFieldSimpleGradient(double *grad_vector) : BField() {
  for (int i=0; i<3; i++) {
    grad[i] = grad_vector[i];
  }
}

/**
 * implementation of getField()
 * - uses position and gradient matrix to compute gradient field
 */
int BFieldSimpleGradient::getField(double *out, BFieldVars &vars) {
  for (int i=0; i<3; i++) {
    out[i] = grad[i] * vars.pos[i];
  }
  return 1;
}

int BFieldSimpleQuadraticGradient::sign_matrix[9] = {1,-1,-1,-1,1,-1,-1,-1,1};

/**
 * implementation of getField()
 * - uses position and gradient matrix to compute gradient field
 * - to satisfy divergenless of the field, minimal derived values
 *   from the 3 parameters as QZxx, QYyy, QZzz are found:
 *   QYxy = QZxz = -1/2 QXxx
 *   QXxy = QZyz = -1/2 QYyy
 *   QXxz = QYyz = -1/2 QXxx
 *
 * The minimal format of the quadratic field to preserve divergenceless
 * \f$\vec{\nabla} \cdot \vec{B} = 0\f$ is then:
 * \f[
  \vec{\Delta B_{Q}} =
    \left(
      \begin{array}{c}
        \begin{matrix} x & 0 & 0 \end{matrix} \\
        \begin{matrix} 0 & y & 0 \end{matrix} \\
        \begin{matrix} 0 & 0 & z \end{matrix}
      \end{array}
    \right) \cdot
    \left(
      \begin{array}{c}
        \begin{matrix}  1 & -1 & -1 \end{matrix} \\
        \begin{matrix} -1 &  1 & -1 \end{matrix} \\
        \begin{matrix} -1 & -1 &  1 \end{matrix}
      \end{array}
    \right) \cdot
    \left(
      \begin{array}{c}
        \begin{matrix} x & 0 & 0 \end{matrix} \\
        \begin{matrix} 0 & y & 0 \end{matrix} \\
        \begin{matrix} 0 & 0 & z \end{matrix}
      \end{array}
    \right) \cdot
    \left(
      \begin{array}{c}
        a \\
	b \\
	c
      \end{array}
    \right) =
    \left(
      \begin{array}{c}
        a\,x^{2} - b\,x\,y - c\,x\,z \\
        b\,y^{2} - a\,x\,y - c\,y\,z \\
        c\,z^{2} - a\,x\,z - b\,y\,z
      \end{array}
    \right)
 * \f]
 */
int BFieldSimpleQuadraticGradient::getField(double *out, BFieldVars &vars) {
  double temp[3]; // a*x, b*y, c*z
  for (int i=0; i<3; i++) {
    temp[i] = grad[i] * vars.pos[i];
    out[i] = 0;
  }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      out[i] += sign_matrix[3*i + j] * temp[j];
    }
    out[i] *= vars.pos[i];
  }
  return 1;
}

/**
 * Destructor implementation
 */
BFieldMotional::~BFieldMotional() {
  if (E_field) {
    delete E_field;
  }
}

/**
 * implementation of getField()
 * - uses velocity and stored Electric field
 */
int BFieldMotional::getField(double *out, BFieldVars &vars) {
  if (!E_field) {
    out[0] = out[1] = out[2] = 0;
    return 1;
  }
  double E0[3];
  E_field->getField(E0, vars);
  Vector::CrossProduct(E0, vars.vel, out);

  const double factor_c2 = 1.11265005605362e-17;

  out[0] *= factor_c2;
  out[1] *= factor_c2;
  out[2] *= factor_c2;
  return 1;
}
