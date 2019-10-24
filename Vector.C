#include "Vector.h"

#include <iostream>
#include <math.h>
using namespace std;

/**
 * Copy vector \f$ \vec{B} = \vec{A} \f$
 */
void Vector::Copy(Double_t *A, Double_t *B) {
  B[0] = A[0];
  B[1] = A[1];
  B[2] = A[2];
}

/**
 * Add vectors \f$ \vec{C} = \vec{A} + \vec{B} \f$
 */
Double_t *Vector::Add(Double_t *A, Double_t *B, Double_t *C) {

  C[0] = A[0] + B[0];
  C[1] = A[1] + B[1];
  C[2] = A[2] + B[2];

  return C;
}

/**
 * Dot product \f$ \vec{A} \cdot \vec{B} \f$
 */
Double_t Vector::DotProduct(Double_t *A, Double_t *B) {
  if (A && B) {
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
  }
  else
    return 0.0;
}

/**
 * Cross product \f$ \vec{C} = \vec{A} \times \vec{B} \f$
 *
 * it is safe to assign as output either of the inputs
 */
Double_t *Vector::CrossProduct(Double_t *A, Double_t *B, Double_t *C) {
  Double_t temp[3];

  if (A && B && C) {
    temp[0] = A[1]*B[2] - A[2]*B[1];
    temp[1] = A[2]*B[0] - A[0]*B[2];
    temp[2] = A[0]*B[1] - A[1]*B[0];
  }
  else {
    return 0;
  }

  C[0] = temp[0];
  C[1] = temp[1];
  C[2] = temp[2];

  return C;
}

/**
 * Cross product \f$ \vec{A} \times \vec{B} \f$
 *
 * creates output vector
 */
Double_t *Vector::CrossProduct(Double_t *A, Double_t *B) {

  Double_t *C = 0;

  if (A && B)
    C = new Double_t[3];

  return CrossProduct(A,B,C);
}

/**
 * Return magnitude of vector \f$ \left| \vec{A} \right| \f$
 */
Double_t Vector::Norm(Double_t *A) {
  return sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
}

/**
 * Return magnitude squared of vector \f$ \left| \vec{A} \right|^{2} \f$
 */
Double_t Vector::Norm2(Double_t *A) {
  return A[0]*A[0] + A[1]*A[1] + A[2]*A[2];
}

/**
 * Normalize vector \f$ \hat{a} = \frac{\vec{A}}{\left| \vec{A} \right|} \f$
 */
Double_t Vector::Normalize(Double_t *A, Double_t *out) {

  if (!A)
     return -1.0;

  Double_t temp = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
  Double_t temp2 = 1/temp;

  out[0] = A[0] * temp2;
  out[1] = A[1] * temp2;
  out[2] = A[2] * temp2;

  return temp;
}

/**
 * Normalize vector in place
 * \f$ \hat{a} = \frac{\vec{A}}{\left| \vec{A} \right|} \f$
 */
Double_t Vector::Normalize(Double_t *A) {
  return Normalize(A, A);
}

/**
 * Normalize vector in place and scale by magnitude
 * \f$ \vec{A}^{\prime} = m\; \frac{\vec{A}}{\left| \vec{A} \right|} \f$
 */
Double_t Vector::Normalize(Double_t *A, Double_t mag) {
  Double_t temp;

  temp = Normalize(A, A);

  A[0] *= mag;
  A[1] *= mag;
  A[2] *= mag;

  return temp;
}

/**
 * Compute angle between vectors A and B
 */
Double_t Vector::Theta(Double_t *A, Double_t *B) {

  Double_t c;
  Double_t temp[] = {0., 0., 0.};

  c = DotProduct(A, B);
  c /= Norm(A);
  c /= Norm(B);

  CrossProduct(A, B, temp);
  CrossProduct(temp, A, temp);

  if (DotProduct(B, temp) >= 0)
    return acos(c);
  else
    return 2. * PI - acos(c);
}

/**
 * Print vector
 */
void Vector::Print(Double_t *v) {
  cout << "{" << v[0]
       << ", " << v[1]
       << ", " << v[2]
       << "}";
}
