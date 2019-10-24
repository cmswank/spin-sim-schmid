#ifndef __VECTOR_H
#define __VECTOR_H

#include "TROOT.h"

#define PI 3.141592653589793238462643383

class Vector {

 public:
  Vector() {};


  static void Copy(Double_t *A, Double_t *B);

  static Double_t *Add(Double_t *A, Double_t *B, Double_t *C);

  static Double_t DotProduct(Double_t *A, Double_t *B);

  static Double_t *CrossProduct(Double_t *A, Double_t *B, Double_t *C);

  static Double_t *CrossProduct(Double_t *A, Double_t *B);

  static Double_t Norm(Double_t *A);

  static Double_t Norm2(Double_t *A);

  static Double_t Normalize(Double_t *A, Double_t *out);

  static Double_t Normalize(Double_t *A);

  static Double_t Normalize(Double_t *A, Double_t mag);

  static Double_t Theta(Double_t *A, Double_t *B);

  static void Print(Double_t *v);
};

#endif
