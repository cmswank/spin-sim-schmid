#ifndef _SCATTERING_H_
#define _SCATTERING_H_

#include "TROOT.h"

class Scattering {


 private:
  Double_t next_scat;

  Int_t type; // the type of scattering:
              // 0 -> mean free path
              // 1 -> mean time between collisions

  Double_t constant;


 public:

  Int_t next();

  Double_t Set(Double_t constant);

  Int_t check(Double_t time);

  Double_t reach(Double_t time);

  Double_t go(Double_t time);

  Scattering(Int_t type = 0);

  static Double_t get_expo(Double_t tau);


};


#endif
