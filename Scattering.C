#include "Scattering.h"

#include "TRandom.h"

#include <math.h>

Scattering::Scattering(Int_t type) {
 
  next_scat = -1;
  this->type = type;

}

Double_t Scattering::Set(Double_t constant) {
  
  this->constant = constant;

  return constant;

}

Int_t Scattering::next() {

  this->next_scat = get_expo(this->constant);

  return 1;

}

/*
 * Check to see if the particle has scattered.
 * Return value: -1 -> Has not scattered yet
 *                0 -> It is at the scattering point
 *                1 -> It should have scattered already
 *
 */
Int_t Scattering::check(Double_t time) {
 
  if (time < next_scat || next_scat < 0.) 
    return -1;
  else if (time > next_scat) 
    return 1;

  return 0;

}

/* 
 * returns the amount of time to the next scatter
 */
Double_t Scattering::reach(Double_t time) {

  if (time >= next_scat)
    return next_scat;

  return time;
}

/*
 * If possible, advance by time
 */
Double_t Scattering::go(Double_t time) {

  if (check(time) >= 0)
    return -1;

  next_scat -= time;

  return time;

}

/*
 * Returns a random value from an exponential distribution
 * with coefficient tau_c
 */
Double_t Scattering::get_expo(Double_t tau) {

  Double_t temp;
  
  temp = - tau * log(gRandom->Rndm());

  return temp;

}


