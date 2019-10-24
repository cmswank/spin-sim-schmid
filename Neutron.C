#include "Neutron.h"

#include <assert.h>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

#include "TRandom3.h"

#include "Vector.h"

ParticleStatus::ParticleStatus() : valid(0), angle(0), bounces(0),
				   propagations(0), backsteps(0), lifetime(0), lifetimeSeconds(0)
{
  position[0] = position[1] = position[2] = 0;
  velocity[0] = velocity[1] = velocity[2] = 0;
  spin[0] = spin[1] = spin[2] = 0;
}

void Neutron::defaultConstructor() {
  // default constructor invoked at the beginning 
  // by all constructors
  boundary = 0;
  field = 0;
  scattering = 0;
  nScat = 0;

  kutta = 0;
  maxError = 1e-3;
  maxAngle = PI / 8.;

  current.valid = 1;
  previous.valid = 0;
  last.valid = 0;
  last.nextTstep = -1.;

  MIN_T_STEP = 1e-12;
}

Neutron::Neutron() : gamma(NEUTRON_GAMMA) {

  defaultConstructor();

}

Neutron::Neutron(Double_t *position, Double_t *velocity) : gamma(NEUTRON_GAMMA) {

  defaultConstructor();

  this->current.position[0] = position[0];
  this->current.position[1] = position[1];
  this->current.position[2] = position[2];

  this->current.velocity[0]= velocity[0];
  this->current.velocity[1]= velocity[1];
  this->current.velocity[2]= velocity[2];

  this->current.spin[0] = 0.0;
  this->current.spin[1] = 0.0;
  this->current.spin[2] = 0.0;



}

Neutron::Neutron(Double_t *position, Double_t *velocity, Double_t *spin) : gamma(NEUTRON_GAMMA) {

  defaultConstructor();

  this->current.position[0] = position[0];
  this->current.position[1] = position[1];
  this->current.position[2] = position[2];

  this->current.velocity[0]= velocity[0];
  this->current.velocity[1]= velocity[1];
  this->current.velocity[2]= velocity[2];
  
  this->current.spin[0] = spin[0];
  this->current.spin[1] = spin[1];
  this->current.spin[2] = spin[2];


}

Double_t *Neutron::GetPosition(Int_t prev) {

  if (prev)
    return previous.position;
  else
    return current.position;

}



Double_t *Neutron::GetVelocity(Int_t prev) {
  
  if (prev)
    return previous.velocity;
  else
    return current.velocity;

}

Double_t *Neutron::GetSpin(Int_t prev) {
  
  if (prev)
    return previous.spin;
  else
    return current.spin;

}

Int_t Neutron::GetBounces() {

  return current.bounces;

}

Int_t Neutron::GetProp() {

  return current.propagations;

}

 
Double_t Neutron::GetLifetime() {
  // lifetime is calculated by adding the seconds to the 
  // partial lifetime (assures higher precision)
  
  return (Double_t)current.lifetimeSeconds + current.lifetime;

}

Double_t Neutron::GetGamma() {
  return gamma;
}

/**
 * Return the gamma factor in units of m/(V*s)
 * where edm is in units of e*m and the E field
 * in units of V/m
 */
Double_t Neutron::GetGammaE() {
  return 2 * edm * EV / HBAR;
}

Double_t Neutron::GetAngle() {
  
  return current.angle;

}

Int_t Neutron::GetLastValid() {
  return last.valid;
}

Double_t Neutron::GetLastTime() {
  if (last.valid)
    return last.time;
  
  return 0.;
}

Int_t Neutron::GetLastBounce() {
  if (last.valid)
    return last.bounce;

  return 0;
}

Double_t Neutron::GetLastAngle() {
  if (last.valid)
    return last.angle;

  return 0;
}

Double_t Neutron::GetLastError() {
  if (last.valid)
    return last.error;

  return 0;
}

	     

void Neutron::SetBoundary(Boundary *b) {

  boundary = b;

  if (boundary->CheckBoundary(current.position) > 0) {
    cout << "WARNING: SetBoundary, Neutron created outside boundary" << endl;
    cout << "position: (" << current.position[0]
	 << ", " << current.position[1]
	 << ", " << current.position[2] << ")" << endl;
  }
}

void Neutron::SetField (ParticleField *f) {
  field = f;

  Double_t *E0 = field->GetE();
  B0dir[0] = E0[0];
  B0dir[1] = E0[1];
  B0dir[2] = E0[2];
  Vector::Normalize(B0dir);
}

void Neutron::SetScatter(Scattering *scatter) {
  this->scattering = scatter;
}

void Neutron::SetPosition(Double_t *pos) {

  current.position[0] = pos[0];
  current.position[1] = pos[1];
  current.position[2] = pos[2];

}

void Neutron::SetVelocity(Double_t *vel) {

  current.velocity[0] = vel[0];
  current.velocity[1] = vel[1];
  current.velocity[2] = vel[2];

}

void Neutron::SetSpin(Double_t *spin) {

  this->current.spin[0] = spin[0];
  this->current.spin[1] = spin[1];
  this->current.spin[2] = spin[2];

}

Int_t Neutron::UndoPropagate() {
  if (!previous.valid)
    return -1;

  current = previous;
  previous.valid = 0;

  return 0;
}

void Neutron::SavePropagate() {
  previous = current;
}

Double_t Neutron::AddTime(Double_t time) {

  Double_t temp;

  current.lifetime += time;
  if (current.lifetime > 1.) {
    temp = (Double_t)floor(current.lifetime);
    current.lifetimeSeconds += (Int_t)temp;
    current.lifetime -= temp;
  }

  return (Double_t)current.lifetimeSeconds + current.lifetime;

}

/*
 * Advance a neutron for an interval of time.
 * bounceOption : 0 - do not stop
 *                1 - stop on bounce
 *                2 - stop on bounce/scattering
 */
Double_t Neutron::Advance(Double_t time, AdvanceOption bounceOption) {

  Double_t step =  time;
  Double_t total = 0;

  while (total < time) {
    if (step > time - total)
      step = time - total;
    Double_t temp = Propagate(step);

    if (last.bounce == -1)
      return -1.;

    if (temp < 0) { // last propagation did not work
      if (UndoPropagate() < 0) {
	cout << "FATAL undo called when non valid" << endl;
	return -1.;
      }

    }
    else {
      total += temp;
      if (bounceOption == ADVANCE_STOP_ON_BOUNCE &&
	  GetLastBounce() >= 1) {
	return total;
      }
    }

  } // end while

  return total;
}

/**
 * Propagate particle at max time step. Adjust step time to account
 * for boundary and scattering and actually propagate particle. Return
 * the adjusted propagation time step. In addition, modify the
 * parameter step to the next collision.
 */
Double_t Neutron::Propagate(Double_t &step) {

  int reflection = 0;

  if (step + current.lifetime <= current.lifetime) {
    cout << "ERROR : Propagate(step) called with step = " << step
	 << " <= 0.0" << endl;
    return -1.;
  }

  SavePropagate();

  last.valid = 1;

  Double_t prop_time = step;

  // Check for scattering
  if (scattering && scattering->check(prop_time) >= 0) {
    prop_time = scattering->reach(prop_time);
  }

  if (boundary) {
    Int_t check = boundary->CheckBoundary(current.position);
    if (check > 0) {
      cout << "Particle escaped! ("
	   << current.position[0] << ", "
	   << current.position[1] << ", "
	   << current.position[2] << ")" << endl;
      last.bounce = -1;
      return -1.;
    }
    if (last.bounce == 2 && check == 0) {
      // if the last step had a scattering event and the particle within
      // the width of the wall, reflect from the wall, just in case the direction
      // of motion is towards the wall
      int check_bounce = boundary->Reflect(current.position, current.velocity);
      if (check_bounce < 0) {
	cout << "Error bouncing right after scattering" << endl;
	last.bounce = -1;
	return -1;
      }
    }
    // check the move step
    Int_t wall_collision = boundary->CheckMove(current.position, current.velocity, prop_time);
    if (wall_collision == 0) {
      reflection = 1;
    }
    else if (wall_collision > 0) {
      reflection = 1;
      Double_t tempTime = boundary->ReachBoundary(current.position, current.velocity);
      if (tempTime < 0) {
	tempTime = 0;
      }
      if (tempTime < prop_time) {
	prop_time = tempTime;
      }
    }
  }

  // do the actual propagation+precession
  if (prop_time > 0.) {
    DoPropagate(prop_time);
  }

  // check for bounces
  last.bounce = 0;
  Int_t bounces = 0;
  if (reflection) {
    bounces = boundary->Reflect(current.position, current.velocity);
    if (bounces < 0) {
      cout << "Error Bouncing " << endl;
      last.bounce = -1;
      return -1.;
    }
    else if (bounces > 0) {
      current.bounces += bounces;
      last.bounce = 1;
      // compute next step time based on boundary
      step = boundary->ReachBoundary(current.position, current.velocity);
    }
  }

  // Check if scattered
  if (scattering && scattering->check(prop_time) == 0) {
    last.bounce = 2;
    scattering->next();
  }
  else if (scattering && scattering->go(prop_time) < 0) {
    cout << "Error Scattering " << endl;
    last.bounce = -1;
    return -1.;
  }
  

  current.propagations++;

  last.time = prop_time;

  last.valid = 1;

  return prop_time;
}

void Neutron::DoPropagate(Double_t time) {
  Double_t tempPos[3];
  Double_t tempVel[3];
  Double_t save_lifetime = current.lifetime;
  Int_t save_lifetime_seconds = current.lifetimeSeconds;
  Double_t temp;

  // Save the position, velocity to which the particle will propagate
  Vector::Copy(current.position, tempPos);
  Vector::Copy(current.velocity, tempVel);
  if (boundary) {
    boundary->Move(tempPos, tempVel, time);
  }
  // Precession
  if (field) {
    temp = PrecessLoop(time);
    assert(temp >= time);
  }

  current.lifetime = save_lifetime;
  current.lifetimeSeconds = save_lifetime_seconds;

  // calculate the new lifetime
  current.lifetime += time;
  if (current.lifetime > 1.) {
    temp = (Double_t)floor(current.lifetime);
    current.lifetimeSeconds += (Int_t)temp;
    current.lifetime -= temp;
  }

  // Now restore the new position and velocity
  Vector::Copy(tempPos, current.position);
  Vector::Copy(tempVel, current.velocity);
}

/*
 * Loop through precession of variable stepsize.
 * Effectively separating propagation from spin precession
 */
Double_t Neutron::PrecessLoop(Double_t time) {

  Double_t t_step = time,
    tempTime = 0.;

  Double_t tempPos[3], tempVel[3];

  // save original position and velocity
  Vector::Copy(current.position, tempPos);
  Vector::Copy(current.velocity, tempVel);

  Int_t steps = 0;

  Int_t i = 0;

  
  while (tempTime < time) {

    if ((kutta == 1 || kutta == 5) && last.nextTstep > MIN_T_STEP)
      t_step = last.nextTstep;

    if (t_step > time - tempTime) t_step = time - tempTime;

    if (t_step < MIN_T_STEP) t_step = MIN_T_STEP;

    if (Precess(t_step) < 0) {
      if (t_step <= last.nextTstep) {
	cout << last.error << " > error_max "
	     << "time: " << t_step
	     << " >> " << last.nextTstep << endl;
      }
      current.backsteps++;
      i++;
    }
    else {
      i=0;
      tempTime += t_step;
      AddTime(t_step);
      // calculate position, velocity from original
      if (boundary) {
	Vector::Copy(tempPos, current.position);
	Vector::Copy(tempVel, current.velocity);
	boundary->Move(current.position, current.velocity, tempTime);
      }
      else {
	current.position[0] = tempPos[0] + current.velocity[0] * tempTime;
	current.position[1] = tempPos[1] + current.velocity[1] * tempTime;
	current.position[2] = tempPos[2] + current.velocity[2] * tempTime;
      }
      steps++;
    }

    if (i > 100) {
      cout << "Precession change size > 100";
      cout << " t_step = " << t_step;
      cout << " last.nextTstep = " << last.nextTstep << endl;
      return -1;
    }
    
  } // end while loop

  if (steps > 0) 
    current.propagations += steps - 1;

  return tempTime;

}

Double_t Neutron::Precess(Double_t time) {

  Double_t newSpin[3];

  Double_t tempB0[3];

  Double_t tempTime;

  tempTime = time;

  if (kutta == 1) {
    last.angle = RungeKutta(time, newSpin);
    if (last.error > maxError)
      tempTime = -1;
  }
  else if (kutta == 5) {
    tempTime = RungeKutta5(time, newSpin);
    if (last.error > maxError)
      tempTime = -1;
  }
  else if (kutta == 0) {
    BFieldVars args(GetLifetime(), current.position, current.velocity);
    field->getField(tempB0, args);
    tempTime = PrecessExact(time,tempB0, current.spin);
    this->current.angle += last.angle;
    return tempTime;
  }
  else if (kutta == 2) 
    last.angle = PrecessApprox(time, newSpin);
  else if (kutta == 6)
    tempTime = PrecessStep(time, newSpin);
  
  // If the precession was not satisfactory, do not store the spin
  if (tempTime < 0)
    return tempTime;
  
  this->current.angle += last.angle;
  
  current.spin[0] = newSpin[0];
  current.spin[1] = newSpin[1];
  current.spin[2] = newSpin[2];
  
  return tempTime;
}

Double_t Neutron::PrecessExact(Double_t time, Double_t *B0, Double_t *spin) {

  Double_t spinold[3], cross[] = {0., 0., 0.};
  Double_t sb, rotation;
  Double_t bmag, bmag2, c, s, A;

  spinold[0] = spin[0];
  spinold[1] = spin[1];
  spinold[2] = spin[2];

  Double_t gamma = GetGamma();
  
  Vector::CrossProduct(spin, B0, cross);
  
  sb = Vector::DotProduct(spin, B0);

  bmag = Vector::Norm(B0);
  bmag2 = bmag * bmag;

  rotation = gamma * bmag * time;

  c = cos(rotation);
  s = sin(rotation);

  A = sb / bmag2;

  spin[0] = B0[0] * A + (spin[0] - B0[0] * A) * c + (cross[0] / bmag) * s;
  spin[1] = B0[1] * A + (spin[1] - B0[1] * A) * c + (cross[1] / bmag) * s;
  spin[2] = B0[2] * A + (spin[2] - B0[2] * A) * c + (cross[2] / bmag) * s;

  if (!this->no_normalize) {
    Vector::Normalize(spin);
  }
  
  rotation = floor(rotation / (2. * PI)) * 2. * PI;
  
  rotation += CalculateAngle(spinold, spin);

  last.angle = rotation;

  last.error = 0.;

  return time;

}



Double_t Neutron::RungeKutta(Double_t time, Double_t *newSpin) {
  
  Double_t tempB0[3], k1[3], k2[3], k3[3], k4[3];
  Double_t *Efield = field->GetE(), add_edm[] = {0., 0., 0.};

  Double_t tempSpin[3], tempPos[3], tempVel[3];
  
  Double_t angle = 1., rotations, partialAngle, tempAngle;

  Int_t use_edm = (edm != 0)? 1 : 0;
  Double_t gamma_edm = GetGammaE();

  Double_t gamma = GetGamma();
  
  // Get k1:
  BFieldVars vars1(GetLifetime(), current.position, current.velocity);
  field->getField(tempB0, vars1);

  Vector::CrossProduct(current.spin, tempB0, k1);
  
  k1[0] *= gamma;
  k1[1] *= gamma;
  k1[2] *= gamma;
  if (use_edm) {
    Efield = field->GetE();
    Vector::CrossProduct(current.spin, Efield, add_edm);
    k1[0] += add_edm[0] * gamma_edm;
    k1[1] += add_edm[1] * gamma_edm;
    k1[2] += add_edm[2] * gamma_edm;
  }

  // Get k2:
  tempSpin[0] = current.spin[0] + (time / 2.) * k1[0];
  tempSpin[1] = current.spin[1] + (time / 2.) * k1[1];
  tempSpin[2] = current.spin[2] + (time / 2.) * k1[2];

  Vector::Copy(current.position, tempPos);
  Vector::Copy(current.velocity, tempVel);
  if (boundary) {
    boundary->Move(tempPos, tempVel, time/2.);
  }
  else {
    tempPos[0] += tempVel[0] * time / 2.;
    tempPos[1] += tempVel[1] * time / 2.;
    tempPos[2] += tempVel[2] * time / 2.;
  }

  BFieldVars vars2(GetLifetime() + time/2., current.position, current.velocity);
  field->getField(tempB0, vars2);
  Vector::CrossProduct(tempSpin, tempB0, k2);

  k2[0] *= gamma;
  k2[1] *= gamma;
  k2[2] *= gamma;
  if (use_edm) {
    Efield = field->GetE();
    Vector::CrossProduct(tempSpin, Efield, add_edm);
    k2[0] += add_edm[0] * gamma_edm;
    k2[1] += add_edm[1] * gamma_edm;
    k2[2] += add_edm[2] * gamma_edm;
  }
  
  
  // Get k3:
  tempSpin[0] = current.spin[0] + (time / 2.) * k2[0];
  tempSpin[1] = current.spin[1] + (time / 2.) * k2[1];
  tempSpin[2] = current.spin[2] + (time / 2.) * k2[2];
  
  Vector::CrossProduct(tempSpin, tempB0, k3);

  k3[0] *= gamma;
  k3[1] *= gamma;
  k3[2] *= gamma;
  if (use_edm) {
    Efield = field->GetE();
    Vector::CrossProduct(tempSpin, Efield, add_edm);
    k3[0] += add_edm[0] * gamma_edm;
    k3[1] += add_edm[1] * gamma_edm;
    k3[2] += add_edm[2] * gamma_edm;
  }

  // Get k4:
  tempSpin[0] = current.spin[0] + time * k3[0];
  tempSpin[1] = current.spin[1] + time * k3[1];
  tempSpin[2] = current.spin[2] + time * k3[2];

  Vector::Copy(current.position, tempPos);
  Vector::Copy(current.velocity, tempVel);
  if (boundary) {
    boundary->Move(tempPos, tempVel, time);
  }
  else {
    tempPos[0] += tempVel[0] * time;
    tempPos[1] += tempVel[1] * time;
    tempPos[2] += tempVel[2] * time;
  }

  BFieldVars vars4(GetLifetime() + time, tempPos, tempVel);
  field->getField(tempB0, vars4);
  Vector::CrossProduct(tempSpin, tempB0, k4);

  k4[0] *= gamma;
  k4[1] *= gamma;
  k4[2] *= gamma;
  if (use_edm) {
    Efield = field->GetE();
    Vector::CrossProduct(tempSpin, Efield, add_edm);
    k4[0] += add_edm[0] * gamma_edm;
    k4[1] += add_edm[1] * gamma_edm;
    k4[2] += add_edm[2] * gamma_edm;
  }
  
  newSpin[0] = current.spin[0] + (time / 6.) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
  newSpin[1] = current.spin[1] + (time / 6.) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
  newSpin[2] = current.spin[2] + (time / 6.) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
  
  // older
  angle = time / 6. * (Vector::Norm(k1) + 2. * Vector::Norm(k2)
		     + 2 * Vector::Norm(k3) + Vector::Norm(k4));
  
  if (angle < 0.) 
    angle = - angle;
  
  if (!this->no_normalize) {
    Vector::Normalize(newSpin);
  }
  
  rotations = floor(angle / (2.0 * PI));
  
  partialAngle = CalculateAngle(current.spin, newSpin);
  
  tempAngle = rotations * 2.0 * PI + partialAngle;

  if (angle >= tempAngle) {
    if (tempAngle + PI  < angle)
      tempAngle = tempAngle + 2.0 * PI;
  }
  else {
    if (angle  < tempAngle - PI)
      tempAngle = tempAngle - 2.0 * PI;
  }
  // See if the error is above the threshold
  //if (tempAngle >= limitAngle || angle >= limitAngle)
  //  last.error = 1;

  if ((angle = gamma * Vector::Norm(tempB0) * time) < 0.) angle = -angle;
  

  if (angle > maxAngle) {
    last.error = 1;
    last.nextTstep = time * maxAngle / angle;
  }
  else {
    last.error = 0;
    last.nextTstep = time * maxAngle / angle;
  }

  return tempAngle;

}

Double_t Neutron::RungeKutta5(Double_t time, Double_t *newSpin) {

  // Cash-Karp Parameters for Embedded Runge-Kutta Method
  // From Numerical Recipes in C++
  const Double_t a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875;
  
  const Double_t b21 = 0.2;
  const Double_t b31 = 0.075, b32 = 0.225;
  const Double_t b41 = 0.3, b42 = -0.9, b43 = 1.2;
  const Double_t b51 = -2.03703703703703692e-1,
    b52 = 2.5, b53 = -2.59259259259259256,
    b54 = 1.29629629629629628;
  const Double_t b61 = 2.94958043981481469e-2,
    b62 = 3.41796875e-1,
    b63 = 4.15943287037037063e-2,
    b64 = 4.00345413773148140e-1,
    b65 = 6.1767578125e-2;
  const Double_t 
    c1 = 9.78835978835978782e-2,
    c2 = 0.,
    c3 = 4.02576489533011284e-1,
    c4 = 2.10437710437710451e-1,
    c5 = 0.,
    c6 = 2.89102202145680387e-1;
  const Double_t 
    d1 = 1.02177372685185189e-1,
    d2 = 0.,
    d3 = 3.83907903439153431e-1,
    d4 = 2.44592737268518517e-1,
    d5 = 1.93219866071428562e-2,
    d6 = 0.25;
  
  static Double_t gamma = GetGamma();

  Int_t use_edm = (edm != 0)? 1 : 0;
  static Double_t gamma_edm = GetGammaE();
  
  Double_t temp;

  Double_t tempSpin[3], tempB0[3], tempPos[3], tempVel[3], tempErr[3], tempNewSpin[3];
  Double_t tempTime, tempAngle;

  Double_t *Efield = field->GetE(), add_edm[] = {0., 0., 0.};

  Double_t Bmag[6];

  Double_t k1[3], k2[3], k3[3], k4[3], k5[3], k6[3];

  // Get k1
  BFieldVars vars1(GetLifetime(), current.position, current.velocity);
  field->getField(tempB0, vars1);

  Vector::CrossProduct(current.spin, tempB0, k1);
  k1[0] *= time * gamma;
  k1[1] *= time * gamma;
  k1[2] *= time * gamma;
  if (use_edm) {
    Efield = field->GetE();
    Vector::CrossProduct(current.spin, Efield, add_edm);
    k1[0] += add_edm[0] * time * gamma_edm;
    k1[1] += add_edm[1] * time * gamma_edm;
    k1[2] += add_edm[2] * time * gamma_edm;
  }

  Bmag[0] = Vector::Norm(tempB0);

  // Get k2
  tempTime = a2 * time;
  Vector::Copy(current.position, tempPos);
  Vector::Copy(current.velocity, tempVel);
  if (boundary) {
    boundary->Move(tempPos, tempVel, tempTime);
  }
  else {
    tempPos[0] += tempTime * tempVel[0];
    tempPos[0] += tempTime * tempVel[0];
    tempPos[0] += tempTime * tempVel[0];
  }
  BFieldVars vars2(GetLifetime() + tempTime, tempPos, tempVel);
  field->getField(tempB0, vars2);

  tempSpin[0] = current.spin[0] + b21 * k1[0];
  tempSpin[1] = current.spin[1] + b21 * k1[1];
  tempSpin[2] = current.spin[2] + b21 * k1[2];

  Vector::CrossProduct(tempSpin, tempB0, k2);
  k2[0] *= time * gamma;
  k2[1] *= time * gamma;
  k2[2] *= time * gamma;
  if (use_edm) {
    Efield = field->GetE();
    Vector::CrossProduct(current.spin, Efield, add_edm);
    k2[0] += add_edm[0] * time * gamma_edm;
    k2[1] += add_edm[1] * time * gamma_edm;
    k2[2] += add_edm[2] * time * gamma_edm;
  }

  Bmag[1] = Vector::Norm(tempB0);

  // Get k3
  tempTime = a3 * time;
  Vector::Copy(current.position, tempPos);
  Vector::Copy(current.velocity, tempVel);
  if (boundary) {
    boundary->Move(tempPos, tempVel, tempTime);
  }
  else {
    tempPos[0] += tempTime * tempVel[0];
    tempPos[0] += tempTime * tempVel[0];
    tempPos[0] += tempTime * tempVel[0];
  }
  BFieldVars vars3(GetLifetime() + tempTime, tempPos, tempVel);
  field->getField(tempB0, vars3);

  tempSpin[0] = current.spin[0] + b31 * k1[0] + b32 * k2[0];
  tempSpin[1] = current.spin[1] + b31 * k1[1] + b32 * k2[1];
  tempSpin[2] = current.spin[2] + b31 * k1[2] + b32 * k2[2];

  Vector::CrossProduct(tempSpin, tempB0, k3);
  k3[0] *= time * gamma;
  k3[1] *= time * gamma;
  k3[2] *= time * gamma;
  if (use_edm) {
    Efield = field->GetE();
    Vector::CrossProduct(current.spin, Efield, add_edm);
    k3[0] += add_edm[0] * time * gamma_edm;
    k3[1] += add_edm[1] * time * gamma_edm;
    k3[2] += add_edm[2] * time * gamma_edm;
  }

  Bmag[2] = Vector::Norm(tempB0);

  // Get k4
  tempTime = a4 * time;
  Vector::Copy(current.position, tempPos);
  Vector::Copy(current.velocity, tempVel);
  if (boundary) {
    boundary->Move(tempPos, tempVel, tempTime);
  }
  else {
    tempPos[0] += tempTime * tempVel[0];
    tempPos[0] += tempTime * tempVel[0];
    tempPos[0] += tempTime * tempVel[0];
  }
  BFieldVars vars4(GetLifetime() + tempTime, tempPos, tempVel);
  field->getField(tempB0, vars4);

  tempSpin[0] = current.spin[0] + b41 * k1[0] + b42 * k2[0] + b43 * k3[0];
  tempSpin[1] = current.spin[1] + b41 * k1[1] + b42 * k2[1] + b43 * k3[1];
  tempSpin[2] = current.spin[2] + b41 * k1[2] + b42 * k2[2] + b43 * k3[2];

  Vector::CrossProduct(tempSpin, tempB0, k4);
  k4[0] *= time * gamma;
  k4[1] *= time * gamma;
  k4[2] *= time * gamma;
  if (use_edm) {
    Efield = field->GetE();
    Vector::CrossProduct(current.spin, Efield, add_edm);
    k4[0] += add_edm[0] * time * gamma_edm;
    k4[1] += add_edm[1] * time * gamma_edm;
    k4[2] += add_edm[2] * time * gamma_edm;
  }

  Bmag[3] = Vector::Norm(tempB0);

  // Get k5
  tempTime = a5 * time;
  Vector::Copy(current.position, tempPos);
  Vector::Copy(current.velocity, tempVel);
  if (boundary) {
    boundary->Move(tempPos, tempVel, tempTime);
  }
  else {
    tempPos[0] += tempTime * tempVel[0];
    tempPos[0] += tempTime * tempVel[0];
    tempPos[0] += tempTime * tempVel[0];
  }
  BFieldVars vars5(GetLifetime() + tempTime, tempPos, tempVel);
  field->getField(tempB0, vars5);

  tempSpin[0] = current.spin[0] + b51 * k1[0] + b52 * k2[0] + b53 * k3[0] + b54 * k4[0];
  tempSpin[1] = current.spin[1] + b51 * k1[1] + b52 * k2[1] + b53 * k3[1] + b54 * k4[1];
  tempSpin[2] = current.spin[2] + b51 * k1[2] + b52 * k2[2] + b53 * k3[2] + b54 * k4[2];

  Vector::CrossProduct(tempSpin, tempB0, k5);
  k5[0] *= time * gamma;
  k5[1] *= time * gamma;
  k5[2] *= time * gamma;
  if (use_edm) {
    Efield = field->GetE();
    Vector::CrossProduct(current.spin, Efield, add_edm);
    k5[0] += add_edm[0] * time * gamma_edm;
    k5[1] += add_edm[1] * time * gamma_edm;
    k5[2] += add_edm[2] * time * gamma_edm;
  }

  Bmag[4] = Vector::Norm(tempB0);

  // Get k6
  tempTime = a6 * time;
  Vector::Copy(current.position, tempPos);
  Vector::Copy(current.velocity, tempVel);
  if (boundary) {
    boundary->Move(tempPos, tempVel, tempTime);
  }
  else {
    tempPos[0] += tempTime * tempVel[0];
    tempPos[0] += tempTime * tempVel[0];
    tempPos[0] += tempTime * tempVel[0];
  }
  BFieldVars vars6(GetLifetime() + tempTime, tempPos, tempVel);
  field->getField(tempB0, vars6);

  tempSpin[0] = current.spin[0] + b61*k1[0] + b62*k2[0] + b63*k3[0] + b64*k4[0] + b65*k5[0];
  tempSpin[1] = current.spin[1] + b61*k1[1] + b62*k2[1] + b63*k3[1] + b64*k4[1] + b65*k5[1];
  tempSpin[2] = current.spin[2] + b61*k1[2] + b62*k2[2] + b63*k3[2] + b64*k4[2] + b65*k5[2];

  Vector::CrossProduct(tempSpin, tempB0, k6);
  k6[0] *= time * gamma;
  k6[1] *= time * gamma;
  k6[2] *= time * gamma;
  if (use_edm) {
    Efield = field->GetE();
    Vector::CrossProduct(current.spin, Efield, add_edm);
    k6[0] += add_edm[0] * time * gamma_edm;
    k6[1] += add_edm[1] * time * gamma_edm;
    k6[2] += add_edm[2] * time * gamma_edm;
  }

  Bmag[5] = Vector::Norm(tempB0);

  // calculate spins:
  // N.B. no normalization occurs at this stage so that the error can be calculated
  // in a more consistent way

  tempNewSpin[0] = current.spin[0] + c1*k1[0] + c2*k2[0] + c3*k3[0] + c4*k4[0] + c5*k5[0] + c6*k6[0];
  tempNewSpin[1] = current.spin[1] + c1*k1[1] + c2*k2[1] + c3*k3[1] + c4*k4[1] + c5*k5[1] + c6*k6[1];
  tempNewSpin[2] = current.spin[2] + c1*k1[2] + c2*k2[2] + c3*k3[2] + c4*k4[2] + c5*k5[2] + c6*k6[2];
  
  tempSpin[0] = current.spin[0] + d1*k1[0] + d2*k2[0] + d3*k3[0] + d4*k4[0] + d5*k5[0] + d6*k6[0];
  tempSpin[1] = current.spin[1] + d1*k1[1] + d2*k2[1] + d3*k3[1] + d4*k4[1] + d5*k5[1] + d6*k6[1];
  tempSpin[2] = current.spin[2] + d1*k1[2] + d2*k2[2] + d3*k3[2] + d4*k4[2] + d5*k5[2] + d6*k6[2];
  
  // error estimations
  tempErr[0] = tempNewSpin[0] - tempSpin[0];
  tempErr[1] = tempNewSpin[1] - tempSpin[1];
  tempErr[2] = tempNewSpin[2] - tempSpin[2];
  
  last.error = Vector::Norm(tempErr);
  
  temp = last.error / maxError;
  if (temp > 1.0) {
    last.nextTstep = 0.9 * time * pow(temp, -0.25);
    if (last.nextTstep < time * 0.1)
      last.nextTstep = time * 0.1;
  }
  else if (temp > 1.89e-4) {
    last.nextTstep = 0.9 * time * pow(temp, -0.2);
    if (last.nextTstep > time * 5.0)
      last.nextTstep = time * 5.0;
  }
  else {
    last.nextTstep = time * 5.0;
  }

  // Calculate the rotation angle
  tempAngle = time * gamma * (c1*Bmag[0] + c2*Bmag[1] + c3*Bmag[2] + c4*Bmag[3] + c5*Bmag[4] + c6*Bmag[5]);
  if (tempAngle < 0.) tempAngle = -tempAngle;

  if (tempAngle > maxAngle) {
    last.error = 1.;
    last.nextTstep = 0.9 * time * maxAngle / tempAngle;
  }
  else if (tempAngle * last.nextTstep / time > maxAngle) {
    // make sure a larger time will not interfere with angle
    last.nextTstep = 0.9 * time * maxAngle / tempAngle;
  }

  if (!this->no_normalize) {
    // necessary but can be done outside of rk
    Vector::Normalize(tempNewSpin);
  }

  newSpin[0] = tempNewSpin[0];
  newSpin[1] = tempNewSpin[1];
  newSpin[2] = tempNewSpin[2];

  last.angle = tempAngle;

  return time;

}

Double_t Neutron::PrecessStep(Double_t time, Double_t *newSpin) {
  Double_t B[3];
  BFieldVars vars(GetLifetime(), current.position, current.velocity);
  field->getField(B, vars);

  Double_t dSpin[3], gamma = GetGamma();

  Vector::CrossProduct(current.spin, B, dSpin);
  dSpin[0] = gamma * dSpin[0];
  dSpin[1] = gamma * dSpin[1];
  dSpin[2] = gamma * dSpin[2];

  newSpin[0] += time * dSpin[0];
  newSpin[1] += time * dSpin[1];
  newSpin[2] += time * dSpin[2];
  return 0.;
}

Double_t Neutron::PrecessApprox(Double_t time, Double_t *newSpin) {

  Double_t B1[3], B2[3];
  Double_t tempPos[3], tempVel[3];

  BFieldVars vars1(GetLifetime(), current.position, current.velocity);
  field->getField(B1, vars1);
  Vector::Copy(current.position, tempPos);
  Vector::Copy(current.velocity, tempVel);
  if (boundary) {
    boundary->Move(tempPos, tempVel, time);
  }
  else {
    tempPos[0] += tempVel[0] * time;
    tempPos[1] += tempVel[1] * time;
    tempPos[2] += tempVel[2] * time;
  }
  BFieldVars vars2(GetLifetime(), tempPos, tempVel);
  field->getField(B2, vars2);

  return PrecessApprox(time, newSpin, B1, B2);
}


Double_t Neutron::PrecessApprox(Double_t time, Double_t *newSpin, 
				Double_t *B1, Double_t *B2) {

  return PrecessApprox(time, current.spin, newSpin, B1, B2);

}

Double_t Neutron::PrecessApprox(Double_t time, 
				Double_t *spin, Double_t *newSpin, 
				Double_t *B1, Double_t *B2) {

  Double_t Bhalf[3];

  Double_t tempSpin[3], spin1[3], spin2[3];
  Double_t angle2, partialAngle, rotations;

  Double_t tempAngle;

  Bhalf[0] = (B1[0] + B2[0]) / 2.0;
  Bhalf[1] = (B1[1] + B2[1]) / 2.0;
  Bhalf[2] = (B1[2] + B2[2]) / 2.0;

  spin1[0] = tempSpin[0] = spin[0];
  spin1[1] = tempSpin[1] = spin[1];
  spin1[2] = tempSpin[2] = spin[2];
  
  PrecessExact(time, B1, spin1);
  PrecessExact(time, B2, tempSpin);
  
  spin1[0] = (spin1[0] + tempSpin[0]) / 2.0;
  spin1[1] = (spin1[1] + tempSpin[1]) / 2.0;
  spin1[2] = (spin1[2] + tempSpin[2]) / 2.0;
  
  if (!this->no_normalize) {
    Vector::Normalize(spin1);
  }

  spin2[0] = tempSpin[0] = spin[0];
  spin2[1] = tempSpin[1] = spin[1];
  spin2[2] = tempSpin[2] = spin[2];
  
  angle2 = PrecessExact(time / 2.0, B1, spin2);
  angle2 += PrecessExact(time / 2.0, Bhalf, tempSpin);
  
  spin2[0] = (spin2[0] + tempSpin[0]) / 2.0;
  spin2[1] = (spin2[1] + tempSpin[1]) / 2.0;
  spin2[2] = (spin2[2] + tempSpin[2]) / 2.0;
  
  if (!this->no_normalize) {
    Vector::Normalize(spin2);
  }

  tempSpin[0] = spin2[0];
  tempSpin[1] = spin2[1];
  tempSpin[2] = spin2[2];

  angle2 += PrecessExact(time / 2.0, Bhalf, spin2);
  angle2 += PrecessExact(time / 2.0, B2, tempSpin);

  angle2 /= 2.0;
  
  newSpin[0] = (spin2[0] + tempSpin[0]) / 2.0;
  newSpin[1] = (spin2[1] + tempSpin[1]) / 2.0;
  newSpin[2] = (spin2[2] + tempSpin[2]) / 2.0;

  if (!this->no_normalize) {
    Vector::Normalize(newSpin);
  }

  tempSpin[0] = spin1[0] - newSpin[0];
  tempSpin[1] = spin1[1] - newSpin[1];
  tempSpin[2] = spin1[2] - newSpin[2];

  last.error = Vector::Norm(tempSpin);
  //cout << " error = " << last.error;
  

  rotations = floor(angle2 / (2.0 * PI));
  
  partialAngle = CalculateAngle(spin, spin2);
  
  tempAngle = rotations * 2.0 * PI + partialAngle;

  if (angle2 >= tempAngle) {
    if (tempAngle + PI  < angle2)
      tempAngle = tempAngle + 2.0 * PI;
  }
  else {
    if (angle2  < tempAngle - PI)
      tempAngle = tempAngle - 2.0 * PI;
  }

  //cout << " (angle = " << tempAngle << ")" << endl;

  return tempAngle;


  }

Double_t Neutron::CalculateAngle(Double_t *oldSpin, Double_t *newSpin) {
  
  Double_t temp1[] = {0., 0., 0.};
  Double_t temp2[] = {0., 0., 0.};
  Double_t c, s;

  // B0 dir is not the right direction: the right direction
  // for inhomogeneous B field would be the field itself.
  // Better if at the end of the propagation and the precession
  Vector::CrossProduct(B0dir, oldSpin, temp1);
  Vector::Normalize(temp1);

  Vector::CrossProduct(B0dir, newSpin, temp2);
  Vector::Normalize(temp2);
  
  c = Vector::DotProduct(temp1, temp2);
  
  Vector::CrossProduct(temp1, temp2, temp1);
  s = Vector::DotProduct(temp1, B0dir);

  if (s > 0.)
    return acos(c);
  else
    return 2. * PI - acos(c);
  
  

}

Helium3::Helium3() : Neutron::Neutron(), vDistribution(0) {
  gamma = HE3_GAMMA;
}

Helium3::Helium3(Double_t *position, Double_t *velocity) : Neutron::Neutron(position, velocity), vDistribution(0) {
  gamma = HE3_GAMMA;
}

Helium3::Helium3(Double_t *position, Double_t *velocity, Double_t *spin) : Neutron::Neutron(position, velocity, spin), vDistribution(0) {
  gamma = HE3_GAMMA;
}

/*
 * Implementation of Advance in He3
 *  -- adds thermalization and scattering
 * Advance a neutron for an interval of time.
 * bounceOption : 0 - do not stop
 *                1 - stop on bounce
 *                2 - stop on bounce/scattering
 */
Double_t Helium3::Advance(Double_t time, AdvanceOption bounceOption) {
  Double_t total = 0.;
  Double_t next_time = GetLifetime() + time;
  while (GetLifetime() < next_time) {
    time = next_time - GetLifetime();
    Double_t step_ret = Neutron::Advance(time, ADVANCE_STOP_ON_BOUNCE);
    if (step_ret < 0) {
      return -1;
    }
    int collision = GetLastBounce();
    if (collision < 0) {
      return -1;
    }
    total += step_ret;

    if (collision > 0) {
      Double_t new_speed = -1;
      if (vDistribution) {
	// Thermalization
	new_speed = vDistribution->GetRandom();
	if (new_speed < 0.) {
	  cout << "W: vDistribution returns a value " << new_speed << " < 0" << ". The previous value will be used instead" << endl;
	  new_speed = -1;
	}
      }
      Double_t *vel_vec = GetVelocity();
      if (collision == 1) {
	if (new_speed > 0) {
	  Vector::Normalize(vel_vec, new_speed);
	}
      }
      else if (collision == 2) {
	nScat++;
	if (new_speed < 0) {
	  new_speed = Vector::Norm(vel_vec);
	}
	gRandom->Sphere(vel_vec[0], vel_vec[1], vel_vec[2], new_speed);
      }
      if (ADVANCE_STOP_ON_BOUNCE) {
	return total;
      }
    }
  }
  return total;
}
