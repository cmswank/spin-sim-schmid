#include "Boundary.h"

#include <iostream>
#include <math.h>
using namespace std;

#include "TRandom.h"

#include "Vector.h"


/**
 * Set value for gravity
 */
void Boundary::SetGravity(Double_t gravity) {
  this->gravity = gravity;
}

/**
 * Set value for diffuse bounces
 */
void Boundary::SetDiffusion(Double_t diff) {
  this->diff = diff;
}

/**
 * Create a particle inside a boundary
 */
void Boundary::CreateInside(Double_t *position) {
  position[0] = position[1] = position[2] = 0;
}

/**
 * Check if the particle is inside (-1), outside (1)
 * or right on the boundary (0)
 */
Int_t Boundary::CheckBoundary(Double_t *position) {
  return -1;
}

/**
 * Check the position of the next move
 */
Int_t Boundary::CheckMove(Double_t *position, Double_t *velocity, Double_t time) {
  Double_t tempPos[3];
  Double_t tempVel[3];
  Vector::Copy(position, tempPos);
  Vector::Copy(velocity, tempVel);
  this->Move(tempPos, tempVel, time);
  return this->CheckBoundary(tempPos);
}

/**
 * Calculate how long until the next wall is reached
 */
Double_t Boundary::ReachBoundary(Double_t *position, Double_t *velocity) {
  return 1e3;
}

/**
 * Reflect from a wall
 */
Int_t Boundary::Reflect (Double_t *position, Double_t *velocity) {
  for (int i=0; i < 3; i++) {
    velocity[i] *= -1;
  }
  return 1;
}

/**
 * Old functionality
 */
Int_t Boundary::RectDiff(Int_t wall, Double_t *velocity) {

  Double_t speed, tempVel[3];

  speed = Vector::Norm(velocity);
  
  if (gRandom->Rndm() < this->diff) {
    
    gRandom->Sphere(tempVel[0], tempVel[1], tempVel[2], 1.0);
    
    switch (wall) {
    case 1:
      if (tempVel[0] < 0.)
	tempVel[0] = -tempVel[0];
      break;
    case 2:
      if (tempVel[0] > 0.)
	tempVel[0] = -tempVel[0];
      break;
    case 3:
      if (tempVel[1] < 0.)
	tempVel[1] = -tempVel[1];
      break;
    case 4:
      if (tempVel[1] > 0.)
	tempVel[1] = -tempVel[1];
      break;
    case 5:
      if (tempVel[2] < 0.)
	tempVel[2] = -tempVel[2];
      break;
    case 6:
      if (tempVel[2] > 0.)
	tempVel[2] = -tempVel[2];
      break;

    }
    
    Vector::Normalize(tempVel, speed);

    velocity[0] = tempVel[0];
    velocity[1] = tempVel[1];
    velocity[2] = tempVel[2];
    

    return 1;

  }


  return 0;  


}

/**
 * Move a particle under the effect of velocity and possibly gravity
 * calculate the next position and velocity.
 * only works for box type geometry
 */
void Boundary::Move(Double_t *position, Double_t *velocity, Double_t time) {
  position[0] += velocity[0]*time;
  position[1] += velocity[1]*time;
  position[2] += velocity[2]*time;
  // add gravity
  if (this->gravity != 0) {
    position[1] += .5*this->gravity*time*time;
    velocity[1] += this->gravity*time;
  }
}

/**
 * Specular bounce off wall
 */
int Boundary::SpecularBounce(Double_t *velocity, Double_t *normal) {
  Double_t dot = Vector::DotProduct(velocity, normal);
  if (dot > 0) {
    return 0;
  }
  for (int i=0; i<3; i++) {
    velocity[i] -= 2*dot*normal[i];
  }
  return 1;
}

/**
 * Diffuse bounce off wall
 */
int Boundary::DiffuseBounce(Double_t *velocity, Double_t *normal) {
  Double_t speed = Vector::Norm(velocity);
  Double_t v_normal = Vector::DotProduct(velocity, normal);
  Double_t costh_i = v_normal / speed;
  if (costh_i > 0) {
    cout << "W: NO BOUNCE, adjusting" << endl;
    costh_i *= -1;
    v_normal *= -1;
    Vector::Normalize(normal, -1);
  }
  if (gRandom->Rndm() < -costh_i * diff) {
    // Steyerl approach
    Double_t th, phi, cosphi, sinphi, sinth, costh;
    // get sin\th cos^2\th distro
    // using inverse y = cos^3 \th
    th = acos(pow(gRandom->Rndm(), .333333333333));
    phi = gRandom->Rndm() * 6.28318530717959; // 2 pi

    Double_t v0[3];
    for (int i=0; i < 3; i++) {
      v0[i] = velocity[i] - v_normal*normal[i];
    }
    Vector::Normalize(v0);

    cosphi = cos(phi);
    sinphi = sin(phi);

    sinth = sin(th);
    costh = cos(th);

    Double_t v1[3];
    Vector::CrossProduct(normal, v0, v1);
    for (int i=0; i < 3; i++) {
      v1[i] = v0[i]*cosphi + v1[i]*sinphi;
      v1[i] *= sinth;
      velocity[i] = v1[i] + costh*normal[i];
    }
    Vector::Normalize(velocity, speed);
    return 1;
  }
  return 0;
}

/**
 * Implementation of CreateInside
 */
void CircleBoundary::CreateInside(Double_t *position) {
  Double_t r = sqrt(gRandom->Rndm()) * (radius - MARGIN);
  gRandom->Circle(*position, *(position+1), r);
  position[2] = 0.;
}

/**
 * Implementation of CheckBoundary
 */
Int_t CircleBoundary::CheckBoundary(Double_t *position) {

  Double_t d = sqrt(pow(position[0], 2) + pow(position[1],2)) - radius;
  if (d > -MARGIN) {
    if (d > 0) {
      return 1; // outside
    }
    return 0; // on margin
  }
  return -1; // inside
}

/**
 * Implementation of ReachBoundary
 */
Double_t CircleBoundary::ReachBoundary(Double_t *position, Double_t *velocity) {

  // solve for line circle intersection
  Double_t a = pow(velocity[0],2) + pow(velocity[1],2);
  Double_t b2 = position[0] * velocity[0] + position[1] * velocity[1];
  Double_t c = pow(position[0],2) + pow(position[1],2) - pow(radius-MARGINH,2);

  Double_t det = pow(b2,2)-a*c;
  if (det < 0) {
    // no solution, return "infinity"
    return 1e30;
  }

  Double_t time = (Double_t)((-b2 + sqrt(det)) / a);
  return time;
}

/**
 * Implementation of Reflect
 */
Int_t CircleBoundary::Reflect (Double_t *position, Double_t *velocity) {

  // check to see if on boundary
  if (this->CheckBoundary(position) != 0) {
    cout << endl << "ERROR: not on boundary" << endl;
    cout << position[0] << ", " << position[1] << ", " << position[2] << endl;
    return -1;
  }

  Double_t posVector[3];
  Double_t temp;

  posVector[0] = position[0];
  posVector[1] = position[1];
  posVector[2] = 0.;

  velocity[2] = 0.;
  Vector::Normalize(posVector);
  // add diffusive bouncing for 2D
  // TODO: this should probably depend on the velocity normal to the wall
  if (this->diff > 0 && gRandom->Rndm() < this->diff) {
    temp = Vector::Norm(velocity);
    gRandom->Circle(velocity[0], velocity[1], temp);
  }
  // specular reflection if velocity opposite to wall normal
  Double_t dotProduct = Vector::DotProduct(posVector, velocity);
  if (dotProduct > 0.) {
    velocity[0] -= 2 * dotProduct * posVector[0];
    velocity[1] -= 2 * dotProduct * posVector[1];
  }
  return 1;
}

/**
 * Box constructor
 */
BoxBoundary::BoxBoundary(Double_t *px, Double_t *py, Double_t *pz) : Boundary::Boundary() {

  this->px[0] = px[0];
  this->px[1] = px[1];
  this->py[0] = py[0];
  this->py[1] = py[1];
  this->pz[0] = pz[0];
  this->pz[1] = pz[1];

}

/**
 * Box constructor
 */
BoxBoundary::BoxBoundary(Double_t *box_low, Double_t *box_high) : Boundary::Boundary() {

  this->px[0] = box_low[0];
  this->px[1] = box_high[0];
  this->py[0] = box_low[1];
  this->py[1] = box_high[1];
  this->pz[0] = box_low[2];
  this->pz[1] = box_high[2];

}


/**
 * Implementation of CreateInside
 */
void BoxBoundary::CreateInside(Double_t *position) {
  position[0] = px[0] + MARGIN + gRandom->Rndm() * (px[1] - px[0] - 2*MARGIN);
  position[1] = py[0] + MARGIN + gRandom->Rndm() * (py[1] - py[0] - 2*MARGIN);
  position[2] = pz[0] + MARGIN + gRandom->Rndm() * (pz[1] - pz[0] - 2*MARGIN);
}

/**
 * Implementation of CheckBoundary
 */
Int_t BoxBoundary::CheckBoundary(Double_t *position) {
  // Modified Jan 26, 2011 to address a bug in which particles were trapped outside the cell
  // now no margin outside of the cell. Particles must stay inside box at all times

  // Check if outside margin box
  if (position[0] < px[0] + MARGIN ||
      position[0] > px[1] - MARGIN ||
      position[1] < py[0] + MARGIN ||
      position[1] > py[1] - MARGIN ||
      position[2] < pz[0] + MARGIN ||
      position[2] > pz[1] - MARGIN) {
    // Check if outside box
    if (position[0] < px[0] ||
	position[0] > px[1] ||
	position[1] < py[0] ||
	position[1] > py[1] ||
	position[2] < pz[0] ||
	position[2] > pz[1]) {
      return 1; // Outside boundary
    }
    return 0; // On boundary
  }
  return -1; // In boundary
}

/**
 * Helper method to reach a single wall
 */
Double_t BoxBoundary::reach_walls(Double_t *walls, Double_t pos, Double_t vel) {
  if (vel < 0) {
    Double_t temp = (walls[0] + MARGINH - pos) / vel;
    if (temp > 0.)
      return temp;
  }
  else {
    Double_t temp = (walls[1] - MARGINH - pos) / vel;
    if (temp > 0.)
      return temp;
  }
  return 1e30;
}

/**
 * Helper method to reach a single wall considering
 * the effect of gravity
 */
Double_t BoxBoundary::reach_walls_gravity(Double_t *walls, Double_t pos, Double_t vel) {
  // solve for time in intersection of parabola with y-plane
  Double_t a;
  Double_t v2 = vel*vel;
  Double_t ret = 1e30;
  Double_t temp;
  // bottom
  a = pos - (walls[0] + MARGINH);
  a = v2 - 2*gravity*a;
  if (a >= 0.) {
    a = sqrt(a);
    // solution with negative velocity (v = g*t)
    temp = (-vel - a) / gravity;
    if (temp > 0. && temp < ret) {
      ret = temp;
    }
  }
  // top
  a = pos - (walls[1] - MARGINH);
  a = v2 - 2*gravity*a;
  if (a >= 0.) {
    a = sqrt(a);
    // solution with positive velocity (v = g*t)
    temp = (-vel + a) / gravity;
    if (temp > 0. && temp < ret) {
      ret = temp;
    }
  }
  return ret;
}

/**
 * Returns time to reach closest boundary wall
 */
Double_t BoxBoundary::ReachBoundary(Double_t *position, Double_t *velocity) {

  Double_t time = 1e50;
  Double_t temp;

  // x component
  temp = reach_walls(px, position[0], velocity[0]);
  if (temp < time) {
    time = temp;
  }
  // y component (no gravity)
  if (gravity == 0) {
    temp = reach_walls(py, position[1], velocity[1]);
  }
  else {
    temp = reach_walls_gravity(py, position[1], velocity[1]);
  }
  if (temp < time) {
    time = temp;
  }
  // z component
  temp = reach_walls(pz, position[2], velocity[2]);
  if (temp < time) {
    time = temp;
  }
  return time;
}

/**
 * Implementation of Reflect
 */
Int_t BoxBoundary::Reflect (Double_t *position, Double_t *velocity) {

  if (this->CheckBoundary(position) != 0) {
    cout << endl << "ERROR: not on boundary" << endl;
    cout << position[0] << ", " << position[1] << ", " << position[2] << endl;
    return -1;
  }
  // get the diffusion
  int bounce_z;
  int bounce_count;
  Double_t depth, max_depth;;

  bounce_count = 0;
  max_depth = -1.;
  bounce_z = -1;
  Double_t norm_dir = 1;

  // Check for closest wall
  // X wall bounce
  if ((depth = (px[0] + MARGIN) - position[0]) > 0) {
    bounce_count++;
    if (velocity[0] < 0. && depth > max_depth) {
      bounce_z = 0;
      norm_dir = 1;
      max_depth = depth;
    }
  }
  else if ((depth = position[0] - (px[1] - MARGIN)) > 0) {
    bounce_count++;
    if (velocity[0] > 0. && depth > max_depth) {
      bounce_z = 0;
      norm_dir = -1;
      max_depth = depth;
    }
  }
  // Y wall bounce
  if ((depth = (py[0] + MARGIN) - position[1]) > 0) {
    bounce_count++;
    if (velocity[1] < 0. && depth > max_depth) {
      bounce_z = 1;
      norm_dir = 1;
      max_depth = depth;
    }
  }
  else if ((depth = position[1] - (py[1] - MARGIN)) > 0) {
    bounce_count++;
    if (velocity[1] > 0. && depth > max_depth) {
      bounce_z = 1;
      norm_dir = -1;
      max_depth = depth;
    }
  }
  // Z wall bounce
  if ((depth = (pz[0] + MARGIN) - position[2]) > 0) {
    bounce_count++;
    if (velocity[2] < 0. && depth > max_depth) {
      bounce_z = 2;
      norm_dir = 1;
      max_depth = depth;
    }
  }
  else if ((depth = position[2] - (pz[1] - MARGIN)) > 0) {
    bounce_count++;
    if (velocity[2] > 0. && depth > max_depth) {
      bounce_z = 2;
      norm_dir = -1;
      max_depth = depth;
    }
  }

  if (bounce_z < 0) {
    // no wall to bounce from
    return 0;
  }

  int reflected = 0;
  if (this->diff > 0) {
    Double_t normal[3] = {0, 0, 0};
    normal[bounce_z] = norm_dir;
    reflected = DiffuseBounce(velocity, normal);
  }

  if (!reflected) {
    // normal bounce
    velocity[bounce_z] = -velocity[bounce_z];
  }

  if (bounce_count > 1) {
    // more than one wall encountered,
    // make sure to not bounce out of the box
    return 1 + this->Reflect(position, velocity);
  }
  return 1;
}

/**
 * Constructor
 */
CylindricalBoundary::CylindricalBoundary(Double_t radius, Double_t z_low, Double_t z_high) : BoxBoundary(), CircleBoundary(radius) {
  pz[0] = z_low;
  pz[1] = z_high;
}

/**
 * Override gravity direction
 */
void CylindricalBoundary::Move(Double_t *position, Double_t *velocity, Double_t time) {
  position[0] += velocity[0]*time;
  position[1] += velocity[1]*time;
  position[2] += velocity[2]*time;
  // add gravity
  if (this->gravity != 0) {
    position[2] += .5*this->gravity*time*time;
    velocity[2] += this->gravity*time;
  }
}

/**
 * Implementation of CreateInside
 */
void CylindricalBoundary::CreateInside(Double_t *position) {
  CircleBoundary::CreateInside(position);
  position[2] = pz[0] + BoxBoundary::MARGIN + gRandom->Rndm() * (pz[1] - pz[0] - 2*BoxBoundary::MARGIN);
}

/**
 * Implementation of CheckBoundary
 */
Int_t CylindricalBoundary::CheckBoundary(Double_t *position) {
  Int_t ret = -1;
  if (position[2] < pz[0] + BoxBoundary::MARGIN || position[2] > pz[1] - BoxBoundary::MARGIN) {
    if (position[2] < pz[0] || position[2] > pz[1])
      return 1;
    ret = 0;
  }
  // now need to check circular boundary no matter what
  Int_t check_circle = CircleBoundary::CheckBoundary(position);
  if (check_circle > ret)
    return check_circle;
  return ret;
}

/**
 * Implementation of ReachBoundary
 */
Double_t CylindricalBoundary::ReachBoundary(Double_t *position, Double_t *velocity) {
  // circle part:
  Double_t ret = CircleBoundary::ReachBoundary(position, velocity);
  // z part
  Double_t temp;
  if (BoxBoundary::gravity != 0) {
    temp = BoxBoundary::reach_walls_gravity(pz, position[2], velocity[2]);
  }
  else {
    temp = BoxBoundary::reach_walls(pz, position[2], velocity[2]);
  }
  if (temp < ret) {
    ret = temp;
  }
  return ret;
}

/**
 * Implementation of Reflect
 */
Int_t CylindricalBoundary::Reflect(Double_t *position, Double_t *velocity) {
  if (CheckBoundary(position) != 0) {
    cout << endl << "ERROR: not on boundary" << endl;
    cout << position[0] << ", " << position[1] << ", " << position[2] << endl;
    return -1;
  }
  Int_t bounce_count = 0;
  Int_t reflection = 0;
  Double_t normal[3] = {0,0,0};
  // Z wall bounce
  if ((pz[0] + BoxBoundary::MARGIN) - position[2] > 0) {
    bounce_count++;
    if (velocity[2] < 0) {
      normal[2] = 1;
      reflection = 1;
    }
  }
  else if (position[2] - (pz[1] - BoxBoundary::MARGIN) > 0) {
    bounce_count++;
    if (velocity[2] > 0) {
      normal[2] = -1;
      reflection = 1;
    }
  }
  Double_t r = sqrt(position[0]*position[0] + position[1]*position[1]);
  if (r > (radius - BoxBoundary::MARGIN)) {
    bounce_count++;
    if (!reflection) {
      normal[0] = -position[0];
      normal[1] = -position[1];
      normal[2] = 0;
      if (Vector::DotProduct(velocity, normal) < 0) {
	Vector::Normalize(normal);
	reflection = 1;
      }
    }
  }
  Int_t reflected = 0;
  if (reflection) {
    if (BoxBoundary::diff > 0) {
      reflected = BoxBoundary::DiffuseBounce(velocity, normal);
    }
    if (!reflected) {
      reflected = BoxBoundary::SpecularBounce(velocity, normal);
    }
  }
  if (reflected && bounce_count > 1) {
    return 1 + Reflect(position, velocity);
  }
  return reflected;
}

/**
 * Implementation of CreateInside
 */
void SphereBoundary::CreateInside(Double_t *position) {
  Double_t r = pow(gRandom->Rndm(),1./3) * radius;
  gRandom->Sphere(*position, *(position+1), *(position+2), r);
}

/**
 * Implementation of CheckBoundary
 */
Int_t SphereBoundary::CheckBoundary(Double_t *position) {
  Double_t d = Vector::Norm(position) - radius;
  if (d > -MARGIN) {
    if (d > 0)
      return 1;
    return 0;
  }
  return -1;
}

/**
 * Implementation of ReachBoundary
 */
Double_t SphereBoundary::ReachBoundary(Double_t *position, Double_t *velocity) {

  Double_t time = 0;

  Double_t a, b2, c;

  // solve for line sphere intersection
  a = Vector::Norm2(velocity);
  b2 = Vector::DotProduct(position, velocity);
  c = Vector::Norm2(position) - pow(radius-MARGINH, 2);

  Double_t dt = (Double_t)((-b2 + sqrt(b2*b2 - a*c)) / a);

  // Make sure to end up on boundary
  Double_t tempPos[3];
  Double_t tempVel[3];
  int check;
  do {
    time += dt;
    Vector::Copy(position, tempPos);
    Vector::Copy(velocity, tempVel);
    Move(tempPos, tempVel, time);
    check = CheckBoundary(tempPos);
    if (check < 0) {
      // recursevly find boundary
      dt = ReachBoundary(tempPos, tempVel);
    }
    else if (check > 0) {
      // simple binary search
      time -= dt;
      dt /= 2;
    }
  } while (check != 0);
  return time;
}

/**
 * Implementation of Reflect
 */
Int_t SphereBoundary::Reflect (Double_t *position, Double_t *velocity) {

  if (this->CheckBoundary(position) != 0) {
    cout << endl << "ERROR: not on boundary" << endl;
    cout << position[0] << ", " << position[1] << ", " << position[2] << endl;
    return -1;
  }

  Double_t posVector[3];

  Double_t temp;

  posVector[0] = position[0];
  posVector[1] = position[1];
  posVector[2] = position[2];

  Vector::Normalize(posVector);
  Double_t dotProduct = Vector::DotProduct(posVector, velocity);

  if (dotProduct < 0) {
    return 0;
  }

  velocity[0] -= 2 * dotProduct * posVector[0];
  velocity[1] -= 2 * dotProduct * posVector[1];
  velocity[2] -= 2 * dotProduct * posVector[2];

  if (diff > 0) {
    // Diffuse if necessary
    // TODO: Fix to use Lambertian
    if (gRandom->Rndm() < this->diff) {
      temp = Vector::Norm(velocity);
      gRandom->Sphere(velocity[0],velocity[1],velocity[2],temp);
    }
  }

  return 1;
}

/**
 * Box 2D constructor
 */
Box2DBoundary::Box2DBoundary(Double_t *box1, Double_t *box2) : BoxBoundary() {
  px[0] = box1[0];
  py[0] = box1[1];

  px[1] = box2[0];
  py[1] = box2[1];
}

/**
 * Implementation of CreateInside
 */
void Box2DBoundary::CreateInside(Double_t *position) {
  position[0] = px[0] + MARGIN + gRandom->Rndm() * (px[1] - px[0] - 2*MARGIN);
  position[1] = py[0] + MARGIN + gRandom->Rndm() * (py[1] - py[0] - 2*MARGIN);
  position[2] = 0.;
}

/**
 * Implementation of CheckBoundary
 */
Int_t Box2DBoundary::CheckBoundary(Double_t *position) {
  if (position[0] < px[0] + MARGIN ||
      position[0] > px[1] - MARGIN ||
      position[1] < py[0] + MARGIN ||
      position[1] > py[1] - MARGIN) {
    if (position[0] < px[0] ||
	position[0] > px[1] ||
	position[1] < py[0] ||
	position[1] > py[1])
      return 1;
    return 0;
  }
  return -1;
}

void Box2DBoundary::Move(Double_t *position, Double_t *velocity, Double_t time) {
  position[0] += velocity[0] * time;
  position[1] += velocity[1] * time;
}


/**
 * Implementation of ReachBoundary
 */
Double_t Box2DBoundary::ReachBoundary(Double_t *position, Double_t *velocity) {

  Double_t times[2];

  times[0] = reach_walls(px, position[0], velocity[0]);
  times[1] = reach_walls(py, position[1], velocity[1]);
  if (times[0] < times[1])
    return times[0];
  return times[1];
}

/**
 * Implementation of Reflect
 */
Int_t Box2DBoundary::Reflect (Double_t *position, Double_t *velocity) {
  if (this->CheckBoundary(position) != 0) {
    cout << endl << "ERROR: not on boundary: "
	 << position[0] << ", " << position[1] << endl;
    return -1;
  }
  int bounce_z;
  int bounce_count;
  Double_t depth, max_depth;;

  bounce_count = 0;
  max_depth = -1.;
  bounce_z = -1;
  Double_t norm_dir = 1;

  // Check for closest wall
  // X wall bounce
  if ((depth = (px[0] + MARGIN) - position[0]) > 0) {
    bounce_count++;
    if (velocity[0] < 0. && depth > max_depth) {
      bounce_z = 0;
      norm_dir = 1;
      max_depth = depth;
    }
  }
  else if ((depth = position[0] - (px[1] - MARGIN)) > 0) {
    bounce_count++;
    if (velocity[0] > 0. && depth > max_depth) {
      bounce_z = 0;
      norm_dir = -1;
      max_depth = depth;
    }
  }
  // Y wall bounce
  if ((depth = (py[0] + MARGIN) - position[1]) > 0) {
    bounce_count++;
    if (velocity[1] < 0. && depth > max_depth) {
      bounce_z = 1;
      norm_dir = 1;
      max_depth = depth;
    }
  }
  else if ((depth = position[1] - (py[1] - MARGIN)) > 0) {
    bounce_count++;
    if (velocity[1] > 0. && depth > max_depth) {
      bounce_z = 1;
      norm_dir = -1;
      max_depth = depth;
    }
  }

  if (bounce_z < 0) {
    // no wall to bounce from
    return 0;
  }

  int reflected = 0;
  if (this->diff > 0) {
    Double_t normal[3] = {0, 0, 0};
    normal[bounce_z] = norm_dir;
    Double_t speed = Vector::Norm(velocity);
    reflected = DiffuseBounce(velocity, normal);
    if (reflected) {
      velocity[2] = 0.;
      Vector::Normalize(velocity, speed);
    }
  }

  if (!reflected) {
    // normal bounce
    velocity[bounce_z] = -velocity[bounce_z];
  }

  if (bounce_count > 1) {
    // more than one wall encountered,
    // make sure to not bounce out of the box
    return 1 + this->Reflect(position, velocity);
  }
  return 1;
}

/**
 * Constructor
 */
Boundary1D::Boundary1D(Double_t wall_low, Double_t wall_high) : BoxBoundary()  {
  px[0] = wall_low;
  px[1] = wall_high;
}

/**
 * Implementation of CreateInside
 */
void Boundary1D::CreateInside(Double_t *position) {
  position[0] = px[0] + MARGIN + gRandom->Rndm() * (px[1] - px[0] - 2*MARGIN);
  position[1] = 0.;
  position[2] = 0.;
}

/**
 * Implementation of Move
 */
void Boundary1D::Move(Double_t *position, Double_t *velocity, Double_t time) {
  position[0] += velocity[0] * time;
}

/**
 * Implementation of CheckBoundary
 */
Int_t Boundary1D::CheckBoundary(Double_t *position) {
  if (position[0] < px[0] + MARGIN || position[0] > px[1] - MARGIN) {
    if (position[0] < px[0] || position[0] > px[1])
      return 1;
    return 0;
  }
  return -1;
}

/**
 * Implementation of ReachBoundary
 */
Double_t Boundary1D::ReachBoundary(Double_t *position, Double_t *velocity) {
  return reach_walls(px, position[0], velocity[0]);
}

/**
 * Implementation of Reflect
 */
Int_t Boundary1D::Reflect(Double_t *position, Double_t *velocity) {
  if (CheckBoundary(position) != 0) {
    cout << endl << "ERROR: not on boundary: (x = " << position[0] << endl;
    return -1;
  }
  velocity[0] = -velocity[0];
  return 1;
}
