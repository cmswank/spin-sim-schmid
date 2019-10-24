#ifndef __BOUNDARY_H
#define __BOUNDARY_H

#include "TROOT.h"

class Boundary {
 protected:
  Double_t MARGIN; // = 1e-6
  Double_t MARGIN2; // = 1e-12; // MARGIN^2
  Double_t MARGINH; // = 5e-7; // 1/2 MARGIN


  Double_t diff;
  Double_t gravity;

  Int_t RectDiff(Int_t wall, Double_t *velocity);

 public:
  Boundary() : MARGIN(1e-6), MARGIN2(1e-12), MARGINH(5e-7),
    diff(-1), gravity(0) {};

  void SetGravity(Double_t gravity);
  void SetDiffusion(Double_t diff);

  virtual void CreateInside(Double_t *position);

  virtual Int_t CheckBoundary(Double_t *position);
  virtual Int_t CheckMove(Double_t *position, Double_t *velocity, Double_t time);

  virtual Double_t ReachBoundary(Double_t *position, Double_t *velocity);


  virtual Int_t Reflect(Double_t *position, Double_t *velocity);

  virtual void Move(Double_t *position, Double_t *velocity, Double_t time);

  virtual int SpecularBounce(Double_t *velocity, Double_t *normal);
  virtual int DiffuseBounce(Double_t *velocity, Double_t *normal);

};

class CircleBoundary : virtual public Boundary {
 public:
  CircleBoundary(Double_t radius) : Boundary(), radius(radius) {};

  virtual void CreateInside(Double_t *position);
  virtual Int_t CheckBoundary(Double_t *position);

  virtual Double_t ReachBoundary(Double_t *position, Double_t *velocity);
  virtual Int_t Reflect(Double_t *position, Double_t *velocity);

 protected:
  Double_t radius;
};

class BoxBoundary : virtual public Boundary {
 protected:
  BoxBoundary() : Boundary() {};
 public:
  BoxBoundary(Double_t *px, Double_t *py, Double_t *pz);
  BoxBoundary(Double_t *box_low, Double_t *box_high);
  virtual void CreateInside(Double_t *position);
  virtual Int_t CheckBoundary(Double_t *position);
  virtual Double_t ReachBoundary(Double_t *position, Double_t *velocity);
  virtual Int_t Reflect(Double_t *position, Double_t *velocity);

 protected:
  Double_t px[2], py[2], pz[2];
  Double_t reach_walls(Double_t *walls, Double_t pos, Double_t vel);
  Double_t reach_walls_gravity(Double_t *walls, Double_t pos, Double_t vel);
};

class CylindricalBoundary : public BoxBoundary, public CircleBoundary {
 public:
  CylindricalBoundary(Double_t radius, Double_t *pz) : BoxBoundary(pz,pz,pz), CircleBoundary(radius) {};
  CylindricalBoundary(Double_t radius, Double_t z_low, Double_t z_high);

  virtual void Move(Double_t *position, Double_t *velocity, Double_t time);
  virtual void CreateInside(Double_t *position);
  virtual Int_t CheckBoundary(Double_t *position);
  virtual Double_t ReachBoundary(Double_t *position, Double_t *velocity);
  virtual Int_t Reflect(Double_t *position, Double_t *velocity);

};

class SphereBoundary : public Boundary {
 public:
  SphereBoundary(Double_t radius) : Boundary::Boundary(), radius(radius) {};

  virtual void CreateInside(Double_t *position);
  virtual Int_t CheckBoundary(Double_t *position);
  virtual Double_t ReachBoundary(Double_t *position, Double_t *velocity);
  virtual Int_t Reflect(Double_t *position, Double_t *velocity);

 protected:
  Double_t radius;
};

class Box2DBoundary : public BoxBoundary {
 public:
  Box2DBoundary(Double_t *box1, Double_t *box2);

  virtual void Move(Double_t *position, Double_t *velocity, Double_t time);
  virtual void CreateInside(Double_t *position);
  virtual Int_t CheckBoundary(Double_t *position);
  virtual Double_t ReachBoundary(Double_t *position, Double_t *velocity);
  virtual Int_t Reflect(Double_t *position, Double_t *velocity);
};

class Boundary1D : public BoxBoundary {
 public:
  Boundary1D(Double_t wall_low, Double_t wall_high);

  virtual void Move(Double_t *position, Double_t *velocity, Double_t time);
  virtual void CreateInside(Double_t *position);
  virtual Int_t CheckBoundary(Double_t *position);
  virtual Double_t ReachBoundary(Double_t *position, Double_t *velocity);
  virtual Int_t Reflect(Double_t *position, Double_t *velocity);

};

#endif

