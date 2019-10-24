#ifndef __BFIELD_H__
#define __BFIELD_H__

/**
 * Class to hold the arguments necessary to specify the field
 */
class BFieldVars {
 public:
  BFieldVars() : t(0) {
    pos[0] = pos[1] = pos[2] = 0;
    vel[0] = vel[1] = vel[2] = 0;
  }
  BFieldVars(double t, double *pos, double *vel);

  double t;
  double pos[3];
  double vel[3];
};

/**
 * Base class for fields
 */
class BField
{
 public:
  BField() {};
  virtual ~BField() {};
  virtual int getField(double *out, BFieldVars &vars) = 0;
};

/**
 * Constant field initialized by a vector
 */
class BFieldConst : public BField
{
 public:
  BFieldConst(double *field = 0);
  virtual int getField(double *out, BFieldVars &vars);
 protected:
  double field_const[3];
};

/**
 * Base class for multiplier factors
 */
class BFactor
{
 public:
  BFactor(double f = 1.) : factor(f) {};
  virtual ~BFactor() {};
  virtual double getFactor(BFieldVars &vars);
 protected:
  double factor;
};

#endif
