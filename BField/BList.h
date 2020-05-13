#ifndef __BLIST_H__
#define __BLIST_H__

#include <vector>

#include "BField.h"

// dummy class helps dumb rootcint
class BList {};

/**
 * Store fields in a vector<BField*> as
 * field superpositions.
 */
class BFieldAdder : public BField
{
 public:
  BFieldAdder(BField *A = 0, BField *B = 0);
  virtual ~BFieldAdder();
  virtual int append(BField *field);
  virtual int getField(double *out, BFieldVars &vars);
  std::vector<BField *> fields;
 protected:
  
};

/**
 * use BFactor to scale BField
 */
class BFieldScaler : public BField
{
 public:
  BFieldScaler(BField *field, BFactor *f = 0) : BField(), base_field(field), factor(f) {};
  virtual ~BFieldScaler();
  virtual int getField(double *out, BFieldVars &vars);
  void SetFactor(BFactor *f);
  
 protected:
  BField *base_field;
  BFactor *factor;
};

#endif
