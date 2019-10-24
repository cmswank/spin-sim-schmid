#include "BFieldROOT.h"

/**
 * Implementation of getFactor
 */
double BTF1_Factor::getFactor(BFieldVars &vars) {
  if (!function)
    return 0;
  return function->Eval(vars.t);
}

/**
 * Implementation of getFactor
 */
double BTFormula_Factor::getFactor(BFieldVars &vars) {
  if (!formula4)
    return 0;
  return formula4->Eval(vars.pos[0], vars.pos[1], vars.pos[2], vars.t);
}

/**
 * Implementation of getFactor for Spin Dressing why are my parameters wrong?
 */
double BTSDFormula_Factor::getFactor(BFieldVars &vars) {
  if (!formula4)
    return 0;
  return (formula4->Eval(vars.pos[0], vars.pos[1], vars.pos[2], vars.t))*(this->BSDF->getFactor(vars));
}
