
#include "../interface/MVA_input_var_tools.h"

Int_t CalcBendFromPattern( Int_t pattern, Double_t eta ) {

  Int_t bend = -88;
  if (pattern < 0)
    return bend;

  if (pattern == 10)
    bend = 0;
  else if ( (pattern % 2) == 0 )
    bend = (10 - pattern) / 2;
  else if ( (pattern % 2) == 1 )
    bend = -1 * (11 - pattern) / 2;

  // Reverse to match dPhi convention
  if (eta > 0)
    bend *= -1;

  return bend;
}
