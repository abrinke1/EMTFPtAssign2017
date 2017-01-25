
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

Int_t MaxOfSix( Int_t d1, Int_t d2, Int_t d3,
		Int_t d4, Int_t d5, Int_t d6) {
  
  Int_t f1 = abs(d1);
  Int_t f2 = abs(d2);
  Int_t f3 = abs(d3);
  Int_t f4 = abs(d4);
  Int_t f5 = abs(d5);
  Int_t f6 = abs(d6);
  
  if (f1 >= f2 && f1 >= f3 && f1 >= f4 && f1 >= f5 && f1 >= f6) return d1;
  if (f2 >= f1 && f2 >= f3 && f2 >= f4 && f2 >= f5 && f2 >= f6) return d2;
  if (f3 >= f1 && f3 >= f2 && f3 >= f4 && f3 >= f5 && f3 >= f6) return d3;
  if (f4 >= f1 && f4 >= f2 && f4 >= f3 && f4 >= f5 && f4 >= f6) return d4;
  if (f5 >= f1 && f5 >= f2 && f5 >= f3 && f5 >= f4 && f5 >= f6) return d5;
  if (f6 >= f1 && f6 >= f2 && f6 >= f3 && f6 >= f4 && f6 >= f5) return d6;
  
  return -88;
}

Int_t MaxOfThree( Int_t d1, Int_t d2, Int_t d3) {
  
  Int_t f1 = abs(d1);
  Int_t f2 = abs(d2);
  Int_t f3 = abs(d3);
  
  if (f1 >= f2 && f1 >= f3) return d1;
  if (f2 >= f1 && f2 >= f3) return d2;
  if (f3 >= f1 && f3 >= f2) return d3;
  
  return -88;
}
