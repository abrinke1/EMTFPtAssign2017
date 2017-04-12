
#include "../interface/PtLutVarCalc.h"
#include "../src/PtAssignmentEngineAux2017.cc"

PtAssignmentEngineAux2017 ENG;


int CalcTrackTheta( const int th1, const int th2, const int th3, const int th4,
		    const int st1_ring2, const int mode, const bool BIT_COMP ) {

  int theta = -99;

  if      ( (mode % 8) / 4 > 0 ) // Has station 2 hit
    theta = th2;
  else if ( (mode % 4) / 2 > 0 ) // Has station 3 hit
    theta = th3;
  else if ( (mode % 2) > 0 ) // Has station 4 hit
    theta = th4;

  assert( theta > 0 );

  if (BIT_COMP) {
    int bits = (mode == 15 ? 4 : 5);
    theta = ENG.getTheta(theta, st1_ring2, bits);
  }

  return theta;
}


void CalcDeltaPhis( int& dPh12, int& dPh13, int& dPh14, int& dPh23, int& dPh24, int& dPh34, int& dPhSign,
		    int& dPhSum4, int& dPhSum4A, int& dPhSum3, int& dPhSum3A, int& outStPh,
		    const int ph1, const int ph2, const int ph3, const int ph4, const int mode, const bool BIT_COMP ) {

  dPh12 = ph2 - ph1;
  dPh13 = ph3 - ph1;
  dPh14 = ph4 - ph1;
  dPh23 = ph3 - ph2;
  dPh24 = ph4 - ph2;
  dPh34 = ph4 - ph3;
  dPhSign = 0;

  // std::cout << "Computing dPhi signs" << std::endl;
  if (mode >= 8) {                   // First hit is station 1
    if      ( (mode % 8) / 4 > 0 )   // Has station 2 hit
      dPhSign = dPh12 / max(1, abs(dPh12));
    else if ( (mode % 4) / 2 > 0 )   // Has station 3 hit
      dPhSign = dPh13 / max(1, abs(dPh13));
    else if ( (mode % 2) > 0 )       // Has station 4 hit
      dPhSign = dPh14 / max(1, abs(dPh14));
  } else if ( (mode % 8) / 4 > 0 ) { // First hit is station 2
    if      ( (mode % 4) / 2 > 0 )   // Has station 3 hit
      dPhSign = dPh23 / max(1, abs(dPh23));
    else if ( (mode % 2) > 0 )       // Has station 4 hit
      dPhSign = dPh24 / max(1, abs(dPh24));
  } else if ( (mode % 4) / 2 > 0 ) { // First hit is station 3
    if      ( (mode % 2) > 0 )       // Has station 4 hit
      dPhSign = dPh34 / max(1, abs(dPh34));
  }

  if (dPhSign == 0)
    dPhSign = 1;

  dPh12 *= dPhSign;
  dPh13 *= dPhSign;
  dPh14 *= dPhSign;
  dPh23 *= dPhSign;
  dPh24 *= dPhSign;
  dPh34 *= dPhSign;

  if (BIT_COMP) {
    // std::cout << "Doing compressed bit calculations" << std::endl;
    dPh12 = ENG.getNLBdPhi(dPh12, 7, 512);
    dPh13 = ENG.getNLBdPhi(dPh13, 7, 512);
    dPh14 = ENG.getNLBdPhi(dPh14, 7, 512);
    dPh23 = ENG.getNLBdPhi(dPh23, 7, 512);
    dPh24 = ENG.getNLBdPhi(dPh24, 7, 512);
    dPh34 = ENG.getNLBdPhi(dPh34, 7, 512);

    if (mode == 15)
      dPh34 = ENG.getNLBdPhi(dPh34, 4, 256);

    if (mode == 7 || mode == 11 || mode > 12) {
      if (mode != 7)
	dPh23 = ENG.getNLBdPhi(dPh23, 5, 256);
      dPh24 = ENG.getNLBdPhi(dPh23, 5, 256);
      dPh34 = ENG.getNLBdPhi(dPh23, 5, 256);
    }
  } // End conditional: if (BIT_COMP)


  // Compute summed quantities
  if (mode == 15) {
    dPhSum4  = dPh12 + dPh13 + dPh14 + dPh23 + dPh24 + dPh34;
    dPhSum4A = abs(dPh12) + abs(dPh13) + abs(dPh14) + abs(dPh23) + abs(dPh24) + abs(dPh34);
    int devSt1 = abs(dPh12) + abs(dPh13) + abs(dPh14);
    int devSt2 = abs(dPh12) + abs(dPh23) + abs(dPh24);
    int devSt3 = abs(dPh13) + abs(dPh23) + abs(dPh34);
    int devSt4 = abs(dPh14) + abs(dPh24) + abs(dPh34);
    
    if      (devSt4 > devSt3 && devSt4 > devSt2 && devSt4 > devSt1)  outStPh = 4;
    else if (devSt3 > devSt4 && devSt3 > devSt2 && devSt3 > devSt1)  outStPh = 3;
    else if (devSt2 > devSt4 && devSt2 > devSt3 && devSt2 > devSt1)  outStPh = 2;
    else if (devSt1 > devSt4 && devSt1 > devSt3 && devSt1 > devSt2)  outStPh = 1;
    else                                                             outStPh = 0;
    
    if      (outStPh == 4) {
      dPhSum3  = dPh12 + dPh13 + dPh23;
      dPhSum3A = abs(dPh12) + abs(dPh13) + abs(dPh23);
    } else if (outStPh == 3) {
      dPhSum3  = dPh12 + dPh14 + dPh24;
      dPhSum3A = abs(dPh12) + abs(dPh14) + abs(dPh24);
    } else if (outStPh == 2) {
      dPhSum3  = dPh13 + dPh14 + dPh34;
      dPhSum3A = abs(dPh13) + abs(dPh14) + abs(dPh34);
    } else {
      dPhSum3  = dPh23 + dPh24 + dPh34;
      dPhSum3A = abs(dPh23) + abs(dPh24) + abs(dPh34);
    }
  }

} // End function: CalcDeltaPhis()



void CalcDeltaThetas( int& dTh12, int& dTh13, int& dTh14, int& dTh23, int& dTh24, int& dTh34,
		      const int th1, const int th2, const int th3, const int th4, const int mode, const bool BIT_COMP ) {
  
  dTh12 = th2 - th1;
  dTh13 = th3 - th1;
  dTh14 = th4 - th1;
  dTh23 = th3 - th2;
  dTh24 = th4 - th2;
  dTh34 = th4 - th3;

  if (BIT_COMP) {
    dTh12 = ENG.getdTheta(dTh12, 3);
    dTh13 = ENG.getdTheta(dTh13, 3);
    dTh14 = ENG.getdTheta(dTh14, 3);
    dTh23 = ENG.getdTheta(dTh23, 3);
    dTh24 = ENG.getdTheta(dTh24, 3);
    dTh34 = ENG.getdTheta(dTh34, 3);

    if (mode == 15)
      dTh14 = ENG.getdTheta(dTh14, 2);
  } // End conditional: if (BIT_COMP)

} // CalcDeltaThetas()



void CalcBends( int& bend1, int& bend2, int& bend3, int& bend4,
		const int pat1, const int pat2, const int pat3, const int pat4,
		const int dPhSign, const int endcap, const int mode, const bool BIT_COMP ) {

  bend1 = CalcBendFromPattern( pat1, endcap );
  bend2 = CalcBendFromPattern( pat2, endcap );
  bend3 = CalcBendFromPattern( pat3, endcap );
  bend4 = CalcBendFromPattern( pat4, endcap );
  
  if (BIT_COMP) {
    int bits = 3;
    if (mode == 7 || mode == 11 || mode > 12)
      bits = 2;

    bend1 = ENG.getCLCT( pat1, endcap, dPhSign, bits );
    bend2 = ENG.getCLCT( pat2, endcap, dPhSign, bits );
    bend3 = ENG.getCLCT( pat3, endcap, dPhSign, bits );
    bend4 = ENG.getCLCT( pat4, endcap, dPhSign, bits );
  } // End conditional: if (BIT_COMP)

} // End function: CalcBends()



int CalcBendFromPattern( const int pattern, const int endcap ) {

  int bend = -99;
  if (pattern < 0)
    return bend;

  if (pattern == 10)
    bend = 0;
  else if ( (pattern % 2) == 0 )
    bend = (10 - pattern) / 2;
  else if ( (pattern % 2) == 1 )
    bend = -1 * (11 - pattern) / 2;

  // Reverse to match dPhi convention
  if (endcap == 1)
    bend *= -1;

  assert( bend != -99 );
  return bend;
}


