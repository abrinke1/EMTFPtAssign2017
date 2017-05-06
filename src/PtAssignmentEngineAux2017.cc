// #include "L1Trigger/L1TMuonEndCap/interface/PtAssignmentEngineAux2017.hh"
#include "../interface/PtAssignmentEngineAux2017.hh"

// // ModeVariables is a 2D arrary indexed by [TrackMode(13 Total Listed Below)][VariableNumber(20 Total Constructed Above)]
// // Variable numbering
// // 0 = dPhi12
// // 1 = dPhi13
// // 2 = dPhi14
// // 3 = dPhi23
// // 4 = dPhi24
// // 5 = dPhi34
// // 6 = dEta12
// // 7 = dEta13
// // 8 = dEta14
// // 9 = dEta23
// // 10 = dEta24
// // 11 = dEta34
// // 12 = CLCT1
// // 13 = CLCT2
// // 14 = CLCT3
// // 15 = CLCT4
// // 16 = CSCID1
// // 17 = CSCID2
// // 18 = CSCID3
// // 19 = CSCID4
// // 20 = FR1
// // 21 = FR2
// // 22 = FR3
// // 23 = FR4

// // Bobby's Scheme3 (or "SchemeC"), with 30 bit compression //
// //3:TrackEta:dPhi12:dEta12:CLCT1:CLCT2:FR1
// //4:Single Station Track Not Possible
// //5:TrackEta:dPhi13:dEta13:CLCT1:CLCT3:FR1
// //6:TrackEta:dPhi23:dEta23:CLCT2:CLCT3:FR2
// //7:TrackEta:dPhi12:dPhi23:dEta13:CLCT1:FR1
// //8:Single Station Track Not Possible
// //9:TrackEta:dPhi14:dEta14:CLCT1:CLCT4:FR1
// //10:TrackEta:dPhi24:dEta24:CLCT2:CLCT4:FR2
// //11:TrackEta:dPhi12:dPhi24:dEta14:CLCT1:FR1
// //12:TrackEta:dPhi34:dEta34:CLCT3:CLCT4:FR3
// //13:TrackEta:dPhi13:dPhi34:dEta14:CLCT1:FR1
// //14:TrackEta:dPhi23:dPhi34:dEta24:CLCT2
// //15:TrackEta:dPhi12:dPhi23:dPhi34:FR1

// static const int ModeVariables_Scheme3[13][6] =
// {
//     {0,6,12,13,20,-999},              // 3
//     {-999,-999,-999,-999,-999,-999},  // 4
//     {1,7,12,14,20,-999},              // 5
//     {3,9,13,14,21,-999},              // 6
//     {0,3,7,12,20,-999},               // 7
//     {-999,-999,-999,-999,-999,-999},  // 8
//     {2,8,12,15,20,-999},              // 9
//     {4,10,13,15,21,-999},             // 10
//     {0,4,8,12,20,-999},               // 11
//     {5,11,14,15,22,-999},             // 12
//     {1,5,8,16,20,-999},               // 13
//     {3,5,10,13,-999,-999},            // 14
//     {0,3,5,20,-999,-999}              // 15
// };


// Arrays that map the integer dPhi --> dPhi-units. 1/60th of a degree per unit; 255 units --> 4.25 degrees, 511 --> 8.52 degrees

// 256 max units----

// For use in dPhi34 in mode 15.  Derived manually from dPhiNLBMap_5bit_256Max for now; should generate algorithmically. - AWB 17.03.17
static const int dPhiNLBMap_4bit_256Max[16] = {0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 20, 25, 31, 46, 68, 136};

// For use in dPhi23, dPhi24, and dPhi34 in 3- and 4-station modes (7, 11, 13, 14, 15), except for dPhi23 in mode 7 and dPhi34 in mode 15
static const int dPhiNLBMap_5bit_256Max[32] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 
					       16, 17, 19, 20, 21, 23, 25, 28, 31, 34, 39, 46, 55, 68, 91, 136};

// 512 max units----

// For use in all dPhiAB (where "A" and "B" are the first two stations in the track) in all modes
static const int dPhiNLBMap_7bit_512Max[128] =  {  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15, 
						  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31, 
						  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47, 
						  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63, 
						  64,  65,  66,  67,  68,  69,  71,  72,  73,  74,  75,  76,  77,  79,  80,  81, 
						  83,  84,  86,  87,  89,  91,  92,  94,  96,  98, 100, 102, 105, 107, 110, 112, 
						 115, 118, 121, 124, 127, 131, 135, 138, 143, 147, 152, 157, 162, 168, 174, 181, 
						 188, 196, 204, 214, 224, 235, 247, 261, 276, 294, 313, 336, 361, 391, 427, 470};

// const int (*PtAssignmentEngineAux2017::getModeVariables() const)[6] {
//   return ModeVariables_Scheme3;
// }

int PtAssignmentEngineAux2017::getNLBdPhi(int dPhi, int bits, int max) const {
  assert( (bits == 4 && max == 256) || 
	  (bits == 5 && max == 256) || 
	  (bits == 7 && max == 512) );

  int dPhi_ = max;
  int sign_ = 1;
  if (dPhi < 0)
    sign_ = -1;
  dPhi = sign_ * dPhi;

  if (max == 256) {
    if (bits == 4) {
      dPhi_ = dPhiNLBMap_4bit_256Max[(1 << bits) - 1];
      for (int edge = 0; edge < (1 << bits) - 1; edge++) {
        if (dPhiNLBMap_4bit_256Max[edge]  <= dPhi && 
	    dPhiNLBMap_4bit_256Max[edge+1] > dPhi) {
          dPhi_ = dPhiNLBMap_4bit_256Max[edge];
          break;
        }
      }
    } // End conditional: if (bits == 4)
    if (bits == 5) {
      dPhi_ = dPhiNLBMap_5bit_256Max[(1 << bits) - 1];
      for (int edge = 0; edge < (1 << bits) - 1; edge++) {
        if (dPhiNLBMap_5bit_256Max[edge]  <= dPhi && 
	    dPhiNLBMap_5bit_256Max[edge+1] > dPhi) {
          dPhi_ = dPhiNLBMap_5bit_256Max[edge];
          break;
        }
      }
    } // End conditional: if (bits == 5)
  } // End conditional: if (max == 256)

  else if (max == 512) {
    if (bits == 7) {
      dPhi_ = dPhiNLBMap_7bit_512Max[(1 << bits) - 1];
      for (int edge = 0; edge < (1 << bits) - 1; edge++) {
        if (dPhiNLBMap_7bit_512Max[edge]  <= dPhi && 
	    dPhiNLBMap_7bit_512Max[edge+1] > dPhi) {
          dPhi_ = dPhiNLBMap_7bit_512Max[edge];
          break;
        }
      }
    } // End conditional: if (bits == 7)
  } // End conditional: else if (max == 512)

  assert( abs(sign_) == 1 && dPhi_ >= 0 && dPhi_ < max);
  return (sign_ * dPhi_);
} // End function: nt PtAssignmentEngineAux2017::getNLBdPhi()


int PtAssignmentEngineAux2017::getNLBdPhiBin(int dPhi, int bits, int max) const {
  assert( (bits == 4 && max == 256) || 
	  (bits == 5 && max == 256) || 
	  (bits == 7 && max == 512) );
  
  int dPhiBin_ = (1 << bits) - 1;
  int sign_ = 1;
  if (dPhi < 0)
    sign_ = -1;
  dPhi = sign_ * dPhi;
  
  if (max == 256) {
    if (bits == 4) {
      for (int edge = 0; edge < (1 << bits) - 1; edge++) {
        if (dPhiNLBMap_4bit_256Max[edge] <= dPhi && 
	    dPhiNLBMap_4bit_256Max[edge+1] > dPhi) {
          dPhiBin_ = edge;
          break;
        }
      }
    } // End conditional: if (bits == 4)
    if (bits == 5) {
      for (int edge = 0; edge < (1 << bits) - 1; edge++) {
        if (dPhiNLBMap_5bit_256Max[edge]  <= dPhi && 
	    dPhiNLBMap_5bit_256Max[edge+1] > dPhi) {
          dPhiBin_ = edge;
          break;
        }
      }
    } // End conditional: if (bits == 5) 
  } // End conditional: if (max == 256)

  else if (max == 512) {
    if (bits == 7) {
      for (int edge = 0; edge < (1 << bits) - 1; edge++) {
        if (dPhiNLBMap_7bit_512Max[edge]  <= dPhi && 
	    dPhiNLBMap_7bit_512Max[edge+1] > dPhi) {
          dPhiBin_ = edge;
          break;
        }
      }
    } // End conditional: if (bits == 7)
  } // End conditional: else if (max == 512)
  
  assert(dPhiBin_ >= 0 && dPhiBin_ < pow(2, bits));
  return (dPhiBin_);
} // End function: int PtAssignmentEngineAux2017::getNLBdPhiBin()


int PtAssignmentEngineAux2017::getdPhiFromBin(int dPhiBin, int bits, int max) const {
  assert( (bits == 4 && max == 256) || 
	  (bits == 5 && max == 256) || 
	  (bits == 7 && max == 512) );
  
  int dPhi_ = (1 << bits) - 1;

  if (dPhiBin > (1 << bits) - 1)
    dPhiBin = (1 << bits) - 1;
  
  if (max == 256) {
    if (bits == 4)
      dPhi_ = dPhiNLBMap_4bit_256Max[dPhiBin];
    if (bits == 5)
      dPhi_ = dPhiNLBMap_5bit_256Max[dPhiBin];
  } // End conditional: if (max == 256)

  else if (max == 512) {
    if (bits == 7)
      dPhi_ = dPhiNLBMap_7bit_512Max[dPhiBin];
  } // End conditional: else if (max == 512)

  assert(dPhi_ >= 0 && dPhi_ < max);
  return (dPhi_);
} // End function: int PtAssignmentEngineAux2017::getdPhiFromBin()


int PtAssignmentEngineAux2017::getCLCT(int clct, int endcap, int dPhiSign, int bits) const {
  assert( clct >= 0 && clct <= 10 && abs(endcap) == 1 && 
	  abs(dPhiSign) == 1 && (bits == 2 || bits == 3) );

  // Convention here: endcap == +/-1, dPhiSign = +/-1.  May need to change to match FW. - AWB 17.03.17
  int clct_ = 0;
  int sign_ = -1 * endcap * dPhiSign;  // CLCT bend is with dPhi in ME-, opposite in ME+

  if (clct < 2) {
    // std::cout << "\n\n*** In endcap " << endcap << ", CLCT = " << clct << std::endl;
    clct = 2;
  }

  // CLCT pattern can be converted into |bend| x sign as follows:
  // |bend| = (10 + (pattern % 2) - pattern) / 2
  //   * 10 --> 0, 9/8 --> 1, 7/6 --> 2, 5/4 --> 3, 3/2 --> 4, 0 indicates RPC hit
  //  sign  = ((pattern % 2) == 1 ? -1 : 1) * (endcap == 1 ? -1 : 1)   
  //   * In ME+, even CLCTs have negative sign, odd CLCTs have positive

  // For use in all 3- and 4-station modes (7, 11, 13, 14, 15)
  // Bends [-4, -3, -2] --> 0, [-1, 0] --> 1, [+1] --> 2, [+2, +3, +4] --> 3
  if (bits == 2) {
    assert(clct >= 2);
    switch (clct) {
    case 10: clct_ = 1;                  break;
    case  9: clct_ = (sign_ > 0 ? 1 : 2); break;
    case  8: clct_ = (sign_ < 0 ? 1 : 2); break;
    case  7: clct_ = (sign_ > 0 ? 0 : 3); break;
    case  6: clct_ = (sign_ < 0 ? 0 : 3); break;
    case  5: clct_ = (sign_ > 0 ? 0 : 3); break;
    case  4: clct_ = (sign_ < 0 ? 0 : 3); break;
    case  3: clct_ = (sign_ > 0 ? 0 : 3); break;
    case  2: clct_ = (sign_ < 0 ? 0 : 3); break;
    default: clct_ = 0;                   break;
    }
  } // End conditional: if (bits == 2)

  // For use in all 2-station modes (3, 5, 6, 9, 10, 12)
  // Bends [isRPC] --> 0, [-4, -3] --> 1, [-2] --> 2, [-1] --> 3, [0] --> 4, [+1] --> 5, [+2] --> 6, [+3, +4] --> 7
  else if (bits == 3) {
    assert(clct >= 2 || clct == 0);
    switch (clct) {
    case 10: clct_ = 4;                   break;
    case  9: clct_ = (sign_ > 0 ? 3 : 5); break;
    case  8: clct_ = (sign_ < 0 ? 3 : 5); break;
    case  7: clct_ = (sign_ > 0 ? 2 : 6); break;
    case  6: clct_ = (sign_ < 0 ? 2 : 6); break;
    case  5: clct_ = (sign_ > 0 ? 1 : 7); break;
    case  4: clct_ = (sign_ < 0 ? 1 : 7); break;
    case  3: clct_ = (sign_ > 0 ? 1 : 7); break;
    case  2: clct_ = (sign_ < 0 ? 1 : 7); break;
    case  0: clct_ = 0;                   break;
    default: clct_ = 0;                   break;
    }
  } // End conditional: else if (bits == 3)

  assert(clct_ >= 0 && clct_ < pow(2, bits));
  return clct_;
} // End function: int PtAssignmentEngineAux2017::getCLCT()


int PtAssignmentEngineAux2017::getdTheta(int dTheta, int bits) const {
  assert( bits == 2 || bits == 3 );

  int dTheta_ = -99;

  // For use in mode 15
  if (bits == 2) {
    if      (abs(dTheta) <= 1)
      dTheta_ = 2;
    else if (abs(dTheta) <= 2)
      dTheta_ = 1;
    else if (dTheta <= -3)
      dTheta_ = 0;
    else
      dTheta_ = 3;
  } // End conditional: if (bits == 2)

  // For use in all 2- and 3-station modes (all modes except 15)
  else if (bits == 3) {
    if      (dTheta <= -4)
      dTheta_ = 0;
    else if (dTheta == -3)
      dTheta_ = 1;
    else if (dTheta == -2)
      dTheta_ = 2;
    else if (dTheta == -1)
      dTheta_ = 3;
    else if (dTheta ==  0)
      dTheta_ = 4;
    else if (dTheta == +1)
      dTheta_ = 5;
    else if (dTheta == +2)
      dTheta_ = 6;
    else
      dTheta_ = 7;
  } // End conditional: if (bits == 3)

  assert(dTheta_ >= 0 && dTheta_ < pow(2, bits));
  return (dTheta_);
} // End function: int PtAssignmentEngineAux2017::getdTheta()


int PtAssignmentEngineAux2017::getTheta(int theta, int st1_ring2, int bits) const {
  assert( theta >= 5 && theta < 128 && 
	  (st1_ring2 == 0 || st1_ring2 == 1) && 
	  (bits == 4 || bits == 5) );

  int theta_ = -99;

  // For use in mode 15
  if (bits == 4) {
    if (st1_ring2 == 0) {
      // Should never fail with dTheta < 4 windows ... should change to using ME1 for track theta - AWB 17.03.17
      if (theta > 58) {
	std::cout << "\n\n*** Bizzare case of mode 15 track with ME1/1 LCT and track theta = " << theta << std::endl;
      }     
      theta_ = (std::min( std::max(theta, 5), 58) - 5) / 9;
    }
    else if (st1_ring2 == 1) {
      // Should rarely fail with dTheta < 4 windows ... should change to using ME1 for track theta - AWB 17.03.17
      if (theta < 44 || theta > 88) {
	std::cout << "\n\n*** Bizzare case of mode 15 track with ME1/2 LCT and track theta = " << theta << std::endl;
      }
      theta_ = ((std::min( std::max(theta, 44), 88) - 44) / 9) + 6;
    }
  } // End conditional: if (bits == 4)

  // For use in all 2- and 3-station modes (all modes except 15)
  else if (bits == 5) {
    if (st1_ring2 == 0) {
      theta_ = (std::max(theta, 1) - 1) / 4;
    }
    else if (st1_ring2 == 1) {
      theta_ = ((std::min(theta, 104) - 1) / 4) + 6;
    }
  } // End conditional: else if (bits == 5)

  assert(theta_ >= 0 && ((bits == 4 && theta_ <= 10) || (bits == 5 && theta_ < pow(2, bits))) );
  return (theta_);
} // End function: int PtAssignmentEngineAux2017::getTheta()


// // Need to re-check / verify this - AWB 17.03.17
// // front-rear LUTs
// // [sector[0]][station 0-4][chamber id]
// // chamber numbers start from 1, so add an extra low bit for invalid chamber = 0
// static const int FRLUT[2][5] = {
//   {0b0000000100100, 0b0000001011010, 0b0101010101010, 0b0010101010100, 0b0010101010100},
//   {0b0000000100100, 0b0000001011010, 0b0111010100100, 0b0000101011010, 0b0000101011010}
// };

// int PtAssignmentEngineAux2017::getFRLUT(int sector, int station, int chamber) const {
//   int bits = FRLUT[(sector-1)%2][station];
//   bool isFront = bits & (1<<chamber);
//   return isFront;
// }


// _____________________________________________________________________________
static const int GMT_eta_from_theta[128] = {
  239, 235, 233, 230, 227, 224, 222, 219, 217, 214, 212, 210, 207, 205, 203, 201,
  199, 197, 195, 193, 191, 189, 187, 186, 184, 182, 180, 179, 177, 176, 174, 172,
  171, 169, 168, 166, 165, 164, 162, 161, 160, 158, 157, 156, 154, 153, 152, 151,
  149, 148, 147, 146, 145, 143, 142, 141, 140, 139, 138, 137, 136, 135, 134, 133,
  132, 131, 130, 129, 128, 127, 126, 125, 124, 123, 122, 121, 120, 119, 118, 117,
  116, 116, 115, 114, 113, 112, 111, 110, 110, 109, 108, 107, 106, 106, 105, 104,
  103, 102, 102, 101, 100,  99,  99,  98,  97,  96,  96,  95,  94,  93,  93,  92,
   91,  91,  90,  89,  89,  88,  87,  87,  86,  85,  84,  84,  83,  83,  82,  81
};

int PtAssignmentEngineAux2017::getGMTPt(float pt) const {
  // compressed pt = pt*2 (scale) + 1 (pt = 0 is empty candidate)
  int gmt_pt = (pt * 2) + 1;
  gmt_pt = (gmt_pt > 511) ? 511 : gmt_pt;
  return gmt_pt;
}

float PtAssignmentEngineAux2017::getPtFromGMTPt(int gmt_pt) const {
  float pt = (gmt_pt <= 0) ?  0 : 0.5 * (gmt_pt-1);
  return pt;
}

int PtAssignmentEngineAux2017::getGMTPhi(int phi) const {
  // convert phi into gmt scale according to DN15-017
  // full scale is -16 to 100, or 116 values, covers range -10 to 62.5 deg
  // my internal ph scale is 0..5000, covers from -22 to 63.333 deg
  // converted to GMT scale it is from -35 to 95
  // bt_phi * 107.01/4096, equivalent to bt_phi * 6849/0x40000
  phi *= 6849;
  phi >>= 18; // divide by 0x40000
  phi -= 35;  // offset of -22 deg
  return phi;
}

int PtAssignmentEngineAux2017::getGMTPhiV2(int phi) const {
  // convert phi into gmt scale according to DN15-017
  phi *= 6991;
  phi >>= 18; // divide by 0x40000
  phi -= 35;  // offset of -22 deg
  return phi;
}

int PtAssignmentEngineAux2017::getGMTEta(int theta, int endcap) const {
  if (theta < 0)
    return 0;
  if (endcap == -1 && theta > 127)
    return -240;
  if (endcap == +1 && theta > 127)
    return 239;

  int eta = GMT_eta_from_theta[theta];
  if (endcap == -1)
    eta = -eta;
  return eta;
}

int PtAssignmentEngineAux2017::getGMTQuality(int mode, int theta) const {
  int quality = 0;
  if (theta > 87) {  // if (eta < 1.2)
    switch (mode) {
    case 15:  quality = 8;  break;
    case 14:  quality = 4;  break;
    case 13:  quality = 4;  break;
    case 12:  quality = 4;  break;
    case 11:  quality = 4;  break;
    default:  quality = 4;  break;
    }
  } else {
    switch (mode) {
    case 15:  quality = 12; break;
    case 14:  quality = 12; break;
    case 13:  quality = 12; break;
    case 12:  quality = 8;  break;
    case 11:  quality = 12; break;
    case 10:  quality = 8;  break;
    case 7:   quality = 8;  break;
    default:  quality = 4;  break;
    }
  }
  quality |= (mode & 3);
  return quality;
}

// std::pair<int,int> PtAssignmentEngineAux2017::getGMTCharge(int mode, const std::vector<int>& phidiffs) const {
//   // -1 = postive physical charge to match pdgId code (i.e. -13 is positive, anti-muon). +1 = negative physical charge.
//   // Also matches DN-2015/017 format for track finder --> uGMT interface format, where 0 indicates positive, 1 negative.
//   int emuCharge = 0;

//   // Note: sign_ph[0] == 1 in firmware actually translates to phidiffs[0] >= 0 (instead of phidiffs[0] > 0 in the old emulator)
//   // The effect needs to be checked

//   switch (mode) {
//   case 15:  // 1-2-3-4
//     if (phidiffs[0] >= 0)                         // 1-2 (should use > 0)
//       emuCharge = 1;
//     else if (phidiffs[0] == 0 && phidiffs[1] < 0) // 1-3
//       emuCharge = 1;
//     else if (phidiffs[1] == 0 && phidiffs[2] < 0) // 1-4
//       emuCharge = 1;
//     else
//       emuCharge = -1;
//     break;

//   default:
//     //emuCharge = -1;
//     emuCharge = 0;
//     break;
//   }

//   int charge = 0;
//   if (emuCharge == 1)
//     charge = 1;

//   int charge_valid = 1;
//   if (emuCharge == 0)
//     charge_valid = 0;
//   return std::make_pair(charge, charge_valid);
// }
