#ifndef L1TMuonEndCap_PtAssignmentEngineAux2017_hh
#define L1TMuonEndCap_PtAssignmentEngineAux2017_hh

#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>


class PtAssignmentEngineAux2017 {
public:
  // // Functions for pT assignment
  // const int (*getModeVariables() const)[6];

  int getNLBdPhi(int dPhi, int bits=7, int max=512) const;

  int getNLBdPhiBin(int dPhi, int bits=7, int max=512) const;

  int getdPhiFromBin(int dPhiBin, int bits=7, int max=512) const;

  int getCLCT(int clct, int endcap, int dPhiSign, int bits=3) const;

  int getdTheta(int dTheta, int bits=3) const;

  int getTheta(int theta, int ring2, int bits=5) const;

  // Need to re-check / verify this - AWB 17.03.17
  // int getFRLUT(int sector, int station, int chamber) const;

  // Functions for GMT quantities
  int getGMTPt(float pt) const;

  float getPtFromGMTPt(int gmt_pt) const;

  int getGMTPhi(int phi) const;
  int getGMTPhiV2(int phi) const;

  int getGMTEta(int theta, int endcap) const;

  int getGMTQuality(int mode, int theta) const;

  std::pair<int,int> getGMTCharge(int mode, const std::vector<int>& phidiffs) const;

};

#endif
