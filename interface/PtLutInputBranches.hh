
// Defines structs for input NTuples - AWB 09.12.16
// Copied from /afs/cern.ch/user/a/abrinke1/JiaFuEmulator/CMSSW_8_0_24/src/EMTFAnalyzer/NTupleMaker/interface/PtLutInputBranches.hh

#ifndef PtLutInputBranches_hh
#define PtLutInputBranches_hh

// GEN particles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// EMTF classes
#include "DataFormatsSep2016/L1TMuon/interface/EMTFHitExtra.h"
#include "DataFormatsSep2016/L1TMuon/interface/EMTFTrackExtra.h"

// Helpful tools
#include "EMTFAnalyzer/NTupleMaker/interface/HelperFunctions.hh"


// Size of branches
const unsigned int N_GEN =  2;
const unsigned int N_HIT = 24;
const unsigned int N_TRK =  4;

// Default fill values
const int   DINT = -999;
const float DFLT = -999.0;

// Structs for output tree
typedef struct {
  int nMuons;
  float pt[N_GEN], eta[N_GEN], theta[N_GEN], phi[N_GEN];
  int charge[N_GEN];
  void Initialize();
  void Fill(unsigned int i, reco::GenParticle genMuon);
} GenMuonBranch;

typedef struct {
  int nHits;
  float eta[N_HIT], theta[N_HIT], phi[N_HIT], phi_loc[N_HIT];
  int eta_int[N_HIT], theta_int[N_HIT], phi_int[N_HIT];
  int endcap[N_HIT], sector[N_HIT], sector_index[N_HIT], station[N_HIT], ring[N_HIT];
  int CSC_ID[N_HIT], chamber[N_HIT], FR[N_HIT], pattern[N_HIT];
  int roll[N_HIT], subsector[N_HIT], isRPC[N_HIT], vetoed[N_HIT];
  void Initialize();
  void Fill(unsigned int i, L1TMuonEndCap::EMTFHitExtra emtfHit);
} EMTFHitBranch;

typedef struct {
  int nTracks;
  float pt[N_TRK], eta[N_TRK], theta[N_TRK], phi[N_TRK], phi_loc[N_TRK];
  int pt_int[N_TRK], eta_int[N_TRK], theta_int[N_TRK], phi_int[N_TRK];
  int endcap[N_TRK], sector[N_TRK], sector_index[N_TRK], mode[N_TRK], charge[N_TRK];

  int nHits[N_TRK], nRPC[N_TRK];
  float hit_eta[N_TRK][4], hit_theta[N_TRK][4], hit_phi[N_TRK][4], hit_phi_loc[N_TRK][4];
  int hit_eta_int[N_TRK][4], hit_theta_int[N_TRK][4], hit_phi_int[N_TRK][4];
  int hit_endcap[N_TRK][4], hit_sector[N_TRK][4], hit_sector_index[N_TRK][4], hit_station[N_TRK][4], hit_ring[N_TRK][4];
  int hit_CSC_ID[N_TRK][4], hit_chamber[N_TRK][4], hit_FR[N_TRK][4], hit_pattern[N_TRK][4];
  int hit_roll[N_TRK][4], hit_subsector[N_TRK][4], hit_isRPC[N_TRK][4], hit_vetoed[N_TRK][4];

  void Initialize();
  void Fill(unsigned int i, L1TMuonEndCap::EMTFTrackExtra emtfTrk);
} EMTFTrackBranch;


#endif /* define PtLutInputBranches_hh */
