
// *** User *** ///
const TString USER = "AWB";  // Settings applied in User.h

// *** Events to process *** //
const int MAX_EVT    =  2000000;  // Number of MC events to process          (4M default for mode 15, 12M available)
const int MAX_TR     =  1000000;  // Number of MC events to use for training (2M default for mode 15)
const int REPORT_EVT =    10000;  // Report every Nth event during processing
const int MAX_ZB_FIL =       50;  // Number of ZeroBias files to include (~200 CSC-only, ~30 with RPC)

/* // // ***** Test settings ***** // // */
/* const int MAX_EVT    =  40000;  // Number of MC events to process */
/* const int MAX_TR     =  20000;  // Number of MC events to use for training */
/* const int REPORT_EVT =   1000;  // Report every Nth event during processing */
/* const int MAX_ZB_FIL =      3;  // Number of ZeroBias files to include */

// *** Track-building settings *** //
const int MODE      = 15;     // Track mode to build - settings applied in Modes.h
const bool USE_RPC  = true;   // Use RPC hits in track-building       
const bool BIT_COMP = true;   // Use bit-compressed versions of input variables
const std::vector<int> CSC_MASK = {};  // Mask CSC LCTs in these stations
const std::vector<int> RPC_MASK = {};  // Mask RPC hits in these stations

// *** Target and event weights *** //
// Choose "pt", "invPt", "logPt", and/or "charge"
const std::vector<TString> TARG_VARS = {"invPt"}; // Default "invPt"
// Choose "no", "invPt", and/or "invPtSq"
const std::vector<TString> EVT_WGTS  = {"invPt"}; // Default "invPt" 

// *** EMTF tracks *** //
const bool REQ_EMTF     = true;  // Require that an EMTF muon be matched to the GEN muon
const bool USE_EMTF_CSC = true;  // Only use CSC LCTs from the original EMTF track
const std::vector<int> EMTF_MODES = {15, 14, 13, 12, 11, 10, 9, 7, 6, 5, 3};  // Modes of saved EMTF tracks

// *** Output data options *** //
const bool SPEC_VARS = true;  // When generating final XMLs, set to "false" to leave out spectators

// *** High-pT muons *** //
const double PTMIN_TR =    1.;  // Minimum GEN pT for training
const double PTMAX_TR =  256.;  // Maximum GEN pT for training
const double PTMAX_TRG = 128.;  // Maximum trigger pT assigned
const bool CLEAN_HI_PT = true;  // Remove showering high-pT mode 15 tracks from training
