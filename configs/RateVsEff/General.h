
// *** User *** ///
const TString USER = "AWB";  // Settings applied in User.h

// *** Events to process *** //
const int MAX_EVT    =      -1;  // Number of events to process
const int REPORT_EVT =  100000;  // Report every Nth event during processing

// *** Muon / track settings *** //
const double ETAMIN  = 1.20;  // Minimum GEN / trigger eta to consider
const double ETAMAX  = 2.40;  // Maximum GEN / trigger eta to consider

// *** Rate plot settings *** //
const double PTMIN   =  0.0;  // Minimum pT
const double PTMAX   = 50.0;  // Maximum pT
const int    PTBINS  =   50;  // Number of bins in rate plot

const int    PTDIVS  =   10;  // Number of trigger pT bins per GEN pT bin, before scaling
const double BIT     = 0.1 * (PTMAX - PTMIN) / (PTBINS * PTDIVS);  // Small pT offset increment
