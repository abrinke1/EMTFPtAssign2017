
/////////////////////////////////////////////////////////
///         Simple macro to read EMTF NTuples         ///
///            Andrew Brinkerhoff 10.12.16            ///
///                                                   ///
///  TChain can be used to read multiple files.       ///
///  NTuple format in interface/PtLutInputBranches.hh ///
/////////////////////////////////////////////////////////

#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

const int MAX_EVT  = 10;   // Max number of events to process
const int PRT_EVT  =  1;   // Print every N events
const bool verbose = true; // Print information about muons, tracks, and hits

void ReadNTupleChain() {

  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);

  // List of input files
  std::vector<TString> in_file_names;
  TString store = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/ntuples";
  TString in_dir = "SingleMu_Pt1To1000_FlatRandomOneOverPt/EMTF_MuGun/170113_165434/0000";
  TString file_name;
  for (int i = 1; i < 99; i++) {
    file_name.Form("%s/%s/EMTF_MC_NTuple_SingleMu_noRPC_%d.root", store.Data(), in_dir.Data(), i);
    std::cout << "Adding file " << file_name.Data() << std::endl;
    in_file_names.push_back(file_name.Data());
    if (i*100000 > MAX_EVT) break; // ~100k events per file
  }

  // Open all input files
  for (UInt_t i = 0; i < in_file_names.size(); i++) {
    if ( !gSystem->AccessPathName(in_file_names.at(i)) )
      file_tmp = TFile::Open( in_file_names.at(i) ); // Check if file exists
    if (!file_tmp) {
      std::cout << "ERROR: could not open data file " << in_file_names.at(i) << std::endl;
      return;
    }
  }

  // Add trees from the input files to the TChain
  std::vector<TChain*> in_chains; 
  // Super-hacky ... but using "GetBranch" with a single chain with multiple files causes a segfault - AWB 19.01.16
  for (UInt_t i = 0; i < in_file_names.size(); i++) {
    TChain *tmp_chain = new TChain("ntuple/tree");
    tmp_chain->Add( in_file_names.at(i) );
    in_chains.push_back(tmp_chain);
  }

  std::cout << "\n******* About to loop over chains *******" << std::endl;
  UInt_t iEvt = 0;
  for (int iCh = 0; iCh < in_chains.size(); iCh++) {
    TChain *in_chain = in_chains.at(iCh);
    
    // Get branches from the chain
    TBranch *muon_br  = in_chain->GetBranch("muon");
    TBranch *hit_br   = in_chain->GetBranch("hit");
    TBranch *track_br = in_chain->GetBranch("track");
    
    std::cout << "\n******* About to enter the event loop for chain " << iCh+1 << " *******" << std::endl;
    for (UInt_t jEvt = 0; jEvt < in_chain->GetEntries(); jEvt++) {
      
      if (iEvt > MAX_EVT) break;
      if ( (iEvt % PRT_EVT) == 0 ) {
	std::cout << "\n*********************" << std::endl;
	std::cout << "Looking at event " << iEvt <<  std::endl;
	std::cout << "*********************" << std::endl;
      }
      iEvt += 1;
      
      in_chain->GetEntry(jEvt);
      
      UInt_t nMuons  = (muon_br->GetLeaf("nMuons"))->GetValue();
      UInt_t nHits   = (hit_br->GetLeaf("nHits"))->GetValue();
      UInt_t nTracks = (track_br->GetLeaf("nTracks"))->GetValue();
      
      if (verbose) std::cout << "\nnMuons = " << nMuons << std::endl;
      if (verbose) std::cout << "nHits = " << nHits << std::endl;
      if (verbose) std::cout << "nTracks = " << nTracks << std::endl;
      
      // Print info for GEN level muons
      for (UInt_t iMu = 0; iMu < nMuons; iMu++) {
	Double_t mu_pt  = (muon_br->GetLeaf("pt"))->GetValue(iMu);
	Double_t mu_eta = (muon_br->GetLeaf("eta"))->GetValue(iMu);
	Double_t mu_phi = (muon_br->GetLeaf("phi"))->GetValue(iMu);
	if (verbose) std::cout << "\nMuon " << iMu+1 << " has pT = " << mu_pt << ", eta = " << mu_eta << ", phi = " << mu_phi << std::endl;
	
	// Print info for hits in the same endcap
	for (UInt_t iHit = 0; iHit < nHits; iHit++) {
	  Int_t hit_station = int( (hit_br->GetLeaf("station"))->GetValue(iHit) );
	  Double_t hit_eta  = (hit_br->GetLeaf("eta"))->GetValue(iHit);
	  Double_t hit_phi  = (hit_br->GetLeaf("phi"))->GetValue(iHit);
	  
	  if (hit_station < 0) continue;
	  if ( (hit_eta > 0) != (mu_eta > 0) ) continue;
	  if (verbose) std::cout << "    * Station " << hit_station << " hit has eta = " << hit_eta << ", phi = " << hit_phi << std::endl;
	} // End loop: for (UInt_t iHit = 0; iHit < nHits; iHit++)
	
	// Print info for EMTF tracks in the same endcap
	for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++) {
	  
	  Int_t trk_sector = int( (track_br->GetLeaf("sector"))->GetValue(iTrk) );
	  Int_t trk_mode   = int( (track_br->GetLeaf("mode"))->GetValue(iTrk) );
	  Double_t trk_pt  = (track_br->GetLeaf("pt"))->GetValue(iTrk);
	  Double_t trk_eta = (track_br->GetLeaf("eta"))->GetValue(iTrk);
	  Double_t trk_phi = (track_br->GetLeaf("phi"))->GetValue(iTrk);
	  
	  if (trk_sector < 0) continue;
	  if ( (trk_eta > 0) != (mu_eta > 0) ) continue;
	  if (verbose) std::cout << "  * Sector " << trk_sector << " track has pT = " << trk_pt << ", eta = " << trk_eta 
				 << ", phi = " << trk_phi << ", mode = " << trk_mode << std::endl;
	  // Print info for hits in the EMTF tracks
	  for (UInt_t iTrkHit = 4*iTrk; iTrkHit < (4*iTrk) + 4; iTrkHit++) {
	    Int_t hit_station = int( (track_br->GetLeaf("hit_station"))->GetValue(iTrkHit) );
	    Double_t hit_eta  = (track_br->GetLeaf("hit_eta"))->GetValue(iTrkHit);
	    Double_t hit_phi  = (track_br->GetLeaf("hit_phi"))->GetValue(iTrkHit);
	    
	    if (hit_station < 0) continue;
	    if (verbose) std::cout << "    * Station " << hit_station << " hit has eta = " << hit_eta << ", phi = " << hit_phi << std::endl;
	  } // End loop: for (UInt_t iTrkHit = iTrk*4; iTrkHit < (iTrk*4) + 4; iTrkHit++)
	} // End loop: for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++)

      } // End loop: for (UInt_t iMu = 0; iMu < nMuons; iMu++)
      
    } // End loop: for (UInt_t jEvt = 0; jEvt < in_chain->GetEntries(); jEvt++)
    std::cout << "\n******* Leaving the event loop for chain " << iCh+1 << " *******" << std::endl;

  } // End loop: for (int iCh = 0; iCh < in_chains.size(); iCh++)
  std::cout << "\n******* Leaving the loop over chains *******" << std::endl;

  std::cout << "\nExiting ReadNTupleChain()\n";
  
} // End void ReadNTupleChain()
  
