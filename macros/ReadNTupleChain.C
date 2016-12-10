
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

void ReadNTupleChain() {

  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);

  // List of input files
  std::vector<TString> in_file_names;
  in_file_names.push_back("/afs/cern.ch/work/a/abrinke1/public/EMTF/Analyzer/ntuples/EMTF_MC_NTuple_SingleMu_noRPC_300k.root");

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
  TChain *in_chain = new TChain("ntuple/tree");
  for (UInt_t i = 0; i < in_file_names.size(); i++)
    in_chain->Add( in_file_names.at(i) );

  // Get branches from the chain
  TBranch *muon_br  = in_chain->GetBranch("muon");
  TBranch *hit_br   = in_chain->GetBranch("hit");
  TBranch *track_br = in_chain->GetBranch("track");


  std::cout << "\n******* About to enter the event loop *******" << std::endl;
  for (UInt_t iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++) {

    if (iEvt > 10) break;
    if ( (iEvt % 1) == 0 ) {
      std::cout << "\n*********************" << std::endl;
      std::cout << "Looking at event " << iEvt <<  std::endl;
      std::cout << "*********************" << std::endl;
    }

    in_chain->GetEntry(iEvt);

    UInt_t nMuons  = (muon_br->GetLeaf("nMuons"))->GetValue();
    UInt_t nHits   = (hit_br->GetLeaf("nHits"))->GetValue();
    UInt_t nTracks = (track_br->GetLeaf("nTracks"))->GetValue();

      std::cout << "\nnMuons = " << nMuons << std::endl;
      std::cout << "nHits = " << nHits << std::endl;
      std::cout << "nTracks = " << nTracks << std::endl;

    
    // Print info for GEN level muons
    for (UInt_t iMu = 0; iMu < nMuons; iMu++) {
      Double_t mu_pt  = (muon_br->GetLeaf("pt"))->GetValue(iMu);
      Double_t mu_eta = (muon_br->GetLeaf("eta"))->GetValue(iMu);
      Double_t mu_phi = (muon_br->GetLeaf("phi"))->GetValue(iMu);
      std::cout << "\nMuon " << iMu+1 << " has pT = " << mu_pt << ", eta = " << mu_eta << ", phi = " << mu_phi << std::endl;

      // Print info for hits in the same endcap
      for (UInt_t iHit = 0; iHit < nHits; iHit++) {
      	Int_t hit_station = int( (hit_br->GetLeaf("station"))->GetValue(iHit) );
      	Double_t hit_eta  = (hit_br->GetLeaf("eta"))->GetValue(iHit);
      	Double_t hit_phi  = (hit_br->GetLeaf("phi"))->GetValue(iHit);
	
      	if (hit_station < 0) continue;
      	if ( (hit_eta > 0) != (mu_eta > 0) ) continue;
      	std::cout << "    * Station " << hit_station << " hit has eta = " << hit_eta << ", phi = " << hit_phi << std::endl;
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
	std::cout << "  * Sector " << trk_sector << " track has pT = " << trk_pt << ", eta = " << trk_eta << ", phi = " << trk_phi << ", mode = " << trk_mode << std::endl;

	// Print info for hits in the EMTF tracks
	for (UInt_t iTrkHit = 4*iTrk; iTrkHit < (4*iTrk) + 4; iTrkHit++) {
	  Int_t hit_station = int( (track_br->GetLeaf("hit_station"))->GetValue(iTrkHit) );
	  Double_t hit_eta  = (track_br->GetLeaf("hit_eta"))->GetValue(iTrkHit);
	  Double_t hit_phi  = (track_br->GetLeaf("hit_phi"))->GetValue(iTrkHit);
	  
	  if (hit_station < 0) continue;
	  std::cout << "    * Station " << hit_station << " hit has eta = " << hit_eta << ", phi = " << hit_phi << std::endl;
	} // End loop: for (UInt_t iTrkHit = iTrk*4; iTrkHit < (iTrk*4) + 4; iTrkHit++)
      } // End loop: for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++)

    } // End loop: for (UInt_t iMu = 0; iMu < nMuons; iMu++)

  } // End loop: for (UInt_t iEvt = 0; iEvt < in_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the event loop *******" << std::endl;

  std::cout << "\nExiting ReadNTupleChain()\n";

} // End void ReadNTupleChain()
