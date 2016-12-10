
/////////////////////////////////////////////////////////
///         Simple macro to read EMTF NTuples         ///
///            Andrew Brinkerhoff 10.12.16            ///
///                                                   ///
///  TTreeReader pre-loads the branches to be read    ///
///  ONLY WORKS WHEN BRANCHES HAVE SIMPLE STRUCTURE   ///
///  e.g. a branch which is an int or float can be 
///  read with TTreeReaderValue; arrays of ints or 
//   floats with TTreeReaderArray.
///  NTuple format in interface/PtLutInputBranches.hh ///
/////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void ReadNTupleReader() {

  TString in_file_name = "/afs/cern.ch/work/a/abrinke1/public/EMTF/Analyzer/hiPt/EMTF_NTuple_highPt200MuonSkim_emtfStage2Digis_2016BCD.root";
  
  TFile *in_file = TFile::Open(in_file_name);
  
  if (in_file == 0) {
    // If we cannot open the file, print an error message and return immediatly
    std::cout << "Error: cannot open " << in_file_name << std::endl;
    return;
  }

  // Create a tree reader (of type Int_t) on the branch "fEventSize"
  TTreeReader reader("ntuple/tree", in_file);

  ////////////////////////////////
  // Pre-book branches in the tree
  ////////////////////////////////

  // RECO muons
  TTreeReaderValue<Int_t>   _numRecoMuons (reader, "numRecoMuons");
  TTreeReaderArray<Float_t> _recoPt       (reader, "recoPt");
  TTreeReaderArray<Float_t> _recoEta      (reader, "recoEta");
  TTreeReaderArray<Float_t> _recoPhi      (reader, "recoPhi");

  // Unpacked EMTF
  TTreeReaderValue<Int_t>   _numUnpTrks  (reader, "numUnpTrks");
  TTreeReaderArray<Float_t> _unp_trkPt   (reader, "unp_trkPt");
  TTreeReaderArray<Float_t> _unp_trkEta  (reader, "unp_trkEta");
  TTreeReaderArray<Float_t> _unp_trkPhi  (reader, "unp_trkPhi");
  TTreeReaderArray<Int_t>   _unp_trkMode (reader, "unp_trkMode");
  
  // Re-emulated EMTF
  TTreeReaderValue<Int_t>   _numTrks  (reader, "numTrks");
  TTreeReaderArray<Float_t> _trkPt    (reader, "trkPt");
  TTreeReaderArray<Float_t> _trkEta   (reader, "trkEta");
  TTreeReaderArray<Float_t> _trkPhi   (reader, "trkPhi");
  TTreeReaderArray<Int_t>   _trkMode  (reader, "trkMode");
  
  // Re-emulated (?) CSCTF
  TTreeReaderValue<Int_t>   _numLegTrks  (reader, "numLegTrks");
  TTreeReaderArray<Float_t> _leg_trkPt   (reader, "leg_trkPt");
  TTreeReaderArray<Float_t> _leg_trkEta  (reader, "leg_trkEta");
  TTreeReaderArray<Float_t> _leg_trkPhi  (reader, "leg_trkPhi");
  TTreeReaderArray<Int_t>   _leg_trkMode (reader, "leg_trkMode");

  
  std::cout << "\n******* About to enter the event loop *******" << std::endl;
  UInt_t iEvt = 0;
  while (reader.Next()) {
    iEvt += 1;
    if (iEvt > 10) break;
    if ( (iEvt % 1) == 0 ) {
      std::cout << "\n*********************" << std::endl;
      std::cout << "Looking at event " << iEvt <<  std::endl;
      std::cout << "*********************" << std::endl;
    }

    UInt_t numRecoMuons = (*_numRecoMuons);
    UInt_t numUnpTrks   = (*_numUnpTrks);
    UInt_t numTrks      = (*_numTrks);
    UInt_t numLegTrks   = (*_numLegTrks);

    std::cout << "\nnumRecoMuons = " << numRecoMuons << std::endl;
    std::cout << "numUnpTrks = " << numUnpTrks << std::endl;
    std::cout << "numTrks = " << numTrks << std::endl;
    std::cout << "numLegTrks = " << numLegTrks << std::endl;

    // Print info for RECO level muons
    for (UInt_t iReco = 0; iReco < numRecoMuons; iReco++) {
      std::cout << "\nRECO muon " << iReco+1 << " has pT = " << _recoPt[iReco] << ", eta = " << _recoEta[iReco] << ", phi = " << _recoPhi[iReco] << std::endl;

      // Print info for unpacked EMTF tracks in the same endcap
      for (UInt_t iUnp = 0; iUnp < numUnpTrks; iUnp++) {
    	if ( (_unp_trkEta[iUnp] > 0) != (_recoEta[iReco] > 0) ) continue;
	std::cout << "  * Unpacked EMTF track has pT = " << _unp_trkPt[iUnp] << ", eta = " << _unp_trkEta[iUnp] 
		  << ", phi = " << _unp_trkPhi[iUnp] << ", mode = " << _unp_trkMode[iUnp] << std::endl;
      }
      
      // Print info for unpacked EMTF tracks in the same endcap
      for (UInt_t iTrk = 0; iTrk < numTrks; iTrk++) {
    	if ( (_trkEta[iTrk] > 0) != (_recoEta[iReco] > 0) ) continue;
	std::cout << "  * Emulated EMTF track has pT = " << _trkPt[iTrk] << ", eta = " << _trkEta[iTrk] 
		  << ", phi = " << _trkPhi[iTrk] << ", mode = " << _trkMode[iTrk] << std::endl;
      }

      // Print info for emulated (?) CSCTF tracks in the same endcap
      for (UInt_t iLeg = 0; iLeg < numLegTrks; iLeg++) {
	// CSCTF track eta not signed in these trees
    	// if ( (_leg_trkEta[iLeg] > 0) != (_recoEta[iReco] > 0) ) continue;
	std::cout << "  * Legacy CSCTF track has pT = " << _leg_trkPt[iLeg] << ", eta = +/-" << _leg_trkEta[iLeg] 
		  << ", phi = " << _leg_trkPhi[iLeg] << ", mode = " << _leg_trkMode[iLeg] << std::endl;
      }

    } // End loop: for (UInt_t iReco = 0; iReco < numRecoMuons; iReco++)

  } // End loop: for (UInt_t iEvt = 0; iEvt < chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the event loop *******" << std::endl;

  std::cout << "\nExiting ReadNTupleReader()\n";

} // End void ReadNTupleReader()
