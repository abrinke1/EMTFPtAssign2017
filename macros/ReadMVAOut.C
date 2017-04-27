
/////////////////////////////////////////////////////////
///          Macro to compute pT resolution           ///
///            Andrew Brinkerhoff 10.12.16            ///
///                                                   ///
/////////////////////////////////////////////////////////

#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "ReadMVAOut.h"

void ReadMVAOut() {

  std::cout<<"foo = "<<foo<<endl;
  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);

  // List of input files
  std::vector<TString> in_file_names;
  //in_file_names.push_back("/afs/cern.ch/user/a/abrinke1/TMVA/EMTFPtAssign2017/PtRegression_AWB_v0_16_12_09.root");
  in_file_names.push_back(str_file_location);

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
  TChain *train_chain = new TChain("dataset/TrainTree");
  TChain *test_chain = new TChain("dataset/TestTree");
  for (UInt_t i = 0; i < in_file_names.size(); i++) {
    train_chain->Add( in_file_names.at(i) );
    test_chain->Add( in_file_names.at(i) );
  }

  // Get branches from the chains
  Float_t EMTF_pt_train_br;
  Float_t GEN_pt_train_br;
  Float_t inv_GEN_pt_train_br;
  Float_t BDTG_default_train_br;

  Float_t EMTF_pt_test_br;
  Float_t GEN_pt_test_br;
  Float_t inv_GEN_pt_test_br;
  Float_t BDTG_default_test_br;

  train_chain->SetBranchAddress("EMTF_pt", &EMTF_pt_train_br);
  train_chain->SetBranchAddress("GEN_pt", &GEN_pt_train_br);
  train_chain->SetBranchAddress("inv_GEN_pt", &inv_GEN_pt_train_br);
  train_chain->SetBranchAddress("BDTG_default", &BDTG_default_train_br);

  test_chain->SetBranchAddress("EMTF_pt", &EMTF_pt_test_br);
  test_chain->SetBranchAddress("GEN_pt", &GEN_pt_test_br);
  test_chain->SetBranchAddress("inv_GEN_pt", &inv_GEN_pt_test_br);
  test_chain->SetBranchAddress("BDTG_default", &BDTG_default_test_br);


  std::cout << "\n******* About to enter the train event loop *******" << std::endl;
  for (UInt_t iEvt = 0; iEvt < train_chain->GetEntries(); iEvt++) {

    if (iEvt > 10) break;
    if ( (iEvt % 1) == 0 )
      std::cout << "\n*** Looking at event " << iEvt << " ***" <<  std::endl;

    train_chain->GetEntry(iEvt);

    std::cout << "EMTF_pt = " << EMTF_pt_train_br << std::endl;
    std::cout << "GEN_pt = " << GEN_pt_train_br << std::endl;
    std::cout << "inv_GEN_pt_train = " << inv_GEN_pt_train_br << std::endl;
    std::cout << "BDTG_default_train = " << BDTG_default_train_br << std::endl;

  } // End loop: for (UInt_t iEvt = 0; iEvt < train_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the train event loop *******" << std::endl;


  std::cout << "\n******* About to enter the test event loop *******" << std::endl;
  for (UInt_t iEvt = 0; iEvt < test_chain->GetEntries(); iEvt++) {

    if (iEvt > 10) break;
    if ( (iEvt % 1) == 0 ) 
      std::cout << "\n*** Looking at event " << iEvt << " ***" <<  std::endl;

    test_chain->GetEntry(iEvt);

    std::cout << "EMTF_pt = " << EMTF_pt_test_br << std::endl;
    std::cout << "GEN_pt = " << GEN_pt_test_br << std::endl;
    std::cout << "inv_GEN_pt_test = " << inv_GEN_pt_test_br << std::endl;
    std::cout << "BDTG_default_test = " << BDTG_default_test_br << std::endl;

  } // End loop: for (UInt_t iEvt = 0; iEvt < test_chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the test event loop *******" << std::endl;

  std::cout << "\nExiting ReadMVAOut()\n";

} // End void ReadMVAOut()
