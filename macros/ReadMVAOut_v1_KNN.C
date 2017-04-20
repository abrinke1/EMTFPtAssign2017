
/////////////////////////////////////////////////////////
///          Macro to compute rate and efficiency     ///
///            Wei Shi         03.24.17               ///
///             Original use for KNN algorithm        ///
/////////////////////////////////////////////////////////

#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include <vector>
#include "stdio.h"
#include "math.h"
#include "TMath.h"
#include "TGraph.h"

void ReadMVAOut_v1_KNN() {

  // Initialize empty file to access each file in the list
  TFile *file_tmp(0);

  // List of input files
  std::vector<TString> in_file_names;
  in_file_names.push_back("/afs/cern.ch/YOURFILEDIRECTORY");

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
  TChain *train_chain = new TChain("f_0x0000110d_0x2_invPt/TrainTree");
  TChain *test_chain = new TChain("f_0x0000110d_0x2_invPt/TestTree");
  for (UInt_t i = 0; i < in_file_names.size(); i++) {
    train_chain->Add( in_file_names.at(i) );
    test_chain->Add( in_file_names.at(i) );
  }
    
    //bins for trigger efficiency
    double trigger_Cut[59]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,55,60,65,70,75,80,85,90,95};//for cut
    double trigger_Cut_scaled[59]={0};//scaled so KNN reach 90% at each trigger cut
    double EMTF_trigger_Cut_scaled[59]={0};//scale EMTF to 90%
    
    double GEN_pT[60]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,55,60,65,70,75,80,85,90,95,100};//for counting purpose only
    double count[59]={0};//total count in each bin
    
    double pass_count[59][59]={0};//KNN predict
    double EMTF_pass_count[59][59]={0};//current EMTF
    double efficiency[59][59]={0};//KNN efficiency
    double EMTF_efficiency[59][59]={0};//EMTF efficiency
    
    //choose special trigger cuts to evaluate how scale factor change
    double trigger_Cut_special[6] = {8,16,24,32,40,48};
    
    double pass_count_special[6][59]={0};//KNN predict
    double EMTF_pass_count_special[6][59]={0};//current EMTF
    double efficiency_special[6][59]={0};//KNN efficiency
    double EMTF_efficiency_special[6][59]={0};//EMTF efficiency

    
  // Get branches from the chains
    Long64_t train_events=0;
    Long64_t train_events_1_8=0;
    Long64_t train_events_8_30=0;
    Long64_t train_events_30_120=0;
    Long64_t train_events_120_1000=0;
    
    Long64_t test_events=0;
    Long64_t test_events_1_8=0;
    Long64_t test_events_8_30=0;
    Long64_t test_events_30_120=0;
    Long64_t test_events_120_1000=0;
    
    //KNN
    Double_t Deviation_test=0;
    Double_t Standard_Deviation_test=0;
    Double_t Deviation_train=0;
    Double_t Standard_Deviation_train=0;
    //1-8GeV
    Double_t Deviation_test_1_8=0;
    Double_t Standard_Deviation_test_1_8=0;
    Double_t Deviation_train_1_8=0;
    Double_t Standard_Deviation_train_1_8=0;
    //8-30GeV
    Double_t Deviation_test_8_30=0;
    Double_t Standard_Deviation_test_8_30=0;
    Double_t Deviation_train_8_30=0;
    Double_t Standard_Deviation_train_8_30=0;
    //30-120GeV
    Double_t Deviation_test_30_120=0;
    Double_t Standard_Deviation_test_30_120=0;
    Double_t Deviation_train_30_120=0;
    Double_t Standard_Deviation_train_30_120=0;
    //120-1000GeV
    Double_t Deviation_test_120_1000=0;
    Double_t Standard_Deviation_test_120_1000=0;
    Double_t Deviation_train_120_1000=0;
    Double_t Standard_Deviation_train_120_1000=0;
    
    //EMTF SD for comparision
    Double_t EMTF_Deviation_test=0;
    Double_t EMTF_Standard_Deviation_test=0;
    Double_t EMTF_Deviation_train=0;
    Double_t EMTF_Standard_Deviation_train=0;
    
    Double_t EMTF_Deviation_test_1_8=0;
    Double_t EMTF_Standard_Deviation_test_1_8=0;
    Double_t EMTF_Deviation_train_1_8=0;
    Double_t EMTF_Standard_Deviation_train_1_8=0;
    
    Double_t EMTF_Deviation_test_8_30=0;
    Double_t EMTF_Standard_Deviation_test_8_30=0;
    Double_t EMTF_Deviation_train_8_30=0;
    Double_t EMTF_Standard_Deviation_train_8_30=0;
    
    Double_t EMTF_Deviation_test_30_120=0;
    Double_t EMTF_Standard_Deviation_test_30_120=0;
    Double_t EMTF_Deviation_train_30_120=0;
    Double_t EMTF_Standard_Deviation_train_30_120=0;
    
    Double_t EMTF_Deviation_test_120_1000=0;
    Double_t EMTF_Standard_Deviation_test_120_1000=0;
    Double_t EMTF_Deviation_train_120_1000=0;
    Double_t EMTF_Standard_Deviation_train_120_1000=0;
    
    Float_t dPhiSum3_train_br;
    Int_t FR_1_train_br;
    Float_t dTh_13_train_br;
    Int_t bend_1_train_br;
    Float_t theta_train_br;
    Float_t dPhi_12_train_br;
    Float_t dPhi_23_train_br;
    Float_t dPhi_34_train_br;
    Float_t EMTF_pt_train_br;
    Float_t GEN_pt_train_br;
    Float_t inv_EMTF_pt_train_br;
    Float_t inv_GEN_pt_train_br;
    Float_t KNN_train_br;
  
    Float_t dPhiSum3_test_br;
    Int_t FR_1_test_br;
    Float_t dTh_13_test_br;
    Int_t bend_1_test_br;
    Float_t theta_test_br;
    Float_t dPhi_12_test_br;
    Float_t dPhi_23_test_br;
    Float_t dPhi_34_test_br;
    Float_t EMTF_pt_test_br;
    Float_t GEN_pt_test_br;
    Float_t inv_EMTF_pt_test_br;
    Float_t inv_GEN_pt_test_br;
    Float_t KNN_test_br;
    //special for identifying ZeroBias data in test
    Float_t GEN_eta_test_br;
    Int_t GEN_charge_test_br;
    
    //train
    //train_chain->SetBranchAddress("dPhiSum3", &dPhiSum3_train_br);
    train_chain->SetBranchAddress("FR_1", &FR_1_train_br);
    //train_chain->SetBranchAddress("dTh_13", &dTh_13_train_br);
    train_chain->SetBranchAddress("bend_1", &bend_1_train_br);
    train_chain->SetBranchAddress("theta", &theta_train_br);
    train_chain->SetBranchAddress("dPhi_12", &dPhi_12_train_br);
    train_chain->SetBranchAddress("dPhi_23", &dPhi_23_train_br);
    //train_chain->SetBranchAddress("dPhi_34", &dPhi_34_train_br);
    train_chain->SetBranchAddress("EMTF_pt", &EMTF_pt_train_br);
    train_chain->SetBranchAddress("GEN_pt", &GEN_pt_train_br);
    train_chain->SetBranchAddress("inv_GEN_pt", &inv_GEN_pt_train_br);
    train_chain->SetBranchAddress("inv_EMTF_pt", &inv_EMTF_pt_train_br);
    train_chain->SetBranchAddress("KNN", &KNN_train_br);

    //test
    //test_chain->SetBranchAddress("dPhiSum3", &dPhiSum3_test_br);
    test_chain->SetBranchAddress("FR_1", &FR_1_test_br);
    //test_chain->SetBranchAddress("dTh_13", &dTh_13_test_br);
    test_chain->SetBranchAddress("bend_1", &bend_1_test_br);
    test_chain->SetBranchAddress("theta", &theta_test_br);
    test_chain->SetBranchAddress("dPhi_12", &dPhi_12_test_br);
    test_chain->SetBranchAddress("dPhi_23", &dPhi_23_test_br);
    //test_chain->SetBranchAddress("dPhi_34", &dPhi_34_test_br);
    test_chain->SetBranchAddress("EMTF_pt", &EMTF_pt_test_br);
    test_chain->SetBranchAddress("GEN_pt", &GEN_pt_test_br);
    test_chain->SetBranchAddress("inv_EMTF_pt", &inv_EMTF_pt_test_br);
    test_chain->SetBranchAddress("inv_GEN_pt", &inv_GEN_pt_test_br);
    test_chain->SetBranchAddress("KNN", &KNN_test_br);
    test_chain->SetBranchAddress("GEN_eta", &GEN_eta_test_br);
    test_chain->SetBranchAddress("GEN_charge", &GEN_charge_test_br);
    
    
    //===================================================
    //book histograms for target/input for train sample
    //===================================================
    
    TH2F *Deviation_dPhiSum3_train = new TH2F("Deviation_dPhiSum3_train","Deviation vs dPhiSum3 train",30,-2,2,100,-0.1,0.1);
    TH2F *Deviation_FR_1_train = new TH2F("Deviation_FR_1_train","Deviation vs FR_1 train",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_dTh_13_train = new TH2F("Deviation_dTh_13_train","Deviation vs dTh_13 train",100,-10,10,100,-0.1,0.1);
    TH2F *Deviation_bend_1_train = new TH2F("Deviation_bend_1_train","Deviation vs bend_1 train",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_theta_train = new TH2F("Deviation_theta_train","Deviation vs theta train",100,0,90,100,-0.1,0.1);
    TH2F *Deviation_dPhi_12_train = new TH2F("Deviation_dPhi_12_train","Deviation vs dPhi_12 train",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_23_train = new TH2F("Deviation_dPhi_23_train","Deviation vs dPhi_23 train",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_34_train = new TH2F("Deviation_dPhi_34_train","Deviation vs dPhi_34 train",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_GEN_pt_train = new TH2F("Deviation_GEN_pt_train","Deviation vs GEN_pt train",100,0,1000,100,-0.1,0.1);
    TH2F *EMTF_Deviation_GEN_pt_train = new TH2F("EMTF_Deviation_GEN_pt_train","EMTF Deviation vs GEN_pt train",100,0,1000,100,-0.1,0.1);
    
    //for GEN pT: 1-8 GeV
    TH2F *Deviation_dPhiSum3_train_1_8 = new TH2F("Deviation_dPhiSum3_train_1_8","Deviation vs dPhiSum3 train 1<pT<=8 GeV",30,-2,2,100,-0.1,0.1);
    TH2F *Deviation_FR_1_train_1_8 = new TH2F("Deviation_FR_1_train_1_8","Deviation vs FR_1 train 1<pT<=8 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_dTh_13_train_1_8 = new TH2F("Deviation_dTh_13_train_1_8","Deviation vs dTh_13 train 1<pT<=8 GeV",100,-10,10,100,-0.1,0.1);
    TH2F *Deviation_bend_1_train_1_8 = new TH2F("Deviation_bend_1_train_1_8","Deviation vs bend_1 train 1<pT<=8 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_theta_train_1_8 = new TH2F("Deviation_theta_train_1_8","Deviation vs theta train 1<pT<=8 GeV",100,0,90,100,-0.1,0.1);
    TH2F *Deviation_dPhi_12_train_1_8 = new TH2F("Deviation_dPhi_12_train_1_8","Deviation vs dPhi_12 train 1<pT<=8 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_23_train_1_8 = new TH2F("Deviation_dPhi_23_train_1_8","Deviation vs dPhi_23 train 1<pT<=8 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_34_train_1_8 = new TH2F("Deviation_dPhi_34_train_1_8","Deviation vs dPhi_34 train 1<pT<=8 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_GEN_pt_train_1_8 = new TH2F("Deviation_GEN_pt_train_1_8","Deviation vs GEN_pt train 1<pT<=8 GeV",100,0,10,100,-0.1,0.1);
    TH2F *EMTF_Deviation_GEN_pt_train_1_8 = new TH2F("EMTF_Deviation_GEN_pt_train_1_8","EMTF Deviation vs GEN_pt train 1<pT<=8 GeV",100,0,10,100,-0.1,0.1);
    
    //for GEN pT: 8-30 GeV
    TH2F *Deviation_dPhiSum3_train_8_30 = new TH2F("Deviation_dPhiSum3_train_8_30","Deviation vs dPhiSum3 train 8<pT<=30 GeV",30,-2,2,100,-0.1,0.1);
    TH2F *Deviation_FR_1_train_8_30 = new TH2F("Deviation_FR_1_train_8_30","Deviation vs FR_1 train 8<pT<=30 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_dTh_13_train_8_30 = new TH2F("Deviation_dTh_13_train_8_30","Deviation vs dTh_13 train 8<pT<=30 GeV",100,-10,10,100,-0.1,0.1);
    TH2F *Deviation_bend_1_train_8_30 = new TH2F("Deviation_bend_1_train_8_30","Deviation vs bend_1 train 8<pT<=30 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_theta_train_8_30 = new TH2F("Deviation_theta_train_8_30","Deviation vs theta train 8<pT<=30 GeV",100,0,90,100,-0.1,0.1);
    TH2F *Deviation_dPhi_12_train_8_30 = new TH2F("Deviation_dPhi_12_train_8_30","Deviation vs dPhi_12 train 8<pT<=30 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_23_train_8_30 = new TH2F("Deviation_dPhi_23_train_8_30","Deviation vs dPhi_23 train 8<pT<=30 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_34_train_8_30 = new TH2F("Deviation_dPhi_34_train_8_30","Deviation vs dPhi_34 train 8<pT<=30 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_GEN_pt_train_8_30 = new TH2F("Deviation_GEN_pt_train_8_30","Deviation vs GEN_pt train 8<pT<=30 GeV",100,8,30,100,-0.1,0.1);
    TH2F *EMTF_Deviation_GEN_pt_train_8_30 = new TH2F("EMTF_Deviation_GEN_pt_train_8_30","EMTF Deviation vs GEN_pt train 8<pT<=30 GeV",100,8,30,100,-0.1,0.1);
    
    //for GEN pT: 30-120 GeV
    TH2F *Deviation_dPhiSum3_train_30_120 = new TH2F("Deviation_dPhiSum3_train_30_120","Deviation vs dPhiSum3 train 30<pT<=120 GeV",30,-2,2,100,-0.1,0.1);
    TH2F *Deviation_FR_1_train_30_120 = new TH2F("Deviation_FR_1_train_30_120","Deviation vs FR_1 train 30<pT<=120 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_dTh_13_train_30_120 = new TH2F("Deviation_dTh_13_train_30_120","Deviation vs dTh_13 train 30<pT<=120 GeV",100,-10,10,100,-0.1,0.1);
    TH2F *Deviation_bend_1_train_30_120 = new TH2F("Deviation_bend_1_train_30_120","Deviation vs bend_1 train 30<pT<=120 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_theta_train_30_120 = new TH2F("Deviation_theta_train_30_120","Deviation vs theta train 30<pT<=120 GeV",100,0,90,100,-0.1,0.1);
    TH2F *Deviation_dPhi_12_train_30_120 = new TH2F("Deviation_dPhi_12_train_30_120","Deviation vs dPhi_12 train 30<pT<=120 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_23_train_30_120 = new TH2F("Deviation_dPhi_23_train_30_120","Deviation vs dPhi_23 train 30<pT<=120 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_34_train_30_120 = new TH2F("Deviation_dPhi_34_train_30_120","Deviation vs dPhi_34 train 30<pT<=120 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_GEN_pt_train_30_120 = new TH2F("Deviation_GEN_pt_train_30_120","Deviation vs GEN_pt train 30<pT<=120 GeV",100,30,120,100,-0.1,0.1);
    TH2F *EMTF_Deviation_GEN_pt_train_30_120 = new TH2F("EMTF_Deviation_GEN_pt_train_30_120","EMTF Deviation vs GEN_pt train 30<pT<=120 GeV",100,30,120,100,-0.1,0.1);
    
    //for GEN pT: 120-1000 GeV
    TH2F *Deviation_dPhiSum3_train_120_1000 = new TH2F("Deviation_dPhiSum3_train_120_1000","Deviation vs dPhiSum3 train 120<pT<=1000 GeV",30,-2,2,100,-0.1,0.1);
    TH2F *Deviation_FR_1_train_120_1000 = new TH2F("Deviation_FR_1_train_120_1000","Deviation vs FR_1 train 120<pT<=1000 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_dTh_13_train_120_1000 = new TH2F("Deviation_dTh_13_train_120_1000","Deviation vs dTh_13 train 120<pT<=1000 GeV",100,-10,10,100,-0.1,0.1);
    TH2F *Deviation_bend_1_train_120_1000 = new TH2F("Deviation_bend_1_train_120_1000","Deviation vs bend_1 train 120<pT<=1000 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_theta_train_120_1000 = new TH2F("Deviation_theta_train_120_1000","Deviation vs theta train 120<pT<=1000 GeV",100,0,90,100,-0.1,0.1);
    TH2F *Deviation_dPhi_12_train_120_1000 = new TH2F("Deviation_dPhi_12_train_120_1000","Deviation vs dPhi_12 train 120<pT<=1000 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_23_train_120_1000 = new TH2F("Deviation_dPhi_23_train_120_1000","Deviation vs dPhi_23 train 120<pT<=1000 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_34_train_120_1000 = new TH2F("Deviation_dPhi_34_train_120_1000","Deviation vs dPhi_34 train 120<pT<=1000 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_GEN_pt_train_120_1000 = new TH2F("Deviation_GEN_pt_train_120_1000","Deviation vs GEN_pt train 120<pT<=1000 GeV",100,120,1000,100,-0.1,0.1);
    TH2F *EMTF_Deviation_GEN_pt_train_120_1000 = new TH2F("EMTF_Deviation_GEN_pt_train_120_1000","EMTF Deviation vs GEN_pt train 120<pT<=1000 GeV",100,120,1000,100,-0.1,0.1);

    //==================================================
    //book histograms for target/input for test sample
    //==================================================
    
    TH2F *Deviation_dPhiSum3_test = new TH2F("Deviation_dPhiSum3_test","Deviation vs dPhiSum3 test",30,-2,2,100,-0.1,0.1);
    TH2F *Deviation_FR_1_test = new TH2F("Deviation_FR_1_test","Deviation vs FR_1 test",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_dTh_13_test = new TH2F("Deviation_dTh_13_test","Deviation vs dTh_13 test",100,-10,10,100,-0.1,0.1);
    TH2F *Deviation_bend_1_test = new TH2F("Deviation_bend_1_test","Deviation vs bend_1 test",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_theta_test = new TH2F("Deviation_theta_test","Deviation vs theta test",100,0,90,100,-0.1,0.1);
    TH2F *Deviation_dPhi_12_test = new TH2F("Deviation_dPhi_12_test","Deviation vs dPhi_12 test",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_23_test = new TH2F("Deviation_dPhi_23_test","Deviation vs dPhi_23 test",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_34_test = new TH2F("Deviation_dPhi_34_test","Deviation vs dPhi_34 test",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_GEN_pt_test = new TH2F("Deviation_GEN_pt_test","Deviation vs GEN_pt test",100,0,1000,100,-0.1,0.1);
    TH2F *EMTF_Deviation_GEN_pt_test = new TH2F("EMTF_Deviation_GEN_pt_test","EMTF Deviation vs GEN_pt test",100,0,1000,100,-0.1,0.1);
    
    //for GEN pT: 1-8 GeV
    TH2F *Deviation_dPhiSum3_test_1_8 = new TH2F("Deviation_dPhiSum3_test_1_8","Deviation vs dPhiSum3 test 1<pT<=8 GeV",30,-2,2,100,-0.1,0.1);
    TH2F *Deviation_FR_1_test_1_8 = new TH2F("Deviation_FR_1_test_1_8","Deviation vs FR_1 test 1<pT<=8 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_dTh_13_test_1_8 = new TH2F("Deviation_dTh_13_test_1_8","Deviation vs dTh_13 test 1<pT<=8 GeV",100,-10,10,100,-0.1,0.1);
    TH2F *Deviation_bend_1_test_1_8 = new TH2F("Deviation_bend_1_test_1_8","Deviation vs bend_1 test 1<pT<=8 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_theta_test_1_8 = new TH2F("Deviation_theta_test_1_8","Deviation vs theta test 1<pT<=8 GeV",100,0,90,100,-0.1,0.1);
    TH2F *Deviation_dPhi_12_test_1_8 = new TH2F("Deviation_dPhi_12_test_1_8","Deviation vs dPhi_12 test 1<pT<=8 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_23_test_1_8 = new TH2F("Deviation_dPhi_23_test_1_8","Deviation vs dPhi_23 test 1<pT<=8 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_34_test_1_8 = new TH2F("Deviation_dPhi_34_test_1_8","Deviation vs dPhi_34 test 1<pT<=8 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_GEN_pt_test_1_8 = new TH2F("Deviation_GEN_pt_test_1_8","Deviation vs GEN_pt test 1<pT<=8 GeV",100,0,10,100,-0.1,0.1);
    TH2F *EMTF_Deviation_GEN_pt_test_1_8 = new TH2F("EMTF_Deviation_GEN_pt_test_1_8","EMTF Deviation vs GEN_pt test 1<pT<=8 GeV",100,0,10,100,-0.1,0.1);
    
    //for GEN pT: 8-30 GeV
    TH2F *Deviation_dPhiSum3_test_8_30 = new TH2F("Deviation_dPhiSum3_test_8_30","Deviation vs dPhiSum3 test 8<pT<=30 GeV",30,-2,2,100,-0.1,0.1);
    TH2F *Deviation_FR_1_test_8_30 = new TH2F("Deviation_FR_1_test_8_30","Deviation vs FR_1 test 8<pT<=30 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_dTh_13_test_8_30 = new TH2F("Deviation_dTh_13_test_8_30","Deviation vs dTh_13 test 8<pT<=30 GeV",100,-10,10,100,-0.1,0.1);
    TH2F *Deviation_bend_1_test_8_30 = new TH2F("Deviation_bend_1_test_8_30","Deviation vs bend_1 test 8<pT<=30 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_theta_test_8_30 = new TH2F("Deviation_theta_test_8_30","Deviation vs theta test 8<pT<=30 GeV",100,0,90,100,-0.1,0.1);
    TH2F *Deviation_dPhi_12_test_8_30 = new TH2F("Deviation_dPhi_12_test_8_30","Deviation vs dPhi_12 test 8<pT<=30 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_23_test_8_30 = new TH2F("Deviation_dPhi_23_test_8_30","Deviation vs dPhi_23 test 8<pT<=30 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_34_test_8_30 = new TH2F("Deviation_dPhi_34_test_8_30","Deviation vs dPhi_34 test 8<pT<=30 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_GEN_pt_test_8_30 = new TH2F("Deviation_GEN_pt_test_8_30","Deviation vs GEN_pt test 8<pT<=30 GeV",100,8,30,100,-0.1,0.1);
    TH2F *EMTF_Deviation_GEN_pt_test_8_30 = new TH2F("EMTF_Deviation_GEN_pt_test_8_30","EMTF Deviation vs GEN_pt test 8<pT<=30 GeV",100,8,30,100,-0.1,0.1);
    
    //for GEN pT: 30-120 GeV
    TH2F *Deviation_dPhiSum3_test_30_120 = new TH2F("Deviation_dPhiSum3_test_30_120","Deviation vs dPhiSum3 test 30<pT<=120 GeV",30,-2,2,100,-0.1,0.1);
    TH2F *Deviation_FR_1_test_30_120 = new TH2F("Deviation_FR_1_test_30_120","Deviation vs FR_1 test 30<pT<=120 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_dTh_13_test_30_120 = new TH2F("Deviation_dTh_13_test_30_120","Deviation vs dTh_13 test 30<pT<=120 GeV",100,-10,10,100,-0.1,0.1);
    TH2F *Deviation_bend_1_test_30_120 = new TH2F("Deviation_bend_1_test_30_120","Deviation vs bend_1 test 30<pT<=120 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_theta_test_30_120 = new TH2F("Deviation_theta_test_30_120","Deviation vs theta test 30<pT<=120 GeV",100,0,90,100,-0.1,0.1);
    TH2F *Deviation_dPhi_12_test_30_120 = new TH2F("Deviation_dPhi_12_test_30_120","Deviation vs dPhi_12 test 30<pT<=120 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_23_test_30_120 = new TH2F("Deviation_dPhi_23_test_30_120","Deviation vs dPhi_23 test 30<pT<=120 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_34_test_30_120 = new TH2F("Deviation_dPhi_34_test_30_120","Deviation vs dPhi_34 test 30<pT<=120 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_GEN_pt_test_30_120 = new TH2F("Deviation_GEN_pt_test_30_120","Deviation vs GEN_pt test 30<pT<=120 GeV",100,30,120,100,-0.1,0.1);
    TH2F *EMTF_Deviation_GEN_pt_test_30_120 = new TH2F("EMTF_Deviation_GEN_pt_test_30_120","EMTF Deviation vs GEN_pt test 30<pT<=120 GeV",100,30,120,100,-0.1,0.1);
    
    //for GEN pT: 120-1000 GeV
    TH2F *Deviation_dPhiSum3_test_120_1000 = new TH2F("Deviation_dPhiSum3_test_120_1000","Deviation vs dPhiSum3 test 120<pT<=1000 GeV",30,-2,2,100,-0.1,0.1);
    TH2F *Deviation_FR_1_test_120_1000 = new TH2F("Deviation_FR_1_test_120_1000","Deviation vs FR_1 test 120<pT<=1000 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_dTh_13_test_120_1000 = new TH2F("Deviation_dTh_13_test_120_1000","Deviation vs dTh_13 test 120<pT<=1000 GeV",100,-10,10,100,-0.1,0.1);
    TH2F *Deviation_bend_1_test_120_1000 = new TH2F("Deviation_bend_1_test_120_1000","Deviation vs bend_1 test 120<pT<=1000 GeV",3,-1,2,100,-0.1,0.1);
    TH2F *Deviation_theta_test_120_1000 = new TH2F("Deviation_theta_test_120_1000","Deviation vs theta test 120<pT<=1000 GeV",100,0,90,100,-0.1,0.1);
    TH2F *Deviation_dPhi_12_test_120_1000 = new TH2F("Deviation_dPhi_12_test_120_1000","Deviation vs dPhi_12 test 120<pT<=1000 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_23_test_120_1000 = new TH2F("Deviation_dPhi_23_test_120_1000","Deviation vs dPhi_23 test 120<pT<=1000 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_dPhi_34_test_120_1000 = new TH2F("Deviation_dPhi_34_test_120_1000","Deviation vs dPhi_34 test 120<pT<=1000 GeV",100,-4,4,100,-0.1,0.1);
    TH2F *Deviation_GEN_pt_test_120_1000 = new TH2F("Deviation_GEN_pt_test_120_1000","Deviation vs GEN_pt test 120<pT<=1000 GeV",100,120,1000,100,-0.1,0.1);
    TH2F *EMTF_Deviation_GEN_pt_test_120_1000 = new TH2F("EMTF_Deviation_GEN_pt_test_120_1000","EMTF Deviation vs GEN_pt test 120<pT<=1000 GeV",100,120,1000,100,-0.1,0.1);
    
    train_events = train_chain->GetEntries();
    test_events = test_chain->GetEntries();
    
  //std::cout << "\n******* About to enter the train event loop *******" << std::endl;
  for (UInt_t iEvt = 0; iEvt < train_chain->GetEntries(); iEvt++) {

  //  if (iEvt > 10) break;
  //  if ( (iEvt % 1) == 0 )
    //  std::cout << "\n*** Looking at event " << iEvt << " ***" <<  std::endl;

    train_chain->GetEntry(iEvt);

    //std::cout << "EMTF_pt = " << EMTF_pt_train_br << std::endl;
    //std::cout << "GEN_pt = " << GEN_pt_train_br << std::endl;
    //std::cout << "inv_GEN_pt_train = " << inv_GEN_pt_train_br << std::endl;
    //std::cout << "KNN_train = " << KNN_train_br << std::endl;
    //std::cout << "dPhi_12 = " << dPhi_12_train_br << std::endl;
      
      //Deviation_dPhiSum3_train->Fill(dPhiSum3_train_br,KNN_train_br-inv_GEN_pt_train_br);
      Deviation_FR_1_train->Fill(FR_1_train_br,KNN_train_br-inv_GEN_pt_train_br);
      //Deviation_dTh_13_train->Fill(dTh_13_train_br,KNN_train_br-inv_GEN_pt_train_br);
      Deviation_bend_1_train->Fill(bend_1_train_br,KNN_train_br-inv_GEN_pt_train_br);
      Deviation_theta_train->Fill(theta_train_br,KNN_train_br-inv_GEN_pt_train_br);
      Deviation_dPhi_12_train->Fill(dPhi_12_train_br,KNN_train_br-inv_GEN_pt_train_br);
      Deviation_dPhi_23_train->Fill(dPhi_23_train_br,KNN_train_br-inv_GEN_pt_train_br);
      //Deviation_dPhi_34_train->Fill(dPhi_34_train_br,KNN_train_br-inv_GEN_pt_train_br);
      Deviation_GEN_pt_train->Fill(GEN_pt_train_br,KNN_train_br-inv_GEN_pt_train_br);
      EMTF_Deviation_GEN_pt_train->Fill(GEN_pt_train_br,inv_EMTF_pt_train_br-inv_GEN_pt_train_br);
      
      Deviation_train += (KNN_train_br-inv_GEN_pt_train_br)*(KNN_train_br-inv_GEN_pt_train_br);
      EMTF_Deviation_train += (inv_EMTF_pt_train_br-inv_GEN_pt_train_br)*(inv_EMTF_pt_train_br-inv_GEN_pt_train_br);
      
      //devide by pT:1-8GeV
      if (GEN_pt_train_br > 1 && GEN_pt_train_br <= 8){
          
          //Deviation_dPhiSum3_train_1_8->Fill(dPhiSum3_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_FR_1_train_1_8->Fill(FR_1_train_br,KNN_train_br-inv_GEN_pt_train_br);
          //Deviation_dTh_13_train_1_8->Fill(dTh_13_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_bend_1_train_1_8->Fill(bend_1_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_theta_train_1_8->Fill(theta_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_dPhi_12_train_1_8->Fill(dPhi_12_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_dPhi_23_train_1_8->Fill(dPhi_23_train_br,KNN_train_br-inv_GEN_pt_train_br);
          //Deviation_dPhi_34_train_1_8->Fill(dPhi_34_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_GEN_pt_train_1_8->Fill(GEN_pt_train_br,KNN_train_br-inv_GEN_pt_train_br);
          EMTF_Deviation_GEN_pt_train_1_8->Fill(GEN_pt_train_br,inv_EMTF_pt_train_br-inv_GEN_pt_train_br);
          
          train_events_1_8 += 1;
          Deviation_train_1_8 += (KNN_train_br-inv_GEN_pt_train_br)*(KNN_train_br-inv_GEN_pt_train_br);
          EMTF_Deviation_train_1_8 += (inv_EMTF_pt_train_br-inv_GEN_pt_train_br)*(inv_EMTF_pt_train_br-inv_GEN_pt_train_br);
      }
      else if (GEN_pt_train_br > 8 && GEN_pt_train_br <= 30){
          
          //Deviation_dPhiSum3_train_8_30->Fill(dPhiSum3_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_FR_1_train_8_30->Fill(FR_1_train_br,KNN_train_br-inv_GEN_pt_train_br);
          //Deviation_dTh_13_train_8_30->Fill(dTh_13_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_bend_1_train_8_30->Fill(bend_1_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_theta_train_8_30->Fill(theta_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_dPhi_12_train_8_30->Fill(dPhi_12_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_dPhi_23_train_8_30->Fill(dPhi_23_train_br,KNN_train_br-inv_GEN_pt_train_br);
          //Deviation_dPhi_34_train_8_30->Fill(dPhi_34_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_GEN_pt_train_8_30->Fill(GEN_pt_train_br,KNN_train_br-inv_GEN_pt_train_br);
          EMTF_Deviation_GEN_pt_train_8_30->Fill(GEN_pt_train_br,inv_EMTF_pt_train_br-inv_GEN_pt_train_br);
          
          train_events_8_30 += 1;
          Deviation_train_8_30 += (KNN_train_br-inv_GEN_pt_train_br)*(KNN_train_br-inv_GEN_pt_train_br);
          EMTF_Deviation_train_8_30 += (inv_EMTF_pt_train_br-inv_GEN_pt_train_br)*(inv_EMTF_pt_train_br-inv_GEN_pt_train_br);
      }
      else if (GEN_pt_train_br > 30 && GEN_pt_train_br <= 120){
          
          //Deviation_dPhiSum3_train_30_120->Fill(dPhiSum3_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_FR_1_train_30_120->Fill(FR_1_train_br,KNN_train_br-inv_GEN_pt_train_br);
          //Deviation_dTh_13_train_30_120->Fill(dTh_13_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_bend_1_train_30_120->Fill(bend_1_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_theta_train_30_120->Fill(theta_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_dPhi_12_train_30_120->Fill(dPhi_12_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_dPhi_23_train_30_120->Fill(dPhi_23_train_br,KNN_train_br-inv_GEN_pt_train_br);
          //Deviation_dPhi_34_train_30_120->Fill(dPhi_34_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_GEN_pt_train_30_120->Fill(GEN_pt_train_br,KNN_train_br-inv_GEN_pt_train_br);
          EMTF_Deviation_GEN_pt_train_30_120->Fill(GEN_pt_train_br,inv_EMTF_pt_train_br-inv_GEN_pt_train_br);
          
          train_events_30_120 += 1;
          Deviation_train_30_120 += (KNN_train_br-inv_GEN_pt_train_br)*(KNN_train_br-inv_GEN_pt_train_br);
          EMTF_Deviation_train_30_120 += (inv_EMTF_pt_train_br-inv_GEN_pt_train_br)*(inv_EMTF_pt_train_br-inv_GEN_pt_train_br);
      }
      else if (GEN_pt_train_br > 120 && GEN_pt_train_br <= 1000){
          
          //Deviation_dPhiSum3_train_120_1000->Fill(dPhiSum3_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_FR_1_train_120_1000->Fill(FR_1_train_br,KNN_train_br-inv_GEN_pt_train_br);
          //Deviation_dTh_13_train_120_1000->Fill(dTh_13_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_bend_1_train_120_1000->Fill(bend_1_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_theta_train_120_1000->Fill(theta_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_dPhi_12_train_120_1000->Fill(dPhi_12_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_dPhi_23_train_120_1000->Fill(dPhi_23_train_br,KNN_train_br-inv_GEN_pt_train_br);
          //Deviation_dPhi_34_train_120_1000->Fill(dPhi_34_train_br,KNN_train_br-inv_GEN_pt_train_br);
          Deviation_GEN_pt_train_120_1000->Fill(GEN_pt_train_br,KNN_train_br-inv_GEN_pt_train_br);
          EMTF_Deviation_GEN_pt_train_120_1000->Fill(GEN_pt_train_br,inv_EMTF_pt_train_br-inv_GEN_pt_train_br);
          
          train_events_120_1000 += 1;
          Deviation_train_120_1000 += (KNN_train_br-inv_GEN_pt_train_br)*(KNN_train_br-inv_GEN_pt_train_br);
          EMTF_Deviation_train_120_1000 += (inv_EMTF_pt_train_br-inv_GEN_pt_train_br)*(inv_EMTF_pt_train_br-inv_GEN_pt_train_br);
      }
      
  } // End loop: for (UInt_t iEvt = 0; iEvt < train_chain->GetEntries(); iEvt++)
    //std::cout << "\n******* Leaving the train event loop *******" << std::endl;
    
  std::cout << "\n******* 1st Enter the test event loop *******" << std::endl;
  for (UInt_t iEvt = 0; iEvt < test_chain->GetEntries(); iEvt++) {

  //  if (iEvt > 10) break;
  //  if ( (iEvt % 1) == 0 ) 
      //std::cout << "\n*** Looking at event " << iEvt << " ***" <<  std::endl;

    test_chain->GetEntry(iEvt);

    //std::cout << "EMTF_pt = " << EMTF_pt_test_br << std::endl;
    //std::cout << "GEN_pt = " << GEN_pt_test_br << std::endl;
    //std::cout << "inv_GEN_pt_test = " << inv_GEN_pt_test_br << std::endl;
    //std::cout << "KNN_test = " << KNN_test_br << std::endl;
    //std::cout << "dPhi_12 = " << dPhi_12_test_br << std::endl;
     
      //Deviation_dPhiSum3_test->Fill(dPhiSum3_test_br,KNN_test_br-inv_GEN_pt_test_br);
      Deviation_FR_1_test->Fill(FR_1_test_br,KNN_test_br-inv_GEN_pt_test_br);
      //Deviation_dTh_13_test->Fill(dTh_13_test_br,KNN_test_br-inv_GEN_pt_test_br);
      Deviation_bend_1_test->Fill(bend_1_test_br,KNN_test_br-inv_GEN_pt_test_br);
      Deviation_theta_test->Fill(theta_test_br,KNN_test_br-inv_GEN_pt_test_br);
      Deviation_dPhi_12_test->Fill(dPhi_12_test_br,KNN_test_br-inv_GEN_pt_test_br);
      Deviation_dPhi_23_test->Fill(dPhi_23_test_br,KNN_test_br-inv_GEN_pt_test_br);
      //Deviation_dPhi_34_test->Fill(dPhi_34_test_br,KNN_test_br-inv_GEN_pt_test_br);
      Deviation_GEN_pt_test->Fill(GEN_pt_test_br,KNN_test_br-inv_GEN_pt_test_br);
      EMTF_Deviation_GEN_pt_test->Fill(GEN_pt_test_br,inv_EMTF_pt_test_br-inv_GEN_pt_test_br);
      
      Deviation_test += (KNN_test_br-inv_GEN_pt_test_br)*(KNN_test_br-inv_GEN_pt_test_br);
      EMTF_Deviation_test += (inv_EMTF_pt_test_br-inv_GEN_pt_test_br)*(inv_EMTF_pt_test_br-inv_GEN_pt_test_br);
      
      //devide by pT:1-8GeV
      if (GEN_pt_test_br > 1 && GEN_pt_test_br <= 8){
          
          //Deviation_dPhiSum3_test_1_8->Fill(dPhiSum3_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_FR_1_test_1_8->Fill(FR_1_test_br,KNN_test_br-inv_GEN_pt_test_br);
          //Deviation_dTh_13_test_1_8->Fill(dTh_13_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_bend_1_test_1_8->Fill(bend_1_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_theta_test_1_8->Fill(theta_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_dPhi_12_test_1_8->Fill(dPhi_12_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_dPhi_23_test_1_8->Fill(dPhi_23_test_br,KNN_test_br-inv_GEN_pt_test_br);
          //Deviation_dPhi_34_test_1_8->Fill(dPhi_34_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_GEN_pt_test_1_8->Fill(GEN_pt_test_br,KNN_test_br-inv_GEN_pt_test_br);
          EMTF_Deviation_GEN_pt_test_1_8->Fill(GEN_pt_test_br,inv_EMTF_pt_test_br-inv_GEN_pt_test_br);
          
          test_events_1_8 += 1;
          Deviation_test_1_8 += (KNN_test_br-inv_GEN_pt_test_br)*(KNN_test_br-inv_GEN_pt_test_br);
          EMTF_Deviation_test_1_8 += (inv_EMTF_pt_test_br-inv_GEN_pt_test_br)*(inv_EMTF_pt_test_br-inv_GEN_pt_test_br);
      }
      else if (GEN_pt_test_br > 8 && GEN_pt_test_br <= 30){
          //Deviation_dPhiSum3_test_8_30->Fill(dPhiSum3_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_FR_1_test_8_30->Fill(FR_1_test_br,KNN_test_br-inv_GEN_pt_test_br);
          //Deviation_dTh_13_test_8_30->Fill(dTh_13_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_bend_1_test_8_30->Fill(bend_1_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_theta_test_8_30->Fill(theta_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_dPhi_12_test_8_30->Fill(dPhi_12_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_dPhi_23_test_8_30->Fill(dPhi_23_test_br,KNN_test_br-inv_GEN_pt_test_br);
          //Deviation_dPhi_34_test_8_30->Fill(dPhi_34_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_GEN_pt_test_8_30->Fill(GEN_pt_test_br,KNN_test_br-inv_GEN_pt_test_br);
          EMTF_Deviation_GEN_pt_test_8_30->Fill(GEN_pt_test_br,inv_EMTF_pt_test_br-inv_GEN_pt_test_br);
          
          test_events_8_30 += 1;
          Deviation_test_8_30 += (KNN_test_br-inv_GEN_pt_test_br)*(KNN_test_br-inv_GEN_pt_test_br);
          EMTF_Deviation_test_8_30 += (inv_EMTF_pt_test_br-inv_GEN_pt_test_br)*(inv_EMTF_pt_test_br-inv_GEN_pt_test_br);
      }
      else if (GEN_pt_test_br > 30 && GEN_pt_test_br <= 120){
          //Deviation_dPhiSum3_test_30_120->Fill(dPhiSum3_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_FR_1_test_30_120->Fill(FR_1_test_br,KNN_test_br-inv_GEN_pt_test_br);
          //Deviation_dTh_13_test_30_120->Fill(dTh_13_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_bend_1_test_30_120->Fill(bend_1_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_theta_test_30_120->Fill(theta_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_dPhi_12_test_30_120->Fill(dPhi_12_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_dPhi_23_test_30_120->Fill(dPhi_23_test_br,KNN_test_br-inv_GEN_pt_test_br);
          //Deviation_dPhi_34_test_30_120->Fill(dPhi_34_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_GEN_pt_test_30_120->Fill(GEN_pt_test_br,KNN_test_br-inv_GEN_pt_test_br);
          EMTF_Deviation_GEN_pt_test_30_120->Fill(GEN_pt_test_br,inv_EMTF_pt_test_br-inv_GEN_pt_test_br);
          
          test_events_30_120 += 1;
          Deviation_test_30_120 += (KNN_test_br-inv_GEN_pt_test_br)*(KNN_test_br-inv_GEN_pt_test_br);
          EMTF_Deviation_test_30_120 += (inv_EMTF_pt_test_br-inv_GEN_pt_test_br)*(inv_EMTF_pt_test_br-inv_GEN_pt_test_br);
      }
      else if (GEN_pt_test_br > 120 && GEN_pt_test_br <= 1000){
          //Deviation_dPhiSum3_test_120_1000->Fill(dPhiSum3_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_FR_1_test_120_1000->Fill(FR_1_test_br,KNN_test_br-inv_GEN_pt_test_br);
          //Deviation_dTh_13_test_120_1000->Fill(dTh_13_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_bend_1_test_120_1000->Fill(bend_1_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_theta_test_120_1000->Fill(theta_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_dPhi_12_test_120_1000->Fill(dPhi_12_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_dPhi_23_test_120_1000->Fill(dPhi_23_test_br,KNN_test_br-inv_GEN_pt_test_br);
          //Deviation_dPhi_34_test_120_1000->Fill(dPhi_34_test_br,KNN_test_br-inv_GEN_pt_test_br);
          Deviation_GEN_pt_test_120_1000->Fill(GEN_pt_test_br,KNN_test_br-inv_GEN_pt_test_br);
          EMTF_Deviation_GEN_pt_test_120_1000->Fill(GEN_pt_test_br,inv_EMTF_pt_test_br-inv_GEN_pt_test_br);
          
          test_events_120_1000 += 1;
          Deviation_test_120_1000 += (KNN_test_br-inv_GEN_pt_test_br)*(KNN_test_br-inv_GEN_pt_test_br);
          EMTF_Deviation_test_120_1000 += (inv_EMTF_pt_test_br-inv_GEN_pt_test_br)*(inv_EMTF_pt_test_br-inv_GEN_pt_test_br);
      }
      
      //test trigger efficiency 59 bins, 59 cuts as well
      for (int Gen_bin=0;Gen_bin<59;Gen_bin++) {
          if (GEN_pt_test_br > GEN_pT[Gen_bin] && GEN_pt_test_br <= GEN_pT[Gen_bin+1]){
              //need to count muons in all bins, don't use break
              count[Gen_bin]++;
              
              //deal with each cut
              for(int Cut_bin=0;Cut_bin<59;Cut_bin++){
                  //KNN is targeting 1/pT, need to convert
                  if (1./KNN_test_br > trigger_Cut[Cut_bin]) {
                      pass_count[Cut_bin][Gen_bin]++;
                  }//end if KNN
                  if (EMTF_pt_test_br > trigger_Cut[Cut_bin]){
                      EMTF_pass_count[Cut_bin][Gen_bin]++;
                  }//end if EMTF
              }//end trigger cuts
              
              //deal with special selected cuts to compare efficiency
              for(int Cut_special_bin=0;Cut_special_bin<6;Cut_special_bin++){
                  //KNN is targeting 1/pT, need to convert
                  if (1./KNN_test_br > trigger_Cut_special[Cut_special_bin]) {
                      pass_count_special[Cut_special_bin][Gen_bin]++;
                  }//end if KNN
                  if (EMTF_pt_test_br > trigger_Cut_special[Cut_special_bin]){
                      EMTF_pass_count_special[Cut_special_bin][Gen_bin]++;
                  }//end if EMTF
             }//end for special cuts
              
          }//end if GEN_pt
      }//end for Gen_bin
    
  } // End loop: for (UInt_t iEvt = 0; iEvt < test_chain->GetEntries(); iEvt++)
  //std::cout << "\n******* Leaving the test event loop *******" << std::endl;
    
    //KNN
    Standard_Deviation_test = Deviation_test/test_events;
    Standard_Deviation_test_1_8 = Deviation_test_1_8/test_events_1_8;
    Standard_Deviation_test_8_30 = Deviation_test_8_30/test_events_8_30;
    Standard_Deviation_test_30_120 = Deviation_test_30_120/test_events_30_120;
    Standard_Deviation_test_120_1000 = Deviation_test_120_1000/test_events_120_1000;
    
    Standard_Deviation_train = Deviation_train/train_events;
    Standard_Deviation_train_1_8 = Deviation_train_1_8/train_events_1_8;
    Standard_Deviation_train_8_30 = Deviation_train_8_30/train_events_8_30;
    Standard_Deviation_train_30_120 = Deviation_train_30_120/train_events_30_120;
    Standard_Deviation_train_120_1000 = Deviation_train_120_1000/train_events_120_1000;
    //EMTF
    EMTF_Standard_Deviation_test = EMTF_Deviation_test/test_events;
    EMTF_Standard_Deviation_test_1_8 = EMTF_Deviation_test_1_8/test_events_1_8;
    EMTF_Standard_Deviation_test_8_30 = EMTF_Deviation_test_8_30/test_events_8_30;
    EMTF_Standard_Deviation_test_30_120 = EMTF_Deviation_test_30_120/test_events_30_120;
    EMTF_Standard_Deviation_test_120_1000 = EMTF_Deviation_test_120_1000/test_events_120_1000;
    
    EMTF_Standard_Deviation_train = EMTF_Deviation_train/train_events;
    EMTF_Standard_Deviation_train_1_8 = EMTF_Deviation_train_1_8/train_events_1_8;
    EMTF_Standard_Deviation_train_8_30 = EMTF_Deviation_train_8_30/train_events_8_30;
    EMTF_Standard_Deviation_train_30_120 = EMTF_Deviation_train_30_120/train_events_30_120;
    EMTF_Standard_Deviation_train_120_1000 = EMTF_Deviation_train_120_1000/train_events_120_1000;
    
    //we're interested in rate of low pT: 1-50 GeV, warning may appear because overflow after 50
    KNN_trigger_Gen_efficiency = new TProfile2D("KNN_trigger_Gen_efficiency","KNN trigger efficiency versus thresholds and GEN pT",49,1,50,49,1,50,0,1);
    EMTF_trigger_Gen_efficiency = new TProfile2D("EMTF_trigger_Gen_efficiency","EMTF trigger efficiency versus thresholds and GEN pT",49,1,50,49,1,50,0,1);
    
    //======================================================
    //calculate trigger efficiency KNN not scaled to 90% yet
    //======================================================
    for (int Gen_bin=0;Gen_bin<59;Gen_bin++) {
        
        for(int Cut_bin=0;Cut_bin<59;Cut_bin++){
            
            efficiency[Cut_bin][Gen_bin]=pass_count[Cut_bin][Gen_bin]/count[Gen_bin];
            EMTF_efficiency[Cut_bin][Gen_bin]=EMTF_pass_count[Cut_bin][Gen_bin]/count[Gen_bin];
            //fill 2D profile for non-scaled KNN
            KNN_trigger_Gen_efficiency->Fill(GEN_pT[Gen_bin],trigger_Cut[Cut_bin],efficiency[Cut_bin][Gen_bin]);
            EMTF_trigger_Gen_efficiency->Fill(GEN_pT[Gen_bin],trigger_Cut[Cut_bin],EMTF_efficiency[Cut_bin][Gen_bin]);
        }//end Cut bin
        
        for(int Cut_special_bin=0;Cut_special_bin<6;Cut_special_bin++){
            efficiency_special[Cut_special_bin][Gen_bin]=pass_count_special[Cut_special_bin][Gen_bin]/count[Gen_bin];
            EMTF_efficiency_special[Cut_special_bin][Gen_bin]=EMTF_pass_count_special[Cut_special_bin][Gen_bin]/count[Gen_bin];
        }//end Cut bin special
        
    }//end Gen bin
    
    //write to output file
    TFile myPlot("/afs/cern.ch/OUTPUTFILEDIRECTORY","RECREATE");
    
    //write 2D non scaled efficiency plot
    KNN_trigger_Gen_efficiency->Write();
    EMTF_trigger_Gen_efficiency->Write();
    
    //book graph for 6 special trigger efficiency,EMTF not scaled
    //initialization
    const Int_t special_cuts=6;
    TGraph *KNN_eff[special_cuts];
    TGraph *EMTF_eff[special_cuts];
    TCanvas *C[special_cuts];
    TMultiGraph *mg[special_cuts];
    for(Int_t i=0; i<special_cuts;i++){
        KNN_eff[i] = new TGraph(59,trigger_Cut,efficiency_special[i]); KNN_eff[i]->SetMarkerStyle(21); KNN_eff[i]->SetMarkerColor(2);//red
        EMTF_eff[i] = new TGraph(59,trigger_Cut,EMTF_efficiency_special[i]); EMTF_eff[i]->SetMarkerStyle(21); EMTF_eff[i]->SetMarkerColor(1);//black
        C[i] = new TCanvas(Form("C%d",i),Form("Efficiency_%d",i),700,500);
        mg[i] = new TMultiGraph();
        C[i]->cd();
        mg[i]->SetTitle(Form("Mode 14 trigger efficiency pT > %f GeV",trigger_Cut_special[i]));
        mg[i]->Add(KNN_eff[i]);
        mg[i]->Add(EMTF_eff[i]);
        mg[i]->Draw();
        mg[i]->Write();
    }
    
    //================================================================================
    //Scale KNN(and EMTF) to 90%:
    //loop over true pT, for KNN, needs to rescale threshold X by multiply a factor[0,1]
    //to achieve 90% eff when true pT = X, these factors vary depend
    //on what the X is. Under each X, loop over change the scale factor; Under each
    //scale factor, loop over all test events to figure out how many pT > X*scale factor,
    //divide by total count in [bin X] to calculate the efficiency. If not close
    //to 90% eff cut, continue incresing KNN assigned pT scale factor until it closes.
    //================================================================================
    //scale factor of threshold of KNN assigned pT to achieve 90% efficiency,
    //it's between 0 and 1. It will be used in rate plot
    
    //scale factor of current EMTF assigned pT to achieve 90% efficiency at trigger cuts,
    //it's b/t[0,1]. It can be used in comparison with KNN.
    double EMTF_scale_Max=1.0;//Current EMTF is mostly sacaled to 85%
    double EMTF_scale_Min=0.0;
    double EMTF_scale[59]={0};
    double EMTF_pass_count_scaled[59][59]={0};
    double EMTF_efficiency_scaled[59][59]={0};
    double EMTF_efficiency_scale=0.90;
    double EMTF_scale_consistency[59]={0};
    EMTF_scale_plot = new TProfile("EMTF_scale_plot","EMTF scale versus thresholds",49,1,50,0,1);
    EMTF_scale_consistency_plot = new TProfile("EMTF_scale_consistency_plot","EMTF scale factor to 90% at thresholds",49,1,50,0,1);
    
    double KNN_scale_Max=1.0;
    double KNN_scale_Min=0.0;
    double KNN_scale[59]={0};
    //after KNN scaled
    double pass_count_scaled[59][59]={0};
    double efficiency_scaled[59][59]={0};
    double efficiency_scale=0.90;
    double KNN_scale_consistency[59]={0};
    KNN_scale_plot = new TProfile("KNN_scale_plot","KNN scale versus thresholds",49,1,50,0,1);
    KNN_scale_consistency_plot = new TProfile("KNN_scale_consistency_plot","KNN scale factor to 90% at thresholds",49,1,50,0,1);
    
    for(int Cut_bin=0;Cut_bin<59;Cut_bin++){
        
        double efficiency_tmp=0;
        double KNN_min=1.1;
        double KNN_scale_tmp=0;
        int flag=0;//stop loop over the scale for this bin when flag=1;
        
        double EMTF_efficiency_tmp=0;
        double EMTF_min=1.1;
        double EMTF_scale_tmp=0;
        int EMTF_flag=0;//stop loop over the scale for this bin when flag=1;
        
        for (int i=0; i<101; i++) {
            
            //initialize # KNN pT> X * this scale factor
            double pass_count_tmp=0;
            double EMTF_pass_count_tmp=0;
            if (flag == 1 && EMTF_flag == 1) break;
            
            //scale KNN assigned pT up to compare rate
            KNN_scale_tmp = KNN_scale_Max - i*(KNN_scale_Max-KNN_scale_Min)/100.0;
            EMTF_scale_tmp = EMTF_scale_Max - i*(EMTF_scale_Max-EMTF_scale_Min)/100.0;
            
            //==================================================
            //Enter test again for scale determine
            //===================================================
            //std::cout << "\n******* Enter the test event loop for Cut bin: "<<Cut_bin<< "; Current scale: "<<KNN_scale_tmp<<" *******" << std::endl;
            for (UInt_t iEvt = 0; iEvt < test_chain->GetEntries(); iEvt++) {
            
            //  if ( (iEvt % 1) == 0 )
            //std::cout << "\n*** Looking at event " << iEvt << " ***" <<  std::endl;
                test_chain->GetEntry(iEvt);
                
                //need to specify events in this bin
                if (GEN_pt_test_br > GEN_pT[Cut_bin] && GEN_pt_test_br <= GEN_pT[Cut_bin+1]){
            //KNN is targeting 1/pT, need to convert
                    if (1./KNN_test_br > trigger_Cut[Cut_bin]*KNN_scale_tmp) {
                        pass_count_tmp++;
                    }//end if KNN
                    if (EMTF_pt_test_br > trigger_Cut[Cut_bin]*EMTF_scale_tmp){
                        EMTF_pass_count_tmp++;
                    }//end if EMTF
                }//end if
                
            }//End loop: for (UInt_t iEvt = 0; iEvt < test_chain->GetEntries(); iEvt++)
            //****Leaving the test event loop ****
            
            //calculate eff under this scale factor
            efficiency_tmp = pass_count_tmp/count[Cut_bin];
            EMTF_efficiency_tmp = EMTF_pass_count_tmp/count[Cut_bin];
            
            if (fabs(efficiency_tmp-efficiency_scale) < KNN_min) {
                KNN_min = fabs(efficiency_tmp-efficiency_scale);
                KNN_scale[Cut_bin] = KNN_scale_tmp;
                KNN_scale_consistency[Cut_bin] = efficiency_tmp;
                //if goes past the scale, stop loop over efficiency scale
                if(efficiency_tmp >= efficiency_scale){
                    flag=1;
                }//end if flag
            }//end if
            
            if (fabs(EMTF_efficiency_tmp-EMTF_efficiency_scale) < EMTF_min) {
                EMTF_min = fabs(EMTF_efficiency_tmp-EMTF_efficiency_scale);
                EMTF_scale[Cut_bin] = EMTF_scale_tmp;
                EMTF_scale_consistency[Cut_bin] = EMTF_efficiency_tmp;
                //if goes past the scale, stop loop over efficiency scale
                if(EMTF_efficiency_tmp >= EMTF_efficiency_scale){
                    EMTF_flag=1;
                }//end if EMTF flag
            }//end if
            
        }//end varying scale fac
        std::cout << "\n******* Leave the test event loop for Cut bin: "<<Cut_bin<< "; KNN scale: "<<KNN_scale[Cut_bin]<<"; EMTF scale: "<<EMTF_scale[Cut_bin]<<" *******" << std::endl;
        
        //see if all 90% at all trigger thresholds
        KNN_scale_consistency_plot->Fill(trigger_Cut[Cut_bin],KNN_scale_consistency[Cut_bin]);
        KNN_scale_plot->Fill(trigger_Cut[Cut_bin],KNN_scale[Cut_bin]);
        trigger_Cut_scaled[Cut_bin] = trigger_Cut[Cut_bin]*KNN_scale[Cut_bin];
        
        EMTF_scale_consistency_plot->Fill(trigger_Cut[Cut_bin],EMTF_scale_consistency[Cut_bin]);
        EMTF_scale_plot->Fill(trigger_Cut[Cut_bin],EMTF_scale[Cut_bin]);
        EMTF_trigger_Cut_scaled[Cut_bin] = trigger_Cut[Cut_bin]*EMTF_scale[Cut_bin];
        
    }//end loop over Cuts
    
    KNN_scale_plot->Write();
    KNN_scale_consistency_plot->Write();
    
    EMTF_scale_plot->Write();
    EMTF_scale_consistency_plot->Write();
    
    //====================================================
    //recount using the new KNN scaled threshold
    //====================================================
    std::cout << "\n******* Enter the test event loop to recount passing trigger events *******" << std::endl;
    for (UInt_t iEvt = 0; iEvt < test_chain->GetEntries(); iEvt++) {
        
        //  if (iEvt > 10) break;
        //  if ( (iEvt % 1) == 0 )
        //std::cout << "\n*** Looking at event " << iEvt << " ***" <<  std::endl;
        
        test_chain->GetEntry(iEvt);
        
        for (int Gen_bin=0;Gen_bin<59;Gen_bin++) {
            
            if (GEN_pt_test_br > GEN_pT[Gen_bin] && GEN_pt_test_br <= GEN_pT[Gen_bin+1]){
            //deal with each cut
            for(int Cut_bin=0;Cut_bin<59;Cut_bin++){
                
                if (1./KNN_test_br > trigger_Cut_scaled[Cut_bin]) {//this scaled cut correspond to original trigger cut when plot 2d
                    pass_count_scaled[Cut_bin][Gen_bin]++;
                }//end if KNN
                if (EMTF_pt_test_br > EMTF_trigger_Cut_scaled[Cut_bin]){
                    EMTF_pass_count_scaled[Cut_bin][Gen_bin]++;
                }//end if EMTF
                
            }//end trigger cuts
            
        }//end if GEN_pt
        }//end Gen bin
        
    }// End loop: for (UInt_t iEvt = 0; iEvt < test_chain->GetEntries(); iEvt++)
    std::cout << "\n******* Finish recount, leaving the test event loop *******" << std::endl;
    
    //=========================================
    //calculate trigger efficiency after rescale
    //=========================================
    KNN_trigger_Gen_efficiency_scaled = new TProfile2D("KNN_trigger_Gen_efficiency_scaled","KNN trigger efficiency versus thresholds and GEN pT SCALED",49,1,50,49,1,50,0,1);
    EMTF_trigger_Gen_efficiency_scaled = new TProfile2D("EMTF_trigger_Gen_efficiency_scaled","EMTF trigger efficiency versus thresholds and GEN pT SCALED",49,1,50,49,1,50,0,1);
    
    //====================================================
    //check efficiency consistent for rate plot
    //====================================================
    double efficiency_threshold=0.90;
    double KNN_cuts[59]={0};
    double EMTF_cuts[59]={0};
    double KNN_efficiency_cuts[59]={0};
    double EMTF_efficiency_cuts[59]={0};
    KNN_efficiency_cuts_consistency = new TProfile("KNN_efficiency_cuts_consistency","KNN cut efficiency versus GEN pT",49,1,50,0,1);
    EMTF_efficiency_cuts_consistency = new TProfile("EMTF_efficiency_cuts_consistency","EMTF cut efficiency versus GEN pT",49,1,50,0,1);
    KNN_cuts_vs_Gen_pT = new TProfile("KNN_cuts_vs_Gen_pT","KNN cuts versus GEN pT",49,1,50,0,50);
    EMTF_cuts_vs_Gen_pT = new TProfile("EMTF_cuts_vs_Gen_pT","EMTF cuts versus GEN pT",49,1,50,0,50);
    
    
    for (int Gen_bin=0;Gen_bin<59;Gen_bin++) {
        
        double KNN_min=1.1;//calculate the difference of 0.95 and each bin, find the closet to 0.95
        double EMTF_min=1.1;
        
        for(int Cut_bin=0;Cut_bin<59;Cut_bin++){
            
            //KNN eff after scaled
            efficiency_scaled[Cut_bin][Gen_bin]=pass_count_scaled[Cut_bin][Gen_bin]/count[Gen_bin];
            EMTF_efficiency_scaled[Cut_bin][Gen_bin]=EMTF_pass_count_scaled[Cut_bin][Gen_bin]/count[Gen_bin];
            
            //===========================================================================
            //fill 2D profile, don't fill trigger cut scaled! fill original trigger cut!
            //===========================================================================
            KNN_trigger_Gen_efficiency_scaled->Fill(GEN_pT[Gen_bin],trigger_Cut[Cut_bin],efficiency_scaled[Cut_bin][Gen_bin]);
            EMTF_trigger_Gen_efficiency_scaled->Fill(GEN_pT[Gen_bin],trigger_Cut[Cut_bin],EMTF_efficiency_scaled[Cut_bin][Gen_bin]);
            
            if (fabs(efficiency_scaled[Cut_bin][Gen_bin]-efficiency_threshold)<KNN_min) {
                KNN_min = fabs(efficiency_scaled[Cut_bin][Gen_bin]-efficiency_threshold);
                KNN_cuts[Gen_bin] = trigger_Cut_scaled[Cut_bin];//need to multiply the corresponding scale factor, trigger_Cut_scaled[Cut_bin] = trigger_Cut[Cut_bin]*KNN_scale[Cut_bin];
                KNN_efficiency_cuts[Gen_bin] = efficiency_scaled[Cut_bin][Gen_bin];
            }//end if
            
            if (fabs(EMTF_efficiency_scaled[Cut_bin][Gen_bin]-efficiency_threshold)<EMTF_min) {
                EMTF_min = fabs(EMTF_efficiency_scaled[Cut_bin][Gen_bin]-efficiency_threshold);
                EMTF_cuts[Gen_bin] = EMTF_trigger_Cut_scaled[Cut_bin];
                EMTF_efficiency_cuts[Gen_bin] = EMTF_efficiency_scaled[Cut_bin][Gen_bin];
            }
        }//end Cut bin
        
    }//end Gen bin
    
    //fill 1D profile
    for (int Gen_bin=0;Gen_bin<59;Gen_bin++) {
        KNN_efficiency_cuts_consistency->Fill(GEN_pT[Gen_bin],KNN_efficiency_cuts[Gen_bin]);
        EMTF_efficiency_cuts_consistency->Fill(GEN_pT[Gen_bin],EMTF_efficiency_cuts[Gen_bin]);
        KNN_cuts_vs_Gen_pT->Fill(GEN_pT[Gen_bin],KNN_cuts[Gen_bin]);
        EMTF_cuts_vs_Gen_pT->Fill(GEN_pT[Gen_bin],EMTF_cuts[Gen_bin]);
    }
    
    KNN_trigger_Gen_efficiency_scaled->Write();
    EMTF_trigger_Gen_efficiency_scaled->Write();
    KNN_efficiency_cuts_consistency->Write();
    EMTF_efficiency_cuts_consistency->Write();
    KNN_cuts_vs_Gen_pT->Write();
    EMTF_cuts_vs_Gen_pT->Write();
    
    //===================================================
    //Enter test again for zero bias data to final rate plot
    //===================================================
    double rate_count[59]={0};
    double EMTF_rate_count[59]={0};
    
    std::cout << "\n******* Enter Zerobias test event to calculate rate *******" << std::endl;
    for (UInt_t iEvt = 0; iEvt < test_chain->GetEntries(); iEvt++) {
        
        //  if (iEvt > 10) break;
        //  if ( (iEvt % 1) == 0 )
        //std::cout << "\n*** Looking at event " << iEvt << " ***" <<  std::endl;
        
        test_chain->GetEntry(iEvt);
    
        //rate reduction using ZeroBias data in test sample, useful: eta==-99(623817),pt==999
        for (int Gen_bin =0;Gen_bin<59;Gen_bin++){
            if (GEN_pt_test_br == 999 && GEN_eta_test_br == -99){//both in zero bias
                if(1./KNN_test_br > KNN_cuts[Gen_bin]){
                    rate_count[Gen_bin]++;
                }
                if(EMTF_pt_test_br > EMTF_cuts[Gen_bin]){
                    EMTF_rate_count[Gen_bin]++;
                }
            }//end if Zerobias
            
        }//end for Gen_bin
    }// End loop: for (UInt_t iEvt = 0; iEvt < test_chain->GetEntries(); iEvt++)
    //std::cout << "\n******* Leaving the test event loop *******" << std::endl;
    
    //=========
    //rate plot
    //=========
    TGraph *KNN_rate = new TGraph(59,trigger_Cut,rate_count); KNN_rate->SetMarkerStyle(21); KNN_rate->SetMarkerColor(2);//red
    TGraph *EMTF_rate = new TGraph(59,trigger_Cut,EMTF_rate_count); EMTF_rate->SetMarkerStyle(21); EMTF_rate->SetMarkerColor(1);//black
    TCanvas *C_rate = new TCanvas("C_rate","Mode 14 rate",700,500);
    TMultiGraph *mg_rate = new TMultiGraph();
    C_rate->cd();
    mg_rate->SetTitle(Form("Mode 14 rate vs %f efficiency cut",efficiency_threshold));
    mg_rate->Add(KNN_rate);
    mg_rate->Add(EMTF_rate);
    mg_rate->Draw();
    mg_rate->Write();
    
    //==================
    //make log rate plot
    //==================
    double rate_count_log[59]={0};
    double EMTF_rate_count_log[59]={0};
    
    for (i=0;i<59;i++){rate_count_log[i] = log(rate_count[i])/log(10);};
    for (i=0;i<59;i++){EMTF_rate_count_log[i] = log(EMTF_rate_count[i])/log(10);};
    
    TGraph *KNN_rate_log = new TGraph(59,trigger_Cut,rate_count_log); KNN_rate_log->SetMarkerStyle(21); KNN_rate_log->SetMarkerColor(2);//red
    TGraph *EMTF_rate_log = new TGraph(59,trigger_Cut,EMTF_rate_count_log); EMTF_rate_log->SetMarkerStyle(21); EMTF_rate_log->SetMarkerColor(1);//black
    TCanvas *C_rate_log = new TCanvas("C_rate_log","Mode 14 log rate",700,500);
    TMultiGraph *mg_rate_log = new TMultiGraph();
    C_rate_log->cd();
    mg_rate_log->SetTitle(Form("Mode 14 log(rate)vs %f efficiency cut",efficiency_threshold));
    mg_rate_log->Add(KNN_rate_log);
    mg_rate_log->Add(EMTF_rate_log);
    mg_rate_log->Draw();
    mg_rate_log->Write();
    
    //==============
    //rate ratio plot
    //==============
    double rate_ratio[59]={0};
    for (i=0;i<59;i++){rate_ratio[i] = rate_count[i]/EMTF_rate_count[i];};
    TGraph *rate_ratio_plot = new TGraph(59,trigger_Cut,rate_ratio); rate_ratio_plot->SetMarkerStyle(21); rate_ratio_plot->SetMarkerColor(1);//black
    TCanvas *C_rate_ratio = new TCanvas("C_rate_ratio","Mode 14 rate ratio: KNN/EMTF",700,500);
    C_rate_ratio->cd();
    rate_ratio_plot->Draw();
    rate_ratio_plot->Write();
    
    //=====
    //test
    //=====
    Deviation_dPhiSum3_test->Write();
    Deviation_FR_1_test->Write();
    Deviation_dTh_13_test->Write();
    Deviation_bend_1_test->Write();
    Deviation_theta_test->Write();
    Deviation_dPhi_12_test->Write();
    Deviation_dPhi_23_test->Write();
    Deviation_dPhi_34_test->Write();
    Deviation_GEN_pt_test->Write();
    EMTF_Deviation_GEN_pt_test->Write();
    
    Deviation_dPhiSum3_test_1_8->Write();
    Deviation_FR_1_test_1_8->Write();
    Deviation_dTh_13_test_1_8->Write();
    Deviation_bend_1_test_1_8->Write();
    Deviation_theta_test_1_8->Write();
    Deviation_dPhi_12_test_1_8->Write();
    Deviation_dPhi_23_test_1_8->Write();
    Deviation_dPhi_34_test_1_8->Write();
    Deviation_GEN_pt_test_1_8->Write();
    EMTF_Deviation_GEN_pt_test_1_8->Write();
    
    Deviation_dPhiSum3_test_8_30->Write();
    Deviation_FR_1_test_8_30->Write();
    Deviation_dTh_13_test_8_30->Write();
    Deviation_bend_1_test_8_30->Write();
    Deviation_theta_test_8_30->Write();
    Deviation_dPhi_12_test_8_30->Write();
    Deviation_dPhi_23_test_8_30->Write();
    Deviation_dPhi_34_test_8_30->Write();
    Deviation_GEN_pt_test_8_30->Write();
    EMTF_Deviation_GEN_pt_test_8_30->Write();
    
    Deviation_dPhiSum3_test_30_120->Write();
    Deviation_FR_1_test_30_120->Write();
    Deviation_dTh_13_test_30_120->Write();
    Deviation_bend_1_test_30_120->Write();
    Deviation_theta_test_30_120->Write();
    Deviation_dPhi_12_test_30_120->Write();
    Deviation_dPhi_23_test_30_120->Write();
    Deviation_dPhi_34_test_30_120->Write();
    Deviation_GEN_pt_test_30_120->Write();
    EMTF_Deviation_GEN_pt_test_30_120->Write();
    
    Deviation_dPhiSum3_test_120_1000->Write();
    Deviation_FR_1_test_120_1000->Write();
    Deviation_dTh_13_test_120_1000->Write();
    Deviation_bend_1_test_120_1000->Write();
    Deviation_theta_test_120_1000->Write();
    Deviation_dPhi_12_test_120_1000->Write();
    Deviation_dPhi_23_test_120_1000->Write();
    Deviation_dPhi_34_test_120_1000->Write();
    Deviation_GEN_pt_test_120_1000->Write();
    EMTF_Deviation_GEN_pt_test_120_1000->Write();
    
    //=====
    //train
    //=====
    Deviation_dPhiSum3_train->Write();
    Deviation_FR_1_train->Write();
    Deviation_dTh_13_train->Write();
    Deviation_bend_1_train->Write();
    Deviation_theta_train->Write();
    Deviation_dPhi_12_train->Write();
    Deviation_dPhi_23_train->Write();
    Deviation_dPhi_34_train->Write();
    Deviation_GEN_pt_train->Write();
    EMTF_Deviation_GEN_pt_train->Write();
    
    Deviation_dPhiSum3_train_1_8->Write();
    Deviation_FR_1_train_1_8->Write();
    Deviation_dTh_13_train_1_8->Write();
    Deviation_bend_1_train_1_8->Write();
    Deviation_theta_train_1_8->Write();
    Deviation_dPhi_12_train_1_8->Write();
    Deviation_dPhi_23_train_1_8->Write();
    Deviation_dPhi_34_train_1_8->Write();
    Deviation_GEN_pt_train_1_8->Write();
    EMTF_Deviation_GEN_pt_train_1_8->Write();
    
    Deviation_dPhiSum3_train_8_30->Write();
    Deviation_FR_1_train_8_30->Write();
    Deviation_dTh_13_train_8_30->Write();
    Deviation_bend_1_train_8_30->Write();
    Deviation_theta_train_8_30->Write();
    Deviation_dPhi_12_train_8_30->Write();
    Deviation_dPhi_23_train_8_30->Write();
    Deviation_dPhi_34_train_8_30->Write();
    Deviation_GEN_pt_train_8_30->Write();
    EMTF_Deviation_GEN_pt_train_8_30->Write();
    
    Deviation_dPhiSum3_train_30_120->Write();
    Deviation_FR_1_train_30_120->Write();
    Deviation_dTh_13_train_30_120->Write();
    Deviation_bend_1_train_30_120->Write();
    Deviation_theta_train_30_120->Write();
    Deviation_dPhi_12_train_30_120->Write();
    Deviation_dPhi_23_train_30_120->Write();
    Deviation_dPhi_34_train_30_120->Write();
    Deviation_GEN_pt_train_30_120->Write();
    EMTF_Deviation_GEN_pt_train_30_120->Write();
    
    Deviation_dPhiSum3_train_120_1000->Write();
    Deviation_FR_1_train_120_1000->Write();
    Deviation_dTh_13_train_120_1000->Write();
    Deviation_bend_1_train_120_1000->Write();
    Deviation_theta_train_120_1000->Write();
    Deviation_dPhi_12_train_120_1000->Write();
    Deviation_dPhi_23_train_120_1000->Write();
    Deviation_dPhi_34_train_120_1000->Write();
    Deviation_GEN_pt_train_120_1000->Write();
    EMTF_Deviation_GEN_pt_train_120_1000->Write();
    
    //KNN
    //std::cout << "D of train: " << Deviation_train << std::endl;
    //std::cout << "D of test: " << Deviation_test << std::endl;
    std::cout << "In the order of train and test: overall train, test; 1-8 train, test; 8-30 train,test; 30-120 train, test; 120-1000 train,test;" << std::endl;
    std::cout << Standard_Deviation_train << std::endl;
    std::cout << Standard_Deviation_test << std::endl;
    std::cout << Standard_Deviation_train_1_8 << std::endl;
    std::cout << Standard_Deviation_test_1_8 << std::endl;
    std::cout << Standard_Deviation_train_8_30 << std::endl;
    std::cout << Standard_Deviation_test_8_30 << std::endl;
    std::cout << Standard_Deviation_train_30_120 << std::endl;
    std::cout << Standard_Deviation_test_30_120 << std::endl;
    std::cout << Standard_Deviation_train_120_1000 << std::endl;
    std::cout << Standard_Deviation_test_120_1000 << std::endl;
    
    //EMTF peformance
    //std::cout << "EMTF D of train: " << EMTF_Deviation_train << std::endl;
    //std::cout << "EMTF D of test: " << EMTF_Deviation_test << std::endl;
    std::cout << "****************" << std::endl;
    std::cout << "train events: " << train_events << std::endl;
    std::cout << "test events: " << test_events << std::endl;
    std::cout << "EMTF SD of train: " << EMTF_Standard_Deviation_train << std::endl;
    std::cout << "EMTF SD of test: " << EMTF_Standard_Deviation_test << std::endl;
    
    std::cout << "***1-8 GeV***" << std::endl;
    std::cout << "train events: " << train_events_1_8 << std::endl;
    std::cout << "test events: " << test_events_1_8 << std::endl;
    std::cout << "EMTF SD of train: " << EMTF_Standard_Deviation_train_1_8 << std::endl;
    std::cout << "EMTF SD of test: " << EMTF_Standard_Deviation_test_1_8 << std::endl;
    
    std::cout << "***8-30 GeV***" << std::endl;
    std::cout << "train events: " << train_events_8_30 << std::endl;
    std::cout << "test events: " << test_events_8_30 << std::endl;
    std::cout << "EMTF SD of train: " << EMTF_Standard_Deviation_train_8_30 << std::endl;
    std::cout << "EMTF SD of test: " << EMTF_Standard_Deviation_test_8_30 << std::endl;
    
    std::cout << "***30-120 GeV***" << std::endl;
    std::cout << "train events: " << train_events_30_120 << std::endl;
    std::cout << "test events: " << test_events_30_120 << std::endl;
    std::cout << "EMTF SD of train: " << EMTF_Standard_Deviation_train_30_120 << std::endl;
    std::cout << "EMTF SD of test: " << EMTF_Standard_Deviation_test_30_120 << std::endl;
    
    std::cout << "***120-1000 GeV***" << std::endl;
    std::cout << "train events: " << train_events_120_1000 << std::endl;
    std::cout << "test events: " << test_events_120_1000 << std::endl;
    std::cout << "EMTF SD of train: " << EMTF_Standard_Deviation_train_120_1000 << std::endl;
    std::cout << "EMTF SD of test: " << EMTF_Standard_Deviation_test_120_1000 << std::endl;
    
    
  //std::cout << "\nExiting ReadMVAOut_v1_KNN()\n";

} // End void ReadMVAOut()
