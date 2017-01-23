
/////////////////////////////////////////////////////////////////////
///  Simulataneous pT regression for multiple factories and MVAs  ///
///                  Andrew Brinkerhoff 23.01.17                  ///
///                                                               ///
///  Adapted from ROOT TMVARegression.C                           ///
///  Run using "root -l PtRegression_AWB_v1.C                     /// 
/////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"

// Extra tools - AWB 07.12.16
#include "interface/MVA_helper.h"
#include "src/MVA_input_var_tools.cc"

// // Class with NTuple branch definitions
// #include "interface/PtLutInputBranchesClass.hh"

const int MAX_EVT =    200000;
const int REPORT_EVT =  10000;
// const int MAX_EVT =   10000;
// const int REPORT_EVT = 1000;

const double PI = 3.14159265359;
const double PT_SCALE = 1.;
// const double PT_SCALE = (1. / 1.25); // EMTF pT was scaled up by ~1.25 in 2016 (for 4-hit tracks, mode 15)
const double BIT = 0.000001; // Tiny value or offset

using namespace TMVA;

void PtRegression_AWB_v1 ( TString myMethodList = "" ) {

   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDEFoam"]         = 0;
   Use["KNN"]             = 0;
   //
   // Linear Discriminant Analysis
   Use["LD"]		  = 0;
   //
   // Function Discriminant analysis
   Use["FDA_GA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   //
   // Neural Network
   Use["MLP"]             = 0;
   Use["DNN"]             = 0;
   //
   // Support Vector Machine
   Use["SVM"]             = 0;
   //
   // Boosted Decision Trees
   Use["BDT"]                     = 0;

   Use["BDTG_default"]            = 0;

   Use["BDTG_AWB"]                = 1;
   Use["BDTG_AWB_lite"]           = 1;

   Use["BDTG_AWB_64_trees"]       = 0;
   Use["BDTG_AWB_250_trees"]      = 0;
   Use["BDTG_AWB_1000_trees"]     = 0;

   Use["BDTG_AWB_4_deep"]         = 0;
   Use["BDTG_AWB_5_deep"]         = 0;
   Use["BDTG_AWB_6_deep"]         = 0;

   Use["BDTG_LeastSq"]            = 0;

   Use["BDTG_Carnes_AbsDev"]      = 0;
   Use["BDTG_Carnes_Huber"]       = 0;
   Use["BDTG_Carnes_LeastSq"]     = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start PtRegression_AWB_v1" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Here the preparation phase begins

   // Create a new root output file
   TString outfileName( "PtRegression_AWB_v1_17_01_23.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Read training and test data (see TMVAClassification for reading ASCII files)
   // load the signal and background event samples from ROOT trees
   TFile *input(0);
   TString store = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/ntuples";
   TString in_dir = "SingleMu_Pt1To1000_FlatRandomOneOverPt/EMTF_MuGun/170113_165434/0000";
   std::vector<TString> in_file_names;
   TString file_name;
   for (int i = 1; i < 99; i++) {
     file_name.Form("%s/%s/EMTF_MC_NTuple_SingleMu_noRPC_%d.root", store.Data(), in_dir.Data(), i);
     std::cout << "Adding file " << file_name.Data() << std::endl;
     in_file_names.push_back(file_name.Data());
     if (i*100000 > MAX_EVT) break; // ~100k events per file
   }

   for (UInt_t i = 0; i < in_file_names.size(); i++) {
     if ( !gSystem->AccessPathName(in_file_names.at(i)) )
       input = TFile::Open( in_file_names.at(i) ); // check if file in local directory exists
     if (!input) {
       std::cout << "ERROR: could not open data file " << in_file_names.at(i) << std::endl;
       exit(1);
     }
   }

   // Add trees from the input files to the TChain
   // Have to use TChain for both SetBranchAddress and GetEntry to work
   std::vector<TChain*> in_chains;
   // Super-hacky ... but using "GetBranch" with a single chain with multiple files causes a segfault - AWB 19.01.16
   for (UInt_t i = 0; i < in_file_names.size(); i++) {
     TChain *tmp_chain = new TChain("ntuple/tree");
     tmp_chain->Add( in_file_names.at(i) );
     in_chains.push_back(tmp_chain);
   }

   //////////////////////////////////////////////////////////////////////////
   ///  Factories: Use different sets of variables, target, weights, etc. ///
   //////////////////////////////////////////////////////////////////////////
   
   TString fact_set = "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression";
   std::vector<TString> var_names; // Holds names of variables for a given factory and permutation
   std::vector<Double_t> var_vals; // Holds values of variables for a given factory and permutation
   TMVA::Factory* nullF = new TMVA::Factory("NULL", outputFile, fact_set); // Placeholder factory
   TMVA::DataLoader* nullL = new TMVA::DataLoader("NULL");                 // Placeholder loader

   // Tuple is defined by the factory and dataloader,  followed by a name, 
   // var name and value vectors, and hex bit masks for input and target variables.
   // Each hex bit represents four variables, e.g. 0x1 would select only the 1st variable, 
   // 0xf the 1st 4, 0xff the 1st 8, 0xa the 2nd and 4th, 0xf1 the 1st and 5th-8th, etc.
   std::vector< std::tuple<TMVA::Factory*, TMVA::DataLoader*, TString, std::vector<TString>, std::vector<Double_t>, int, int> > factories;

   // Original set of training variables
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x0000011d_0x2", 
   					 var_names, var_vals, 0x0000011d, 0x2) );
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x0000011d_0x4", 
   					 var_names, var_vals, 0x0000011d, 0x4) );
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x0000011d_0x2_invPt", 
   					 var_names, var_vals, 0x0000011d, 0x2) );
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x0000011d_0x4_invPt", 
   					 var_names, var_vals, 0x0000011d, 0x4) );
   
   // Original set of training variables + combinations
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01fd_0x2", 
   					 var_names, var_vals, 0x001f01fd, 0x2) );
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01fd_0x4", 
   					 var_names, var_vals, 0x001f01fd, 0x4) );
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01fd_0x2_invPt", 
   					 var_names, var_vals, 0x001f01fd, 0x2) );
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01fd_0x4_invPt", 
   					 var_names, var_vals, 0x001f01fd, 0x4) );

   // Original set of training variables + combinations + St 1 ring, FR bits, bending
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x001fffff_0x2", 
   					 var_names, var_vals, 0x001fffff, 0x2) );
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x001fffff_0x4", 
   					 var_names, var_vals, 0x001fffff, 0x4) );
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x001fffff_0x2_invPt", 
   					 var_names, var_vals, 0x001fffff, 0x2) );
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x001fffff_0x4_invPt", 
   					 var_names, var_vals, 0x001fffff, 0x4) );

   // Initialize factories and dataloaders
   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
     std::get<0>(factories.at(iFact)) = new TMVA::Factory( std::get<2>(factories.at(iFact)), outputFile, fact_set );
     std::get<1>(factories.at(iFact)) = new TMVA::DataLoader( std::get<2>(factories.at(iFact)) );
   }

   // Defined in interface/MVA_helper.h
   // MVA_var(TString name, TString descr, TString unit, TString type, Double_t def_val)
   std::vector<MVA_var> in_vars;   // All input variables
   std::vector<MVA_var> targ_vars; // All target variables (should only define 1, unless using MLP)
   std::vector<MVA_var> spec_vars; // All spectator variables
   std::vector<MVA_var> all_vars;  // All variables
   
   /////////////////////////////////////////////////////////
   ///  Input variables: used in BDT to estimate the pT  ///
   /////////////////////////////////////////////////////////
   
   in_vars.push_back( MVA_var( "theta",     "Track #theta",      "deg", 'F', -88 ) ); // 0x0000 0001  * 2016 variable
   in_vars.push_back( MVA_var( "St1_ring",  "St 1 LCT ring",     "int", 'I', -88 ) ); // 0x0000 0002  
   in_vars.push_back( MVA_var( "dPhi_12",   "#phi(2) - #phi(1)", "deg", 'F', -88 ) ); // 0x0000 0004  * 2016 variable
   in_vars.push_back( MVA_var( "dPhi_23",   "#phi(3) - #phi(2)", "deg", 'F', -88 ) ); // 0x0000 0008  * 2016 variable

   in_vars.push_back( MVA_var( "dPhi_34",   "#phi(4) - #phi(3)", "deg", 'F', -88 ) ); // 0x0000 0010  * 2016 variable
   in_vars.push_back( MVA_var( "dPhi_13",   "#phi(3) - #phi(1)", "deg", 'F', -88 ) ); // 0x0000 0020  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "dPhi_14",   "#phi(4) - #phi(1)", "deg", 'F', -88 ) ); // 0x0000 0040  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "dPhi_24",   "#phi(4) - #phi(2)", "deg", 'F', -88 ) ); // 0x0000 0080  * Derivable from 2016 var

   in_vars.push_back( MVA_var( "FR_1",      "St 1 LCT F/R",      "int", 'I', -88 ) ); // 0x0000 0100  * 2016 variable
   in_vars.push_back( MVA_var( "FR_2",      "St 2 LCT F/R",      "int", 'I', -88 ) ); // 0x0000 0200
   in_vars.push_back( MVA_var( "FR_3",      "St 3 LCT F/R",      "int", 'I', -88 ) ); // 0x0000 0400
   in_vars.push_back( MVA_var( "FR_4",      "St 4 LCT F/R",      "int", 'I', -88 ) ); // 0x0000 0800

   in_vars.push_back( MVA_var( "bend_1",    "St 1 LCT bending",  "int", 'I', -88 ) ); // 0x0000 1000
   in_vars.push_back( MVA_var( "bend_2",    "St 2 LCT bending",  "int", 'I', -88 ) ); // 0x0000 2000
   in_vars.push_back( MVA_var( "bend_3",    "St 3 LCT bending",  "int", 'I', -88 ) ); // 0x0000 4000
   in_vars.push_back( MVA_var( "bend_4",    "St 4 LCT bending",  "int", 'I', -88 ) ); // 0x0000 8000

   in_vars.push_back( MVA_var( "dPhiSum4",  "#Sigmad#phi (6)",   "deg", 'F', -88 ) ); // 0x0001 0000  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "dPhiSum4A", "#Sigma|d#phi| (6)", "deg", 'F', -88 ) ); // 0x0002 0000  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "dPhiSum3",  "#Sigmad#phi (3)",   "deg", 'F', -88 ) ); // 0x0004 0000  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "dPhiSum3A", "#Sigma|d#phi| (3)", "deg", 'F', -88 ) ); // 0x0008 0000  * Derivable from 2016 var
   
   in_vars.push_back( MVA_var( "outStPhi",  "#phi outlier St",   "int", 'I', -88 ) ); // 0x0010 0000  * Derivable from 2016 var
   // // Not yet computed - AWB 23.01.17
   // in_vars.push_back( MVA_var( "outStTh",   "#theta outlier St", "int", 'I', -88 ) ); // 0x0020 0000
   // in_vars.push_back( MVA_var( "dThMax4",   "Max d#theta (6)",   "deg", 'F', -88 ) ); // 0x0040 0000
   // in_vars.push_back( MVA_var( "dThMax3",   "Max d#theta (3)",   "deg", 'F', -88 ) ); // 0x0080 0000

   // in_vars.push_back( MVA_var( "dThSum4",   "#Sigmad#phi (6)",   "deg", 'F', -88 ) ); // 0x0100 0000
   // in_vars.push_back( MVA_var( "dThSum4A",  "#Sigma|d#phi| (6)", "deg", 'F', -88 ) ); // 0x0200 0000
   // in_vars.push_back( MVA_var( "dThSum3",   "#Sigmad#phi (3)",   "deg", 'F', -88 ) ); // 0x0400 0000
   // in_vars.push_back( MVA_var( "dThSum3A",  "#Sigma|d#phi| (3)", "deg", 'F', -88 ) ); // 0x0800 0000

   ////////////////////////////////////////////////////////////
   //  Target variable: true muon pT, or 1/pT, or log2(pT)  ///
   ////////////////////////////////////////////////////////////
   
   targ_vars.push_back( MVA_var( "GEN_pt",       "GEN p_{T}",               "GeV",      'F', -99 ) ); // 0x1
   targ_vars.push_back( MVA_var( "inv_GEN_pt",   "1 / GEN muon p_{T}",      "GeV^{-1}", 'F', -99 ) ); // 0x2
   targ_vars.push_back( MVA_var( "log2_GEN_pt",  "log_{2}(GEN muon p_{T})", "GeV",      'F', -99 ) ); // 0x4

   /////////////////////////////////////////////////////////////////////////////
   ///  Spectator variables: not used in training, but saved in output tree  ///
   /////////////////////////////////////////////////////////////////////////////
   
   spec_vars.push_back( MVA_var( "GEN_pt",       "GEN p_{T}",               "GeV",      'F', -77 ) ); // 0x01
   spec_vars.push_back( MVA_var( "EMTF_pt",      "EMTF p_{T}",              "GeV",      'F', -77 ) ); // 0x02
   spec_vars.push_back( MVA_var( "inv_GEN_pt",   "1 / GEN muon p_{T}",      "GeV^{-1}", 'F', -77 ) ); // 0x04
   spec_vars.push_back( MVA_var( "inv_EMTF_pt",  "1 / EMTF p_{T}",          "GeV^{-1}", 'F', -77 ) ); // 0x08
   spec_vars.push_back( MVA_var( "log2_GEN_pt",  "log_{2}(GEN muon p_{T})", "GeV",      'F', -77 ) ); // 0x10
   spec_vars.push_back( MVA_var( "log2_EMTF_pt", "log_{2}(EMTF p_{T})",     "GeV",      'F', -77 ) ); // 0x20
   spec_vars.push_back( MVA_var( "GEN_eta",      "GEN #eta",                "",         'F', -77 ) ); // 0x40
   spec_vars.push_back( MVA_var( "EMTF_eta",     "EMTF #eta",               "",         'F', -77 ) ); // 0x80
   spec_vars.push_back( MVA_var( "GEN_charge",   "GEN charge",              "",         'I', -77 ) ); // 0x40
   spec_vars.push_back( MVA_var( "EMTF_charge",  "EMTF charge",             "",         'I', -77 ) ); // 0x80
   spec_vars.push_back( MVA_var( "dPhi_12_sign", "#phi(2) - #phi(1) sign",  "",         'I', -77 ) ); // 0x80

   assert( in_vars.size() > 0 );   // You need at least one input variable
   assert( targ_vars.size() > 0 ); // You need at least one target variable
   // Order is important: input variables first, then target, then specator
   all_vars.insert( all_vars.end(), in_vars.begin(), in_vars.end() );
   all_vars.insert( all_vars.end(), targ_vars.begin(), targ_vars.end() );
   all_vars.insert( all_vars.end(), spec_vars.begin(), spec_vars.end() );

   // Fill each factory with the correct set of variables
   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
     std::cout << "\n*** Factory " << std::get<2>(factories.at(iFact)) << " variables ***" << std::endl;
       
     std::cout << "*** Input ***" << std::endl;
     for (UInt_t i = 0; i < in_vars.size(); i++) {
       if ( 0x1 & (std::get<5>(factories.at(iFact)) >> i) ) { // Hex bit mask for in_vars
	 MVA_var v = in_vars.at(i);
	 std::cout << v.name << std::endl;
	 std::get<1>(factories.at(iFact))->AddVariable( v.name, v.descr, v.unit, v.type ); // Add var to dataloader 
	 std::get<3>(factories.at(iFact)).push_back( v.name );    // Add to vector of var names
	 std::get<4>(factories.at(iFact)).push_back( v.def_val ); // Add to vector of var values
       }
     }

     TString targ_str = ""; // Save name of target variable
     std::cout << "*** Target ***" << std::endl;
     for (UInt_t i = 0; i < targ_vars.size(); i++) {
       if ( 0x1 & (std::get<6>(factories.at(iFact)) >> i) ) { // Hex bit mask for targ_vars
	 MVA_var v = targ_vars.at(i);
	 std::cout << v.name << std::endl;
	 targ_str = v.name;
	 std::get<1>(factories.at(iFact))->AddTarget( v.name, v.descr, v.unit, v.type );
	 std::get<3>(factories.at(iFact)).push_back( v.name );
	 std::get<4>(factories.at(iFact)).push_back( v.def_val );
       }
     }

     std::cout << "*** Spectator ***" << std::endl;
     for (UInt_t i = 0; i < spec_vars.size(); i++) {
       MVA_var v = spec_vars.at(i);
       if (v.name == targ_str) continue; // Don't add target variable
       std::cout << v.name << std::endl;
       std::get<1>(factories.at(iFact))->AddSpectator( v.name, v.descr, v.unit, v.type );
       std::get<3>(factories.at(iFact)).push_back( v.name );
       std::get<4>(factories.at(iFact)).push_back( v.def_val );
     }
   } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++)


   std::cout << "\n******* About to loop over chains *******" << std::endl;
   UInt_t iEvt = 0;
   UInt_t nTrain = 0;
   UInt_t nTest  = 0;
   for (int iCh = 0; iCh < in_chains.size(); iCh++) {
     TChain *in_chain = in_chains.at(iCh);
     
     // Get branches from the chain
     TBranch *muon_br  = in_chain->GetBranch("muon");
     TBranch *hit_br   = in_chain->GetBranch("hit");
     TBranch *track_br = in_chain->GetBranch("track");
     
     std::cout << "\n******* About to enter the event loop for chain " << iCh+1 << " *******" << std::endl;
     
     for (UInt_t jEvt = 0; jEvt < in_chain->GetEntries(); jEvt++) {
       if (iEvt > MAX_EVT) break;
       if ( (iEvt % REPORT_EVT) == 0 )
	 std::cout << "Looking at event " << iEvt << std::endl;
       in_chain->GetEntry(jEvt);
       
       UInt_t nMuons  = (muon_br->GetLeaf("nMuons"))->GetValue();
       UInt_t nTracks = (track_br->GetLeaf("nTracks"))->GetValue();
       // std::cout << "There are " << nMuons << " GEN muons and " << nTracks << " EMTF tracks\n" << std::endl;
       
       for (UInt_t iMu = 0; iMu < nMuons; iMu++) {
	 Double_t mu_pt  = (muon_br->GetLeaf("pt"))->GetValue(iMu);
	 Double_t mu_eta = (muon_br->GetLeaf("eta"))->GetValue(iMu);
	 Double_t mu_phi = (muon_br->GetLeaf("phi"))->GetValue(iMu);
	 Int_t mu_charge = (muon_br->GetLeaf("charge"))->GetValue(iMu);
	 if ( fabs( mu_eta ) < 1.10 ) continue;
	 if ( fabs( mu_eta ) > 2.50 ) continue;
	 // if ( mu_pt < 10 ) continue;
	 // std::cout << "Muon " << iMu+1 << " has pt = " << mu_pt << ", eta = " << mu_eta << ", phi = " << mu_phi << std::endl;
	 
	 for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++) {
	   Double_t trk_pt    = (track_br->GetLeaf("pt"))->GetValue(iTrk);
	   Double_t trk_eta   = (track_br->GetLeaf("eta"))->GetValue(iTrk);
	   trk_pt *= PT_SCALE;
	   Double_t trk_theta = (track_br->GetLeaf("theta"))->GetValue(iTrk);
	   Double_t trk_phi   = (track_br->GetLeaf("phi"))->GetValue(iTrk);
	   Int_t trk_charge   = (track_br->GetLeaf("charge"))->GetValue(iTrk);
	   Int_t trk_mode     = (track_br->GetLeaf("mode"))->GetValue(iTrk);
	   if ( trk_mode != 15 ) continue;
	   if ( ( mu_eta > 0 ) != ( trk_eta > 0 ) ) continue;
	   // std::cout << "  * Track " << iTrk+1 << " has pt = " << trk_pt << ", eta = " << trk_eta 
	   // 	   << ", phi = " << mu_phi << ", mode = " << trk_mode << std::endl;
	   
	   Int_t st1_ring = ((int) (track_br->GetLeaf("hit_ring"))->GetValue(4*iTrk + 0)) % 3; // Ring 4 --> Ring 1 (ME1/1a)

	   Double_t phi1 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk + 0);
	   Double_t phi2 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk + 1);
	   Double_t phi3 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk + 2);
	   Double_t phi4 = (track_br->GetLeaf("hit_phi"))->GetValue(4*iTrk + 3);
	   // std::cout << "  * Hits have phi = " << phi1 << ", " << phi2 << ", " << phi3 << ", " << phi4 << std::endl;
	   
	   Double_t dPhi12 = acos( cos( (phi2 - phi1)*(PI/180.) ) );
	   dPhi12 *= ( sin( (phi2 - phi1)*(PI/180.) ) / max( BIT, abs( sin( (phi2 - phi1)*(PI/180.) ) ) ) );
	   dPhi12 *= (180./PI);
	   Double_t dPhi13 = acos( cos( (phi3 - phi1)*(PI/180.) ) );
	   dPhi13 *= ( sin( (phi3 - phi1)*(PI/180.) ) / max( BIT, abs( sin( (phi3 - phi1)*(PI/180.) ) ) ) );
	   dPhi13 *= (180./PI);
	   Double_t dPhi14 = acos( cos( (phi4 - phi1)*(PI/180.) ) );
	   dPhi14 *= ( sin( (phi4 - phi1)*(PI/180.) ) / max( BIT, abs( sin( (phi4 - phi1)*(PI/180.) ) ) ) );
	   dPhi14 *= (180./PI);
	   Double_t dPhi23 = acos( cos( (phi3 - phi2)*(PI/180.) ) );
	   dPhi23 *= ( sin( (phi3 - phi2)*(PI/180.) ) / max( BIT, abs( sin( (phi3 - phi2)*(PI/180.) ) ) ) );
	   dPhi23 *= (180./PI);
	   Double_t dPhi24 = acos( cos( (phi4 - phi2)*(PI/180.) ) );
	   dPhi24 *= ( sin( (phi4 - phi2)*(PI/180.) ) / max( BIT, abs( sin( (phi4 - phi2)*(PI/180.) ) ) ) );
	   dPhi24 *= (180./PI);
	   Double_t dPhi34 = acos( cos( (phi4 - phi3)*(PI/180.) ) );
	   dPhi34 *= ( sin( (phi4 - phi3)*(PI/180.) ) / max( BIT, abs( sin( (phi4 - phi3)*(PI/180.) ) ) ) );
	   dPhi34 *= (180./PI);
	   
	   // Define all dPhi values relative to dPhi12
	   Int_t dPhi12_sign = ( (dPhi12 < 0) ? -1. : 1. );
	   dPhi13 *= dPhi12_sign;
	   dPhi14 *= dPhi12_sign;
	   dPhi23 *= dPhi12_sign;
	   dPhi24 *= dPhi12_sign;
	   dPhi34 *= dPhi12_sign;
	   dPhi12  = fabs(dPhi12);
	   
	   Double_t dPhi_sum_4  = dPhi12 + dPhi13 + dPhi14 + dPhi23 + dPhi24 + dPhi34;
	   Double_t dPhi_sum_4A = fabs(dPhi12) + fabs(dPhi13) + fabs(dPhi14) + fabs(dPhi23) + fabs(dPhi24) + fabs(dPhi34);
	   Double_t dev_st1 = fabs(dPhi12) + fabs(dPhi13) + fabs(dPhi14);
	   Double_t dev_st2 = fabs(dPhi12) + fabs(dPhi23) + fabs(dPhi24);
	   Double_t dev_st3 = fabs(dPhi13) + fabs(dPhi23) + fabs(dPhi34);
	   Double_t dev_st4 = fabs(dPhi14) + fabs(dPhi24) + fabs(dPhi34);
	   
	   Int_t out_st_phi = -88;
	   if      (dev_st4 > dev_st3 && dev_st4 > dev_st2 && dev_st4 > dev_st1)  out_st_phi = 4;
	   else if (dev_st3 > dev_st4 && dev_st3 > dev_st2 && dev_st3 > dev_st1)  out_st_phi = 3;
	   else if (dev_st2 > dev_st4 && dev_st2 > dev_st3 && dev_st2 > dev_st1)  out_st_phi = 2;
	   else if (dev_st1 > dev_st4 && dev_st1 > dev_st3 && dev_st1 > dev_st2)  out_st_phi = 1;
	   else                                                                   out_st_phi = 0;
	   
	   Double_t dPhi_sum_3  = -88;
	   Double_t dPhi_sum_3A = -88;
	   if      (out_st_phi == 4) {
	     dPhi_sum_3  = dPhi12 + dPhi13 + dPhi23;
	     dPhi_sum_3A = fabs(dPhi12) + fabs(dPhi13) + fabs(dPhi23);
	   } else if (out_st_phi == 3) {
	     dPhi_sum_3  = dPhi12 + dPhi14 + dPhi24;
	     dPhi_sum_3A = fabs(dPhi12) + fabs(dPhi14) + fabs(dPhi24);
	   } else if (out_st_phi == 2) {
	     dPhi_sum_3  = dPhi13 + dPhi14 + dPhi34;
	     dPhi_sum_3A = fabs(dPhi13) + fabs(dPhi14) + fabs(dPhi34);
	   } else {
	     dPhi_sum_3  = dPhi23 + dPhi24 + dPhi34;
	     dPhi_sum_3A = fabs(dPhi23) + fabs(dPhi24) + fabs(dPhi34);
	   }
	   
	   Int_t FR1 = (track_br->GetLeaf("hit_FR"))->GetValue(4*iTrk + 0);
	   Int_t FR2 = (track_br->GetLeaf("hit_FR"))->GetValue(4*iTrk + 1);
	   Int_t FR3 = (track_br->GetLeaf("hit_FR"))->GetValue(4*iTrk + 2);
	   Int_t FR4 = (track_br->GetLeaf("hit_FR"))->GetValue(4*iTrk + 3);
	   
	   Int_t bend1 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 0), trk_eta ) * dPhi12_sign;
	   Int_t bend2 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 1), trk_eta ) * dPhi12_sign;
	   Int_t bend3 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 2), trk_eta ) * dPhi12_sign;
	   Int_t bend4 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 3), trk_eta ) * dPhi12_sign;
	   
	   /////////////////////////////////////////////////////
	   ///  Loop over factories and set variable values  ///
	   /////////////////////////////////////////////////////
	   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
	     
	     // Set vars equal to default vector of variables for this factory
	     var_names = std::get<3>(factories.at(iFact));
	     var_vals = std::get<4>(factories.at(iFact));
	     
	     // Fill all variables
	     for (UInt_t iVar = 0; iVar < var_names.size(); iVar++) {
	       TString vName = var_names.at(iVar);
	       
	       /////////////////////////
	       ///  Input variables  ///
	       /////////////////////////
	       
	       if ( vName == "theta" )
		 var_vals.at(iVar) = fabs(trk_theta);
	       if ( vName == "St1_ring" )
		 var_vals.at(iVar) = st1_ring;
	       if ( vName == "dPhi_12" )
		 var_vals.at(iVar) = dPhi12;
	       if ( vName == "dPhi_12_sign" )
		 var_vals.at(iVar) = dPhi12_sign;
	       if ( vName == "dPhi_13" )
		 var_vals.at(iVar) = dPhi13;
	       if ( vName == "dPhi_14" )
		 var_vals.at(iVar) = dPhi14;
	       if ( vName == "dPhi_23" )
		 var_vals.at(iVar) = dPhi23;
	       if ( vName == "dPhi_24" )
		 var_vals.at(iVar) = dPhi24;
	       if ( vName == "dPhi_34" )
		 var_vals.at(iVar) = dPhi34;
	       
	       if ( vName == "FR_1" )
		 var_vals.at(iVar) = FR1;
	       if ( vName == "FR_2" )
		 var_vals.at(iVar) = FR2;
	       if ( vName == "FR_3" )
		 var_vals.at(iVar) = FR3;
	       if ( vName == "FR_4" )
		 var_vals.at(iVar) = FR4;
	       
	       if ( vName == "bend_1" )
		 var_vals.at(iVar) = bend1;
	       if ( vName == "bend_2" )
		 var_vals.at(iVar) = bend2;
	       if ( vName == "bend_3" )
		 var_vals.at(iVar) = bend3;
	       if ( vName == "bend_4" )
		 var_vals.at(iVar) = bend4;
	       
	       if ( vName == "dPhiSum4" )
		 var_vals.at(iVar) = dPhi_sum_4;
	       if ( vName == "dPhiSum4A" )
		 var_vals.at(iVar) = dPhi_sum_4A;
	       if ( vName == "dPhiSum3" )
		 var_vals.at(iVar) = dPhi_sum_3;
	       if ( vName == "dPhiSum3A" )
		 var_vals.at(iVar) = dPhi_sum_3A;
	       if ( vName == "outStPhi" )
		 var_vals.at(iVar) = out_st_phi;
	       
	       
	       ////////////////////////////////////////
	       ///  Target and spectator variables  ///
	       ////////////////////////////////////////
	       
	       if ( vName == "GEN_pt" )
		 var_vals.at(iVar) = mu_pt;
	       if ( vName == "EMTF_pt" )
		 var_vals.at(iVar) = trk_pt;
	       if ( vName == "inv_GEN_pt" )
		 var_vals.at(iVar) = 1. / mu_pt;
	       if ( vName == "inv_EMTF_pt" )
		 var_vals.at(iVar) = 1. / trk_pt;
	       if ( vName == "log2_GEN_pt" )
		 var_vals.at(iVar) = log2(mu_pt);
	       if ( vName == "log2_EMTF_pt" )
		 var_vals.at(iVar) = log2(trk_pt);
	       if ( vName == "GEN_eta" )
		 var_vals.at(iVar) = mu_eta;
	       if ( vName == "EMTF_eta" )
		 var_vals.at(iVar) = trk_eta;
	       if ( vName == "GEN_charge" )
		 var_vals.at(iVar) = mu_charge;
	       if ( vName == "EMTF_charge" )
		 var_vals.at(iVar) = trk_charge;
	       
	     } // End loop: for (UInt_t iVar = 0; iVar < var_names.size(); iVar++)
	     
	     // Unweighted distribution: flat in eta and 1/pT
	     Double_t evt_weight = 1.0;
	     // Weight by 1/pT so overall distribution is (1/pT)^2
	     if ( std::get<2>(factories.at(iFact)).Contains("_invPt") )
	       evt_weight = 1. / mu_pt;
	     
	     // Load values into event
	     if ( (iEvt % 2) == 0 ) {
	       std::get<1>(factories.at(iFact))->AddTrainingEvent( "Regression", var_vals, evt_weight );
	       if (iFact == 0) nTrain += 1;
	     }
	     else {
	       std::get<1>(factories.at(iFact))->AddTestEvent( "Regression", var_vals, evt_weight );
	       if (iFact == 0) nTest += 1;
	     }
	     
	   } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++) 
	   
	 } // End loop: for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++)
       } // End loop: for (UInt_t iMu = 0; iMu < nMuons; iMu++)
       iEvt += 1;
     } // End loop: for (UInt_t jEvt = 0; jEvt < in_chain->GetEntries(); jEvt++)
   } // End loop: for (int iCh = 0; iCh < in_chains.size(); iCh++) {

   std::cout << "******* Made it out of the event loop *******" << std::endl;
   
   string NTr;
   string NTe;

   ostringstream convertTr;
   convertTr << nTrain;
   NTr = convertTr.str();
   ostringstream convertTe;
   convertTe << nTest;
   NTe = convertTe.str();

   string numTrainStr = "nTrain_Regression="+NTr+":nTest_Regression="+NTe+":";
   std::cout << "NTr: " << NTr << ", NTe: " << NTe << std::endl;

   // // global event weights per tree (see below for setting event-wise weights)
   // Double_t regWeight  = 1.0;

   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
     
     TMVA::Factory* factX = std::get<0>(factories.at(iFact));
     TMVA::DataLoader* loadX = std::get<1>(factories.at(iFact));
     
     // // You can add an arbitrary number of regression trees
     // loadX->AddRegressionTree( regTree, regWeight );
     
     // // This would set individual event weights (the variables defined in the
     // // expression need to exist in the original TTree)
     // loadX->SetWeightExpression( "var1", "Regression" );
     loadX->SetWeightExpression( 1.0 );
     
     // // Apply additional cuts on the signal and background samples (can be different)
     // TCut mycut = "( abs(muon.eta[0]) > 1.25 && abs(muon.eta[1]) < 2.4 )"; // && track.mode[0] == 15 )"; 
     
     // Set nTest_Regression to 0 to tell the DataLoader to use all remaining events in the trees after training for testing:
     loadX->PrepareTrainingAndTestTree( "", numTrainStr+"SplitMode=Random:NormMode=NumEvents:!V" );   
     // loadX->PrepareTrainingAndTestTree( mycut, "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );
     
     // If no numbers of events are given, half of the events in the tree are used
     // for training, and the other half for testing:
     //
     //     loadX->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
     
     // Book MVA methods
     //
     // Please lookup the various method configuration options in the corresponding cxx files, eg:
     // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
     // it is possible to preset ranges in the option string in which the cut optimisation should be done:
     // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable
     
     // Linear discriminant
     if (Use["LD"])
       factX->BookMethod( loadX,  TMVA::Types::kLD, "LD",
			  "!H:!V:VarTransform=None" );
     
     // Neural network (MLP)
     if (Use["MLP"])
       factX->BookMethod( loadX,  TMVA::Types::kMLP, "MLP", (string)
			  "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:"+
			  "TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:"+
			  "ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );
     
     if (Use["DNN"])
       {
	 
	 // TString layoutString ("Layout=TANH|(N+100)*2,LINEAR");
	 // TString layoutString ("Layout=SOFTSIGN|100,SOFTSIGN|50,SOFTSIGN|20,LINEAR");
	 // TString layoutString ("Layout=RELU|300,RELU|100,RELU|30,RELU|10,LINEAR");
	 // TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|30,SOFTSIGN|20,SOFTSIGN|10,LINEAR");
	 // TString layoutString ("Layout=TANH|50,TANH|30,TANH|20,TANH|10,LINEAR");
	 // TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|20,LINEAR");
	 // TString layoutString ("Layout=TANH|100,TANH|30,LINEAR");

	 TString layoutString ("Layout=TANH|100,LINEAR");
	 
	 TString training0 ( (string) "LearningRate=1e-5,Momentum=0.5,Repetitions=1,"+
			     "ConvergenceSteps=500,BatchSize=50,TestRepetitions=7,WeightDecay=0.01,"+
			     "Regularization=NONE,DropConfig=0.5+0.5+0.5+0.5,DropRepetitions=2");
	 TString training1 ( (string) "LearningRate=1e-5,Momentum=0.9,Repetitions=1,"+
			     "ConvergenceSteps=170,BatchSize=30,TestRepetitions=7,WeightDecay=0.01,"+
			     "Regularization=L2,DropConfig=0.1+0.1+0.1,DropRepetitions=1");
	 TString training2 ( (string) "LearningRate=1e-5,Momentum=0.3,Repetitions=1,ConvergenceSteps=150,"+
			     "BatchSize=40,TestRepetitions=7,WeightDecay=0.01,Regularization=NONE");
	 TString training3 ( (string) "LearningRate=1e-6,Momentum=0.1,Repetitions=1,ConvergenceSteps=500,"+
			     "BatchSize=100,TestRepetitions=7,WeightDecay=0.0001,Regularization=NONE");
	 
	 TString trainingStrategyString ("TrainingStrategy=");
	 trainingStrategyString += training0 + "|" + training1 + "|" + training2 + "|" + training3;
	 
	 
	 // TString trainingStrategyString ( (string) "TrainingStrategy=LearningRate=1e-1,Momentum=0.3,"+
	 // 				  "Repetitions=3,ConvergenceSteps=20,BatchSize=30,TestRepetitions=7,"+
	 // 				  "WeightDecay=0.0,L1=false,DropFraction=0.0,DropRepetitions=5");
	 
	 TString nnOptions ("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G:WeightInitialization=XAVIERUNIFORM");
	 // TString nnOptions ("!H:V:VarTransform=Normalize:ErrorStrategy=CHECKGRADIENTS");
	 nnOptions.Append (":"); nnOptions.Append (layoutString);
	 nnOptions.Append (":"); nnOptions.Append (trainingStrategyString);
	 
	 factX->BookMethod(loadX, TMVA::Types::kDNN, "DNN", nnOptions ); // NN
       }


     // Support Vector Machine
     if (Use["SVM"])
       factX->BookMethod( loadX,  TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
     
     // Boosted Decision Trees
     if (Use["BDT"])
       factX->BookMethod( loadX,  TMVA::Types::kBDT, "BDT", (string)
			  "!H:!V:NTrees=100:MinNodeSize=1.0%:BoostType=AdaBoostR2:SeparationType=RegressionVariance"+
			  ":nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );
     
     // Default TMVA settings
     if (Use["BDTG_default"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_default", (string)
			  "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			  "BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );

     // AWB settings - AbsoluteDeviation
     if (Use["BDTG_AWB"]) // Optimized settings
       // Actually, would prefer NTrees=400 and MinNodeSize=0.001 ... but causing "all events went to the same branch" issue
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB", (string)
			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.01:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );

     if (Use["BDTG_AWB_lite"]) // Fast, simple BDT
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_lite", (string)
			  "!H:!V:NTrees=40::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.01:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );

     if (Use["BDTG_AWB_64_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_64_trees", (string)
			  "!H:!V:NTrees=64::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_250_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_250_trees", (string)
			  "!H:!V:NTrees=250::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_1000_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_1000_trees", (string)
			  "!H:!V:NTrees=1000::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     
     if (Use["BDTG_AWB_4_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_4_deep", (string)
			  "!H:!V:NTrees=128::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=4:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_5_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_5_deep", (string)
			  "!H:!V:NTrees=128::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_6_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_6_deep", (string)
			  "!H:!V:NTrees=128::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=6:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     
     // Default TMVA settings with LeastSquares loss function
     if (Use["BDTG_LeastSq"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_LeastSq", (string)
			  "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:"+
			  "BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:"+
			  "RegressionLossFunctionBDTG=LeastSquares");
     
     // Factory settings from Andrew Carnes ... what do they do? - AWB 04.01.17
     if (Use["BDTG_Carnes_AbsDev"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_Carnes_AbsDev", (string)
 			  "!H:!V:NTrees=64::BoostType=Grad:Shrinkage=0.3:nCuts=99999:MaxDepth=4:MinNodeSize=0.001:"+
			  "NegWeightTreatment=IgnoreNegWeightsInTraining:PruneMethod=NoPruning:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     
     if (Use["BDTG_Carnes_Huber"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_Carnes_Huber", (string)
 			  "!H:!V:NTrees=64::BoostType=Grad:Shrinkage=0.3:nCuts=99999:MaxDepth=4:MinNodeSize=0.001:"+
			  "NegWeightTreatment=IgnoreNegWeightsInTraining:PruneMethod=NoPruning:"+
			  "RegressionLossFunctionBDTG=Huber" );
     
     if (Use["BDTG_Carnes_LeastSq"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_Carnes_LeastSq", (string)
 			  "!H:!V:NTrees=64::BoostType=Grad:Shrinkage=0.3:nCuts=99999:MaxDepth=4:MinNodeSize=0.001:"+
			  "NegWeightTreatment=IgnoreNegWeightsInTraining:PruneMethod=NoPruning:"+
			  "RegressionLossFunctionBDTG=LeastSquares" );
     
     
     // --------------------------------------------------------------------------------------------------
     
     // Now you can tell the factory to train, test, and evaluate the MVAs
     
     // Train MVAs using the set of training events
     factX->TrainAllMethods();
     
     // Evaluate all MVAs using the set of test events
     factX->TestAllMethods();
     
     // Evaluate and compare performance of all configured MVAs
     factX->EvaluateAllMethods();
     
     // --------------------------------------------------------------
     
   } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++)
   
   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;

   // delete factory;
   // delete dataloader;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVARegGui( outfileName );
}


int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   PtRegression_AWB_v1(methodList);
   return 0;
}

