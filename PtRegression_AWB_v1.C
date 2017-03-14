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


// Typical settings with ZeroBias, without RPC
const int MAX_EVT    = 5000000; // Maximum number of MuGun events to process
const int MAX_TR     =  500000; // Maximum number of MuGun training events
const int REPORT_EVT =  100000;

// // // Typical settings without ZeroBias, with RPC
// const int MAX_EVT    = 30000000; // Maximum number of MuGun events to process
// const int MAX_TR     = 20000000; // Maximum number of MuGun training events
// const int REPORT_EVT =   100000;

// // Typical test settings
// const int MAX_EVT    = 20000;
// const int MAX_TR     = 10000;
// const int REPORT_EVT =  1000;

const double PI = 3.14159265359;
const double PT_SCALE = 1.;
// const double PT_SCALE = (1. / 1.25); // EMTF pT was scaled up by ~1.25 in 2016 (for 4-hit tracks, mode 15)
const double BIT = 0.000001; // Tiny value or offset

const double PTMIN  =    1.; // Minimum GEN pT
const double PTMAX  =  256.; // Maximum GEN pT
const double ETAMIN =   1.0; // Minimum GEN |eta|
const double ETAMAX =   2.5; // Maximum GEN |eta|

const double PTMAX_TRG = 128.;
const bool RPC_STUDY   = false;
const bool CLEAN_HI_PT = true;

const std::vector<int> MODES_CSC = {15};

// const std::vector<int> MODES_CSC = {7, 15};
// const std::vector<int> MODES_CSC = {11, 15};
const std::vector<int> MODES_RPC = {15};

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
   Use["BDTG_AWB_Sq"]             = 0;
   Use["BDTG_AWB_lite"]           = 0;

   Use["BDTG_AWB_50_trees"]       = 0;
   Use["BDTG_AWB_100_trees"]      = 0;
   Use["BDTG_AWB_200_trees"]      = 0;
   Use["BDTG_AWB_400_trees"]      = 0;
   Use["BDTG_AWB_800_trees"]      = 0;

   Use["BDTG_AWB_3_deep"]         = 0;
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
   std::vector<TString> mlist;
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      mlist = gTools().SplitString( myMethodList, ',' );
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
   TString out_dir = "/afs/cern.ch/work/a/abrinke1/public/EMTF/PtAssign2017";
   // out_dir = ".";
   TString out_file_name;
   out_file_name.Form( "%s/PtRegression_AWB_v1_17_03_14_mode_15_opt_bends_dTh_pt_1_256_clean.root", out_dir.Data() );
   // out_file_name.Form( "%s/PtRegression_AWB_v1_17_02_14_mode_7_RPC.root", out_dir.Data() );
   // out_file_name.Form( "%s/PtRegression_AWB_v1_17_02_14_mode_11_RPC.root", out_dir.Data() );
   // out_file_name.Form( "%s/PtRegression_AWB_v1_17_02_14_mode_19_RPC.root", out_dir.Data() );
   TFile* out_file = TFile::Open( out_file_name, "RECREATE" );

   // Read training and test data (see TMVAClassification for reading ASCII files)
   // load the signal and background event samples from ROOT trees
   TFile *input(0);
   TString store = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/ntuples";

   std::vector<TString> in_file_names;
   TString in_file_name;

   TString in_dirs[4] = { "ZeroBiasIsolatedBunch0/Slim/170130_224405/0000",
   			  "ZeroBiasIsolatedBunch1/Slim/170130_175144/0000",
   			  "ZeroBiasIsolatedBunch4/Slim/170130_175005/0000",
   			  "ZeroBiasIsolatedBunch5/Slim/170130_174947/0000" };
   for (int i = 0; i < 4; i++) {
     for (int j = 1; j < 50; j++) {
       in_file_name.Form("%s/%s/tuple_%d.root", store.Data(), in_dirs[i].Data(), j);
       std::cout << "Adding file " << in_file_name.Data() << std::endl;
       in_file_names.push_back(in_file_name.Data());
     }
   }

   // Load files with RPC hits
   TString in_dir_RPC = "SingleMu_Pt1To1000_FlatRandomOneOverPt/RPC/170213_173255/0000";
   Int_t nRPC_in = 0;
   for (int i = 1; i < 99; i++) {
     if (!RPC_STUDY) continue;
     in_file_name.Form("%s/%s/tuple_%d.root", store.Data(), in_dir_RPC.Data(), i);
     std::cout << "Adding file " << in_file_name.Data() << std::endl;
     in_file_names.push_back(in_file_name.Data());
     nRPC_in += 1;
     if (i*200000 > MAX_EVT) break; // ~100k events per file; half of events should have RPC hits
   }

   // Load files without RPC hits
   TString in_dir_CSC = "SingleMu_Pt1To1000_FlatRandomOneOverPt/EMTF_MuGun/170113_165434/0000";
   for (int i = 1; i < 99; i++) {
     in_file_name.Form("%s/%s/EMTF_MC_NTuple_SingleMu_noRPC_%d.root", store.Data(), in_dir_CSC.Data(), i);
     std::cout << "Adding file " << in_file_name.Data() << std::endl;
     in_file_names.push_back(in_file_name.Data());
     if ((i + nRPC_in)*100000 > MAX_EVT) break; // ~100k events per file
   }

   for (UInt_t i = 0; i < in_file_names.size(); i++) {
     if ( !gSystem->AccessPathName(in_file_names.at(i)) )
       input = TFile::Open( in_file_names.at(i) ); // check if file in local directory exists
     if (!input) {
       std::cout << "ERROR: could not open data file " << in_file_names.at(i) << std::endl;
       in_file_names.erase( in_file_names.begin()+i );
       i -= 1;
       // exit(1);
     }
   }

   // Add trees from the input files to the TChain
   // Have to use TChain for both SetBranchAddress and GetEntry to work
   std::vector<TChain*> in_chains;
   // Super-hacky ... but using "GetBranch" with a single chain with multiple files causes a segfault - AWB 19.01.16
   Int_t nChains_RPC = -99;
   for (UInt_t i = 0; i < in_file_names.size(); i++) {
     TChain *tmp_chain = new TChain("ntuple/tree");
     tmp_chain->Add( in_file_names.at(i) );
     in_chains.push_back(tmp_chain);
     if (i == nRPC_in - 1 && nChains_RPC < 0) 
       nChains_RPC = i + 1;
   }

   //////////////////////////////////////////////////////////////////////////
   ///  Factories: Use different sets of variables, target, weights, etc. ///
   //////////////////////////////////////////////////////////////////////////
   
   TString fact_set = "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression";
   std::vector<TString> var_names; // Holds names of variables for a given factory and permutation
   std::vector<Double_t> var_vals; // Holds values of variables for a given factory and permutation
   TMVA::Factory* nullF = new TMVA::Factory("NULL", out_file, fact_set); // Placeholder factory
   TMVA::DataLoader* nullL = new TMVA::DataLoader("NULL");                 // Placeholder loader

   // Tuple is defined by the factory and dataloader,  followed by a name, 
   // var name and value vectors, and hex bit masks for input and target variables.
   // Each hex bit represents four variables, e.g. 0x1 would select only the 1st variable, 
   // 0xf the 1st 4, 0xff the 1st 8, 0xa the 2nd and 4th, 0xf1 the 1st and 5th-8th, etc.
   std::vector< std::tuple<TMVA::Factory*, TMVA::DataLoader*, TString, std::vector<TString>, std::vector<Double_t>, int, int> > factories;

   // Optimized mode 15 - AWB 26.01.17
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01ff_0x4_invPt",    // dPhi12, 23, 34, theta, combs, FR1, St1 ring
   					 var_names, var_vals, 0x001f01ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01ff_0x4_invPtSq",  // dPhi12, 23, 34, theta, combs, FR1, St1 ring
   // 					 var_names, var_vals, 0x001f01ff, 0x4) );

   // Optimized mode 15 + dTh14 + bend1 - AWB 26.01.17
   factories.push_back( std::make_tuple( nullF, nullL, "f_0x401f11ff_0x4_invPt",    // dPhi12, 23, 34, theta, combs, FR1, St1 ring, dTh14, bend1
   					 var_names, var_vals, 0x401f11ff, 0x4) );

   // Optimized mode 15 + dTh14 + bend1 + bendMax3 - AWB 14.03.17
   factories.push_back( std::make_tuple( nullF, nullL, "f_0xc01f11ff_0x4_invPt",    // dPhi12, 23, 34, theta, combs, FR1, St1 ring, dTh14, bend1/max3
   					 var_names, var_vals, 0xc01f11ff, 0x4) );

   // // Optimized mode 19 - AWB 26.01.17
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0xf001f01ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, FR1, St1 ring, dPhi CSC-RPC
   // 					 var_names, var_vals, 0xf001f01ff, 0x4) );

   // // Original set of training variables
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x0000011d_0x2", 
   // 					 var_names, var_vals, 0x0000011d, 0x2) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x0000011d_0x4", 
   // 					 var_names, var_vals, 0x0000011d, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x0000011d_0x2_invPt", 
   // 					 var_names, var_vals, 0x0000011d, 0x2) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x0000011d_0x4_invPt", 
   // 					 var_names, var_vals, 0x0000011d, 0x4) );
   
   // // Original set of training variables + combinations
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01fd_0x2", 
   // 					 var_names, var_vals, 0x001f01fd, 0x2) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01fd_0x4", 
   // 					 var_names, var_vals, 0x001f01fd, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01fd_0x2_invPt", 
   // 					 var_names, var_vals, 0x001f01fd, 0x2) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01fd_0x4_invPt", 
   // 					 var_names, var_vals, 0x001f01fd, 0x4) );

   // // Original set of training variables + combinations + St 1 ring, FR bits, bending
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001fffff_0x2", 
   // 					 var_names, var_vals, 0x001fffff, 0x2) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001fffff_0x4", 
   // 					 var_names, var_vals, 0x001fffff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001fffff_0x2_invPt", 
   // 					 var_names, var_vals, 0x001fffff, 0x2) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001fffff_0x4_invPt", 
   // 					 var_names, var_vals, 0x001fffff, 0x4) );

   // // Various sets of variables
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x00000004_0x4_invPt",  // dPhi12
   // 					 var_names, var_vals, 0x00000004, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x00000005_0x4_invPt",  // dPhi12, theta
   // 					 var_names, var_vals, 0x00000005, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x0000000d_0x4_invPt",  // dPhi12, 23, theta 
   // 					 var_names, var_vals, 0x0000000d, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x00000085_0x4_invPt",  // dPhi12, 24, theta
   // 					 var_names, var_vals, 0x00000085, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x0000001d_0x4_invPt",  // dPhi12, 23, 34, theta 
   // 					 var_names, var_vals, 0x0000001d, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f00fd_0x4_invPt",  // dPhi12, 23, 34, theta, combs
   // 					 var_names, var_vals, 0x001f00fd, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f0ffd_0x4_invPt",  // dPhi12, 23, 34, theta, combs, FRs
   // 					 var_names, var_vals, 0x001f0ffd, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f0fff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, FRs, St1 ring
   // 					 var_names, var_vals, 0x001f0fff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001fffff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, FRs, St1 ring, bends
   // 					 var_names, var_vals, 0x001fffff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x8fff0fff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, FRs, St1 ring, dThetas
   // 					 var_names, var_vals, 0x8fff0fff, 0x4) );

   // // FRs
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f00ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring
   // 					 var_names, var_vals, 0x001f00ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, FR1
   // 					 var_names, var_vals, 0x001f01ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f03ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, FR1/2
   // 					 var_names, var_vals, 0x001f03ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f05ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, FR1/3
   // 					 var_names, var_vals, 0x001f05ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f09ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, FR1/4
   // 					 var_names, var_vals, 0x001f09ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f0fff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, all FRs
   // 					 var_names, var_vals, 0x001f0fff, 0x4) );

   // // Bends
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f01ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, FR1
   // 					 var_names, var_vals, 0x001f01ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f11ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, FR1, bend 1
   // 					 var_names, var_vals, 0x001f11ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f31ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, FR1, bend 1/2
   // 					 var_names, var_vals, 0x001f31ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f51ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, FR1, bend 1/3
   // 					 var_names, var_vals, 0x001f51ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001f91ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, FR1, bend 1/3
   // 					 var_names, var_vals, 0x001f91ff, 0x4) );
   // factories.push_back( std::make_tuple( nullF, nullL, "f_0x001ff1ff_0x4_invPt",  // dPhi12, 23, 34, theta, combs, St1 ring, FR1, all bends
   // 					 var_names, var_vals, 0x001ff1ff, 0x4) );


   // Initialize factories and dataloaders
   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
     std::get<0>(factories.at(iFact)) = new TMVA::Factory( std::get<2>(factories.at(iFact)), out_file, fact_set );
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
   
   in_vars.push_back( MVA_var( "theta",     "Track #theta",        "int", 'I', -88 ) ); // 0x0000 0001  * 2016 variable
   in_vars.push_back( MVA_var( "St1_ring",  "St 1 LCT ring",       "int", 'I', -88 ) ); // 0x0000 0002  
   in_vars.push_back( MVA_var( "dPhi_12",   "#phi(2) - #phi(1)",   "int", 'I', -88 ) ); // 0x0000 0004  * 2016 variable
   in_vars.push_back( MVA_var( "dPhi_23",   "#phi(3) - #phi(2)",   "int", 'I', -88 ) ); // 0x0000 0008  * 2016 variable

   in_vars.push_back( MVA_var( "dPhi_34",   "#phi(4) - #phi(3)",   "int", 'I', -88 ) ); // 0x0000 0010  * 2016 variable
   in_vars.push_back( MVA_var( "dPhi_13",   "#phi(3) - #phi(1)",   "int", 'I', -88 ) ); // 0x0000 0020  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "dPhi_14",   "#phi(4) - #phi(1)",   "int", 'I', -88 ) ); // 0x0000 0040  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "dPhi_24",   "#phi(4) - #phi(2)",   "int", 'I', -88 ) ); // 0x0000 0080  * Derivable from 2016 var

   in_vars.push_back( MVA_var( "FR_1",      "St 1 LCT F/R",        "int", 'I', -88 ) ); // 0x0000 0100  * 2016 variable
   in_vars.push_back( MVA_var( "FR_2",      "St 2 LCT F/R",        "int", 'I', -88 ) ); // 0x0000 0200
   in_vars.push_back( MVA_var( "FR_3",      "St 3 LCT F/R",        "int", 'I', -88 ) ); // 0x0000 0400
   in_vars.push_back( MVA_var( "FR_4",      "St 4 LCT F/R",        "int", 'I', -88 ) ); // 0x0000 0800

   in_vars.push_back( MVA_var( "bend_1",    "St 1 LCT bending",    "int", 'I', -88 ) ); // 0x0000 1000
   in_vars.push_back( MVA_var( "bend_2",    "St 2 LCT bending",    "int", 'I', -88 ) ); // 0x0000 2000
   in_vars.push_back( MVA_var( "bend_3",    "St 3 LCT bending",    "int", 'I', -88 ) ); // 0x0000 4000
   in_vars.push_back( MVA_var( "bend_4",    "St 4 LCT bending",    "int", 'I', -88 ) ); // 0x0000 8000

   in_vars.push_back( MVA_var( "dPhiSum4",  "#Sigmad#phi (6)",     "int", 'I', -88 ) ); // 0x0001 0000  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "dPhiSum4A", "#Sigma|d#phi| (6)",   "int", 'I', -88 ) ); // 0x0002 0000  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "dPhiSum3",  "#Sigmad#phi (3)",     "int", 'I', -88 ) ); // 0x0004 0000  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "dPhiSum3A", "#Sigma|d#phi| (3)",   "int", 'I', -88 ) ); // 0x0008 0000  * Derivable from 2016 var
   
   in_vars.push_back( MVA_var( "outStPhi",  "#phi outlier St",     "int", 'I', -88 ) ); // 0x0010 0000  * Derivable from 2016 var
   in_vars.push_back( MVA_var( "outStTh",   "#theta outlier St",   "int", 'I', -88 ) ); // 0x0020 0000
   in_vars.push_back( MVA_var( "dThMax4",   "Max d#theta (6)",     "int", 'I', -88 ) ); // 0x0040 0000
   in_vars.push_back( MVA_var( "dThMax3",   "Max d#theta (3)",     "int", 'I', -88 ) ); // 0x0080 0000

   in_vars.push_back( MVA_var( "dThSum4",   "#Sigmad#theta (6)",   "int", 'I', -88 ) ); // 0x0100 0000
   in_vars.push_back( MVA_var( "dThSum4A",  "#Sigma|d#theta| (6)", "int", 'I', -88 ) ); // 0x0200 0000
   in_vars.push_back( MVA_var( "dThSum3",   "#Sigmad#theta (3)",   "int", 'I', -88 ) ); // 0x0400 0000
   in_vars.push_back( MVA_var( "dThSum3A",  "#Sigma|d#theta| (3)", "int", 'I', -88 ) ); // 0x0800 0000

   in_vars.push_back( MVA_var( "dTh_12",  "#theta(2) - #theta(1)", "int", 'I', -88 ) ); // 0x1000 0000
   in_vars.push_back( MVA_var( "dTh_13",  "#theta(3) - #theta(1)", "int", 'I', -88 ) ); // 0x2000 0000
   in_vars.push_back( MVA_var( "dTh_14",  "#theta(4) - #theta(1)", "int", 'I', -88 ) ); // 0x4000 0000
   in_vars.push_back( MVA_var( "bendMax3", "Max bend, st 2-3-4",   "int", 'I', -88 ) ); // 0x8000 0000

   // in_vars.push_back( MVA_var( "dPhi_11", "CSC-RPC #Delta#phi(1)", "int", 'I', -88 ) ); // 0x1 0000 0000
   // in_vars.push_back( MVA_var( "dPhi_22", "CSC-RPC #Delta#phi(2)", "int", 'I', -88 ) ); // 0x2 0000 0000
   // in_vars.push_back( MVA_var( "dPhi_33", "CSC-RPC #Delta#phi(3)", "int", 'I', -88 ) ); // 0x4 0000 0000
   // in_vars.push_back( MVA_var( "dPhi_44", "CSC-RPC #Delta#phi(4)", "int", 'I', -88 ) ); // 0x8 0000 0000


   ////////////////////////////////////////////////////////////
   //  Target variable: true muon pT, or 1/pT, or log2(pT)  ///
   ////////////////////////////////////////////////////////////
   
   targ_vars.push_back( MVA_var( "GEN_pt",       "GEN p_{T}",               "GeV",      'F', -99 ) ); // 0x1
   targ_vars.push_back( MVA_var( "inv_GEN_pt",   "1 / GEN muon p_{T}",      "GeV^{-1}", 'F', -99 ) ); // 0x2
   targ_vars.push_back( MVA_var( "log2_GEN_pt",  "log_{2}(GEN muon p_{T})", "GeV",      'F', -99 ) ); // 0x4

   /////////////////////////////////////////////////////////////////////////////
   ///  Spectator variables: not used in training, but saved in output tree  ///
   /////////////////////////////////////////////////////////////////////////////
   
   spec_vars.push_back( MVA_var( "GEN_pt",       "GEN p_{T}",               "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_pt",      "EMTF p_{T}",              "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "inv_GEN_pt",   "1 / GEN muon p_{T}",      "GeV^{-1}", 'F', -77 ) );
   spec_vars.push_back( MVA_var( "inv_EMTF_pt",  "1 / EMTF p_{T}",          "GeV^{-1}", 'F', -77 ) );
   spec_vars.push_back( MVA_var( "log2_GEN_pt",  "log_{2}(GEN muon p_{T})", "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "log2_EMTF_pt", "log_{2}(EMTF p_{T})",     "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "GEN_eta",      "GEN #eta",                "",         'F', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_eta",     "EMTF #eta",               "",         'F', -77 ) );
   spec_vars.push_back( MVA_var( "GEN_charge",   "GEN charge",              "",         'I', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_charge",  "EMTF charge",             "",         'I', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_mode",    "EMTF mode",               "",         'I', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_hasRPC",  "EMTF has RPC hit",        "",         'I', -77 ) );
   spec_vars.push_back( MVA_var( "dPhi_12_sign", "#phi(2) - #phi(1) sign",  "",         'I', -77 ) );

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
   UInt_t iEvtZB = 0;
   UInt_t nTrain = 0;
   UInt_t nTest  = 0;

   UInt_t nAlready = 0;
   UInt_t nFound   = 0;
   UInt_t nMissing = 0;

   for (int iCh = 0; iCh < in_chains.size(); iCh++) {
     TChain *in_chain = in_chains.at(iCh);
     
     // Get branches from the chain
     TBranch *muon_br  = in_chain->GetBranch("muon");
     TBranch *hit_br   = in_chain->GetBranch("hit");
     TBranch *track_br = in_chain->GetBranch("track");
     
     std::cout << "\n******* About to enter the event loop for chain " << iCh+1 << " *******" << std::endl;
     
     for (UInt_t jEvt = 0; jEvt < in_chain->GetEntries(); jEvt++) {
       if (iEvt > MAX_EVT) break;
       // if ( (iEvt % REPORT_EVT) == 0 )
       // 	 std::cout << "Looking at event " << iEvt << std::endl;
       in_chain->GetEntry(jEvt);
       
       UInt_t nMuons   = (muon_br->GetLeaf("nMuons"))->GetValue();
       UInt_t nTracks  = (track_br->GetLeaf("nTracks"))->GetValue();
       UInt_t nHits    = (track_br->GetLeaf("nHits"))->GetValue();
       Bool_t isMC     = (nMuons > 0);
       Bool_t trainEvt = true;  // Can use the event for training
       if (not isMC)  // Process ZeroBias anyway
	 nMuons = nTracks;
       // std::cout << "There are " << nMuons << " GEN muons and " << nTracks << " EMTF tracks\n" << std::endl;
       
       if ( ( (iEvt % REPORT_EVT) == 0 && isMC) || (iEvtZB > 0 && (iEvtZB % 100000) == 0) )
	 std::cout << "Looking at event " << iEvt << " (" << iEvtZB << ")" << std::endl;
       
       for (UInt_t iMu = 0; iMu < nMuons; iMu++) {
	 Double_t mu_pt ;
	 Double_t mu_eta;
	 Double_t mu_phi;
	 Int_t mu_charge;
	 if (isMC) {
	   mu_pt     = (muon_br->GetLeaf("pt"))->GetValue(iMu);
	   mu_eta    = (muon_br->GetLeaf("eta"))->GetValue(iMu);
	   mu_phi    = (muon_br->GetLeaf("phi"))->GetValue(iMu);
	   mu_charge = (muon_br->GetLeaf("charge"))->GetValue(iMu);
	 } else {
	   mu_pt     = 999.;
	   mu_eta    = -99.;
	   mu_phi    = -99.;
	   mu_charge = -99;
	 }

	 if ( mu_pt < PTMIN && isMC ) trainEvt = false;
	 if ( mu_pt > PTMAX && isMC ) trainEvt = false;
	 if ( fabs( mu_eta ) < ETAMIN && isMC ) continue;
	 if ( fabs( mu_eta ) > ETAMAX && isMC ) continue;
	 // std::cout << "Muon " << iMu+1 << " has pt = " << mu_pt << ", eta = " << mu_eta << ", phi = " << mu_phi << std::endl;
	 
	 for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++) {
	   Double_t trk_pt    = (track_br->GetLeaf("pt"))->GetValue(iTrk);
	   Double_t trk_eta   = (track_br->GetLeaf("eta"))->GetValue(iTrk);
	   trk_pt *= PT_SCALE;
	   Int_t trk_theta_int = (track_br->GetLeaf("theta_int"))->GetValue(iTrk);
	   Double_t trk_phi    = (track_br->GetLeaf("phi"))->GetValue(iTrk);
	   Int_t trk_phi_int   = (track_br->GetLeaf("phi_int"))->GetValue(iTrk);
	   Int_t trk_charge    = (track_br->GetLeaf("charge"))->GetValue(iTrk);
	   Int_t trk_mode      = (track_br->GetLeaf("mode"))->GetValue(iTrk);

	   if ( ( mu_eta > 0 ) != ( trk_eta > 0 ) && isMC ) continue;
	   if (trk_theta_int == 0) continue;  // Not sure what causes this buggy condition - AWB 14.02.17
	   Bool_t goodMode = false;
	   // std::cout << "iCh = " << iCh << ", RPC_STUDY = " << RPC_STUDY << ", trk_mode = " << trk_mode << std::endl;
	   for (UInt_t iMode = 0; iMode < MODES_CSC.size(); iMode++)
	     if ( (iCh >= nChains_RPC || !RPC_STUDY) && trk_mode == MODES_CSC.at(iMode))
	       goodMode = true;
	   for (UInt_t iMode = 0; iMode < MODES_RPC.size(); iMode++)
	     if ( (iCh <  nChains_RPC &&  RPC_STUDY) && trk_mode == MODES_RPC.at(iMode))
	       goodMode = true;
	   if (not goodMode) continue;
	   // std::cout << "  * Track " << iTrk+1 << " has pt = " << trk_pt << ", eta = " << trk_eta 
	   // 	   << ", phi = " << mu_phi << ", mode = " << trk_mode << std::endl;

	   // Specifically looking at station 1 or 2 for now
	   Int_t has_RPC = (track_br->GetLeaf("hit_isRPC"))->GetValue(4*iTrk + 0);
	   // Int_t has_RPC = (track_br->GetLeaf("hit_isRPC"))->GetValue(4*iTrk + 1);
	   Int_t hit_theta_int = -9999;
	   Int_t hit_phi_int = -9999;
	   if ( iCh < nChains_RPC && RPC_STUDY ) {
	     if (has_RPC == 1) {
	       hit_theta_int = (track_br->GetLeaf("hit_theta_int"))->GetValue(4*iTrk + 0);
	       hit_phi_int   = (track_br->GetLeaf("hit_phi_int"))->GetValue(4*iTrk + 0);
	       // hit_theta_int = (track_br->GetLeaf("hit_theta_int"))->GetValue(4*iTrk + 1);
	       // hit_phi_int   = (track_br->GetLeaf("hit_phi_int"))->GetValue(4*iTrk + 1);
	       nAlready += 1;
	     } else {
	       for (UInt_t iHit = 0; iHit < nHits; iHit++) {
		 if ( (hit_br->GetLeaf("isRPC"))->GetValue(iHit) == 1 &&
		      (hit_br->GetLeaf("station"))->GetValue(iHit) == 1 &&
		      // (hit_br->GetLeaf("station"))->GetValue(iHit) == 2 &&
		      (hit_br->GetLeaf("sector_index"))->GetValue(iHit) == (track_br->GetLeaf("sector_index"))->GetValue(iTrk) && 
		      abs( (hit_br->GetLeaf("theta_int"))->GetValue(iHit) - trk_theta_int ) < 8 &&
		      abs( (hit_br->GetLeaf("phi_int"))->GetValue(iHit) - trk_phi_int ) < 512 &&
		      abs( (hit_br->GetLeaf("phi_int"))->GetValue(iHit) - trk_phi_int ) < abs(hit_phi_int - trk_phi_int) ) {
		   hit_theta_int = (hit_br->GetLeaf("theta_int"))->GetValue(iHit);
		   hit_phi_int   = (hit_br->GetLeaf("phi_int"))->GetValue(iHit);
		   has_RPC = 1;
		 }
	       }
	     }
	     if (has_RPC != 1) {
	       nMissing += 1;
	       continue;
	     } else {
	       nFound += 1;
	     }
	   } else if (has_RPC != 1) {
	     has_RPC = 0;
	   }


	   Int_t st1_ring = ((int) (track_br->GetLeaf("hit_ring"))->GetValue(4*iTrk + 0)) % 3; // Ring 4 --> Ring 1 (ME1/1a)

	   Int_t phi1 = (track_br->GetLeaf("hit_phi_int"))->GetValue(4*iTrk + 0);
	   Int_t phi2 = (track_br->GetLeaf("hit_phi_int"))->GetValue(4*iTrk + 1);
	   Int_t phi3 = (track_br->GetLeaf("hit_phi_int"))->GetValue(4*iTrk + 2);
	   Int_t phi4 = (track_br->GetLeaf("hit_phi_int"))->GetValue(4*iTrk + 3);
	   // std::cout << "  * Hits have phi = " << phi1 << ", " << phi2 << ", " << phi3 << ", " << phi4 << std::endl;
	   
	   Int_t th1 = (track_br->GetLeaf("hit_theta_int"))->GetValue(4*iTrk + 0);
	   Int_t th2 = (track_br->GetLeaf("hit_theta_int"))->GetValue(4*iTrk + 1);
	   Int_t th3 = (track_br->GetLeaf("hit_theta_int"))->GetValue(4*iTrk + 2);
	   Int_t th4 = (track_br->GetLeaf("hit_theta_int"))->GetValue(4*iTrk + 3);
	   // std::cout << "  * Hits have theta = " << th1 << ", " << th2 << ", " << th3 << ", " << th4 << std::endl;

	   if ( iCh < nChains_RPC && RPC_STUDY ) phi1 = hit_phi_int;
	   if ( iCh < nChains_RPC && RPC_STUDY ) th1 = hit_theta_int;
	   // if ( iCh < nChains_RPC && RPC_STUDY ) phi2 = hit_phi_int;
	   // if ( iCh < nChains_RPC && RPC_STUDY ) th2 = hit_theta_int;
	   
	   Int_t dPhi12 = phi2 - phi1;
	   Int_t dPhi13 = phi3 - phi1;
	   Int_t dPhi14 = phi4 - phi1;
	   Int_t dPhi23 = phi3 - phi2;
	   Int_t dPhi24 = phi4 - phi2;
	   Int_t dPhi34 = phi4 - phi3;

	   Int_t dTh12 = th2 - th1;
	   Int_t dTh13 = th3 - th1;
	   Int_t dTh14 = th4 - th1;
	   Int_t dTh23 = th3 - th2;
	   Int_t dTh24 = th4 - th2;
	   Int_t dTh34 = th4 - th3;

	   if (trk_mode == 7 && has_RPC != 1) {
	     dPhi12 = 0;
	     dPhi13 = dPhi23;
	     dPhi14 = dPhi34;
	   
	     dTh12 = 0;
	     dTh13 = dTh23;
	     dTh14 = dTh24;
	   }
	   // if (trk_mode == 11 && has_RPC != 1) {
	   //   dPhi12 = dPhi13 / 2;
	   //   dPhi23 = dPhi13 - dPhi12;
	   //   dPhi24 = dPhi23 + dPhi34;
	   
	   //   dTh12 = dTh13 / 2;
	   //   dTh23 = dTh13 - dTh12;
	   //   dTh24 = dTh23 + dTh34;
	   // }
	   
	   
	   // Define all dPhi values relative to dPhi12
	   Int_t dPhi12_sign = ( (dPhi12 < 0) ? -1. : 1. );
	   dPhi13 *= dPhi12_sign;
	   dPhi14 *= dPhi12_sign;
	   dPhi23 *= dPhi12_sign;
	   dPhi24 *= dPhi12_sign;
	   dPhi34 *= dPhi12_sign;
	   dPhi12  = abs(dPhi12);
	   
	   Int_t dPhi_sum_4  = dPhi12 + dPhi13 + dPhi14 + dPhi23 + dPhi24 + dPhi34;
	   Int_t dPhi_sum_4A = abs(dPhi12) + abs(dPhi13) + abs(dPhi14) + abs(dPhi23) + abs(dPhi24) + abs(dPhi34);
	   Int_t pDev_st1 = abs(dPhi12) + abs(dPhi13) + abs(dPhi14);
	   Int_t pDev_st2 = abs(dPhi12) + abs(dPhi23) + abs(dPhi24);
	   Int_t pDev_st3 = abs(dPhi13) + abs(dPhi23) + abs(dPhi34);
	   Int_t pDev_st4 = abs(dPhi14) + abs(dPhi24) + abs(dPhi34);
	   
	   Int_t dTh_sum_4  = dTh12 + dTh13 + dTh14 + dTh23 + dTh24 + dTh34;
	   Int_t dTh_sum_4A = abs(dTh12) + abs(dTh13) + abs(dTh14) + abs(dTh23) + abs(dTh24) + abs(dTh34);
	   Int_t tDev_st1 = abs(dTh12) + abs(dTh13) + abs(dTh14);
	   Int_t tDev_st2 = abs(dTh12) + abs(dTh23) + abs(dTh24);
	   Int_t tDev_st3 = abs(dTh13) + abs(dTh23) + abs(dTh34);
	   Int_t tDev_st4 = abs(dTh14) + abs(dTh24) + abs(dTh34);
	   
	   Int_t out_st_phi = -88;
	   if      (pDev_st4 > pDev_st3 && pDev_st4 > pDev_st2 && pDev_st4 > pDev_st1)  out_st_phi = 4;
	   else if (pDev_st3 > pDev_st4 && pDev_st3 > pDev_st2 && pDev_st3 > pDev_st1)  out_st_phi = 3;
	   else if (pDev_st2 > pDev_st4 && pDev_st2 > pDev_st3 && pDev_st2 > pDev_st1)  out_st_phi = 2;
	   else if (pDev_st1 > pDev_st4 && pDev_st1 > pDev_st3 && pDev_st1 > pDev_st2)  out_st_phi = 1;
	   else                                                                         out_st_phi = 0;
	   
	   Int_t out_st_th = -88;
	   if      (tDev_st4 > tDev_st3 && tDev_st4 > tDev_st2 && tDev_st4 > tDev_st1)  out_st_th = 4;
	   else if (tDev_st3 > tDev_st4 && tDev_st3 > tDev_st2 && tDev_st3 > tDev_st1)  out_st_th = 3;
	   else if (tDev_st2 > tDev_st4 && tDev_st2 > tDev_st3 && tDev_st2 > tDev_st1)  out_st_th = 2;
	   else if (tDev_st1 > tDev_st4 && tDev_st1 > tDev_st3 && tDev_st1 > tDev_st2)  out_st_th = 1;
	   else                                                                         out_st_th = 0;
	   
 	   Int_t dPhi_sum_3  = -88;
	   Int_t dPhi_sum_3A = -88;
	   if      (out_st_phi == 4) {
	     dPhi_sum_3  = dPhi12 + dPhi13 + dPhi23;
	     dPhi_sum_3A = abs(dPhi12) + abs(dPhi13) + abs(dPhi23);
	   } else if (out_st_phi == 3) {
	     dPhi_sum_3  = dPhi12 + dPhi14 + dPhi24;
	     dPhi_sum_3A = abs(dPhi12) + abs(dPhi14) + abs(dPhi24);
	   } else if (out_st_phi == 2) {
	     dPhi_sum_3  = dPhi13 + dPhi14 + dPhi34;
	     dPhi_sum_3A = abs(dPhi13) + abs(dPhi14) + abs(dPhi34);
	   } else {
	     dPhi_sum_3  = dPhi23 + dPhi24 + dPhi34;
	     dPhi_sum_3A = abs(dPhi23) + abs(dPhi24) + abs(dPhi34);
	   }

	   // Clean out showering muons with outlier station 1, or >= 2 outlier stations
	   if (isMC && log2(mu_pt) > 6 && CLEAN_HI_PT)
	     if ( dPhi_sum_4A >= std::max(40., 332. - 40*log2(mu_pt)) )
	       if ( out_st_phi < 2 || dPhi_sum_3A >= std::max(24., 174. - 20*log2(mu_pt)) )
		 trainEvt = false;

	   Int_t dTh_max_4 = MaxOfSix( dTh12, dTh13, dTh14, dTh23, dTh24, dTh34 );
	   Int_t dTh_max_3  = -88;
 	   Int_t dTh_sum_3  = -88;
	   Int_t dTh_sum_3A = -88;
	   if      (out_st_th == 4) {
	     dTh_sum_3  = dTh12 + dTh13 + dTh23;
	     dTh_sum_3A = abs(dTh12) + abs(dTh13) + abs(dTh23);
	     dTh_max_3  = MaxOfThree(dTh12, dTh13, dTh23);
	   } else if (out_st_th == 3) {
	     dTh_sum_3  = dTh12 + dTh14 + dTh24;
	     dTh_sum_3A = abs(dTh12) + abs(dTh14) + abs(dTh24);
	     dTh_max_3	 = MaxOfThree(dTh12, dTh14, dTh24);
	   } else if (out_st_th == 2) {
	     dTh_sum_3  = dTh13 + dTh14 + dTh34;
	     dTh_sum_3A = abs(dTh13) + abs(dTh14) + abs(dTh34);
	     dTh_max_3	 = MaxOfThree(dTh13, dTh14, dTh34);
	   } else {
	     dTh_sum_3  = dTh23 + dTh24 + dTh34;
	     dTh_sum_3A = abs(dTh23) + abs(dTh24) + abs(dTh34);
	     dTh_max_3	 = MaxOfThree(dTh23, dTh24, dTh34);
	   }
	   
	   Int_t FR1 = (track_br->GetLeaf("hit_FR"))->GetValue(4*iTrk + 0);
	   Int_t FR2 = (track_br->GetLeaf("hit_FR"))->GetValue(4*iTrk + 1);
	   Int_t FR3 = (track_br->GetLeaf("hit_FR"))->GetValue(4*iTrk + 2);
	   Int_t FR4 = (track_br->GetLeaf("hit_FR"))->GetValue(4*iTrk + 3);
	   
	   Int_t bend1 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 0), trk_eta ) * dPhi12_sign;
	   Int_t bend2 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 1), trk_eta ) * dPhi12_sign;
	   Int_t bend3 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 2), trk_eta ) * dPhi12_sign;
	   Int_t bend4 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 3), trk_eta ) * dPhi12_sign;

	   Int_t bendMax_3 = 0;
	   if (abs(bend2) > abs(bendMax_3)) bendMax_3 = bend2;
	   if (abs(bend3) > abs(bendMax_3)) bendMax_3 = bend3;
	   if (abs(bend4) > abs(bendMax_3)) bendMax_3 = bend4;
	   
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
		 var_vals.at(iVar) = abs(trk_theta_int);
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
	       if ( vName == "outStTh" )
		 var_vals.at(iVar) = out_st_th;
	       if ( vName == "dThMax4" )
		 var_vals.at(iVar) = dTh_max_4;
	       if ( vName == "dThMax3" )
		 var_vals.at(iVar) = dTh_max_3;

	       if ( vName == "dThSum4" )
		 var_vals.at(iVar) = dTh_sum_4;
	       if ( vName == "dThSum4A" )
		 var_vals.at(iVar) = dTh_sum_4A;
	       if ( vName == "dThSum3" )
		 var_vals.at(iVar) = dTh_sum_3;
	       if ( vName == "dThSum3A" )
		 var_vals.at(iVar) = dTh_sum_3A;

	       if ( vName == "dTh_12" )
		 var_vals.at(iVar) = dTh12;
	       if ( vName == "dTh_13" )
		 var_vals.at(iVar) = dTh13;
	       if ( vName == "dTh_14" )
		 var_vals.at(iVar) = dTh14;
	       if ( vName == "bendMax3" )
		 var_vals.at(iVar) = bendMax_3;

	       // if ( vName == "dPhi_11" )
	       // 	 var_vals.at(iVar) = dPhi11;
	       // if ( vName == "dPhi_22" )
	       // 	 var_vals.at(iVar) = dPhi22;
	       // if ( vName == "dPhi_33" )
	       // 	 var_vals.at(iVar) = dPhi33;
	       // if ( vName == "dPhi_44" )
	       // 	 var_vals.at(iVar) = dPhi44;

	       
	       ////////////////////////////////////////
	       ///  Target and spectator variables  ///
	       ////////////////////////////////////////
	       
	       if ( vName == "GEN_pt" )
		 var_vals.at(iVar) = min(mu_pt, PTMAX_TRG);
	       if ( vName == "EMTF_pt" )
		 var_vals.at(iVar) = trk_pt;
	       if ( vName == "inv_GEN_pt" )
		 var_vals.at(iVar) = 1. / min(mu_pt, PTMAX_TRG);
	       if ( vName == "inv_EMTF_pt" )
		 var_vals.at(iVar) = 1. / trk_pt;
	       if ( vName == "log2_GEN_pt" )
		 var_vals.at(iVar) = log2(min(mu_pt, PTMAX_TRG));
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
	       if ( vName == "EMTF_mode" )
		 var_vals.at(iVar) = trk_mode;
	       if ( vName == "EMTF_hasRPC" )
		 var_vals.at(iVar) = has_RPC;
	       
	     } // End loop: for (UInt_t iVar = 0; iVar < var_names.size(); iVar++)
	     
	     // Unweighted distribution: flat in eta and 1/pT
	     Double_t evt_weight = 1.0;
	     // Weight by 1/pT so overall distribution is (1/pT)^2
	     if ( std::get<2>(factories.at(iFact)).Contains("_invPt") )
	       evt_weight = 1. / mu_pt;
	     if ( std::get<2>(factories.at(iFact)).Contains("_invPtSq") )
	       evt_weight = 1. / pow(mu_pt, 2);
	     
	     // Load values into event
	     if ( (iEvt % 2) == 0 && isMC && trainEvt && nTrain < (MAX_TR - (iFact == 0)) && (!RPC_STUDY || iCh < nChains_RPC) ) {
	       std::get<1>(factories.at(iFact))->AddTrainingEvent( "Regression", var_vals, evt_weight );
	       if (iFact == 0) nTrain += 1;
	       // std::cout << "Added train event " << nTrain << std::endl;
	     }
	     else {
	       std::get<1>(factories.at(iFact))->AddTestEvent( "Regression", var_vals, evt_weight );
	       if (iFact == 0) nTest += 1;
	       // std::cout << "Added test event " << nTest << std::endl;
	     }
	     
	   } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++) 
	   
	 } // End loop: for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++)
       } // End loop: for (UInt_t iMu = 0; iMu < nMuons; iMu++)
       if (isMC) iEvt += 1;
       else iEvtZB += 1;
     } // End loop: for (UInt_t jEvt = 0; jEvt < in_chain->GetEntries(); jEvt++)
   } // End loop: for (int iCh = 0; iCh < in_chains.size(); iCh++) {

   std::cout << "******* Made it out of the event loop *******" << std::endl;

   std::cout << "nAlready = " << nAlready << ", nFound = " << nFound << ", nMissing = " << nMissing << std::endl;
   
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
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB", (string)
			  "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.000001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     // AWB settings - LeastSquares
     if (Use["BDTG_AWB_Sq"]) // Optimized settings
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_Sq", (string)
			  "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.000001:"+
			  "RegressionLossFunctionBDTG=LeastSquares" );

     if (Use["BDTG_AWB_lite"]) // Fast, simple BDT
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_lite", (string)
			  "!H:!V:NTrees=40::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.01:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );

     if (Use["BDTG_AWB_50_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_50_trees", (string)
			  "!H:!V:NTrees=50::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_100_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_100_trees", (string)
			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_200_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_200_trees", (string)
			  "!H:!V:NTrees=200::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_400_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_400_trees", (string)
			  "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_800_trees"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_800_trees", (string)
			  "!H:!V:NTrees=800::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=3:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     
     if (Use["BDTG_AWB_3_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_3_deep", (string)
			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=4:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_4_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_4_deep", (string)
			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=4:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_5_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_5_deep", (string)
			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.001:"+
			  "RegressionLossFunctionBDTG=AbsoluteDeviation" );
     if (Use["BDTG_AWB_6_deep"])
       factX->BookMethod( loadX, TMVA::Types::kBDT, "BDTG_AWB_6_deep", (string)
			  "!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=6:MinNodeSize=0.001:"+
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
     
     // // Evaluate and compare performance of all configured MVAs
     // factX->EvaluateAllMethods();

     // Instead of "EvaluateAllMethods()", just write out the training and testing trees
     // Skip unnecessary evaluatioh histograms, which take time on large datasets 
     // Code gleaned from original "EvaluateAllMethods()" function in tmva/tmva/src/Factory.cxx - AWB 31.01.17
     if ( factX->fMethodsMap.empty() )
       std::cout << "factX->fMethodsMap is empty" << std::endl;
     
     std::map<TString, std::vector<IMethod*>*>::iterator itrMap;
     for (itrMap = factX->fMethodsMap.begin(); itrMap != factX->fMethodsMap.end(); itrMap++) {
       
       std::vector<IMethod*> *methods = itrMap->second;
       std::list<TString> datasets;
       Int_t nmeth_used[2] = {int(mlist.size()), 1};
       
       for (Int_t k = 0; k < 2; k++) {
     	 for (Int_t i = 0; i < nmeth_used[k]; i++) {
     	   MethodBase* theMethod = dynamic_cast<MethodBase*>((*methods)[i]);
     	   if (theMethod == 0) {
     	     std::cout << "For k = " << k << ", i = " << i << ", no valid method" << std::endl;
     	     continue;
     	   }
     	   TDirectory* RootBaseDir = (TDirectory*) out_file;
     	   RootBaseDir->cd( std::get<2>(factories.at(iFact)) );
     	   if ( std::find( datasets.begin(), datasets.end(), std::get<2>(factories.at(iFact)) ) == datasets.end() ) {
     	     theMethod->Data()->GetTree(Types::kTesting)->Write( "", TObject::kOverwrite );
     	     theMethod->Data()->GetTree(Types::kTraining)->Write( "", TObject::kOverwrite );
     	     datasets.push_back( std::get<2>(factories.at(iFact)) );
     	   }
     	 } // End loop: for (Int_t i = 0; i < nmeth_used[k]; i++)
       } // End loop: for (Int_t k = 0; k < 2; k++)
     } // End loop: for (itrMap = factX->fMethodsMap.begin(); itrMap != factX->fMethodsMap.end(); itrMap++) 

     // --------------------------------------------------------------
     
   } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++)
   
   // Save the output
   out_file->Close();

   std::cout << "==> Wrote root file: " << out_file->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;

   // delete factory;
   // delete dataloader;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVARegGui( out_file_name );
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

