/////////////////////////////////////////////////////////////////////
///  Simulataneous pT regression for multiple factories and MVAs  ///
///                  Andrew Brinkerhoff 23.01.17                  ///
///                                                               ///
///  Adapted from ROOT TMVARegression.C                           ///
///  Run using "root -l PtRegression_AWB_v2.C                     /// 
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
#include "src/TrackBuilder.cc"
#include "src/PtLutVarCalc.cc"


////////////////////////////////
///  Configuration settings  ///
////////////////////////////////

// *** Number of events *** //
// // Typical settings with ZeroBias, without RPC
// const int MAX_EVT    = 5000000; // Maximum number of MuGun events to process
// const int MAX_TR     =  500000; // Maximum number of MuGun training events
// const int REPORT_EVT =  100000;

// // Typical settings with ZeroBias, with RPC
// const int MAX_EVT    = 10000000; // Maximum number of MuGun events to process
// const int MAX_TR     =  5000000; // Maximum number of MuGun training events
// const int REPORT_EVT =   100000;

// // // Typical settings without ZeroBias, with RPC
// const int MAX_EVT    = 30000000; // Maximum number of MuGun events to process
// const int MAX_TR     = 20000000; // Maximum number of MuGun training events
// const int REPORT_EVT =   100000;

// Typical test settings
const int MAX_EVT    = 2000;
const int MAX_TR     = 1000;
const int REPORT_EVT =  100;

// *** Constants *** //
const double PI  = 3.14159265359;
const double BIT = 0.000001; // Tiny value or offset

// *** Muon kinematics *** //
const double PTMIN  =      1.; // Minimum GEN pT
const double PTMAX  =   1000.; // Maximum GEN pT
const double ETAMIN =     1.0; // Minimum GEN |eta|
const double ETAMAX =     2.5; // Maximum GEN |eta|

// *** Track-building settings *** //
const int  MODE     = 13;    // Track mode to build
const bool USE_RPC  = false; // Use input files with RPC hits 
const int  MAX_RPC  = 0;     // Maximum number of RPC hits per track
const int  MIN_CSC  = 2;     // Minimum number of CSC LCTs per track
const int  MAX_DPH  = 1024;  // Maximum dPhi between hits for track-building (excludes maximum)
const int  MAX_DTH  = 8;     // Maximum dTheta between hits for track-building (includes maximum)
const bool BIT_COMP = false; // Use bit-compressed versions of input variables
const std::vector<int> CSC_MASK = {}; // Mask CSC LCTs in these stations  
const std::vector<int> RPC_MASK = {}; // Mask RPC hits in these stations  

// *** High-pT muons *** //
const double PTMAX_TR =  256.;  // Maximum GEN pT for training
const double PTMAX_TRG = 128.;  // Maximum trigger pT assigned
const bool CLEAN_HI_PT = false; // Remove showering high-pT mode 15 tracks from training

// *** EMTF tracks *** //
const bool REQ_EMTF = false; // Require that an EMTF muon be matched to the GEN muon 
const std::vector<int> EMTF_MODES = {15, 14, 13, 12, 11, 9}; // Modes of saved EMTF tracks

// *** Output data options *** //
const bool SPEC_VARS = true;


//////////////////////////////////
///  Main executable function  ///
//////////////////////////////////

using namespace TMVA;

void PtRegression_AWB_v2 ( TString myMethodList = "" ) {

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

   Use["BDTG_AWB"]                = 0;
   Use["BDTG_AWB_Sq"]             = 0;
   Use["BDTG_AWB_lite"]           = 1;

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
   std::cout << "==> Start PtRegression_AWB_v2" << std::endl;

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
   out_dir = ".";
   TString out_file_name;
   out_file_name.Form( "%s/PtRegression_AWB_v2_17_04_07_test.root", out_dir.Data() );
   TFile* out_file = TFile::Open( out_file_name, "RECREATE" );

   // Read training and test data (see TMVAClassification for reading ASCII files)
   // load the signal and background event samples from ROOT trees
   TFile *input(0);
   TString store = "root://eoscms.cern.ch//store/user/abrinke1/EMTF/Emulator/ntuples";

   std::vector<TString> in_file_names;
   TString in_file_name;

   int nZB_in = 0;
   // if (USE_RPC) {
   //   TString in_dir = "ZeroBiasIsolatedBunch0/Slim_RPC/170213_174254/0000";
   //   for (int j = 1; j < 50; j++) {
   //     in_file_name.Form("%s/%s/tuple_%d.root", store.Data(), in_dir.Data(), j);
   //     std::cout << "Adding file " << in_file_name.Data() << std::endl;
   //     in_file_names.push_back(in_file_name.Data());
   //     nZB_in += 1;
   //   }
   // } else {
   //   TString in_dirs[4] = { "ZeroBiasIsolatedBunch0/Slim/170130_224405/0000",
   // 			    "ZeroBiasIsolatedBunch1/Slim/170130_175144/0000",
   // 			    "ZeroBiasIsolatedBunch4/Slim/170130_175005/0000",
   // 			    "ZeroBiasIsolatedBunch5/Slim/170130_174947/0000" };
     
   //   for (int i = 0; i < 4; i++) {
   //     for (int j = 1; j < 50; j++) {
   // 	 in_file_name.Form("%s/%s/tuple_%d.root", store.Data(), in_dirs[i].Data(), j);
   // 	 std::cout << "Adding file " << in_file_name.Data() << std::endl;
   // 	 in_file_names.push_back(in_file_name.Data());
   // 	 nZB_in += 1;
   //     }
   //   }
   // }


   // Load files with RPC hits
   TString in_dir_RPC = "SingleMu_Pt1To1000_FlatRandomOneOverPt/RPC/170213_173255/0000";
   int nRPC_in = 0;
   for (int i = 1; i < 99; i++) {
     if (!USE_RPC) continue;
     in_file_name.Form("%s/%s/tuple_%d.root", store.Data(), in_dir_RPC.Data(), i);
     std::cout << "Adding file " << in_file_name.Data() << std::endl;
     in_file_names.push_back(in_file_name.Data());
     nRPC_in += 1;
     if (i*200000 > MAX_EVT) break; // ~100k events per file; half of events should have RPC hits
   }
   
   // Load files without RPC hits
   TString in_dir_CSC = "SingleMu_Pt1To1000_FlatRandomOneOverPt/EMTF_MuGun/170113_165434/0000";
   for (int i = 1; i < 99; i++) {
     if (USE_RPC) continue;
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
       if (i < nZB_in) 
	 nZB_in -= 1;
       else if (i < (nZB_in + nRPC_in)) 
	 nRPC_in -= 1;
       i -= 1;
       // exit(1);
     }
   }

   // Add trees from the input files to the TChain
   // Have to use TChain for both SetBranchAddress and GetEntry to work
   std::vector<TChain*> in_chains;
   // Super-hacky ... but using "GetBranch" with a single chain with multiple files causes a segfault - AWB 19.01.16
   int nChains_ZB  = -99;
   int nChains_RPC = -99;
   for (UInt_t i = 0; i < in_file_names.size(); i++) {
     TChain *tmp_chain = new TChain("ntuple/tree");
     tmp_chain->Add( in_file_names.at(i) );
     in_chains.push_back(tmp_chain);
     if ( i == (nZB_in - 1) && nChains_ZB < 0)
       nChains_ZB = i + 1;
     if ( i == (nZB_in + nRPC_in - 1) && nChains_RPC < 0) 
       nChains_RPC = i + 1 - nChains_ZB;
   }
   // assert(nChains_ZB > 0 && nChains_RPC > 0);
   // assert(nChains_RPC > 0);

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

   if        (MODE == 15) {
     // BASELINE mode 15 - dPhi12/23/34 + combos, theta, FR1, St1 ring, dTh14, bend1, RPC 1/2/3/4
     factories.push_back( std::make_tuple( nullF, nullL, "f_0x041f11ff_0x4_invPt",    
					   var_names, var_vals, 0x041f11ff, 0x4) );
   } else if (MODE == 13) {
     // BASELINE mode 13 - dPhi12, 24, theta, FR1/2, St1 ring, dTh14, bend1
     factories.push_back( std::make_tuple( nullF, nullL, "f_0x100013c7_0x4_invPt",
					   var_names, var_vals, 0x100013c7, 0x4) );
   }

     // // BASELINE mode 7 - AWB 18.03.17
     // factories.push_back( std::make_tuple( nullF, nullL, "f_0x40002299_0x4_invPt",    // dPhi23, 34, theta, FR2, dTh24, bend2
     // 					 var_names, var_vals, 0x40002299, 0x4) );
     
     // // BASELINE mode 11 - AWB 18.03.17
     // factories.push_back( std::make_tuple( nullF, nullL, "f_0x10001573_0x4_invPt",    // dPhi13, 34, theta, FR1/3, St1 ring, dTh14, bend1
     // 					 var_names, var_vals, 0x10001573, 0x4) );
     
     // // BASELINE mode 14 - AWB 18.03.17
     // factories.push_back( std::make_tuple( nullF, nullL, "f_0x2000132f_0x4_invPt",    // dPhi12, 23, theta, FR1/2, St1 ring, dTh13, bend1
     // 					 var_names, var_vals, 0x2000132f, 0x4) );





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
   
   in_vars.push_back( MVA_var( "theta",     "Track #theta",          "int", 'I', -88 ) ); // 0x0000 0001
   in_vars.push_back( MVA_var( "St1_ring2", "St 1 hit in ring 2",    "int", 'I', -88 ) ); // 0x0000 0002  
   in_vars.push_back( MVA_var( "dPhi_12",   "#phi(2) - #phi(1)",     "int", 'I', -88 ) ); // 0x0000 0004
   in_vars.push_back( MVA_var( "dPhi_23",   "#phi(3) - #phi(2)",     "int", 'I', -88 ) ); // 0x0000 0008

   in_vars.push_back( MVA_var( "dPhi_34",   "#phi(4) - #phi(3)",     "int", 'I', -88 ) ); // 0x0000 0010
   in_vars.push_back( MVA_var( "dPhi_13",   "#phi(3) - #phi(1)",     "int", 'I', -88 ) ); // 0x0000 0020
   in_vars.push_back( MVA_var( "dPhi_14",   "#phi(4) - #phi(1)",     "int", 'I', -88 ) ); // 0x0000 0040
   in_vars.push_back( MVA_var( "dPhi_24",   "#phi(4) - #phi(2)",     "int", 'I', -88 ) ); // 0x0000 0080

   in_vars.push_back( MVA_var( "FR_1",      "St 1 LCT F/R",          "int", 'I', -88 ) ); // 0x0000 0100
   in_vars.push_back( MVA_var( "FR_2",      "St 2 LCT F/R",          "int", 'I', -88 ) ); // 0x0000 0200
   in_vars.push_back( MVA_var( "FR_3",      "St 3 LCT F/R",          "int", 'I', -88 ) ); // 0x0000 0400
   in_vars.push_back( MVA_var( "FR_4",      "St 4 LCT F/R",          "int", 'I', -88 ) ); // 0x0000 0800

   in_vars.push_back( MVA_var( "bend_1",    "St 1 LCT bending",      "int", 'I', -88 ) ); // 0x0000 1000
   in_vars.push_back( MVA_var( "bend_2",    "St 2 LCT bending",      "int", 'I', -88 ) ); // 0x0000 2000
   in_vars.push_back( MVA_var( "bend_3",    "St 3 LCT bending",      "int", 'I', -88 ) ); // 0x0000 4000
   in_vars.push_back( MVA_var( "bend_4",    "St 4 LCT bending",      "int", 'I', -88 ) ); // 0x0000 8000

   in_vars.push_back( MVA_var( "dPhiSum4",  "#Sigmad#phi (6)",       "int", 'I', -88 ) ); // 0x0001 0000
   in_vars.push_back( MVA_var( "dPhiSum4A", "#Sigma|d#phi| (6)",     "int", 'I', -88 ) ); // 0x0002 0000
   in_vars.push_back( MVA_var( "dPhiSum3",  "#Sigmad#phi (3)",       "int", 'I', -88 ) ); // 0x0004 0000
   in_vars.push_back( MVA_var( "dPhiSum3A", "#Sigma|d#phi| (3)",     "int", 'I', -88 ) ); // 0x0008 0000
   
   in_vars.push_back( MVA_var( "outStPhi",  "#phi outlier St",       "int", 'I', -88 ) ); // 0x0010 0000
   in_vars.push_back( MVA_var( "filler",    "Filler",                "int", 'I', -88 ) ); // 0x0020 0000
   in_vars.push_back( MVA_var( "dTh_12",    "#theta(2) - #theta(1)", "int", 'I', -88 ) ); // 0x0040 0004
   in_vars.push_back( MVA_var( "dTh_23",    "#theta(3) - #theta(2)", "int", 'I', -88 ) ); // 0x0080 0008

   in_vars.push_back( MVA_var( "dTh_34",    "#theta(4) - #theta(3)", "int", 'I', -88 ) ); // 0x0100 0000
   in_vars.push_back( MVA_var( "dTh_13",    "#theta(3) - #theta(1)", "int", 'I', -88 ) ); // 0x0200 0000
   in_vars.push_back( MVA_var( "dTh_14",    "#theta(4) - #theta(1)", "int", 'I', -88 ) ); // 0x0400 0000
   in_vars.push_back( MVA_var( "dTh_24",    "#theta(4) - #theta(2)", "int", 'I', -88 ) ); // 0x0800 0000

   // in_vars.push_back( MVA_var( "RPC_1",     "St 1 hit is RPC",       "int", 'I', -88 ) ); // 0x1000 0000
   // in_vars.push_back( MVA_var( "RPC_2",     "St 2 hit is RPC",       "int", 'I', -88 ) ); // 0x2000 0000
   // in_vars.push_back( MVA_var( "RPC_3",     "St 3 hit is RPC",       "int", 'I', -88 ) ); // 0x4000 0000
   // in_vars.push_back( MVA_var( "RPC_4",     "St 4 hit is RPC",       "int", 'I', -88 ) ); // 0x8000 0000


   ////////////////////////////////////////////////////////////
   //  Target variable: true muon pT, or 1/pT, or log2(pT)  ///
   ////////////////////////////////////////////////////////////
   
   targ_vars.push_back( MVA_var( "GEN_pt_trg",      "GEN p_{T} for trigger",               "GeV",      'F', -99 ) ); // 0x1
   targ_vars.push_back( MVA_var( "inv_GEN_pt_trg",  "1 / GEN muon p_{T} for trigger",      "GeV^{-1}", 'F', -99 ) ); // 0x2
   targ_vars.push_back( MVA_var( "log2_GEN_pt_trg", "log_{2}(GEN muon p_{T} for trigger)", "GeV",      'F', -99 ) ); // 0x4

   /////////////////////////////////////////////////////////////////////////////
   ///  Spectator variables: not used in training, but saved in output tree  ///
   /////////////////////////////////////////////////////////////////////////////
   
   spec_vars.push_back( MVA_var( "GEN_pt",       "GEN p_{T}",               "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_pt",      "EMTF p_{T}",              "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "inv_GEN_pt",   "1 / GEN muon p_{T}",      "GeV^{-1}", 'F', -77 ) );
   spec_vars.push_back( MVA_var( "inv_EMTF_pt",  "1 / EMTF p_{T}",          "GeV^{-1}", 'F', -77 ) );
   spec_vars.push_back( MVA_var( "log2_GEN_pt",  "log_{2}(GEN muon p_{T})", "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "log2_EMTF_pt", "log_{2}(EMTF p_{T})",     "GeV",      'F', -77 ) );

   spec_vars.push_back( MVA_var( "GEN_eta",       "GEN #eta",                "",         'F', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_eta",      "EMTF #eta",               "",         'F', -77 ) );
   spec_vars.push_back( MVA_var( "TRK_eta",       "Track #eta",              "",         'F', -77 ) );
   spec_vars.push_back( MVA_var( "GEN_phi",       "GEN #phi",                "",         'F', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_phi",      "EMTF #phi",               "",         'F', -77 ) );
   spec_vars.push_back( MVA_var( "TRK_phi",       "Track #phi",              "",         'F', -77 ) );

   spec_vars.push_back( MVA_var( "GEN_charge",    "GEN charge",              "",         'I', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_charge",   "EMTF charge",             "",         'I', -77 ) );

   spec_vars.push_back( MVA_var( "EMTF_mode",     "EMTF mode",               "",         'I', -77 ) );
   spec_vars.push_back( MVA_var( "TRK_mode",      "Track mode",              "",         'I', -77 ) );
   spec_vars.push_back( MVA_var( "SHARED_mode",   "EMTF-track shared mode",  "",         'I', -77 ) );
   spec_vars.push_back( MVA_var( "TRK_mode_CSC",  "Track CSC-only mode",     "",         'I', -77 ) );
   spec_vars.push_back( MVA_var( "TRK_mode_RPC",  "Track RPC-only mode",     "",         'I', -77 ) );

   spec_vars.push_back( MVA_var( "dPhi_sign",    "#phi(B) - #phi(A) sign",  "",         'I', -77 ) );


   assert( in_vars.size() > 0 );   // You need at least one input variable
   assert( targ_vars.size() > 0 ); // You need at least one target variable
   // Order is important: input variables first, then target, then specator
   all_vars.insert( all_vars.end(), in_vars.begin(), in_vars.end() );
   all_vars.insert( all_vars.end(), targ_vars.begin(), targ_vars.end() );
   if (SPEC_VARS) all_vars.insert( all_vars.end(), spec_vars.begin(), spec_vars.end() );


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

   for (int iCh = 0; iCh < in_chains.size(); iCh++) {
     TChain *in_chain = in_chains.at(iCh);
     
     // Get branches from the chain
     TBranch *muon_br  = in_chain->GetBranch("muon");
     TBranch *hit_br   = in_chain->GetBranch("hit");
     TBranch *trk_br = in_chain->GetBranch("track");
     
     std::cout << "\n******* About to enter the event loop for chain " << iCh+1 << " *******" << std::endl;
     
     for (UInt_t jEvt = 0; jEvt < in_chain->GetEntries(); jEvt++) {
       if (iEvt > MAX_EVT) break;

       in_chain->GetEntry(jEvt);
       
       UInt_t nMuons = (muon_br->GetLeaf("nMuons"))->GetValue();
       UInt_t nHits  = (hit_br->GetLeaf("nHits"))->GetValue();
       UInt_t nTrks  = (trk_br->GetLeaf("nTracks"))->GetValue();
       Bool_t isMC     = (nMuons > 0);
       Bool_t trainEvt = true;  // Can use the event for training
       if (not isMC)  // Process ZeroBias anyway
	 nMuons = nTrks;
       // std::cout << "There are " << nMuons << " GEN muons and " << nTrks << " EMTF tracks\n" << std::endl;

       if ( ( (iEvt % REPORT_EVT) == 0 && isMC) || (iEvtZB > 0 && (iEvtZB % REPORT_EVT) == 0) )
	 std::cout << "Looking at MC event " << iEvt << " (ZB = " << iEvtZB << ")" << std::endl;
       
       for (UInt_t iMu = 0; iMu < nMuons; iMu++) {
	 double mu_pt  = 999.;
	 double mu_eta = -99.;
	 double mu_phi = -99.;
	 int mu_charge = -99;
	 if (isMC) {
	   mu_pt     = (muon_br->GetLeaf("pt"))->GetValue(iMu);
	   mu_eta    = (muon_br->GetLeaf("eta"))->GetValue(iMu);
	   mu_phi    = (muon_br->GetLeaf("phi"))->GetValue(iMu);
	   mu_charge = (muon_br->GetLeaf("charge"))->GetValue(iMu);
	 }

	 if ( mu_pt < PTMIN && isMC ) continue;
	 if ( mu_pt > PTMAX && isMC ) continue;
	 if ( fabs( mu_eta ) < ETAMIN && isMC ) continue;
	 if ( fabs( mu_eta ) > ETAMAX && isMC ) continue;
	 if ( mu_pt > PTMAX_TR && isMC ) trainEvt = false;
	 // std::cout << "\nMuon " << iMu+1 << " has pt = " << mu_pt << ", eta = " << mu_eta << ", phi = " << mu_phi << std::endl;


	 // Find the relevant EMTF track
	 double emtf_pt    = 999.;
	 double emtf_eta   = -99.;
	 double emtf_phi   = -99.;
	 int emtf_charge   = -99;
	 int emtf_mode     = -99;
	 int emtf_sect_idx = -99;
	 std::vector<int> emtf_ph = {-99, -99, -99, -99};
	 std::vector<int> emtf_th = {-99, -99, -99, -99};

	 for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

	   // Require same endcap
	   emtf_eta  = (trk_br->GetLeaf("eta"))->GetValue(iTrk);
	   if ((emtf_eta > 0) != (mu_eta > 0)) {
	     emtf_eta = -99.;
	     continue;
	   }

	   // Require valid mode
	   emtf_mode = (trk_br->GetLeaf("mode"))->GetValue(iTrk);
	   bool good_emtf_mode = false;
	   for (UInt_t jMode = 0; jMode < EMTF_MODES.size(); jMode++) {
	     if (emtf_mode == EMTF_MODES.at(jMode))
	       good_emtf_mode = true;
	   }
	   if (!good_emtf_mode) {
	     emtf_mode = -99;
	     continue;
	   }

	   emtf_pt       = (trk_br->GetLeaf("pt"))->GetValue(iTrk);
	   emtf_phi      = (trk_br->GetLeaf("phi"))->GetValue(iTrk);;
	   emtf_charge   = (trk_br->GetLeaf("charge"))->GetValue(iTrk);
	   emtf_sect_idx = (trk_br->GetLeaf("sector_index"))->GetValue(iTrk);
	   
	   for (int ii = 0; ii < 4; ii++) {
	     if (emtf_mode % int(pow(2, 4 - ii)) / pow(2, 3 - ii) > 0) {
	       emtf_ph.at(ii) = (trk_br->GetLeaf("hit_phi_int"))->GetValue(iTrk*4 + ii);
	       emtf_th.at(ii) = (trk_br->GetLeaf("hit_theta_int"))->GetValue(iTrk*4 + ii);
	     }
	   }
	   
	   break; // Only one EMTF track per GEN muon considered

	 } // End loop: for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++)

	 if (REQ_EMTF && emtf_mode < 0)
	   continue;

	 // std::cout << "  * EMTF track has sector_index = " << emtf_sect_idx
	 // 	   << ", eta = " << emtf_eta << ", phi = " << emtf_phi << std::endl;
	 // std::cout << "    - St. 1: theta = " << emtf_ph.at(0) << ", phi = " << emtf_th.at(0) << std::endl;
	 // std::cout << "    - St. 2: theta = " << emtf_ph.at(1) << ", phi = " << emtf_th.at(1) << std::endl;
	 // std::cout << "    - St. 3: theta = " << emtf_ph.at(2) << ", phi = " << emtf_th.at(2) << std::endl;
	 // std::cout << "    - St. 4: theta = " << emtf_ph.at(3) << ", phi = " << emtf_th.at(3) << std::endl;


	 //////////////////////////////////////////
	 ///  Build tracks from available hits  ///
	 //////////////////////////////////////////

	 std::vector< std::vector< std::vector<int> > > id; // All hit index values, by sector and station
	 std::vector< std::vector< std::vector<int> > > ph; // All full-precision integer phi values
	 std::vector< std::vector< std::vector<int> > > th; // All full-precision integer theta values
	 std::vector< std::vector< std::vector<int> > > dt; // All detector values (0 for none, 1 for CSC, 2 for RPC)

	 // Fill vectors with empty sectors / stations
	 for (int ii = 0; ii < 12; ii++) {
	   id.push_back({{}});
	   ph.push_back({{}});
	   th.push_back({{}});
	   dt.push_back({{}});
	   for (int jj = 0; jj < 4; jj++) {
	     id.at(ii).push_back({});
	     ph.at(ii).push_back({});
	     th.at(ii).push_back({});
	     dt.at(ii).push_back({});
	   }
	 }

	 
	 // Loop over all hits
	 for (UInt_t iHit = 0; iHit < nHits; iHit++) {
	   if ( (mu_eta > 0) != ((hit_br->GetLeaf("eta"))->GetValue(iHit) > 0) )
	     continue;
	   int iSc = (hit_br->GetLeaf("sector_index"))->GetValue(iHit) - 1;
	   int iSt = (hit_br->GetLeaf("station"))     ->GetValue(iHit) - 1;
	   id.at(iSc).at(iSt).push_back( iHit );
	   ph.at(iSc).at(iSt).push_back( (hit_br->GetLeaf("phi_int"))  ->GetValue(iHit) );
	   th.at(iSc).at(iSt).push_back( (hit_br->GetLeaf("theta_int"))->GetValue(iHit) );
	   dt.at(iSc).at(iSt).push_back( (hit_br->GetLeaf("isRPC"))    ->GetValue(iHit) ? 2 : 1 );
	 }

	 // Remove masked hits
	 std::vector<std::tuple<int, int, int>> to_erase;
	 for (int ii = 0; ii < 12; ii++) { // Loop over sectors
	   for (int jj = 0; jj < 4; jj++) { // Loop over stations
	     for (int kk = 0; kk < dt.at(ii).at(jj).size(); kk++) { // Loop over hits
	       for (int ll = 0; ll < CSC_MASK.size(); ll++)
		 if (jj+1 == CSC_MASK.at(ll) && dt.at(ii).at(jj).at(kk) == 1)
		   to_erase.push_back(std::make_tuple(ii, jj, kk));
	       for (int ll = 0; ll < RPC_MASK.size(); ll++)
		 if (jj+1 == RPC_MASK.at(ll) && dt.at(ii).at(jj).at(kk) == 2)
		   to_erase.push_back(std::make_tuple(ii, jj, kk));
	     }
	   }
	 }
	 for (int ii = int(to_erase.size()) - 1; ii >= 0; ii--) {
	   int iSc = std::get<0>(to_erase.at(ii));
	   int iSt = std::get<1>(to_erase.at(ii));
	   int iHt = std::get<2>(to_erase.at(ii));
	   id.at(iSc).at(iSt).erase(id.at(iSc).at(iSt).begin() + iHt);
	   ph.at(iSc).at(iSt).erase(ph.at(iSc).at(iSt).begin() + iHt);
	   th.at(iSc).at(iSt).erase(th.at(iSc).at(iSt).begin() + iHt);
	   dt.at(iSc).at(iSt).erase(dt.at(iSc).at(iSt).begin() + iHt);
	 }
	 
	 // Vector indices of hits in each track; 4 in each, stations 1-2-3-4 
	 std::vector< std::vector<int> > all_trk_hits;
	 // Vector of mode, CSC mode, RPC mode, and sumAbsDPhi in each track
	 std::vector< std::vector<int> > all_trk_modes;
      
	 // Build tracks for the specified mode
	 std::vector< std::vector<int> > trk_hits;
	 std::vector< std::vector<int> > trk_modes;
	 
	 BuildTracks( trk_hits, trk_modes, id, ph, th, dt, MODE, MAX_RPC, MIN_CSC, MAX_DPH, MAX_DTH );
	 all_trk_hits.insert( all_trk_hits.end(), trk_hits.begin(), trk_hits.end() );
	 all_trk_modes.insert( all_trk_modes.end(), trk_modes.begin(), trk_modes.end() );

	 // std::cout << "  * Built " << all_trk_hits.size() << " tracks out of " << nHits << " hits" << std::endl;
	 assert(all_trk_modes.size() == all_trk_hits.size());
	 
	 ///////////////////////////////
	 ///  Loop over built tracks ///
	 ///////////////////////////////
	 
	 for (UInt_t iTrk = 0; iTrk < all_trk_hits.size(); iTrk++) {
	   
	   std::vector<int> trk_hits  = all_trk_hits.at(iTrk);
	   std::vector<int> trk_modes = all_trk_modes.at(iTrk);
	   assert(trk_hits.size() == 4);
	   assert(trk_modes.size() == 5);

	   int i1 = trk_hits.at(0);
	   int i2 = trk_hits.at(1);
	   int i3 = trk_hits.at(2);
	   int i4 = trk_hits.at(3);

	   int mode     = trk_modes.at(0);
	   int mode_CSC = trk_modes.at(1);
	   int mode_RPC = trk_modes.at(2);
	   int shared_mode = 0;
	   assert(mode == MODE);

	   // std::cout << "\n    - i1 = " << i1 <<", i2 = " << i2<< ", i3 = " <<i3 << ", i4 = "<< i4 << std::endl;

	   // Properties of hits
	   int ph1 = (i1 >= 0 ? (hit_br->GetLeaf("phi_int"))->GetValue(i1) : -99);
	   int ph2 = (i2 >= 0 ? (hit_br->GetLeaf("phi_int"))->GetValue(i2) : -99);
	   int ph3 = (i3 >= 0 ? (hit_br->GetLeaf("phi_int"))->GetValue(i3) : -99);
	   int ph4 = (i4 >= 0 ? (hit_br->GetLeaf("phi_int"))->GetValue(i4) : -99);

	   // std::cout << "    - ph1 = " << ph1 << ", ph2 = " << ph2 << ", ph3 = " << ph3 << ", ph4 = " << ph4 << std::endl;
	   
	   int th1 = (i1 >= 0 ? (hit_br->GetLeaf("theta_int"))->GetValue(i1) : -99);
	   int th2 = (i2 >= 0 ? (hit_br->GetLeaf("theta_int"))->GetValue(i2) : -99);
	   int th3 = (i3 >= 0 ? (hit_br->GetLeaf("theta_int"))->GetValue(i3) : -99);
	   int th4 = (i4 >= 0 ? (hit_br->GetLeaf("theta_int"))->GetValue(i4) : -99);

	   // std::cout << "    - th1 = " << th1 << ", th2 = " << th2 << ", th3 = " << th3 << ", th4 = " << th4 << std::endl;

	   int pat1 = (i1 >= 0 ? (hit_br->GetLeaf("pattern"))->GetValue(i1) : -99);
	   int pat2 = (i2 >= 0 ? (hit_br->GetLeaf("pattern"))->GetValue(i2) : -99);
	   int pat3 = (i3 >= 0 ? (hit_br->GetLeaf("pattern"))->GetValue(i3) : -99);
	   int pat4 = (i4 >= 0 ? (hit_br->GetLeaf("pattern"))->GetValue(i4) : -99);

	   int st1_ring2 = (i1 >= 0 ? ((hit_br->GetLeaf("ring"))->GetValue(i1) == 2) : -99);

	   double eta;
	   double phi;
	   int endcap;
	   if      (i2 >= 0) { eta = (hit_br->GetLeaf("eta"))->GetValue(i2); phi = (hit_br->GetLeaf("phi"))->GetValue(i2); }
	   else if (i3 >= 0) { eta = (hit_br->GetLeaf("eta"))->GetValue(i3); phi = (hit_br->GetLeaf("phi"))->GetValue(i3); }
	   else if (i4 >= 0) { eta = (hit_br->GetLeaf("eta"))->GetValue(i4); phi = (hit_br->GetLeaf("phi"))->GetValue(i4); }
	   else if (i1 >= 0) { eta = (hit_br->GetLeaf("eta"))->GetValue(i1); phi = (hit_br->GetLeaf("phi"))->GetValue(i1); }
	   endcap = (eta > 0 ? +1 : -1);

	   // Check which hits match between EMTF track and built track
	   if (i1 >= 0 && ph1 == emtf_ph.at(0) && th1 == emtf_th.at(0)) shared_mode += 8;
	   if (i2 >= 0 && ph2 == emtf_ph.at(1) && th2 == emtf_th.at(1)) shared_mode += 4;
	   if (i3 >= 0 && ph3 == emtf_ph.at(2) && th3 == emtf_th.at(2)) shared_mode += 2;
	   if (i4 >= 0 && ph4 == emtf_ph.at(3) && th4 == emtf_th.at(3)) shared_mode += 1;

	   // Variables to go into BDT
	   int theta;
	   int dPh12, dPh13, dPh14, dPh23, dPh24, dPh34, dPhSign;
	   int dPhSum4, dPhSum4A, dPhSum3, dPhSum3A, outStPh;
	   int dTh12, dTh13, dTh14, dTh23, dTh24, dTh34;
	   int FR1, FR2, FR3, FR4;
	   int bend1, bend2, bend3, bend4;
	   int RPC1, RPC2, RPC3, RPC4;

	   // std::cout << "    - Computing theta" << std::endl;
	   theta = CalcTrackTheta( th1, th2, th3, th4, st1_ring2, mode, BIT_COMP );
	   
	   // std::cout << "    - Computing dPhis" << std::endl;
	   CalcDeltaPhis( dPh12, dPh13, dPh14, dPh23, dPh24, dPh34, dPhSign,
			  dPhSum4, dPhSum4A, dPhSum3, dPhSum3A, outStPh,
			  ph1, ph2, ph3, ph4, mode, BIT_COMP );
	   
	   // std::cout << "    - Computing dThetas" << std::endl;
	   CalcDeltaThetas( dTh12, dTh13, dTh14, dTh23, dTh24, dTh34,
			    th1, th2, th3, th4, mode, BIT_COMP );

	   // std::cout << "    - Computing FRs" << std::endl;
	   FR1 = (i1 >= 0 ? (hit_br->GetLeaf("FR"))->GetValue(i1) : -99);
	   FR2 = (i2 >= 0 ? (hit_br->GetLeaf("FR"))->GetValue(i2) : -99);
	   FR3 = (i3 >= 0 ? (hit_br->GetLeaf("FR"))->GetValue(i3) : -99);
	   FR4 = (i4 >= 0 ? (hit_br->GetLeaf("FR"))->GetValue(i4) : -99);

	   // std::cout << "    - Computing bend" << std::endl;
	   CalcBends( bend1, bend2, bend3, bend4,
		      pat1, pat2, pat3, pat4, 
		      dPhSign, endcap, mode, BIT_COMP );

	   // std::cout << "    - Computing RPCs" << std::endl;
	   RPC1 = (i1 >= 0 ? ((hit_br->GetLeaf("isRPC"))->GetValue(i1) == 1 ? 1 : 0) : -99);
	   RPC2 = (i2 >= 0 ? ((hit_br->GetLeaf("isRPC"))->GetValue(i2) == 1 ? 1 : 0) : -99);
	   RPC3 = (i3 >= 0 ? ((hit_br->GetLeaf("isRPC"))->GetValue(i3) == 1 ? 1 : 0) : -99);
	   RPC4 = (i4 >= 0 ? ((hit_br->GetLeaf("isRPC"))->GetValue(i4) == 1 ? 1 : 0) : -99);
	   
	   // Clean out showering muons with outlier station 1, or >= 2 outlier stations
	   if (isMC && log2(mu_pt) > 6 && CLEAN_HI_PT)
	     if ( dPhSum4A >= std::max(40., 332. - 40*log2(mu_pt)) )
	       if ( outStPh < 2 || dPhSum3A >= std::max(24., 174. - 20*log2(mu_pt)) )
		 trainEvt = false;
	   

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
		 var_vals.at(iVar) = theta;
	       if ( vName == "St1_ring2" )
		 var_vals.at(iVar) = st1_ring2;

	       if ( vName == "dPhi_12" )
		 var_vals.at(iVar) = dPh12;
	       if ( vName == "dPhi_13" )
		 var_vals.at(iVar) = dPh13;
	       if ( vName == "dPhi_14" )
		 var_vals.at(iVar) = dPh14;
	       if ( vName == "dPhi_23" )
		 var_vals.at(iVar) = dPh23;
	       if ( vName == "dPhi_24" )
		 var_vals.at(iVar) = dPh24;
	       if ( vName == "dPhi_34" )
		 var_vals.at(iVar) = dPh34;
	       
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
		 var_vals.at(iVar) = dPhSum4;
	       if ( vName == "dPhiSum4A" )
		 var_vals.at(iVar) = dPhSum4A;
	       if ( vName == "dPhiSum3" )
		 var_vals.at(iVar) = dPhSum3;
	       if ( vName == "dPhiSum3A" )
		 var_vals.at(iVar) = dPhSum3A;
	       if ( vName == "outStPhi" )
		 var_vals.at(iVar) = outStPh;

	       if ( vName == "dTh_12" )
		 var_vals.at(iVar) = dTh12;
	       if ( vName == "dTh_13" )
		 var_vals.at(iVar) = dTh13;
	       if ( vName == "dTh_14" )
		 var_vals.at(iVar) = dTh14;
	       if ( vName == "dTh_23" )
		 var_vals.at(iVar) = dTh23;
	       if ( vName == "dTh_24" )
		 var_vals.at(iVar) = dTh24;
	       if ( vName == "dTh_34" )
		 var_vals.at(iVar) = dTh34;

	       if ( vName == "RPC_1" )
		 var_vals.at(iVar) = RPC1;
	       if ( vName == "RPC_2" )
		 var_vals.at(iVar) = RPC2;
	       if ( vName == "RPC_3" )
		 var_vals.at(iVar) = RPC3;
	       if ( vName == "RPC_4" )
		 var_vals.at(iVar) = RPC4;

	       
	       //////////////////////////////
	       ///  Target and variables  ///
	       //////////////////////////////

	       if ( vName == "GEN_pt_trg" )
		 var_vals.at(iVar) = min(mu_pt, PTMAX_TRG);
	       if ( vName == "inv_GEN_pt_trg" )
		 var_vals.at(iVar) = 1. / min(mu_pt, PTMAX_TRG);
	       if ( vName == "log2_GEN_pt_trg" )
		 var_vals.at(iVar) = log2(min(mu_pt, PTMAX_TRG));

	       /////////////////////////////
	       ///  Spectator variables  ///
	       /////////////////////////////

	       if ( vName == "GEN_pt" )
		 var_vals.at(iVar) = mu_pt;
	       if ( vName == "EMTF_pt" )
	       	 var_vals.at(iVar) = emtf_pt;
	       if ( vName == "inv_GEN_pt" )
		 var_vals.at(iVar) = 1. / mu_pt;
	       if ( vName == "inv_EMTF_pt" )
	       	 var_vals.at(iVar) = 1. / emtf_pt;
	       if ( vName == "log2_GEN_pt" )
		 var_vals.at(iVar) = log2(mu_pt);
	       if ( vName == "log2_EMTF_pt" )
	       	 var_vals.at(iVar) = (emtf_pt > 0 ? log2(emtf_pt) : -99);

	       if ( vName == "GEN_eta" )
		 var_vals.at(iVar) = mu_eta;
	       if ( vName == "EMTF_eta" )
	       	 var_vals.at(iVar) = emtf_eta;
	       if ( vName == "TRK_eta" )
	       	 var_vals.at(iVar) = eta;
	       if ( vName == "GEN_phi" )
		 var_vals.at(iVar) = mu_phi;
	       if ( vName == "EMTF_phi" )
	       	 var_vals.at(iVar) = emtf_phi;
	       if ( vName == "TRK_phi" )
	       	 var_vals.at(iVar) = phi;

	       if ( vName == "GEN_charge" )
		 var_vals.at(iVar) = mu_charge;
	       if ( vName == "EMTF_charge" )
	       	 var_vals.at(iVar) = emtf_charge;

	       if ( vName == "EMTF_mode" )
	       	 var_vals.at(iVar) = emtf_mode;
	       if ( vName == "TRK_mode" )
	       	 var_vals.at(iVar) = mode;
	       if ( vName == "SHARED_mode" )
	       	 var_vals.at(iVar) = shared_mode;
	       if ( vName == "TRK_mode_CSC" )
	       	 var_vals.at(iVar) = mode_CSC;
	       if ( vName == "TRK_mode_RPC" )
	       	 var_vals.at(iVar) = mode_RPC;

	       if ( vName == "dPhi_sign" )
		 var_vals.at(iVar) = dPhSign;

	       
	     } // End loop: for (UInt_t iVar = 0; iVar < var_names.size(); iVar++)
	     
	     // Unweighted distribution: flat in eta and 1/pT
	     Double_t evt_weight = 1.0;
	     // Weight by 1/pT so overall distribution is (1/pT)^2
	     if ( std::get<2>(factories.at(iFact)).Contains("_invPt") )
	       evt_weight = 1. / mu_pt;
	     if ( std::get<2>(factories.at(iFact)).Contains("_invPtSq") )
	       evt_weight = 1. / pow(mu_pt, 2);


	     // Load values into event
	     if ( (iEvt % 2) == 0 && isMC && trainEvt && nTrain < (MAX_TR - (iFact == 0)) && (!USE_RPC || iCh < (nChains_ZB + nChains_RPC)) ) { 
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
   PtRegression_AWB_v2(methodList);
   return 0;
}

