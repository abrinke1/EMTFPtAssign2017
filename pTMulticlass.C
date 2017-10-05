//====================================================================
//=     This macro is based on
//=     https://root.cern.ch/doc/v608/TMVAMulticlass_8C_source.html
//=     and PtRegression_Apr_2017.C in this repository.
//=     It is a recipe for the EMTF pT training and
//=     testing of the TMVA multiclass classification.
//=
//=     Author: Wei Shi
//====================================================================

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAMultiClassGui.h"

// Extra tools
#include "interface/MVA_helper.h"
#include "src/TrackBuilder.cc"
#include "src/PtLutVarCalc.cc"

// Configuration settings
#include "configs/pTMulticlass/Standard.h" // Settings that are not likely to change
#include "configs/pTMulticlass/General.h"  // General settings relevant for all modes
#include "configs/pTMulticlass/User.h"     // Specific settings for each user
#include "configs/pTMulticlass/Modes.h"    // Specific settigns for each mode

using namespace TMVA;
void pTMulticlass( TString myMethodList = "" ){
    
    // This loads the library
    TMVA::Tools::Instance();
    
    // Default MVA methods to be trained + tested
    std::map<std::string,int> Use;
    Use["MLP"]             = 1;
    Use["BDTG"]            = 1;
    
    std::cout << std::endl;
    std::cout << "==> Start pTMulticlass" << std::endl;
    
    // Select methods (don't look at this code - not of interest)
    std::vector<TString> mlist;
    if (myMethodList != "") {
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
        
        mlist = TMVA::gTools().SplitString( myMethodList, ',' );
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
    
    // Configure settings for this mode and user
    pTMulticlass_cfg::ConfigureMode( MODE );
    pTMulticlass_cfg::ConfigureUser( USER );

    // Create a new root output file
    TString out_file_str;
    TString bit_str = (BIT_COMP ? "bitCompr" : "noBitCompr");
    TString RPC_str = (USE_RPC  ? "RPC"      : "noRPC");

    out_file_str.Form( "%s/%s_MODE_%d_%s_%s.root", OUT_DIR_NAME.Data(), OUT_FILE_NAME.Data(), MODE, bit_str.Data(), RPC_str.Data() );
    TFile* out_file = TFile::Open( out_file_str, "RECREATE" );
    
    // Read training and test data
    // load the signal and background event samples from ROOT trees
    TFile *input(0);
    
    std::vector<TString> in_file_names;
    TString in_file_name;
    
    int nZB_in = 0;
    if (USE_RPC) {
      TString in_dir = "ZeroBiasIsolatedBunch0/Slim_RPC/170213_174254/0000";
      for (int j = 20; j < 50; j++) {  // First 19 files are empty
        if (nZB_in >= MAX_ZB_FIL) break;
        in_file_name.Form("%s/%s/tuple_%d.root", EOS_DIR_NAME.Data(), in_dir.Data(), j);
        std::cout << "Adding file " << in_file_name.Data() << std::endl;
        in_file_names.push_back(in_file_name.Data());
        nZB_in += 1;
      }
    } else {
      TString in_dirs[4] = { 
                "ZeroBiasIsolatedBunch0/Slim/170130_224405/0000",
   			    "ZeroBiasIsolatedBunch1/Slim/170130_175144/0000",
   			    "ZeroBiasIsolatedBunch4/Slim/170130_175005/0000",
   			    "ZeroBiasIsolatedBunch5/Slim/170130_174947/0000"};     
      for (int i = 0; i < 4; i++) {
        for (int j = 1; j < 50; j++) {
            if (nZB_in >= MAX_ZB_FIL) break;
   	        in_file_name.Form("%s/%s/tuple_%d.root", EOS_DIR_NAME.Data(), in_dirs[i].Data(), j);
   	        std::cout << "Adding file " << in_file_name.Data() << std::endl;
   	        in_file_names.push_back(in_file_name.Data());
   	        nZB_in += 1;
        }
      }
    }//end else

    // Load files with RPC hits
    TString in_dir_RPC = "SingleMu_Pt1To1000_FlatRandomOneOverPt/RPC/170213_173255/0000";
    for (int i = 1; i < 99; i++) {
      if (!USE_RPC) continue;
      in_file_name.Form("%s/%s/tuple_%d.root", EOS_DIR_NAME.Data(), in_dir_RPC.Data(), i);
      std::cout << "Adding file " << in_file_name.Data() << std::endl;
      in_file_names.push_back(in_file_name.Data());
      if (i*100000 > MAX_EVT) break; // ~100k events per file
    }
   
    // Load files without RPC hits
    TString in_dir_CSC = "SingleMu_Pt1To1000_FlatRandomOneOverPt/EMTF_MuGun/170113_165434/0000";
    for (int i = 1; i < 99; i++) {
      if (USE_RPC) continue;
      in_file_name.Form("%s/%s/EMTF_MC_NTuple_SingleMu_noRPC_%d.root", EOS_DIR_NAME.Data(), in_dir_CSC.Data(), i);
      std::cout << "Adding file " << in_file_name.Data() << std::endl;
      in_file_names.push_back(in_file_name.Data());
      if (i*100000 > MAX_EVT) break; // ~100k events per file
    }
 
    for (UInt_t i = 0; i < in_file_names.size(); i++) {
      if ( !gSystem->AccessPathName(in_file_names.at(i)) )
        input = TFile::Open( in_file_names.at(i) ); // check if file in local directory exists
      if (!input) {
        std::cout << "ERROR: could not open data file " << in_file_names.at(i) << std::endl;
        in_file_names.erase( in_file_names.begin()+i );
        if (i < nZB_in) 
            nZB_in -= 1;
        i -= 1;
      }
    }
    
    // Add trees from the input files to the TChain
    std::vector<TChain*> in_chains;
    int nChains_ZB  = -99;
    for (UInt_t i = 0; i < in_file_names.size(); i++) {
      TChain *tmp_chain = new TChain("ntuple/tree");
      tmp_chain->Add( in_file_names.at(i) );
      in_chains.push_back(tmp_chain);
      if ( i == (nZB_in - 1) && nChains_ZB < 0)
          nChains_ZB = i + 1;
    }
    
    TString fact_set = "!V:!Silent:Color:DrawProgressBar:AnalysisType=multiclass";
    std::vector<TString> var_names; // Holds names of variables for a given factory and permutation
    std::vector<Double_t> var_vals; // Holds values of variables for a given factory and permutation
    TMVA::Factory* nullF = new TMVA::Factory("NULL", out_file, fact_set); // Placeholder factory
    TMVA::DataLoader* nullL = new TMVA::DataLoader("NULL");               // Placeholder loader
    
    // Tuple is defined by the factory and dataloader,  followed by a name, 
    // var name and value vectors, and hex bit masks for input variables.
    // Each hex bit represents four variables, e.g. 0x1 would select only the 1st variable, 
    // 0xf the 1st 4, 0xff the 1st 8, 0xa the 2nd and 4th, 0xf1 the 1st and 5th-8th, etc.
    std::vector< std::tuple<TMVA::Factory*, TMVA::DataLoader*, TString, std::vector<TString>, std::vector<Double_t>, int> > factories;
    
    for (int iWgt = 0; iWgt < EVT_WGTS.size(); iWgt++) {
       TString factName;  // "Wgt" components not arbitrary - correspond to specific options later on
       factName.Form( "f_MODE_%d_%sWgt_%s_%s", MODE, EVT_WGTS.at(iWgt).Data(), bit_str.Data(), RPC_str.Data() );

       // 4-station tracks
       if        (MODE == 15) {
	 // BASELINE mode 15 - dPhi12/23/34 + combos, theta, FR1, St1 ring, dTh14, bend1, RPC 1/2/3/4
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0xf41f11ff) );
       } 

       // 3-station tracks
       else if   (MODE == 14) {
	 // BASELINE mode 14 - dPhi12/23/13, theta, FR1/2, St1 ring, dTh13, bend1, RPC 1/2/3
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0x7200132f) );
       } else if (MODE == 13) {
	 // BASELINE mode 13 - dPhi12/24/14, theta, FR1/2, St1 ring, dTh14, bend1, RPC 1/2/4
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0xb40013c7) );
       } else if (MODE == 11) {
	 // BASELINE mode 11 - dPhi13/34/14, theta, FR1/3, St1 ring, dTh14, bend1, RPC 1/3/4
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0xd4001573) );
       } else if (MODE ==  7) {
	 // BASELINE mode  7 - dPhi23/34/24, theta, FR2, dTh24, bend2, RPC 2/3/4
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0xe8002299) );
       } 

       // 2-station tracks
       else if   (MODE == 12) {
	 // BASELINE mode 12 - dPhi12, theta, FR1/2, St1 ring, dTh12, bend1/2, RPC 1/2
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0x30403307) );
       } else if (MODE == 10) {
	 // BASELINE mode 10 - dPhi13, theta, FR1/3, St1 ring, dTh13, bend1/3, RPC 1/3
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0x52005523) );
       } else if (MODE ==  9) {
	 // BASELINE mode  9 - dPhi14, theta, FR1/4, St1 ring, dTh14, bend1/4, RPC 1/4
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0x94009943) );
       } else if (MODE ==  6) {
	 // BASELINE mode  6 - dPhi23, theta, FR2/3, dTh23, bend2/3, RPC 2/3
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0x60806609) );
       } else if (MODE ==  5) {
	 // BASELINE mode  5 - dPhi24, theta, FR2/4, dTh24, bend2/4, RPC 2/4
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0xa800aa81) );
       } else if (MODE ==  3) {
	 // BASELINE mode  3 - dPhi34, theta, FR3/4, dTh12, bend3/4, RPC 3/4
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0xc100cc11) );
       }
       else if (MODE == 0) {
	 // Null track, for testing EMTF performance
	 factories.push_back( std::make_tuple( nullF, nullL, factName, var_names, var_vals, 0xc0000001) );
       }
    }//end evt weight loop
	
    // Initialize factories and dataloaders
    for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
	    std::get<0>(factories.at(iFact)) = new TMVA::Factory( std::get<2>(factories.at(iFact)), out_file, fact_set );
	    std::get<1>(factories.at(iFact)) = new TMVA::DataLoader( std::get<2>(factories.at(iFact)) );
    }
	
    // Defined in interface/MVA_helper.h
    // MVA_var(TString name, TString descr, TString unit, TString type, Double_t def_val)
    std::vector<MVA_var> in_vars;   // All input variables
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
   in_vars.push_back( MVA_var( "dTh_12",    "#theta(2) - #theta(1)", "int", 'I', -88 ) ); // 0x0040 0000
   in_vars.push_back( MVA_var( "dTh_23",    "#theta(3) - #theta(2)", "int", 'I', -88 ) ); // 0x0080 0000

   in_vars.push_back( MVA_var( "dTh_34",    "#theta(4) - #theta(3)", "int", 'I', -88 ) ); // 0x0100 0000
   in_vars.push_back( MVA_var( "dTh_13",    "#theta(3) - #theta(1)", "int", 'I', -88 ) ); // 0x0200 0000
   in_vars.push_back( MVA_var( "dTh_14",    "#theta(4) - #theta(1)", "int", 'I', -88 ) ); // 0x0400 0000
   in_vars.push_back( MVA_var( "dTh_24",    "#theta(4) - #theta(2)", "int", 'I', -88 ) ); // 0x0800 0000

   if (USE_RPC) {
     in_vars.push_back( MVA_var( "RPC_1",     "St 1 hit is RPC",       "int", 'I', -88 ) ); // 0x1000 0000
     in_vars.push_back( MVA_var( "RPC_2",     "St 2 hit is RPC",       "int", 'I', -88 ) ); // 0x2000 0000
     in_vars.push_back( MVA_var( "RPC_3",     "St 3 hit is RPC",       "int", 'I', -88 ) ); // 0x4000 0000
     in_vars.push_back( MVA_var( "RPC_4",     "St 4 hit is RPC",       "int", 'I', -88 ) ); // 0x8000 0000
   }
    
   /////////////////////////////////////////////////////////////////////////////
   ///  Spectator variables: not used in training, but saved in output tree  ///
   /////////////////////////////////////////////////////////////////////////////
   
   spec_vars.push_back( MVA_var( "GEN_pt",       "GEN p_{T}",               "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_pt",      "EMTF p_{T}",              "GeV",      'F', -77 ) );
   // spec_vars.push_back( MVA_var( "inv_GEN_pt",   "1 / GEN muon p_{T}",      "GeV^{-1}", 'F', -77 ) );
   // spec_vars.push_back( MVA_var( "inv_EMTF_pt",  "1 / EMTF p_{T}",          "GeV^{-1}", 'F', -77 ) );
   // spec_vars.push_back( MVA_var( "log2_GEN_pt",  "log_{2}(GEN muon p_{T})", "GeV",      'F', -77 ) );
   // spec_vars.push_back( MVA_var( "log2_EMTF_pt", "log_{2}(EMTF p_{T})",     "GeV",      'F', -77 ) );

   spec_vars.push_back( MVA_var( "GEN_eta",       "GEN #eta",                "", 'F', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_eta",      "EMTF #eta",               "", 'F', -77 ) );
   spec_vars.push_back( MVA_var( "TRK_eta",       "Track #eta",              "", 'F', -77 ) );
   // spec_vars.push_back( MVA_var( "GEN_phi",       "GEN #phi",                "", 'F', -77 ) );
   // spec_vars.push_back( MVA_var( "EMTF_phi",      "EMTF #phi",               "", 'F', -77 ) );
   // spec_vars.push_back( MVA_var( "TRK_phi",       "Track #phi",              "", 'F', -77 ) );
   spec_vars.push_back( MVA_var( "GEN_charge",    "GEN charge",              "", 'I', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_charge",   "EMTF charge",             "", 'I', -77 ) );

   spec_vars.push_back( MVA_var( "EMTF_mode",     "EMTF mode",                   "", 'I', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_mode_CSC", "EMTF CSC-only mode",          "", 'I', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_mode_RPC", "EMTF RPC-only",               "", 'I', -77 ) );
   spec_vars.push_back( MVA_var( "TRK_mode",      "Track mode",                  "", 'I', -77 ) );
   spec_vars.push_back( MVA_var( "TRK_mode_CSC",  "Track CSC-only mode",         "", 'I', -77 ) );
   spec_vars.push_back( MVA_var( "TRK_mode_RPC",  "Track RPC-only mode",         "", 'I', -77 ) );
   // spec_vars.push_back( MVA_var( "SHRD_mode",     "EMTF-track shared mode",      "", 'I', -77 ) );
   // spec_vars.push_back( MVA_var( "SHRD_mode_CSC", "EMTF-track shared CSC mode",  "", 'I', -77 ) );
   // spec_vars.push_back( MVA_var( "SHRD_mode_RPC", "EMTF-track shared RPC mode",  "", 'I', -77 ) );

   spec_vars.push_back( MVA_var( "dPhi_sign",  "#phi(B) - #phi(A) sign",    "", 'I', -77 ) );
   // spec_vars.push_back( MVA_var( "nTRK",       "Number of tracks built",    "", 'I', -77 ) );
   spec_vars.push_back( MVA_var( "evt_weight", "Event weight for training", "", 'F', -77 ) );
	
    assert( in_vars.size() > 0 );   // You need at least one input variable
    // Order is important: input variables first, then target, then specator
    all_vars.insert( all_vars.end(), in_vars.begin(), in_vars.end() );
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
    
    std::cout << "*** Spectator ***" << std::endl;
     for (UInt_t i = 0; i < spec_vars.size(); i++) {
       MVA_var v = spec_vars.at(i);
       std::cout << v.name << std::endl;
       std::get<1>(factories.at(iFact))->AddSpectator( v.name, v.descr, v.unit, v.type );
       std::get<3>(factories.at(iFact)).push_back( v.name );
       std::get<4>(factories.at(iFact)).push_back( v.def_val );
     }
   } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++)

   std::cout << "\n******* About to loop over chains *******" << std::endl;
   
   UInt_t iEvt = 0;
   UInt_t iEvtZB = 0;

   for (int iCh = 0; iCh < in_chains.size(); iCh++) {
     TChain *in_chain = in_chains.at(iCh);
     
     // Get branches from the chain
     TBranch *muon_br  = in_chain->GetBranch("muon");
     TBranch *hit_br   = in_chain->GetBranch("hit");
     TBranch *trk_br = in_chain->GetBranch("track");
     
     std::cout << "\n******* About to enter the event loop for chain " << iCh+1 << " *******" << std::endl;
     
     for (UInt_t jEvt = 0; jEvt < in_chain->GetEntries(); jEvt++) {
       if (iEvt > MAX_EVT) break; 
       //use all MC events

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
	 std::cout << "Looking at MC event " << iEvt << " (ZeroBias event " << iEvtZB << ")" << std::endl;
       
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

	 if ( isMC && (mu_pt < PTMIN || mu_pt > PTMAX) ) continue;
	 if ( isMC && (fabs( mu_eta ) < ETAMIN || fabs( mu_eta ) > ETAMAX) ) continue;
	 if ( isMC && (mu_pt < PTMIN_TR || mu_pt > PTMAX_TR) ) trainEvt = false;
	 // std::cout << "\nMuon " << iMu+1 << " has pt = " << mu_pt << ", eta = " << mu_eta << ", phi = " << mu_phi << std::endl;


	 // Find the relevant EMTF track
	 double emtf_pt    = 999.;
	 double emtf_eta   = -99.;
	 int emtf_eta_int  = -99;
	 double emtf_phi   = -99.;
	 int emtf_charge   = -99;
	 int emtf_mode     = -99;
	 int emtf_mode_CSC = -99;
	 int emtf_mode_RPC = -99;
	 int emtf_sect_idx = -99;
	 std::array<int, 4> emtf_id = {-99, -99, -99, -99};
	 std::array<int, 4> emtf_ph = {-99, -99, -99, -99};
	 std::array<int, 4> emtf_th = {-99, -99, -99, -99};
	 std::array<int, 4> emtf_dt = {-99, -99, -99, -99};

	 for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

	   // Require same endcap
	   emtf_eta     = (trk_br->GetLeaf("eta"))->GetValue(iTrk);
	   emtf_eta_int = (trk_br->GetLeaf("eta_int"))->GetValue(iTrk);
	   if (isMC && (emtf_eta > 0) != (mu_eta > 0)) {
	     emtf_eta = -99.;
	     continue;
	   }

	   emtf_mode = (trk_br->GetLeaf("mode"))->GetValue(iTrk);
	   emtf_mode_CSC = 0;
	   emtf_mode_RPC = 0;

	   // Require valid mode
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
	     if (emtf_mode < 0)
	       continue;
	     if ( (emtf_mode % int(pow(2, 4 - ii))) / int(pow(2, 3 - ii)) > 0) {
	       emtf_id.at(ii) = iTrk*4 + ii;
	       emtf_ph.at(ii) = (trk_br->GetLeaf("hit_phi_int"))->GetValue(iTrk*4 + ii);
	       emtf_th.at(ii) = (trk_br->GetLeaf("hit_theta_int"))->GetValue(iTrk*4 + ii);
	       emtf_dt.at(ii) = ( (trk_br->GetLeaf("hit_isRPC"))->GetValue(iTrk*4 + ii) == 1 ? 2 : 1);
	       if (emtf_dt.at(ii) == 1)
		 emtf_mode_CSC += int(pow(2, 3 - ii));
	       else if (emtf_dt.at(ii) == 2)
		 emtf_mode_RPC += int(pow(2, 3 - ii));
	     }
	   }

	   if (emtf_mode_CSC + emtf_mode_RPC != emtf_mode) {
	     std::cout << "\n\n*** Super-bizzare case where EMTF mode = " << emtf_mode  << ", but CSC mode = " 
		       << emtf_mode_CSC << " and RPC mode = " << emtf_mode_RPC << " ***" << std::endl;
	     std::cout << "  - Rare bug in EMTF emulator - skipping.\n\n" << std::endl;
	     continue;
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

	 std::array< std::array< std::vector<int>, 4>, 12> id; // All hit index values, by sector and station
	 std::array< std::array< std::vector<int>, 4>, 12> ph; // All full-precision integer phi values
	 std::array< std::array< std::vector<int>, 4>, 12> th; // All full-precision integer theta values
	 std::array< std::array< std::vector<int>, 4>, 12> dt; // All detector values (0 for none, 1 for CSC, 2 for RPC)
	 std::vector<bool> emtf_found = {false, false, false, false}; // Check if hits in EMTF track were found in hits

	 // Fill hits with LCTs from the EMTF track, rather than all the LCTs in the event
	 if (USE_EMTF_CSC && emtf_mode > 0) {
	   for (int ii = 0; ii < 12; ii++) {
	     for (int jj = 0; jj < 4; jj++) {
	       if ( (trk_br->GetLeaf("hit_sector_index"))->GetValue(emtf_id.at(jj)) == ii+1 &&
		    emtf_dt.at(jj) == 1 ) {
		 id.at(ii).at(jj).push_back( emtf_id.at(jj) );
		 ph.at(ii).at(jj).push_back( emtf_ph.at(jj) );
		 th.at(ii).at(jj).push_back( emtf_th.at(jj) );
		 dt.at(ii).at(jj).push_back( emtf_dt.at(jj) );
		 // std::cout << "In sector " << ii+1 << ", station " << jj+1 << ", adding hit with "
		 // 	   << "phi = " << emtf_ph.at(jj) << ", theta = " << emtf_th.at(jj) << std::endl;
	       }
	     }
	   } // End loop over stations
	 } // End loop over sector indices

	 
	 // Loop over all hits
	 for (UInt_t iHit = 0; iHit < nHits; iHit++) {
	   if ( (mu_eta > 0) != ((hit_br->GetLeaf("eta"))->GetValue(iHit) > 0) )
	     continue;
	   int iSc = (hit_br->GetLeaf("sector_index"))->GetValue(iHit) - 1;
	   int iSt = (hit_br->GetLeaf("station"))     ->GetValue(iHit) - 1;
	   int iPh = (hit_br->GetLeaf("phi_int"))     ->GetValue(iHit);
	   int iTh = (hit_br->GetLeaf("theta_int"))   ->GetValue(iHit);
	   int iDt = (hit_br->GetLeaf("isRPC"))       ->GetValue(iHit) ? 2 : 1;

	   if (USE_EMTF_CSC && iDt == 1) {
	     if (id.at(iSc).at(iSt).size() > 0) {
	       assert( USE_RPC || id.at(iSc).at(iSt).size() == 1 ); // There should only be one LCT per station
	       if ( ph.at(iSc).at(iSt).at(0) == iPh &&
		    th.at(iSc).at(iSt).at(0) == iTh &&
		    dt.at(iSc).at(iSt).at(0) == iDt ) {
		 id.at(iSc).at(iSt).at(0) = iHit; // Change the index to the hit_br index
		 emtf_found.at(iSt) = true;       // Hit in EMTF track was found in general collection
	       }
	     }
	     continue; // Only look at CSC LCTs if they were included in the EMTF track
	   }

	   id.at(iSc).at(iSt).push_back( iHit );
	   ph.at(iSc).at(iSt).push_back( iPh );
	   th.at(iSc).at(iSt).push_back( iTh );
	   dt.at(iSc).at(iSt).push_back( iDt );
	 }

	 bool found_all_EMTF_LCTs = true;
	 for (int ii = 0; ii < 4; ii++) {
	   if (emtf_dt.at(ii) == 1 && !emtf_found.at(ii))
	     found_all_EMTF_LCTs = false;
	 }
	 if (USE_EMTF_CSC && !found_all_EMTF_LCTs) {
	   // std::cout << "\n  * Rare case where not all LCTs in EMTF track were in the hit collection\n" << std::endl;
	   continue;
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
	 
	 // Array indices of hits in each track; 4 in each, stations 1-2-3-4 
	 std::vector< std::array<int, 4> > all_trk_hits;
	 // Array of mode, CSC mode, RPC mode, sumAbsDPhi, and sumAbsDTheta in each track
	 std::vector< std::array<int, 5> > all_trk_modes;

	 if (MODE > 0) {
	   // Build tracks for the specified mode
	   BuildTracks( all_trk_hits, all_trk_modes, id, ph, th, dt, MODE, MAX_RPC, MIN_CSC, MAX_DPH, MAX_DTH );
	   // std::cout << "  * Built " << all_trk_hits.size() << " tracks out of " << nHits << " hits" << std::endl;
	   assert(all_trk_modes.size() == all_trk_hits.size());
	 } else { 
	   // Skip track building, just store EMTF info
	   all_trk_hits.push_back({-99, -99, -99, -99});
	   all_trk_modes.push_back({0, 0, 0, 0, 0});
	 }
	 
	 ///////////////////////////////
	 ///  Loop over built tracks ///
	 ///////////////////////////////
	 
	 for (UInt_t iTrk = 0; iTrk < all_trk_hits.size(); iTrk++) {
	   
	   std::array<int, 4> trk_hits  = all_trk_hits.at(iTrk);
	   std::array<int, 5> trk_modes = all_trk_modes.at(iTrk);

	   int i1 = trk_hits.at(0);
	   int i2 = trk_hits.at(1);
	   int i3 = trk_hits.at(2);
	   int i4 = trk_hits.at(3);

	   int mode     = trk_modes.at(0);
	   int mode_CSC = trk_modes.at(1);
	   int mode_RPC = trk_modes.at(2);
	   int shared_mode     = 0;
	   int shared_mode_CSC = 0;
	   int shared_mode_RPC = 0;
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

	   int st1_ring2 = (i1 >= 0 ? ((hit_br->GetLeaf("ring"))->GetValue(i1) == 2 || (hit_br->GetLeaf("ring"))->GetValue(i1) == 3) : 0);

	   double eta;
	   double phi;
	   int endcap;
	   if      (i2 >= 0) { eta = (hit_br->GetLeaf("eta"))->GetValue(i2); phi = (hit_br->GetLeaf("phi"))->GetValue(i2); }
	   else if (i3 >= 0) { eta = (hit_br->GetLeaf("eta"))->GetValue(i3); phi = (hit_br->GetLeaf("phi"))->GetValue(i3); }
	   else if (i4 >= 0) { eta = (hit_br->GetLeaf("eta"))->GetValue(i4); phi = (hit_br->GetLeaf("phi"))->GetValue(i4); }
	   else if (i1 >= 0) { eta = (hit_br->GetLeaf("eta"))->GetValue(i1); phi = (hit_br->GetLeaf("phi"))->GetValue(i1); }
	   endcap = (eta > 0 ? +1 : -1);

	   // Check which hits match between EMTF track and built track
	   if (i1 >= 0 && ph1 == emtf_ph.at(0) && th1 == emtf_th.at(0)) {
	     shared_mode     += 8;
	     shared_mode_CSC += 8 * ((hit_br->GetLeaf("isRPC"))->GetValue(i1) == 0);
	     shared_mode_RPC += 8 * ((hit_br->GetLeaf("isRPC"))->GetValue(i1) == 1);
	   }
	   if (i2 >= 0 && ph2 == emtf_ph.at(1) && th2 == emtf_th.at(1)) {
	     shared_mode     += 4;
	     shared_mode_CSC += 4 * ((hit_br->GetLeaf("isRPC"))->GetValue(i2) == 0);
	     shared_mode_RPC += 4 * ((hit_br->GetLeaf("isRPC"))->GetValue(i2) == 1);
	   }
	   if (i3 >= 0 && ph3 == emtf_ph.at(2) && th3 == emtf_th.at(2)) {
	     shared_mode     += 2;
	     shared_mode_CSC += 2 * ((hit_br->GetLeaf("isRPC"))->GetValue(i3) == 0);
	     shared_mode_RPC += 2 * ((hit_br->GetLeaf("isRPC"))->GetValue(i3) == 1);
	   }
	   if (i4 >= 0 && ph4 == emtf_ph.at(3) && th4 == emtf_th.at(3)) {
	     shared_mode     += 1;
	     shared_mode_CSC += 1 * ((hit_br->GetLeaf("isRPC"))->GetValue(i4) == 0);
	     shared_mode_RPC += 1 * ((hit_br->GetLeaf("isRPC"))->GetValue(i4) == 1);
	   }


	   // Variables to go into BDT
	   int theta;
	   int dPh12, dPh13, dPh14, dPh23, dPh24, dPh34, dPhSign;
	   int dPhSum4, dPhSum4A, dPhSum3, dPhSum3A, outStPh;
	   int dTh12, dTh13, dTh14, dTh23, dTh24, dTh34;
	   int FR1, FR2, FR3, FR4;
	   int bend1, bend2, bend3, bend4;
	   int RPC1, RPC2, RPC3, RPC4;

	   // Extra variables for FR computation
	   int ring1, cham1, cham2, cham3, cham4;

	   if (MODE == 0) {
	     theta = emtf_eta_int;
	     goto EMTF_ONLY;
	   }

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

	   // // FR bit directly out of the NTuples
	   // FR1 = (i1 >= 0 ? (hit_br->GetLeaf("FR"))->GetValue(i1) : -99);
	   // FR2 = (i2 >= 0 ? (hit_br->GetLeaf("FR"))->GetValue(i2) : -99);
	   // FR3 = (i3 >= 0 ? (hit_br->GetLeaf("FR"))->GetValue(i3) : -99);
	   // FR4 = (i4 >= 0 ? (hit_br->GetLeaf("FR"))->GetValue(i4) : -99);

	   // In firmware, RPC 'FR' bit set according to FR of corresponding CSC chamber
	   ring1 = (i1 >= 0 ? (hit_br->GetLeaf("ring"))   ->GetValue(i1) : -99);
	   cham1 = (i1 >= 0 ? (hit_br->GetLeaf("chamber"))->GetValue(i1) : -99);
	   cham2 = (i2 >= 0 ? (hit_br->GetLeaf("chamber"))->GetValue(i2) : -99);
	   cham3 = (i3 >= 0 ? (hit_br->GetLeaf("chamber"))->GetValue(i3) : -99);
	   cham4 = (i4 >= 0 ? (hit_br->GetLeaf("chamber"))->GetValue(i4) : -99);

	   FR1 = (i1 >= 0 ? (cham1 % 2 == 0) : -99);  // Odd chambers are bolted to the iron,
	   FR2 = (i2 >= 0 ? (cham2 % 2 == 0) : -99);  // which faces forwared in stations 1 & 2,
	   FR3 = (i3 >= 0 ? (cham3 % 2 == 1) : -99);  // backwards in 3 & 4
	   FR4 = (i4 >= 0 ? (cham4 % 2 == 1) : -99);
	   if (ring1 == 3) FR1 = 0;                   // In ME1/3 chambers are non-overlapping

	   // std::cout << "    - Computing bend" << std::endl;
	   CalcBends( bend1, bend2, bend3, bend4,
		      pat1, pat2, pat3, pat4, 
		      dPhSign, endcap, mode, BIT_COMP );

	   // std::cout << "    - Computing RPCs" << std::endl;
	   RPC1 = (i1 >= 0 ? ((hit_br->GetLeaf("isRPC"))->GetValue(i1) == 1 ? 1 : 0) : -99);
	   RPC2 = (i2 >= 0 ? ((hit_br->GetLeaf("isRPC"))->GetValue(i2) == 1 ? 1 : 0) : -99);
	   RPC3 = (i3 >= 0 ? ((hit_br->GetLeaf("isRPC"))->GetValue(i3) == 1 ? 1 : 0) : -99);
	   RPC4 = (i4 >= 0 ? ((hit_br->GetLeaf("isRPC"))->GetValue(i4) == 1 ? 1 : 0) : -99);

	   CalcRPCs( RPC1, RPC2, RPC3, RPC4, mode, st1_ring2, theta, BIT_COMP );
	   
	   // Clean out showering muons with outlier station 1, or >= 2 outlier stations
	   if (isMC && log2(mu_pt) > 6 && CLEAN_HI_PT && MODE == 15)
	     if ( dPhSum4A >= fmax(40., 332. - 40*log2(mu_pt)) )
	       if ( outStPh < 2 || dPhSum3A >= fmax(24., 174. - 20*log2(mu_pt)) )
		 trainEvt = false;

	 EMTF_ONLY: // Skip track building, just store EMTF info

	   /////////////////////////////////////////////////////
	   ///  Loop over factories and set variable values  ///
	   /////////////////////////////////////////////////////
	   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
	     
	     // Set vars equal to default vector of variables for this factory
	     var_names = std::get<3>(factories.at(iFact));
	     var_vals = std::get<4>(factories.at(iFact));

	     // Unweighted distribution: flat in eta and 1/pT
	     Double_t evt_weight = 1.0;
	     
	     // Weight by 1/pT or (1/pT)^2 so overall distribution is (1/pT)^2 or (1/pT)^3
	     //add more weights forms to study the effect
	     if      ( std::get<2>(factories.at(iFact)).Contains("_Pt0p5Wgt") )
	       evt_weight = pow(mu_pt,0.5);
	     else if ( std::get<2>(factories.at(iFact)).Contains("_log2PtWgt") )
	       evt_weight = log2(mu_pt + BIT);
	     else if ( std::get<2>(factories.at(iFact)).Contains("_PtWgt") )
	       evt_weight = mu_pt;  
	     else if ( std::get<2>(factories.at(iFact)).Contains("_PtSqWgt") )
	       evt_weight = pow(mu_pt, 2);
	     else if ( std::get<2>(factories.at(iFact)).Contains("_invPt0p5Wgt") )
	       evt_weight = 1. / pow(mu_pt, 0.5);
             else if ( std::get<2>(factories.at(iFact)).Contains("_invlog2PtWgt") )
	       evt_weight = 1. / log2(mu_pt + BIT); //mu_pt+ BIT offset in case of zero weight
	     else if ( std::get<2>(factories.at(iFact)).Contains("_invPtWgt") )
	       evt_weight = 1. / mu_pt;
	     else if ( std::get<2>(factories.at(iFact)).Contains("_invPt1p5Wgt") )
	       evt_weight = 1. / pow(mu_pt, 1.5);
	     else if ( std::get<2>(factories.at(iFact)).Contains("_invPtSqWgt") )
	       evt_weight = 1. / pow(mu_pt, 2);
	     else if ( std::get<2>(factories.at(iFact)).Contains("_invPt2p5Wgt") )
	       evt_weight = 1. / pow(mu_pt, 2.5);
             else if ( std::get<2>(factories.at(iFact)).Contains("_invPtCubWgt") )
	       evt_weight = 1. / pow(mu_pt, 3);
	     else if ( std::get<2>(factories.at(iFact)).Contains("_invPtQuadWgt") )
	       evt_weight = 1. / pow(mu_pt, 4);
	     else
	       assert( std::get<2>(factories.at(iFact)).Contains("_noWgt") );
	     
	     // Weight by number of tracks in the event
	     evt_weight *= (1. / all_trk_hits.size());
	     
	     // De-weight tracks with one or more RPC hits
	     evt_weight *= (1. / pow( 4, ((RPC1 == 1) + (RPC2 == 1) + (RPC3 == 1) + (RPC4 == 1)) ) );

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
	       if ( vName == "EMTF_mode_CSC" )
	       	 var_vals.at(iVar) = emtf_mode_CSC;
	       if ( vName == "EMTF_mode_RPC" )
	       	 var_vals.at(iVar) = emtf_mode_RPC;
	       if ( vName == "TRK_mode" )
	       	 var_vals.at(iVar) = mode;
	       if ( vName == "TRK_mode_CSC" )
	       	 var_vals.at(iVar) = mode_CSC;
	       if ( vName == "TRK_mode_RPC" )
	       	 var_vals.at(iVar) = mode_RPC;
	       if ( vName == "SHRD_mode" )
	       	 var_vals.at(iVar) = shared_mode;
	       if ( vName == "SHRD_mode_CSC" )
	       	 var_vals.at(iVar) = shared_mode_CSC;
	       if ( vName == "SHRD_mode_RPC" )
	       	 var_vals.at(iVar) = shared_mode_RPC;

	       if ( vName == "dPhi_sign" )
		 var_vals.at(iVar) = dPhSign;
	       if ( vName == "nTRK" )
		 var_vals.at(iVar) = all_trk_hits.size();
	       if ( vName == "evt_weight" )
		 var_vals.at(iVar) = evt_weight;

	       
	     } // End loop: for (UInt_t iVar = 0; iVar < var_names.size(); iVar++)
	     
	     // Load values into event
	     if ( (iEvt % 2) == 0 && isMC && trainEvt ) {
	       if (mu_pt > 30) {
		 std::get<1>(factories.at(iFact))->AddTrainingEvent( "class1", var_vals, evt_weight );
	       } else if (mu_pt<30 && mu_pt>20) {
		 std::get<1>(factories.at(iFact))->AddTrainingEvent( "class2", var_vals, evt_weight );
	       } else {
		 std::get<1>(factories.at(iFact))->AddTrainingEvent( "class3", var_vals, evt_weight );
	       } 
	     }
	     else {
	       if (mu_pt > 30) {
		 std::get<1>(factories.at(iFact))->AddTestEvent( "class1", var_vals, evt_weight );
	       } else if (mu_pt<30 && mu_pt>20){
		 std::get<1>(factories.at(iFact))->AddTestEvent( "class2", var_vals, evt_weight );
	       } else {
		 std::get<1>(factories.at(iFact))->AddTestEvent( "class3", var_vals, evt_weight );  
	       } 
	     }//end if else
	     
	   } // End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++) 
	   
	 } // End loop: for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++)
       } // End loop: for (UInt_t iMu = 0; iMu < nMuons; iMu++)
	     
       if (isMC) iEvt += 1;
       else iEvtZB += 1;
	     
     } // End loop: for (UInt_t jEvt = 0; jEvt < in_chain->GetEntries(); jEvt++)
   } // End loop: for (int iCh = 0; iCh < in_chains.size(); iCh++) {

   std::cout << "******* Made it out of the event loop *******" << std::endl;

   for (UInt_t iFact = 0; iFact < factories.size(); iFact++) {
     
     TMVA::Factory* factX = std::get<0>(factories.at(iFact));
     TMVA::DataLoader* loadX = std::get<1>(factories.at(iFact));
         
     // Tell the dataloader how to use the training and testing events
     loadX->PrepareTrainingAndTestTree( "", "NormMode=None:!V" );   
	   
    // Book MVA methods
    if (Use["MLP"]) // Multilayer Perceptron
        factX->BookMethod( loadX,  TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:NCycles=500:HiddenLayers=N,5:Sampling=0.7:SamplingEpoch=250:SamplingImportance=0.8:TestRate=5:LearningRate=0.01:DecayRate=0.005:EstimatorType=MSE:ConvergenceTests=10:UseRegulator=True:VarTransform=N,P");
    if (Use["BDTG"]) // Gradient Boosted Decision Trees
        factX->BookMethod( loadX,  TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.000001:RegressionLossFunctionBDTG=LeastSquares");
    
    // Train MVAs using the set of training events
    factX->TrainAllMethods();
    // Evaluate all MVAs using the set of test events
    factX->TestAllMethods();
    // Evaluate and compare performance of all configured MVAs
    //factX->EvaluateAllMethods();
	  
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
	   
   }// End loop: for (UInt_t iFact = 0; iFact < factories.size(); iFact++)
	
    // Save the output
    out_file->Close();
    std::cout << "==> Wrote root file: " << out_file->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;
    //delete factory;
    //delete dataloader;
    
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVAMultiClassGui( out_file_str );
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
    pTMulticlass(methodList);
    return 0;
}
