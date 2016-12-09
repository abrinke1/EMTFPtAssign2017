/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides examples for the training and testing of the
/// TMVA classifiers.
///
/// As input data is used a toy-MC sample consisting of four Gaussian-distributed
/// and linearly correlated input variables.
///
/// The methods to be used can be switched on and off by means of booleans, or
/// via the prompt command, for example:
///
///     root -l PtRegression_AWB_v0.C\(\"LD,MLP\"\)
///
/// (note that the backslashes are mandatory)
/// If no method given, a default set is used.
///
/// The output file "TMVAReg.root" can be analysed with the use of dedicated
/// macros (simply say: root -l <macro.C>), which can be conveniently
/// invoked through a GUI that will appear at the end of the run of this macro.
/// - Project   : TMVA - a Root-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: PtRegression_AWB_v0
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker

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

const int MAX_EVT =   10000;
const int REPORT_EVT = 1000;

const double PI = 3.14159265359;
const double PT_SCALE = (1. / 1.25); // EMTF pT was scaled up by ~1.25 in 2016 (for 4-hit tracks, mode 15)
const double BIT = 0.000001; // Tiny value or offset

using namespace TMVA;

void PtRegression_AWB_v0( TString myMethodList = "" )
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   //     mylinux~> root -l PtRegression_AWB_v0.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //

   //---------------------------------------------------------------
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
   Use["BDT"]                 = 0;
   Use["BDTG_default"]        = 1;
   Use["BDTG_Carnes_AbsDev"]  = 0;
   Use["BDTG_Carnes_Huber"]   = 0;
   Use["BDTG_Carnes_LeastSq"] = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start PtRegression_AWB_v0" << std::endl;

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
   TString outfileName( "PtRegression_AWB_v0_16_12_09.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory will
   // then run the performance analysis for you.
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );

   ////////////////////////////////////////
   // Advanced data loading, event by event
   ////////////////////////////////////////
   TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

   // Read training and test data (see TMVAClassification for reading ASCII files)
   // load the signal and background event samples from ROOT trees
   TFile *input(0);
   std::vector<TString> fnames;
   fnames.push_back("/afs/cern.ch/work/a/abrinke1/public/EMTF/Analyzer/ntuples/EMTF_MC_NTuple_SingleMu_noRPC_300k.root");

   for (UInt_t i = 0; i < fnames.size(); i++) {
     if ( !gSystem->AccessPathName(fnames.at(i)) )
       input = TFile::Open( fnames.at(i) ); // check if file in local directory exists

     if (!input) {
       std::cout << "ERROR: could not open data file " << fnames.at(i) << std::endl;
       exit(1);
     }
   }

   // Have to use TChain for both SetBranchAddress and GetEntry to work
   TChain *regChain = new TChain("ntuple/tree");
   for (UInt_t i = 0; i < fnames.size(); i++) {
     regChain->Add( fnames.at(i) );
     std::cout << "--- TMVARegression           : Using input file: " << fnames.at(i) << std::endl;
   }

   TBranch *muon_br  = regChain->GetBranch("muon");
   TBranch *hit_br   = regChain->GetBranch("hit");
   TBranch *track_br = regChain->GetBranch("track");


   std::vector<MVA_var> in_vars;   // All input variables
   std::vector<MVA_var> targ_vars; // All target variables (should only define 1, unless using MLP)
   std::vector<MVA_var> spec_vars; // All spectator variables
   std::vector<MVA_var> all_vars;  // All variables

   // Defined in interface/MVA_helper.h
   // MVA_var(TString name, TString descr, TString unit, TString type, Double_t def_val)

   /////////////////////////////////////////////////////////
   ///  Input variables: used in BDT to estimate the pT  ///
   /////////////////////////////////////////////////////////
   
   in_vars.push_back( MVA_var( "theta",     "Track #theta",      "deg", 'F', -88 ) );
   in_vars.push_back( MVA_var( "dPhi_12",   "#phi(2) - #phi(1)", "deg", 'F', -88 ) ); 
   in_vars.push_back( MVA_var( "dPhi_23",   "#phi(3) - #phi(2)", "deg", 'F', -88 ) ); 
   in_vars.push_back( MVA_var( "dPhi_34",   "#phi(4) - #phi(3)", "deg", 'F', -88 ) );
   // in_vars.push_back( MVA_var( "dPhi_13",   "#phi(3) - #phi(1)", "deg", 'F', -88 ) );
   // in_vars.push_back( MVA_var( "dPhi_14",   "#phi(4) - #phi(1)", "deg", 'F', -88 ) );
   // in_vars.push_back( MVA_var( "dPhi_24",   "#phi(4) - #phi(2)", "deg", 'F', -88 ) );

   // in_vars.push_back( MVA_var( "bend_1",    "St 1 LCT bending",  "int", 'I', -88 ) );
   // in_vars.push_back( MVA_var( "bend_2",    "St 2 LCT bending",  "int", 'I', -88 ) );
   // in_vars.push_back( MVA_var( "bend_3",    "St 3 LCT bending",  "int", 'I', -88 ) );
   // in_vars.push_back( MVA_var( "bend_4",    "St 4 LCT bending",  "int", 'I', -88 ) );

   // in_vars.push_back( MVA_var( "dPhiSum4",  "#Sigmad#phi (6)",   "deg", 'F', -88 ) );
   // in_vars.push_back( MVA_var( "dPhiSum4A", "#Sigma|d#phi| (6)", "deg", 'F', -88 ) );
   // in_vars.push_back( MVA_var( "dPhiSum3",  "#Sigmad#phi (3)",   "deg", 'F', -88 ) );
   // in_vars.push_back( MVA_var( "dPhiSum3A", "#Sigma|d#phi| (3)", "deg", 'F', -88 ) );
   // in_vars.push_back( MVA_var( "outSt",     "Outlier station",   "int", 'I', -88 ) );

   ////////////////////////////////////////////////////////////
   //  Target variable: true muon pT, or 1/pT, or log2(pT)  ///
   ////////////////////////////////////////////////////////////
   
   // targ_vars.push_back( MVA_var( "GEN_pt",       "GEN p_{T}",               "GeV",      'F', -99 ) );
   targ_vars.push_back( MVA_var( "inv_GEN_pt",   "1 / GEN muon p_{T}",      "GeV^{-1}", 'F', -99 ) );
   // targ_vars.push_back( MVA_var( "log2_GEN_pt",  "log_{2}(GEN muon p_{T})", "GeV",      'F', -99 ) );

   /////////////////////////////////////////////////////////////////////////////
   ///  Spectator variables: not used in training, but saved in output tree  ///
   /////////////////////////////////////////////////////////////////////////////
   
   spec_vars.push_back( MVA_var( "GEN_pt",       "GEN p_{T}",               "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "EMTF_pt",      "EMTF p_{T}",              "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "inv_GEN_pt",   "1 / GEN muon p_{T}",      "GeV^{-1}", 'F', -77 ) );
   spec_vars.push_back( MVA_var( "inv_EMTF_pt",  "1 / EMTF p_{T}",          "GeV^{-1}", 'F', -77 ) );
   spec_vars.push_back( MVA_var( "log2_GEN_pt",  "log_{2}(GEN muon p_{T})", "GeV",      'F', -77 ) );
   spec_vars.push_back( MVA_var( "log2_EMTF_pt", "log_{2}(EMTF p_{T})",     "GeV",      'F', -77 ) );
   for (UInt_t i = 0; i < targ_vars.size(); i++)
     for (UInt_t j = 0; j < spec_vars.size(); j++)
       if (spec_vars.at(j).name == targ_vars.at(i).name)
	 spec_vars.erase(spec_vars.begin() + j);
	 
   assert( in_vars.size() > 0 );    // You need at least one input variable
   assert( targ_vars.size() == 1 ); // You should only define one target variable
   // Order is important: input variables first, then target, then specator
   all_vars.insert( all_vars.end(), in_vars.begin(), in_vars.end() );
   all_vars.insert( all_vars.end(), targ_vars.begin(), targ_vars.end() );
   all_vars.insert( all_vars.end(), spec_vars.begin(), spec_vars.end() );

   // Fill factory with the correct set of variables
   std::cout << "\n*** Input variables ***" << std::endl;
   for (UInt_t i = 0; i < in_vars.size(); i++) {
     MVA_var v = in_vars.at(i);
     std::cout << v.name << std::endl;
     dataloader->AddVariable( v.name, v.descr, v.unit, v.type ); }
   std::cout << "\n*** Target variables ***" << std::endl;
   for (UInt_t i = 0; i < targ_vars.size(); i++) {
     MVA_var v = targ_vars.at(i);
     std::cout << v.name << std::endl;
     dataloader->AddTarget( v.name, v.descr, v.unit, v.type ); }
   std::cout << "\n*** Spectator variables ***" << std::endl;
   for (UInt_t i = 0; i < spec_vars.size(); i++) {
     MVA_var v = spec_vars.at(i);
     std::cout << v.name << std::endl;
     dataloader->AddSpectator( v.name, v.descr, v.unit, v.type ); }

   std::vector<Double_t> var_vals; // Holds values of variables


   std::cout << "\n******* About to enter the event loop *******" << std::endl;
   UInt_t nTrain = 0;
   UInt_t nTest  = 0;
   for (UInt_t iEvt = 0; iEvt < regChain->GetEntries(); iEvt++) {
     if (iEvt > MAX_EVT) break;
     if ( (iEvt % REPORT_EVT) == 0 )
       std::cout << "Looking at event " << iEvt << std::endl;
     regChain->GetEntry(iEvt);

     UInt_t nMuons  = (muon_br->GetLeaf("nMuons"))->GetValue();
     UInt_t nTracks = (track_br->GetLeaf("nTracks"))->GetValue();
     // std::cout << "There are " << nMuons << " GEN muons and " << nTracks << " EMTF tracks\n" << std::endl;
     
     for (UInt_t iMu = 0; iMu < nMuons; iMu++) {
       Double_t mu_pt  = (muon_br->GetLeaf("pt"))->GetValue(iMu);
       Double_t mu_eta = (muon_br->GetLeaf("eta"))->GetValue(iMu);
       Double_t mu_phi = (muon_br->GetLeaf("phi"))->GetValue(iMu);
       if ( fabs( mu_eta ) < 1.25 ) continue;
       if ( fabs( mu_eta ) > 2.35 ) continue;
       // if ( mu_pt < 10 ) continue;
       // std::cout << "Muon " << iMu+1 << " has pt = " << mu_pt << ", eta = " << mu_eta << ", phi = " << mu_phi << std::endl;
       
       for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++) {
	 Double_t trk_pt    = (track_br->GetLeaf("pt"))->GetValue(iTrk);
	 Double_t trk_eta   = (track_br->GetLeaf("eta"))->GetValue(iTrk);
	 trk_pt *= PT_SCALE;
	 Double_t trk_theta = (track_br->GetLeaf("theta"))->GetValue(iTrk);
	 Double_t trk_phi   = (track_br->GetLeaf("phi"))->GetValue(iTrk);
	 int trk_mode       = (track_br->GetLeaf("mode"))->GetValue(iTrk);
	 if ( trk_mode != 15 ) continue;
	 if ( ( mu_eta > 0 ) != ( trk_eta > 0 ) ) continue;
	 // std::cout << "  * Track " << iTrk+1 << " has pt = " << trk_pt << ", eta = " << trk_eta 
	 // 	   << ", phi = " << mu_phi << ", mode = " << trk_mode << std::endl;
	 
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
	 
	 Double_t dPhi_sum_4  = dPhi12 + dPhi13 + dPhi14 + dPhi23 + dPhi24 + dPhi34;
	 Double_t dPhi_sum_4A = fabs(dPhi12) + fabs(dPhi13) + fabs(dPhi14) + fabs(dPhi23) + fabs(dPhi24) + fabs(dPhi34);
	 Double_t dev_st1 = fabs(dPhi12) + fabs(dPhi13) + fabs(dPhi14);
	 Double_t dev_st2 = fabs(dPhi12) + fabs(dPhi23) + fabs(dPhi24);
	 Double_t dev_st3 = fabs(dPhi13) + fabs(dPhi23) + fabs(dPhi34);
	 Double_t dev_st4 = fabs(dPhi14) + fabs(dPhi24) + fabs(dPhi34);

	 Int_t out_st = -88;
	 if      (dev_st4 > dev_st3 && dev_st4 > dev_st2 && dev_st4 > dev_st1)  out_st = 4;
	 else if (dev_st3 > dev_st4 && dev_st3 > dev_st2 && dev_st3 > dev_st1)  out_st = 3;
	 else if (dev_st2 > dev_st4 && dev_st2 > dev_st3 && dev_st2 > dev_st1)  out_st = 2;
	 else if (dev_st1 > dev_st4 && dev_st1 > dev_st3 && dev_st1 > dev_st2)  out_st = 1;
	 else                                                                   out_st = 0;

	 Double_t dPhi_sum_3  = -88;
	 Double_t dPhi_sum_3A = -88;
	 if      (out_st == 4) {
	   dPhi_sum_3  = dPhi12 + dPhi13 + dPhi23;
	   dPhi_sum_3A = fabs(dPhi12) + fabs(dPhi13) + fabs(dPhi23);
	 } else if (out_st == 3) {
	   dPhi_sum_3  = dPhi12 + dPhi14 + dPhi24;
	   dPhi_sum_3A = fabs(dPhi12) + fabs(dPhi14) + fabs(dPhi24);
	 } else if (out_st == 2) {
	   dPhi_sum_3  = dPhi13 + dPhi14 + dPhi34;
	   dPhi_sum_3A = fabs(dPhi13) + fabs(dPhi14) + fabs(dPhi34);
	 } else {
	   dPhi_sum_3  = dPhi23 + dPhi24 + dPhi34;
	   dPhi_sum_3A = fabs(dPhi23) + fabs(dPhi24) + fabs(dPhi34);
	 }

	 Int_t bend1 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 0), trk_eta );
	 Int_t bend2 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 1), trk_eta );
	 Int_t bend3 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 2), trk_eta );
	 Int_t bend4 = CalcBendFromPattern( (track_br->GetLeaf("hit_pattern"))->GetValue(4*iTrk + 3), trk_eta );


	 // Initialize all variables
	 var_vals.clear();
	 for (UInt_t i = 0; i < all_vars.size(); i++)
	   var_vals.push_back ( all_vars.at(i).def_val );
	 
	 // Fill all variables
	 for (UInt_t iVar = 0; iVar < all_vars.size(); iVar++) {
	   TString vName = all_vars.at(iVar).name;

	   /////////////////////////
	   ///  Input variables  ///
	   /////////////////////////

	   if ( vName == "theta" )
	     var_vals.at(iVar) = fabs(trk_theta);
	   if ( vName == "dPhi_12" )
	     var_vals.at(iVar) = dPhi12;
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
	   if ( vName == "outSt" )
	     var_vals.at(iVar) = out_st;


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

	 } // End loop: for (UInt_t iVar = 0; iVar < all_vars.size(); iVar++)

	 Double_t evt_weight = 1.0;
	 // // Weight by 1/pT so overall distribution is (1/pT)^2
	 // Double_t evt_weight = 1. / mu_pt;

	 // Load values into event
	 if ( (iEvt % 2) == 0 ) {
	   dataloader->AddTrainingEvent( "Regression", var_vals, evt_weight );
	   nTrain += 1;
	 }
	 else {
	   dataloader->AddTestEvent( "Regression", var_vals, evt_weight );
	   nTest += 1;
	 }

       } // End loop: for (UInt_t iTrk = 0; iTrk < nTracks; iTrk++)
     } // End loop: for (UInt_t iMu = 0; iMu < nMuons; iMu++)

   } // End loop: for (UInt_t i = 0; i < regChain->GetEntries(); i++)

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

   // // You can add an arbitrary number of regression trees
   // dataloader->AddRegressionTree( regTree, regWeight );

   // // This would set individual event weights (the variables defined in the
   // // expression need to exist in the original TTree)
   // dataloader->SetWeightExpression( "var1", "Regression" );
   dataloader->SetWeightExpression( 1.0 );

   // // Apply additional cuts on the signal and background samples (can be different)
   // TCut mycut = "( abs(muon.eta[0]) > 1.25 && abs(muon.eta[1]) < 2.4 )"; // && track.mode[0] == 15 )"; 

   // Set nTest_Regression to 0 to tell the DataLoader to use all remaining events in the trees after training for testing:
   dataloader->PrepareTrainingAndTestTree( "", numTrainStr+"SplitMode=Random:NormMode=NumEvents:!V" );   
   // dataloader->PrepareTrainingAndTestTree( mycut, "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // If no numbers of events are given, half of the events in the tree are used
   // for training, and the other half for testing:
   //
   //     dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );

   // Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // PDE - RS method
   if (Use["PDERS"])
      factory->BookMethod( dataloader,  TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=40:NEventsMax=60:VarTransform=None" );
   // And the options strings for the MinMax and RMS methods, respectively:
   //
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );

   if (Use["PDEFoam"])
       factory->BookMethod( dataloader,  TMVA::Types::kPDEFoam, "PDEFoam",
 			    "!H:!V:MultiTargetRegression=F:TargetSelection=Mpv:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Compress=T:Kernel=None:Nmin=10:VarTransform=None" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( dataloader,  TMVA::Types::kKNN, "KNN",
                           "nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // Linear discriminant
   if (Use["LD"])
      factory->BookMethod( dataloader,  TMVA::Types::kLD, "LD",
                           "!H:!V:VarTransform=None" );

 	// Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_MC",
                          "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=MC:SampleSize=100000:Sigma=0.1:VarTransform=D" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options) .. the formula of this example is good for parabolas
      factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_GA",
                           "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:PopSize=100:Cycles=3:Steps=30:Trim=True:SaveBestGen=1:VarTransform=Norm" );

   if (Use["FDA_MT"])
      factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_MT",
                           "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( dataloader,  TMVA::Types::kFDA, "FDA_GAMT",
                           "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   // Neural network (MLP)
   if (Use["MLP"])
      factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );

   if (Use["DNN"])
   {
   /*
       TString layoutString ("Layout=TANH|(N+100)*2,LINEAR");
       TString layoutString ("Layout=SOFTSIGN|100,SOFTSIGN|50,SOFTSIGN|20,LINEAR");
       TString layoutString ("Layout=RELU|300,RELU|100,RELU|30,RELU|10,LINEAR");
       TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|30,SOFTSIGN|20,SOFTSIGN|10,LINEAR");
       TString layoutString ("Layout=TANH|50,TANH|30,TANH|20,TANH|10,LINEAR");
       TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|20,LINEAR");
       TString layoutString ("Layout=TANH|100,TANH|30,LINEAR");
    */
       TString layoutString ("Layout=TANH|100,LINEAR");

       TString training0 ("LearningRate=1e-5,Momentum=0.5,Repetitions=1,ConvergenceSteps=500,BatchSize=50,TestRepetitions=7,WeightDecay=0.01,Regularization=NONE,DropConfig=0.5+0.5+0.5+0.5,DropRepetitions=2");
       TString training1 ("LearningRate=1e-5,Momentum=0.9,Repetitions=1,ConvergenceSteps=170,BatchSize=30,TestRepetitions=7,WeightDecay=0.01,Regularization=L2,DropConfig=0.1+0.1+0.1,DropRepetitions=1");
       TString training2 ("LearningRate=1e-5,Momentum=0.3,Repetitions=1,ConvergenceSteps=150,BatchSize=40,TestRepetitions=7,WeightDecay=0.01,Regularization=NONE");
       TString training3 ("LearningRate=1e-6,Momentum=0.1,Repetitions=1,ConvergenceSteps=500,BatchSize=100,TestRepetitions=7,WeightDecay=0.0001,Regularization=NONE");

       TString trainingStrategyString ("TrainingStrategy=");
       trainingStrategyString += training0 + "|" + training1 + "|" + training2 + "|" + training3;


 //       TString trainingStrategyString ("TrainingStrategy=LearningRate=1e-1,Momentum=0.3,Repetitions=3,ConvergenceSteps=20,BatchSize=30,TestRepetitions=7,WeightDecay=0.0,L1=false,DropFraction=0.0,DropRepetitions=5");

       TString nnOptions ("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G:WeightInitialization=XAVIERUNIFORM");
 //       TString nnOptions ("!H:V:VarTransform=Normalize:ErrorStrategy=CHECKGRADIENTS");
       nnOptions.Append (":"); nnOptions.Append (layoutString);
       nnOptions.Append (":"); nnOptions.Append (trainingStrategyString);

       factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN", nnOptions ); // NN
   }



   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( dataloader,  TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   // Boosted Decision Trees
   if (Use["BDT"])
     factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=100:MinNodeSize=1.0%:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );

     // Default TMVA settings ... error with MaxDepth? - AWB 07.12.16
   if (Use["BDTG_default"])
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG_default",
 			  "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:MaxDepth=4" );
   // A. Carnes option #1
   if (Use["BDTG_Carnes_AbsDev"])
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG_Carnes_AbsDev",
 			  "!H:!V:NTrees=64::BoostType=Grad:Shrinkage=0.3:nCuts=99999:MaxDepth=4:MinNodeSize=0.001:NegWeightTreatment=IgnoreNegWeightsInTraining:PruneMethod=NoPruning:RegressionLossFunctionBDTG=AbsoluteDeviation" );

   if (Use["BDTG_Carnes_Huber"])
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG_Carnes_Huber",
 			  "!H:!V:NTrees=64::BoostType=Grad:Shrinkage=0.3:nCuts=99999:MaxDepth=4:MinNodeSize=0.001:NegWeightTreatment=IgnoreNegWeightsInTraining:PruneMethod=NoPruning:RegressionLossFunctionBDTG=Huber" );

   if (Use["BDTG_Carnes_LeastSq"])
     factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG_Carnes_LeastSq",
 			  "!H:!V:NTrees=64::BoostType=Grad:Shrinkage=0.3:nCuts=99999:MaxDepth=4:MinNodeSize=0.001:NegWeightTreatment=IgnoreNegWeightsInTraining:PruneMethod=NoPruning:RegressionLossFunctionBDTG=LeastSquares" );

   // --------------------------------------------------------------------------------------------------

   // Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;

   delete factory;
   delete dataloader;

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
   PtRegression_AWB_v0(methodList);
   return 0;
}

