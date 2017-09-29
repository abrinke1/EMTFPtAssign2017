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
    Use["SVM"]             = 0;
    Use["MLP"]             = 0;
    Use["BDTG"]            = 1;
    Use["DNN"]             = 0;
    
    std::cout << std::endl;
    std::cout << "==> Start pTMulticlass" << std::endl;
    
    // Select methods (don't look at this code - not of interest)
    if (myMethodList != "") {
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;
        
        std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
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
    
    TString fact_set = "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=multiclass";
    std::vector<TString> var_names; // Holds names of variables for a given factory and permutation
    std::vector<Double_t> var_vals; // Holds values of variables for a given factory and permutation
    TMVA::Factory* nullF = new TMVA::Factory("NULL", out_file, fact_set); // Placeholder factory
    TMVA::DataLoader* nullL = new TMVA::DataLoader("NULL");               // Placeholder loader
    
    // Tuple is defined by the factory and dataloader,  followed by a name, 
    // var name and value vectors, and hex bit masks for input variables.
    // Each hex bit represents four variables, e.g. 0x1 would select only the 1st variable, 
    // 0xf the 1st 4, 0xff the 1st 8, 0xa the 2nd and 4th, 0xf1 the 1st and 5th-8th, etc.
    std::vector< std::tuple<TMVA::Factory*, TMVA::DataLoader*, TString, std::vector<TString>, std::vector<Double_t>, int> > factories;
    
    
    TTree *signalTree  = (TTree*)input->Get("TreeS");
    TTree *background0 = (TTree*)input->Get("TreeB0");
    TTree *background1 = (TTree*)input->Get("TreeB1");
    TTree *background2 = (TTree*)input->Get("TreeB2");
    gROOT->cd( outfileName+TString(":/") );
    dataloader->AddTree    (signalTree,"Signal");
    dataloader->AddTree    (background0,"bg0");
    dataloader->AddTree    (background1,"bg1");
    dataloader->AddTree    (background2,"bg2");
    
    dataloader->PrepareTrainingAndTestTree( "", "SplitMode=Random:NormMode=NumEvents:!V" );
    
    if (Use["BDTG"]) // gradient boosted decision trees
        factory->BookMethod( dataloader,  TMVA::Types::kBDT, "BDTG", "!H:!V:NTrees=400::BoostType=Grad:Shrinkage=0.1:nCuts=1000:MaxDepth=5:MinNodeSize=0.000001:"+
			  "RegressionLossFunctionBDTG=LeastSquares");
    if (Use["MLP"]) // neural network
        factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=tanh:NCycles=1000:HiddenLayers=N+5,5:TestRate=5:EstimatorType=MSE");
    if (Use["DNN"]) {
        TString layoutString ("Layout=TANH|100,TANH|50,TANH|10,LINEAR");
        TString training0 ("LearningRate=1e-1, Momentum=0.5, Repetitions=1, ConvergenceSteps=10,"
                           " BatchSize=256, TestRepetitions=10, Multithreading=True");
        TString training1 ("LearningRate=1e-2, Momentum=0.0, Repetitions=1, ConvergenceSteps=10,"
                           " BatchSize=256, TestRepetitions=7, Multithreading=True");
        TString trainingStrategyString ("TrainingStrategy=");
        trainingStrategyString += training0 + "|" + training1;
        TString nnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                           "WeightInitialization=XAVIERUNIFORM:Architecture=STANDARD");
        nnOptions.Append (":"); nnOptions.Append (layoutString);
        nnOptions.Append (":"); nnOptions.Append (trainingStrategyString);
        factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN", nnOptions );
    }
    
    // Train MVAs using the set of training events
    factory->TrainAllMethods();
    // Evaluate all MVAs using the set of test events
    factory->TestAllMethods();
    // Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();
    
    // Save the output
    outputFile->Close();
    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;
    delete factory;
    delete dataloader;
    
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVAMultiClassGui( outfileName );
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
