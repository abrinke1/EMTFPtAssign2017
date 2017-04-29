
// *** Default user settings *** //
TString IN_DIR_NAME   = ".";          // Directory for input ROOT files
TString OUT_DIR_NAME  = "plots";      // Directory for output ROOT file
TString OUT_FILE_NAME = "RateVsEff";  // Name base for output ROOT file
std::vector<PtAlgo> ALGOS    = {};    // Vector of factory-MVA-mode sets for evaluation
std::vector<int>    EFF_CUTS = {};    // Vector of efficiency thresholds (%)
std::vector<int>    TURN_ONS = {};    // Vector of pT cuts for turn-on curves

namespace RateVsEff_cfg {

  inline void ConfigureUser( const TString USER ) {
    
    std::cout << "\nConfiguring RateVsEff code for user " << USER << std::endl;
    
    if (USER == "AWB") {

      IN_DIR_NAME = "/afs/cern.ch/work/a/abrinke1/public/EMTF/PtAssign2017/files/";
      EFF_CUTS    = {90};
      TURN_ONS    = {8, 16, 24};

      PtAlgo EMTF15;  // 2016 EMTF pT algorithm, mode 15
      EMTF15.in_file_name = "PtRegression_Apr_2017_04_20_MODE_15_bitCompr.root";
      EMTF15.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_bitCompr";
      EMTF15.MVA_name     = "EMTF_pt";
      EMTF15.unique_ID    = "EMTF_mode_15_15_0";
      EMTF15.alias        = "Mode 15";
      EMTF15.modes        = {15};
      EMTF15.modes_CSC    = {15};
      EMTF15.modes_RPC    = {0};  // No RPC hits allowed
      EMTF15.trg_pt_scale = 1./1.4;  // Had been scaled up by 1.4 from original regression 
      EMTF15.color        = 1;  // kBlack 

      PtAlgo MVA1;
      MVA1.in_file_name = "PtRegression_Apr_2017_04_20_MODE_15_noBitCompr.root";
      MVA1.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_noBitCompr";
      MVA1.MVA_name     = "BDTG_AWB";
      MVA1.unique_ID    = "BDT_CSC_only_uncompr_mode_15_15_0";
      MVA1.alias        = "Mode 15 (15/0) CSC-only BDT, uncompressed";
      MVA1.modes        = {15};
      MVA1.modes_CSC    = {15};
      MVA1.modes_RPC    = {0};  // No RPC hits allowed
      MVA1.color        = 4;  // kBlue 

      PtAlgo MVA2;
      MVA2.in_file_name = "PtRegression_Apr_2017_04_20_MODE_15_bitCompr.root";
      MVA2.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_bitCompr";
      MVA2.MVA_name     = "BDTG_AWB";
      MVA2.unique_ID    = "BDT_CSC_only_compr_mode_15_15_0";
      MVA2.alias        = "Mode 15 (15/0) CSC-only BDT, compressed";
      MVA2.modes        = {15};
      MVA2.modes_CSC    = {15};
      MVA2.modes_RPC    = {0};  // No RPC hits allowed
      MVA2.color        = 3;  // kGreen 

      PtAlgo MVA3;
      MVA3.in_file_name = "PtRegression_Apr_2017_04_28_MODE_15_bitCompr_RPC.root";
      MVA3.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_bitCompr_RPC";
      MVA3.MVA_name     = "BDTG_AWB";
      MVA3.unique_ID    = "BDT_CSC_RPC_compr_mode_15_15_0";
      MVA3.alias        = "Mode 15 (15/0) CSC+RPC BDT, compressed";
      MVA3.modes        = {15};
      MVA3.modes_CSC    = {15};
      MVA3.modes_RPC    = {0};  // No RPC hits allowed
      MVA3.color        = 2;  // kRed 

      PtAlgo MVA4;
      MVA4.in_file_name = "PtRegression_Apr_2017_04_28_MODE_15_bitCompr_RPC.root";
      MVA4.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_bitCompr_RPC";
      MVA4.MVA_name     = "BDTG_AWB";
      MVA4.unique_ID    = "BDT_CSC_RPC_compr_mode_15_14_1";
      MVA4.alias        = "Mode 15 (14/1) CSC+RPC BDT, compressed";
      MVA4.modes        = {15};
      MVA4.modes_CSC    = {14};
      MVA4.modes_RPC    = {1};
      MVA4.color        = 2;  // kRed 

      PtAlgo MVA5;
      MVA5.in_file_name = "PtRegression_Apr_2017_04_28_MODE_15_bitCompr_RPC.root";
      MVA5.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_bitCompr_RPC";
      MVA5.MVA_name     = "BDTG_AWB";
      MVA5.unique_ID    = "BDT_CSC_RPC_compr_mode_15_13_2";
      MVA5.alias        = "Mode 15 (13/2) CSC+RPC BDT, compressed";
      MVA5.modes        = {15};
      MVA5.modes_CSC    = {13};
      MVA5.modes_RPC    = {2};
      MVA5.color        = 2;  // kRed 

      PtAlgo MVA6;
      MVA6.in_file_name = "PtRegression_Apr_2017_04_28_MODE_15_bitCompr_RPC.root";
      MVA6.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_bitCompr_RPC";
      MVA6.MVA_name     = "BDTG_AWB";
      MVA6.unique_ID    = "BDT_CSC_RPC_compr_mode_15_11_4";
      MVA6.alias        = "Mode 15 (11/4) CSC+RPC BDT, compressed";
      MVA6.modes        = {15};
      MVA6.modes_CSC    = {11};
      MVA6.modes_RPC    = {4};
      MVA6.color        = 2;  // kRed 

      PtAlgo MVA7;
      MVA7.in_file_name = "PtRegression_Apr_2017_04_28_MODE_15_bitCompr_RPC.root";
      MVA7.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_bitCompr_RPC";
      MVA7.MVA_name     = "BDTG_AWB";
      MVA7.unique_ID    = "BDT_CSC_RPC_compr_mode_15_7_8";
      MVA7.alias        = "Mode 15 (7/8) CSC+RPC BDT, compressed";
      MVA7.modes        = {15};
      MVA7.modes_CSC    = {7};
      MVA7.modes_RPC    = {8};
      MVA7.color        = 2;  // kRed 

      ALGOS.push_back(EMTF15);  // First algo is always the standard comparison algo
      ALGOS.push_back(MVA1);
      ALGOS.push_back(MVA2);
      /* ALGOS.push_back(MVA3); */
      /* ALGOS.push_back(MVA4); */
      /* ALGOS.push_back(MVA5); */
      /* ALGOS.push_back(MVA6); */
      /* ALGOS.push_back(MVA7); */

    } // End conditional: if (USER == "AWB")

    if (USER == "JTR") {
    } // End conditional: if (USER == "JTR")
    
  } // End function: inline void ConfigureUser()    

} // End namespace RateVsEff_cfg
