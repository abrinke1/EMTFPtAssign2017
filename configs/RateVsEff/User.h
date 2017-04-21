
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

      IN_DIR_NAME = "/afs/cern.ch/work/a/abrinke1/public/EMTF/PtAssign2017";
      EFF_CUTS    = {90};
      TURN_ONS    = {8, 16, 24};

      PtAlgo EMTF;  // 2016 EMTF pT algorithm
      EMTF.in_file_name = "PtRegression_Apr_2017_04_20_MODE_15_bitCompr.root";
      EMTF.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_bitCompr";
      EMTF.MVA_name     = "EMTF_pt";
      EMTF.modes        = {15};
      EMTF.modes_CSC    = {15};
      EMTF.modes_RPC    = {0};  // No RPC hits allowed
      EMTF.trg_pt_scale = 1./1.4;  // Had been scaled up by 1.4 from original regression 
      EMTF.color        = 1;  // kBlack 

      PtAlgo MVA1;
      MVA1.in_file_name = "PtRegression_Apr_2017_04_20_MODE_15_noBitCompr.root";
      MVA1.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_noBitCompr";
      MVA1.MVA_name     = "BDTG_AWB";
      MVA1.modes        = {15};
      MVA1.modes_CSC    = {15};
      MVA1.modes_RPC    = {0};  // No RPC hits allowed
      MVA1.color        = 4;  // kBlue 

      PtAlgo MVA2;
      MVA2.in_file_name = "PtRegression_Apr_2017_04_20_MODE_15_bitCompr.root";
      MVA2.fact_name    = "f_MODE_15_logPtTarg_invPtWgt_bitCompr";
      MVA2.MVA_name     = "BDTG_AWB";
      MVA2.modes        = {15};
      MVA2.modes_CSC    = {15};
      MVA2.modes_RPC    = {0};  // No RPC hits allowed
      MVA2.color        = 2;  // kRed 

      ALGOS.push_back(EMTF);  // First algo is always the standard comparison algo
      ALGOS.push_back(MVA1);
      ALGOS.push_back(MVA2);

    } // End conditional: if (USER == "AWB")

    if (USER == "JTR") {
    } // End conditional: if (USER == "JTR")
    
  } // End function: inline void ConfigureUser()    

} // End namespace RateVsEff_cfg
