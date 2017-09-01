
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

      const int MODE = 14;  // Specify one mode in particular to look at

      IN_DIR_NAME   = "/afs/cern.ch/work/a/abrinke1/public/EMTF/PtAssign2017/files";
      TString out_str;
      out_str.Form("RateVsEff_mode_%d_eta_1p2_2p5", MODE);
      OUT_FILE_NAME = out_str;

      EFF_CUTS    = {90};
      TURN_ONS    = {8, 16, 24};

      TString in_str;
      TString fact_str;      
      TString ID_str;
      TString alias_str;

      PtAlgo EMTF15;  // 2016 EMTF pT algorithm, mode 15
      EMTF15.in_file_name = "PtRegression_Apr_2017_05_06_MODE_0_noBitCompr_noRPC.root";
      EMTF15.fact_name    = "f_MODE_0_logPtTarg_noWgt_noBitCompr_noRPC";
      EMTF15.MVA_name     = "EMTF_pt";
      EMTF15.unique_ID    = "EMTF15";
      EMTF15.alias        = "EMTF mode 15";
      EMTF15.modes        = {15};
      EMTF15.modes_CSC    = {15};
      EMTF15.modes_RPC    = {0};  // No RPC hits allowed
      EMTF15.trg_pt_scale = 1./1.4;  // Had been scaled up by 1.4 from original regression 
      EMTF15.color        = 1;  // kBlack 

      PtAlgo EMTF = EMTF15;
      ID_str.Form("EMTF%d", MODE);
      alias_str.Form("EMTF mode %d", MODE);
      EMTF.unique_ID    = ID_str;
      EMTF.alias        = alias_str;
      EMTF.modes        = {MODE};
      EMTF.modes_CSC    = {MODE};
      EMTF.color        = 880;  // kViolet



      // Mode 15, invPt pT target, invPt weight
      PtAlgo BDT15_invPt_invPt_Sq;
      BDT15_invPt_invPt_Sq.in_file_name = "PtRegression_Apr_2017_05_09_invPtTarg_invPtWgt_MODE_15_bitCompr_noRPC.root";
      BDT15_invPt_invPt_Sq.fact_name    = "f_MODE_15_invPtTarg_invPtWgt_bitCompr_noRPC";
      BDT15_invPt_invPt_Sq.MVA_name     = "BDTG_AWB";
      BDT15_invPt_invPt_Sq.unique_ID    = "BDT_15_invPt_Sq";
      BDT15_invPt_invPt_Sq.alias        = "invPt target, LeastSq loss";
      BDT15_invPt_invPt_Sq.modes        = {15};
      BDT15_invPt_invPt_Sq.modes_CSC    = {15};
      BDT15_invPt_invPt_Sq.modes_RPC    = {0};
      BDT15_invPt_invPt_Sq.color        = 840;  // kTeal


      // logPt target, invPt weight
      PtAlgo BDT_logPt_Abs;
      in_str.Form("PtRegression_Apr_2017_05_06_MODE_%d_bitCompr_RPC.root", MODE);
      fact_str.Form("f_MODE_%d_logPtTarg_invPtWgt_bitCompr_RPC", MODE);
      /* in_str.Form("PtRegression_Apr_2017_05_06_MODE_%d_bitCompr_noRPC.root", MODE); */
      /* fact_str.Form("f_MODE_%d_logPtTarg_invPtWgt_bitCompr_noRPC", MODE); */
      ID_str.Form("BDT%d_logPt_Abs", MODE);

      BDT_logPt_Abs.in_file_name = in_str;
      BDT_logPt_Abs.fact_name    = fact_str;
      BDT_logPt_Abs.MVA_name     = "BDTG_AWB";
      BDT_logPt_Abs.unique_ID    = ID_str;
      BDT_logPt_Abs.alias        = "logPt target, AbsDev loss";
      BDT_logPt_Abs.modes        = {MODE};
      BDT_logPt_Abs.modes_CSC    = {MODE};
      BDT_logPt_Abs.modes_RPC    = {0};
      BDT_logPt_Abs.color        = 4;  // kBlue


      // invPt target, invPt weight
      PtAlgo BDT_invPt_Abs;
      in_str.Form("PtRegression_Apr_2017_05_10_invPtTarg_invPtWgt_MODE_%d_bitCompr_RPC.root", MODE);
      fact_str.Form("f_MODE_%d_invPtTarg_invPtWgt_bitCompr_RPC", MODE);
      ID_str.Form("BDT%d_invPt_Abs", MODE);

      BDT_invPt_Abs.in_file_name = in_str;
      BDT_invPt_Abs.fact_name    = fact_str;
      BDT_invPt_Abs.MVA_name     = "BDTG_AWB";
      BDT_invPt_Abs.unique_ID    = ID_str;
      BDT_invPt_Abs.alias        = "invPt target, AbsDev loss";
      BDT_invPt_Abs.modes        = {MODE};
      BDT_invPt_Abs.modes_CSC    = {MODE};
      BDT_invPt_Abs.modes_RPC    = {0};
      BDT_invPt_Abs.color        = 2;  // kRed

      PtAlgo BDT_invPt_Hub = BDT_invPt_Abs;
      BDT_invPt_Hub.MVA_name     = "BDTG_AWB_Hub";
      BDT_invPt_Hub.unique_ID    = ID_str.ReplaceAll("_Abs", "_Hub");
      BDT_invPt_Hub.alias        = "invPt target, Huber loss";
      BDT_invPt_Hub.color        = 3;  // kGreen

      PtAlgo BDT_invPt_Sq = BDT_invPt_Abs;
      BDT_invPt_Sq.MVA_name     = "BDTG_AWB_Sq";
      BDT_invPt_Sq.unique_ID    = ID_str.ReplaceAll("_Hub", "_Sq");
      BDT_invPt_Sq.alias        = "invPt target, LeastSq loss, CSC-only";
      BDT_invPt_Sq.color        = 4;  // kBlue

      PtAlgo BDT_invPt_Sq_match = BDT_invPt_Sq;
      BDT_invPt_Sq_match.unique_ID  = ID_str.ReplaceAll("_Sq", "_Sq_match");
      BDT_invPt_Sq_match.alias      = "invPt target, LeastSq loss, CSC-only, matching same-mode EMTF track";
      BDT_invPt_Sq_match.color      = 2;  // kRed
      BDT_invPt_Sq_match.match_EMTF = 1;  // Require EMTF track in same event with same mode / CSC mode

      PtAlgo BDT_invPt_Sq_RPC = BDT_invPt_Sq;
      BDT_invPt_Sq_RPC.unique_ID    = ID_str.ReplaceAll("_Sq", "_Sq_RPC");
      BDT_invPt_Sq_RPC.alias        = "invPt target, LeastSq loss, CSC + 1 RPC";
      BDT_invPt_Sq_RPC.modes_CSC    = {MODE-1, MODE-2, MODE-4, MODE-8};
      BDT_invPt_Sq_RPC.modes_RPC    = {1, 2, 4, 8};
      BDT_invPt_Sq_RPC.color        = 2;  // kRed


      ALGOS.push_back(EMTF15);  // First algo is always the standard comparison algo
      ALGOS.push_back(EMTF);

      ALGOS.push_back(BDT15_invPt_invPt_Sq);

      ALGOS.push_back(BDT_logPt_Abs);
      ALGOS.push_back(BDT_invPt_Abs);
      ALGOS.push_back(BDT_invPt_Hub);
      ALGOS.push_back(BDT_invPt_Sq);
      ALGOS.push_back(BDT_invPt_Sq_match);
      ALGOS.push_back(BDT_invPt_Sq_RPC);

    } // End conditional: if (USER == "AWB")

    if (USER == "JTR") {
    } // End conditional: if (USER == "JTR")
    
  } // End function: inline void ConfigureUser()    

} // End namespace RateVsEff_cfg
