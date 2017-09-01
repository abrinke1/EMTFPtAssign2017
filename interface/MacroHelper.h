
Float_t GetMedian(const TH1D* hist);

void WeightByResScore(TH1D& hist, const Float_t med_ratio);

Float_t GetResScore(const TH1D* hist, const Float_t med_ratio);

Float_t GetResScoreErr(const TH1D* hist, const Float_t med_ratio);


// Defines single input variable to MVA in TMVA macro
class MVA_var {

 public:

  // Default constructor
  MVA_var( TString _name = "", TString _descr = "", TString _unit = "", 
	   char _type = 'F', Double_t _def_val = -99 ) {
    name    = _name;
    descr   = _descr;
    unit    = _unit;
    type    = _type;
    def_val = _def_val;
  } // End default constructor MVA_var()

  TString name;
  TString descr;
  TString unit;
  char type;
  Double_t def_val;
  
}; // End class MVA_var


// Defines factory-MVA-mode for performance evaluation
class PtAlgo {

 public:
  
  // Default constructor
  PtAlgo() {

    in_file_name = "";
    
    TFile*  file_tmp(0);
    TChain* chain_tmp_1(0);
    TChain* chain_tmp_2(0);

    in_file      = file_tmp;
    train_tree   = chain_tmp_1;
    test_tree    = chain_tmp_2;

    fact_name = "";
    MVA_name  = "";
    unique_ID = "";
    alias     = "";
    
    modes      = {};
    modes_CSC  = {};
    modes_RPC  = {};
    match_EMTF = 0;
    
    trg_pt_scale = 1.0;
    color        = 1;
    
    h_ZB_rates = {};
    h_turn_ons = {};
    h_MC_etas  = {};
    h_ZB_etas  = {};

  } // End default constructor PtAlgo()
  
  // Input data, from TMVA training output
  TString in_file_name;
  TFile*  in_file;
  TChain* train_tree;
  TChain* test_tree;
  
  // Algorithm to test
  TString fact_name;
  TString MVA_name;
  TString unique_ID;
  TString alias;
  float   MVA_val;
  
  // Track modes to test
  std::vector<int> modes;
  std::vector<int> modes_CSC;
  std::vector<int> modes_RPC;
  bool match_EMTF; // Require that mode and mode_CSC match the EMTF track in the event
  
  double trg_pt_scale;  // Uniform initial scaling of trigger pT   
  int color;            // 1D histogram color

  // Train / test pairs of histograms with trigger and GEN pT, in MC and ZeroBias
  // "Counts" hisotgrams are normalized to unit area, "rate" histograms normalized to 1 in first bin

  std::pair<TH2D*, TH2D*> h_trg_vs_GEN_pt;      // 2D counts vs. trigger pT vs. GEN pT
  std::pair<TH2D*, TH2D*> h_trg_eff_vs_GEN_pt;  // 2D efficiency vs. trigger pT vs. GEN pT
  
  TH1D*                                  h_ZB_count;      // 1D ZeroBias counts vs. trigger pT
  TH2D*                                  h_ZB_count_eta;  // 2D ZeroBias counts vs. trigger pT and eta
  std::vector< std::pair<TH1D*, TH1D*> > h_ZB_rates;      // 1D ZeroBias rate vs. trigger pT for multiple efficiency thresholds
  std::vector< std::pair<TH1D*, TH1D*> > h_ZB_rates_eta;  // 2D ZeroBias rate vs. trigger pT and eta for multiple efficiency thresholds
  std::vector< std::pair<TH1D*, TH1D*> > h_pt_scales;     // 1D scaling factor vs. initial trigger pT needed to achieve XX% efficiency
  std::vector< std::pair<TH1D*, TH1D*> > h_turn_ons;      // 1D turn-on efficiency plots for multiple pT thresholds
  std::vector< std::pair<TH1D*, TH1D*> > h_MC_etas;       // 1D MC rate vs. muon eta for multiple pT thresholds
  std::vector< std::pair<TH1D*, TH1D*> > h_ZB_etas;       // 1D ZeroBias rate vs. muon eta for multiple pT thresholds
  TH1D*                                  h_charge_eff;    // 1D EMTF correct-charge efficiency vs. pT 
};
