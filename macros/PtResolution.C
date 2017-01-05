
/////////////////////////////////////////////////////////
///          Macro to plot pT resolution and          ///
///             compute resolution "score"            ///   
///            Andrew Brinkerhoff 10.12.16            ///
///                                                   ///
/// * Resolution plot is log2(TRG pT / GEN pT)        ///
/// * Resolution score is the average event weight:   ///
///   MAX[ (TRG pT / GEN pT), (TRG pT / GEN pT) ]     ///
/// * In the score calculation, TRG pT is scaled so   ///
///   that the median (TRG pT / GEN pT) value = 1.0   ///
/////////////////////////////////////////////////////////

#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

#include <iomanip>  // std::cout formatting

#include "../interface/PtResolution.h"  // Function declarations
#include "../src/MacroHelper.C"             // Helpful common functions (GetMedian, GetResScore, etc.)


const int    PRTEVT  =      10000;  // When processing file, print every X events
const int    MAXEVT  =         -1;  // Maximum number of events to process
const int     NBINS  =       1600;  // Number of bins in log2(TRG pT / GEN pT) plot
const int     XMIN   =         -8;  // Minimum log2(TRG pT / GEN pT) value in plot
const int     XMAX   =          8;  // Maximum log2(TRG pT / GEN pT) value in plot  
const int     REBIN  =         10;  // Rebin log2(TRG pT / GEN pT) by X for display purposes
const double  PTMIN  =         1.;  // Minimum GEN and TRG pT values used
const double  PTMAX  =       256.;  // Maximum GEN and TRG pT values used
const TString TRGPT  =      "inv";  // Target trigger pT: "inv" for 1/pT, "log2" for log2(pT), "" for pT

const double  EMTFSC =        1.0;  // Scaling of EMTF pT to yield median ratio of 1.0 (optional)
const double  BDTGSC =        1.0;  // Scaling of BDTG pT to yield median ratio of 1.0 (optional)
// const double  EMTFSC = (1./1.245);  // Scaling of EMTF pT to yield median ratio of 1.0 (optional)
// const double  BDTGSC = (1./1.020);  // Scaling of BDTG pT to yield median ratio of 1.0 (optional)

// const TString EVTWGT = "invLog10";  // Weight events: "invLog10" --> 1/(1 + log10(pT))
// const TString EVTWGT = "invLog2";   // Weight events: "invLog2" --> 1/(1 + log2(pT))
const TString EVTWGT =  "invSqrt";  // Weight events: "invSqrt"  --> 1/(sqrt(pT))
// const TString EVTWGT =    "invPt";  // Weight events: "invPt"  --> 1/pT
// const TString EVTWGT =         "";  // Don't weight events


///////////////////////////////////////
///  Main function: PtResolution()  /// 
///////////////////////////////////////

void PtResolution() {

  // Initialize empty file to access each file in the list
  TFile* file_tmp(0);

  // List of input files
  std::vector<TString> in_file_names;
  // in_file_names.push_back("/afs/cern.ch/user/a/abrinke1/TMVA/EMTFPtAssign2017/PtRegression_AWB_v0_17_01_03_no_wgt_inv_pt_orig_vars.root");
  in_file_names.push_back("/afs/cern.ch/user/a/abrinke1/TMVA/EMTFPtAssign2017/PtRegression_AWB_v0_17_01_05_no_wgt_inv_pt_orig_deriv_vars_10k.root");

  // Open all input files
  for (int i = 0; i < in_file_names.size(); i++) {
    if ( !gSystem->AccessPathName(in_file_names.at(i)) )
      file_tmp = TFile::Open( in_file_names.at(i) ); // Check if file exists
    if (!file_tmp) {
      std::cout << "ERROR: could not open data file " << in_file_names.at(i) << std::endl;
      return;
    }
  }

  // Add trees from the input files to the TChain
  TChain* train_chain = new TChain("dataset/TrainTree");
  TChain* test_chain = new TChain("dataset/TestTree");
  for (int i = 0; i < in_file_names.size(); i++) {
    train_chain->Add( in_file_names.at(i) );
    test_chain->Add( in_file_names.at(i) );
  }

  TString out_file_name = "plots/PtResolution_AWB_v0_17_01_05.root";
  TFile *out_file = new TFile(out_file_name, "recreate");


  //////////////////////////////////////////////////
  // Set up histograms to fill and scores to compute
  //////////////////////////////////////////////////

  // Tuple for trigger pT values: name, branch name, branch value, type of value, median ratio to GEN pT
  std::vector<std::tuple<TString, TString, float, TString, float>> pt_br;
  pt_br.push_back( std::make_tuple("EMTF_train",                   "EMTF_pt",                -99.,    "", 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_default_train",           "BDTG_default",           -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_default_20k_trees_train", "BDTG_default_20k_trees", -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_default_5_deep_train",    "BDTG_default_5_deep",    -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_LeastSq_train",           "BDTG_LeastSq",           -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_LeastSq_20k_trees_train", "BDTG_LeastSq_20k_trees", -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_LeastSq_5_deep_train",    "BDTG_LeastSq_5_deep",    -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_Carnes_Huber_train",      "BDTG_Carnes_Huber",      -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_Carnes_AbsDev_train",     "BDTG_Carnes_AbsDev" ,    -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_Carnes_LeastSq_train",    "BDTG_Carnes_LeastSq",    -99., TRGPT, 1.0) );

  pt_br.push_back( std::make_tuple("EMTF_test",                   "EMTF_pt",                -99.,    "", 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_default_test",           "BDTG_default",           -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_default_20k_trees_test", "BDTG_default_20k_trees", -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_default_5_deep_test",    "BDTG_default_5_deep",    -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_LeastSq_test",           "BDTG_LeastSq",           -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_LeastSq_20k_trees_test", "BDTG_LeastSq_20k_trees", -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_LeastSq_5_deep_test",    "BDTG_LeastSq_5_deep",    -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_Carnes_Huber_test",      "BDTG_Carnes_Huber",      -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_Carnes_AbsDev_test",     "BDTG_Carnes_AbsDev" ,    -99., TRGPT, 1.0) );
  pt_br.push_back( std::make_tuple("BDTG_Carnes_LeastSq_test",    "BDTG_Carnes_LeastSq",    -99., TRGPT, 1.0) );


  // Set bins in GEN pT and eta
  std::vector<std::tuple<TString, float, float, TString>> pt_bins;
  std::vector<std::tuple<TString, float, float, TString>> eta_bins;

  pt_bins.push_back( std::make_tuple("all",      1.,   1000., "    ALL    ") );
  pt_bins.push_back( std::make_tuple("1_8",      1.,      8., "[  1,    8]") );
  pt_bins.push_back( std::make_tuple("8_30",     8.,     30., "[  8,   30]") );
  // pt_bins.push_back( std::make_tuple("1_4",      1.,      4., "[  1,    4]") );
  // pt_bins.push_back( std::make_tuple("4_8",      4.,      8., "[  4,    8]") );
  // pt_bins.push_back( std::make_tuple("8_15",     8.,     15., "[  8,   15]") );
  // pt_bins.push_back( std::make_tuple("15_30",   15.,     30., "[ 15,   30]") );
  pt_bins.push_back( std::make_tuple("30_120",   30.,   120., "[ 30,  120]") );
  pt_bins.push_back( std::make_tuple("120_1000", 120., 1000., "[120, 1000]") );

  eta_bins.push_back( std::make_tuple("all",       1.2,  2.4 , "     ALL    ") );
  eta_bins.push_back( std::make_tuple("1p2_1p55",  1.2,  1.55, "[1.2 , 1.55]") );
  eta_bins.push_back( std::make_tuple("1p55_1p85", 1.55, 1.85, "[1.55, 1.85]") );
  eta_bins.push_back( std::make_tuple("1p85_2p1",  1.85, 2.1 , "[1.85, 2.1 ]") );
  eta_bins.push_back( std::make_tuple("2p1_2p4",   2.1,  2.4 , "[2.1 , 2.4 ]") );

  // 1D resolution plots
  std::vector< std::vector< std::vector<TH1D*> > > h_res;
  std::vector< std::vector<TH1D*> > empty1a;
  std::vector<TH1D*> empty1b;

  // 2D resolution score, vs. pT / eta
  std::vector<TH2D*> h_score;

  // Book 1D pT resolution histograms for each trigger and pT/eta bin
  for (int iBr = 0; iBr < pt_br.size(); iBr++) {
    h_res.push_back( empty1a );
    TString br_name = std::get<0>(pt_br.at(iBr));
    BookScoreHist( pt_br, iBr, pt_bins, eta_bins, br_name, h_score );

    for (int iPt = 0; iPt < pt_bins.size(); iPt++) {
      h_res.at(iBr).push_back( empty1b );

      for (int iEta = 0; iEta < eta_bins.size(); iEta++) {
	BookResHist( pt_br, iBr, pt_bins, iPt, eta_bins, iEta, br_name, h_res );
      }
    }
  }

  /////////////////////////////////////////////
  // Loop over events from train and test trees
  /////////////////////////////////////////////
  LoopOverEvents( train_chain, "train", pt_bins, eta_bins, pt_br, h_res);
  LoopOverEvents( test_chain,  "test",  pt_bins, eta_bins, pt_br, h_res);

  //////////////////////////////////////////
  // Write out histograms and compute scores
  //////////////////////////////////////////

  out_file->cd();

  // Compute median ratios to GEN values, store in pt_br
  StoreMedians( pt_bins, eta_bins, h_res, pt_br );

  // Initialize vectors and of scores
  int num_pt_bins = pt_bins.size();
  int num_trg     = (pt_br.size() / 2) - 1;  // Only non-EMTF branches; divide by 2 for train/test

  std::vector<float> EMTF_vec_train     ; for (int i = 0; i < num_pt_bins; i++) EMTF_vec_train.    push_back(-99);
  std::vector<float> EMTF_vec_train_err ; for (int i = 0; i < num_pt_bins; i++) EMTF_vec_train_err.push_back(-99);
  std::vector<float> EMTF_vec_test      ; for (int i = 0; i < num_pt_bins; i++) EMTF_vec_test.     push_back(-99);
  std::vector<float> EMTF_vec_test_err  ; for (int i = 0; i < num_pt_bins; i++) EMTF_vec_test_err. push_back( -99);
  std::vector<float> ratio_vec_train    ; for (int i = 0; i < num_pt_bins * num_trg; i++) ratio_vec_train    .push_back(-99);
  std::vector<float> ratio_vec_train_err; for (int i = 0; i < num_pt_bins * num_trg; i++) ratio_vec_train_err.push_back(-99);
  std::vector<float> ratio_vec_test     ; for (int i = 0; i < num_pt_bins * num_trg; i++) ratio_vec_test     .push_back(-99);
  std::vector<float> ratio_vec_test_err ; for (int i = 0; i < num_pt_bins * num_trg; i++) ratio_vec_test_err .push_back(-99);

  // Write ratio histograms and print out scores
  FillScores( pt_br, pt_bins, eta_bins, num_pt_bins, num_trg, 
	      EMTF_vec_train,  EMTF_vec_train_err,  EMTF_vec_test,  EMTF_vec_test_err,
	      ratio_vec_train, ratio_vec_train_err, ratio_vec_test, ratio_vec_test_err,
	      h_res, h_score );

  // Fill 2D graphs of resolution ratio to EMTF
  FillRatioGraphs( pt_br, pt_bins, eta_bins, num_pt_bins, num_trg,
		   EMTF_vec_train,  EMTF_vec_train_err,  EMTF_vec_test,  EMTF_vec_test_err,
		   ratio_vec_train, ratio_vec_train_err, ratio_vec_test, ratio_vec_test_err );

  // Print out ratios of scores from different triggers
  PrintRatios( pt_br, pt_bins, eta_bins, h_res );
  
  out_file->Close();

  std::cout << "\nExiting PtResolution()\n";

} // End void PtResolution()


//////////////////////////////////////////
///  Functions used in PtResolution()  /// 
//////////////////////////////////////////

void BookScoreHist( const std::vector<std::tuple<TString, TString, float, TString, float>> pt_br, const int iBr,
		    const std::vector<std::tuple<TString, float, float, TString>> pt_bins, 
		    const std::vector<std::tuple<TString, float, float, TString>> eta_bins, 
		    TString& br_name, std::vector<TH2D*>& h_score ) {

  TString h_score_str = "h_score_"+std::get<0>(pt_br.at(iBr));
  h_score.push_back( new TH2D(h_score_str, h_score_str, pt_bins.size(), 0, pt_bins.size(), eta_bins.size(), 0, eta_bins.size()) );
  h_score.at(iBr)->Sumw2();
  br_name.ReplaceAll("_", " ");
  br_name.ReplaceAll("train", "(train)");
  br_name.ReplaceAll("test",  "(test)");
  h_score.at(iBr)->SetTitle(br_name+" p_{T} resolution score");
  h_score.at(iBr)->GetXaxis()->SetTitle("p_{T} range (GeV)");
  h_score.at(iBr)->GetYaxis()->SetTitle("|#eta| range");
  
} // End void BookScoreHist()

void BookResHist( const std::vector<std::tuple<TString, TString, float, TString, float>> pt_br, const int iBr,
		  const std::vector<std::tuple<TString, float, float, TString>> pt_bins, const int iPt,
		  const std::vector<std::tuple<TString, float, float, TString>> eta_bins, const int iEta,
		  TString& br_name, std::vector< std::vector< std::vector<TH1D*> > >& h_res ) {
  
  TString h_res_str = "h_res_"+std::get<0>(pt_br.at(iBr))+"_pt_"+std::get<0>(pt_bins.at(iPt))+"_eta_"+std::get<0>(eta_bins.at(iEta));
  h_res.at(iBr).at(iPt).push_back( new TH1D(h_res_str, h_res_str, NBINS, XMIN, XMAX) );
  h_res.at(iBr).at(iPt).at(iEta)->Sumw2();
  
  std::ostringstream ss_title;
  ss_title << br_name << " p_{T} resolution, "
	   << std::get<1>(pt_bins.at(iPt))   << " < p_{T} < "  << std::get<2>(pt_bins.at(iPt)) << " GeV, "
	   << std::get<1>(eta_bins.at(iEta)) << " < |#eta| < " << std::get<2>(eta_bins.at(iEta));
  h_res.at(iBr).at(iPt).at(iEta)->SetTitle( (TString) ss_title.str() );
  h_res.at(iBr).at(iPt).at(iEta)->GetXaxis()->SetTitle("log_{2} trigger p_{T} / GEN p_{T}");
  // std::cout << h_res.at(iBr).at(iPt).at(iEta)->GetName() << std::endl;	
  
} // End void BookResHist()

void LoopOverEvents( TChain* chain, const TString opt_str, 
		     const std::vector<std::tuple<TString, float, float, TString>> pt_bins, 
		     const std::vector<std::tuple<TString, float, float, TString>> eta_bins, 
		     std::vector<std::tuple<TString, TString, float, TString, float>>& pt_br, 
		     std::vector< std::vector< std::vector<TH1D*> > >& h_res ) {
  
  // Get GEN branches from the chains
  float GEN_pt_br;
  float GEN_eta_br;
  chain->SetBranchAddress("GEN_pt", &GEN_pt_br);
  chain->SetBranchAddress("GEN_eta", &GEN_eta_br);

  // Get trigger branches from the chains
  for (int iBr = 0; iBr < pt_br.size(); iBr++) {
    if ( std::get<0>(pt_br.at(iBr)).Contains("train") && opt_str.Contains("train") )
      chain->SetBranchAddress( std::get<1>(pt_br.at(iBr)), &(std::get<2>(pt_br.at(iBr))) );
    if ( std::get<0>(pt_br.at(iBr)).Contains("test")  && opt_str.Contains("test") )
      chain->SetBranchAddress( std::get<1>(pt_br.at(iBr)), &(std::get<2>(pt_br.at(iBr))) );
  }

  std::cout << "\n******* About to enter the " << opt_str << " event loop *******" << std::endl;
  for (int iEvt = 0; iEvt < chain->GetEntries(); iEvt++) {
    
    if (iEvt > MAXEVT && MAXEVT > 0) break;
    if ( (iEvt % PRTEVT) == 0 ) std::cout << "*** Looking at event " << iEvt << " ***" << std::endl;
    
    chain->GetEntry(iEvt);

    // Loop over different trigger pT computations
    for (int iBr = 0; iBr < pt_br.size(); iBr++) {
      if ( !(std::get<0>(pt_br.at(iBr)).Contains(opt_str)) ) continue;

      // Access trigger pT
      double TRG_pt = std::get<2>(pt_br.at(iBr));
      if ( std::get<3>(pt_br.at(iBr)) == "inv" )
	TRG_pt = 1. / max(0.001, TRG_pt); // Protect against negative 1/pT values
      if ( std::get<3>(pt_br.at(iBr)) == "log2" )
	TRG_pt = pow(2, TRG_pt);
      double GEN_pt = double(GEN_pt_br);

      // Scale trigger pT to shift median TRG pT / GEN pT value to 1.0 (optional)
      if (std::get<0>(pt_br.at(iBr)).Contains("EMTF")) TRG_pt *= EMTFSC;
      if (std::get<0>(pt_br.at(iBr)).Contains("BDTG")) TRG_pt *= BDTGSC;

      // Loop over pT and eta bins
      for (int iPt = 0; iPt < pt_bins.size(); iPt++) {
	if ( GEN_pt_br < std::get<1>(pt_bins.at(iPt)) || GEN_pt_br > std::get<2>(pt_bins.at(iPt)) ) continue;
	for (int iEta = 0; iEta < eta_bins.size(); iEta++) {
	  if ( fabs(GEN_eta_br) < std::get<1>(eta_bins.at(iEta)) || fabs(GEN_eta_br) > std::get<2>(eta_bins.at(iEta)) ) continue;
	  
	  // Enforce pT range, fill resolution hisogram
	  double evt_wgt = 1.0;
	  if (EVTWGT == "invLog10") evt_wgt = 1. / (1. + log10(GEN_pt));
	  if (EVTWGT == "invLog2")  evt_wgt = 1. / (1. + log2(GEN_pt));
	  if (EVTWGT == "invSqrt")  evt_wgt = 1. / sqrt(GEN_pt);
	  if (EVTWGT == "invPt")    evt_wgt = 1. / GEN_pt;

	  TRG_pt = min( PTMAX, max( PTMIN, TRG_pt ) );
	  GEN_pt = min( PTMAX, max( PTMIN, GEN_pt ) );

	  h_res.at(iBr).at(iPt).at(iEta)->Fill( log2( TRG_pt / GEN_pt ), evt_wgt );
	  
	} // End loop: for (int iEta = 0; iEta < eta_bins.size(); iEta++)
      } // End loop: for (int iPt = 0; iPt < pt_bins.size(); iPt++)
    } // End loop: for (int iBr = 0; iBr < pt_br.size(); iBr++)

  } // End loop: for (int iEvt = 0; iEvt < chain->GetEntries(); iEvt++)
  std::cout << "******* Leaving the " << opt_str << " event loop *******" << std::endl;

} // End void LoopOverEvents()

void StoreMedians( const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
		   const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
		   const std::vector< std::vector< std::vector<TH1D*> > > h_res,
		   std::vector<std::tuple<TString, TString, float, TString, float>>& pt_br ) {

  std::cout << "\n" << std::string(33,'*') << "\n" << "Median trigger pT / GEN pT values" << "\n" << std::string(33,'*') << std::endl;
  for (int iBr = 0; iBr < pt_br.size(); iBr++) {
    for (int iPt = 0; iPt < pt_bins.size(); iPt++) {
      if ( !(std::get<0>(pt_bins.at(iPt)).Contains("all")) ) continue;
      for (int iEta = 0; iEta < eta_bins.size(); iEta++) {
	if ( !(std::get<0>(eta_bins.at(iEta)).Contains("all")) ) continue;
	std::get<4>(pt_br.at(iBr)) = pow(2, GetMedian( h_res.at(iBr).at(iPt).at(iEta) ));  // Store median ratio in pt_br
	std::cout << std::setw(35) << std::left << std::get<0>(pt_br.at(iBr)) 
		  << std::fixed << std::setprecision(3) << std::get<4>(pt_br.at(iBr)) << std::endl;
      }
    }
  }

} // End void StoreMedians()

void FillScores( const std::vector<std::tuple<TString, TString, float, TString, float>> pt_br, 
		 const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
		 const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
		 const int num_pt_bins, const int num_trg,
		 std::vector<float>& EMTF_vec_train,  std::vector<float>& EMTF_vec_train_err,  
		 std::vector<float>& EMTF_vec_test,   std::vector<float>& EMTF_vec_test_err,
		 std::vector<float>& ratio_vec_train, std::vector<float>& ratio_vec_train_err, 
		 std::vector<float>& ratio_vec_test,  std::vector<float>& ratio_vec_test_err,
		 const std::vector< std::vector< std::vector<TH1D*> > > h_res,
		 std::vector<TH2D*>& h_score ) {

  for (int iBr = 0; iBr < pt_br.size(); iBr++) {
    
    // Print table header
    if ( std::get<0>(pt_br.at(iBr)).Contains("test") )
      std::cout << "\n" << std::string(45,'*') << "\n" << std::get<0>(pt_br.at(iBr)) 
		<< " resolution score" << "\n" << std::string(45,'*') << std::endl;
    
    for (int iPt = 0; iPt < pt_bins.size(); iPt++) {
      
      if ( std::get<0>(pt_br.at(iBr)).Contains("test") ) {
  	if (iPt == 0) {  // Print top row with eta bins
  	  std::cout << "              ";
  	  for (int iEta = 0; iEta < eta_bins.size(); iEta++) 
  	    std::cout << "  " << std::get<3>(eta_bins.at(iEta));
  	  std::cout << std::endl;
  	}
	// Print first column with pT bins
	std::cout << std::get<3>(pt_bins.at(iPt));
      }

      // Print pT resolution score values for each pT/eta bin
      for (int iEta = 0; iEta < eta_bins.size(); iEta++) {
	
  	if ( std::get<0>(pt_br.at(iBr)).Contains("test") ) {
	  // Print score
  	  std::cout << "          " << std::fixed << std::setprecision(2) << GetResScore( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
	  // Print score with uncertainty
  	  // std::cout << "   " << std::fixed << std::setprecision(2) << GetResScore( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) )
  	  // 	    << " +/- " << GetResScoreErr( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
	  // // Print median ratio to GEN
  	  // std::cout << "          " << std::fixed << std::setprecision(2) << GetMedian( h_res.at(iBr).at(iPt).at(iEta) );
	  // // Print number of events
  	  // std::cout << std::setw(14) << std::right << std::fixed << std::setprecision(0) << h_res.at(iBr).at(iPt).at(iEta)->Integral();
  	}
	
	// Write out resolution histogram
  	TH1D* h_tmp = (TH1D*) ( (TH1D*) h_res.at(iBr).at(iPt).at(iEta)->Clone() )->Rebin(REBIN);  // Don't rebin h_res
  	h_tmp->SetLineWidth(2);
	if ( std::get<0>(pt_br.at(iBr)).Contains("EMTF") )
	  h_tmp->SetLineColor(kBlack);
	else
	  h_tmp->SetLineColor(kRed);
  	h_tmp->Write();
	
  	TH1D* h_tmp_wgt = (TH1D*) h_res.at(iBr).at(iPt).at(iEta)->Clone();
	WeightByResScore( (*h_tmp_wgt), std::get<4>(pt_br.at(iBr)) );
	h_tmp_wgt = (TH1D*) h_tmp_wgt->Rebin(REBIN);
  	h_tmp_wgt->SetLineWidth(2);
	if ( std::get<0>(pt_br.at(iBr)).Contains("EMTF") )
	  h_tmp_wgt->SetLineColor(kBlack);
	else
	  h_tmp_wgt->SetLineColor(kRed);
	h_tmp_wgt->SetName( (TString) h_tmp->GetName()+"_wgt" );
	h_tmp_wgt->SetTitle( (TString) h_tmp->GetTitle()+" (weighted)");
  	h_tmp_wgt->Write();
	
	// Fill score histogram
	h_score.at(iBr)->SetBinContent(iPt+1, iEta+1, GetResScore(    h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) ) );
	h_score.at(iBr)->SetBinError  (iPt+1, iEta+1, GetResScoreErr( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) ) );
	h_score.at(iBr)->GetXaxis()->SetBinLabel( iPt+1,  std::get<3>(pt_bins.at(iPt))   );
	h_score.at(iBr)->GetYaxis()->SetBinLabel( iEta+1, std::get<3>(eta_bins.at(iEta)) );
	
	// Fill vectors of scores
	if ( std::get<0>(eta_bins.at(iEta)).Contains("all") ) {
	  if ( std::get<0>(pt_br.at(iBr)).Contains("train") ) {
	    if ( std::get<0>(pt_br.at(iBr)).Contains("EMTF") ) {
	      EMTF_vec_train    .at(iPt) = GetResScore   ( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
	      EMTF_vec_train_err.at(iPt) = GetResScoreErr( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
	    } else {
	      ratio_vec_train    .at((iBr - 1)*num_pt_bins + iPt) = GetResScore   ( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
	      ratio_vec_train_err.at((iBr - 1)*num_pt_bins + iPt) = GetResScoreErr( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
	    }
	  } else {
	      if ( std::get<0>(pt_br.at(iBr)).Contains("EMTF") ) {
		EMTF_vec_test    .at(iPt) = GetResScore   ( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
		EMTF_vec_test_err.at(iPt) = GetResScoreErr( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
	      } else {
		ratio_vec_test    .at((iBr - 2 - num_trg)*num_pt_bins + iPt) = GetResScore   ( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
		ratio_vec_test_err.at((iBr - 2 - num_trg)*num_pt_bins + iPt) = GetResScoreErr( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
	      }
	  }
	  }
	
      } // End loop: for (int iEta = 0; iEta < eta_bins.size(); iEta++)
      if ( std::get<0>(pt_br.at(iBr)).Contains("test") )
  	std::cout << std::endl;
    } // End loop: for (int iPt = 0; iPt < pt_bins.size(); iPt++)
    
    // Fill score histogram
    h_score.at(iBr)->Write();
    
  } // End loop: for (int iBr = 0; iBr < pt_br.size(); iBr++)

} // End void FillScores()

void FillRatioGraphs( const std::vector<std::tuple<TString, TString, float, TString, float>> pt_br, 
		      const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
		      const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
		      const int num_pt_bins, const int num_trg,
		      const std::vector<float> EMTF_vec_train,  const std::vector<float> EMTF_vec_train_err,  
		      const std::vector<float> EMTF_vec_test,   const std::vector<float> EMTF_vec_test_err,
		      const std::vector<float> ratio_vec_train, const std::vector<float> ratio_vec_train_err, 
		      const std::vector<float> ratio_vec_test,  const std::vector<float> ratio_vec_test_err ) {

  // Initialize arrays of scores
  float EMTF_arr_train     [num_pt_bins]; for (int i = 0; i < num_pt_bins; i++) EMTF_arr_train    [i] = EMTF_vec_train.at(i);
  float EMTF_arr_train_err [num_pt_bins]; for (int i = 0; i < num_pt_bins; i++) EMTF_arr_train_err[i] = EMTF_vec_train_err.at(i);
  float EMTF_arr_test      [num_pt_bins]; for (int i = 0; i < num_pt_bins; i++) EMTF_arr_test     [i] = EMTF_vec_test.at(i);
  float EMTF_arr_test_err  [num_pt_bins]; for (int i = 0; i < num_pt_bins; i++) EMTF_arr_test_err [i] = EMTF_vec_test_err.at(i);
  float ratio_arr_train    [num_pt_bins * num_trg]; for (int i = 0; i < num_pt_bins * num_trg; i++) ratio_arr_train    [i] = ratio_vec_train.at(i);
  float ratio_arr_train_err[num_pt_bins * num_trg]; for (int i = 0; i < num_pt_bins * num_trg; i++) ratio_arr_train_err[i] = ratio_vec_train_err.at(i);
  float ratio_arr_test     [num_pt_bins * num_trg]; for (int i = 0; i < num_pt_bins * num_trg; i++) ratio_arr_test     [i] = ratio_vec_test.at(i);
  float ratio_arr_test_err [num_pt_bins * num_trg]; for (int i = 0; i < num_pt_bins * num_trg; i++) ratio_arr_test_err [i] = ratio_vec_test_err.at(i);

  for (int iRat = 0; iRat < num_pt_bins * num_trg; iRat++) {
    int jRat = iRat % num_pt_bins;
    ratio_arr_train_err[iRat] = (ratio_arr_train[iRat] / EMTF_arr_train[jRat]) * sqrt( pow(ratio_arr_train_err[iRat] / ratio_arr_train[iRat], 2) +
										       pow(EMTF_arr_train_err [jRat] / EMTF_arr_train [jRat], 2) );
    ratio_arr_train    [iRat] =  ratio_arr_train[iRat] / EMTF_arr_train[jRat];
    ratio_arr_test_err [iRat] = (ratio_arr_test [iRat] / EMTF_arr_test [jRat]) * sqrt( pow(ratio_arr_test_err [iRat] / ratio_arr_test [iRat], 2) +
										       pow(EMTF_arr_test_err  [jRat] / EMTF_arr_test  [jRat], 2) );
    ratio_arr_test     [iRat] =  ratio_arr_test [iRat] / EMTF_arr_test [jRat];
  }

  float x_arr_train[num_pt_bins * num_trg];  
  float x_arr_test [num_pt_bins * num_trg];
  float x_arr_zero [num_pt_bins * num_trg];
  for (int i = 0; i < num_pt_bins * num_trg; i++) {
    x_arr_train[i] = i*1.0 - 0.1;
    x_arr_test [i] = i*1.0 + 0.1;
    x_arr_zero [i] = 0.0;
  }
  TGraphErrors* h_ratio_train = new TGraphErrors(num_pt_bins * num_trg, x_arr_train, ratio_arr_train, x_arr_zero, ratio_arr_train_err);
  TGraphErrors* h_ratio_test  = new TGraphErrors(num_pt_bins * num_trg, x_arr_test,  ratio_arr_test,  x_arr_zero, ratio_arr_test_err );

  h_ratio_train->SetName("h_ratio_train");
  h_ratio_train->SetTitle("Ratio of p_{T} resolution scores");
  h_ratio_train->GetYaxis()->SetTitle("score / EMTF score");
  h_ratio_train->SetMarkerSize(1);
  h_ratio_train->SetMarkerStyle(20);
  h_ratio_train->SetMarkerColor(1);  // kBlack
  h_ratio_train->SetLineWidth(2);
  h_ratio_train->SetLineColor(1);  // kBlack
  for (int i = 0; i < num_pt_bins * num_trg; i++) {
    float bin_val = h_ratio_train->GetXaxis()->FindBin(i*1.0);
    if ( (i % num_pt_bins) == 0 ) h_ratio_train->GetXaxis()->SetBinLabel( bin_val, std::get<1>(pt_br.at  ( (i / num_pt_bins) + 1 )) );
    else                          h_ratio_train->GetXaxis()->SetBinLabel( bin_val, std::get<3>(pt_bins.at( (i % num_pt_bins)     )) );
  }
  h_ratio_train->Write();

  h_ratio_test->SetName("h_ratio_test");
  h_ratio_test->SetTitle("Ratio of p_{T} resolution scores");
  h_ratio_test->GetYaxis()->SetTitle("score / EMTF score");
  h_ratio_test->SetMarkerSize(1);
  h_ratio_test->SetMarkerStyle(21);
  h_ratio_test->SetMarkerColor(2);  // kRed
  h_ratio_test->SetLineWidth(2);
  h_ratio_test->SetLineColor(2);  // kRed
  for (int i = 0; i < num_pt_bins * num_trg; i++) {
    float bin_val = h_ratio_test->GetXaxis()->FindBin(i*1.0);
    if ( (i % num_pt_bins) == 0 ) h_ratio_test->GetXaxis()->SetBinLabel( bin_val, std::get<1>(pt_br.at  ( (i / num_pt_bins) + 1 )) );
    else                          h_ratio_test->GetXaxis()->SetBinLabel( bin_val, std::get<3>(pt_bins.at( (i % num_pt_bins)     )) );
  }
  h_ratio_test->Write();

} // End FillRatioGraphs()

void PrintRatios( const std::vector<std::tuple<TString, TString, float, TString, float>> pt_br,
		  const std::vector<std::tuple<TString, float, float, TString>> pt_bins, 
		  const std::vector<std::tuple<TString, float, float, TString>> eta_bins, 
		  const std::vector< std::vector< std::vector<TH1D*> > > h_res ) {

  for (int iBr = 0; iBr < pt_br.size(); iBr++) {  // Denominator trigger
    for (int jBr = iBr+1; jBr < pt_br.size(); jBr++) {  // Numerator trigger
      // Don't print unnecessary tables
      if ( ( std::get<0>(pt_br.at(iBr)).Contains("test")  == std::get<0>(pt_br.at(jBr)).Contains("test") ) ==
  	   ( std::get<0>(pt_br.at(iBr)).Contains("EMTF")  == std::get<0>(pt_br.at(jBr)).Contains("EMTF") ) ) continue;
      if ( ( std::get<0>(pt_br.at(iBr)).Contains("test")  != std::get<0>(pt_br.at(jBr)).Contains("test") ) &&
  	   ( std::get<1>(pt_br.at(iBr))                   != std::get<1>(pt_br.at(jBr))                  ) ) continue;
      if (   std::get<0>(pt_br.at(iBr)).Contains("train") && std::get<0>(pt_br.at(jBr)).Contains("train")  ) continue;

      // Print table header
      std::cout << "\n" << std::string(70,'*') << "\n" << std::get<0>(pt_br.at(jBr)) << " / "
		<< std::get<0>(pt_br.at(iBr)) << " score ratio" << "\n" << std::string(70,'*') << std::endl;
      for (int iPt = 0; iPt < pt_bins.size(); iPt++) {
  	if (iPt == 0) {  // Print top row with eta bins
  	  std::cout << "            ";
  	  for (int iEta = 0; iEta < eta_bins.size(); iEta++) 
  	    std::cout << "     " << std::get<3>(eta_bins.at(iEta));
  	  std::cout << std::endl;
  	}
	// Print first column with pT bins
  	std::cout << std::get<3>(pt_bins.at(iPt)) << "  ";
	// Print pT resolution score ratios for each pT/eta bin
  	for (int iEta = 0; iEta < eta_bins.size(); iEta++) {
  	  Float_t score_num     = GetResScore(    h_res.at(jBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
  	  Float_t score_den     = GetResScore(    h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
  	  Float_t score_err_num = GetResScoreErr( h_res.at(jBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
  	  Float_t score_err_den = GetResScoreErr( h_res.at(iBr).at(iPt).at(iEta), std::get<4>(pt_br.at(iBr)) );
  	  Float_t ratio         = score_num / score_den;
  	  Float_t ratio_err     = ratio * sqrt( pow(score_err_num/score_num, 2) + pow(score_err_den/score_den, 2) );

  	  // std::cout << "            " << std::fixed << std::setprecision(2) << ratio;
  	  std::cout << "    " << std::fixed << std::setprecision(2) << ratio << " +/- " << ratio_err;
  	}
  	std::cout << std::endl;
      }
    }
  }

} // End void PrintRatios()


