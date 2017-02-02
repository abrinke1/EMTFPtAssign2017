
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
#include "../src/MacroHelper.C"         // Helpful common functions (GetMedian, GetResScore, etc.)


const int    PRTEVT  =     100000;  // When processing file, print every X events
const int    MAXEVT  =         -1;  // Maximum number of events to process
const int     NBINS  =       1600;  // Number of bins in log2(TRG pT / GEN pT) plot
const int     XMIN   =         -8;  // Minimum log2(TRG pT / GEN pT) value in plot
const int     XMAX   =          8;  // Maximum log2(TRG pT / GEN pT) value in plot  
const int     REBIN  =         10;  // Rebin log2(TRG pT / GEN pT) by X for display purposes
const double  PTMIN  =         1.;  // Minimum GEN and TRG pT values used
const double  PTMAX  =       256.;  // Maximum GEN and TRG pT values used

// const double  EMTFSC =        1.0;  // Scaling of EMTF pT to yield median ratio of 1.0 (optional)
const double  BDTGSC =        1.0;  // Scaling of BDTG pT to yield median ratio of 1.0 (optional)
const double  EMTFSC = (1./1.25);  // Scaling of EMTF pT to yield median ratio of 1.0 (optional)
// const double  BDTGSC = (1./1.02);  // Scaling of BDTG pT to yield median ratio of 1.0 (optional)

// const TString EVTWGT = "invLog10";  // Weight events: "invLog10" --> 1/(1 + log10(pT))
// const TString EVTWGT = "invLog2";   // Weight events: "invLog2" --> 1/(1 + log2(pT))
const TString EVTWGT =  "invSqrt";  // Weight events: "invSqrt"  --> 1/(sqrt(pT))
// const TString EVTWGT =    "invPt";  // Weight events: "invPt"  --> 1/pT
// const TString EVTWGT =         "";  // Don't weight events


///////////////////////////////////////
///  Main function: PtResolution()  /// 
///////////////////////////////////////

void PtResolution() {

  // Initialize empty file and chain to access each file and chain in the list
  TFile* file_tmp(0);
  TChain* ch_tmp(0);

  // List of input files
  // TString in_dir = "/afs/cern.ch/user/a/abrinke1/TMVA/EMTFPtAssign2017/";
  TString in_dir = "/afs/cern.ch/work/a/abrinke1/public/EMTF/PtAssign2017/";
  std::vector<TString> in_file_names;
  // in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_01_23_400_trees_0p002_node.root");
  // in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_01_24_vars.root");
  // in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_01_24_bends.root");
  // in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_01_26_NTrees.root");
  // in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_01_26_MaxDepth.root");
  in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_01_26_mode_15_opt.root");

  // Open all input files
  for (int i = 0; i < in_file_names.size(); i++) {
    if ( !gSystem->AccessPathName(in_file_names.at(i)) )
      file_tmp = TFile::Open( in_file_names.at(i) ); // Check if file exists
    if (!file_tmp) {
      std::cout << "ERROR: could not open data file " << in_file_names.at(i) << std::endl;
      return;
    }
  }

  // Tuple for factory directories: directory name, target pT, train chain, test chain
  // Target trigger pT: "inv" for 1/pT, "log2" for log2(pT), "" for pT
  std::vector<std::tuple<TString, TString, TChain*, TChain*>> facts;

  // // Original EMTF variables + combinations, invPt target, no weighting
  // facts.push_back( std::make_tuple("f_0x001f01fd_0x2",       "inv",  ch_tmp, ch_tmp) );

  // Optimized mode 15
  facts.push_back( std::make_tuple("f_0x001f01ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );

  // // Training options (target pT, weighting, number of variables)
  // facts.push_back( std::make_tuple("f_0x0000011d_0x2"      , "inv",  ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x0000011d_0x4"      , "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f01fd_0x2"      , "inv",  ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f01fd_0x4"      , "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001fffff_0x2"      , "inv",  ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001fffff_0x4"      , "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x0000011d_0x2_invPt", "inv",  ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x0000011d_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f01fd_0x2_invPt", "inv",  ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f01fd_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001fffff_0x2_invPt", "inv",  ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001fffff_0x4_invPt", "log2", ch_tmp, ch_tmp) );

  // // Different sets of variables
  // facts.push_back( std::make_tuple("f_0x00000004_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x00000005_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x0000000d_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x00000085_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x0000001d_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f00fd_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f0ffd_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f0fff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001fffff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x8fff0fff_0x4_invPt", "log2", ch_tmp, ch_tmp) );

  // // Different FR combinations
  // facts.push_back( std::make_tuple("f_0x001f00ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f01ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f03ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f05ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f09ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f0fff_0x4_invPt", "log2", ch_tmp, ch_tmp) );

  // // Different bend combinations
  // facts.push_back( std::make_tuple("f_0x001f01ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f11ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f31ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f51ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f91ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001ff1ff_0x4_invPt", "log2", ch_tmp, ch_tmp) );


  // Add trees from the input files to the TChain
  for (int iFt = 0; iFt < facts.size(); iFt++) {
    std::get<2>(facts.at(iFt)) = new TChain(std::get<0>(facts.at(iFt))+"/TrainTree");
    std::get<3>(facts.at(iFt)) = new TChain(std::get<0>(facts.at(iFt))+"/TestTree");
    for (int i = 0; i < in_file_names.size(); i++) {
      std::get<2>(facts.at(iFt))->Add( in_file_names.at(i) );
      std::get<3>(facts.at(iFt))->Add( in_file_names.at(i) );
    }
  }

  // TString out_file_name = "plots/PtResolution_AWB_v1_17_01_23_400_trees_0p002_node.root";
  // TString out_file_name = "plots/PtResolution_AWB_v1_17_01_24_vars_best.root";
  // TString out_file_name = "plots/PtResolution_AWB_v1_17_01_24_bends.root";
  // TString out_file_name = "plots/PtResolution_AWB_v1_17_01_26_NTrees.root";
  // TString out_file_name = "plots/PtResolution_AWB_v1_17_01_26_MaxDepth.root";
  TString out_file_name = "plots/PtResolution_AWB_v1_17_01_26_mode_15_opt.root";
  TFile *out_file = new TFile(out_file_name, "recreate");


  //////////////////////////////////////////////////
  // Set up histograms to fill and scores to compute
  //////////////////////////////////////////////////

  // Pair each factory with a vector of MVAs
  std::vector< std::pair< std::tuple<TString, TString, TChain*, TChain*>, std::vector<std::tuple<TString, float, float, float> > > > fMVAs;
  
  for (int iFt = 0; iFt < facts.size(); iFt++) {
    // Tuple for trigger pT values from different MVAs: branch name, branch value, median ratio to GEN pT (train/test)
    std::vector<std::tuple<TString, float, float, float>> MVAs;
    
    MVAs.push_back( std::make_tuple("EMTF_pt",        -99., 1., 1.) );
    
    MVAs.push_back( std::make_tuple("BDTG_AWB",       -99., 1., 1.) );
    // // MVAs.push_back( std::make_tuple("BDTG_AWB_lite",  -99., 1., 1.) );

    // MVAs.push_back( std::make_tuple("BDTG_AWB_50_trees",  -99., 1., 1.) );
    // MVAs.push_back( std::make_tuple("BDTG_AWB_100_trees", -99., 1., 1.) );
    // MVAs.push_back( std::make_tuple("BDTG_AWB_200_trees", -99., 1., 1.) );
    // MVAs.push_back( std::make_tuple("BDTG_AWB_400_trees", -99., 1., 1.) );
    // MVAs.push_back( std::make_tuple("BDTG_AWB_800_trees", -99., 1., 1.) );

    // MVAs.push_back( std::make_tuple("BDTG_AWB_3_deep",    -99., 1., 1.) );
    // MVAs.push_back( std::make_tuple("BDTG_AWB_4_deep",    -99., 1., 1.) );
    // MVAs.push_back( std::make_tuple("BDTG_AWB_5_deep",    -99., 1., 1.) );
    // MVAs.push_back( std::make_tuple("BDTG_AWB_6_deep",    -99., 1., 1.) );


    fMVAs.push_back( std::make_pair(facts.at(iFt), MVAs) );
  } // End loop: for (int iCh = 0; iCh < facts.size(); iCh++)

  std::vector<TString> fact_names;
  std::vector<TString> MVA_names;
  for (int iFt = 0; iFt < facts.size(); iFt++)
    fact_names.push_back( std::get<0>(facts.at(iFt)) );
  for (int iMVA = 0; iMVA < (fMVAs.at(0).second).size(); iMVA++)
    MVA_names.push_back( std::get<0>(fMVAs.at(0).second.at(iMVA)) );
  
  // Set bins in GEN pT and eta
  std::vector<std::tuple<TString, float, float, TString>> pt_bins;
  std::vector<std::tuple<TString, float, float, TString>> eta_bins;
  
  pt_bins.push_back( std::make_tuple("all",      1.,   1000., "    ALL    ") );

  // pt_bins.push_back( std::make_tuple("1_8",      1.,      8., "[  1,    8]") );
  // pt_bins.push_back( std::make_tuple("8_30",     8.,     30., "[  8,   30]") );
  // pt_bins.push_back( std::make_tuple("30_120",   30.,   120., "[ 30,  120]") );
  // pt_bins.push_back( std::make_tuple("120_1000", 120., 1000., "[120, 1000]") );

  pt_bins.push_back( std::make_tuple("1_4",        1.,    4., "[  1,    4]") );
  pt_bins.push_back( std::make_tuple("4_8",        4.,    8., "[  4,    8]") );
  pt_bins.push_back( std::make_tuple("8_15",       8.,   15., "[  8,   15]") );
  pt_bins.push_back( std::make_tuple("15_30",     15.,   30., "[ 15,   30]") );
  pt_bins.push_back( std::make_tuple("30_60",     30.,   60., "[ 30,   60]") );
  pt_bins.push_back( std::make_tuple("60_120",    60.,  120., "[ 60,  120]") );
  pt_bins.push_back( std::make_tuple("120_250",  120.,  250., "[120,  250]") );
  pt_bins.push_back( std::make_tuple("250_1000", 250., 1000., "[250, 1000]") );

  eta_bins.push_back( std::make_tuple("all",       1.2,  2.4 , "     ALL    ") );
  eta_bins.push_back( std::make_tuple("1p2_1p55",  1.2,  1.55, "[1.2 , 1.55]") );
  eta_bins.push_back( std::make_tuple("1p55_1p85", 1.55, 1.85, "[1.55, 1.85]") );
  eta_bins.push_back( std::make_tuple("1p85_2p1",  1.85, 2.1 , "[1.85, 2.1 ]") );
  eta_bins.push_back( std::make_tuple("2p1_2p4",   2.1,  2.4 , "[2.1 , 2.4 ]") );

  // 1D resolution plots
  std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > > h_res;
  std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > empty1a;
  std::vector< std::vector< std::pair<TH1D*, TH1D*> > > empty1b;
  std::vector< std::pair<TH1D*, TH1D*> > empty1c;

  // 2D resolution score, vs. pT / eta
  std::vector< std::vector< std::pair<TH2D*, TH2D*> > > h_score;
  std::vector< std::pair<TH2D*, TH2D*> > empty2;

  // Book 1D pT resolution histograms for each trigger and pT/eta bin
  for (int iFM = 0; iFM < fMVAs.size(); iFM++) {
    h_score.push_back( empty2 );
    h_res.push_back( empty1a );
    for (int iMVA = 0; iMVA < fMVAs.at(iFM).second.size(); iMVA++) {
      h_res.at(iFM).push_back( empty1b );
      TString ft_name  = fact_names.at(iFM);
      TString mva_name = std::get<0>(fMVAs.at(iFM).second.at(iMVA));

      BookScoreHist( iFM, iMVA, pt_bins, eta_bins, ft_name, mva_name+"_train", h_score );
      // std::cout << "Booked score hists: " << (h_score.back().back().first)->GetName() << std::endl;
      // std::cout << "               and  " << (h_score.back().back().second)->GetName() << std::endl;
      
      for (int iPt = 0; iPt < pt_bins.size(); iPt++) {
	h_res.at(iFM).at(iMVA).push_back( empty1c );
	
	for (int iEta = 0; iEta < eta_bins.size(); iEta++) {
	  BookResHist( iFM, iMVA, iPt, iEta, pt_bins, eta_bins, ft_name, mva_name, h_res );
	  // std::cout << "  * Booked res hists: " << (h_res.back().back().back().back().first)->GetName() << std::endl;
	  // std::cout << "                 and  " << (h_res.back().back().back().back().second)->GetName() << std::endl;
	}
      }
    }
  }
  
  /////////////////////////////////////////////
  // Loop over events from train and test trees
  /////////////////////////////////////////////
  for (int iFM = 0; iFM < fMVAs.size(); iFM++) {
    TString ft_name = fact_names.at(iFM);
    TString trgPt   = std::get<1>(fMVAs.at(iFM).first);
    LoopOverEvents( std::get<2>(fMVAs.at(iFM).first), ft_name, trgPt, iFM,
  		    "train", pt_bins, eta_bins, fMVAs.at(iFM).second, h_res);
    LoopOverEvents( std::get<3>(fMVAs.at(iFM).first), ft_name, trgPt, iFM,
  		    "test",  pt_bins, eta_bins, fMVAs.at(iFM).second, h_res);
  }
    
  //////////////////////////////////////////
  // Write out histograms and compute scores
  //////////////////////////////////////////

  out_file->cd();

  // Compute median ratios to GEN values, store in MVAs
  for (int iFM = 0; iFM < fMVAs.size(); iFM++) {
    TString ft_name = fact_names.at(iFM);
    StoreMedians( ft_name, iFM, pt_bins, eta_bins, h_res, fMVAs.at(iFM).second );
  }

  // Initialize vectors and of scores
  int num_FM  = fMVAs.size();
  int num_pt  = pt_bins.size();
  int num_MVA = (fMVAs.at(0).second).size() - 1; // Not including EMTF

  std::vector<float> EMTF_vec_tr     ; for (int i = 0; i < num_FM * num_pt; i++) EMTF_vec_tr.    push_back(-99);
  std::vector<float> EMTF_vec_tr_err ; for (int i = 0; i < num_FM * num_pt; i++) EMTF_vec_tr_err.push_back(-99);
  std::vector<float> EMTF_vec_te     ; for (int i = 0; i < num_FM * num_pt; i++) EMTF_vec_te.    push_back(-99);
  std::vector<float> EMTF_vec_te_err ; for (int i = 0; i < num_FM * num_pt; i++) EMTF_vec_te_err.push_back( -99);
  std::vector<float> ratio_vec_tr    ; for (int i = 0; i < num_FM * num_pt * num_MVA; i++) ratio_vec_tr    .push_back(-99);
  std::vector<float> ratio_vec_tr_err; for (int i = 0; i < num_FM * num_pt * num_MVA; i++) ratio_vec_tr_err.push_back(-99);
  std::vector<float> ratio_vec_te    ; for (int i = 0; i < num_FM * num_pt * num_MVA; i++) ratio_vec_te    .push_back(-99);
  std::vector<float> ratio_vec_te_err; for (int i = 0; i < num_FM * num_pt * num_MVA; i++) ratio_vec_te_err.push_back(-99);

  // Write ratio histograms and print out scores
  for (int iFM = 0; iFM < fMVAs.size(); iFM++) {
    TString ft_name  = fact_names.at(iFM);
    FillScores( ft_name, iFM, fMVAs.at(iFM).second, pt_bins, eta_bins, num_FM, num_pt, num_MVA, 
		EMTF_vec_tr,  EMTF_vec_tr_err,  EMTF_vec_te,  EMTF_vec_te_err,
		ratio_vec_tr, ratio_vec_tr_err, ratio_vec_te, ratio_vec_te_err,
		h_res, h_score );
  }

  // Fill 2D graphs of resolution ratio to EMTF
  FillRatioGraphs( fact_names, MVA_names, pt_bins, eta_bins, num_FM, num_pt, num_MVA,
		   EMTF_vec_tr,  EMTF_vec_tr_err,  EMTF_vec_te,  EMTF_vec_te_err,
		   ratio_vec_tr, ratio_vec_tr_err, ratio_vec_te, ratio_vec_te_err );

  // Print out ratios of scores from different triggers
  for (int iFM = 0; iFM < fMVAs.size(); iFM++) {
    TString ft_name  = fact_names.at(iFM);
    PrintRatios( ft_name, iFM, fMVAs.at(iFM).second, pt_bins, eta_bins, h_res );
  }
  
  out_file->Close();

  std::cout << "\nExiting PtResolution()\n";

} // End void PtResolution()


//////////////////////////////////////////
///  Functions used in PtResolution()  /// 
//////////////////////////////////////////

void BookScoreHist( const int iFM, const int iMVA,
		    const std::vector<std::tuple<TString, float, float, TString>> pt_bins, 
		    const std::vector<std::tuple<TString, float, float, TString>> eta_bins, 
		    const TString ft_name, const TString mva_name, 
		    std::vector< std::vector< std::pair<TH2D*, TH2D*> > >& h_score ) {

  TString h_score_str_tr = "h_score_"+ft_name+"_"+mva_name+"_train";
  TString h_score_str_te = "h_score_"+ft_name+"_"+mva_name+"_test";
  h_score.at(iFM).push_back( std::make_pair( new TH2D( h_score_str_tr, h_score_str_tr, 
						       pt_bins.size(), 0, pt_bins.size(), 
						       eta_bins.size(), 0, eta_bins.size() ),
					     new TH2D( h_score_str_te,  h_score_str_te,  
						       pt_bins.size(), 0, pt_bins.size(), 
						       eta_bins.size(), 0, eta_bins.size() ) ) );
  h_score.at(iFM).at(iMVA).first ->Sumw2();
  h_score.at(iFM).at(iMVA).second->Sumw2();
  
  h_score_str_tr.ReplaceAll("h_score_", "");
  h_score_str_te.ReplaceAll("h_score_", "");
  h_score_str_tr.ReplaceAll("_", " ");
  h_score_str_te.ReplaceAll("_", " ");
  h_score_str_tr.ReplaceAll("train", "(train)");
  h_score_str_te.ReplaceAll("test",  "(test)");
  
  h_score.at(iFM).at(iMVA).first ->SetTitle(h_score_str_tr+" p_{T} resolution score");
  h_score.at(iFM).at(iMVA).second->SetTitle(h_score_str_te+" p_{T} resolution score");
  h_score.at(iFM).at(iMVA).first ->GetXaxis()->SetTitle("p_{T} range (GeV)");
  h_score.at(iFM).at(iMVA).second->GetXaxis()->SetTitle("p_{T} range (GeV)");
  h_score.at(iFM).at(iMVA).first ->GetYaxis()->SetTitle("|#eta| range");
  h_score.at(iFM).at(iMVA).second->GetYaxis()->SetTitle("|#eta| range");
  
} // End void BookScoreHist()

void BookResHist( const int iFM, const int iMVA, const int iPt, const int iEta,
		      const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
		  const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
		  const TString ft_name, const TString mva_name, 
		  std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > >& h_res ) {
  
  TString h_res_str_tr = "h_res_"+ft_name+"_"+mva_name+"_train_pt_"+std::get<0>(pt_bins.at(iPt))+"_eta_"+std::get<0>(eta_bins.at(iEta));
  TString h_res_str_te = "h_res_"+ft_name+"_"+mva_name+"_test_pt_"+std::get<0>(pt_bins.at(iPt))+"_eta_"+std::get<0>(eta_bins.at(iEta));
  h_res.at(iFM).at(iMVA).at(iPt).push_back( std::make_pair( new TH1D(h_res_str_tr, h_res_str_tr, NBINS, XMIN, XMAX),
							    new TH1D(h_res_str_te, h_res_str_te, NBINS, XMIN, XMAX) ) );
  h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first ->Sumw2();
  h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second->Sumw2();
  
  h_res_str_tr.ReplaceAll("h_res_", "");
  h_res_str_te.ReplaceAll("h_res_", "");
  h_res_str_tr.ReplaceAll("_", " ");
  h_res_str_te.ReplaceAll("_", " ");
  h_res_str_tr.ReplaceAll("train", "(train)");
  h_res_str_te.ReplaceAll("test",  "(test)");

  std::ostringstream ss_title_tr;
  std::ostringstream ss_title_te;
  ss_title_tr << h_res_str_tr << " p_{T} resolution, "
	      << std::get<1>(pt_bins.at(iPt))   << " < p_{T} < "  << std::get<2>(pt_bins.at(iPt)) << " GeV, "
	      << std::get<1>(eta_bins.at(iEta)) << " < |#eta| < " << std::get<2>(eta_bins.at(iEta));
  ss_title_te << h_res_str_te << " p_{T} resolution, "
	      << std::get<1>(pt_bins.at(iPt))   << " < p_{T} < "  << std::get<2>(pt_bins.at(iPt)) << " GeV, "
	      << std::get<1>(eta_bins.at(iEta)) << " < |#eta| < " << std::get<2>(eta_bins.at(iEta));

  h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first ->SetTitle( (TString) ss_title_tr.str() );
  h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second->SetTitle( (TString) ss_title_te.str() );
  h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first ->GetXaxis()->SetTitle("log_{2} trigger p_{T} / GEN p_{T}");
  h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second->GetXaxis()->SetTitle("log_{2} trigger p_{T} / GEN p_{T}");
  
} // End void BookResHist()

void LoopOverEvents( TChain* chain, const TString ft_name, const TString trgPt, const int iFM, const TString tr_te,
		     const std::vector<std::tuple<TString, float, float, TString>> pt_bins, 
		     const std::vector<std::tuple<TString, float, float, TString>> eta_bins, 
		     std::vector<std::tuple<TString, float, float, float>>& MVAs,
		     std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > >& h_res ) {
  
  // Get GEN branches from the factories
  float GEN_pt_br;
  float GEN_eta_br;
  chain->SetBranchAddress("GEN_pt", &GEN_pt_br);
  chain->SetBranchAddress("GEN_eta", &GEN_eta_br);

  // Get trigger branches from the factories
  for (int iMVA = 0; iMVA < MVAs.size(); iMVA++)
    chain->SetBranchAddress( std::get<0>(MVAs.at(iMVA)), &(std::get<1>(MVAs.at(iMVA))) );

  std::cout << "\n******* About to enter the " << ft_name << " " << tr_te << " event loop *******" << std::endl;
  for (int iEvt = 0; iEvt < chain->GetEntries(); iEvt++) {
    
    if (iEvt > MAXEVT && MAXEVT > 0) break;
    if ( (iEvt % PRTEVT) == 0 ) std::cout << "*** Looking at event " << iEvt << " ***" << std::endl;
    
    chain->GetEntry(iEvt);
    
    // Loop over different trigger pT computations
    for (int iMVA = 0; iMVA < MVAs.size(); iMVA++) {

      // Access trigger pT
      double TRG_pt = std::get<1>(MVAs.at(iMVA));
      if ( not std::get<0>(MVAs.at(iMVA)).Contains("EMTF") ) {
	if ( trgPt == "inv" )
	  TRG_pt = 1. / max(0.001, TRG_pt); // Protect against negative 1/pT values
	if ( trgPt == "log2" )
	  TRG_pt = pow(2, TRG_pt);
      }
      double GEN_pt = double(GEN_pt_br);

      // Scale trigger pT to shift median TRG pT / GEN pT value to 1.0 (optional)
      if (std::get<0>(MVAs.at(iMVA)).Contains("EMTF")) TRG_pt *= EMTFSC;
      if (std::get<0>(MVAs.at(iMVA)).Contains("BDTG")) TRG_pt *= BDTGSC;

      // std::cout << "GEN_pt = " << GEN_pt << ", TRG_pt = " << TRG_pt << std::endl;

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

	  // std::cout << "Filling " << std::get<0>(MVAs.at(iMVA)) << " histogram: TRG_pt = " << TRG_pt << ", GEN_pt = " << GEN_pt << std::endl;
	  if (tr_te.Contains("train"))
	    h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first->Fill( log2( TRG_pt / GEN_pt ), evt_wgt );
	  else if (tr_te.Contains("test"))
	    h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second->Fill( log2( TRG_pt / GEN_pt ), evt_wgt );
	  else {
	    std::cout << "tr_te = " << tr_te << ", not train or test. Exiting." << std::endl;
	    break;
	  }
	  
	} // End loop: for (int iEta = 0; iEta < eta_bins.size(); iEta++)
      } // End loop: for (int iPt = 0; iPt < pt_bins.size(); iPt++)
    } // End loop: for (int iMVA = 0; iMVA < MVAs.size(); iMVA++)

  } // End loop: for (int iEvt = 0; iEvt < chain->GetEntries(); iEvt++)
  std::cout << "******* Leaving the " << ft_name << " " << tr_te << " event loop *******" << std::endl;

} // End void LoopOverEvents()

void StoreMedians( const TString ft_name, const int iFM, 
		   const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
		   const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
		   const std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > > h_res,
		   std::vector<std::tuple<TString, float, float, float>>& MVAs ) {

  std::cout << "\n" << std::string(33,'*') << "\n" << "Median trigger pT / GEN pT values" << "\n" << std::string(33,'*') << std::endl;
  for (int iMVA = 0; iMVA < MVAs.size(); iMVA++) {
    for (int iPt = 0; iPt < pt_bins.size(); iPt++) {
      if ( !(std::get<0>(pt_bins.at(iPt)).Contains("all")) ) continue;
      for (int iEta = 0; iEta < eta_bins.size(); iEta++) {
	if ( !(std::get<0>(eta_bins.at(iEta)).Contains("all")) ) continue;

	std::get<2>(MVAs.at(iMVA)) = pow(2, GetMedian( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first ));  // Store median ratio in MVAs
	std::cout << std::setw(45) << std::left << ft_name
		  << std::setw(35) << std::left << std::get<0>(MVAs.at(iMVA))+" train"
		  << std::fixed << std::setprecision(3) << std::get<2>(MVAs.at(iMVA)) << std::endl;
	std::get<3>(MVAs.at(iMVA)) = pow(2, GetMedian( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second ));  // Store median ratio in MVAs
	std::cout << std::setw(45) << std::left << ft_name
		  << std::setw(35) << std::left << std::get<0>(MVAs.at(iMVA))+" test"
		  << std::fixed << std::setprecision(3) << std::get<3>(MVAs.at(iMVA)) << std::endl;

      }
    }
  }

} // End void StoreMedians()

void FillScores( const TString ft_name, const int iFM, 
		 const std::vector<std::tuple<TString, float, float, float>> MVAs, 
		 const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
		 const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
		 const int num_FM, const int num_pt, const int num_MVA,
		 std::vector<float>& EMTF_vec_tr,  std::vector<float>& EMTF_vec_tr_err,  
		 std::vector<float>& EMTF_vec_te,  std::vector<float>& EMTF_vec_te_err,
		 std::vector<float>& ratio_vec_tr, std::vector<float>& ratio_vec_tr_err, 
		 std::vector<float>& ratio_vec_te, std::vector<float>& ratio_vec_te_err,
		 const std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > > h_res,
		 std::vector< std::vector< std::pair<TH2D*, TH2D*> > >& h_score ) {

  for (int iMVA = 0; iMVA < MVAs.size(); iMVA++) {
    
    // Print table header
    std::cout << "\n" << std::string(45,'*') << "\n" << ft_name << " " 
	      << std::get<0>(MVAs.at(iMVA)) << " resolution score" << "\n" << std::string(45,'*') << std::endl;
    
    for (int iPt = 0; iPt < pt_bins.size(); iPt++) {
      
      if (iPt == 0) {  // Print top row with eta bins
	std::cout << "              ";
	for (int iEta = 0; iEta < eta_bins.size(); iEta++) 
	  std::cout << "  " << std::get<3>(eta_bins.at(iEta));
	std::cout << std::endl;
      }
      // Print first column with pT bins
      std::cout << std::get<3>(pt_bins.at(iPt));

      // Print pT resolution score values for each pT/eta bin
      for (int iEta = 0; iEta < eta_bins.size(); iEta++) {
	
	// Print score
	std::cout << "          " << std::fixed << std::setprecision(2) << GetResScore( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) );
	// Print score with uncertainty
	// std::cout << "   " << std::fixed << std::setprecision(2) << GetResScore( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) )
	// 	    << " +/- " << GetResScoreErr( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) );
	// // Print median ratio to GEN
	// std::cout << "          " << std::fixed << std::setprecision(2) << GetMedian( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second );
	// // Print number of events
	// std::cout << std::setw(14) << std::right << std::fixed << std::setprecision(0) << h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second->Integral();
	
	// Write out resolution histogram
  	TH1D* h_tmp_tr = (TH1D*) ( (TH1D*) h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first ->Clone() )->Rebin(REBIN);  // Don't rebin h_res
  	TH1D* h_tmp_te = (TH1D*) ( (TH1D*) h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second->Clone() )->Rebin(REBIN);  // Don't rebin h_res
  	h_tmp_tr->SetLineWidth(2);
  	h_tmp_te->SetLineWidth(2);
	if ( std::get<0>(MVAs.at(iMVA)).Contains("EMTF") ) {
	  if (iFM == 0) { // Only write the EMTF histogram from one factory
	    h_tmp_tr->SetName( ((TString) h_tmp_tr->GetName()).ReplaceAll(ft_name+"_", "") );
	    h_tmp_te->SetName( ((TString) h_tmp_te->GetName()).ReplaceAll(ft_name+"_", "") );
	    h_tmp_tr->SetLineColor(kBlack);
	    h_tmp_te->SetLineColor(kBlack);
	    h_tmp_tr->Write();
	    h_tmp_te->Write();
	  }
	} else {
	  h_tmp_tr->SetLineColor(kRed);
	  h_tmp_te->SetLineColor(kRed);
	  h_tmp_tr->Write();
	  h_tmp_te->Write();
	}
	
  	TH1D* h_tmp_wgt_tr = (TH1D*) h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first ->Clone();
  	TH1D* h_tmp_wgt_te = (TH1D*) h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second->Clone();
	WeightByResScore( (*h_tmp_wgt_tr), std::get<2>(MVAs.at(iMVA)) );
	WeightByResScore( (*h_tmp_wgt_te), std::get<3>(MVAs.at(iMVA)) );
	h_tmp_wgt_tr = (TH1D*) h_tmp_wgt_tr->Rebin(REBIN);
	h_tmp_wgt_te = (TH1D*) h_tmp_wgt_te->Rebin(REBIN);
  	h_tmp_wgt_tr->SetLineWidth(2);
  	h_tmp_wgt_te->SetLineWidth(2);
	if ( std::get<0>(MVAs.at(iMVA)).Contains("EMTF") ) {
	  if (iFM == 0) { // Only write the EMTF histogram from one factory
	    h_tmp_wgt_tr->SetName( ((TString) h_tmp_tr->GetName()).ReplaceAll(ft_name+"_", "")+"_wgt" );
	    h_tmp_wgt_te->SetName( ((TString) h_tmp_te->GetName()).ReplaceAll(ft_name+"_", "")+"_wgt" );
	    // h_tmp_wgt_tr->SetTitle( (TString) h_tmp_tr->GetTitle().ReplaceAll(ft_name+" ", "")+" (weighted)");
	    // h_tmp_wgt_te->SetTitle( (TString) h_tmp_te->GetTitle().ReplaceAll(ft_name+" ", "")+" (weighted)");
	    h_tmp_wgt_tr->SetLineColor(kBlack);
	    h_tmp_wgt_te->SetLineColor(kBlack);
	    h_tmp_wgt_tr->Write();
	    h_tmp_wgt_te->Write();
	  }
	} else {
	  h_tmp_wgt_tr->SetLineColor(kRed);
	  h_tmp_wgt_te->SetLineColor(kRed);
	  h_tmp_wgt_tr->SetName( (TString) h_tmp_tr->GetName()+"_wgt" );
	  h_tmp_wgt_te->SetName( (TString) h_tmp_te->GetName()+"_wgt" );
	  h_tmp_wgt_tr->SetTitle( (TString) h_tmp_tr->GetTitle()+" (weighted)");
	  h_tmp_wgt_te->SetTitle( (TString) h_tmp_te->GetTitle()+" (weighted)");
	  h_tmp_wgt_tr->Write();
	  h_tmp_wgt_te->Write();
	}
	
	// Fill score histogram
	h_score.at(iFM).at(iMVA).first ->SetBinContent(iPt+1, iEta+1, GetResScore(    h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first , std::get<2>(MVAs.at(iMVA)) ) );
	h_score.at(iFM).at(iMVA).second->SetBinContent(iPt+1, iEta+1, GetResScore(    h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) ) );
	h_score.at(iFM).at(iMVA).first ->SetBinError  (iPt+1, iEta+1, GetResScoreErr( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first , std::get<2>(MVAs.at(iMVA)) ) );
	h_score.at(iFM).at(iMVA).second->SetBinError  (iPt+1, iEta+1, GetResScoreErr( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) ) );
	h_score.at(iFM).at(iMVA).first ->GetXaxis()->SetBinLabel( iPt+1,  std::get<3>(pt_bins.at(iPt))   );
	h_score.at(iFM).at(iMVA).second->GetXaxis()->SetBinLabel( iPt+1,  std::get<3>(pt_bins.at(iPt))   );
	h_score.at(iFM).at(iMVA).first ->GetYaxis()->SetBinLabel( iEta+1, std::get<3>(eta_bins.at(iEta)) );
	h_score.at(iFM).at(iMVA).second->GetYaxis()->SetBinLabel( iEta+1, std::get<3>(eta_bins.at(iEta)) );
	
	// Fill vectors of scores
	if ( std::get<0>(eta_bins.at(iEta)).Contains("all") ) {
	  if ( std::get<0>(MVAs.at(iMVA)).Contains("EMTF") ) {
	    EMTF_vec_tr    .at(iFM*num_pt + iPt) = GetResScore   ( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<2>(MVAs.at(iMVA)) );
	    EMTF_vec_tr_err.at(iFM*num_pt + iPt) = GetResScoreErr( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<2>(MVAs.at(iMVA)) );
	    EMTF_vec_te     .at(iFM*num_pt + iPt) = GetResScore   ( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) );
	    EMTF_vec_te_err .at(iFM*num_pt + iPt) = GetResScoreErr( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) );
	  } else {
	    ratio_vec_tr    .at(iFM*num_pt*num_MVA + (iMVA-1)*num_pt + iPt) = GetResScore   ( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first,  std::get<2>(MVAs.at(iMVA)) );
	    ratio_vec_tr_err.at(iFM*num_pt*num_MVA + (iMVA-1)*num_pt + iPt) = GetResScoreErr( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).first,  std::get<2>(MVAs.at(iMVA)) );
	    ratio_vec_te     .at(iFM*num_pt*num_MVA + (iMVA-1)*num_pt + iPt) = GetResScore   ( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) );
	    ratio_vec_te_err .at(iFM*num_pt*num_MVA + (iMVA-1)*num_pt + iPt) = GetResScoreErr( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) );
	  }
	}
	
      } // End loop: for (int iEta = 0; iEta < eta_bins.size(); iEta++)
      std::cout << std::endl;
    } // End loop: for (int iPt = 0; iPt < pt_bins.size(); iPt++)
    
    // Fill score histogram
    h_score.at(iFM).at(iMVA).first ->Write();
    h_score.at(iFM).at(iMVA).second->Write();
    
  } // End loop: for (int iMVA = 0; iMVA < MVAs.size(); iMVA++)
  
} // End void FillScores()

void FillRatioGraphs( const std::vector<TString> fact_names, const std::vector<TString> MVA_names,
		      const std::vector<std::tuple<TString, float, float, TString>> pt_bins,
		      const std::vector<std::tuple<TString, float, float, TString>> eta_bins,
		      const int num_FM, const int num_pt, const int num_MVA,
		      const std::vector<float> EMTF_vec_tr,  const std::vector<float> EMTF_vec_tr_err,  
		      const std::vector<float> EMTF_vec_te,   const std::vector<float> EMTF_vec_te_err,
		      const std::vector<float> ratio_vec_tr, const std::vector<float> ratio_vec_tr_err, 
		      const std::vector<float> ratio_vec_te,  const std::vector<float> ratio_vec_te_err ) {

  // Initialize arrays of scores
  float EMTF_arr_tr     [num_FM * num_pt]; for (int i = 0; i < num_FM * num_pt; i++) EMTF_arr_tr    [i] = EMTF_vec_tr.at(i);
  float EMTF_arr_tr_err [num_FM * num_pt]; for (int i = 0; i < num_FM * num_pt; i++) EMTF_arr_tr_err[i] = EMTF_vec_tr_err.at(i);
  float EMTF_arr_te     [num_FM * num_pt]; for (int i = 0; i < num_FM * num_pt; i++) EMTF_arr_te    [i] = EMTF_vec_te.at(i);
  float EMTF_arr_te_err [num_FM * num_pt]; for (int i = 0; i < num_FM * num_pt; i++) EMTF_arr_te_err[i] = EMTF_vec_te_err.at(i);
  float ratio_arr_tr    [num_FM * num_pt * num_MVA]; for (int i = 0; i < num_FM * num_pt * num_MVA; i++) ratio_arr_tr    [i] = ratio_vec_tr.at(i);
  float ratio_arr_tr_err[num_FM * num_pt * num_MVA]; for (int i = 0; i < num_FM * num_pt * num_MVA; i++) ratio_arr_tr_err[i] = ratio_vec_tr_err.at(i);
  float ratio_arr_te    [num_FM * num_pt * num_MVA]; for (int i = 0; i < num_FM * num_pt * num_MVA; i++) ratio_arr_te    [i] = ratio_vec_te.at(i);
  float ratio_arr_te_err[num_FM * num_pt * num_MVA]; for (int i = 0; i < num_FM * num_pt * num_MVA; i++) ratio_arr_te_err[i] = ratio_vec_te_err.at(i);

  float EMTF_arr_tr_pt     [num_pt][num_FM]; for (int i = 0; i < num_pt; i++) for (int j = 0; j < num_FM; j++) EMTF_arr_tr_pt    [i][j] = EMTF_vec_tr.at(j*num_pt + i);
  float EMTF_arr_tr_err_pt [num_pt][num_FM]; for (int i = 0; i < num_pt; i++) for (int j = 0; j < num_FM; j++) EMTF_arr_tr_err_pt[i][j] = EMTF_vec_tr_err.at(j*num_pt + i);
  float EMTF_arr_te_pt     [num_pt][num_FM]; for (int i = 0; i < num_pt; i++) for (int j = 0; j < num_FM; j++) EMTF_arr_te_pt    [i][j] = EMTF_vec_te.at(j*num_pt + i);
  float EMTF_arr_te_err_pt [num_pt][num_FM]; for (int i = 0; i < num_pt; i++) for (int j = 0; j < num_FM; j++) EMTF_arr_te_err_pt[i][j] = EMTF_vec_te_err.at(j*num_pt + i);
  float ratio_arr_tr_pt    [num_pt][num_FM * num_MVA]; for (int i = 0; i < num_pt; i++) for (int j = 0; j < num_FM * num_MVA; j++) ratio_arr_tr_pt    [i][j] = ratio_vec_tr.at(j*num_pt + i);
  float ratio_arr_tr_err_pt[num_pt][num_FM * num_MVA]; for (int i = 0; i < num_pt; i++) for (int j = 0; j < num_FM * num_MVA; j++) ratio_arr_tr_err_pt[i][j] = ratio_vec_tr_err.at(j*num_pt + i);
  float ratio_arr_te_pt    [num_pt][num_FM * num_MVA]; for (int i = 0; i < num_pt; i++) for (int j = 0; j < num_FM * num_MVA; j++) ratio_arr_te_pt    [i][j] = ratio_vec_te.at(j*num_pt + i);
  float ratio_arr_te_err_pt[num_pt][num_FM * num_MVA]; for (int i = 0; i < num_pt; i++) for (int j = 0; j < num_FM * num_MVA; j++) ratio_arr_te_err_pt[i][j] = ratio_vec_te_err.at(j*num_pt + i);

  for (int iRat = 0; iRat < num_FM * num_pt * num_MVA; iRat++) {
    int jRat = iRat % (num_FM * num_pt);
    ratio_arr_tr_err[iRat] = (ratio_arr_tr[iRat] / EMTF_arr_tr[jRat]) * sqrt( pow(ratio_arr_tr_err[iRat] / ratio_arr_tr[iRat], 2) +
									      pow(EMTF_arr_tr_err [jRat] / EMTF_arr_tr [jRat], 2) );
    ratio_arr_tr    [iRat] =  ratio_arr_tr[iRat] / EMTF_arr_tr[jRat];
    ratio_arr_te_err[iRat] = (ratio_arr_te[iRat] / EMTF_arr_te[jRat]) * sqrt( pow(ratio_arr_te_err[iRat] / ratio_arr_te[iRat], 2) +
									      pow(EMTF_arr_te_err [jRat] / EMTF_arr_te [jRat], 2) );
    ratio_arr_te    [iRat] =  ratio_arr_te[iRat] / EMTF_arr_te[jRat];
  }

  for (int iPt = 0; iPt < num_pt; iPt++) {
    for (int iRat = 0; iRat < num_FM * num_MVA; iRat++) {
      int jRat = iRat % num_FM;
      ratio_arr_tr_err_pt[iPt][iRat] = (ratio_arr_tr_pt[iPt][iRat] / EMTF_arr_tr_pt[iPt][jRat]) * sqrt( pow(ratio_arr_tr_err_pt[iPt][iRat] / ratio_arr_tr_pt[iPt][iRat], 2) +
												   pow(EMTF_arr_tr_err_pt [iPt][jRat] / EMTF_arr_tr_pt [iPt][jRat], 2) );
      ratio_arr_tr_pt    [iPt][iRat] =  ratio_arr_tr_pt[iPt][iRat] / EMTF_arr_tr_pt[iPt][jRat];
      ratio_arr_te_err_pt[iPt][iRat] = (ratio_arr_te_pt[iPt][iRat] / EMTF_arr_te_pt[iPt][jRat]) * sqrt( pow(ratio_arr_te_err_pt[iPt][iRat] / ratio_arr_te_pt[iPt][iRat], 2) +
												   pow(EMTF_arr_te_err_pt [iPt][jRat] / EMTF_arr_te_pt [iPt][jRat], 2) );
      ratio_arr_te_pt    [iPt][iRat] =  ratio_arr_te_pt[iPt][iRat] / EMTF_arr_te_pt[iPt][jRat];
    }
  }
    
  float x_arr_tr[num_FM * num_pt * num_MVA];  
  float x_arr_te[num_FM * num_pt * num_MVA];
  float x_arr_zero [num_FM * num_pt * num_MVA];
  for (int i = 0; i < num_FM * num_pt * num_MVA; i++) {
    x_arr_tr[i]    = i*1.0 - 0.1;
    x_arr_te[i]    = i*1.0 + 0.1;
    x_arr_zero [i] = 0.0;
  }
  TGraphErrors* h_ratio_tr = new TGraphErrors(num_FM * num_pt * num_MVA, x_arr_tr, ratio_arr_tr, x_arr_zero, ratio_arr_tr_err);
  TGraphErrors* h_ratio_te = new TGraphErrors(num_FM * num_pt * num_MVA, x_arr_te, ratio_arr_te, x_arr_zero, ratio_arr_te_err);

  TGraphErrors* h_ratio_tr_pt[num_pt];
  TGraphErrors* h_ratio_te_pt[num_pt];
  for (int i = 0; i < num_pt; i++) {
    
    float x_arr_tr_pt[num_FM * num_MVA];  
    float x_arr_te_pt[num_FM * num_MVA];
    float x_arr_zero_pt[num_FM * num_MVA];

    for (int j = 0; j < num_FM * num_MVA; j++) {
      x_arr_tr_pt[j]   = j*1.0 - 0.1;
      x_arr_te_pt[j]   = j*1.0 + 0.1;
      x_arr_zero_pt[j] = 0.0;
    }

    h_ratio_tr_pt[i] = new TGraphErrors(num_FM * num_MVA, x_arr_tr_pt, ratio_arr_tr_pt[i], x_arr_zero_pt, ratio_arr_tr_err_pt[i]);
    h_ratio_te_pt[i] = new TGraphErrors(num_FM * num_MVA, x_arr_te_pt, ratio_arr_te_pt[i], x_arr_zero_pt, ratio_arr_te_err_pt[i]);
  }

  h_ratio_tr->SetName("h_ratio_train");
  h_ratio_tr->SetTitle("Ratio of p_{T} resolution scores");
  h_ratio_tr->GetYaxis()->SetTitle("score / EMTF score");
  h_ratio_tr->SetMarkerSize(1);
  h_ratio_tr->SetMarkerStyle(20);
  h_ratio_tr->SetMarkerColor(kBlack);
  h_ratio_tr->SetLineWidth(2);
  h_ratio_tr->SetLineColor(kBlack);
  for (int i = 0; i < num_FM * num_pt * num_MVA; i++) {
    float bin_val = h_ratio_tr->GetXaxis()->FindBin(i*1.0);
    if ( (i % num_pt) == 0 )
      h_ratio_tr->GetXaxis()->SetBinLabel( bin_val, fact_names.at(i / (num_pt*num_MVA))+" "+MVA_names.at( ((i % (num_MVA*num_pt)) / num_pt) + 1 ) );
    else                     h_ratio_tr->GetXaxis()->SetBinLabel( bin_val, std::get<3>(pt_bins.at(i % num_pt)) );
  }
  h_ratio_tr->Write();

  h_ratio_te->SetName("h_ratio_test");
  h_ratio_te->SetTitle("Ratio of p_{T} resolution scores");
  h_ratio_te->GetYaxis()->SetTitle("score / EMTF score");
  h_ratio_te->SetMarkerSize(1);
  h_ratio_te->SetMarkerStyle(21);
  h_ratio_te->SetMarkerColor(kRed);
  h_ratio_te->SetLineWidth(2);
  h_ratio_te->SetLineColor(kRed);
  for (int i = 0; i < num_FM * num_pt * num_MVA; i++) {
    float bin_val = h_ratio_te->GetXaxis()->FindBin(i*1.0);
    if ( (i % num_pt) == 0 ) h_ratio_te->GetXaxis()->SetBinLabel( bin_val, fact_names.at(i / (num_pt*num_MVA))+" "+MVA_names.at( ((i % (num_MVA*num_pt)) / num_pt) + 1 ) );
    else                     h_ratio_te->GetXaxis()->SetBinLabel( bin_val, std::get<3>(pt_bins.at(i % num_pt)) );
  }
  h_ratio_te->Write();


  for (int i = 0; i < num_pt; i++) {

    std::ostringstream ss_title;
    ss_title << "Ratio of p_{T} resolution scores, " << std::get<1>(pt_bins.at(i)) << " < p_{T} < " << std::get<2>(pt_bins.at(i)) << " GeV";
    
    h_ratio_tr_pt[i]->SetName("h_ratio_train_pt_"+std::get<0>(pt_bins.at(i)));
    h_ratio_tr_pt[i]->SetTitle( (TString) ss_title.str() );
    h_ratio_tr_pt[i]->GetYaxis()->SetTitle("score / EMTF score");
    h_ratio_tr_pt[i]->SetMarkerSize(1);
    h_ratio_tr_pt[i]->SetMarkerStyle(20);
    h_ratio_tr_pt[i]->SetMarkerColor(kBlack);
    h_ratio_tr_pt[i]->SetLineWidth(2);
    h_ratio_tr_pt[i]->SetLineColor(kBlack);
    for (int j = 0; j < num_FM * num_MVA; j++) {
      float bin_val = h_ratio_tr_pt[i]->GetXaxis()->FindBin(j*1.0);
      h_ratio_tr_pt[i]->GetXaxis()->SetBinLabel( bin_val, fact_names.at(j / (num_MVA))+" "+MVA_names.at( (j % (num_MVA)) + 1 ) );
    }
    h_ratio_tr_pt[i]->Write();
    
    h_ratio_te_pt[i]->SetName("h_ratio_test_pt_"+std::get<0>(pt_bins.at(i)));
    h_ratio_te_pt[i]->SetTitle( (TString) ss_title.str() );
    h_ratio_te_pt[i]->GetYaxis()->SetTitle("score / EMTF score");
    h_ratio_te_pt[i]->SetMarkerSize(1);
    h_ratio_te_pt[i]->SetMarkerStyle(21);
    h_ratio_te_pt[i]->SetMarkerColor(kRed);
    h_ratio_te_pt[i]->SetLineWidth(2);
    h_ratio_te_pt[i]->SetLineColor(kRed);
    for (int j = 0; j < num_FM * num_MVA; j++) {
      float bin_val = h_ratio_te_pt[i]->GetXaxis()->FindBin(j*1.0);
      h_ratio_te_pt[i]->GetXaxis()->SetBinLabel( bin_val, fact_names.at(j / (num_MVA))+" "+MVA_names.at( (j % (num_MVA)) + 1 ) );
    }
    h_ratio_te_pt[i]->Write();
    
  }
  

} // End FillRatioGraphs()

void PrintRatios( const TString ft_name, const int iFM,
		  const std::vector<std::tuple<TString, float, float, float>> MVAs,
		  const std::vector<std::tuple<TString, float, float, TString>> pt_bins, 
		  const std::vector<std::tuple<TString, float, float, TString>> eta_bins, 
		  const std::vector< std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > > h_res ) {

  for (int iMVA = 1; iMVA < MVAs.size(); iMVA++) {  // Numerator trigger
    // Print table header
    std::cout << "\n" << std::string(70,'*') << "\n" << std::get<0>(MVAs.at(iMVA)) << " (test) / "
	      << std::get<0>(MVAs.at(0)) << " score ratio" << "\n" << std::string(70,'*') << std::endl;
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
	Float_t score_num     = GetResScore(    h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) );
	Float_t score_den     = GetResScore(    h_res.at(iFM).at(0)   .at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) );
	Float_t score_err_num = GetResScoreErr( h_res.at(iFM).at(iMVA).at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) );
	Float_t score_err_den = GetResScoreErr( h_res.at(iFM).at(0)   .at(iPt).at(iEta).second, std::get<3>(MVAs.at(iMVA)) );
	Float_t ratio         = score_num / score_den;
	Float_t ratio_err     = ratio * sqrt( pow(score_err_num/score_num, 2) + pow(score_err_den/score_den, 2) );
	
	// std::cout << "            " << std::fixed << std::setprecision(2) << ratio;
	std::cout << "    " << std::fixed << std::setprecision(2) << ratio << " +/- " << ratio_err;
      } // End loop: for (int iEta = 0; iEta < eta_bins.size(); iEta++)
      std::cout << std::endl;
    } // End loop: for (int iPt = 0; iPt < pt_bins.size(); iPt++)
  } // End loop: for (int iMVA = 1; iMVA < MVAs.size(); iMVA++)

} // End void PrintRatios()


