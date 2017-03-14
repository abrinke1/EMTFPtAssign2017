
////////////////////////////////////////////////////////
///         Macro to plot the rate reduction         ///
///         at the 95% efficiency threshold          ///
///         for different MVAs and modes             ///
///         Andrew Brinkerhoff 30.01.17              ///
///                                                  ///
////////////////////////////////////////////////////////

#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

#include <iomanip>  // std::cout formatting

#include "../interface/RateVsEff.h"  // Function declarations
#include "../src/MacroHelper.C"      // Helpful common functions (GetMedian, GetResScore, etc.)


const int    PRTEVT  =    100000;  // When processing file, print every X events
const int    MAXEVT  =        -1;  // Maximum number of events to process
const double BIT     =      0.01;  // Shift to trigger pT to get it off bin edge
const double PTMIN   =       0.0;
const double PTMAX   =      50.0;
const int    PTBINS  =        50;
const double ETAMIN  =      1.24;  // Minimum GEN / trigger eta to consider
const double ETAMAX  =      2.40;  // Maximum GEN / trigger eta to consider


////////////////////////////////////
///  Main function: RateVsEff()  ///
////////////////////////////////////

void RateVsEff() {

  // Initialize empty file and chain to access each file and chain in the list
  TFile* file_tmp(0);
  TChain* ch_tmp(0);

  // List of input files                                                                                                                       
  TString in_dir = "/afs/cern.ch/work/a/abrinke1/public/EMTF/PtAssign2017/";
  // in_dir = "./";
  std::vector<TString> in_file_names;
  // in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_01_31_mode_15_opt.root");
  // in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_02_02_mode_15_opt_clean.root");
  // in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_02_09_mode_15_Sq.root");
  // in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_02_24_mode_15_opt_pt_1_256_clean.root");
  // in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_02_24_mode_15_opt_bend1_pt_1_256_clean.root");
  in_file_names.push_back(in_dir+"PtRegression_AWB_v1_17_03_14_mode_15_opt_bends_dTh_pt_1_256_clean.root");

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

  // // Original EMTF variables, invPt target, no weighting
  // facts.push_back( std::make_tuple("f_0x0000011d_0x2",         "inv",  ch_tmp, ch_tmp) );
  // Optimized mode 15
  facts.push_back( std::make_tuple("f_0x001f01ff_0x4_invPt",   "log2", ch_tmp, ch_tmp) );
  // facts.push_back( std::make_tuple("f_0x001f01ff_0x4_invPtSq", "log2", ch_tmp, ch_tmp) );
  // // Optimized mode 15 + bend1
  // facts.push_back( std::make_tuple("f_0x001f11ff_0x4_invPt",   "log2", ch_tmp, ch_tmp) );
  // Optimized mode 15 + dTheta14 + bend1
  facts.push_back( std::make_tuple("f_0x401f11ff_0x4_invPt",   "log2", ch_tmp, ch_tmp) );
  // Optimized mode 15 + dTheta14 + bend1 + bendMax3
  facts.push_back( std::make_tuple("f_0xc01f11ff_0x4_invPt",   "log2", ch_tmp, ch_tmp) );

  
  // Add trees from the input files to the TChain
  for (int iFt = 0; iFt < facts.size(); iFt++) {
    std::get<2>(facts.at(iFt)) = new TChain(std::get<0>(facts.at(iFt))+"/TrainTree");
    std::get<3>(facts.at(iFt)) = new TChain(std::get<0>(facts.at(iFt))+"/TestTree");
    for (int i = 0; i < in_file_names.size(); i++) {
      std::get<2>(facts.at(iFt))->Add( in_file_names.at(i) );
      std::get<3>(facts.at(iFt))->Add( in_file_names.at(i) );
    }
  }

  TString out_file_name = "plots/RateVsEff.root";
  TFile *out_file = new TFile(out_file_name, "recreate");

  
  // Pair each factory with a vector of MVAs
  std::vector< std::pair< std::tuple<TString, TString, TChain*, TChain*>, std::vector<std::tuple<TString, float> > > > fMVAs;
  
  for (int iFt = 0; iFt < facts.size(); iFt++) {
    // Tuple for trigger pT values from different MVAs: branch name, branch value
    std::vector<std::tuple<TString, float>> MVAs;

    MVAs.push_back( std::make_tuple("EMTF_pt",        -99.) );
    
    MVAs.push_back( std::make_tuple("BDTG_AWB",       -99.) );
    // MVAs.push_back( std::make_tuple("BDTG_AWB_Sq",    -99.) );


    fMVAs.push_back( std::make_pair(facts.at(iFt), MVAs) );
  } // End loop: for (int iCh = 0; iCh < facts.size(); iCh++)

  std::vector<TString> fact_names;
  std::vector<TString> MVA_names;
  for (int iFt = 0; iFt < facts.size(); iFt++)
    fact_names.push_back( std::get<0>(facts.at(iFt)) );
  for (int iMVA = 0; iMVA < (fMVAs.at(0).second).size(); iMVA++)
    MVA_names.push_back( std::get<0>(fMVAs.at(0).second.at(iMVA)) );

  // Set efficiency thresholds
  std::vector<int> eff_cuts;
  // eff_cuts.push_back(85);
  eff_cuts.push_back(90);

  // 2D events: trigger vs. GEN pT
  std::vector< std::vector<std::pair<TH2D*, TH2D*> > > h_pt;
  // 2D efficiency: trigger vs. GEN pT
  std::vector< std::vector<std::pair<TH2D*, TH2D*> > > h_eff;
  // 1D counts: vs. threshold pT
  std::vector< std::vector<TH1D*> > h_count;
  // 1D rate: vs. XX% threshold pT
  std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > > h_rate;
  // 1D turn-on efficiency: trigger pT = 20 GeV
  std::vector< std::vector<TH1D*> > h_turn_on;
  // 1D scaling histogram to scale by pT^2
  TH1D* h_scale_pt = new TH1D( "h_scale_pt", "h_scale_pt", PTBINS, PTMIN, PTMAX );
  for (int i = 1; i <= PTBINS; i++) {
    float low_edge = max(PTMIN + (i - 1)*(PTMAX - PTMIN)/PTBINS, 1.0);
    float hi_edge = max(PTMIN + i*(PTMAX - PTMIN)/PTBINS, 1.0);
    h_scale_pt->SetBinContent( i, pow(2/(1/low_edge + 1/hi_edge), 2) );
    h_scale_pt->SetBinError(i, 0.);
  }

  std::vector< std::vector< std::pair<TH1D*, TH1D*> > > empty1;
  std::vector< std::pair<TH1D*, TH1D*> > empty1a;
  std::vector<TH1D*> empty1b;
  std::vector< std::pair<TH2D*, TH2D*> > empty2;

  // Book histograms for each trigger
  for (int iFM = 0; iFM < fMVAs.size(); iFM++) {
    h_pt.push_back( empty2 );
    h_eff.push_back( empty2 );
    h_count.push_back( empty1b );
    h_turn_on.push_back( empty1b );
    h_rate.push_back( empty1 );
    for (int iMVA = 0; iMVA < fMVAs.at(iFM).second.size(); iMVA++) {
      h_rate.at(iFM).push_back( empty1a );
      TString ft_name  = fact_names.at(iFM);
      TString mva_name = std::get<0>(fMVAs.at(iFM).second.at(iMVA));

      BookPtHist( iFM, iMVA, ft_name, mva_name, h_pt );
      BookEffHist( iFM, iMVA, ft_name, mva_name, h_eff );
      BookCountHist( iFM, iMVA, ft_name, mva_name, h_count );
      BookTurnOnHist( iFM, iMVA, ft_name, mva_name, h_turn_on );
      for (int iEff = 0; iEff < eff_cuts.size(); iEff++)
	BookRateHist( iFM, iMVA, iEff, ft_name, mva_name, eff_cuts.at(iEff), h_rate );
    }
  }
    

  for (int iFM = 0; iFM < fMVAs.size(); iFM++) {
    TString ft_name = fact_names.at(iFM);
    TString trgPt   = std::get<1>(fMVAs.at(iFM).first);
    
    // Fill 2D pT and counts vs. pT histograms
    LoopOverEvents( std::get<2>(fMVAs.at(iFM).first), ft_name, trgPt, iFM,
		    "train", fMVAs.at(iFM).second, h_pt, h_count);
    LoopOverEvents( std::get<3>(fMVAs.at(iFM).first), ft_name, trgPt, iFM,
		    "test", fMVAs.at(iFM).second, h_pt, h_count);
    
    for (int iMVA = 0; iMVA < fMVAs.at(iFM).second.size(); iMVA++) {
      
      // Compute 2D efficiency histograms and rate at efficiency threshold histograms
      for (int iBin = 1; iBin <= PTBINS; iBin++) { // Loop over GEN pT bins (x-axis)

	// Efficiency from previous bins
	float preff1[4] = {0, 0, 0, 0};
	float preff2[4] = {0, 0, 0, 0};
	
	for (int jBin = PTBINS; jBin >= 1; jBin--) { // Loop over trigger pT bins (y-axis)
	  float num1 = h_pt.at(iFM).at(iMVA).first ->Integral(iBin, iBin, jBin, PTBINS); // Events passing trigger pT cut
	  float den1 = h_pt.at(iFM).at(iMVA).first ->Integral(iBin, iBin,    1, PTBINS); // All events in GEN pT bin
	  float eff1 = num1 / den1;
	  float num2 = h_pt.at(iFM).at(iMVA).second->Integral(iBin, iBin, jBin, PTBINS);
	  float den2 = h_pt.at(iFM).at(iMVA).second->Integral(iBin, iBin,    1, PTBINS);
	  float eff2 = num2 / den2;

	  h_eff.at(iFM).at(iMVA).first ->SetBinContent(iBin, jBin, eff1);
	  h_eff.at(iFM).at(iMVA).first ->SetBinError  (iBin, jBin, eff1 * sqrt( (1/num1) + (1/den1) ) );
	  h_eff.at(iFM).at(iMVA).second->SetBinContent(iBin, jBin, eff2);
	  h_eff.at(iFM).at(iMVA).second->SetBinError  (iBin, jBin, eff2 * sqrt( (1/num2) + (1/den2) ) );

	  if ( h_eff.at(iFM).at(iMVA).first->GetYaxis()->GetBinLowEdge(jBin) > 19.9 &&
	       h_eff.at(iFM).at(iMVA).first->GetYaxis()->GetBinLowEdge(jBin) < 20.1 ) {
	    h_turn_on.at(iFM).at(iMVA)->SetBinContent(iBin, eff2);
	    h_turn_on.at(iFM).at(iMVA)->SetBinError(iBin, eff2 * sqrt( (1/num2) + (1/den2) ) );
	  }

	  for (int iEff = 0; iEff < eff_cuts.size(); iEff++) {
	    float thresh = eff_cuts.at(iEff)*0.01;
	    bool aboveT1 = (eff1 > thresh || preff1[0] > thresh || preff1[1] > thresh || preff1[2] > thresh || preff1[3] > thresh);
	    bool aboveT2 = (eff2 > thresh || preff2[0] > thresh || preff2[1] > thresh || preff2[2] > thresh || preff2[3] > thresh);
	    bool belowT1 = (eff1 < thresh || preff1[0] < thresh || preff1[1] < thresh || preff1[2] < thresh || preff1[3] < thresh);
	    bool belowT2 = (eff2 < thresh || preff2[0] < thresh || preff2[1] < thresh || preff2[2] < thresh || preff2[3] < thresh);

	    // Find the closest efficiency over 5 EMTF pT cuts
	    float close1 = fabs(eff1 - thresh);
	    float close2 = fabs(eff2 - thresh);
	    int  iClose1 = 0;
	    int  iClose2 = 0;
	    for (int iPre = 0; iPre < 5; iPre++) {
	      if (fabs(preff1[iPre] - thresh) < close1) {
		close1  = fabs(preff1[iPre] - thresh);
		iClose1 = iPre+1;
	      }
	      if (fabs(preff2[iPre] - thresh) < close2) {
		close2  = fabs(preff2[iPre] - thresh);
		iClose2 = iPre+1;
	      }
	    }

	    // bool noRate1 = h_rate.at(iFM).at(iMVA).at(iEff).first ->GetBinContent(iBin) < 0;
	    // bool noRate2 = h_rate.at(iFM).at(iMVA).at(iEff).second->GetBinContent(iBin) < 0;
	    
	    if (aboveT1 && belowT1) // && noRate1?
	      h_rate.at(iFM).at(iMVA).at(iEff).first ->SetBinContent( iBin, h_count.at(iFM).at(iMVA)->Integral(jBin+iClose1, PTBINS) );
	    if (aboveT2 && belowT2) // && noRate2?
	      h_rate.at(iFM).at(iMVA).at(iEff).second->SetBinContent( iBin, h_count.at(iFM).at(iMVA)->Integral(jBin+iClose2, PTBINS) );
	  } // End loop: for (int iEff = 0; iEff < eff_cuts.size(); iEff++)

	  // Push back previous efficiencies
	  preff1[3] = preff1[2];  preff1[2] = preff1[1];  preff1[1] = preff1[0];  preff1[0] = eff1;
	  preff2[3] = preff2[2];  preff2[2] = preff2[1];  preff2[1] = preff2[0];  preff2[0] = eff2;

	} // End loop: for (int jBin = PTBINS; jBin >= 1; jBin--)
      } // End loop: for (int iBin = 1; iBin <= PTBINS; iBin++)

    } // End loop: for (int iMVA = 0; iMVA < fMVAs.at(iFM).second.size(); iMVA++)
  } // End loop: for (int iFM = 0; iFM < fMVAs.size(); iFM++)


  out_file->cd();
  for (int iFM = 0; iFM < fMVAs.size(); iFM++) {
    for (int iMVA = 0; iMVA < fMVAs.at(iFM).second.size(); iMVA++) {
      h_pt.at(iFM).at(iMVA).first ->Write();
      h_pt.at(iFM).at(iMVA).second->Write();
      h_eff.at(iFM).at(iMVA).first ->Write();
      h_eff.at(iFM).at(iMVA).second->Write();
      h_count.at(iFM).at(iMVA)->Write();
      h_turn_on.at(iFM).at(iMVA)->Write();
      for (int iEff = 0; iEff < eff_cuts.size(); iEff++) {
	h_rate.at(iFM).at(iMVA).at(iEff).first ->Write();
	h_rate.at(iFM).at(iMVA).at(iEff).second->Write();

	h_rate.at(iFM).at(iMVA).at(iEff).first ->SetName( ((TString) h_rate.at(iFM).at(iMVA).at(iEff).first ->GetName()).ReplaceAll("train", "train_wgt") );
	h_rate.at(iFM).at(iMVA).at(iEff).second->SetName( ((TString) h_rate.at(iFM).at(iMVA).at(iEff).second->GetName()).ReplaceAll("test", "test_wgt") );
	h_rate.at(iFM).at(iMVA).at(iEff).first ->Multiply( h_scale_pt );
	h_rate.at(iFM).at(iMVA).at(iEff).second->Multiply( h_scale_pt );
	h_rate.at(iFM).at(iMVA).at(iEff).first ->Write();
	h_rate.at(iFM).at(iMVA).at(iEff).second->Write();

	if (iMVA > 0) {
	  h_rate.at(iFM).at(iMVA).at(iEff).first ->SetName( ((TString) h_rate.at(iFM).at(iMVA).at(iEff).first ->GetName()).ReplaceAll("train_wgt", "train") );
	  h_rate.at(iFM).at(iMVA).at(iEff).second->SetName( ((TString) h_rate.at(iFM).at(iMVA).at(iEff).second->GetName()).ReplaceAll("test_wgt", "test") );
	  h_rate.at(iFM).at(iMVA).at(iEff).first ->SetName( ((TString) h_rate.at(iFM).at(iMVA).at(iEff).first ->GetName()).ReplaceAll("rate", "rate_ratio") );
	  h_rate.at(iFM).at(iMVA).at(iEff).second->SetName( ((TString) h_rate.at(iFM).at(iMVA).at(iEff).second->GetName()).ReplaceAll("rate", "rate_ratio") );
	  h_rate.at(iFM).at(iMVA).at(iEff).first ->Divide( h_rate.at(iFM).at(0).at(iEff).first );
	  h_rate.at(iFM).at(iMVA).at(iEff).second->Divide( h_rate.at(iFM).at(0).at(iEff).second);
	  h_rate.at(iFM).at(iMVA).at(iEff).first ->Write();
	  h_rate.at(iFM).at(iMVA).at(iEff).second->Write();
	}
	
      }
    }
  }

  
  out_file->Close();

  std::cout << "\nExiting RateVsEff()" << std::endl;

} // End void RateVsEff()


void BookPtHist( const int iFM, const int iMVA, const TString ft_name, const TString mva_name,
		 std::vector< std::vector<std::pair<TH2D*, TH2D*> > >& h_pt ) {
  
  TString h_pt_str_tr = "h_pt_"+ft_name+"_"+mva_name+"_train";
  TString h_pt_str_te = "h_pt_"+ft_name+"_"+mva_name+"_test";
  h_pt.at(iFM).push_back( std::make_pair( new TH2D( h_pt_str_tr, h_pt_str_tr,
						    PTBINS, PTMIN, PTMAX, PTBINS, PTMIN, PTMAX ),
					  new TH2D( h_pt_str_te,  h_pt_str_te,
						    PTBINS, PTMIN, PTMAX, PTBINS, PTMIN, PTMAX ) ) );

  h_pt.at(iFM).at(iMVA).first ->Sumw2();
  h_pt.at(iFM).at(iMVA).second->Sumw2();
  
  h_pt_str_tr.ReplaceAll("h_pt_", "");
  h_pt_str_te.ReplaceAll("h_pt_", "");
  h_pt_str_tr.ReplaceAll("_", " ");
  h_pt_str_te.ReplaceAll("_", " ");
  h_pt_str_tr.ReplaceAll("train", "(train)");
  h_pt_str_te.ReplaceAll("test",  "(test)");
  
  h_pt.at(iFM).at(iMVA).first ->SetTitle(h_pt_str_tr+" vs. GEN p_{T}");
  h_pt.at(iFM).at(iMVA).second->SetTitle(h_pt_str_te+" vs. GEN p_{T}");
  h_pt.at(iFM).at(iMVA).first ->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  h_pt.at(iFM).at(iMVA).second->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  h_pt.at(iFM).at(iMVA).first ->GetYaxis()->SetTitle(ft_name+" "+mva_name+" (train) p_{T} (GeV)");
  h_pt.at(iFM).at(iMVA).second->GetYaxis()->SetTitle(ft_name+" "+mva_name+" (test) p_{T} (GeV)");

} // End BookPtHist()
  

void BookEffHist( const int iFM, const int iMVA, const TString ft_name, const TString mva_name,
		 std::vector< std::vector<std::pair<TH2D*, TH2D*> > >& h_eff ) {
  
  TString h_eff_str_tr = "h_eff_"+ft_name+"_"+mva_name+"_train";
  TString h_eff_str_te = "h_eff_"+ft_name+"_"+mva_name+"_test";
  h_eff.at(iFM).push_back( std::make_pair( new TH2D( h_eff_str_tr, h_eff_str_tr,
						     PTBINS, PTMIN, PTMAX, PTBINS, PTMIN, PTMAX ),
					   new TH2D( h_eff_str_te,  h_eff_str_te,
						     PTBINS, PTMIN, PTMAX, PTBINS, PTMIN, PTMAX ) ) );
  
  h_eff.at(iFM).at(iMVA).first ->Sumw2();
  h_eff.at(iFM).at(iMVA).second->Sumw2();
  
  h_eff_str_tr.ReplaceAll("h_eff_", "");
  h_eff_str_te.ReplaceAll("h_eff_", "");
  h_eff_str_tr.ReplaceAll("_", " ");
  h_eff_str_te.ReplaceAll("_", " ");
  h_eff_str_tr.ReplaceAll("train", "(train)");
  h_eff_str_te.ReplaceAll("test",  "(test)");
  
  h_eff.at(iFM).at(iMVA).first ->SetTitle(h_eff_str_tr+" vs. GEN p_{T} efficiency");
  h_eff.at(iFM).at(iMVA).second->SetTitle(h_eff_str_te+" vs. GEN p_{T} efficiency");
  h_eff.at(iFM).at(iMVA).first ->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  h_eff.at(iFM).at(iMVA).second->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  h_eff.at(iFM).at(iMVA).first ->GetYaxis()->SetTitle(ft_name+" "+mva_name+" (train) p_{T} (GeV)");
  h_eff.at(iFM).at(iMVA).second->GetYaxis()->SetTitle(ft_name+" "+mva_name+" (test) p_{T} (GeV)");
  
} // End BookEffHist()
  
		 
void BookCountHist( const int iFM, const int iMVA, const TString ft_name, const TString mva_name,
		    std::vector< std::vector<TH1D*> >& h_count ) {
  
  TString h_count_str = "h_count_"+ft_name+"_"+mva_name;
  h_count.at(iFM).push_back( new TH1D( h_count_str, h_count_str, PTBINS, PTMIN, PTMAX ) );
  
  h_count.at(iFM).at(iMVA)->Sumw2();
  h_count.at(iFM).at(iMVA)->SetLineWidth(2);
  h_count.at(iFM).at(iMVA)->SetLineColor( (iMVA == 0 ? 1 : iFM+2) );
  
} // End BookCountHist()

		 
void BookTurnOnHist( const int iFM, const int iMVA, const TString ft_name, const TString mva_name,
		    std::vector< std::vector<TH1D*> >& h_turn_on ) {
  
  TString h_turn_on_str = "h_turn_on_"+ft_name+"_"+mva_name;
  h_turn_on.at(iFM).push_back( new TH1D( h_turn_on_str, h_turn_on_str, PTBINS, PTMIN, PTMAX ) );
  
  h_turn_on.at(iFM).at(iMVA)->Sumw2();
  h_turn_on.at(iFM).at(iMVA)->SetLineWidth(2);
  h_turn_on.at(iFM).at(iMVA)->SetLineColor( (iMVA == 0 ? 1 : iFM+2) );
  
} // End BookTurnOnHist()

		 
void BookRateHist( const int iFM, const int iMVA, const int iEff,
		   const TString ft_name, const TString mva_name, const int eff_cut,
		   std::vector< std::vector< std::vector< std::pair<TH1D*, TH1D*> > > >& h_rate ) {
  
  TString h_rate_str_tr;
  TString h_rate_str_te;
  h_rate_str_tr.Form("h_rate_%s_%s_%d_train", ft_name.Data(), mva_name.Data(), eff_cut);
  h_rate_str_te.Form("h_rate_%s_%s_%d_test",  ft_name.Data(), mva_name.Data(), eff_cut);
  
  h_rate.at(iFM).at(iMVA).push_back( std::make_pair( new TH1D( h_rate_str_tr, h_rate_str_tr,
							       PTBINS, PTMIN, PTMAX ),
						     new TH1D( h_rate_str_te,  h_rate_str_te,
							       PTBINS, PTMIN, PTMAX ) ) );
  
  h_rate.at(iFM).at(iMVA).at(iEff).first ->Sumw2();
  h_rate.at(iFM).at(iMVA).at(iEff).second->Sumw2();
  h_rate_str_tr.Form("%s %s (train) vs. %d%s efficiency p_{T} cut (GeV)", ft_name.Data(), mva_name.Data(), eff_cut, "%");
  h_rate_str_te.Form("%s %s (test) vs. %d%s efficiency p_{T} cut (GeV)",  ft_name.Data(), mva_name.Data(), eff_cut, "%");
  
  h_rate.at(iFM).at(iMVA).at(iEff).first ->SetTitle(h_rate_str_tr);
  h_rate.at(iFM).at(iMVA).at(iEff).second->SetTitle(h_rate_str_te);
  h_rate.at(iFM).at(iMVA).at(iEff).first ->GetXaxis()->SetTitle("p_{T} threshold (GeV)");
  h_rate.at(iFM).at(iMVA).at(iEff).second->GetXaxis()->SetTitle("p_{T} threshold (GeV)");
  h_rate.at(iFM).at(iMVA).at(iEff).first ->GetYaxis()->SetTitle(ft_name+" "+mva_name+" (train) rate");
  h_rate.at(iFM).at(iMVA).at(iEff).second->GetYaxis()->SetTitle(ft_name+" "+mva_name+" (test) rate");

  for (int i = 1; i <= PTBINS; i++) {
    h_rate.at(iFM).at(iMVA).at(iEff).first ->SetBinContent(i, -99);
    h_rate.at(iFM).at(iMVA).at(iEff).second->SetBinContent(i, -99);
  }
  
  h_rate.at(iFM).at(iMVA).at(iEff).first ->SetLineWidth(2);
  h_rate.at(iFM).at(iMVA).at(iEff).second->SetLineWidth(2);
  h_rate.at(iFM).at(iMVA).at(iEff).first ->SetLineColor( (iMVA == 0 ? 1 : iFM+2) );
  h_rate.at(iFM).at(iMVA).at(iEff).second->SetLineColor( (iMVA == 0 ? 1 : iFM+2) );
  
} // End BookRateHist()

		 
void LoopOverEvents( TChain* chain, const TString ft_name, const TString trgPt, const int iFM, const TString tr_te,
		     std::vector<std::tuple<TString, float>>& MVAs,
		     std::vector< std::vector< std::pair<TH2D*, TH2D*> > >& h_pt,
		     std::vector< std::vector<TH1D*> >& h_count ) {

  // Get GEN branches from the factories
  float GEN_pt_br;
  float GEN_eta_br;
  float EMTF_eta_br;
  chain->SetBranchAddress("GEN_pt", &GEN_pt_br);
  chain->SetBranchAddress("GEN_eta", &GEN_eta_br);
  chain->SetBranchAddress("EMTF_eta", &EMTF_eta_br);

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
      TRG_pt += BIT;

      double GEN_pt = double(GEN_pt_br);
      double GEN_eta = double(GEN_eta_br);
      double EMTF_eta = double(EMTF_eta_br);

      if ( fabs(EMTF_eta) < ETAMIN ) continue;
      if ( fabs(EMTF_eta) > ETAMAX ) continue;

      // Fill counts from ZeroBias events
      if (GEN_eta < -10 && tr_te.Contains("test")) {
	h_count.at(iFM).at(iMVA)->Fill( TRG_pt );
      }
      
      if ( TRG_pt < 0 ) continue;
      if ( fabs(GEN_eta) < ETAMIN ) continue;
      if ( fabs(GEN_eta) > ETAMAX ) continue;
      TRG_pt = min(PTMAX - BIT, TRG_pt);
      GEN_pt = min(PTMAX - BIT, GEN_pt);

      if (tr_te.Contains("train"))
	h_pt.at(iFM).at(iMVA).first ->Fill( GEN_pt, TRG_pt );
      else if (tr_te.Contains("test"))
	h_pt.at(iFM).at(iMVA).second->Fill( GEN_pt, TRG_pt );
      else {
	std::cout << "tr_te = " << tr_te << ", not train or test. Exiting." << std::endl;
	break;
      }

    } // End loop: for (int iMVA = 0; iMVA < MVAs.size(); iMVA++)
  } // End loop: for (int iEvt = 0; iEvt < chain->GetEntries(); iEvt++)
  std::cout << "******* Leaving the " << ft_name << " " << tr_te << " event loop *******" << std::endl;

} // End function: void LoopOverEvents()

