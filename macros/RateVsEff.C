
////////////////////////////////////////////////////////
///         Macro to plot the rate reduction         ///
///         at the 90% efficiency threshold          ///
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

#include "../interface/RateVsEff.h"   // Function declarations

#include "../configs/RateVsEff/Standard.h"  // Settings that are not likely to change
#include "../configs/RateVsEff/General.h"   // General settings
#include "../configs/RateVsEff/User.h"      // Specific PtAlgos to consider

////////////////////////////////////
///  Main function: RateVsEff()  ///
////////////////////////////////////

void RateVsEff() {

  // Configure settings for this user
  RateVsEff_cfg::ConfigureUser( USER );

  // ALGOS: PtAlgo to test from User.h
  for (int i = 0; i < ALGOS.size(); i++) {
    PtAlgo algo = ALGOS.at(i);
    if ( !gSystem->AccessPathName( IN_DIR_NAME+"/"+algo.in_file_name ) )
      algo.in_file = TFile::Open( IN_DIR_NAME+"/"+algo.in_file_name ); // Check if file exists
    if ( !algo.in_file ) {
      std::cout << "ERROR: could not open data file " << IN_DIR_NAME+"/"+algo.in_file_name << std::endl;
      return;
    }
    ALGOS.at(i) = algo; // Update ALGOS
  }

  // Add trees from the input files to the TChain
  for (int i = 0; i < ALGOS.size(); i++) {
    PtAlgo algo = ALGOS.at(i);
    algo.train_tree = new TChain(algo.fact_name+"/TrainTree");
    algo.test_tree  = new TChain(algo.fact_name+"/TestTree");
    algo.train_tree->Add( IN_DIR_NAME+"/"+algo.in_file_name );
    algo.test_tree ->Add( IN_DIR_NAME+"/"+algo.in_file_name );
    ALGOS.at(i) = algo; // Update ALGOS
  }

  TString out_file_name = OUT_DIR_NAME+"/"+OUT_FILE_NAME+".root";
  TFile *out_file = new TFile(out_file_name, "recreate");


  for (int i = 0; i < ALGOS.size(); i++) {
    PtAlgo algo = ALGOS.at(i);
    BookPtHist( algo );
    BookEffHist( algo );
    BookCountHist( algo );
    for (int iPt = 0; iPt < TURN_ONS.size(); iPt++)
      BookTurnOnHist( algo, iPt, TURN_ONS.at(iPt) );
    for (int iEff = 0; iEff < EFF_CUTS.size(); iEff++)
      BookRateHist( algo, iEff, EFF_CUTS.at(iEff) );
    ALGOS.at(i) = algo; // Update ALGOS
  }


  for (int i = 0; i < ALGOS.size(); i++) {
    PtAlgo algo = ALGOS.at(i);

    // Fill 2D pT and counts vs. pT histograms
    LoopOverEvents( algo, "train" );
    LoopOverEvents( algo, "test" );

    // Compute 2D efficiency histograms and rate at efficiency threshold histograms
    for (int iBin = 1; iBin <= PTBINS; iBin++) { // Loop over GEN pT bins (x-axis)

      // Efficiency from previous bins
      float preff1[4] = {0, 0, 0, 0};
      float preff2[4] = {0, 0, 0, 0};
      
      for (int jBin = PTBINS; jBin >= 1; jBin--) { // Loop over trigger pT bins (y-axis)
	float num1 = algo.h_trg_vs_GEN_pt.first ->Integral(iBin, iBin, jBin, PTBINS); // Events passing trigger pT cut
	float den1 = algo.h_trg_vs_GEN_pt.first ->Integral(iBin, iBin,    1, PTBINS); // All events in GEN pT bin
	float eff1 = num1 / den1;
	float num2 = algo.h_trg_vs_GEN_pt.second->Integral(iBin, iBin, jBin, PTBINS);
	float den2 = algo.h_trg_vs_GEN_pt.second->Integral(iBin, iBin,    1, PTBINS);
	float eff2 = num2 / den2;

	algo.h_trg_vs_GEN_pt_eff.first ->SetBinContent(iBin, jBin, eff1);
	algo.h_trg_vs_GEN_pt_eff.first ->SetBinError  (iBin, jBin, eff1 * sqrt( (1/num1) + (1/den1) ) );
	algo.h_trg_vs_GEN_pt_eff.second->SetBinContent(iBin, jBin, eff2);
	algo.h_trg_vs_GEN_pt_eff.second->SetBinError  (iBin, jBin, eff2 * sqrt( (1/num2) + (1/den2) ) );

	for (int iPt = 0; iPt < TURN_ONS.size(); iPt++) {
	  float pt_cut = TURN_ONS.at(iPt)*1.0;
	  
	  if ( algo.h_trg_vs_GEN_pt_eff.first->GetYaxis()->GetBinLowEdge(jBin) > pt_cut - 2*BIT &&
	       algo.h_trg_vs_GEN_pt_eff.first->GetYaxis()->GetBinLowEdge(jBin) < pt_cut + 2*BIT ) {
	    algo.h_turn_ons.at(iPt)->SetBinContent(iBin, eff2);
	    algo.h_turn_ons.at(iPt)->SetBinError(iBin, eff2 * sqrt( (1/num2) + (1/den2) ) );
	  }
	} // End loop: for (int iPt = 0; iPt < TURN_ONS.size(); iPt++)
	
	for (int iEff = 0; iEff < EFF_CUTS.size(); iEff++) {
	  float thresh = EFF_CUTS.at(iEff)*0.01;
	  bool aboveT1 = (eff1 > thresh || preff1[0] > thresh || preff1[1] > thresh || preff1[2] > thresh || preff1[3] > thresh);
	  bool aboveT2 = (eff2 > thresh || preff2[0] > thresh || preff2[1] > thresh || preff2[2] > thresh || preff2[3] > thresh);
	  bool belowT1 = (eff1 < thresh || preff1[0] < thresh || preff1[1] < thresh || preff1[2] < thresh || preff1[3] < thresh);
	  bool belowT2 = (eff2 < thresh || preff2[0] < thresh || preff2[1] < thresh || preff2[2] < thresh || preff2[3] < thresh);
	  
	  // Find the closest efficiency over 5 trigger pT cuts
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
	  
	  // bool noRate1 = algo.h_ZB_rates.at(iEff).first ->GetBinContent(iBin) < 0;
	  // bool noRate2 = algo.h_ZB_rates.at(iEff).second->GetBinContent(iBin) < 0;
	  
	  if (aboveT1 && belowT1) // && noRate1?
	    algo.h_ZB_rates.at(iEff).first ->SetBinContent( iBin, algo.h_ZB_count->Integral(jBin+iClose1, PTBINS) );
	  if (aboveT2 && belowT2) // && noRate2?
	    algo.h_ZB_rates.at(iEff).second->SetBinContent( iBin, algo.h_ZB_count->Integral(jBin+iClose2, PTBINS) );
	} // End loop: for (int iEff = 0; iEff < EFF_CUTS.size(); iEff++)
	
	// Push back previous efficiencies
	preff1[3] = preff1[2];  preff1[2] = preff1[1];  preff1[1] = preff1[0];  preff1[0] = eff1;
	preff2[3] = preff2[2];  preff2[2] = preff2[1];  preff2[1] = preff2[0];  preff2[0] = eff2;
	
      } // End loop: for (int jBin = PTBINS; jBin >= 1; jBin--)
    } // End loop: for (int iBin = 1; iBin <= PTBINS; iBin++)
    
    ALGOS.at(i) = algo; // Update ALGOS
  } // End loop: for (int i = 0; i < ALGOS.size(); i++)

  out_file->cd();

  for (int i = 0; i < ALGOS.size(); i++) {
    PtAlgo algo = ALGOS.at(i);

    algo.h_trg_vs_GEN_pt.first ->Write();
    algo.h_trg_vs_GEN_pt.second->Write();
    algo.h_trg_vs_GEN_pt_eff.first ->Write();
    algo.h_trg_vs_GEN_pt_eff.second->Write();
    algo.h_ZB_count->Write();

    for (int j = 0; j < TURN_ONS.size(); j++)
      algo.h_turn_ons.at(j)->Write();
    
    for (int iEff = 0; iEff < EFF_CUTS.size(); iEff++) {
      algo.h_ZB_rates.at(iEff).first ->Write();
      algo.h_ZB_rates.at(iEff).second->Write();
      
      // Compute ratios w.r.t. the first algorithm
      if (i == 0) continue;

      algo.h_ZB_rates.at(iEff).first ->SetName( ((TString) algo.h_ZB_rates.at(iEff).first ->GetName()).ReplaceAll("rate", "rate_ratio") );
      algo.h_ZB_rates.at(iEff).second->SetName( ((TString) algo.h_ZB_rates.at(iEff).second->GetName()).ReplaceAll("rate", "rate_ratio") );
      algo.h_ZB_rates.at(iEff).first ->Divide( ALGOS.at(0).h_ZB_rates.at(iEff).first );
      algo.h_ZB_rates.at(iEff).second->Divide( ALGOS.at(0).h_ZB_rates.at(iEff).second);
      algo.h_ZB_rates.at(iEff).first ->Write();
      algo.h_ZB_rates.at(iEff).second->Write();
      
    } // End loop: for (int iEff = 0; iEff < EFF_CUTS.size(); iEff++)
    ALGOS.at(i) = algo; // Update ALGOS
  } // End loop: for (int i = 0; i < ALGOS.size(); i++) 
  
  out_file->Close();
  
  std::cout << "\nExiting RateVsEff()" << std::endl;
  
} // End void RateVsEff()


void BookPtHist( PtAlgo& algo ) {
  
  TString h_pt_str_tr = "h_pt_"+algo.fact_name+"_"+algo.MVA_name+"_train";
  TString h_pt_str_te = "h_pt_"+algo.fact_name+"_"+algo.MVA_name+"_test";
  algo.h_trg_vs_GEN_pt = std::make_pair( new TH2D( h_pt_str_tr, h_pt_str_tr,
						   PTBINS, PTMIN, PTMAX, PTBINS, PTMIN, PTMAX ),
					 new TH2D( h_pt_str_te,  h_pt_str_te,
						   PTBINS, PTMIN, PTMAX, PTBINS, PTMIN, PTMAX ) );

  algo.h_trg_vs_GEN_pt.first ->Sumw2();
  algo.h_trg_vs_GEN_pt.second->Sumw2();
  
  h_pt_str_tr.ReplaceAll("h_pt_", "");
  h_pt_str_te.ReplaceAll("h_pt_", "");
  h_pt_str_tr.ReplaceAll("_", " ");
  h_pt_str_te.ReplaceAll("_", " ");
  h_pt_str_tr.ReplaceAll("train", "(train)");
  h_pt_str_te.ReplaceAll("test",  "(test)");
  
  algo.h_trg_vs_GEN_pt.first ->SetTitle(h_pt_str_tr+" vs. GEN p_{T}");
  algo.h_trg_vs_GEN_pt.second->SetTitle(h_pt_str_te+" vs. GEN p_{T}");
  algo.h_trg_vs_GEN_pt.first ->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  algo.h_trg_vs_GEN_pt.second->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  algo.h_trg_vs_GEN_pt.first ->GetYaxis()->SetTitle(algo.fact_name+" "+algo.MVA_name+" (train) p_{T} (GeV)");
  algo.h_trg_vs_GEN_pt.second->GetYaxis()->SetTitle(algo.fact_name+" "+algo.MVA_name+" (test) p_{T} (GeV)");

} // End BookPtHist()
  

void BookEffHist( PtAlgo& algo ) {
  
  TString h_eff_str_tr = "h_eff_"+algo.fact_name+"_"+algo.MVA_name+"_train";
  TString h_eff_str_te = "h_eff_"+algo.fact_name+"_"+algo.MVA_name+"_test";
  algo.h_trg_vs_GEN_pt_eff =  std::make_pair( new TH2D( h_eff_str_tr, h_eff_str_tr,
							PTBINS, PTMIN, PTMAX, PTBINS, PTMIN, PTMAX ),
					      new TH2D( h_eff_str_te,  h_eff_str_te,
							PTBINS, PTMIN, PTMAX, PTBINS, PTMIN, PTMAX ) );
  
  algo.h_trg_vs_GEN_pt_eff.first ->Sumw2();
  algo.h_trg_vs_GEN_pt_eff.second->Sumw2();
  
  h_eff_str_tr.ReplaceAll("h_eff_", "");
  h_eff_str_te.ReplaceAll("h_eff_", "");
  h_eff_str_tr.ReplaceAll("_", " ");
  h_eff_str_te.ReplaceAll("_", " ");
  h_eff_str_tr.ReplaceAll("train", "(train)");
  h_eff_str_te.ReplaceAll("test",  "(test)");
  
  algo.h_trg_vs_GEN_pt_eff.first ->SetTitle(h_eff_str_tr+" vs. GEN p_{T} efficiency");
  algo.h_trg_vs_GEN_pt_eff.second->SetTitle(h_eff_str_te+" vs. GEN p_{T} efficiency");
  algo.h_trg_vs_GEN_pt_eff.first ->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  algo.h_trg_vs_GEN_pt_eff.second->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  algo.h_trg_vs_GEN_pt_eff.first ->GetYaxis()->SetTitle(algo.fact_name+" "+algo.MVA_name+" (train) p_{T} (GeV)");
  algo.h_trg_vs_GEN_pt_eff.second->GetYaxis()->SetTitle(algo.fact_name+" "+algo.MVA_name+" (test) p_{T} (GeV)");
  
} // End BookEffHist()
  
		 
void BookCountHist( PtAlgo& algo ) {
  
  TString h_ZB_count_str = "h_ZB_count_"+algo.fact_name+"_"+algo.MVA_name;
  algo.h_ZB_count = new TH1D( h_ZB_count_str, h_ZB_count_str, PTBINS, PTMIN, PTMAX );
  
  algo.h_ZB_count->Sumw2();
  algo.h_ZB_count->SetLineWidth(2);
  algo.h_ZB_count->SetLineColor(algo.color);
  
} // End BookCountHist()

		 
void BookTurnOnHist( PtAlgo& algo, const int iPt, const int pt_cut ) {

  TString h_turn_on_str;
  h_turn_on_str.Form("h_turn_on_%s_%s_%d", algo.fact_name.Data(), algo.MVA_name.Data(), pt_cut);
  algo.h_turn_ons.push_back( new TH1D( h_turn_on_str, h_turn_on_str, PTBINS, PTMIN, PTMAX ) );
  
  algo.h_turn_ons.at(iPt)->Sumw2();
  algo.h_turn_ons.at(iPt)->SetLineWidth(2);
  algo.h_turn_ons.at(iPt)->SetLineColor(algo.color);
  
} // End BookTurnOnHist()

		 
void BookRateHist( PtAlgo& algo, const int iEff, const int eff_cut ) {
  
  TString h_ZB_rate_str_tr;
  TString h_ZB_rate_str_te;
  h_ZB_rate_str_tr.Form("h_ZB_rates_%s_%s_%d_train", algo.fact_name.Data(), algo.MVA_name.Data(), eff_cut);
  h_ZB_rate_str_te.Form("h_ZB_rates_%s_%s_%d_test",  algo.fact_name.Data(), algo.MVA_name.Data(), eff_cut);
  
  algo.h_ZB_rates.push_back( std::make_pair( new TH1D( h_ZB_rate_str_tr, h_ZB_rate_str_tr,
						       PTBINS, PTMIN, PTMAX ),
					     new TH1D( h_ZB_rate_str_te,  h_ZB_rate_str_te,
						       PTBINS, PTMIN, PTMAX ) ) );
  
  algo.h_ZB_rates.at(iEff).first ->Sumw2();
  algo.h_ZB_rates.at(iEff).second->Sumw2();
  h_ZB_rate_str_tr.Form("%s %s (train) vs. %d%s efficiency p_{T} cut (GeV)", algo.fact_name.Data(), algo.MVA_name.Data(), eff_cut, "%");
  h_ZB_rate_str_te.Form("%s %s (test) vs. %d%s efficiency p_{T} cut (GeV)",  algo.fact_name.Data(), algo.MVA_name.Data(), eff_cut, "%");
  
  algo.h_ZB_rates.at(iEff).first ->SetTitle(h_ZB_rate_str_tr);
  algo.h_ZB_rates.at(iEff).second->SetTitle(h_ZB_rate_str_te);
  algo.h_ZB_rates.at(iEff).first ->GetXaxis()->SetTitle("p_{T} threshold (GeV)");
  algo.h_ZB_rates.at(iEff).second->GetXaxis()->SetTitle("p_{T} threshold (GeV)");
  algo.h_ZB_rates.at(iEff).first ->GetYaxis()->SetTitle(algo.fact_name+" "+algo.MVA_name+" (train) rate");
  algo.h_ZB_rates.at(iEff).second->GetYaxis()->SetTitle(algo.fact_name+" "+algo.MVA_name+" (test) rate");

  for (int i = 1; i <= PTBINS; i++) {
    algo.h_ZB_rates.at(iEff).first ->SetBinContent(i, -99);
    algo.h_ZB_rates.at(iEff).second->SetBinContent(i, -99);
  }
  
  algo.h_ZB_rates.at(iEff).first ->SetLineWidth(2);
  algo.h_ZB_rates.at(iEff).second->SetLineWidth(2);
  algo.h_ZB_rates.at(iEff).first ->SetLineColor(algo.color);
  algo.h_ZB_rates.at(iEff).second->SetLineColor(algo.color);
  
} // End BookRateHist()

		 
void LoopOverEvents( PtAlgo& algo, const TString tr_te ) {

  // Get GEN branches from the factories
  float GEN_pt_br;
  float GEN_eta_br;
  float TRK_eta_br;
  float TRK_mode_br;
  float TRK_mode_CSC_br;
  float TRK_mode_RPC_br;

  TChain* chain(0);
  if (tr_te == "train")
    chain = algo.train_tree;
  else if (tr_te == "test")
    chain = algo.test_tree;
  else {
    std::cout << "tr_te = " << tr_te << ", not train or test. Exiting." << std::endl;
    return;
  }

  bool isEMTF = algo.MVA_name.Contains("EMTF");

  chain->SetBranchAddress("GEN_pt", &GEN_pt_br);
  chain->SetBranchAddress("GEN_eta", &GEN_eta_br);
  if (isEMTF) {
    chain->SetBranchAddress("EMTF_eta", &TRK_eta_br);
    chain->SetBranchAddress("EMTF_mode", &TRK_mode_br);
    chain->SetBranchAddress("EMTF_mode_CSC", &TRK_mode_CSC_br);
    chain->SetBranchAddress("EMTF_mode_RPC", &TRK_mode_RPC_br);
  } else {
    chain->SetBranchAddress("TRK_eta", &TRK_eta_br);
    chain->SetBranchAddress("TRK_mode", &TRK_mode_br);
    chain->SetBranchAddress("TRK_mode_CSC", &TRK_mode_CSC_br);
    chain->SetBranchAddress("TRK_mode_RPC", &TRK_mode_RPC_br);
  }

  // Get trigger branch from the factory
  chain->SetBranchAddress( algo.MVA_name, &(algo.MVA_val) );
  
  std::cout << "\n******* About to enter the " << algo.fact_name << " " << tr_te << " event loop *******" << std::endl;
  for (int iEvt = 0; iEvt < chain->GetEntries(); iEvt++) {
    
    if (iEvt > MAX_EVT && MAX_EVT > 0) break;
    if ( (iEvt % REPORT_EVT) == 0 ) std::cout << "*** Looking at event " << iEvt << " ***" << std::endl;
    
    chain->GetEntry(iEvt);
    
    // // Loop over different trigger pT computations (should maybe keep? - AWB 21.04.17)
    // for (int iMVA = 0; iMVA < MVAs.size(); iMVA++) {
    
    // Access trigger pT
    double TRG_pt = algo.MVA_val;
    if ( not isEMTF ) {
      if ( algo.fact_name.Contains("invPtTarg") )
	TRG_pt = 1. / max(0.001, TRG_pt); // Protect against negative 1/pT values
      if ( algo.fact_name.Contains("logPtTarg") )
	TRG_pt = pow(2, TRG_pt);
    }
    TRG_pt *= algo.trg_pt_scale;
    TRG_pt += BIT;

    double GEN_pt  = double(GEN_pt_br);
    double GEN_eta = double(GEN_eta_br);
    double TRK_eta = double(TRK_eta_br);
    int TRK_mode     =  int(TRK_mode_br);
    int TRK_mode_CSC =  int(TRK_mode_CSC_br);
    int TRK_mode_RPC =  int(TRK_mode_RPC_br);
    
    if ( fabs(TRK_eta) < ETAMIN ) continue;
    if ( fabs(TRK_eta) > ETAMAX ) continue;

    bool good_mode     = false;
    bool good_mode_CSC = false;
    bool good_mode_RPC = false;
    
    for (int iMode = 0; iMode < algo.modes.size(); iMode++)
      if (algo.modes.at(iMode) == TRK_mode)
	good_mode = true;
    for (int iMode = 0; iMode < algo.modes_CSC.size(); iMode++)
      if (algo.modes_CSC.at(iMode) == TRK_mode_CSC)
	good_mode_CSC = true;
    for (int iMode = 0; iMode < algo.modes_RPC.size(); iMode++)
      if (algo.modes_RPC.at(iMode) == TRK_mode_RPC)
	good_mode_RPC = true;

    // Only use events with the proper modes
    if (not (good_mode && good_mode_CSC && good_mode_RPC) ) continue;

    // Fill counts from ZeroBias events
    if (GEN_eta < -10 && tr_te.Contains("test")) {
      algo.h_ZB_count->Fill( TRG_pt );
    }
      
    if ( TRG_pt < 0 ) continue;
    if ( fabs(GEN_eta) < ETAMIN ) continue;
    if ( fabs(GEN_eta) > ETAMAX ) continue;
    TRG_pt = min(PTMAX - BIT, TRG_pt);
    GEN_pt = min(PTMAX - BIT, GEN_pt);
    
    if (tr_te == "train")
      algo.h_trg_vs_GEN_pt.first ->Fill( GEN_pt, TRG_pt );
    else if (tr_te == "test")
      algo.h_trg_vs_GEN_pt.second->Fill( GEN_pt, TRG_pt );
    
  } // End loop: for (int iEvt = 0; iEvt < chain->GetEntries(); iEvt++)
  std::cout << "******* Leaving the " << algo.fact_name << " " << tr_te << " event loop *******" << std::endl;

} // End function: void LoopOverEvents()

