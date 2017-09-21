
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
    if (algo.unique_ID == "") algo.unique_ID = algo.fact_name+"_"+algo.MVA_name;
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
    BookChargeHist( algo );
    for (int iPt = 0; iPt < TURN_ONS.size(); iPt++) {
      BookTurnOnHist( algo, iPt, TURN_ONS.at(iPt) );
      BookEtaHist( algo, iPt, TURN_ONS.at(iPt) );
    }
    BookEtaHist( algo, TURN_ONS.size(), 0 );
    for (int iEff = 0; iEff < EFF_CUTS.size(); iEff++) {
      BookPtScaleHist( algo, iEff, EFF_CUTS.at(iEff) );
      BookRateHist( algo, iEff, EFF_CUTS.at(iEff) );
    }
    ALGOS.at(i) = algo; // Update ALGOS
  }


  for (int i = 0; i < ALGOS.size(); i++) {
    PtAlgo algo = ALGOS.at(i);

    // Fill 2D pT and counts vs. pT histograms
    LoopOverEvents( algo, "train" );
    LoopOverEvents( algo, "test" );

    // Compute 2D efficiency histograms and rate at efficiency threshold histograms
    int jBinPrev1 = 1; // Save threshold for previous (lower-pT) bin
    int jBinPrev2 = 1;
    for (int iBin = 1; iBin <= PTBINS; iBin++) { // Loop over GEN pT bins (x-axis)
      int   nBinsEx = 0;  // Number of bins on either side to consider, to add stats
      float iBinMin = PTMIN + (iBin - 1) * (PTMAX - PTMIN) / PTBINS;
      float iBinMax = PTMIN +  iBin      * (PTMAX - PTMIN) / PTBINS;
      iBinMin = fmax(iBinMin, BIT);

      float den1 = fmax(algo.h_trg_vs_GEN_pt.first ->Integral(iBin, iBin, 1, PTBINS*PTDIVS), BIT); // All events in GEN pT bin
      float den2 = fmax(algo.h_trg_vs_GEN_pt.second->Integral(iBin, iBin, 1, PTBINS*PTDIVS), BIT);

      // Get enough stats in denominator to achieve ~2% uncertainty
      while(den1 < pow(1./ERRMIN, 2) && den2 < pow(1./ERRMIN, 2) && iBin-nBinsEx > 1 && iBin+nBinsEx < PTBINS) {
	nBinsEx += 1;
	iBinMin = PTMIN + (iBin-nBinsEx - 1) * (PTMAX - PTMIN) / PTBINS;
	iBinMax = PTMIN + (iBin+nBinsEx    ) * (PTMAX - PTMIN) / PTBINS;
	iBinMin = fmax(iBinMin-nBinsEx, BIT);
	den1 = fmax(algo.h_trg_vs_GEN_pt.first ->Integral(iBin-nBinsEx, iBin+nBinsEx, 1, PTBINS*PTDIVS), BIT); // All events in GEN pT bin
	den2 = fmax(algo.h_trg_vs_GEN_pt.second->Integral(iBin-nBinsEx, iBin+nBinsEx, 1, PTBINS*PTDIVS), BIT);
      }


      // Fill EMTF charge-efficiency plot
      if (algo.MVA_name.Contains("EMTF")) {
	algo.h_charge_eff->SetBinContent(iBin, algo.h_charge_eff->GetBinContent(iBin) / (den1 + den2));
	algo.h_charge_eff->SetBinError  (iBin, 1. / sqrt(den1 + den2));
      }

      // Fill trigger pT scaling and ZeroBias rate plots
      for (int jBin = PTBINS*PTDIVS; jBin >= 1; jBin--) { // Loop over trigger pT bins (y-axis)
	float jBinMin = PTMIN + (jBin - 1) * (PTMAX - PTMIN) / (PTBINS*PTDIVS);
	float jBinMax = PTMIN +  jBin      * (PTMAX - PTMIN) / (PTBINS*PTDIVS);
	jBinMin = fmax(jBinMin, BIT);

	if (jBinMin > iBinMin + BIT) continue; // Don't set threshold higher than nominal pT

	float num1 = algo.h_trg_vs_GEN_pt.first ->Integral(iBin-nBinsEx, iBin+nBinsEx, jBin, PTBINS*PTDIVS); // Events passing trigger pT cut
	float eff1 = num1 / den1;
	float num2 = algo.h_trg_vs_GEN_pt.second->Integral(iBin-nBinsEx, iBin+nBinsEx, jBin, PTBINS*PTDIVS);
	float eff2 = num2 / den2;

	algo.h_trg_eff_vs_GEN_pt.first ->SetBinContent(iBin, jBin, eff1);
	algo.h_trg_eff_vs_GEN_pt.first ->SetBinError  (iBin, jBin, eff1 / sqrt(den1) );
	algo.h_trg_eff_vs_GEN_pt.second->SetBinContent(iBin, jBin, eff2);
	algo.h_trg_eff_vs_GEN_pt.second->SetBinError  (iBin, jBin, eff2 / sqrt(den2) );

	// Loop over efficiency-at-turn-on working points
	for (int iEff = 0; iEff < EFF_CUTS.size(); iEff++) {
	  float thresh = EFF_CUTS.at(iEff)*0.01;

	  // Fill if efficiency reaches working point
	  if ( (eff1 >= thresh || jBin == 1) && algo.h_pt_scales.at(iEff).first ->GetBinContent(iBin) < 0) {
	    algo.h_ZB_rates.at(iEff).first  ->SetBinContent( iBin, algo.h_ZB_count->Integral(jBin, PTBINS*PTDIVS) / (2*nBinsEx + 1) );
	    algo.h_pt_scales.at(iEff).first ->SetBinContent( iBin, iBinMin / jBinMin );
	    jBinPrev1 = jBin;
	  }
	  if ( (eff2 >= thresh || jBin == 1) && algo.h_pt_scales.at(iEff).second->GetBinContent(iBin) < 0) {
	    algo.h_ZB_rates.at(iEff).second ->SetBinContent( iBin, algo.h_ZB_count->Integral(jBin, PTBINS*PTDIVS) / (2*nBinsEx + 1) );
	    algo.h_pt_scales.at(iEff).second->SetBinContent( iBin, iBinMin / jBinMin );
	    jBinPrev2 = jBin;
	  }

	} // End loop: for (int iEff = 0; iEff < EFF_CUTS.size(); iEff++)
	
      } // End loop: for (int jBin = PTBINS*PTDIVS; jBin >= 1; jBin--)


      // Loop over pT thresholds for turn-on curves
      for (int iPt = 0; iPt < TURN_ONS.size(); iPt++) {
	assert(algo.h_pt_scales.at(0).first ->GetBinContent(iBin) > 0);
	assert(algo.h_pt_scales.at(0).second->GetBinContent(iBin) > 0);

	float pt_cut1 = 1.0 * TURN_ONS.at(iPt) / algo.h_pt_scales.at(0).first ->GetBinContent(iBin);
	float pt_cut2 = 1.0 * TURN_ONS.at(iPt) / algo.h_pt_scales.at(0).second->GetBinContent(iBin);
	
	int jBin1 = int( PTBINS*PTDIVS * fmin(1.0, fmax(0.0, (pt_cut1 - PTMIN) / (PTMAX - PTMIN))) );
	int jBin2 = int( PTBINS*PTDIVS * fmin(1.0, fmax(0.0, (pt_cut2 - PTMIN) / (PTMAX - PTMIN))) );

	float num1 = algo.h_trg_vs_GEN_pt.first ->Integral(iBin-nBinsEx, iBin+nBinsEx, jBin1, PTBINS*PTDIVS); // Events passing trigger pT cut
	float eff1 = num1 / den1;
	float num2 = algo.h_trg_vs_GEN_pt.second->Integral(iBin-nBinsEx, iBin+nBinsEx, jBin2, PTBINS*PTDIVS);
	float eff2 = num2 / den2;

	algo.h_turn_ons.at(iPt).first ->SetBinContent(iBin, eff1);
	algo.h_turn_ons.at(iPt).first ->SetBinError(iBin, eff1 / sqrt(den1) );
	algo.h_turn_ons.at(iPt).second->SetBinContent(iBin, eff2);
	algo.h_turn_ons.at(iPt).second->SetBinError(iBin, eff2 / sqrt(den2) );

	for (int iEta = 1; iEta <= ETABINS; iEta++) {
	  algo.h_ZB_etas.at(iPt).first ->SetBinContent( iEta, algo.h_ZB_count_eta->Integral(jBin1, PTBINS*PTDIVS, iEta, iEta) / (2*nBinsEx + 1) );
	  algo.h_ZB_etas.at(iPt).second->SetBinContent( iEta, algo.h_ZB_count_eta->Integral(jBin2, PTBINS*PTDIVS, iEta, iEta) / (2*nBinsEx + 1) );

	  if (iPt == TURN_ONS.size() - 1) {  // Also fill inclusive eta plot
	    algo.h_ZB_etas.at(iPt+1).first ->SetBinContent( iEta, algo.h_ZB_count_eta->Integral(1, PTBINS*PTDIVS, iEta, iEta) / (2*nBinsEx + 1) );
	    algo.h_ZB_etas.at(iPt+1).second->SetBinContent( iEta, algo.h_ZB_count_eta->Integral(1, PTBINS*PTDIVS, iEta, iEta) / (2*nBinsEx + 1) );
	  }
	}

      } // End loop: for (int iPt = 0; iPt < TURN_ONS.size(); iPt++)

    } // End loop: for (int iBin = 1; iBin <= PTBINS; iBin++)
    
    ALGOS.at(i) = algo; // Update ALGOS
  } // End loop: for (int i = 0; i < ALGOS.size(); i++)
  
  out_file->cd();

  for (int i = 0; i < ALGOS.size(); i++) {
    PtAlgo algo = ALGOS.at(i);

    algo.h_trg_vs_GEN_pt.first ->Write();
    algo.h_trg_vs_GEN_pt.second->Write();
    algo.h_trg_eff_vs_GEN_pt.first ->Write();
    algo.h_trg_eff_vs_GEN_pt.second->Write();
    algo.h_ZB_count->Write();
    algo.h_ZB_count_eta->Write();
    if (algo.MVA_name.Contains("EMTF"))
      algo.h_charge_eff->Write();

    for (int j = 0; j < TURN_ONS.size(); j++) {
      algo.h_turn_ons.at(j).first ->Write();
      algo.h_turn_ons.at(j).second->Write();
      algo.h_ZB_etas.at(j).first ->Write();
      algo.h_ZB_etas.at(j).second->Write();
    }
    algo.h_ZB_etas.at(TURN_ONS.size()).first ->Write();
    algo.h_ZB_etas.at(TURN_ONS.size()).second->Write();
    
    for (int iEff = 0; iEff < EFF_CUTS.size(); iEff++) {
      algo.h_pt_scales.at(iEff).first ->Write();
      algo.h_pt_scales.at(iEff).second->Write();

      assert(algo.h_ZB_rates.at(iEff).first ->GetBinContent(1) == algo.h_ZB_count->Integral());
      assert(algo.h_ZB_rates.at(iEff).first ->GetBinContent(1) == algo.h_ZB_count->Integral());
      algo.h_ZB_rates.at(iEff).first ->Scale(1.0 / algo.h_ZB_rates.at(iEff).first ->GetBinContent(1));
      algo.h_ZB_rates.at(iEff).second->Scale(1.0 / algo.h_ZB_rates.at(iEff).second->GetBinContent(1));

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
  
  TString h_pt_str_tr = "h_trg_vs_GEN_pt_"+algo.unique_ID+"_train";
  TString h_pt_str_te = "h_trg_vs_GEN_pt_"+algo.unique_ID+"_test";
  algo.h_trg_vs_GEN_pt = std::make_pair( new TH2D( h_pt_str_tr, h_pt_str_tr,
						   PTBINS, PTMIN, PTMAX, PTBINS*PTDIVS, PTMIN, PTMAX ),
					 new TH2D( h_pt_str_te,  h_pt_str_te,
						   PTBINS, PTMIN, PTMAX, PTBINS*PTDIVS, PTMIN, PTMAX ) );

  algo.h_trg_vs_GEN_pt.first ->Sumw2();
  algo.h_trg_vs_GEN_pt.second->Sumw2();
  
  algo.h_trg_vs_GEN_pt.first ->SetTitle(algo.alias+" (train) vs. GEN p_{T}");
  algo.h_trg_vs_GEN_pt.second->SetTitle(algo.alias+" (test) vs. GEN p_{T}");
  algo.h_trg_vs_GEN_pt.first ->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  algo.h_trg_vs_GEN_pt.second->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  algo.h_trg_vs_GEN_pt.first ->GetYaxis()->SetTitle("Unscaled trigger p_{T} (GeV)");
  algo.h_trg_vs_GEN_pt.second->GetYaxis()->SetTitle("Unscaled trigger p_{T} (GeV)");

} // End BookPtHist()
  

void BookEffHist( PtAlgo& algo ) {
  
  TString h_eff_str_tr = "h_trg_eff_vs_GEN_pt_"+algo.unique_ID+"_train";
  TString h_eff_str_te = "h_trg_eff_vs_GEN_pt_"+algo.unique_ID+"_test";
  algo.h_trg_eff_vs_GEN_pt =  std::make_pair( new TH2D( h_eff_str_tr, h_eff_str_tr,
							PTBINS, PTMIN, PTMAX, PTBINS*PTDIVS, PTMIN, PTMAX ),
					      new TH2D( h_eff_str_te,  h_eff_str_te,
							PTBINS, PTMIN, PTMAX, PTBINS*PTDIVS, PTMIN, PTMAX ) );
  
  algo.h_trg_eff_vs_GEN_pt.first ->Sumw2();
  algo.h_trg_eff_vs_GEN_pt.second->Sumw2();
  
  algo.h_trg_eff_vs_GEN_pt.first ->SetTitle(algo.alias+" (train) vs. GEN p_{T} efficiency");
  algo.h_trg_eff_vs_GEN_pt.second->SetTitle(algo.alias+" (test) vs. GEN p_{T} efficiency");
  algo.h_trg_eff_vs_GEN_pt.first ->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  algo.h_trg_eff_vs_GEN_pt.second->GetXaxis()->SetTitle("GEN p_{T} (GeV)");
  algo.h_trg_eff_vs_GEN_pt.first ->GetYaxis()->SetTitle("Unscaled trigger p_{T} (GeV)");
  algo.h_trg_eff_vs_GEN_pt.second->GetYaxis()->SetTitle("Unscaled trigger p_{T} (GeV)");
  
} // End BookEffHist()
  
		 
void BookCountHist( PtAlgo& algo ) {
  
  TString h_ZB_count_str = "h_ZB_count_"+algo.unique_ID;
  TString h_ZB_count_eta_str = "h_ZB_count_eta_"+algo.unique_ID;
  algo.h_ZB_count     = new TH1D( h_ZB_count_str,     h_ZB_count_str,     PTBINS*PTDIVS, PTMIN, PTMAX );
  algo.h_ZB_count_eta = new TH2D( h_ZB_count_eta_str, h_ZB_count_eta_str, PTBINS*PTDIVS, PTMIN, PTMAX, ETABINS, ETAMIN, ETAMAX );
  
  algo.h_ZB_count->Sumw2();
  algo.h_ZB_count->SetLineWidth(2);
  algo.h_ZB_count->SetLineColor(algo.color);

  algo.h_ZB_count_eta->Sumw2();
  algo.h_ZB_count_eta->SetLineWidth(2);
  algo.h_ZB_count_eta->SetLineColor(algo.color);

  algo.h_ZB_count->SetTitle(algo.alias+" events vs. p_{T} threshold (GeV)") ;
  algo.h_ZB_count->GetXaxis()->SetTitle("Unscaled trigger p_{T} threshold (GeV)") ;
  algo.h_ZB_count->GetYaxis()->SetTitle("ZeroBias events") ;
  
  algo.h_ZB_count_eta->SetTitle(algo.alias+" events vs. p_{T} threshold (GeV) vs. muon |#eta|") ;
  algo.h_ZB_count_eta->GetXaxis()->SetTitle("Unscaled trigger p_{T} threshold (GeV)") ;
  algo.h_ZB_count_eta->GetYaxis()->SetTitle("Muon |#eta|");
  
} // End BookCountHist()

		 
void BookChargeHist( PtAlgo& algo ) {

  if (!algo.MVA_name.Contains("EMTF")) return; // Only plot for EMTF charge assignment
  
  TString h_charge_eff_str = "h_charge_eff_"+algo.unique_ID;
  algo.h_charge_eff = new TH1D( h_charge_eff_str, h_charge_eff_str, PTBINS, PTMIN, PTMAX );
  
  algo.h_charge_eff->Sumw2();
  algo.h_charge_eff->SetLineWidth(2);
  algo.h_charge_eff->SetLineColor(algo.color);

  algo.h_charge_eff->SetTitle(algo.alias+" correct-charge efficiency vs. GEN p_{T} (GeV)") ;
  algo.h_charge_eff->GetXaxis()->SetTitle("GEN muon p_{T} (GeV)") ;
  algo.h_charge_eff->GetYaxis()->SetTitle("Efficiency") ;
  
} // End BookChargeHist()

		 
void BookTurnOnHist( PtAlgo& algo, const int iPt, const int pt_cut ) {

  TString h_turn_on_str_tr;
  TString h_turn_on_str_te;
  h_turn_on_str_tr.Form("h_turn_on_%s_%d_train", algo.unique_ID.Data(), pt_cut);
  h_turn_on_str_te.Form("h_turn_on_%s_%d_test",  algo.unique_ID.Data(), pt_cut);
  algo.h_turn_ons.push_back( std::make_pair( new TH1D( h_turn_on_str_tr, h_turn_on_str_tr, PTBINS, PTMIN, PTMAX ),
					     new TH1D( h_turn_on_str_te, h_turn_on_str_te, PTBINS, PTMIN, PTMAX ) ) );
  
  algo.h_turn_ons.at(iPt).first ->Sumw2();
  algo.h_turn_ons.at(iPt).second->Sumw2();
  algo.h_turn_ons.at(iPt).first ->SetLineWidth(2);
  algo.h_turn_ons.at(iPt).second->SetLineWidth(2);
  algo.h_turn_ons.at(iPt).first ->SetLineColor(algo.color);
  algo.h_turn_ons.at(iPt).second->SetLineColor(algo.color);

  h_turn_on_str_tr.Form("%s efficiency at %d GeV threshold (train)", algo.alias.Data(), pt_cut);
  h_turn_on_str_te.Form("%s efficiency at %d GeV threshold (test)",  algo.alias.Data(), pt_cut);
  
  algo.h_turn_ons.at(iPt).first ->SetTitle(h_turn_on_str_tr);
  algo.h_turn_ons.at(iPt).second->SetTitle(h_turn_on_str_te);
  algo.h_turn_ons.at(iPt).first ->GetXaxis()->SetTitle("GEN p_{T} (GeV)") ;
  algo.h_turn_ons.at(iPt).second->GetXaxis()->SetTitle("GEN p_{T} (GeV)") ;
  algo.h_turn_ons.at(iPt).first ->GetXaxis()->SetTitle("Efficiency") ;
  algo.h_turn_ons.at(iPt).second->GetXaxis()->SetTitle("Efficiency") ;
  
} // End BookTurnOnHist()

		 
void BookEtaHist( PtAlgo& algo, const int iPt, const int pt_cut ) {

  TString h_ZB_eta_str_tr;
  TString h_ZB_eta_str_te;
  h_ZB_eta_str_tr.Form("h_ZB_eta_%s_%d_train", algo.unique_ID.Data(), pt_cut);
  h_ZB_eta_str_te.Form("h_ZB_eta_%s_%d_test",  algo.unique_ID.Data(), pt_cut);
  algo.h_ZB_etas.push_back( std::make_pair( new TH1D( h_ZB_eta_str_tr, h_ZB_eta_str_tr, ETABINS, ETAMIN, ETAMAX ),
					    new TH1D( h_ZB_eta_str_te, h_ZB_eta_str_te, ETABINS, ETAMIN, ETAMAX ) ) );
  
  algo.h_ZB_etas.at(iPt).first ->Sumw2();
  algo.h_ZB_etas.at(iPt).second->Sumw2();
  algo.h_ZB_etas.at(iPt).first ->SetLineWidth(2);
  algo.h_ZB_etas.at(iPt).second->SetLineWidth(2);
  algo.h_ZB_etas.at(iPt).first ->SetLineColor(algo.color);
  algo.h_ZB_etas.at(iPt).second->SetLineColor(algo.color);

  h_ZB_eta_str_tr.Form("%s ZeroBias |#eta| at %d GeV threshold (train)", algo.alias.Data(), pt_cut);
  h_ZB_eta_str_te.Form("%s ZeroBias |#eta| at %d GeV threshold (test)",  algo.alias.Data(), pt_cut);
  
  algo.h_ZB_etas.at(iPt).first ->SetTitle(h_ZB_eta_str_tr);
  algo.h_ZB_etas.at(iPt).second->SetTitle(h_ZB_eta_str_te);
  algo.h_ZB_etas.at(iPt).first ->GetXaxis()->SetTitle("Muon |#eta|") ;
  algo.h_ZB_etas.at(iPt).second->GetXaxis()->SetTitle("Muon |#eta|") ;
  algo.h_ZB_etas.at(iPt).first ->GetXaxis()->SetTitle("Counts") ;
  algo.h_ZB_etas.at(iPt).second->GetXaxis()->SetTitle("Counts") ;
  
} // End BookEtaHist()

		 
void BookPtScaleHist( PtAlgo& algo, const int iEff, const int eff_cut ) {
  
  TString h_pt_scale_str_tr;
  TString h_pt_scale_str_te;
  h_pt_scale_str_tr.Form("h_pt_scale_%s_%d_train", algo.unique_ID.Data(), eff_cut);
  h_pt_scale_str_te.Form("h_pt_scale_%s_%d_test",  algo.unique_ID.Data(), eff_cut);
  
  algo.h_pt_scales.push_back( std::make_pair( new TH1D( h_pt_scale_str_tr, h_pt_scale_str_tr,
							PTBINS, PTMIN, PTMAX ),
					      new TH1D( h_pt_scale_str_te,  h_pt_scale_str_te,
							PTBINS, PTMIN, PTMAX ) ) );
  
  algo.h_pt_scales.at(iEff).first ->Sumw2();
  algo.h_pt_scales.at(iEff).second->Sumw2();
  algo.h_pt_scales.at(iEff).first ->SetLineWidth(2);
  algo.h_pt_scales.at(iEff).second->SetLineWidth(2);
  algo.h_pt_scales.at(iEff).first ->SetLineColor(algo.color);
  algo.h_pt_scales.at(iEff).second->SetLineColor(algo.color);
  
  h_pt_scale_str_tr.Form("%s (train) %d%s efficiency scaling vs. p_{T} (GeV)", algo.alias.Data(), eff_cut, "%");
  h_pt_scale_str_te.Form("%s (test) %d%s efficiency scaling vs. p_{T} (GeV)",  algo.alias.Data(), eff_cut, "%");
  
  algo.h_pt_scales.at(iEff).first ->SetTitle(h_pt_scale_str_tr);
  algo.h_pt_scales.at(iEff).second->SetTitle(h_pt_scale_str_te);
  algo.h_pt_scales.at(iEff).first ->GetXaxis()->SetTitle("p_{T} threshold (GeV)");
  algo.h_pt_scales.at(iEff).second->GetXaxis()->SetTitle("p_{T} threshold (GeV)");
  algo.h_pt_scales.at(iEff).first ->GetYaxis()->SetTitle("Trigger p_{T} scaling");
  algo.h_pt_scales.at(iEff).second->GetYaxis()->SetTitle("Trigger p_{T} scaling");
  
  for (int i = 1; i <= PTBINS; i++) {
    algo.h_pt_scales.at(iEff).first ->SetBinContent(i, -99);
    algo.h_pt_scales.at(iEff).second->SetBinContent(i, -99);
  }

} // End BookPtScaleHist()

		 
void BookRateHist( PtAlgo& algo, const int iEff, const int eff_cut ) {
  
  TString h_ZB_rate_str_tr;
  TString h_ZB_rate_str_te;
  h_ZB_rate_str_tr.Form("h_ZB_rate_%s_%d_train", algo.unique_ID.Data(), eff_cut);
  h_ZB_rate_str_te.Form("h_ZB_rate_%s_%d_test",  algo.unique_ID.Data(), eff_cut);
  
  algo.h_ZB_rates.push_back( std::make_pair( new TH1D( h_ZB_rate_str_tr, h_ZB_rate_str_tr,
						       PTBINS, PTMIN, PTMAX ),
					     new TH1D( h_ZB_rate_str_te,  h_ZB_rate_str_te,
						       PTBINS, PTMIN, PTMAX ) ) );
  
  algo.h_ZB_rates.at(iEff).first ->Sumw2();
  algo.h_ZB_rates.at(iEff).second->Sumw2();
  h_ZB_rate_str_tr.Form("%s rate (train) vs. %d%s efficiency p_{T} cut (GeV)", algo.alias.Data(), eff_cut, "%");
  h_ZB_rate_str_te.Form("%s rate (test) vs. %d%s efficiency p_{T} cut (GeV)",  algo.alias.Data(), eff_cut, "%");
  
  algo.h_ZB_rates.at(iEff).first ->SetTitle(h_ZB_rate_str_tr);
  algo.h_ZB_rates.at(iEff).second->SetTitle(h_ZB_rate_str_te);
  algo.h_ZB_rates.at(iEff).first ->GetXaxis()->SetTitle("p_{T} threshold (GeV)");
  algo.h_ZB_rates.at(iEff).second->GetXaxis()->SetTitle("p_{T} threshold (GeV)");
  algo.h_ZB_rates.at(iEff).first ->GetYaxis()->SetTitle("Fraction of inclusive rate");
  algo.h_ZB_rates.at(iEff).second->GetYaxis()->SetTitle("Fraction of inclusive rate");

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
  float GEN_charge_br;
  float EMTF_charge_br;
  float EMTF_mode_br;
  float EMTF_mode_CSC_br;
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
  chain->SetBranchAddress("GEN_charge", &GEN_charge_br);
  chain->SetBranchAddress("EMTF_charge", &EMTF_charge_br);
  chain->SetBranchAddress("EMTF_mode", &EMTF_mode_br);
  chain->SetBranchAddress("EMTF_mode_CSC", &EMTF_mode_CSC_br);
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
  
  std::cout << "\n******* About to enter the " << algo.fact_name << " (" << algo.unique_ID << ") " << tr_te << " event loop *******" << std::endl;
  for (int iEvt = 0; iEvt < chain->GetEntries(); iEvt++) {
    
    if (iEvt > MAX_EVT && MAX_EVT > 0) break;
    if ( (iEvt % REPORT_EVT) == 0 ) std::cout << "*** Looking at event " << iEvt << " ***" << std::endl;
    
    chain->GetEntry(iEvt);
    
    // // Loop over different trigger pT computations (should maybe keep? - AWB 21.04.17)
    // for (int iMVA = 0; iMVA < MVAs.size(); iMVA++) {
    
    // Access trigger pT
    double TRG_pt = algo.MVA_val;
    if ( not isEMTF ) {
      if ( algo.fact_name.Contains("ptTarg") )
	TRG_pt = TRG_pt;
      if ( algo.fact_name.Contains("invPtTarg") )
	TRG_pt = 1. / fmax(0.001, TRG_pt); // Protect against negative 1/pT values
      if ( algo.fact_name.Contains("logPtTarg") )
	TRG_pt = pow(2, TRG_pt);
      if ( algo.fact_name.Contains("sqrtPtTarg") )
	TRG_pt = pow(TRG_pt, 2);
    }
    TRG_pt *= algo.trg_pt_scale;
    TRG_pt += BIT; // Small value to offset EMTF trigger pT right on 0.0/0.5 GeV boundaries 

    double GEN_pt     = double(GEN_pt_br);
    double GEN_eta    = double(GEN_eta_br);
    double TRK_eta    = double(TRK_eta_br);
    int GEN_charge    = int(GEN_charge_br);
    int EMTF_charge   = int(EMTF_charge_br);
    int EMTF_mode     = int(EMTF_mode_br);
    int EMTF_mode_CSC = int(EMTF_mode_CSC_br);
    int TRK_mode      = int(TRK_mode_br);
    int TRK_mode_CSC  = int(TRK_mode_CSC_br);
    int TRK_mode_RPC  = int(TRK_mode_RPC_br);

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

    // Impose range on trigger pT
    TRG_pt = fmin(PTMAX - BIT, fmax(PTMIN + BIT, TRG_pt));
    assert(TRG_pt > PTMIN);
    assert(TRG_pt < PTMAX);

    // Fill counts from ZeroBias events
    if (GEN_eta < -10 && tr_te.Contains("test")) {
      if ( isEMTF || (!algo.match_EMTF) || (TRK_mode == EMTF_mode && TRK_mode_CSC == EMTF_mode_CSC) ) {
	algo.h_ZB_count    ->Fill( TRG_pt );
	algo.h_ZB_count_eta->Fill( TRG_pt, fabs(TRK_eta) );
      }
    }

    if ( GEN_pt < PTMIN + BIT ) continue;
    if ( GEN_pt > PTMAX - BIT ) continue;
    if ( fabs(GEN_eta) < ETAMIN ) continue;
    if ( fabs(GEN_eta) > ETAMAX ) continue;

    if (tr_te == "train")
      algo.h_trg_vs_GEN_pt.first ->Fill( GEN_pt, TRG_pt );
    else if (tr_te == "test")
      algo.h_trg_vs_GEN_pt.second->Fill( GEN_pt, TRG_pt );

    if (isEMTF && EMTF_charge == GEN_charge)
      algo.h_charge_eff->Fill( GEN_pt );
    
  } // End loop: for (int iEvt = 0; iEvt < chain->GetEntries(); iEvt++)
  std::cout << "\n******* Leaving the " << algo.fact_name << " (" << algo.unique_ID << ") " << tr_te << " event loop *******" << std::endl;

} // End function: void LoopOverEvents()

