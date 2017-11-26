#include <iostream>
#include <iomanip>
using namespace std;
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TAttFill.h"
#include "TCanvas.h"
#include <vector>
#include "stdio.h"
#include <stdlib.h>
#include "math.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "THStack.h"
#include "TFitResultPtr.h"

//****************************************************************************************
//*Study the pT classifier performance and compare with BDT regression normally used at P5
//*including ROC cureve, signal-to-background ratio, etc
//*
//*Wei Shi @ Nov 11, 2017 CERN Geneva
//****************************************************************************************
//*Reference Table
//*=======================================================================================
//*         Test MC Events   #                     |  GEN < PT_CUT    |    GEN >= PT_CUT
//*======================================================================================
//*MC True: Signal (class1)                        |       NO         |         S
//*---------------------------------------------------------------------------------------
//*MC True: Background (class2)                    |       B          |         NO
//*=======================================================================================
//*Cut: class1 >= a (Signal)                       |       S1         |         S2
//*---------------------------------------------------------------------------------------
//*Cut: class1 < a (Background)                    |       B1         |         B2
//*=======================================================================================
//*In binary classifier, two classes cut "class1>=a && class2<=b" is redundant. Since it 
//*only make sense to have b<=1-a, a>=1-b, this implies "a+b=1", ie, cut on one class is 
//*sufficienct. As "a" decreases, efficiency increase as well as rate(monotonic increase)
//*Use Eff_REF as the stop point for further decreasing "a";
//*True Postive Rate: S2/(S2+B2), i.e. S2/S, signal efficiency/plateau trigger efficiency
//*False Positive Rate: S1/(S1+B1)
//*S=S2+B2
//*B=S1+B1
void ClassifierROC()
{
        //USER modify here ONLY//
        //================================================================
        Int_t PT_CUT = 32;//the classifier trained on this cut
        Float_t EFF_REF = 0.95;//the eff beyond which classifier cut stops
        Int_t Bins=10;//bins on class cut
        Int_t lxplus=1;//machine: lxplus(1) or bonner(0)?
        //================================================================
        TString Cluster="";
        if(lxplus==1){
                Cluster = "/afs/cern.ch/work/w/wshi/public/TMVA_2017/TMVA/";
        }
        else{
                Cluster = "/home/ws13/TMVA/TMVA/";//bonner
        }
        TString fileName = Cluster + "EMTFPtAssign2017/pTMulticlass_MODE_15_bitCompr_RPC_" + Form("%d", PT_CUT) + ".root";
        TString directoryName = "f_MODE_15_noWgt_bitCompr_RPC/TestTree";
        TFile* myFile = new TFile(fileName);
        TTree* myTree = (TTree*) myFile->Get(directoryName);
        
        cout<<"Accessing file:"<<fileName<<endl;
        cout<<"Accessing directory:"<<directoryName<<endl;
        
        TBranch *GEN_pt_br = myTree->GetBranch("GEN_pt");
        TBranch *GEN_charge_br = myTree->GetBranch("GEN_charge");
        TBranch *BDTG_br = myTree->GetBranch("BDTG");
        TBranch *TRK_mode_RPC_br = myTree->GetBranch("TRK_mode_RPC");
        
        double a=1.0;
        double b=0.0;//b is defined but not used in the cut, b==1-a;
        double BIT=0.000001;
        Long64_t MinRATE=9999;
        double OptA=a;//best cut with min rate while high efficiency(>reference eff)
        double OptB=b;
        Int_t fill=0;//only fill 2 classes topology one time
        
        auto ROC = new TProfile("ROC","ROC Curve",100,0,1,0,1);
        auto EFFvsCUTs = new TProfile("Efficiency","Signal Efficiency vs Cuts",Bins,0,1,0,1);
        TString RATEvsCUTsTitle="";
        RATEvsCUTsTitle = RATEvsCUTsTitle + "RATE vs Cuts (Eff > "+Form("%0.2lf", EFF_REF) + ")";
        auto RATEvsCUTs = new TProfile("RATE", RATEvsCUTsTitle, Bins, 0, 1, 0, 10000);
        TH2F *Topology = new TH2F("Topology", "Class2 vs Class1", 100, 0, 1, 100, 0, 1);
        
        /* //debug plots
        TH1F *SUM = new TH1F("SUM", "SUM", 100, 0, 2);
        TH1F *CLASSONE = new TH1F("CLASSONE", "CLASSONE", 200, 0, 1);
        TH1F *CLASSTWO = new TH1F("CLASSTWO", "CLASSTWO", 200, 0, 1);
        TH1F *CUTA = new TH1F("CUTA", "CUTA", 200, 0, 1);
        TH1F *CUTB = new TH1F("CUTB", "CUTB", 200, 0, 1);
        TH1F *CHARGE = new TH1F("CHARGE", "CHARGE", 200, -100, 100);
        TH1F *PT = new TH1F("PT", "PT", 1000, 0, 1000);
        */
        
        Long64_t numEvents = myTree->GetEntries();
        cout<<">>>>>>>>>>>>>>>>>>>>>"<<endl;
        cout<<numEvents<<" events to process..."<<endl;
      
        //loop over cut on class1
        for(int i = 1; i < Bins; i++){
          
          a = (Bins-i)*1.0/Bins;//update cut on class1
          b = 1.0 - a;//store b, b is not used in cut
          
          /*
          Long64_t Z=0;
          Long64_t C1=0;
          Long64_t C2=0;
          Long64_t C3=0;
          Long64_t C4=0;
          Long64_t S=0;
          Long64_t B=0;
          */
                
          Long64_t S1=0;
          Long64_t S2=0;
          Long64_t B1=0;
          Long64_t B2=0;
          double TPR=-1.0;
          double FPR=-1.0;
          Long64_t RATE=0;
          
          for(Long64_t iEntry = 0; iEntry <numEvents; iEntry++){
              
              myTree->GetEntry(iEntry);
                    
              //access leaves under branch
              Float_t GEN_pt = (GEN_pt_br->GetLeaf("GEN_pt"))->GetValue();
              Float_t GEN_charge = (GEN_charge_br->GetLeaf("GEN_charge"))->GetValue();
              Float_t BDTG_class1 = (BDTG_br->GetLeaf("class1"))->GetValue();
              Float_t BDTG_class2 = (BDTG_br->GetLeaf("class2"))->GetValue();//not used in the cut
                
              if(fill==0){
                Topology->Fill(BDTG_class1,BDTG_class2);//sanity check off diagnoal: class2=1-class1;
              }
                  
              /* 
              //@@@Debug if accesses classes right
              SUM->Fill(BDTG_class1+BDTG_class2);
              CLASSONE->Fill(BDTG_class1);
              CLASSTWO->Fill(BDTG_class2);
              CUTA->Fill(a);
              CUTB->Fill(b);
              CHARGE->Fill(GEN_charge);
              PT->Fill(GEN_pt);
                   
              //ZB events
              if(GEN_charge < -2){Z=Z+1;}
                    
              //MC events
              if(GEN_charge > -2 && GEN_pt >= PT_CUT){S=S+1;}
              if(GEN_charge > -2 && GEN_pt < PT_CUT){B=B+1;}
              if(BDTG_class1 >= a && BDTG_class2 < b){C1=C1+1;}
              if(BDTG_class1 < a || BDTG_class2 >= b){C2=C2+1;}
              if(GEN_charge > -2 && BDTG_class1 >= a && BDTG_class2 < b){C3=C3+1;}
              if(GEN_charge > -2 && (BDTG_class1 < a || BDTG_class2 >= b)){C4=C4+1;}
              */
                    
              if(GEN_charge > -2 && GEN_pt >= PT_CUT && BDTG_class1 >= a){S2=S2+1;}
              if(GEN_charge > -2 && GEN_pt < PT_CUT && BDTG_class1 >= a){S1=S1+1;}
              if(GEN_charge > -2 && GEN_pt >= PT_CUT && BDTG_class1 < a){B2=B2+1;}
              if(GEN_charge > -2 && GEN_pt < PT_CUT && BDTG_class1 < a){B1=B1+1;}
            
            }//end loop over events
            
            fill=1;
            //Fill ROC curve
            TPR=1.0*S2/(S2+B2);
            FPR=1.0*S1/(S1+B1);
            ROC->Fill(FPR,TPR);
                  
            //Fill Signal efficiency vs cut
            EFFvsCUTs->Fill(a,TPR);
                  
            /*
            //@@@debug 
            cout<<">>>>>>>>>>>>>>>>>>>>>"<<endl;
            cout<<"a: "<<a<<" b: "<<b<<endl;
            cout<<"Z: "<<Z<<endl;
            cout<<"S: "<<S<<endl;
            cout<<"B: "<<B<<endl;
            cout<<"S+B+Z: "<<S+B+Z<<endl;
            cout<<"C1: "<<C1<<endl;
            cout<<"C2: "<<C2<<endl;
            cout<<"C1+C2: "<<C1+C2<<endl;
            cout<<"S1: "<<S1<<endl;
            cout<<"S2: "<<S2<<endl;
            cout<<"B1: "<<B1<<endl;
            cout<<"B2: "<<B2<<endl;
            cout<<"S2+B2: "<<S2+B2<<endl;
            cout<<"S1+B1: "<<S1+B1<<endl;
            cout<<"C3: "<<C3<<endl;
            cout<<"C4: "<<C4<<endl;
            cout<<"S1+S2: "<<S1+S2<<endl;
            cout<<"B1+B2: "<<B1+B2<<endl;
            cout<<"TPR: "<<TPR<<endl;
            cout<<"FPR: "<<FPR<<endl;
            //end debug
            */
                  
            //calculate ratr once signal eff higher than EFF_REF
            if(TPR >= EFF_REF){
                    
              for(Long64_t iEntry = 0; iEntry <numEvents; iEntry++){
              
                myTree->GetEntry(iEntry);
                Float_t GEN_pt = (GEN_pt_br->GetLeaf("GEN_pt"))->GetValue();//not used in rate counts
                Float_t GEN_charge = (GEN_charge_br->GetLeaf("GEN_charge"))->GetValue();
                Float_t BDTG_class1 = (BDTG_br->GetLeaf("class1"))->GetValue();
                Float_t BDTG_class2 = (BDTG_br->GetLeaf("class2"))->GetValue();//not used in cuts
                  
                //ZB events
                if(GEN_charge < -2 && BDTG_class1 >= a){RATE=RATE+1;}//after cut
              
              }//end loop over events for rate
              
              //keep note of rate 
              if(RATE < MinRATE){
                MinRATE=RATE;
                OptA=a;
                OptB=1-OptA;
              }
                    
            }//end if TPR higher than reference
             
            cout<<"a:"<<a<<" (b:"<<b<<") TPR:"<<TPR<<" FPR:"<<FPR<<" RATE:"<<RATE<<" S1:"<<S1<<" S2:"<<S2<<" B1:"<<B1<<" B2:"<<B2<<endl;
                  
            //fill rate vs cuts only for eff > ref_eff
            RATEvsCUTs->Fill(a,RATE);
          
        }//end loop over cut on class1     
        
        //==========================================
        //compare eff b/t regression and classifier
        //==========================================
        //regression 
        TString RegfileName = Cluster + "EMTFPtAssign2017/Regression.root";
        TString RegdirectoryName = "f_MODE_15_invPtTarg_invPtWgt_bitCompr_RPC/TestTree";
        TFile* myRegFile = new TFile(RegfileName);
        TTree* myRegTree = (TTree*) myRegFile->Get(RegdirectoryName);
        
        Long64_t RegnumEvents = myRegTree->GetEntries();
        
        TBranch *RegGEN_pt_br = myRegTree->GetBranch("GEN_pt");
        TBranch *RegGEN_charge_br = myRegTree->GetBranch("GEN_charge");
        TBranch *RegBDTG_br = myRegTree->GetBranch("BDTG_AWB_Sq");
        TBranch *RegTRK_mode_RPC_br = myRegTree->GetBranch("TRK_mode_RPC");

        //GEN pt distribution
        TH1F *RegCSConlyMC = new TH1F("RegCSConlyMC", "RegCSConlyMC", 50, 0, 10);
        TH1F *RegCSConlyMCCut = new TH1F("RegCSConlyMCCut", "RegCSConlyMCCut", 50, 0, 10);
        TH1F *CSConlyMC = new TH1F("CSConlyMC", "CSConlyMC", 50, 0, 10);
        TH1F *CSConlyMCCut = new TH1F("CSConlyMCCut", "CSConlyMCCut", 50, 0, 10);
        
        for(Long64_t iEntry = 0; iEntry <RegnumEvents; iEntry++){
              
                myRegTree->GetEntry(iEntry);
                Float_t GEN_pt = (RegGEN_pt_br->GetLeaf("GEN_pt"))->GetValue();
                Float_t GEN_charge = (RegGEN_charge_br->GetLeaf("GEN_charge"))->GetValue();
                Float_t BDTG = (RegBDTG_br->GetLeaf("inv_GEN_pt_trg"))->GetValue();
                Float_t TRK_mode_RPC = (RegTRK_mode_RPC_br->GetLeaf("TRK_mode_RPC"))->GetValue();
                  
                //CSC-only GEN pT distributions
                if(GEN_charge > -2 && TRK_mode_RPC == 0){
                        RegCSConlyMC->Fill(TMath::Log2(GEN_pt));
                        if(1./BDTG > 16){
                                RegCSConlyMCCut->Fill(TMath::Log2(GEN_pt));
                        }
                }
                
        }//end loop over Regression events
        
        //classifier with OptA and OptB
        for(Long64_t iEntry = 0; iEntry <numEvents; iEntry++){
              
                myTree->GetEntry(iEntry);
                Float_t GEN_pt = (GEN_pt_br->GetLeaf("GEN_pt"))->GetValue();
                Float_t GEN_charge = (GEN_charge_br->GetLeaf("GEN_charge"))->GetValue();
                Float_t BDTG_class1 = (BDTG_br->GetLeaf("class1"))->GetValue();
                Float_t BDTG_class2 = (BDTG_br->GetLeaf("class2"))->GetValue();//not used here
                Float_t TRK_mode_RPC = (TRK_mode_RPC_br->GetLeaf("TRK_mode_RPC"))->GetValue();
                  
                //CSC-only GEN pT distributions
                if(GEN_charge > -2 && TRK_mode_RPC == 0){
                        CSConlyMC->Fill(TMath::Log2(GEN_pt));
                        if(BDTG_class1 >= OptA){
                                CSConlyMCCut->Fill(TMath::Log2(GEN_pt));
                        }
                }
        }//end loop over events for rate
        
        cout<<">>>>>>>>>>>>>>>>>>>>>"<<endl;
        cout<<"OptA:"<<OptA<<" (OptB:"<<OptB<<") MinRATE:"<<MinRATE<<endl;
        
        //write to output file
        TString outFile = Cluster + "EMTFPtAssign2017/ClassifierROC_" + Form("%d", PT_CUT) + ".root";
        TFile myPlot(outFile,"RECREATE");
        
        ROC->GetXaxis()->SetTitle("FPR");
        ROC->GetYaxis()->SetTitle("TPR");
        ROC->Write();
        EFFvsCUTs->GetXaxis()->SetTitle("class1 cut");
        EFFvsCUTs->GetYaxis()->SetTitle("efficiency");
        EFFvsCUTs->Write();
        RATEvsCUTs->GetXaxis()->SetTitle("class1 cut");
        RATEvsCUTs->GetYaxis()->SetTitle("ZeroBias rate");
        RATEvsCUTs->Write();
        Topology->GetXaxis()->SetTitle("class1");
        Topology->GetYaxis()->SetTitle("class2");
        Topology->Write();
        
        /*
        SUM->Write();
        CLASSONE->Write();
        CLASSTWO->Write();
        CUTA->Write();
        CUTB->Write();
        CHARGE->Write();
        PT->Write();
        */
        
        TCanvas *C1=new TCanvas("C1","C1",700,500);
        THStack *CSConlyGENpt = new THStack("CSConlyGENpt","CSC only GEN pt: Regression vs Classifier");
        
        C1->cd();
        RegCSConlyMC->SetLineColor(1);//black
        RegCSConlyMC->SetLineStyle(1);//solid
        RegCSConlyMC->SetLineWidth(2);
        gStyle->SetOptStat(0);
        
        RegCSConlyMCCut->SetLineColor(2);//red
        RegCSConlyMCCut->SetLineStyle(1);
        RegCSConlyMCCut->SetLineWidth(2);
        gStyle->SetOptStat(0);
        
        CSConlyMC->SetLineColor(1);//black
        CSConlyMC->SetLineStyle(2);//dash
        CSConlyMC->SetLineWidth(2);
        gStyle->SetOptStat(0);
        
        CSConlyMCCut->SetLineColor(2);//red
        CSConlyMCCut->SetLineStyle(2);//dash
        CSConlyMCCut->SetLineWidth(2);
        gStyle->SetOptStat(0);
        
        CSConlyGENpt->Add(RegCSConlyMC);
        CSConlyGENpt->Add(RegCSConlyMCCut);
        CSConlyGENpt->Add(CSConlyMC);
        CSConlyGENpt->Add(CSConlyMCCut);
        CSConlyGENpt->Draw("nostack");
        CSConlyGENpt->GetXaxis()->SetTitle("log2(GEN pT)");
        C1->Modified();
        
        TLegend* L1 = new TLegend(0.1,0.7,0.7,0.9);
        TString ClassifierL1="";
        ClassifierL1 = ClassifierL1 + "Classifier: GEN pT(class1 >= " + Form("%0.4lf", OptA) + ")";
        L1->AddEntry(RegCSConlyMC, "Regression: GEN pT");
        L1->AddEntry(RegCSConlyMCCut,"Regression: GEN pT(trigger pT > 16 GeV)");
        L1->AddEntry(CSConlyMC, "Classifier: GEN pT");
        L1->AddEntry(CSConlyMCCut, ClassifierL1);
        L1->SetFillStyle(0);
        L1->SetBorderSize(0);
        L1->Draw(); 
        C1->Write();
        
        //divide histograms for eff
        TCanvas *C2=new TCanvas("C2","C2",700,500);
        THStack *CSConlyEff = new THStack("CSConlyEff","CSC only Efficiency: Regression vs Classifier");
        RegCSConlyMCCut->Divide(RegCSConlyMC);
        CSConlyMCCut->Divide(CSConlyMC);
        CSConlyEff->Add(RegCSConlyMCCut);
        CSConlyEff->Add(CSConlyMCCut);
        CSConlyEff->Draw("nostack");
        CSConlyEff->GetXaxis()->SetTitle("log2(GEN pT)");
        CSConlyEff->GetYaxis()->SetTitle("efficiency");
        C2->Modified();
        
        TLegend* L2 = new TLegend(0.1,0.7,0.7,0.9);
        TString ClassifierL2="";
        ClassifierL2 = ClassifierL2 + "Classifier: class1 >= " + Form("%0.4lf", OptA);
        L2->AddEntry(RegCSConlyMCCut,"Regression: trigger pT > 16 GeV");
        L2->AddEntry(CSConlyMCCut, ClassifierL2);
        L2->SetFillStyle(0);
        L2->SetBorderSize(0);
        L2->Draw(); 
        C2->Write();
        
        myPlot.Close();
          
}//end macro
