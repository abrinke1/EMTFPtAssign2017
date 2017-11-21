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

//=***************************************************************************************
//=Study the pT classifier performance and compare with BDT regression normally used at P5
//=including ROC cureve, signal-to-background ratio, etc
//=
//=Wei Shi @ Nov 11, 2017 CERN Geneva
//=***************************************************************************************
//=Reference Table
//========================================================================================
//=         Test MC Events   #                     |  GEN < PT_CUT    |    GEN >= PT_CUT
//========================================================================================
//=MC True: Signal (class1)                        |       NO         |         S
//=---------------------------------------------------------------------------------------
//=MC True: Background (class2)                    |       B          |         NO
//========================================================================================
//=Cut: class1 >= a && class2 < b (Signal)         |       S1         |         S2
//=---------------------------------------------------------------------------------------
//=Cut: complementary cut to signal (Background)   |       B1         |         B2
//========================================================================================
//=True Postive Rate: S2/(S2+B2), i.e. S2/S, signal efficiency/plateau trigger efficiency
//=False Positive Rate: S1/(S1+B1)
//S=S2+B2
//B=S1+B1
void ClassifierROC()
{
        //USER modify here ONLY//
        //====================================================
        Int_t PT_CUT = 32;//the classifier trained on this cut
        TString pt_cut = "32";//string format of PT_CUT
        Float_t EFF_REF = 0.95;//compare rate with BDT Regression
        TString eff_ref = "0.95";//string format of EFF_REF
        Int_t Bins=10;//bins on class cut
        //===================================================
        
        TString fileName="/home/ws13/TMVA/TMVA/EMTFPtAssign2017/pTMulticlass_MODE_15_bitCompr_RPC_"+ pt_cut +".root";
        TString directoryName="f_MODE_15_noWgt_bitCompr_RPC/TestTree";
        TFile* myFile = new TFile(fileName);
        TTree* myTree = (TTree*) myFile->Get(directoryName);
        
        cout<<"Accessing file:"<<fileName<<endl;
        
        Float_t GEN_pt;
        Float_t GEN_charge;
        Float_t BDTG_class1;
        Float_t BDTG_class2;
        double a=0.0;
        double b=1.0;//use b <= 1-a
        double BIT=0.00001;
      
        myTree->SetBranchAddress("GEN_pt",&GEN_pt);
        myTree->SetBranchAddress("GEN_charge",&GEN_charge);
        myTree->SetBranchAddress("BDTG",&BDTG_class2);//first is allocating class2 in the branch
        myTree->SetBranchAddress("BDTG",&BDTG_class1);
        cout<<"Accessing directory:"<<directoryName<<endl;
        
        auto ROC = new TProfile("ROC","ROC Curve",100,0,1,0,1);
        auto EFFvsCUTs = new TProfile2D("Efficiency","Signal Efficiency vs Cuts",Bins,0,1,Bins,0,1,0,1);
        auto RATEvsCUTs = new TProfile2D("RATE","RATE vs Cuts (Eff > " + eff_ref +")",Bins,0,1,Bins,0,1,0,10000);
        /*
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
          
          a = i*1.0/Bins;//update cut on class1
          b = 1.0 - a;//update cut on class2
                
          //loop over cut on class2
          while(b > BIT){
            
            Long64_t Z=0;
            Long64_t C1=0;
            Long64_t C2=0;
            Long64_t C3=0;
            Long64_t C4=0;
            Long64_t S=0;
            Long64_t B=0;
            Long64_t S1=0;
            Long64_t S2=0;
            Long64_t B1=0;
            Long64_t B2=0;
            double TPR=-1.0;
            double FPR=-1.0;
            Long64_t RATE=0;
            
            for(Long64_t iEntry = 0; iEntry <numEvents; iEntry++){
              
              myTree->GetEntry(iEntry);
              /*      
              //@@@Debug if accesses classes right
              SUM->Fill(BDTG_class1+BDTG_class2);
              CLASSONE->Fill(BDTG_class1);
              CLASSTWO->Fill(BDTG_class2);
              CUTA->Fill(a);
              CUTB->Fill(b);
              CHARGE->Fill(GEN_charge);
              PT->Fill(GEN_pt);
              */
                    
              //ZB events
              //if(GEN_charge < -2){Z=Z+1;}
                    
              //MC events
              //if(GEN_charge > -2 && GEN_pt >= PT_CUT){S=S+1;}
              //if(GEN_charge > -2 && GEN_pt < PT_CUT){B=B+1;}
              //if(BDTG_class1 >= a && BDTG_class2 < b){C1=C1+1;}
              //if(BDTG_class1 < a || BDTG_class2 >= b){C2=C2+1;}
              //if(GEN_charge > -2 && BDTG_class1 >= a && BDTG_class2 < b){C3=C3+1;}
              //if(GEN_charge > -2 && (BDTG_class1 < a || BDTG_class2 >= b)){C4=C4+1;}
              if(GEN_charge > -2 && GEN_pt >= PT_CUT && BDTG_class1 >= a && BDTG_class2 < b){S2=S2+1;}
              if(GEN_charge > -2 && GEN_pt < PT_CUT && BDTG_class1 >= a && BDTG_class2 < b){S1=S1+1;}
              if(GEN_charge > -2 && GEN_pt >= PT_CUT && (BDTG_class1 < a || BDTG_class2 >= b)){B2=B2+1;}
              if(GEN_charge > -2 && GEN_pt < PT_CUT && (BDTG_class1 < a || BDTG_class2 >= b)){B1=B1+1;}
            
            }//end loop over events
                  
            //Fill ROC curve
            TPR=1.0*S2/(S2+B2);
            FPR=1.0*S1/(S1+B1);
            ROC->Fill(FPR,TPR);
                  
            //Fill Signal efficiency vs cut
            EFFvsCUTs->Fill(a,b,TPR);
                  
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
                  
            //keep note of the rate whenever efficiency is higher than EFF_REF
            if(TPR >= EFF_REF){
                    
              for(Long64_t iEntry = 0; iEntry <numEvents; iEntry++){
              
                myTree->GetEntry(iEntry);
                  
                //ZB events
                if(GEN_charge < -2 && BDTG_class1 >= a && BDTG_class2 < b){RATE=RATE+1;}//after cut
              
              }//end loop over events for rate
                    
            }//end if TPR higher than reference
             
            cout<<"RATE: "<<RATE<<endl;
                  
            //fill rate vs cuts
            RATEvsCUTs->Fill(a,b,RATE);
            b = b - 1.0/Bins;//update b
                  
          }//end while for class 2 cut
          
        }//end loop over cut on class1     
         
        //write to output file
        TFile myPlot("/home/ws13/TMVA/TMVA/EMTFPtAssign2017/ClassifierROC_" + pt_cut + ".root","RECREATE");
        ROC->Write();
        EFFvsCUTs->Write();
        RATEvsCUTs->Write();
        /*
        SUM->Write();
        CLASSONE->Write();
        CLASSTWO->Write();
        CUTA->Write();
        CUTB->Write();
        CHARGE->Write();
        PT->Write();
        */
        myPlot.Close();
          
}//end macro
