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
        int PT_CUT = 32;//the classifier trained on this cut
        TString pt_cut = "32";//string format of PT_CUT
        Float_t EFF_REF = 0.95;//compare rate with BDT Regression
        TString eff_ref = "0.95";//string format of EFF_REF
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
        Float_t a=0.0;
        Float_t b=0.0;
      
        myTree->SetBranchAddress("GEN_pt",&GEN_pt);
        myTree->SetBranchAddress("GEN_charge",&GEN_charge);
        myTree->SetBranchAddress("BDTG",&BDTG_class1);
        myTree->SetBranchAddress("BDTG",&BDTG_class2);
        cout<<"Accessing directory:"<<directoryName<<endl;
        
        auto ROC = new TProfile("ROC","ROC Curve",100,0,1,0,1);
        auto EFFvsCUTs = new TProfile2D("Efficiency","Signal Efficiency vs Cuts",100,0,1,100,0,1,0,1);
        auto RATEvsCUTs = new TProfile2D("RATE","RATE vs Cuts (Eff > " + eff_ref +")",100,0,1,100,0,1,0,1000);
  
        Long64_t numEvents = myTree->GetEntries();
        cout<<">>>>>>>>>>>>>>>>>>>>>"<<endl;
        cout<<numEvents<<" events to process..."<<endl;
      
        //loop over cut on class1
        for(int i = 1; i < 10; i++){
          
          
          //loop over cut on class2
          for(int j = 1; j < 10; j++){
            
            Float_t S1=0;
            Float_t S2=0;
            Float_t B1=0;
            Float_t B2=0;
            Float_t TPR=-1.0;
            Float_t FPR=-1.0;
            Float_t RATE=0;
            a = 0.1*i;//update cut on class1
            b = 0.1*j;//update cut on class2
            
            for(Long64_t iEntry = 0; iEntry <numEvents; iEntry++){
              
              myTree->GetEntry(iEntry);
                     
              //MC events
              if(GEN_charge > -2){
                
                //predict signal
                if(BDTG_class1 >= a && BDTG_class2 < b){
                  if(GEN_pt >= PT_CUT){
                    S2++;
                  }
                  else{
                    S1++;
                  }
                }
                
                //predict bkg
                else{
                
                  if(GEN_pt >= PT_CUT){
                    B2++;
                  }
                  else{
                    B1++;
                  }
                
                }//end if else prediction
              
              }//end if MC
            
            }//end loop over events
              
            //Fill ROC curve
            TPR=S2/(S2+B2);
            FPR=S1/(S1+B1);
            ROC->Fill(FPR,TPR);
                
            //Fill Signal efficiency vs cut
            EFFvsCUTs->Fill(a,b,TPR);
            
            //keep note of the rate whenever efficiency is higher than EFF_REF
            if(TPR >= EFF_REF){
                    
              for(Long64_t iEntry = 0; iEntry <numEvents; iEntry++){
              
                myTree->GetEntry(iEntry);
                     
                //ZB events
                if(GEN_charge < -2){
                  //rate
                  if(BDTG_class1 >= a && BDTG_class2 < b){
                    RATE++;
                  }//after cut
                
                }//end ZB
              
              }//end loop over events for rate
                    
            }//end if TPR higher than reference
                  
            //fill rate vs cuts
            RATEvsCUTs->Fill(a,b,RATE);

          }//end loop over cut on class2
          
        }//end loop over cut on class1     
         
        //write to output file
        TFile myPlot("/home/ws13/TMVA/TMVA/EMTFPtAssign2017/ClassifierROC_" + pt_cut + ".root","RECREATE");
        ROC->Write();
        EFFvsCUTs->Write();
        RATEvsCUTs->Write();
        
        myPlot.Close();
          
}//end macro
