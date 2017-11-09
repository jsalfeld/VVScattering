#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include "TLorentzVector.h"

//#include "MitAnalysisRunII/panda/macros/80x/auxiliar_old.h"
#include "MitAnalysisRunII/panda/macros/80x/auxiliar.h"

void testTrigger(){

const int nTrgBinTestPt1  = 5; Float_t xTrgBinTestPt1[nTrgBinTestPt1+1]   = {25,30,35,45,50,10000};
const int nTrgBinTestPt2  = 7; Float_t xTrgBinTestPt2[nTrgBinTestPt2+1]   = {10,15,25,30,35,45,50,10000};
const int nTrgBinTestEta1 = 9; Float_t xTrgBinTestEta1[nTrgBinTestEta1+1] = {-2.5,-2.3,-2.2,-2.1,-1.7,1.7,2.1,2.2,2.3,2.5};
const int nTrgBinTestEta2 = 9; Float_t xTrgBinTestEta2[nTrgBinTestEta2+1] = {-2.5,-2.3,-2.2,-2.1,-1.7,1.7,2.1,2.2,2.3,2.5};

TH1D *hDTrgBinTestPt1  = new TH1D(Form("hDTrgBinTestPt1"),  Form("hDTrgBinTestPt1"),  nTrgBinTestPt1,  xTrgBinTestPt1);
TH1D *hDTrgBinTestPt2  = new TH1D(Form("hDTrgBinTestPt2"),  Form("hDTrgBinTestPt2"),  nTrgBinTestPt2,  xTrgBinTestPt2);
TH1D *hDTrgBinTestEta1 = new TH1D(Form("hDTrgBinTestEta1"), Form("hDTrgBinTestEta1"), nTrgBinTestEta1, xTrgBinTestEta1);
TH1D *hDTrgBinTestEta2 = new TH1D(Form("hDTrgBinTestEta2"), Form("hDTrgBinTestEta2"), nTrgBinTestEta2, xTrgBinTestEta2);

double trgEff[2][2][nTrgBinPt1][nTrgBinPt2][nTrgBinEta1][nTrgBinEta2];
initialize_trgEff(trgEff);

for(int pdgId1=11; pdgId1<=13; pdgId1+=2){
  for(int pdgId2=11; pdgId2<=13; pdgId2+=2){
    for(int npt1=0; npt1<nTrgBinTestPt1; npt1++){
      for(int npt2=0; npt2<nTrgBinTestPt2; npt2++){
        for(int neta1=0; neta1<nTrgBinTestEta1; neta1++){
          for(int neta2=0; neta2<nTrgBinTestEta2; neta2++){
	    if(hDTrgBinTestPt2->GetBinCenter(npt2+1) > hDTrgBinTestPt1->GetBinCenter(npt1+1)) continue;
	    int npdgId1 = TMath::Min(pdgId1-11,1);
	    int npdgId2 = TMath::Min(pdgId2-11,1);
            //printf("trgEff[%d][%d][%d][%d][%d][%d] = %f;\n",npdgId1,npdgId2,npt1,npt2,neta1,neta2,trigger_sf(hDTrgBinTestPt1->GetBinCenter(npt1+1),hDTrgBinTestEta1->GetBinCenter(neta1+1),pdgId1,hDTrgBinTestPt2->GetBinCenter(npt2+1),hDTrgBinTestEta2->GetBinCenter(neta2+1),pdgId2));
            printf("trgEff[%d][%d][%d][%d][%d][%d] = %f;\n",npdgId1,npdgId2,npt1,npt2,neta1,neta2,trigger_sf(trgEff,hDTrgBinTestPt1->GetBinCenter(npt1+1),hDTrgBinTestEta1->GetBinCenter(neta1+1),pdgId1,hDTrgBinTestPt2->GetBinCenter(npt2+1),hDTrgBinTestEta2->GetBinCenter(neta2+1),pdgId2));
          }
	}
      }
    }
  }
}
}
