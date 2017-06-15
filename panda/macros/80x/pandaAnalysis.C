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

#include "MitAnalysisRunII/panda/macros/80x/pandaFlat.C"
#include "MitAnalysisRunII/data/80x/RoccoR.cc"

enum TriggerBits {
    kMETTrig	   =(1<<0),
    kSinglePhoTrig =(1<<1),
    kMuEGTrig	   =(1<<2),
    kMuMuTrig	   =(1<<3),
    kMuTrig	   =(1<<4),
    kEGEGTrig	   =(1<<5),
    kEGTrig	   =(1<<6)
};

enum SelectionBit {
 kLoose   =(1<<0),
 kFake    =(1<<1),
 kMedium  =(1<<2),
 kTight   =(1<<3)
};

const double mass_el = 0.000510998928;
const double mass_mu = 0.10566;
void pandaAnalysis()
{
  TString dirPathRM = TString(gSystem->Getenv("CMSSW_BASE")) + "/src/MitAnalysisRunII/data/80x/rcdata.2016.v3";
  RoccoR rmcor(dirPathRM.Data());
  double lumi = 35.8;
  //*******************************************************
  //Input Files
  //*******************************************************
  TString filesPath = "/data/t3home000/ceballos/panda/v_004_0/";
  vector<TString> infileName_;  
  vector<Int_t> infileCat_;  
  infileName_.push_back(Form("%sdata.root",filesPath.Data()));                 infileCat_.push_back(0);
  infileName_.push_back(Form("%sDYJetsToLL_M-50_LO.root" ,filesPath.Data()));  infileCat_.push_back(1);
  ///*
  infileName_.push_back(Form("%sDYJetsToLL_M-10to50.root" ,filesPath.Data())); infileCat_.push_back(1);
  infileName_.push_back(Form("%sqqWW.root" ,filesPath.Data()));                infileCat_.push_back(2);
  infileName_.push_back(Form("%sggWW.root" ,filesPath.Data()));                infileCat_.push_back(2);
  infileName_.push_back(Form("%sZZ.root" ,filesPath.Data()));                  infileCat_.push_back(3);
  infileName_.push_back(Form("%sWZ.root" ,filesPath.Data()));                  infileCat_.push_back(4);
  infileName_.push_back(Form("%sWGstar.root" ,filesPath.Data()));              infileCat_.push_back(4);
  infileName_.push_back(Form("%sVVV.root" ,filesPath.Data()));                 infileCat_.push_back(5);
  infileName_.push_back(Form("%sTTV.root" ,filesPath.Data()));                 infileCat_.push_back(5);
  infileName_.push_back(Form("%sTT2L.root" ,filesPath.Data()));                infileCat_.push_back(6);
  infileName_.push_back(Form("%sTW.root" ,filesPath.Data()));                  infileCat_.push_back(6);
  infileName_.push_back(Form("%sVG.root" ,filesPath.Data()));                  infileCat_.push_back(7);
  infileName_.push_back(Form("%sH125.root" ,filesPath.Data()));                infileCat_.push_back(8);
  //*/
  const int nBinPt = 64; Float_t xbinsPt[nBinPt+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                       10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                       20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                                       30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                                       40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 
						       90, 95,100,120,140,160,180,200,250,300,
						      350,400,450,500,1000};

  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 10;
  const int histBins = 9;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {".Data", "DY", "WW", "ZZ", "WZ", "VVV", "Top", "VG", "H125"};
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  1) {nBinPlot =  60; xminPlot = 91.1876-15; xmaxPlot = 91.1876+15;}
    else if(thePlot >=  2 && thePlot <=  3) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot = 200;}
    else if(thePlot >=  4 && thePlot <=  5) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot =2000;}
    else if(thePlot >=  6 && thePlot <=  7) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot =  20;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Reset();histos->Clear();
  }

  TH2D* histoPtRecGen[2];
  histoPtRecGen[0] = new TH2D("histoPtRecGen_0", "histoPtRecGen_0", nBinPt, xbinsPt, nBinPt, xbinsPt);
  histoPtRecGen[1] = new TH2D("histoPtRecGen_1", "histoPtRecGen_1", nBinPt, xbinsPt, nBinPt, xbinsPt);
  TH1D* histoPtRecDA[2];
  histoPtRecDA[0] = new TH1D("histoPtRecDA_0", "histoPtRecDA_0", nBinPt, xbinsPt);
  histoPtRecDA[1] = new TH1D("histoPtRecDA_1", "histoPtRecDA_1", nBinPt, xbinsPt);
  TH1D* histoPtRecDY[2];
  histoPtRecDY[0] = new TH1D("histoPtRecDY_0", "histoPtRecDY_0", nBinPt, xbinsPt);
  histoPtRecDY[1] = new TH1D("histoPtRecDY_1", "histoPtRecDY_1", nBinPt, xbinsPt);

  for(UInt_t ifile=0; ifile<infileName_.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infileName_[ifile].Data());
    TFile *the_input_file = TFile::Open(infileName_[ifile].Data());
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
    
    pandaFlat thePandaFlat(the_input_tree);
    Long64_t nentries = thePandaFlat.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = thePandaFlat.LoadTree(jentry);
      if (ientry < 0) break;
      nb = thePandaFlat.fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%1000000 == 0) printf("--- reading event %8lld (%8lld) of %8lld\n",jentry,ientry,nentries);
      
      bool passTrigger = (thePandaFlat.trigger & kMuEGTrig) == kMuEGTrig || (thePandaFlat.trigger & kMuMuTrig) == kMuMuTrig ||
                         (thePandaFlat.trigger & kMuTrig)   == kMuTrig   || (thePandaFlat.trigger & kEGEGTrig) == kEGEGTrig ||
			 (thePandaFlat.trigger & kEGTrig)   == kEGTrig;
      if(passTrigger == false) continue;
      if(thePandaFlat.metFilter == 0) continue;
      
      if(thePandaFlat.nLooseLep != 2) continue;
      if(thePandaFlat.looseLep1Pt <= 25 || thePandaFlat.looseLep2Pt <= 25) continue;
      bool passLepId = ((thePandaFlat.looseLep1SelBit & kMedium) == kMedium) && ((thePandaFlat.looseLep2SelBit & kMedium) == kMedium);
      if(passLepId == false) continue;

      int lepType = -1;
      if     (abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 0;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 1;
      else if(abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 2;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 2;
      else if(abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 3;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 4;
      else if(abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 5;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 5;
      else printf("Impossible dilepton combination: %d %d\n",thePandaFlat.looseLep1PdgId,thePandaFlat.looseLep2PdgId);

      if(lepType > 1) continue;
      double thePDGMass[2] = {mass_mu, mass_mu};
      if     (abs(lepType) == 1) {thePDGMass[0] = mass_el; thePDGMass[1] = mass_el;}
      else if(abs(lepType) == 2 && abs(thePandaFlat.looseLep1PdgId)==11) {thePDGMass[0] = mass_el;}
      else if(abs(lepType) == 2 && abs(thePandaFlat.looseLep2PdgId)==11) {thePDGMass[1] = mass_el;}
      TLorentzVector v1,v2;
      v1.SetPtEtaPhiM(thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      v2.SetPtEtaPhiM(thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);

      if(TMath::Abs((v1+v2).M()-91.1876) >= 15) continue;

      double ZRecPt = (v1+v2).Pt();
      double ZGenPt = 0; bool passFid = false;
      if(thePandaFlat.looseGenLep1PdgId != 0 && thePandaFlat.looseGenLep2PdgId != 0 &&
         thePandaFlat.genLep1Pt > 25 && TMath::Abs(thePandaFlat.genLep1Eta) < 2.5 &&
	 thePandaFlat.genLep2Pt > 25 && TMath::Abs(thePandaFlat.genLep2Eta) < 2.5){
        TLorentzVector vGen1,vGen2;
        vGen1.SetPtEtaPhiM(thePandaFlat.genLep1Pt,thePandaFlat.genLep1Eta,thePandaFlat.genLep1Phi,thePDGMass[0]);
        vGen2.SetPtEtaPhiM(thePandaFlat.genLep2Pt,thePandaFlat.genLep2Eta,thePandaFlat.genLep2Phi,thePDGMass[1]);
	if(TMath::Abs((vGen1+vGen2).M()-91.1876) < 15.0) {
	  ZGenPt = (vGen1+vGen2).Pt(); passFid = true;
	}
      }

      double totalWeight = 1.0;
      int theCategory = infileCat_[ifile];
      if(theCategory != 0){
        totalWeight = thePandaFlat.normalizedWeight * lumi * thePandaFlat.sf_pu *
	              thePandaFlat.sf_trk1 * thePandaFlat.sf_medium1 *
		      thePandaFlat.sf_trk2 * thePandaFlat.sf_medium2;
      }
      histo[lepType+0][theCategory]->Fill((v1+v2).M(),totalWeight);
      histo[lepType+2][theCategory]->Fill(TMath::Min(ZRecPt, 199.999),totalWeight);
      histo[lepType+4][theCategory]->Fill(TMath::Min(ZRecPt,1999.999),totalWeight);
      histo[lepType+6][theCategory]->Fill(TMath::Abs(ZRecPt-ZGenPt),totalWeight);
      
     double dataSF = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 0, 0);
     double mcSF = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, gRandom->Rndm(), gRandom->Rndm(), 0, 0);
     printf("%f %f\n",dataSF,mcSF);

      if     (theCategory == 1){
        histoPtRecDY[lepType]->Fill(TMath::Min(ZRecPt,xbinsPt[nBinPt]-0.001),totalWeight);
      }
      else if(theCategory == 0){
        histoPtRecDA[lepType]->Fill(TMath::Min(ZRecPt,xbinsPt[nBinPt]-0.001),totalWeight);
      }
      else {
        histoPtRecDA[lepType]->Fill(TMath::Min(ZRecPt,xbinsPt[nBinPt]-0.001),-1.0*totalWeight);
      }

      if(theCategory == 1 && passFid == true){
        histoPtRecGen[lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
      }
    }

  }

  char output[200];
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    sprintf(output,"histozll_%d.root",thePlot);	
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    double totBck = 0;
    for(int i=1; i<=8; i++) totBck =totBck + histo[thePlot][i]->GetSumOfWeights();
    printf("(%d) %f (%f+%f+%f+%f+%f+%f+%f+%f)=%f\n",thePlot,histo[thePlot][0]->GetSumOfWeights(),
    histo[thePlot][1]->GetSumOfWeights(),histo[thePlot][2]->GetSumOfWeights(),
    histo[thePlot][3]->GetSumOfWeights(),histo[thePlot][4]->GetSumOfWeights(),
    histo[thePlot][5]->GetSumOfWeights(),histo[thePlot][6]->GetSumOfWeights(),
    histo[thePlot][6]->GetSumOfWeights(),histo[thePlot][8]->GetSumOfWeights(),
    totBck);
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }

  sprintf(output,"histozllPtRecGen.root"); 
  TFile* outFilePlotsA = new TFile(output,"recreate");
  outFilePlotsA->cd();
  histoPtRecGen[0]->Write();
  histoPtRecGen[1]->Write();
  outFilePlotsA->Close();

  sprintf(output,"histozllPtRec.root"); 
  TFile* outFilePlotsB = new TFile(output,"recreate");
  outFilePlotsB->cd();
  histoPtRecDA[0]->Write();
  histoPtRecDA[1]->Write();
  histoPtRecDY[0]->Write();
  histoPtRecDY[1]->Write();
  outFilePlotsB->Close();
}
