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
#include "TRandom.h"
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
void pandaAnalysis(int whichDY = 0, bool isMIT=false)
{
  TString dirPathRM = TString(gSystem->Getenv("CMSSW_BASE")) + "/src/MitAnalysisRunII/data/80x/rcdata.2016.v3";
  RoccoR rmcor(dirPathRM.Data());
  double lumi = 35.8;
  //*******************************************************
  //Input Files
  //*******************************************************
  TString filesPath = "/data/t3home000/ceballos/panda/v_004_0/";
  if(isMIT == false) filesPath = "/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/";
  vector<TString> infileName_;
  vector<Int_t> infileCat_;

  infileName_.push_back(Form("%sdata.root",filesPath.Data()));                 infileCat_.push_back(0);
  if     (whichDY == 0)
 {infileName_.push_back(Form("%sDYJetsToLL_M-50_LO.root" ,filesPath.Data()));  infileCat_.push_back(1);}
  else if(whichDY == 1)
 {infileName_.push_back(Form("%sDYJetsToLL_M-50_NLO.root",filesPath.Data()));  infileCat_.push_back(1);}
  else if(whichDY == 2)
 {infileName_.push_back(Form("%sDYJetsToLL_POWHEG.root",filesPath.Data()));    infileCat_.push_back(1);}
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

  const int nBinPt = 57; Float_t xbinsPt[nBinPt+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                       10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                       20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                                       30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                                       40, 50, 60, 70, 80, 90,100,125,150,175, 
						      200,250,300,350,400,450,500,1000};
  const int nBinPt2 = 50; Float_t xbinsPt2[nBinPt2+1]; for(int i=0; i<=nBinPt2;i++) xbinsPt2[i] = i * 1.0;

  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 20;
  const int histBins = 9;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {".Data", "DY", "WW", "ZZ", "WZ", "VVV", "Top", "VG", "H125"};
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  1) {nBinPlot = 120; xminPlot = 91.1876-15; xmaxPlot = 91.1876+15;}
    else if(thePlot >=  2 && thePlot <=  3) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot = 200;}
    else if(thePlot >=  4 && thePlot <=  5) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot =2000;}
    else if(thePlot >=  6 && thePlot <=  7) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot =  20;}
    else if(thePlot >=  8 && thePlot <= 13) {nBinPlot = 120; xminPlot = 91.1876-15; xmaxPlot = 91.1876+15;}
    else if(thePlot >= 14 && thePlot <= 17) {nBinPlot = 100; xminPlot =-50.0; xmaxPlot = 50;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Reset();histos->Clear();
  }

  TH2D* histoPtRecGen[2];
  TH1D* histoPtRecDA[2];
  TH1D* histoPtRecDY[2];
  TH1D* histoPtRecDA_MonRes[2];
  TH1D* histoPtRecDY_MonRes[2];
  TH1D* histoPtRecDA_PDF[2];
  TH1D* histoPtRecDY_PDF[2];
  TH1D* histoPtRecDA_QCD[2];
  TH1D* histoPtRecDY_QCD[2];
  TH1D* histoPtRecDA_QCDPart[2][6];
  TH1D* histoPtRecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRecGen[i] = new TH2D(Form("histoPtRecGen_%d",i), Form("histoPtRecGen_%d",i), nBinPt, xbinsPt, nBinPt, xbinsPt);
   histoPtRecDA[i]  = new TH1D(Form("histoPtRecDA_%d",i),  Form("histoPtRecDA_%d",i),  nBinPt, xbinsPt);
   histoPtRecDY[i]  = new TH1D(Form("histoPtRecDY_%d",i),  Form("histoPtRecDY_%d",i),  nBinPt, xbinsPt);

   histoPtRecDA_MonRes[i] = new TH1D(Form("histoPtRecDA_MonRes_%d",i), Form("histoPtRecDA_MonRes_%d",i), nBinPt, xbinsPt);
   histoPtRecDY_MonRes[i] = new TH1D(Form("histoPtRecDY_MonRes_%d",i), Form("histoPtRecDY_MonRes_%d",i), nBinPt, xbinsPt);

   histoPtRecDA_PDF[i] = new TH1D(Form("histoPtRecDA_PDF_%d",i), Form("histoPtRecDA_PDF_%d",i), nBinPt, xbinsPt);
   histoPtRecDY_PDF[i] = new TH1D(Form("histoPtRecDY_PDF_%d",i), Form("histoPtRecDY_PDF_%d",i), nBinPt, xbinsPt);

   histoPtRecDA_QCD[i] = new TH1D(Form("histoPtRecDA_QCD_%d",i), Form("histoPtRecDA_QCD_%d",i), nBinPt, xbinsPt);
   histoPtRecDY_QCD[i] = new TH1D(Form("histoPtRecDY_QCD_%d",i), Form("histoPtRecDY_QCD_%d",i), nBinPt, xbinsPt);

   for(int j=0; j<6; j++){
      histoPtRecDA_QCDPart[i][j] = new TH1D(Form("histoPtRecDA_QCDPart_%d_%d",i,j), Form("histoPtRecDA_QCDPart_%d_%d",i,j), nBinPt, xbinsPt);
      histoPtRecDY_QCDPart[i][j] = new TH1D(Form("histoPtRecDY_QCDPart_%d_%d",i,j), Form("histoPtRecDY_QCDPart_%d_%d",i,j), nBinPt, xbinsPt);
    }
  }

  TH2D* histoPt2RecGen[2];
  TH1D* histoPt2RecDA[2];
  TH1D* histoPt2RecDY[2];
  TH1D* histoPt2RecDA_MonRes[2];
  TH1D* histoPt2RecDY_MonRes[2];
  TH1D* histoPt2RecDA_PDF[2];
  TH1D* histoPt2RecDY_PDF[2];
  TH1D* histoPt2RecDA_QCD[2];
  TH1D* histoPt2RecDY_QCD[2];
  TH1D* histoPt2RecDA_QCDPart[2][6];
  TH1D* histoPt2RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPt2RecGen[i] = new TH2D(Form("histoPt2RecGen_%d",i), Form("histoPt2RecGen_%d",i), nBinPt2, xbinsPt2, nBinPt2, xbinsPt2);
   histoPt2RecDA[i]  = new TH1D(Form("histoPt2RecDA_%d",i),  Form("histoPt2RecDA_%d",i),  nBinPt2, xbinsPt2);
   histoPt2RecDY[i]  = new TH1D(Form("histoPt2RecDY_%d",i),  Form("histoPt2RecDY_%d",i),  nBinPt2, xbinsPt2);

   histoPt2RecDA_MonRes[i] = new TH1D(Form("histoPt2RecDA_MonRes_%d",i), Form("histoPt2RecDA_MonRes_%d",i), nBinPt2, xbinsPt2);
   histoPt2RecDY_MonRes[i] = new TH1D(Form("histoPt2RecDY_MonRes_%d",i), Form("histoPt2RecDY_MonRes_%d",i), nBinPt2, xbinsPt2);

   histoPt2RecDA_PDF[i] = new TH1D(Form("histoPt2RecDA_PDF_%d",i), Form("histoPt2RecDA_PDF_%d",i), nBinPt2, xbinsPt2);
   histoPt2RecDY_PDF[i] = new TH1D(Form("histoPt2RecDY_PDF_%d",i), Form("histoPt2RecDY_PDF_%d",i), nBinPt2, xbinsPt2);

   histoPt2RecDA_QCD[i] = new TH1D(Form("histoPt2RecDA_QCD_%d",i), Form("histoPt2RecDA_QCD_%d",i), nBinPt2, xbinsPt2);
   histoPt2RecDY_QCD[i] = new TH1D(Form("histoPt2RecDY_QCD_%d",i), Form("histoPt2RecDY_QCD_%d",i), nBinPt2, xbinsPt2);

   for(int j=0; j<6; j++){
      histoPt2RecDA_QCDPart[i][j] = new TH1D(Form("histoPt2RecDA_QCDPart_%d_%d",i,j), Form("histoPt2RecDA_QCDPart_%d_%d",i,j), nBinPt2, xbinsPt2);
      histoPt2RecDY_QCDPart[i][j] = new TH1D(Form("histoPt2RecDY_QCDPart_%d_%d",i,j), Form("histoPt2RecDY_QCDPart_%d_%d",i,j), nBinPt2, xbinsPt2);
    }
  }

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
      

      int theCategory = infileCat_[ifile];
      float lepPtSF[2] = {1,1};
      float lepPtSFSyst[2] = {1,1};
      if(theCategory == 0) { // Data
        if(abs(thePandaFlat.looseLep1PdgId)==13) {
	  double rnd[2] = {gRandom->Rndm(), gRandom->Rndm()};
          lepPtSF[0] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 0, 0);
          double muonPtSystA[100],muonPtSystB[1],muonPtSystC[5],muonPtSystD[5];
          double muonPtSystAll[4] = {lepPtSF[0], lepPtSF[0], lepPtSF[0], lepPtSF[0]};

          for(int i=0; i<100; i++)
            muonPtSystA[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 1, i);
          for(int i=0; i<1; i++)
            muonPtSystB[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 2, i);
          for(int i=0; i<5; i++)
            muonPtSystC[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 4, i);
          for(int i=0; i<5; i++)
            muonPtSystD[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 5, i);

          for(int i=0; i<100; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystA[i]) > TMath::Abs(lepPtSF[0]-muonPtSystAll[0])) muonPtSystAll[0] = TMath::Abs(lepPtSF[0]-muonPtSystA[i]) + lepPtSF[0];
          for(int i=0; i<1  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystB[i]) > TMath::Abs(lepPtSF[0]-muonPtSystAll[1])) muonPtSystAll[1] = TMath::Abs(lepPtSF[0]-muonPtSystB[i]) + lepPtSF[0];
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystC[i]) > TMath::Abs(lepPtSF[0]-muonPtSystAll[2])) muonPtSystAll[2] = TMath::Abs(lepPtSF[0]-muonPtSystC[i]) + lepPtSF[0];
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystD[i]) > TMath::Abs(lepPtSF[0]-muonPtSystAll[3])) muonPtSystAll[3] = TMath::Abs(lepPtSF[0]-muonPtSystD[i]) + lepPtSF[0];
          lepPtSFSyst[0] = sqrt(muonPtSystAll[0]*muonPtSystAll[0]/lepPtSF[0]/lepPtSF[0]+
                                muonPtSystAll[1]*muonPtSystAll[1]/lepPtSF[0]/lepPtSF[0]+
                                muonPtSystAll[2]*muonPtSystAll[2]/lepPtSF[0]/lepPtSF[0]+
                                muonPtSystAll[3]*muonPtSystAll[3]/lepPtSF[0]/lepPtSF[0]) * lepPtSF[0]/2.;

        }
        if(abs(thePandaFlat.looseLep2PdgId)==13) {
	  double rnd[2] = {gRandom->Rndm(), gRandom->Rndm()};
	  lepPtSF[1] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 0, 0);
          double muonPtSystA[100],muonPtSystB[1],muonPtSystC[5],muonPtSystD[5];
          double muonPtSystAll[4] = {lepPtSF[0], lepPtSF[0], lepPtSF[0], lepPtSF[0]};

          for(int i=0; i<100; i++)
            muonPtSystA[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 1, i);
          for(int i=0; i<1; i++)
            muonPtSystB[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 2, i);
          for(int i=0; i<5; i++)
            muonPtSystC[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 4, i);
          for(int i=0; i<5; i++)
            muonPtSystD[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 5, i);

          for(int i=0; i<100; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystA[i]) > TMath::Abs(lepPtSF[1]-muonPtSystAll[0])) muonPtSystAll[0] = TMath::Abs(lepPtSF[1]-muonPtSystA[i]) + lepPtSF[1];
          for(int i=0; i<1  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystB[i]) > TMath::Abs(lepPtSF[1]-muonPtSystAll[1])) muonPtSystAll[1] = TMath::Abs(lepPtSF[1]-muonPtSystB[i]) + lepPtSF[1];
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystC[i]) > TMath::Abs(lepPtSF[1]-muonPtSystAll[2])) muonPtSystAll[2] = TMath::Abs(lepPtSF[1]-muonPtSystC[i]) + lepPtSF[1];
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystD[i]) > TMath::Abs(lepPtSF[1]-muonPtSystAll[3])) muonPtSystAll[3] = TMath::Abs(lepPtSF[1]-muonPtSystD[i]) + lepPtSF[1];
          lepPtSFSyst[1] = sqrt(muonPtSystAll[0]*muonPtSystAll[0]/lepPtSF[1]/lepPtSF[1]+
                                muonPtSystAll[1]*muonPtSystAll[1]/lepPtSF[1]/lepPtSF[1]+
                                muonPtSystAll[2]*muonPtSystAll[2]/lepPtSF[1]/lepPtSF[1]+
                                muonPtSystAll[3]*muonPtSystAll[3]/lepPtSF[1]/lepPtSF[1]) * lepPtSF[1]/2.;

        }
      } else { // MC
        if(abs(thePandaFlat.looseLep1PdgId)==13) {
	  double rnd[2] = {gRandom->Rndm(), gRandom->Rndm()};
          lepPtSF[0] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, rnd[0], rnd[1], 0, 0);
          double muonPtSyst[11];
          double muonPtSystA[100],muonPtSystB[1],muonPtSystC[5],muonPtSystD[5];
          double muonPtSystAll[4] = {lepPtSF[0], lepPtSF[0], lepPtSF[0], lepPtSF[0]};

          for(int i=0; i<100; i++)
            muonPtSystA[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, rnd[0], rnd[1], 1, i);

          for(int i=0; i<1; i++)
            muonPtSystB[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, rnd[0], rnd[1], 2, i);
          for(int i=0; i<5; i++)
            muonPtSystC[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, rnd[0], rnd[1], 4, i);
          for(int i=0; i<5; i++)
            muonPtSystD[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, rnd[0], rnd[1], 5, i);

          for(int i=0; i<100; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystA[i]) > TMath::Abs(lepPtSF[0]-muonPtSystAll[0])) muonPtSystAll[0] = TMath::Abs(lepPtSF[0]-muonPtSystA[i]) + lepPtSF[0];
          for(int i=0; i<1  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystB[i]) > TMath::Abs(lepPtSF[0]-muonPtSystAll[1])) muonPtSystAll[1] = TMath::Abs(lepPtSF[0]-muonPtSystB[i]) + lepPtSF[0];
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystC[i]) > TMath::Abs(lepPtSF[0]-muonPtSystAll[2])) muonPtSystAll[2] = TMath::Abs(lepPtSF[0]-muonPtSystC[i]) + lepPtSF[0];
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystD[i]) > TMath::Abs(lepPtSF[0]-muonPtSystAll[3])) muonPtSystAll[3] = TMath::Abs(lepPtSF[0]-muonPtSystD[i]) + lepPtSF[0];
          lepPtSFSyst[0] = sqrt(muonPtSystAll[0]*muonPtSystAll[0]/lepPtSF[0]/lepPtSF[0]+
                                muonPtSystAll[1]*muonPtSystAll[1]/lepPtSF[0]/lepPtSF[0]+
                                muonPtSystAll[2]*muonPtSystAll[2]/lepPtSF[0]/lepPtSF[0]+
                                muonPtSystAll[3]*muonPtSystAll[3]/lepPtSF[0]/lepPtSF[0]) * lepPtSF[0]/2.;

        }
        if(abs(thePandaFlat.looseLep2PdgId)==13) {
	  double rnd[2] = {gRandom->Rndm(), gRandom->Rndm()};
	  lepPtSF[1] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 10, rnd[0], rnd[1], 0, 0);
          double muonPtSyst[11];
          double muonPtSystA[100],muonPtSystB[1],muonPtSystC[5],muonPtSystD[5];
          double muonPtSystAll[4] = {lepPtSF[1], lepPtSF[1], lepPtSF[1], lepPtSF[1]};

          for(int i=0; i<100; i++)
            muonPtSystA[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 10, rnd[0], rnd[1], 1, i);

          for(int i=0; i<1; i++)
            muonPtSystB[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 10, rnd[0], rnd[1], 2, i);
          for(int i=0; i<5; i++)
            muonPtSystC[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 10, rnd[0], rnd[1], 4, i);
          for(int i=0; i<5; i++)
            muonPtSystD[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 10, rnd[0], rnd[1], 5, i);

          for(int i=0; i<100; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystA[i]) > TMath::Abs(lepPtSF[1]-muonPtSystAll[0])) muonPtSystAll[0] = TMath::Abs(lepPtSF[1]-muonPtSystA[i]) + lepPtSF[1];
          for(int i=0; i<1  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystB[i]) > TMath::Abs(lepPtSF[1]-muonPtSystAll[1])) muonPtSystAll[1] = TMath::Abs(lepPtSF[1]-muonPtSystB[i]) + lepPtSF[1];
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystC[i]) > TMath::Abs(lepPtSF[1]-muonPtSystAll[2])) muonPtSystAll[2] = TMath::Abs(lepPtSF[1]-muonPtSystC[i]) + lepPtSF[1];
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystD[i]) > TMath::Abs(lepPtSF[1]-muonPtSystAll[3])) muonPtSystAll[3] = TMath::Abs(lepPtSF[1]-muonPtSystD[i]) + lepPtSF[1];
          lepPtSFSyst[1] = sqrt(muonPtSystAll[0]*muonPtSystAll[0]/lepPtSF[1]/lepPtSF[1]+
                                muonPtSystAll[1]*muonPtSystAll[1]/lepPtSF[1]/lepPtSF[1]+
                                muonPtSystAll[2]*muonPtSystAll[2]/lepPtSF[1]/lepPtSF[1]+
                                muonPtSystAll[3]*muonPtSystAll[3]/lepPtSF[1]/lepPtSF[1]) * lepPtSF[1]/2.;

        }
        if(abs(thePandaFlat.looseLep1PdgId)==11) {
          if(TMath::Abs(thePandaFlat.looseLep1Eta) <  1.5) lepPtSF[0] = gRandom->Gaus(1.000,0.013);
          else                                             lepPtSF[0] = gRandom->Gaus(0.993,0.025);
        }
        if(abs(thePandaFlat.looseLep2PdgId)==11) {
          if(TMath::Abs(thePandaFlat.looseLep2Eta) <  1.5) lepPtSF[1] = gRandom->Gaus(1.000,0.013);
          else                                             lepPtSF[1] = gRandom->Gaus(0.993,0.025);
        }
      }

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
      v1.SetPtEtaPhiM(thePandaFlat.looseLep1Pt*lepPtSF[0],thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      v2.SetPtEtaPhiM(thePandaFlat.looseLep2Pt*lepPtSF[1],thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      TLorentzVector vMonRes1,vMonRes2;
      vMonRes1.SetPtEtaPhiM(thePandaFlat.looseLep1Pt*sqrt(lepPtSFSyst[0]*lepPtSFSyst[0]*lepPtSFSyst[0]),thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      vMonRes2.SetPtEtaPhiM(thePandaFlat.looseLep2Pt*sqrt(lepPtSFSyst[1]*lepPtSFSyst[1]*lepPtSFSyst[1]),thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      bool passSel = TMath::Abs((v1+v2).M()-91.1876) < 15 && v1.Pt() > 25 && v2.Pt() > 25;
      bool passSystSel[1] = {TMath::Abs((vMonRes1+vMonRes2).M()-91.1876) < 15 && vMonRes1.Pt() > 25 && vMonRes2.Pt() > 25};

      double ZRecPt = (v1+v2).Pt();

      double ZRecSystPt[1] = {(vMonRes1+vMonRes2).Pt()};
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
      if(theCategory != 0){
        totalWeight = thePandaFlat.normalizedWeight * lumi * thePandaFlat.sf_pu *
	              thePandaFlat.sf_trk1 * thePandaFlat.sf_medium1 *
		      thePandaFlat.sf_trk2 * thePandaFlat.sf_medium2;
      }
      
      if(passSel){
	histo[lepType+0][theCategory]->Fill((v1+v2).M(),totalWeight);
	histo[lepType+2][theCategory]->Fill(TMath::Min(ZRecPt, 199.999),totalWeight);
	histo[lepType+4][theCategory]->Fill(TMath::Min(ZRecPt,1999.999),totalWeight);
	histo[lepType+6][theCategory]->Fill(TMath::Abs(ZRecPt-ZGenPt),totalWeight);
        if     (TMath::Abs(thePandaFlat.looseLep1Eta) <  1.5 && TMath::Abs(thePandaFlat.looseLep2Eta) <  1.5){
	  histo[lepType+ 8][theCategory]->Fill((v1+v2).M(),totalWeight);
        }
        else if(TMath::Abs(thePandaFlat.looseLep1Eta) >= 1.5 && TMath::Abs(thePandaFlat.looseLep2Eta) >= 1.5){
	  histo[lepType+12][theCategory]->Fill((v1+v2).M(),totalWeight);
        }
        else {
	  histo[lepType+10][theCategory]->Fill((v1+v2).M(),totalWeight);
        }
	histo[lepType+14][theCategory]->Fill((v1+v2).Px(),totalWeight);
	histo[lepType+16][theCategory]->Fill((v1+v2).Py(),totalWeight);

        double maxQCDscale = 1.0;
        if(theCategory != 0 && thePandaFlat.scale[0] != -1){
          maxQCDscale = (TMath::Abs(1+thePandaFlat.scale[0])+TMath::Abs(1+thePandaFlat.scale[1])+TMath::Abs(1+thePandaFlat.scale[2])+
        		 TMath::Abs(1+thePandaFlat.scale[3])+TMath::Abs(1+thePandaFlat.scale[4])+TMath::Abs(1+thePandaFlat.scale[4]))/6.0;
        }

	if     (theCategory == 1){
          // Pt
          histoPtRecDY        [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRecDY_PDF    [lepType]   ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPtRecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          // Pt2
          histoPt2RecDY        [lepType]   ->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDY_PDF    [lepType]   ->Fill(ZRecPt*ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPt2RecDY_QCDPart[lepType][0]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPt2RecDY_QCDPart[lepType][1]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPt2RecDY_QCDPart[lepType][2]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPt2RecDY_QCDPart[lepType][3]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPt2RecDY_QCDPart[lepType][4]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPt2RecDY_QCDPart[lepType][5]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
	}
	else if(theCategory == 0){
          // Pt
          histoPtRecDA[lepType]           ->Fill(ZRecPt,totalWeight);
          histoPtRecDA_PDF    [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
          // Pt2
          histoPt2RecDA[lepType]           ->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_PDF    [lepType]   ->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][0]->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][1]->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][2]->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][3]->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][4]->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][5]->Fill(ZRecPt*ZRecPt,totalWeight);
	}
	else {
          // Pt
          histoPtRecDA        [lepType]     ->Fill(ZRecPt,-1.0*totalWeight);
          histoPtRecDA_PDF    [lepType]     ->Fill(ZRecPt,-1.0*totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoPtRecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoPtRecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight);
          }
          // Pt2
          histoPt2RecDA        [lepType]     ->Fill(ZRecPt*ZRecPt,-1.0*totalWeight);
          histoPt2RecDA_PDF    [lepType]     ->Fill(ZRecPt*ZRecPt,-1.0*totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoPt2RecDA_QCDPart[lepType][0]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPt2RecDA_QCDPart[lepType][1]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPt2RecDA_QCDPart[lepType][2]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPt2RecDA_QCDPart[lepType][3]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPt2RecDA_QCDPart[lepType][4]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPt2RecDA_QCDPart[lepType][5]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoPt2RecDA_QCDPart[lepType][0]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight);
            histoPt2RecDA_QCDPart[lepType][1]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight);
            histoPt2RecDA_QCDPart[lepType][2]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight);
            histoPt2RecDA_QCDPart[lepType][3]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight);
            histoPt2RecDA_QCDPart[lepType][4]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight);
            histoPt2RecDA_QCDPart[lepType][5]->Fill(ZRecPt*ZRecPt,-1.0*totalWeight);
          }
	}

	if(theCategory == 1 && passFid == true){
          histoPtRecGen [lepType]->Fill(ZRecPt       ,ZGenPt       ,totalWeight);
          histoPt2RecGen[lepType]->Fill(ZRecPt*ZRecPt,ZGenPt*ZGenPt,totalWeight);
	}
      }
      
      if(passSystSel[0]){
	if     (theCategory == 1){
          histoPtRecDY_MonRes [lepType]->Fill(ZRecSystPt[0],totalWeight);
          histoPt2RecDY_MonRes[lepType]->Fill(ZRecSystPt[0]*ZRecSystPt[0],totalWeight);
	}
	else if(theCategory == 0){
          histoPtRecDA_MonRes[lepType] ->Fill(ZRecSystPt[0],totalWeight);
          histoPt2RecDA_MonRes[lepType]->Fill(ZRecSystPt[0]*ZRecSystPt[0],totalWeight);
	}
	else {
          histoPtRecDA_MonRes[lepType] ->Fill(ZRecSystPt[0],-1.0*totalWeight);
          histoPt2RecDA_MonRes[lepType]->Fill(ZRecSystPt[0]*ZRecSystPt[0],-1.0*totalWeight);
	}
      }
    } // end event loop
  } // end samples loop

  for(int ntype=0; ntype<2; ntype++){
    printf("QCD(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRecDA_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRecDA_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRecDA_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRecDA_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRecDA_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRecDA_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRecDA[ntype]->GetSumOfWeights(),
           histoPtRecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinPt; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRecDA_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRecDA[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecDA[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecDA[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRecDA[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRecDA[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[0]/histoPtRecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRecDA_QCD[ntype]->SetBinContent(nb, histoPtRecDA[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRecDY_QCD[ntype]->SetBinContent(nb, histoPtRecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }
  for(int ntype=0; ntype<2; ntype++){
    printf("QCD2(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPt2RecDA_QCDPart[ntype][0]->GetSumOfWeights(),histoPt2RecDA_QCDPart[ntype][1]->GetSumOfWeights(),histoPt2RecDA_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPt2RecDA_QCDPart[ntype][3]->GetSumOfWeights(),histoPt2RecDA_QCDPart[ntype][4]->GetSumOfWeights(),histoPt2RecDA_QCDPart[ntype][5]->GetSumOfWeights(),histoPt2RecDA[ntype]->GetSumOfWeights(),
           histoPt2RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPt2RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPt2RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPt2RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPt2RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPt2RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPt2RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinPt2; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPt2RecDA_QCDPart[ntype][0]->GetBinContent(nb)-histoPt2RecDA[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPt2RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPt2RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPt2RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPt2RecDA[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPt2RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPt2RecDA[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPt2RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPt2RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPt2RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPt2RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPt2RecDA[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPt2RecDA[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPt2RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[0]/histoPt2RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPt2RecDA_QCD[ntype]->SetBinContent(nb, histoPt2RecDA[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPt2RecDY_QCD[ntype]->SetBinContent(nb, histoPt2RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
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

  histoPtRecDA[0]->Write();
  histoPtRecDA[1]->Write();
  histoPtRecDY[0]->Write();
  histoPtRecDY[1]->Write();
  histoPtRecDA_MonRes[0]->Write();
  histoPtRecDA_MonRes[1]->Write();
  histoPtRecDY_MonRes[0]->Write();
  histoPtRecDY_MonRes[1]->Write();
  histoPtRecDA_PDF[0]->Write();
  histoPtRecDA_PDF[1]->Write();
  histoPtRecDY_PDF[0]->Write();
  histoPtRecDY_PDF[1]->Write();
  histoPtRecDA_QCD[0]->Write();
  histoPtRecDA_QCD[1]->Write();
  histoPtRecDY_QCD[0]->Write();
  histoPtRecDY_QCD[1]->Write();

  outFilePlotsA->Close();

  sprintf(output,"histozllPt2RecGen.root"); 
  TFile* outFilePlotsB = new TFile(output,"recreate");
  outFilePlotsB->cd();

  histoPt2RecGen[0]->Write();
  histoPt2RecGen[1]->Write();

  histoPt2RecDA[0]->Write();
  histoPt2RecDA[1]->Write();
  histoPt2RecDY[0]->Write();
  histoPt2RecDY[1]->Write();
  histoPt2RecDA_MonRes[0]->Write();
  histoPt2RecDA_MonRes[1]->Write();
  histoPt2RecDY_MonRes[0]->Write();
  histoPt2RecDY_MonRes[1]->Write();
  histoPt2RecDA_PDF[0]->Write();
  histoPt2RecDA_PDF[1]->Write();
  histoPt2RecDY_PDF[0]->Write();
  histoPt2RecDY_PDF[1]->Write();
  histoPt2RecDA_QCD[0]->Write();
  histoPt2RecDA_QCD[1]->Write();
  histoPt2RecDY_QCD[0]->Write();
  histoPt2RecDY_QCD[1]->Write();

  outFilePlotsB->Close();
}
