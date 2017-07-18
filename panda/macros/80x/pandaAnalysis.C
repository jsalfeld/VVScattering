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
  double k_eff = -0.5 * sqrt(20285930./12446486.);
  //*******************************************************
  //Input Files
  //*******************************************************
  TString filesPath = "/data/t3home000/ceballos/panda/v_004_0/";
  if(isMIT == false) filesPath = "/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/";
  vector<TString> infileName_;
  vector<Int_t> infileCat_;

  infileName_.push_back(Form("%sdata.root",filesPath.Data()));                 infileCat_.push_back(0);
  //infileName_.push_back(Form("%sqqWW.root" ,filesPath.Data()));                infileCat_.push_back(1);
  //infileName_.push_back(Form("%sggWW.root" ,filesPath.Data()));                infileCat_.push_back(1);
  if     (whichDY == 0)
 {infileName_.push_back(Form("%sDYJetsToLL_M-50_LO.root" ,filesPath.Data()));  infileCat_.push_back(2);}
  else if(whichDY == 1)
 {infileName_.push_back(Form("%sDYJetsToLL_M-50_NLO.root",filesPath.Data()));  infileCat_.push_back(2);}
  else if(whichDY == 2)
 {infileName_.push_back(Form("%sDYJetsToLL_POWHEG.root",filesPath.Data()));    infileCat_.push_back(2);}
  infileName_.push_back(Form("%sDYJetsToLL_M-10to50.root" ,filesPath.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sTT2L.root" ,filesPath.Data()));                infileCat_.push_back(3);
  //infileName_.push_back(Form("%sTW.root" ,filesPath.Data()));                  infileCat_.push_back(3);
  infileName_.push_back(Form("%sZZ.root" ,filesPath.Data()));                  infileCat_.push_back(4);
  infileName_.push_back(Form("%sWZ.root" ,filesPath.Data()));                  infileCat_.push_back(4);
  infileName_.push_back(Form("%sVVV.root" ,filesPath.Data()));                 infileCat_.push_back(4);
  infileName_.push_back(Form("%sTTV.root" ,filesPath.Data()));                 infileCat_.push_back(4);
  //infileName_.push_back(Form("%sWGstar.root" ,filesPath.Data()));              infileCat_.push_back(4);
  //infileName_.push_back(Form("%sVG.root" ,filesPath.Data()));                  infileCat_.push_back(6);
  //infileName_.push_back(Form("%sH125.root" ,filesPath.Data()));                infileCat_.push_back(7);

  const int nBinPt = 57; Float_t xbinsPt[nBinPt+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                       10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                       20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                                       30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                                       40, 50, 60, 70, 80, 90,100,125,150,175, 
						      200,250,300,350,400,450,500,1000};

  const int nBinRecoPt =114; Float_t xbinsRecoPt[nBinRecoPt+1] = 
                                                      {  0,0.5,   1,1.5,   2,2.5,   3,3.5,   4,4.5,   5,5.5,   6,6.5,   7,7.5,   8,8.5,   9,9.5,
                                                       10,10.5, 11,11.5, 12,12.5, 13,13.5, 14,14.5, 15,15.5, 16,16.5, 17,17.5, 18,18.5, 19,19.5,
                                                       20,20.5, 21,21.5, 22,22.5, 23,23.5, 24,24.5, 25,25.5, 26,26.5, 27,27.5, 28,28.5, 29,29.5,
                                                       30,30.5, 31,31.5, 32,32.5, 33,33.5, 34,34.5, 35,35.5, 36,36.5, 37,37.5, 38,38.5, 39,39.5,
                                                         40,45,   50,55,   60,65,   70,75,   80,85,   90,95, 100,115, 125,140, 150,160, 175,190, 
						       200,225, 250,275, 300,325, 350,375, 400,425, 450,475, 500,750,    1000};

  const int nBinPt2 = 50; Float_t xbinsPt2[nBinPt2+1]; for(int i=0; i<=nBinPt2;i++) xbinsPt2[i] = i * 1.0;

  const int nBinRap = 24; Float_t xbinsRap[nBinRap+1]; for(int i=0; i<=nBinRap;i++) xbinsRap[i] = i * 0.1;

  const int nBinPtRap0 = 36; Float_t xbinsPtRap0[nBinPtRap0+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                		   10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                		   20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
                                                		  100,150,200,300,400,500,1000};
  const int nBinRecoPtRap0 = 72; Float_t xbinsRecoPtRap0[nBinRecoPtRap0+1] = 
                                                                 {  0,0.5,   1,1.5,   2,2.5,   3,3.5,   4,4.5,   5,5.5,   6,6.5,   7,7.5,   8,8.5,   9,9.5,
                                                                  10,10.5, 11,11.5, 12,12.5, 13,13.5, 14,14.5, 15,15.5, 16,16.5, 17,17.5, 18,18.5, 19,19.5,
                                                		    20,23,  25,28,    30,33,   35,38,   40,45,   50,55,   60,65,   70,75,   80,85,   90,95,
                                                		  100,125,150,175,200,250,   300,350, 400,450, 500,750,1000};
  const int nBinPtRap1 = 36; Float_t xbinsPtRap1[nBinPtRap1+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                		   10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                		   20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
                                                		  100,150,200,300,400,500,1000};
  const int nBinRecoPtRap1 = 72; Float_t xbinsRecoPtRap1[nBinRecoPtRap1+1] = 
                                                                 {  0,0.5,   1,1.5,   2,2.5,   3,3.5,   4,4.5,   5,5.5,   6,6.5,   7,7.5,   8,8.5,   9,9.5,
                                                                  10,10.5, 11,11.5, 12,12.5, 13,13.5, 14,14.5, 15,15.5, 16,16.5, 17,17.5, 18,18.5, 19,19.5,
                                                		    20,23,  25,28,    30,33,   35,38,   40,45,   50,55,   60,65,   70,75,   80,85,   90,95,
                                                		  100,125,150,175,200,250,   300,350, 400,450, 500,750,1000};
  const int nBinPtRap2 = 36; Float_t xbinsPtRap2[nBinPtRap2+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                		   10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                		   20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
                                                		  100,150,200,300,400,500,1000};
  const int nBinRecoPtRap2 = 72; Float_t xbinsRecoPtRap2[nBinRecoPtRap2+1] = 
                                                                 {  0,0.5,   1,1.5,   2,2.5,   3,3.5,   4,4.5,   5,5.5,   6,6.5,   7,7.5,   8,8.5,   9,9.5,
                                                                  10,10.5, 11,11.5, 12,12.5, 13,13.5, 14,14.5, 15,15.5, 16,16.5, 17,17.5, 18,18.5, 19,19.5,
                                                		    20,23,  25,28,    30,33,   35,38,   40,45,   50,55,   60,65,   70,75,   80,85,   90,95,
                                                		  100,125,150,175,200,250,   300,350, 400,450, 500,750,1000};
  const int nBinPtRap3 = 36; Float_t xbinsPtRap3[nBinPtRap3+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                		   10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                		   20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
                                                		  100,150,200,300,400,500,1000};
  const int nBinRecoPtRap3 = 72; Float_t xbinsRecoPtRap3[nBinRecoPtRap3+1] = 
                                                                 {  0,0.5,   1,1.5,   2,2.5,   3,3.5,   4,4.5,   5,5.5,   6,6.5,   7,7.5,   8,8.5,   9,9.5,
                                                                  10,10.5, 11,11.5, 12,12.5, 13,13.5, 14,14.5, 15,15.5, 16,16.5, 17,17.5, 18,18.5, 19,19.5,
                                                		    20,23,  25,28,    30,33,   35,38,   40,45,   50,55,   60,65,   70,75,   80,85,   90,95,
                                                		  100,125,150,175,200,250,   300,350, 400,450, 500,750,1000};
  const int nBinPtRap4 = 36; Float_t xbinsPtRap4[nBinPtRap4+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                                  10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                                  20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
                                                                 100,150,200,300,400,500,1000};
  const int nBinRecoPtRap4 = 72; Float_t xbinsRecoPtRap4[nBinRecoPtRap4+1] = 
                                                                 {  0,0.5,   1,1.5,   2,2.5,   3,3.5,   4,4.5,   5,5.5,   6,6.5,   7,7.5,   8,8.5,   9,9.5,
                                                                  10,10.5, 11,11.5, 12,12.5, 13,13.5, 14,14.5, 15,15.5, 16,16.5, 17,17.5, 18,18.5, 19,19.5,
                                                		    20,23,  25,28,    30,33,   35,38,   40,45,   50,55,   60,65,   70,75,   80,85,   90,95,
                                                		  100,125,150,175,200,250,   300,350, 400,450, 500,750,1000};

  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 20;
  const int histBins = 9;
  TH1D* histo[allPlots][histBins];
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  1) {nBinPlot = 120; xminPlot = 91.1876-15; xmaxPlot = 91.1876+15;}
    else if(thePlot >=  2 && thePlot <=  3) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot = 200;}
    else if(thePlot >=  4 && thePlot <=  5) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot =2000;}
    else if(thePlot >=  6 && thePlot <=  7) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot =  20;}
    else if(thePlot >=  8 && thePlot <= 13) {nBinPlot = 120; xminPlot = 91.1876-15; xmaxPlot = 91.1876+15;}
    else if(thePlot >= 14 && thePlot <= 17) {nBinPlot = 100; xminPlot =-50.0; xmaxPlot = 50;}
    else if(thePlot >= 18 && thePlot <= 19) {nBinPlot = 100; xminPlot =  0.0; xmaxPlot = 2.5;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Reset();histos->Clear();
  }

  TH2D* histoPtRecGen[2];
  TH2D* histoPtRecGen_LepEff[2];
  TH1D* histoPtRecDA[2];
  TH1D* histoPtRecDY[2];
  TH1D* histoPtRecDA_MomRes[2];
  TH1D* histoPtRecDY_MomRes[2];
  TH1D* histoPtRecDA_LepEff[2];
  TH1D* histoPtRecDY_LepEff[2];
  TH1D* histoPtRecDA_PDF[2];
  TH1D* histoPtRecDY_PDF[2];
  TH1D* histoPtRecDA_QCD[2];
  TH1D* histoPtRecDY_QCD[2];
  TH1D* histoPtRecDA_QCDPart[2][6];
  TH1D* histoPtRecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRecGen[i] = new TH2D(Form("histoPtRecGen_%d",i), Form("histoPtRecGen_%d",i), nBinRecoPt, xbinsRecoPt, nBinPt, xbinsPt);
   histoPtRecGen_LepEff[i] = new TH2D(Form("histoPtRecGen_LepEff_%d",i), Form("histoPtRecGen_LepEff_%d",i), nBinRecoPt, xbinsRecoPt, nBinPt, xbinsPt);
   histoPtRecDA[i]  = new TH1D(Form("histoPtRecDA_%d",i),  Form("histoPtRecDA_%d",i),  nBinRecoPt, xbinsRecoPt);
   histoPtRecDY[i]  = new TH1D(Form("histoPtRecDY_%d",i),  Form("histoPtRecDY_%d",i),  nBinRecoPt, xbinsRecoPt);

   histoPtRecDA_MomRes[i] = new TH1D(Form("histoPtRecDA_MomRes_%d",i), Form("histoPtRecDA_MomRes_%d",i), nBinRecoPt, xbinsRecoPt);
   histoPtRecDY_MomRes[i] = new TH1D(Form("histoPtRecDY_MomRes_%d",i), Form("histoPtRecDY_MomRes_%d",i), nBinRecoPt, xbinsRecoPt);

   histoPtRecDA_LepEff[i] = new TH1D(Form("histoPtRecDA_LepEff_%d",i), Form("histoPtRecDA_LepEff_%d",i), nBinRecoPt, xbinsRecoPt);
   histoPtRecDY_LepEff[i] = new TH1D(Form("histoPtRecDY_LepEff_%d",i), Form("histoPtRecDY_LepEff_%d",i), nBinRecoPt, xbinsRecoPt);

   histoPtRecDA_PDF[i] = new TH1D(Form("histoPtRecDA_PDF_%d",i), Form("histoPtRecDA_PDF_%d",i), nBinRecoPt, xbinsRecoPt);
   histoPtRecDY_PDF[i] = new TH1D(Form("histoPtRecDY_PDF_%d",i), Form("histoPtRecDY_PDF_%d",i), nBinRecoPt, xbinsRecoPt);

   histoPtRecDA_QCD[i] = new TH1D(Form("histoPtRecDA_QCD_%d",i), Form("histoPtRecDA_QCD_%d",i), nBinRecoPt, xbinsRecoPt);
   histoPtRecDY_QCD[i] = new TH1D(Form("histoPtRecDY_QCD_%d",i), Form("histoPtRecDY_QCD_%d",i), nBinRecoPt, xbinsRecoPt);

   for(int j=0; j<6; j++){
      histoPtRecDA_QCDPart[i][j] = new TH1D(Form("histoPtRecDA_QCDPart_%d_%d",i,j), Form("histoPtRecDA_QCDPart_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
      histoPtRecDY_QCDPart[i][j] = new TH1D(Form("histoPtRecDY_QCDPart_%d_%d",i,j), Form("histoPtRecDY_QCDPart_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
    }
  }

  TH2D* histoPt2RecGen[2];
  TH2D* histoPt2RecGen_LepEff[2];
  TH1D* histoPt2RecDA[2];
  TH1D* histoPt2RecDY[2];
  TH1D* histoPt2RecDA_MomRes[2];
  TH1D* histoPt2RecDY_MomRes[2];
  TH1D* histoPt2RecDA_LepEff[2];
  TH1D* histoPt2RecDY_LepEff[2];
  TH1D* histoPt2RecDA_PDF[2];
  TH1D* histoPt2RecDY_PDF[2];
  TH1D* histoPt2RecDA_QCD[2];
  TH1D* histoPt2RecDY_QCD[2];
  TH1D* histoPt2RecDA_QCDPart[2][6];
  TH1D* histoPt2RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPt2RecGen[i] = new TH2D(Form("histoPt2RecGen_%d",i), Form("histoPt2RecGen_%d",i), nBinPt2, xbinsPt2, nBinPt2, xbinsPt2);
   histoPt2RecGen_LepEff[i] = new TH2D(Form("histoPt2RecGen_LepEff_%d",i), Form("histoPt2RecGen_LepEff_%d",i), nBinPt2, xbinsPt2, nBinPt2, xbinsPt2);
   histoPt2RecDA[i]  = new TH1D(Form("histoPt2RecDA_%d",i),  Form("histoPt2RecDA_%d",i),  nBinPt2, xbinsPt2);
   histoPt2RecDY[i]  = new TH1D(Form("histoPt2RecDY_%d",i),  Form("histoPt2RecDY_%d",i),  nBinPt2, xbinsPt2);

   histoPt2RecDA_MomRes[i] = new TH1D(Form("histoPt2RecDA_MomRes_%d",i), Form("histoPt2RecDA_MomRes_%d",i), nBinPt2, xbinsPt2);
   histoPt2RecDY_MomRes[i] = new TH1D(Form("histoPt2RecDY_MomRes_%d",i), Form("histoPt2RecDY_MomRes_%d",i), nBinPt2, xbinsPt2);

   histoPt2RecDA_LepEff[i] = new TH1D(Form("histoPt2RecDA_LepEff_%d",i), Form("histoPt2RecDA_LepEff_%d",i), nBinPt2, xbinsPt2);
   histoPt2RecDY_LepEff[i] = new TH1D(Form("histoPt2RecDY_LepEff_%d",i), Form("histoPt2RecDY_LepEff_%d",i), nBinPt2, xbinsPt2);

   histoPt2RecDA_PDF[i] = new TH1D(Form("histoPt2RecDA_PDF_%d",i), Form("histoPt2RecDA_PDF_%d",i), nBinPt2, xbinsPt2);
   histoPt2RecDY_PDF[i] = new TH1D(Form("histoPt2RecDY_PDF_%d",i), Form("histoPt2RecDY_PDF_%d",i), nBinPt2, xbinsPt2);

   histoPt2RecDA_QCD[i] = new TH1D(Form("histoPt2RecDA_QCD_%d",i), Form("histoPt2RecDA_QCD_%d",i), nBinPt2, xbinsPt2);
   histoPt2RecDY_QCD[i] = new TH1D(Form("histoPt2RecDY_QCD_%d",i), Form("histoPt2RecDY_QCD_%d",i), nBinPt2, xbinsPt2);

   for(int j=0; j<6; j++){
      histoPt2RecDA_QCDPart[i][j] = new TH1D(Form("histoPt2RecDA_QCDPart_%d_%d",i,j), Form("histoPt2RecDA_QCDPart_%d_%d",i,j), nBinPt2, xbinsPt2);
      histoPt2RecDY_QCDPart[i][j] = new TH1D(Form("histoPt2RecDY_QCDPart_%d_%d",i,j), Form("histoPt2RecDY_QCDPart_%d_%d",i,j), nBinPt2, xbinsPt2);
    }
  }

  TH2D* histoRapRecGen[2];
  TH2D* histoRapRecGen_LepEff[2];
  TH1D* histoRapRecDA[2];
  TH1D* histoRapRecDY[2];
  TH1D* histoRapRecDA_MomRes[2];
  TH1D* histoRapRecDY_MomRes[2];
  TH1D* histoRapRecDA_LepEff[2];
  TH1D* histoRapRecDY_LepEff[2];
  TH1D* histoRapRecDA_PDF[2];
  TH1D* histoRapRecDY_PDF[2];
  TH1D* histoRapRecDA_QCD[2];
  TH1D* histoRapRecDY_QCD[2];
  TH1D* histoRapRecDA_QCDPart[2][6];
  TH1D* histoRapRecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoRapRecGen[i] = new TH2D(Form("histoRapRecGen_%d",i), Form("histoRapRecGen_%d",i), nBinRap, xbinsRap, nBinRap, xbinsRap);
   histoRapRecGen_LepEff[i] = new TH2D(Form("histoRapRecGen_LepEff_%d",i), Form("histoRapRecGen_LepEff_%d",i), nBinRap, xbinsRap, nBinRap, xbinsRap);
   histoRapRecDA[i]  = new TH1D(Form("histoRapRecDA_%d",i),  Form("histoRapRecDA_%d",i),  nBinRap, xbinsRap);
   histoRapRecDY[i]  = new TH1D(Form("histoRapRecDY_%d",i),  Form("histoRapRecDY_%d",i),  nBinRap, xbinsRap);

   histoRapRecDA_MomRes[i] = new TH1D(Form("histoRapRecDA_MomRes_%d",i), Form("histoRapRecDA_MomRes_%d",i), nBinRap, xbinsRap);
   histoRapRecDY_MomRes[i] = new TH1D(Form("histoRapRecDY_MomRes_%d",i), Form("histoRapRecDY_MomRes_%d",i), nBinRap, xbinsRap);

   histoRapRecDA_LepEff[i] = new TH1D(Form("histoRapRecDA_LepEff_%d",i), Form("histoRapRecDA_LepEff_%d",i), nBinRap, xbinsRap);
   histoRapRecDY_LepEff[i] = new TH1D(Form("histoRapRecDY_LepEff_%d",i), Form("histoRapRecDY_LepEff_%d",i), nBinRap, xbinsRap);

   histoRapRecDA_PDF[i] = new TH1D(Form("histoRapRecDA_PDF_%d",i), Form("histoRapRecDA_PDF_%d",i), nBinRap, xbinsRap);
   histoRapRecDY_PDF[i] = new TH1D(Form("histoRapRecDY_PDF_%d",i), Form("histoRapRecDY_PDF_%d",i), nBinRap, xbinsRap);

   histoRapRecDA_QCD[i] = new TH1D(Form("histoRapRecDA_QCD_%d",i), Form("histoRapRecDA_QCD_%d",i), nBinRap, xbinsRap);
   histoRapRecDY_QCD[i] = new TH1D(Form("histoRapRecDY_QCD_%d",i), Form("histoRapRecDY_QCD_%d",i), nBinRap, xbinsRap);

   for(int j=0; j<6; j++){
      histoRapRecDA_QCDPart[i][j] = new TH1D(Form("histoRapRecDA_QCDPart_%d_%d",i,j), Form("histoRapRecDA_QCDPart_%d_%d",i,j), nBinRap, xbinsRap);
      histoRapRecDY_QCDPart[i][j] = new TH1D(Form("histoRapRecDY_QCDPart_%d_%d",i,j), Form("histoRapRecDY_QCDPart_%d_%d",i,j), nBinRap, xbinsRap);
    }
  }

  TH2D* histoPtRap0RecGen[2];
  TH2D* histoPtRap0RecGen_LepEff[2];
  TH1D* histoPtRap0RecDA[2];
  TH1D* histoPtRap0RecDY[2];
  TH1D* histoPtRap0RecDA_MomRes[2];
  TH1D* histoPtRap0RecDY_MomRes[2];
  TH1D* histoPtRap0RecDA_LepEff[2];
  TH1D* histoPtRap0RecDY_LepEff[2];
  TH1D* histoPtRap0RecDA_PDF[2];
  TH1D* histoPtRap0RecDY_PDF[2];
  TH1D* histoPtRap0RecDA_QCD[2];
  TH1D* histoPtRap0RecDY_QCD[2];
  TH1D* histoPtRap0RecDA_QCDPart[2][6];
  TH1D* histoPtRap0RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRap0RecGen[i] = new TH2D(Form("histoPtRap0RecGen_%d",i), Form("histoPtRap0RecGen_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0, nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecGen_LepEff[i] = new TH2D(Form("histoPtRap0RecGen_LepEff_%d",i), Form("histoPtRap0RecGen_LepEff_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0, nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecDA[i]  = new TH1D(Form("histoPtRap0RecDA_%d",i),  Form("histoPtRap0RecDA_%d",i),  nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecDY[i]  = new TH1D(Form("histoPtRap0RecDY_%d",i),  Form("histoPtRap0RecDY_%d",i),  nBinRecoPtRap0, xbinsRecoPtRap0);

   histoPtRap0RecDA_MomRes[i] = new TH1D(Form("histoPtRap0RecDA_MomRes_%d",i), Form("histoPtRap0RecDA_MomRes_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecDY_MomRes[i] = new TH1D(Form("histoPtRap0RecDY_MomRes_%d",i), Form("histoPtRap0RecDY_MomRes_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);

   histoPtRap0RecDA_LepEff[i] = new TH1D(Form("histoPtRap0RecDA_LepEff_%d",i), Form("histoPtRap0RecDA_LepEff_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecDY_LepEff[i] = new TH1D(Form("histoPtRap0RecDY_LepEff_%d",i), Form("histoPtRap0RecDY_LepEff_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);

   histoPtRap0RecDA_PDF[i] = new TH1D(Form("histoPtRap0RecDA_PDF_%d",i), Form("histoPtRap0RecDA_PDF_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecDY_PDF[i] = new TH1D(Form("histoPtRap0RecDY_PDF_%d",i), Form("histoPtRap0RecDY_PDF_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);

   histoPtRap0RecDA_QCD[i] = new TH1D(Form("histoPtRap0RecDA_QCD_%d",i), Form("histoPtRap0RecDA_QCD_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecDY_QCD[i] = new TH1D(Form("histoPtRap0RecDY_QCD_%d",i), Form("histoPtRap0RecDY_QCD_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);

   for(int j=0; j<6; j++){
      histoPtRap0RecDA_QCDPart[i][j] = new TH1D(Form("histoPtRap0RecDA_QCDPart_%d_%d",i,j), Form("histoPtRap0RecDA_QCDPart_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
      histoPtRap0RecDY_QCDPart[i][j] = new TH1D(Form("histoPtRap0RecDY_QCDPart_%d_%d",i,j), Form("histoPtRap0RecDY_QCDPart_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
    }
  }

  TH2D* histoPtRap1RecGen[2];
  TH2D* histoPtRap1RecGen_LepEff[2];
  TH1D* histoPtRap1RecDA[2];
  TH1D* histoPtRap1RecDY[2];
  TH1D* histoPtRap1RecDA_MomRes[2];
  TH1D* histoPtRap1RecDY_MomRes[2];
  TH1D* histoPtRap1RecDA_LepEff[2];
  TH1D* histoPtRap1RecDY_LepEff[2];
  TH1D* histoPtRap1RecDA_PDF[2];
  TH1D* histoPtRap1RecDY_PDF[2];
  TH1D* histoPtRap1RecDA_QCD[2];
  TH1D* histoPtRap1RecDY_QCD[2];
  TH1D* histoPtRap1RecDA_QCDPart[2][6];
  TH1D* histoPtRap1RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRap1RecGen[i] = new TH2D(Form("histoPtRap1RecGen_%d",i), Form("histoPtRap1RecGen_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1, nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecGen_LepEff[i] = new TH2D(Form("histoPtRap1RecGen_LepEff_%d",i), Form("histoPtRap1RecGen_LepEff_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1, nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecDA[i]  = new TH1D(Form("histoPtRap1RecDA_%d",i),  Form("histoPtRap1RecDA_%d",i),  nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecDY[i]  = new TH1D(Form("histoPtRap1RecDY_%d",i),  Form("histoPtRap1RecDY_%d",i),  nBinRecoPtRap1, xbinsRecoPtRap1);

   histoPtRap1RecDA_MomRes[i] = new TH1D(Form("histoPtRap1RecDA_MomRes_%d",i), Form("histoPtRap1RecDA_MomRes_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecDY_MomRes[i] = new TH1D(Form("histoPtRap1RecDY_MomRes_%d",i), Form("histoPtRap1RecDY_MomRes_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);

   histoPtRap1RecDA_LepEff[i] = new TH1D(Form("histoPtRap1RecDA_LepEff_%d",i), Form("histoPtRap1RecDA_LepEff_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecDY_LepEff[i] = new TH1D(Form("histoPtRap1RecDY_LepEff_%d",i), Form("histoPtRap1RecDY_LepEff_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);

   histoPtRap1RecDA_PDF[i] = new TH1D(Form("histoPtRap1RecDA_PDF_%d",i), Form("histoPtRap1RecDA_PDF_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecDY_PDF[i] = new TH1D(Form("histoPtRap1RecDY_PDF_%d",i), Form("histoPtRap1RecDY_PDF_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);

   histoPtRap1RecDA_QCD[i] = new TH1D(Form("histoPtRap1RecDA_QCD_%d",i), Form("histoPtRap1RecDA_QCD_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecDY_QCD[i] = new TH1D(Form("histoPtRap1RecDY_QCD_%d",i), Form("histoPtRap1RecDY_QCD_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);

   for(int j=0; j<6; j++){
      histoPtRap1RecDA_QCDPart[i][j] = new TH1D(Form("histoPtRap1RecDA_QCDPart_%d_%d",i,j), Form("histoPtRap1RecDA_QCDPart_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
      histoPtRap1RecDY_QCDPart[i][j] = new TH1D(Form("histoPtRap1RecDY_QCDPart_%d_%d",i,j), Form("histoPtRap1RecDY_QCDPart_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
    }
  }

  TH2D* histoPtRap2RecGen[2];
  TH2D* histoPtRap2RecGen_LepEff[2];
  TH1D* histoPtRap2RecDA[2];
  TH1D* histoPtRap2RecDY[2];
  TH1D* histoPtRap2RecDA_MomRes[2];
  TH1D* histoPtRap2RecDY_MomRes[2];
  TH1D* histoPtRap2RecDA_LepEff[2];
  TH1D* histoPtRap2RecDY_LepEff[2];
  TH1D* histoPtRap2RecDA_PDF[2];
  TH1D* histoPtRap2RecDY_PDF[2];
  TH1D* histoPtRap2RecDA_QCD[2];
  TH1D* histoPtRap2RecDY_QCD[2];
  TH1D* histoPtRap2RecDA_QCDPart[2][6];
  TH1D* histoPtRap2RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRap2RecGen[i] = new TH2D(Form("histoPtRap2RecGen_%d",i), Form("histoPtRap2RecGen_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2, nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecGen_LepEff[i] = new TH2D(Form("histoPtRap2RecGen_LepEff_%d",i), Form("histoPtRap2RecGen_LepEff_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2, nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecDA[i]  = new TH1D(Form("histoPtRap2RecDA_%d",i),  Form("histoPtRap2RecDA_%d",i),  nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecDY[i]  = new TH1D(Form("histoPtRap2RecDY_%d",i),  Form("histoPtRap2RecDY_%d",i),  nBinRecoPtRap2, xbinsRecoPtRap2);

   histoPtRap2RecDA_MomRes[i] = new TH1D(Form("histoPtRap2RecDA_MomRes_%d",i), Form("histoPtRap2RecDA_MomRes_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecDY_MomRes[i] = new TH1D(Form("histoPtRap2RecDY_MomRes_%d",i), Form("histoPtRap2RecDY_MomRes_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);

   histoPtRap2RecDA_LepEff[i] = new TH1D(Form("histoPtRap2RecDA_LepEff_%d",i), Form("histoPtRap2RecDA_LepEff_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecDY_LepEff[i] = new TH1D(Form("histoPtRap2RecDY_LepEff_%d",i), Form("histoPtRap2RecDY_LepEff_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);

   histoPtRap2RecDA_PDF[i] = new TH1D(Form("histoPtRap2RecDA_PDF_%d",i), Form("histoPtRap2RecDA_PDF_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecDY_PDF[i] = new TH1D(Form("histoPtRap2RecDY_PDF_%d",i), Form("histoPtRap2RecDY_PDF_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);

   histoPtRap2RecDA_QCD[i] = new TH1D(Form("histoPtRap2RecDA_QCD_%d",i), Form("histoPtRap2RecDA_QCD_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecDY_QCD[i] = new TH1D(Form("histoPtRap2RecDY_QCD_%d",i), Form("histoPtRap2RecDY_QCD_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);

   for(int j=0; j<6; j++){
      histoPtRap2RecDA_QCDPart[i][j] = new TH1D(Form("histoPtRap2RecDA_QCDPart_%d_%d",i,j), Form("histoPtRap2RecDA_QCDPart_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
      histoPtRap2RecDY_QCDPart[i][j] = new TH1D(Form("histoPtRap2RecDY_QCDPart_%d_%d",i,j), Form("histoPtRap2RecDY_QCDPart_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
    }
  }

  TH2D* histoPtRap3RecGen[2];
  TH2D* histoPtRap3RecGen_LepEff[2];
  TH1D* histoPtRap3RecDA[2];
  TH1D* histoPtRap3RecDY[2];
  TH1D* histoPtRap3RecDA_MomRes[2];
  TH1D* histoPtRap3RecDY_MomRes[2];
  TH1D* histoPtRap3RecDA_LepEff[2];
  TH1D* histoPtRap3RecDY_LepEff[2];
  TH1D* histoPtRap3RecDA_PDF[2];
  TH1D* histoPtRap3RecDY_PDF[2];
  TH1D* histoPtRap3RecDA_QCD[2];
  TH1D* histoPtRap3RecDY_QCD[2];
  TH1D* histoPtRap3RecDA_QCDPart[2][6];
  TH1D* histoPtRap3RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRap3RecGen[i] = new TH2D(Form("histoPtRap3RecGen_%d",i), Form("histoPtRap3RecGen_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3, nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecGen_LepEff[i] = new TH2D(Form("histoPtRap3RecGen_LepEff_%d",i), Form("histoPtRap3RecGen_LepEff_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3, nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecDA[i]  = new TH1D(Form("histoPtRap3RecDA_%d",i),  Form("histoPtRap3RecDA_%d",i),  nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecDY[i]  = new TH1D(Form("histoPtRap3RecDY_%d",i),  Form("histoPtRap3RecDY_%d",i),  nBinRecoPtRap3, xbinsRecoPtRap3);

   histoPtRap3RecDA_MomRes[i] = new TH1D(Form("histoPtRap3RecDA_MomRes_%d",i), Form("histoPtRap3RecDA_MomRes_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecDY_MomRes[i] = new TH1D(Form("histoPtRap3RecDY_MomRes_%d",i), Form("histoPtRap3RecDY_MomRes_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);

   histoPtRap3RecDA_LepEff[i] = new TH1D(Form("histoPtRap3RecDA_LepEff_%d",i), Form("histoPtRap3RecDA_LepEff_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecDY_LepEff[i] = new TH1D(Form("histoPtRap3RecDY_LepEff_%d",i), Form("histoPtRap3RecDY_LepEff_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);

   histoPtRap3RecDA_PDF[i] = new TH1D(Form("histoPtRap3RecDA_PDF_%d",i), Form("histoPtRap3RecDA_PDF_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecDY_PDF[i] = new TH1D(Form("histoPtRap3RecDY_PDF_%d",i), Form("histoPtRap3RecDY_PDF_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);

   histoPtRap3RecDA_QCD[i] = new TH1D(Form("histoPtRap3RecDA_QCD_%d",i), Form("histoPtRap3RecDA_QCD_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecDY_QCD[i] = new TH1D(Form("histoPtRap3RecDY_QCD_%d",i), Form("histoPtRap3RecDY_QCD_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);

   for(int j=0; j<6; j++){
      histoPtRap3RecDA_QCDPart[i][j] = new TH1D(Form("histoPtRap3RecDA_QCDPart_%d_%d",i,j), Form("histoPtRap3RecDA_QCDPart_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
      histoPtRap3RecDY_QCDPart[i][j] = new TH1D(Form("histoPtRap3RecDY_QCDPart_%d_%d",i,j), Form("histoPtRap3RecDY_QCDPart_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
    }
  }

  TH2D* histoPtRap4RecGen[2];
  TH2D* histoPtRap4RecGen_LepEff[2];
  TH1D* histoPtRap4RecDA[2];
  TH1D* histoPtRap4RecDY[2];
  TH1D* histoPtRap4RecDA_MomRes[2];
  TH1D* histoPtRap4RecDY_MomRes[2];
  TH1D* histoPtRap4RecDA_LepEff[2];
  TH1D* histoPtRap4RecDY_LepEff[2];
  TH1D* histoPtRap4RecDA_PDF[2];
  TH1D* histoPtRap4RecDY_PDF[2];
  TH1D* histoPtRap4RecDA_QCD[2];
  TH1D* histoPtRap4RecDY_QCD[2];
  TH1D* histoPtRap4RecDA_QCDPart[2][6];
  TH1D* histoPtRap4RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRap4RecGen[i] = new TH2D(Form("histoPtRap4RecGen_%d",i), Form("histoPtRap4RecGen_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4, nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecGen_LepEff[i] = new TH2D(Form("histoPtRap4RecGen_LepEff_%d",i), Form("histoPtRap4RecGen_LepEff_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4, nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecDA[i]  = new TH1D(Form("histoPtRap4RecDA_%d",i),  Form("histoPtRap4RecDA_%d",i),  nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecDY[i]  = new TH1D(Form("histoPtRap4RecDY_%d",i),  Form("histoPtRap4RecDY_%d",i),  nBinRecoPtRap4, xbinsRecoPtRap4);

   histoPtRap4RecDA_MomRes[i] = new TH1D(Form("histoPtRap4RecDA_MomRes_%d",i), Form("histoPtRap4RecDA_MomRes_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecDY_MomRes[i] = new TH1D(Form("histoPtRap4RecDY_MomRes_%d",i), Form("histoPtRap4RecDY_MomRes_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);

   histoPtRap4RecDA_LepEff[i] = new TH1D(Form("histoPtRap4RecDA_LepEff_%d",i), Form("histoPtRap4RecDA_LepEff_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecDY_LepEff[i] = new TH1D(Form("histoPtRap4RecDY_LepEff_%d",i), Form("histoPtRap4RecDY_LepEff_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);

   histoPtRap4RecDA_PDF[i] = new TH1D(Form("histoPtRap4RecDA_PDF_%d",i), Form("histoPtRap4RecDA_PDF_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecDY_PDF[i] = new TH1D(Form("histoPtRap4RecDY_PDF_%d",i), Form("histoPtRap4RecDY_PDF_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);

   histoPtRap4RecDA_QCD[i] = new TH1D(Form("histoPtRap4RecDA_QCD_%d",i), Form("histoPtRap4RecDA_QCD_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecDY_QCD[i] = new TH1D(Form("histoPtRap4RecDY_QCD_%d",i), Form("histoPtRap4RecDY_QCD_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);

   for(int j=0; j<6; j++){
      histoPtRap4RecDA_QCDPart[i][j] = new TH1D(Form("histoPtRap4RecDA_QCDPart_%d_%d",i,j), Form("histoPtRap4RecDA_QCDPart_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
      histoPtRap4RecDY_QCDPart[i][j] = new TH1D(Form("histoPtRap4RecDY_QCDPart_%d_%d",i,j), Form("histoPtRap4RecDY_QCDPart_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
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

      if(lepType >= 3 || (lepType == 2 && theCategory != 0)) continue;
      
      if(lepType == 2 && theCategory == 0) theCategory = 5; // using data e-mu events to estimate non-resonant background

      double thePDGMass[2] = {mass_mu, mass_mu};
      if     (abs(lepType) == 1) {thePDGMass[0] = mass_el; thePDGMass[1] = mass_el;}
      else if(abs(lepType) == 2 && abs(thePandaFlat.looseLep1PdgId)==11) {thePDGMass[0] = mass_el;}
      else if(abs(lepType) == 2 && abs(thePandaFlat.looseLep2PdgId)==11) {thePDGMass[1] = mass_el;}
      TLorentzVector v1,v2;
      v1.SetPtEtaPhiM(thePandaFlat.looseLep1Pt*lepPtSF[0],thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      v2.SetPtEtaPhiM(thePandaFlat.looseLep2Pt*lepPtSF[1],thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      TLorentzVector vMomRes1,vMomRes2;
      vMomRes1.SetPtEtaPhiM(thePandaFlat.looseLep1Pt*sqrt(lepPtSFSyst[0]*lepPtSFSyst[0]*lepPtSFSyst[0]),thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      vMomRes2.SetPtEtaPhiM(thePandaFlat.looseLep2Pt*sqrt(lepPtSFSyst[1]*lepPtSFSyst[1]*lepPtSFSyst[1]),thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);

      bool passSel = TMath::Abs((v1+v2).M()-91.1876) < 15 && v1.Pt() > 25 && v2.Pt() > 25;
      bool passSystSel[1] = {TMath::Abs((vMomRes1+vMomRes2).M()-91.1876) < 15 && vMomRes1.Pt() > 25 && vMomRes2.Pt() > 25};

      double ZRecPt  = (v1+v2).Pt();
      double ZRecRap = TMath::Abs((v1+v2).Rapidity());

      double ZRecSystPt[1] = {(vMomRes1+vMomRes2).Pt()};
      double ZGenPt = 0; double ZGenRap = 0; bool passPtFid = false; bool passRapFid = false; bool passPtRapFid[5] = {false, false, false, false, false}; 
      if(thePandaFlat.looseGenLep1PdgId != 0 && thePandaFlat.looseGenLep2PdgId != 0 &&
         thePandaFlat.genLep1Pt > 25 && TMath::Abs(thePandaFlat.genLep1Eta) < 2.5 &&
	 thePandaFlat.genLep2Pt > 25 && TMath::Abs(thePandaFlat.genLep2Eta) < 2.5){
        TLorentzVector vGen1,vGen2;
        vGen1.SetPtEtaPhiM(thePandaFlat.genLep1Pt,thePandaFlat.genLep1Eta,thePandaFlat.genLep1Phi,thePDGMass[0]);
        vGen2.SetPtEtaPhiM(thePandaFlat.genLep2Pt,thePandaFlat.genLep2Eta,thePandaFlat.genLep2Phi,thePDGMass[1]);
	if(TMath::Abs((vGen1+vGen2).M()-91.1876) < 15.0) {
	  ZGenPt = (vGen1+vGen2).Pt(); passPtFid = true;
          ZGenRap = TMath::Abs((vGen1+vGen2).Rapidity());
          if(ZGenRap < 2.4) { passRapFid = true;}
          if     (ZGenRap < 0.5) passPtRapFid[0] = true;
          else if(ZGenRap < 1.0) passPtRapFid[1] = true;
          else if(ZGenRap < 1.5) passPtRapFid[2] = true;
          else if(ZGenRap < 2.0) passPtRapFid[3] = true;
          else if(ZGenRap < 2.4) passPtRapFid[4] = true;
	}
      }

      double totalWeight = 1.0;
      double totalWeightLepEff = 1.0;
      if(theCategory != 0 && theCategory != 5){
        totalWeight = thePandaFlat.normalizedWeight * lumi * thePandaFlat.sf_pu *
	              thePandaFlat.sf_trk1 * thePandaFlat.sf_medium1 *
		      thePandaFlat.sf_trk2 * thePandaFlat.sf_medium2;
        totalWeightLepEff = thePandaFlat.normalizedWeight * lumi * thePandaFlat.sf_pu *
	              thePandaFlat.sf_trk1 * thePandaFlat.sf_medium1 * (1.0+thePandaFlat.sf_unc1) *
		      thePandaFlat.sf_trk2 * thePandaFlat.sf_medium2 * (1.0+thePandaFlat.sf_unc2);
      }
      
      if(passSel){

        if(theCategory != 5){
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
	  histo[lepType+18][theCategory]->Fill(ZRecRap,totalWeight);
        }
        else {
          for(int theLepType = 0; theLepType<2; theLepType++) {
            double theKeff = -k_eff;
            if(theLepType == 1) theKeff = -1.0/k_eff;
 	    histo[theLepType+0][theCategory]->Fill((v1+v2).M(),totalWeight*theKeff);
	    histo[theLepType+2][theCategory]->Fill(TMath::Min(ZRecPt, 199.999),totalWeight*theKeff);
	    histo[theLepType+4][theCategory]->Fill(TMath::Min(ZRecPt,1999.999),totalWeight*theKeff);
	    histo[theLepType+6][theCategory]->Fill(TMath::Abs(ZRecPt-ZGenPt),totalWeight*theKeff);
            if     (TMath::Abs(thePandaFlat.looseLep1Eta) <  1.5 && TMath::Abs(thePandaFlat.looseLep2Eta) <  1.5){
	      histo[theLepType+ 8][theCategory]->Fill((v1+v2).M(),totalWeight*theKeff);
            }
            else if(TMath::Abs(thePandaFlat.looseLep1Eta) >= 1.5 && TMath::Abs(thePandaFlat.looseLep2Eta) >= 1.5){
	      histo[theLepType+12][theCategory]->Fill((v1+v2).M(),totalWeight*theKeff);
            }
            else {
	      histo[theLepType+10][theCategory]->Fill((v1+v2).M(),totalWeight*theKeff);
            }
	    histo[theLepType+14][theCategory]->Fill((v1+v2).Px(),totalWeight*theKeff);
	    histo[theLepType+16][theCategory]->Fill((v1+v2).Py(),totalWeight*theKeff);
	    histo[theLepType+18][theCategory]->Fill(ZRecRap,totalWeight*theKeff);
          }
        }

        double maxQCDscale = 1.0;
        if(theCategory != 0 && theCategory != 5 && thePandaFlat.scale[0] != -1){
          maxQCDscale = (TMath::Abs(1+thePandaFlat.scale[0])+TMath::Abs(1+thePandaFlat.scale[1])+TMath::Abs(1+thePandaFlat.scale[2])+
        		 TMath::Abs(1+thePandaFlat.scale[3])+TMath::Abs(1+thePandaFlat.scale[4])+TMath::Abs(1+thePandaFlat.scale[4]))/6.0;
        }

	if     (theCategory == 2){ // DY
          // Pt
          histoPtRecDY        [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRecDY_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRecDY_PDF    [lepType]   ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPtRecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          // Pt2
          histoPt2RecDY        [lepType]   ->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDY_LepEff [lepType]   ->Fill(ZRecPt*ZRecPt,totalWeightLepEff);
          histoPt2RecDY_PDF    [lepType]   ->Fill(ZRecPt*ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPt2RecDY_QCDPart[lepType][0]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPt2RecDY_QCDPart[lepType][1]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPt2RecDY_QCDPart[lepType][2]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPt2RecDY_QCDPart[lepType][3]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPt2RecDY_QCDPart[lepType][4]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPt2RecDY_QCDPart[lepType][5]->Fill(ZRecPt*ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          // Rap
          if(ZRecRap < 2.4){
          histoRapRecDY        [lepType]   ->Fill(ZRecRap,totalWeight);
          histoRapRecDY_LepEff [lepType]   ->Fill(ZRecRap,totalWeightLepEff);
          histoRapRecDY_PDF    [lepType]   ->Fill(ZRecRap,totalWeight*thePandaFlat.pdfUp);
          histoRapRecDY_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoRapRecDY_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoRapRecDY_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoRapRecDY_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoRapRecDY_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoRapRecDY_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
          // PtRap0
          if      (ZRecRap < 0.5){
          histoPtRap0RecDY        [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRap0RecDY_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRap0RecDY_PDF    [lepType]   ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPtRap0RecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPtRap0RecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPtRap0RecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPtRap0RecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPtRap0RecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPtRap0RecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
          else if(ZRecRap < 1.0){
          histoPtRap1RecDY        [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRap1RecDY_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRap1RecDY_PDF    [lepType]   ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPtRap1RecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPtRap1RecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPtRap1RecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPtRap1RecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPtRap1RecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPtRap1RecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
          else if(ZRecRap < 1.5){
          histoPtRap2RecDY        [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRap2RecDY_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRap2RecDY_PDF    [lepType]   ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPtRap2RecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPtRap2RecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPtRap2RecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPtRap2RecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPtRap2RecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPtRap2RecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
          else if(ZRecRap < 2.0){
          histoPtRap3RecDY        [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRap3RecDY_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRap3RecDY_PDF    [lepType]   ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPtRap3RecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPtRap3RecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPtRap3RecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPtRap3RecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPtRap3RecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPtRap3RecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
          else if(ZRecRap < 2.4){
          histoPtRap4RecDY        [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRap4RecDY_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRap4RecDY_PDF    [lepType]   ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPtRap4RecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPtRap4RecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPtRap4RecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPtRap4RecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPtRap4RecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPtRap4RecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
	}
	else if(theCategory == 0){ // Data
          // Pt
          histoPtRecDA[lepType]           ->Fill(ZRecPt,totalWeight);
          histoPtRecDA_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRecDA_PDF    [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
          histoPtRecDA_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
          // Pt2
          histoPt2RecDA[lepType]           ->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_LepEff [lepType]   ->Fill(ZRecPt*ZRecPt,totalWeightLepEff);
          histoPt2RecDA_PDF    [lepType]   ->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][0]->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][1]->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][2]->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][3]->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][4]->Fill(ZRecPt*ZRecPt,totalWeight);
          histoPt2RecDA_QCDPart[lepType][5]->Fill(ZRecPt*ZRecPt,totalWeight);
          // Rap
          if(ZRecRap < 2.4){
          histoRapRecDA[lepType]           ->Fill(ZRecRap,totalWeight);
          histoRapRecDA_LepEff [lepType]   ->Fill(ZRecRap,totalWeightLepEff);
          histoRapRecDA_PDF    [lepType]   ->Fill(ZRecRap,totalWeight);
          histoRapRecDA_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight);
          histoRapRecDA_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight);
          histoRapRecDA_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight);
          histoRapRecDA_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight);
          histoRapRecDA_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight);
          histoRapRecDA_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight);
          }
          // PtRap0
          if     (ZRecRap < 0.5){
          histoPtRap0RecDA[lepType]           ->Fill(ZRecPt,totalWeight);
          histoPtRap0RecDA_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRap0RecDA_PDF    [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRap0RecDA_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
          histoPtRap0RecDA_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
          histoPtRap0RecDA_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
          histoPtRap0RecDA_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
          histoPtRap0RecDA_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
          histoPtRap0RecDA_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
          }
          else if(ZRecRap < 1.0){
          histoPtRap1RecDA[lepType]           ->Fill(ZRecPt,totalWeight);
          histoPtRap1RecDA_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRap1RecDA_PDF    [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRap1RecDA_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
          histoPtRap1RecDA_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
          histoPtRap1RecDA_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
          histoPtRap1RecDA_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
          histoPtRap1RecDA_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
          histoPtRap1RecDA_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
          }
          else if(ZRecRap < 1.5){
          histoPtRap2RecDA[lepType]           ->Fill(ZRecPt,totalWeight);
          histoPtRap2RecDA_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRap2RecDA_PDF    [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRap2RecDA_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
          histoPtRap2RecDA_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
          histoPtRap2RecDA_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
          histoPtRap2RecDA_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
          histoPtRap2RecDA_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
          histoPtRap2RecDA_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
          }
          else if(ZRecRap < 2.0){
          histoPtRap3RecDA[lepType]           ->Fill(ZRecPt,totalWeight);
          histoPtRap3RecDA_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRap3RecDA_PDF    [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRap3RecDA_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
          histoPtRap3RecDA_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
          histoPtRap3RecDA_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
          histoPtRap3RecDA_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
          histoPtRap3RecDA_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
          histoPtRap3RecDA_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
          }
          else if(ZRecRap < 2.4){
          histoPtRap4RecDA[lepType]           ->Fill(ZRecPt,totalWeight);
          histoPtRap4RecDA_LepEff [lepType]   ->Fill(ZRecPt,totalWeightLepEff);
          histoPtRap4RecDA_PDF    [lepType]   ->Fill(ZRecPt,totalWeight);
          histoPtRap4RecDA_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
          histoPtRap4RecDA_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
          histoPtRap4RecDA_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
          histoPtRap4RecDA_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
          histoPtRap4RecDA_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
          histoPtRap4RecDA_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
          }
	}
	else if(theCategory == 5){ // e-mu Data
          for(int theLepType = 0; theLepType<2; theLepType++) {
            double theKeff = k_eff;
            if(theLepType == 1) theKeff = 1.0/k_eff;
            // Pt
            histoPtRecDA[theLepType]           ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRecDA_LepEff [theLepType]   ->Fill(ZRecPt,totalWeightLepEff*theKeff);
            histoPtRecDA_PDF    [theLepType]   ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRecDA_QCDPart[theLepType][0]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRecDA_QCDPart[theLepType][1]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRecDA_QCDPart[theLepType][2]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRecDA_QCDPart[theLepType][3]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRecDA_QCDPart[theLepType][4]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRecDA_QCDPart[theLepType][5]->Fill(ZRecPt,totalWeight*theKeff);
            // Pt2
            histoPt2RecDA[theLepType]           ->Fill(ZRecPt*ZRecPt,totalWeight*theKeff);
            histoPt2RecDA_LepEff [theLepType]   ->Fill(ZRecPt*ZRecPt,totalWeightLepEff*theKeff);
            histoPt2RecDA_PDF    [theLepType]   ->Fill(ZRecPt*ZRecPt,totalWeight*theKeff);
            histoPt2RecDA_QCDPart[theLepType][0]->Fill(ZRecPt*ZRecPt,totalWeight*theKeff);
            histoPt2RecDA_QCDPart[theLepType][1]->Fill(ZRecPt*ZRecPt,totalWeight*theKeff);
            histoPt2RecDA_QCDPart[theLepType][2]->Fill(ZRecPt*ZRecPt,totalWeight*theKeff);
            histoPt2RecDA_QCDPart[theLepType][3]->Fill(ZRecPt*ZRecPt,totalWeight*theKeff);
            histoPt2RecDA_QCDPart[theLepType][4]->Fill(ZRecPt*ZRecPt,totalWeight*theKeff);
            histoPt2RecDA_QCDPart[theLepType][5]->Fill(ZRecPt*ZRecPt,totalWeight*theKeff);
            // Rap
            if(ZRecRap < 2.4){
            histoRapRecDA[theLepType]           ->Fill(ZRecRap,totalWeight*theKeff);
            histoRapRecDA_LepEff [theLepType]   ->Fill(ZRecRap,totalWeightLepEff*theKeff);
            histoRapRecDA_PDF    [theLepType]   ->Fill(ZRecRap,totalWeight*theKeff);
            histoRapRecDA_QCDPart[theLepType][0]->Fill(ZRecRap,totalWeight*theKeff);
            histoRapRecDA_QCDPart[theLepType][1]->Fill(ZRecRap,totalWeight*theKeff);
            histoRapRecDA_QCDPart[theLepType][2]->Fill(ZRecRap,totalWeight*theKeff);
            histoRapRecDA_QCDPart[theLepType][3]->Fill(ZRecRap,totalWeight*theKeff);
            histoRapRecDA_QCDPart[theLepType][4]->Fill(ZRecRap,totalWeight*theKeff);
            histoRapRecDA_QCDPart[theLepType][5]->Fill(ZRecRap,totalWeight*theKeff);
            }
            // PtRap0
            if     (ZRecRap < 0.5){
            histoPtRap0RecDA[theLepType]           ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap0RecDA_LepEff [theLepType]   ->Fill(ZRecPt,totalWeightLepEff*theKeff);
            histoPtRap0RecDA_PDF    [theLepType]   ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap0RecDA_QCDPart[theLepType][0]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap0RecDA_QCDPart[theLepType][1]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap0RecDA_QCDPart[theLepType][2]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap0RecDA_QCDPart[theLepType][3]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap0RecDA_QCDPart[theLepType][4]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap0RecDA_QCDPart[theLepType][5]->Fill(ZRecPt,totalWeight*theKeff);
            }
            else if(ZRecRap < 1.0){
            histoPtRap1RecDA[theLepType]           ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap1RecDA_LepEff [theLepType]   ->Fill(ZRecPt,totalWeightLepEff*theKeff);
            histoPtRap1RecDA_PDF    [theLepType]   ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap1RecDA_QCDPart[theLepType][0]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap1RecDA_QCDPart[theLepType][1]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap1RecDA_QCDPart[theLepType][2]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap1RecDA_QCDPart[theLepType][3]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap1RecDA_QCDPart[theLepType][4]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap1RecDA_QCDPart[theLepType][5]->Fill(ZRecPt,totalWeight*theKeff);
            }
            else if(ZRecRap < 1.5){
            histoPtRap2RecDA[theLepType]           ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap2RecDA_LepEff [theLepType]   ->Fill(ZRecPt,totalWeightLepEff*theKeff);
            histoPtRap2RecDA_PDF    [theLepType]   ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap2RecDA_QCDPart[theLepType][0]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap2RecDA_QCDPart[theLepType][1]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap2RecDA_QCDPart[theLepType][2]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap2RecDA_QCDPart[theLepType][3]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap2RecDA_QCDPart[theLepType][4]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap2RecDA_QCDPart[theLepType][5]->Fill(ZRecPt,totalWeight*theKeff);
            }
            else if(ZRecRap < 2.0){
            histoPtRap3RecDA[theLepType]           ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap3RecDA_LepEff [theLepType]   ->Fill(ZRecPt,totalWeightLepEff*theKeff);
            histoPtRap3RecDA_PDF    [theLepType]   ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap3RecDA_QCDPart[theLepType][0]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap3RecDA_QCDPart[theLepType][1]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap3RecDA_QCDPart[theLepType][2]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap3RecDA_QCDPart[theLepType][3]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap3RecDA_QCDPart[theLepType][4]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap3RecDA_QCDPart[theLepType][5]->Fill(ZRecPt,totalWeight*theKeff);
            }
            else if(ZRecRap < 2.4){
            histoPtRap4RecDA[theLepType]           ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap4RecDA_LepEff [theLepType]   ->Fill(ZRecPt,totalWeightLepEff*theKeff);
            histoPtRap4RecDA_PDF    [theLepType]   ->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap4RecDA_QCDPart[theLepType][0]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap4RecDA_QCDPart[theLepType][1]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap4RecDA_QCDPart[theLepType][2]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap4RecDA_QCDPart[theLepType][3]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap4RecDA_QCDPart[theLepType][4]->Fill(ZRecPt,totalWeight*theKeff);
            histoPtRap4RecDA_QCDPart[theLepType][5]->Fill(ZRecPt,totalWeight*theKeff);
            }
          } // loop over muons and electrons
	}
	else { // MC
          // Pt
          histoPtRecDA        [lepType]     ->Fill(ZRecPt,-1.0*totalWeight);
          histoPtRecDA_LepEff [lepType]     ->Fill(ZRecPt,-1.0*totalWeightLepEff);
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
          histoPt2RecDA_LepEff [lepType]     ->Fill(ZRecPt*ZRecPt,-1.0*totalWeightLepEff);
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
          // Rap
          if(ZRecRap < 2.4){
          histoRapRecDA        [lepType]     ->Fill(ZRecRap,-1.0*totalWeight);
          histoRapRecDA_LepEff [lepType]     ->Fill(ZRecRap,-1.0*totalWeightLepEff);
          histoRapRecDA_PDF    [lepType]     ->Fill(ZRecRap,-1.0*totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoRapRecDA_QCDPart[lepType][0]->Fill(ZRecRap,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoRapRecDA_QCDPart[lepType][1]->Fill(ZRecRap,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoRapRecDA_QCDPart[lepType][2]->Fill(ZRecRap,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoRapRecDA_QCDPart[lepType][3]->Fill(ZRecRap,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoRapRecDA_QCDPart[lepType][4]->Fill(ZRecRap,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoRapRecDA_QCDPart[lepType][5]->Fill(ZRecRap,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoRapRecDA_QCDPart[lepType][0]->Fill(ZRecRap,-1.0*totalWeight);
            histoRapRecDA_QCDPart[lepType][1]->Fill(ZRecRap,-1.0*totalWeight);
            histoRapRecDA_QCDPart[lepType][2]->Fill(ZRecRap,-1.0*totalWeight);
            histoRapRecDA_QCDPart[lepType][3]->Fill(ZRecRap,-1.0*totalWeight);
            histoRapRecDA_QCDPart[lepType][4]->Fill(ZRecRap,-1.0*totalWeight);
            histoRapRecDA_QCDPart[lepType][5]->Fill(ZRecRap,-1.0*totalWeight);
          }
          }
          // PtRap0
          if     (ZRecRap < 0.5){
          histoPtRap0RecDA        [lepType]     ->Fill(ZRecPt,-1.0*totalWeight);
          histoPtRap0RecDA_LepEff [lepType]     ->Fill(ZRecPt,-1.0*totalWeightLepEff);
          histoPtRap0RecDA_PDF    [lepType]     ->Fill(ZRecPt,-1.0*totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoPtRap0RecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRap0RecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRap0RecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRap0RecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRap0RecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRap0RecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoPtRap0RecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap0RecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap0RecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap0RecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap0RecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap0RecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight);
          }
          }
          else if(ZRecRap < 1.0){
          histoPtRap1RecDA        [lepType]     ->Fill(ZRecPt,-1.0*totalWeight);
          histoPtRap1RecDA_LepEff [lepType]     ->Fill(ZRecPt,-1.0*totalWeightLepEff);
          histoPtRap1RecDA_PDF    [lepType]     ->Fill(ZRecPt,-1.0*totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoPtRap1RecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRap1RecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRap1RecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRap1RecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRap1RecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRap1RecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoPtRap1RecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap1RecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap1RecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap1RecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap1RecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap1RecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight);
          }
          }
          else if(ZRecRap < 1.5){
          histoPtRap2RecDA        [lepType]     ->Fill(ZRecPt,-1.0*totalWeight);
          histoPtRap2RecDA_LepEff [lepType]     ->Fill(ZRecPt,-1.0*totalWeightLepEff);
          histoPtRap2RecDA_PDF    [lepType]     ->Fill(ZRecPt,-1.0*totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoPtRap2RecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRap2RecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRap2RecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRap2RecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRap2RecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRap2RecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoPtRap2RecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap2RecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap2RecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap2RecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap2RecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap2RecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight);
          }
          }
          else if(ZRecRap < 2.0){
          histoPtRap3RecDA        [lepType]     ->Fill(ZRecPt,-1.0*totalWeight);
          histoPtRap3RecDA_LepEff [lepType]     ->Fill(ZRecPt,-1.0*totalWeightLepEff);
          histoPtRap3RecDA_PDF    [lepType]     ->Fill(ZRecPt,-1.0*totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoPtRap3RecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRap3RecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRap3RecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRap3RecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRap3RecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRap3RecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoPtRap3RecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap3RecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap3RecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap3RecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap3RecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap3RecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight);
          }
          }
          else if(ZRecRap < 2.4){
          histoPtRap4RecDA        [lepType]     ->Fill(ZRecPt,-1.0*totalWeight);
          histoPtRap4RecDA_LepEff [lepType]     ->Fill(ZRecPt,-1.0*totalWeightLepEff);
          histoPtRap4RecDA_PDF    [lepType]     ->Fill(ZRecPt,-1.0*totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoPtRap4RecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRap4RecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRap4RecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRap4RecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRap4RecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRap4RecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoPtRap4RecDA_QCDPart[lepType][0]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap4RecDA_QCDPart[lepType][1]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap4RecDA_QCDPart[lepType][2]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap4RecDA_QCDPart[lepType][3]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap4RecDA_QCDPart[lepType][4]->Fill(ZRecPt,-1.0*totalWeight);
            histoPtRap4RecDA_QCDPart[lepType][5]->Fill(ZRecPt,-1.0*totalWeight);
          }
          }
	}

	if(theCategory == 2 && passPtFid == true){
          histoPtRecGen       [lepType] ->Fill(ZRecPt       ,ZGenPt       ,totalWeight);
          histoPtRecGen_LepEff[lepType] ->Fill(ZRecPt       ,ZGenPt       ,totalWeightLepEff);
          histoPt2RecGen      [lepType] ->Fill(ZRecPt*ZRecPt,ZGenPt*ZGenPt,totalWeight);
          histoPt2RecGen_LepEff[lepType]->Fill(ZRecPt*ZRecPt,ZGenPt*ZGenPt,totalWeightLepEff);
	}
	if(theCategory == 2 && passRapFid == true && ZRecRap < 2.4){
          histoRapRecGen       [lepType]->Fill(ZRecRap,ZGenRap,totalWeight);
          histoRapRecGen_LepEff[lepType]->Fill(ZRecRap,ZGenRap,totalWeightLepEff);
	}
	if(theCategory == 2 && passPtRapFid[0] == true && ZRecRap >= 0.0 && ZRecRap < 0.5){
          histoPtRap0RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          histoPtRap0RecGen_LepEff[lepType]->Fill(ZRecPt,ZGenPt,totalWeightLepEff);
	}
	if(theCategory == 2 && passPtRapFid[1] == true && ZRecRap >= 0.5 && ZRecRap < 1.0){
          histoPtRap1RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          histoPtRap1RecGen_LepEff[lepType]->Fill(ZRecPt,ZGenPt,totalWeightLepEff);
	}
	if(theCategory == 2 && passPtRapFid[2] == true && ZRecRap >= 1.0 && ZRecRap < 1.5){
          histoPtRap2RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          histoPtRap2RecGen_LepEff[lepType]->Fill(ZRecPt,ZGenPt,totalWeightLepEff);
	}
	if(theCategory == 2 && passPtRapFid[3] == true && ZRecRap >= 1.5 && ZRecRap < 2.0){
          histoPtRap3RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          histoPtRap3RecGen_LepEff[lepType]->Fill(ZRecPt,ZGenPt,totalWeightLepEff);
	}
	if(theCategory == 2 && passPtRapFid[4] == true && ZRecRap >= 2.0 && ZRecRap < 2.4){
          histoPtRap4RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          histoPtRap4RecGen_LepEff[lepType]->Fill(ZRecPt,ZGenPt,totalWeightLepEff);
	}
      }
      
      if(passSystSel[0]){
	if     (theCategory == 2){
          histoPtRecDY_MomRes [lepType]->Fill(ZRecSystPt[0],totalWeight);
          histoPt2RecDY_MomRes[lepType]->Fill(ZRecSystPt[0]*ZRecSystPt[0],totalWeight);
          if(ZRecRap < 2.4) histoRapRecDY_MomRes[lepType]->Fill(ZRecRap,totalWeight);
          if     (ZRecRap < 0.5) histoPtRap0RecDY_MomRes[lepType]->Fill(ZRecSystPt[0],totalWeight);
          else if(ZRecRap < 1.0) histoPtRap1RecDY_MomRes[lepType]->Fill(ZRecSystPt[0],totalWeight);
          else if(ZRecRap < 1.5) histoPtRap2RecDY_MomRes[lepType]->Fill(ZRecSystPt[0],totalWeight);
          else if(ZRecRap < 2.0) histoPtRap3RecDY_MomRes[lepType]->Fill(ZRecSystPt[0],totalWeight);
          else if(ZRecRap < 2.4) histoPtRap4RecDY_MomRes[lepType]->Fill(ZRecSystPt[0],totalWeight);
	}
	else if(theCategory == 0){
          histoPtRecDA_MomRes[lepType] ->Fill(ZRecSystPt[0],totalWeight);
          histoPt2RecDA_MomRes[lepType]->Fill(ZRecSystPt[0]*ZRecSystPt[0],totalWeight);
          if(ZRecRap < 2.4) histoRapRecDA_MomRes[lepType]->Fill(ZRecRap,totalWeight);
          if     (ZRecRap < 0.5) histoPtRap0RecDA_MomRes[lepType]->Fill(ZRecSystPt[0],totalWeight);
          else if(ZRecRap < 1.0) histoPtRap1RecDA_MomRes[lepType]->Fill(ZRecSystPt[0],totalWeight);
          else if(ZRecRap < 1.5) histoPtRap2RecDA_MomRes[lepType]->Fill(ZRecSystPt[0],totalWeight);
          else if(ZRecRap < 2.0) histoPtRap3RecDA_MomRes[lepType]->Fill(ZRecSystPt[0],totalWeight);
          else if(ZRecRap < 2.4) histoPtRap4RecDA_MomRes[lepType]->Fill(ZRecSystPt[0],totalWeight);
	}
	else if(theCategory == 5){
          for(int theLepType = 0; theLepType<2; theLepType++) {
            double theKeff = k_eff;
            if(theLepType == 1) theKeff = 1.0/k_eff;
            histoPtRecDA_MomRes[theLepType] ->Fill(ZRecSystPt[0],totalWeight*theKeff);
            histoPt2RecDA_MomRes[theLepType]->Fill(ZRecSystPt[0]*ZRecSystPt[0],totalWeight*theKeff);
            if(ZRecRap < 2.4) histoRapRecDA_MomRes[theLepType]->Fill(ZRecRap,totalWeight*theKeff);
            if     (ZRecRap < 0.5) histoPtRap0RecDA_MomRes[theLepType]->Fill(ZRecSystPt[0],totalWeight*theKeff);
            else if(ZRecRap < 1.0) histoPtRap1RecDA_MomRes[theLepType]->Fill(ZRecSystPt[0],totalWeight*theKeff);
            else if(ZRecRap < 1.5) histoPtRap2RecDA_MomRes[theLepType]->Fill(ZRecSystPt[0],totalWeight*theKeff);
            else if(ZRecRap < 2.0) histoPtRap3RecDA_MomRes[theLepType]->Fill(ZRecSystPt[0],totalWeight*theKeff);
            else if(ZRecRap < 2.4) histoPtRap4RecDA_MomRes[theLepType]->Fill(ZRecSystPt[0],totalWeight*theKeff);
          }
	}
	else {
          histoPtRecDA_MomRes[lepType] ->Fill(ZRecSystPt[0],-1.0*totalWeight);
          histoPt2RecDA_MomRes[lepType]->Fill(ZRecSystPt[0]*ZRecSystPt[0],-1.0*totalWeight);
          if(ZRecRap < 2.4) histoRapRecDA_MomRes[lepType]->Fill(ZRecRap,-1.0*totalWeight);
          if     (ZRecRap < 0.5) histoPtRap0RecDA_MomRes[lepType]->Fill(ZRecSystPt[0],-1.0*totalWeight);
          else if(ZRecRap < 1.0) histoPtRap1RecDA_MomRes[lepType]->Fill(ZRecSystPt[0],-1.0*totalWeight);
          else if(ZRecRap < 1.5) histoPtRap2RecDA_MomRes[lepType]->Fill(ZRecSystPt[0],-1.0*totalWeight);
          else if(ZRecRap < 2.0) histoPtRap3RecDA_MomRes[lepType]->Fill(ZRecSystPt[0],-1.0*totalWeight);
          else if(ZRecRap < 2.4) histoPtRap4RecDA_MomRes[lepType]->Fill(ZRecSystPt[0],-1.0*totalWeight);
	}
      }
    } // end event loop
  } // end samples loop

  // Pt
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPt(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRecDA_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRecDA_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRecDA_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRecDA_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRecDA_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRecDA_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRecDA[ntype]->GetSumOfWeights(),
           histoPtRecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPt+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRecDA_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRecDA[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecDA[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecDA[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRecDA[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRecDA[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRecDA_QCD[ntype]->SetBinContent(nb, histoPtRecDA[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRecDY_QCD[ntype]->SetBinContent(nb, histoPtRecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // Pt2
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPt2(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPt2RecDA_QCDPart[ntype][0]->GetSumOfWeights(),histoPt2RecDA_QCDPart[ntype][1]->GetSumOfWeights(),histoPt2RecDA_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPt2RecDA_QCDPart[ntype][3]->GetSumOfWeights(),histoPt2RecDA_QCDPart[ntype][4]->GetSumOfWeights(),histoPt2RecDA_QCDPart[ntype][5]->GetSumOfWeights(),histoPt2RecDA[ntype]->GetSumOfWeights(),
           histoPt2RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPt2RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPt2RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPt2RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPt2RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPt2RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPt2RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinPt2+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPt2RecDA_QCDPart[ntype][0]->GetBinContent(nb)-histoPt2RecDA[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPt2RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPt2RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPt2RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPt2RecDA[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPt2RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPt2RecDA[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPt2RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPt2RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPt2RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPt2RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPt2RecDA[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPt2RecDA[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPt2RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPt2RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPt2RecDA_QCD[ntype]->SetBinContent(nb, histoPt2RecDA[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPt2RecDY_QCD[ntype]->SetBinContent(nb, histoPt2RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // Rap
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDRap(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoRapRecDA_QCDPart[ntype][0]->GetSumOfWeights(),histoRapRecDA_QCDPart[ntype][1]->GetSumOfWeights(),histoRapRecDA_QCDPart[ntype][2]->GetSumOfWeights(),
           histoRapRecDA_QCDPart[ntype][3]->GetSumOfWeights(),histoRapRecDA_QCDPart[ntype][4]->GetSumOfWeights(),histoRapRecDA_QCDPart[ntype][5]->GetSumOfWeights(),histoRapRecDA[ntype]->GetSumOfWeights(),
           histoRapRecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoRapRecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoRapRecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoRapRecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoRapRecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoRapRecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoRapRecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRap+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoRapRecDA_QCDPart[ntype][0]->GetBinContent(nb)-histoRapRecDA[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoRapRecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoRapRecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoRapRecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapRecDA[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoRapRecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapRecDA[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoRapRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapRecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoRapRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapRecDY[ntype]->GetBinContent(nb));
      }
      if(histoRapRecDA[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoRapRecDA[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoRapRecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoRapRecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoRapRecDA_QCD[ntype]->SetBinContent(nb, histoRapRecDA[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoRapRecDY_QCD[ntype]->SetBinContent(nb, histoRapRecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PtRap0
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPtRap0(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRap0RecDA_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap0RecDA_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap0RecDA_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap0RecDA_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap0RecDA_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap0RecDA_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap0RecDA[ntype]->GetSumOfWeights(),
           histoPtRap0RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap0RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap0RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap0RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap0RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap0RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap0RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPtRap0+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRap0RecDA_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap0RecDA[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRap0RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap0RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRap0RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap0RecDA[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRap0RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap0RecDA[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRap0RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap0RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRap0RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap0RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRap0RecDA[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRap0RecDA[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRap0RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRap0RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRap0RecDA_QCD[ntype]->SetBinContent(nb, histoPtRap0RecDA[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRap0RecDY_QCD[ntype]->SetBinContent(nb, histoPtRap0RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PtRap1
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPtRap1(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRap1RecDA_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap1RecDA_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap1RecDA_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap1RecDA_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap1RecDA_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap1RecDA_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap1RecDA[ntype]->GetSumOfWeights(),
           histoPtRap1RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap1RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap1RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap1RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap1RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap1RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap1RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPtRap1+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRap1RecDA_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap1RecDA[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRap1RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap1RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRap1RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap1RecDA[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRap1RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap1RecDA[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRap1RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap1RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRap1RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap1RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRap1RecDA[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRap1RecDA[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRap1RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRap1RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRap1RecDA_QCD[ntype]->SetBinContent(nb, histoPtRap1RecDA[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRap1RecDY_QCD[ntype]->SetBinContent(nb, histoPtRap1RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PtRap2
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPtRap2(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRap2RecDA_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap2RecDA_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap2RecDA_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap2RecDA_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap2RecDA_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap2RecDA_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap2RecDA[ntype]->GetSumOfWeights(),
           histoPtRap2RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap2RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap2RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap2RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap2RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap2RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap2RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPtRap2+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRap2RecDA_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap2RecDA[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRap2RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap2RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRap2RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap2RecDA[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRap2RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap2RecDA[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRap2RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap2RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRap2RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap2RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRap2RecDA[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRap2RecDA[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRap2RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRap2RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRap2RecDA_QCD[ntype]->SetBinContent(nb, histoPtRap2RecDA[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRap2RecDY_QCD[ntype]->SetBinContent(nb, histoPtRap2RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PtRap3
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPtRap3(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRap3RecDA_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap3RecDA_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap3RecDA_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap3RecDA_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap3RecDA_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap3RecDA_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap3RecDA[ntype]->GetSumOfWeights(),
           histoPtRap3RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap3RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap3RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap3RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap3RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap3RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap3RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPtRap3+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRap3RecDA_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap3RecDA[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRap3RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap3RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRap3RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap3RecDA[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRap3RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap3RecDA[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRap3RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap3RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRap3RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap3RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRap3RecDA[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRap3RecDA[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRap3RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRap3RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRap3RecDA_QCD[ntype]->SetBinContent(nb, histoPtRap3RecDA[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRap3RecDY_QCD[ntype]->SetBinContent(nb, histoPtRap3RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PtRap4
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPtRap4(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRap4RecDA_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap4RecDA_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap4RecDA_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap4RecDA_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap4RecDA_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap4RecDA_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap4RecDA[ntype]->GetSumOfWeights(),
           histoPtRap4RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap4RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap4RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap4RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap4RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap4RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap4RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPtRap4+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRap4RecDA_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap4RecDA[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRap4RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap4RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRap4RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap4RecDA[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRap4RecDA_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap4RecDA[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRap4RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap4RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRap4RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap4RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRap4RecDA[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRap4RecDA[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRap4RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRap4RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRap4RecDA_QCD[ntype]->SetBinContent(nb, histoPtRap4RecDA[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRap4RecDY_QCD[ntype]->SetBinContent(nb, histoPtRap4RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  char output[200];
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    sprintf(output,"histoDY%dzll_%d.root",whichDY,thePlot);	
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

  {
  sprintf(output,"histoDY%dzllPtRecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  histoPtRecGen[0]->Write();
  histoPtRecGen[1]->Write();
  histoPtRecGen_LepEff[0]->Write();
  histoPtRecGen_LepEff[1]->Write();

  histoPtRecDA[0]->Write();
  histoPtRecDA[1]->Write();
  histoPtRecDY[0]->Write();
  histoPtRecDY[1]->Write();
  histoPtRecDA_MomRes[0]->Write();
  histoPtRecDA_MomRes[1]->Write();
  histoPtRecDY_MomRes[0]->Write();
  histoPtRecDY_MomRes[1]->Write();
  histoPtRecDA_LepEff[0]->Write();
  histoPtRecDA_LepEff[1]->Write();
  histoPtRecDY_LepEff[0]->Write();
  histoPtRecDY_LepEff[1]->Write();
  histoPtRecDA_PDF[0]->Write();
  histoPtRecDA_PDF[1]->Write();
  histoPtRecDY_PDF[0]->Write();
  histoPtRecDY_PDF[1]->Write();
  histoPtRecDA_QCD[0]->Write();
  histoPtRecDA_QCD[1]->Write();
  histoPtRecDY_QCD[0]->Write();
  histoPtRecDY_QCD[1]->Write();
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPt2RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();

  histoPt2RecGen[0]->Write();
  histoPt2RecGen[1]->Write();
  histoPt2RecGen_LepEff[0]->Write();
  histoPt2RecGen_LepEff[1]->Write();

  histoPt2RecDA[0]->Write();
  histoPt2RecDA[1]->Write();
  histoPt2RecDY[0]->Write();
  histoPt2RecDY[1]->Write();
  histoPt2RecDA_MomRes[0]->Write();
  histoPt2RecDA_MomRes[1]->Write();
  histoPt2RecDY_MomRes[0]->Write();
  histoPt2RecDY_MomRes[1]->Write();
  histoPt2RecDA_LepEff[0]->Write();
  histoPt2RecDA_LepEff[1]->Write();
  histoPt2RecDY_LepEff[0]->Write();
  histoPt2RecDY_LepEff[1]->Write();
  histoPt2RecDA_PDF[0]->Write();
  histoPt2RecDA_PDF[1]->Write();
  histoPt2RecDY_PDF[0]->Write();
  histoPt2RecDY_PDF[1]->Write();
  histoPt2RecDA_QCD[0]->Write();
  histoPt2RecDA_QCD[1]->Write();
  histoPt2RecDY_QCD[0]->Write();
  histoPt2RecDY_QCD[1]->Write();
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllRapRecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  histoRapRecGen[0]->Write();
  histoRapRecGen[1]->Write();
  histoRapRecGen_LepEff[0]->Write();
  histoRapRecGen_LepEff[1]->Write();

  histoRapRecDA[0]->Write();
  histoRapRecDA[1]->Write();
  histoRapRecDY[0]->Write();
  histoRapRecDY[1]->Write();
  histoRapRecDA_MomRes[0]->Write();
  histoRapRecDA_MomRes[1]->Write();
  histoRapRecDY_MomRes[0]->Write();
  histoRapRecDY_MomRes[1]->Write();
  histoRapRecDA_LepEff[0]->Write();
  histoRapRecDA_LepEff[1]->Write();
  histoRapRecDY_LepEff[0]->Write();
  histoRapRecDY_LepEff[1]->Write();
  histoRapRecDA_PDF[0]->Write();
  histoRapRecDA_PDF[1]->Write();
  histoRapRecDY_PDF[0]->Write();
  histoRapRecDY_PDF[1]->Write();
  histoRapRecDA_QCD[0]->Write();
  histoRapRecDA_QCD[1]->Write();
  histoRapRecDY_QCD[0]->Write();
  histoRapRecDY_QCD[1]->Write();
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap0RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  histoPtRap0RecGen[0]->Write();
  histoPtRap0RecGen[1]->Write();
  histoPtRap0RecGen_LepEff[0]->Write();
  histoPtRap0RecGen_LepEff[1]->Write();

  histoPtRap0RecDA[0]->Write();
  histoPtRap0RecDA[1]->Write();
  histoPtRap0RecDY[0]->Write();
  histoPtRap0RecDY[1]->Write();
  histoPtRap0RecDA_MomRes[0]->Write();
  histoPtRap0RecDA_MomRes[1]->Write();
  histoPtRap0RecDY_MomRes[0]->Write();
  histoPtRap0RecDY_MomRes[1]->Write();
  histoPtRap0RecDA_LepEff[0]->Write();
  histoPtRap0RecDA_LepEff[1]->Write();
  histoPtRap0RecDY_LepEff[0]->Write();
  histoPtRap0RecDY_LepEff[1]->Write();
  histoPtRap0RecDA_PDF[0]->Write();
  histoPtRap0RecDA_PDF[1]->Write();
  histoPtRap0RecDY_PDF[0]->Write();
  histoPtRap0RecDY_PDF[1]->Write();
  histoPtRap0RecDA_QCD[0]->Write();
  histoPtRap0RecDA_QCD[1]->Write();
  histoPtRap0RecDY_QCD[0]->Write();
  histoPtRap0RecDY_QCD[1]->Write();
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap1RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  histoPtRap1RecGen[0]->Write();
  histoPtRap1RecGen[1]->Write();
  histoPtRap1RecGen_LepEff[0]->Write();
  histoPtRap1RecGen_LepEff[1]->Write();

  histoPtRap1RecDA[0]->Write();
  histoPtRap1RecDA[1]->Write();
  histoPtRap1RecDY[0]->Write();
  histoPtRap1RecDY[1]->Write();
  histoPtRap1RecDA_MomRes[0]->Write();
  histoPtRap1RecDA_MomRes[1]->Write();
  histoPtRap1RecDY_MomRes[0]->Write();
  histoPtRap1RecDY_MomRes[1]->Write();
  histoPtRap1RecDA_LepEff[0]->Write();
  histoPtRap1RecDA_LepEff[1]->Write();
  histoPtRap1RecDY_LepEff[0]->Write();
  histoPtRap1RecDY_LepEff[1]->Write();
  histoPtRap1RecDA_PDF[0]->Write();
  histoPtRap1RecDA_PDF[1]->Write();
  histoPtRap1RecDY_PDF[0]->Write();
  histoPtRap1RecDY_PDF[1]->Write();
  histoPtRap1RecDA_QCD[0]->Write();
  histoPtRap1RecDA_QCD[1]->Write();
  histoPtRap1RecDY_QCD[0]->Write();
  histoPtRap1RecDY_QCD[1]->Write();
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap2RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  histoPtRap2RecGen[0]->Write();
  histoPtRap2RecGen[1]->Write();
  histoPtRap2RecGen_LepEff[0]->Write();
  histoPtRap2RecGen_LepEff[1]->Write();

  histoPtRap2RecDA[0]->Write();
  histoPtRap2RecDA[1]->Write();
  histoPtRap2RecDY[0]->Write();
  histoPtRap2RecDY[1]->Write();
  histoPtRap2RecDA_MomRes[0]->Write();
  histoPtRap2RecDA_MomRes[1]->Write();
  histoPtRap2RecDY_MomRes[0]->Write();
  histoPtRap2RecDY_MomRes[1]->Write();
  histoPtRap2RecDA_LepEff[0]->Write();
  histoPtRap2RecDA_LepEff[1]->Write();
  histoPtRap2RecDY_LepEff[0]->Write();
  histoPtRap2RecDY_LepEff[1]->Write();
  histoPtRap2RecDA_PDF[0]->Write();
  histoPtRap2RecDA_PDF[1]->Write();
  histoPtRap2RecDY_PDF[0]->Write();
  histoPtRap2RecDY_PDF[1]->Write();
  histoPtRap2RecDA_QCD[0]->Write();
  histoPtRap2RecDA_QCD[1]->Write();
  histoPtRap2RecDY_QCD[0]->Write();
  histoPtRap2RecDY_QCD[1]->Write();
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap3RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  histoPtRap3RecGen[0]->Write();
  histoPtRap3RecGen[1]->Write();
  histoPtRap3RecGen_LepEff[0]->Write();
  histoPtRap3RecGen_LepEff[1]->Write();

  histoPtRap3RecDA[0]->Write();
  histoPtRap3RecDA[1]->Write();
  histoPtRap3RecDY[0]->Write();
  histoPtRap3RecDY[1]->Write();
  histoPtRap3RecDA_MomRes[0]->Write();
  histoPtRap3RecDA_MomRes[1]->Write();
  histoPtRap3RecDY_MomRes[0]->Write();
  histoPtRap3RecDY_MomRes[1]->Write();
  histoPtRap3RecDA_LepEff[0]->Write();
  histoPtRap3RecDA_LepEff[1]->Write();
  histoPtRap3RecDY_LepEff[0]->Write();
  histoPtRap3RecDY_LepEff[1]->Write();
  histoPtRap3RecDA_PDF[0]->Write();
  histoPtRap3RecDA_PDF[1]->Write();
  histoPtRap3RecDY_PDF[0]->Write();
  histoPtRap3RecDY_PDF[1]->Write();
  histoPtRap3RecDA_QCD[0]->Write();
  histoPtRap3RecDA_QCD[1]->Write();
  histoPtRap3RecDY_QCD[0]->Write();
  histoPtRap3RecDY_QCD[1]->Write();
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap4RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  histoPtRap4RecGen[0]->Write();
  histoPtRap4RecGen[1]->Write();
  histoPtRap4RecGen_LepEff[0]->Write();
  histoPtRap4RecGen_LepEff[1]->Write();

  histoPtRap4RecDA[0]->Write();
  histoPtRap4RecDA[1]->Write();
  histoPtRap4RecDY[0]->Write();
  histoPtRap4RecDY[1]->Write();
  histoPtRap4RecDA_MomRes[0]->Write();
  histoPtRap4RecDA_MomRes[1]->Write();
  histoPtRap4RecDY_MomRes[0]->Write();
  histoPtRap4RecDY_MomRes[1]->Write();
  histoPtRap4RecDA_LepEff[0]->Write();
  histoPtRap4RecDA_LepEff[1]->Write();
  histoPtRap4RecDY_LepEff[0]->Write();
  histoPtRap4RecDY_LepEff[1]->Write();
  histoPtRap4RecDA_PDF[0]->Write();
  histoPtRap4RecDA_PDF[1]->Write();
  histoPtRap4RecDY_PDF[0]->Write();
  histoPtRap4RecDY_PDF[1]->Write();
  histoPtRap4RecDA_QCD[0]->Write();
  histoPtRap4RecDA_QCD[1]->Write();
  histoPtRap4RecDY_QCD[0]->Write();
  histoPtRap4RecDY_QCD[1]->Write();
  outFilePlots->Close();
  }

}
