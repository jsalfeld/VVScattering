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
#include "MitAnalysisRunII/macros/80x/helicity.h"
#include "MitAnalysisRunII/panda/macros/80x/auxiliar.h"

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
 kTight   =(1<<3),
 kDxyz    =(1<<4)
};

const double mass_el = 0.000510998928;
const double mass_mu = 0.10566;
void pandaAnalysis(int whichDY = 0, int whichAnaFlow = 0, bool isMIT=true)
{

  TString isNoDYName = "";
  if(whichAnaFlow == 1) isNoDYName = "_nody";

  TString dirPathRM = TString(gSystem->Getenv("CMSSW_BASE")) + "/src/MitAnalysisRunII/data/80x/rcdata.2016.v3";
  RoccoR rmcor(dirPathRM.Data());
  double lumi = 35.8;
  double k_eff = 0.5 * sqrt(20374493./12953378.);
  //*******************************************************
  //Input Files
  //*******************************************************
  TString filesPath    = "/data/t3home000/ceballos/panda/v_005_0/";
  if(isMIT == false) filesPath = "/afs/cern.ch/work/c/ceballos/public/samples/panda/v_005_0/";
  vector<TString> infileName_;
  vector<Int_t> infileCat_;
  
  if(whichAnaFlow == 0 || whichAnaFlow == 1){
    infileName_.push_back(Form("%sdata.root",filesPath.Data()));		 infileCat_.push_back(0);
    //infileName_.push_back(Form("%sqqWW.root" ,filesPath.Data()));		   infileCat_.push_back(1);
    //infileName_.push_back(Form("%sggWW.root" ,filesPath.Data()));		   infileCat_.push_back(1);
    infileName_.push_back(Form("%sDYJetsToLL_M-10to50.root" ,filesPath.Data())); infileCat_.push_back(2);
    //infileName_.push_back(Form("%sTT2L.root" ,filesPath.Data()));		   infileCat_.push_back(3);
    //infileName_.push_back(Form("%sTW.root" ,filesPath.Data()));		   infileCat_.push_back(3);

    infileName_.push_back(Form("%sqqZZ.root" ,filesPath.Data()));		 infileCat_.push_back(4);
    infileName_.push_back(Form("%sggZZ.root" ,filesPath.Data()));		 infileCat_.push_back(4);
    infileName_.push_back(Form("%sWZ.root" ,filesPath.Data())); 		 infileCat_.push_back(4);
    infileName_.push_back(Form("%sVVV.root" ,filesPath.Data()));		 infileCat_.push_back(4);
    infileName_.push_back(Form("%sTTV.root" ,filesPath.Data()));		 infileCat_.push_back(4);
    //infileName_.push_back(Form("%sWGstar.root" ,filesPath.Data()));		   infileCat_.push_back(4);
    //infileName_.push_back(Form("%sVG.root" ,filesPath.Data()));		   infileCat_.push_back(6);
    //infileName_.push_back(Form("%sH125.root" ,filesPath.Data()));		   infileCat_.push_back(7);
  }

  if(whichAnaFlow == 0 || whichAnaFlow == 2){
    if     (whichDY == 0)
   {infileName_.push_back(Form("%sDYJetsToLL_M-50_LO_Pt000To050.root" ,filesPath.Data()));  infileCat_.push_back(2);
    infileName_.push_back(Form("%sDYJetsToLL_M-50_LO_Pt100to200.root" ,filesPath.Data()));  infileCat_.push_back(2);
    infileName_.push_back(Form("%sDYJetsToLL_M-50_LO_Pt200toInf.root" ,filesPath.Data()));  infileCat_.push_back(2);
    }
    else if(whichDY == 1)
   {infileName_.push_back(Form("%sDYJetsToLL_M-50_NLO.root",filesPath.Data()));  infileCat_.push_back(2);}
    else if(whichDY == 2)
   {infileName_.push_back(Form("%sDYJetsToLL_POWHEG.root",filesPath.Data()));    infileCat_.push_back(2);}
    else if(whichDY == 3)
   {infileName_.push_back(Form("%sDYJetsToLL_Pt0To50.root",filesPath.Data()));   infileCat_.push_back(2);
    infileName_.push_back(Form("%sDYJetsToLL_Pt50To100.root",filesPath.Data())); infileCat_.push_back(2);
    infileName_.push_back(Form("%sDYJetsToLL_Pt100To250.root",filesPath.Data()));infileCat_.push_back(2);
    infileName_.push_back(Form("%sDYJetsToLL_Pt250To400.root",filesPath.Data()));infileCat_.push_back(2);
    infileName_.push_back(Form("%sDYJetsToLL_Pt400To650.root",filesPath.Data()));infileCat_.push_back(2);
    infileName_.push_back(Form("%sDYJetsToLL_Pt650ToInf.root",filesPath.Data()));infileCat_.push_back(2);
    }
  }

  const double dileptonPtCut = 0.0;
  const int nBinTot = 1; Float_t xbinsTot[nBinTot+1] = {0,1};
  const int nBinPt = 36; Float_t xbinsPt[nBinPt+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,500,800,1500};
  //const int nBinPt = 5; Float_t xbinsPt[nBinPt+1] = {250,300,400,500,800,1500};
  const int nBinRap = 12; Float_t xbinsRap[nBinRap+1] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4};
  const int nBinPhiStar = 34; Float_t xbinsPhiStar[nBinPhiStar+1] = {
                                         1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3,
                                         1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2,
                                         1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1,
                                            1,    3,    5,    7,   10,   20,   30,   50};
  const int nBinPtRap0 = 34; Float_t xbinsPtRap0[nBinPtRap0+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,1500};
  const int nBinPtRap1 = 34; Float_t xbinsPtRap1[nBinPtRap1+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,1500};
  const int nBinPtRap2 = 34; Float_t xbinsPtRap2[nBinPtRap2+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,1500};
  const int nBinPtRap3 = 34; Float_t xbinsPtRap3[nBinPtRap3+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,1500};
  const int nBinPtRap4 = 34; Float_t xbinsPtRap4[nBinPtRap4+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,1500};

  const int nBinRecoTot = 1; Float_t xbinsRecoTot[nBinRecoTot+1] = {0,1};
  const int nBinRecoPt     = 72; Float_t xbinsRecoPt[nBinRecoPt+1]	 = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
  									    10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
  									    28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
  									   220, 235,250, 275,300, 350,400,450,500,650,800,1150,1500};
  //const int nBinRecoPt     = 10; Float_t xbinsRecoPt[nBinRecoPt+1]         = {250,275,300,350,400,450,500,650,800,1150,1500};

  const int nBinRecoRap = 24; Float_t xbinsRecoRap[nBinRecoRap+1] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4};

  const int nBinRecoPhiStar = 61; Float_t xbinsRecoPhiStar[nBinRecoPhiStar+1] = {1e-4,
                                         1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3,
                                         1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2,
                                         1e-1, 1.5e-1, 2e-1, 2.5e-1, 3e-1, 3.5e-1, 4e-1, 4.5e-1,
                                         5e-1, 5.5e-1, 6e-1, 6.5e-1, 7e-1, 7.5e-1, 8e-1, 8.5e-1,
                                         9e-1, 9.5e-1,
                                         1,    1.5,    2,    2.5,    3,    3.5,    4,	 4.5,
                                         5,    5.5,    6,    6.5,    7,    7.5,    8,	 9,
                                         10,    15,   20,     25,   30,     35,   40,   45, 50};

  const int nBinRecoPtRap0 = 70; Float_t xbinsRecoPtRap0[nBinRecoPtRap0+1] = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300, 350,400,450,500,750,1500};

  const int nBinRecoPtRap1 = 70; Float_t xbinsRecoPtRap1[nBinRecoPtRap1+1] = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300, 350,400,450,500,750,1500};

  const int nBinRecoPtRap2 = 70; Float_t xbinsRecoPtRap2[nBinRecoPtRap2+1] = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300, 350,400,450,500,750,1500};

  const int nBinRecoPtRap3 = 70; Float_t xbinsRecoPtRap3[nBinRecoPtRap3+1] = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300, 350,400,450,500,750,1500};

  const int nBinRecoPtRap4 = 70; Float_t xbinsRecoPtRap4[nBinRecoPtRap4+1] = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300, 350,400,450,500,750,1500};

  TFile *fLepton_Eta_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_eta_sf_37ifb.root"));
  TH1D* scalefactors_Muon_Eta = (TH1D*)fLepton_Eta_SF->Get("scalefactors_Muon_Eta"); scalefactors_Muon_Eta->SetDirectory(0);
  TH1D* scalefactors_Electron_Eta = (TH1D*)fLepton_Eta_SF->Get("scalefactors_Electron_Eta"); scalefactors_Electron_Eta->SetDirectory(0);
  fLepton_Eta_SF->Close();

  TFile *fLepton_SF_mu_central = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_dylan_MediumIdOnly.root"));
  TH2D* scalefactors_Medium_Muon = (TH2D*)fLepton_SF_mu_central->Get("scalefactors_Medium_Muon"); scalefactors_Medium_Muon->SetDirectory(0);
  fLepton_SF_mu_central->Close();

  TFile *fLepton_SF_el_central = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_dylan_MediumIdOnly.root"));
  TH2D* scalefactors_Medium_Electron = (TH2D*)fLepton_SF_el_central->Get("scalefactors_Medium_Electron"); scalefactors_Medium_Electron->SetDirectory(0);
  fLepton_SF_el_central->Close();

  TFile *fLepton_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_dylan_MediumIdOnly.root"));

  TH2D* scalefactors_Medium_Muon_stat_error_hi      = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_stat_error_hi");      scalefactors_Medium_Muon_stat_error_hi     ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_signalFsrTNP	    = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_signalFsrTNP");       scalefactors_Medium_Muon_signalFsrTNP      ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_signalResTNP	    = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_signalResTNP");       scalefactors_Medium_Muon_signalResTNP      ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_bkgModelTNP	    = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_bkgModelTNP");        scalefactors_Medium_Muon_bkgModelTNP       ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_tagBiasTNP	    = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_tagBiasTNP");         scalefactors_Medium_Muon_tagBiasTNP        ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_generatorChoiceTNP = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_generatorChoiceTNP"); scalefactors_Medium_Muon_generatorChoiceTNP->SetDirectory(0);
  //const int nMuSFBins = 0;
  //TH2D* scalefactors_Medium_Muon_stat_error_hi_bins[nMuSFBins];
  //for(int nj=0; nj<nMuSFBins; nj++) {scalefactors_Medium_Muon_stat_error_hi_bins[nj] =(TH2D*)fLepton_SF->Get(Form("scalefactors_Medium_Muon_stat_error_hi_bins_%d",nj));scalefactors_Medium_Muon_stat_error_hi_bins[nj]->SetDirectory(0);}

  TH2D* scalefactors_Medium_Electron_stat_error_hi      = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_stat_error_hi");      scalefactors_Medium_Electron_stat_error_hi     ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_signalFsrTNP       = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_signalFsrTNP");       scalefactors_Medium_Electron_signalFsrTNP      ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_signalResTNP       = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_signalResTNP");       scalefactors_Medium_Electron_signalResTNP      ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_bkgModelTNP        = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_bkgModelTNP");        scalefactors_Medium_Electron_bkgModelTNP       ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_tagBiasTNP         = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_tagBiasTNP");         scalefactors_Medium_Electron_tagBiasTNP        ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_generatorChoiceTNP = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_generatorChoiceTNP"); scalefactors_Medium_Electron_generatorChoiceTNP->SetDirectory(0);
  //const int nElSFBins = 0;
  //TH2D* scalefactors_Medium_Electron_stat_error_hi_bins[nElSFBins];
  //for(int nj=0; nj<nElSFBins; nj++) {scalefactors_Medium_Electron_stat_error_hi_bins[nj] = (TH2D*)fLepton_SF->Get(Form("scalefactors_Medium_Electron_stat_error_hi_bins_%d",nj)); scalefactors_Medium_Electron_stat_error_hi_bins[nj]->SetDirectory(0);}

  fLepton_SF->Close();

  double getMaxPtForSFs[4] = {scalefactors_Medium_Muon                  ->GetYaxis()->GetBinCenter(scalefactors_Medium_Muon		     ->GetNbinsY()),
                              scalefactors_Medium_Electron              ->GetYaxis()->GetBinCenter(scalefactors_Medium_Electron		     ->GetNbinsY()),
                              scalefactors_Medium_Muon_stat_error_hi    ->GetYaxis()->GetBinCenter(scalefactors_Medium_Muon_stat_error_hi    ->GetNbinsY()),
		              scalefactors_Medium_Electron_stat_error_hi->GetYaxis()->GetBinCenter(scalefactors_Medium_Electron_stat_error_hi->GetNbinsY())
		              };

  printf("getMaxPtForSFs mu central: %f, el central: %f, mu syst: %f, el syst: %f\n",getMaxPtForSFs[0],getMaxPtForSFs[1],getMaxPtForSFs[2],getMaxPtForSFs[3]);

  const int nBinEtaPlot = 26; Float_t xbinsEtaPlot[nBinEtaPlot+1] = {-2.4,-2.3,-2.2,-2.0,-1.8,-1.63,-1.566,-1.4442,-1.2,-1.0,-0.6,-0.4,-0.2,0.0,
                                                                     0.2,0.4,0.6,1.0,1.2,1.4442,1.566,1.63,1.8,2.0,2.2,2.3,2.4};

  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 62;
  const int histBins = 9;
  TH1D* histo[allPlots][histBins];
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  1) {nBinPlot = 120; xminPlot = 91.1876-15; xmaxPlot = 91.1876+15;}
    else if(thePlot >=  2 && thePlot <=  3) {nBinPlot = 500; xminPlot =  0.0; xmaxPlot =500;}
    else if(thePlot >=  4 && thePlot <=  5) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot =  20;}
    else if(thePlot >=  6 && thePlot <= 11) {nBinPlot = 120; xminPlot = 91.1876-15; xmaxPlot = 91.1876+15;}
    else if(thePlot >= 12 && thePlot <= 15) {nBinPlot = 100; xminPlot =-50.0; xmaxPlot = 50;}
    else if(thePlot >= 16 && thePlot <= 17) {nBinPlot =  96; xminPlot = -2.4; xmaxPlot = 2.4;}
    else if(thePlot >= 18 && thePlot <= 19) {nBinPlot = 100; xminPlot = -1.0; xmaxPlot = 1.0;}
    else if(thePlot >= 22 && thePlot <= 59) {nBinPlot =  96; xminPlot = -2.4; xmaxPlot = 2.4;}
    else if(thePlot >= 60 && thePlot <= 61) {nBinPlot =  48; xminPlot =  0.0; xmaxPlot = 2.4;}
    TH1D* histos;
    if     (thePlot >= 24 && thePlot <= 29){
      histos = new TH1D("histos", "histos", nBinEtaPlot, xbinsEtaPlot);
    }    
    else if(thePlot != 20 && thePlot != 21){
      histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    } 
    else {
      histos = new TH1D("histos", "histos", nBinPhiStar, xbinsPhiStar);
    }
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Reset();histos->Clear();
  }

  const int nRecNuisances = 1;
  const int nEffNuisances = 7;
  const int nMomNuisances = 5;

  TH2D* histoTotRecGen[2];
  TH2D* histoTotRecGen_RecEff[2][nRecNuisances];
  TH2D* histoTotRecGen_LepEff[2][nEffNuisances];
  TH2D* histoTotRecGen_MomRes[2][nMomNuisances];
  TH1D* histoTotRecDA[2];
  TH1D* histoTotRecDY[2];
  TH1D* histoTotRecEM[2];
  TH1D* histoTotRecVV[2];
  TH1D* histoTotRecDY_RecEff[2][nRecNuisances];
  TH1D* histoTotRecVV_RecEff[2][nEffNuisances];
  TH1D* histoTotRecDY_LepEff[2][nEffNuisances];
  TH1D* histoTotRecVV_LepEff[2][nEffNuisances];
  TH1D* histoTotRecDA_MomRes[2][nMomNuisances];
  TH1D* histoTotRecDY_MomRes[2][nMomNuisances];
  TH1D* histoTotRecEM_MomRes[2][nMomNuisances];
  TH1D* histoTotRecVV_MomRes[2][nMomNuisances];
  TH1D* histoTotRecVV_PDF[2];
  TH1D* histoTotRecDY_PDF[2];
  TH1D* histoTotRecVV_QCD[2];
  TH1D* histoTotRecDY_QCD[2];
  TH1D* histoTotRecVV_QCDPart[2][6];
  TH1D* histoTotRecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoTotRecGen[i]        = new TH2D(Form("histoTotRecGen_%d",i),        Form("histoTotRecGen_%d",i),        nBinRecoTot, xbinsRecoTot, nBinTot, xbinsTot);
   for(int j=0; j<nRecNuisances; j++){
     histoTotRecGen_RecEff[i][j] = new TH2D(Form("histoTotRecGen_RecEff_%d_%d",i,j), Form("histoTotRecGen_RecEff_%d_%d",i,j), nBinRecoTot, xbinsRecoTot, nBinTot, xbinsTot);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoTotRecGen_LepEff[i][j] = new TH2D(Form("histoTotRecGen_LepEff_%d_%d",i,j), Form("histoTotRecGen_LepEff_%d_%d",i,j), nBinRecoTot, xbinsRecoTot, nBinTot, xbinsTot);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoTotRecGen_MomRes[i][j] = new TH2D(Form("histoTotRecGen_MomRes_%d_%d",i,j), Form("histoTotRecGen_MomRes_%d_%d",i,j), nBinRecoTot, xbinsRecoTot, nBinTot, xbinsTot);
   }
   histoTotRecDA[i]  = new TH1D(Form("histoTotRecDA_%d",i),  Form("histoTotRecDA_%d",i),  nBinRecoTot, xbinsRecoTot);
   histoTotRecDY[i]  = new TH1D(Form("histoTotRecDY_%d",i),  Form("histoTotRecDY_%d",i),  nBinRecoTot, xbinsRecoTot);
   histoTotRecEM[i]  = new TH1D(Form("histoTotRecEM_%d",i),  Form("histoTotRecEM_%d",i),  nBinRecoTot, xbinsRecoTot);
   histoTotRecVV[i]  = new TH1D(Form("histoTotRecVV_%d",i),  Form("histoTotRecVV_%d",i),  nBinRecoTot, xbinsRecoTot);

   for(int j=0; j<nRecNuisances; j++){
     histoTotRecDY_RecEff[i][j] = new TH1D(Form("histoTotRecDY_RecEff_%d_%d",i,j), Form("histoTotRecDY_RecEff_%d_%d",i,j), nBinRecoTot, xbinsRecoTot);
     histoTotRecVV_RecEff[i][j] = new TH1D(Form("histoTotRecVV_RecEff_%d_%d",i,j), Form("histoTotRecVV_RecEff_%d_%d",i,j), nBinRecoTot, xbinsRecoTot);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoTotRecDY_LepEff[i][j] = new TH1D(Form("histoTotRecDY_LepEff_%d_%d",i,j), Form("histoTotRecDY_LepEff_%d_%d",i,j), nBinRecoTot, xbinsRecoTot);
     histoTotRecVV_LepEff[i][j] = new TH1D(Form("histoTotRecVV_LepEff_%d_%d",i,j), Form("histoTotRecVV_LepEff_%d_%d",i,j), nBinRecoTot, xbinsRecoTot);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoTotRecDA_MomRes[i][j] = new TH1D(Form("histoTotRecDA_MomRes_%d_%d",i,j), Form("histoTotRecDA_MomRes_%d_%d",i,j), nBinRecoTot, xbinsRecoTot);
     histoTotRecDY_MomRes[i][j] = new TH1D(Form("histoTotRecDY_MomRes_%d_%d",i,j), Form("histoTotRecDY_MomRes_%d_%d",i,j), nBinRecoTot, xbinsRecoTot);
     histoTotRecEM_MomRes[i][j] = new TH1D(Form("histoTotRecEM_MomRes_%d_%d",i,j), Form("histoTotRecEM_MomRes_%d_%d",i,j), nBinRecoTot, xbinsRecoTot);
     histoTotRecVV_MomRes[i][j] = new TH1D(Form("histoTotRecVV_MomRes_%d_%d",i,j), Form("histoTotRecVV_MomRes_%d_%d",i,j), nBinRecoTot, xbinsRecoTot);
   }
   histoTotRecVV_PDF[i] = new TH1D(Form("histoTotRecVV_PDF_%d",i), Form("histoTotRecVV_PDF_%d",i), nBinRecoTot, xbinsRecoTot);
   histoTotRecDY_PDF[i] = new TH1D(Form("histoTotRecDY_PDF_%d",i), Form("histoTotRecDY_PDF_%d",i), nBinRecoTot, xbinsRecoTot);

   histoTotRecVV_QCD[i] = new TH1D(Form("histoTotRecVV_QCD_%d",i), Form("histoTotRecVV_QCD_%d",i), nBinRecoTot, xbinsRecoTot);
   histoTotRecDY_QCD[i] = new TH1D(Form("histoTotRecDY_QCD_%d",i), Form("histoTotRecDY_QCD_%d",i), nBinRecoTot, xbinsRecoTot);

   for(int j=0; j<6; j++){
      histoTotRecVV_QCDPart[i][j] = new TH1D(Form("histoTotRecVV_QCDPart_%d_%d",i,j), Form("histoTotRecVV_QCDPart_%d_%d",i,j), nBinRecoTot, xbinsRecoTot);
      histoTotRecDY_QCDPart[i][j] = new TH1D(Form("histoTotRecDY_QCDPart_%d_%d",i,j), Form("histoTotRecDY_QCDPart_%d_%d",i,j), nBinRecoTot, xbinsRecoTot);
    }
  }

  TH2D* histoPtRecGen[2];
  TH2D* histoPtRecGen_RecEff[2][nRecNuisances];
  TH2D* histoPtRecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRecDA[2];
  TH1D* histoPtRecDY[2];
  TH1D* histoPtRecEM[2];
  TH1D* histoPtRecVV[2];
  TH1D* histoPtRecDY_RecEff[2][nRecNuisances];
  TH1D* histoPtRecVV_RecEff[2][nEffNuisances];
  TH1D* histoPtRecDY_LepEff[2][nEffNuisances];
  TH1D* histoPtRecVV_LepEff[2][nEffNuisances];
  TH1D* histoPtRecDA_MomRes[2][nMomNuisances];
  TH1D* histoPtRecDY_MomRes[2][nMomNuisances];
  TH1D* histoPtRecEM_MomRes[2][nMomNuisances];
  TH1D* histoPtRecVV_MomRes[2][nMomNuisances];
  TH1D* histoPtRecVV_PDF[2];
  TH1D* histoPtRecDY_PDF[2];
  TH1D* histoPtRecVV_QCD[2];
  TH1D* histoPtRecDY_QCD[2];
  TH1D* histoPtRecVV_QCDPart[2][6];
  TH1D* histoPtRecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRecGen[i]        = new TH2D(Form("histoPtRecGen_%d",i),        Form("histoPtRecGen_%d",i),        nBinRecoPt, xbinsRecoPt, nBinPt, xbinsPt);
   for(int j=0; j<nRecNuisances; j++){
     histoPtRecGen_RecEff[i][j] = new TH2D(Form("histoPtRecGen_RecEff_%d_%d",i,j), Form("histoPtRecGen_RecEff_%d_%d",i,j), nBinRecoPt, xbinsRecoPt, nBinPt, xbinsPt);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRecGen_LepEff[i][j] = new TH2D(Form("histoPtRecGen_LepEff_%d_%d",i,j), Form("histoPtRecGen_LepEff_%d_%d",i,j), nBinRecoPt, xbinsRecoPt, nBinPt, xbinsPt);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRecGen_MomRes[i][j] = new TH2D(Form("histoPtRecGen_MomRes_%d_%d",i,j), Form("histoPtRecGen_MomRes_%d_%d",i,j), nBinRecoPt, xbinsRecoPt, nBinPt, xbinsPt);
   }
   histoPtRecDA[i]  = new TH1D(Form("histoPtRecDA_%d",i),  Form("histoPtRecDA_%d",i),  nBinRecoPt, xbinsRecoPt);
   histoPtRecDY[i]  = new TH1D(Form("histoPtRecDY_%d",i),  Form("histoPtRecDY_%d",i),  nBinRecoPt, xbinsRecoPt);
   histoPtRecEM[i]  = new TH1D(Form("histoPtRecEM_%d",i),  Form("histoPtRecEM_%d",i),  nBinRecoPt, xbinsRecoPt);
   histoPtRecVV[i]  = new TH1D(Form("histoPtRecVV_%d",i),  Form("histoPtRecVV_%d",i),  nBinRecoPt, xbinsRecoPt);

   for(int j=0; j<nRecNuisances; j++){
     histoPtRecDY_RecEff[i][j] = new TH1D(Form("histoPtRecDY_RecEff_%d_%d",i,j), Form("histoPtRecDY_RecEff_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
     histoPtRecVV_RecEff[i][j] = new TH1D(Form("histoPtRecVV_RecEff_%d_%d",i,j), Form("histoPtRecVV_RecEff_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRecDY_LepEff[i][j] = new TH1D(Form("histoPtRecDY_LepEff_%d_%d",i,j), Form("histoPtRecDY_LepEff_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
     histoPtRecVV_LepEff[i][j] = new TH1D(Form("histoPtRecVV_LepEff_%d_%d",i,j), Form("histoPtRecVV_LepEff_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRecDA_MomRes[i][j] = new TH1D(Form("histoPtRecDA_MomRes_%d_%d",i,j), Form("histoPtRecDA_MomRes_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
     histoPtRecDY_MomRes[i][j] = new TH1D(Form("histoPtRecDY_MomRes_%d_%d",i,j), Form("histoPtRecDY_MomRes_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
     histoPtRecEM_MomRes[i][j] = new TH1D(Form("histoPtRecEM_MomRes_%d_%d",i,j), Form("histoPtRecEM_MomRes_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
     histoPtRecVV_MomRes[i][j] = new TH1D(Form("histoPtRecVV_MomRes_%d_%d",i,j), Form("histoPtRecVV_MomRes_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
   }
   histoPtRecVV_PDF[i] = new TH1D(Form("histoPtRecVV_PDF_%d",i), Form("histoPtRecVV_PDF_%d",i), nBinRecoPt, xbinsRecoPt);
   histoPtRecDY_PDF[i] = new TH1D(Form("histoPtRecDY_PDF_%d",i), Form("histoPtRecDY_PDF_%d",i), nBinRecoPt, xbinsRecoPt);

   histoPtRecVV_QCD[i] = new TH1D(Form("histoPtRecVV_QCD_%d",i), Form("histoPtRecVV_QCD_%d",i), nBinRecoPt, xbinsRecoPt);
   histoPtRecDY_QCD[i] = new TH1D(Form("histoPtRecDY_QCD_%d",i), Form("histoPtRecDY_QCD_%d",i), nBinRecoPt, xbinsRecoPt);

   for(int j=0; j<6; j++){
      histoPtRecVV_QCDPart[i][j] = new TH1D(Form("histoPtRecVV_QCDPart_%d_%d",i,j), Form("histoPtRecVV_QCDPart_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
      histoPtRecDY_QCDPart[i][j] = new TH1D(Form("histoPtRecDY_QCDPart_%d_%d",i,j), Form("histoPtRecDY_QCDPart_%d_%d",i,j), nBinRecoPt, xbinsRecoPt);
    }
  }

  TH2D* histoRapRecGen[2];
  TH2D* histoRapRecGen_RecEff[2][nRecNuisances];
  TH2D* histoRapRecGen_LepEff[2][nEffNuisances];
  TH2D* histoRapRecGen_MomRes[2][nMomNuisances];
  TH1D* histoRapRecDA[2];
  TH1D* histoRapRecDY[2];
  TH1D* histoRapRecEM[2];
  TH1D* histoRapRecVV[2];
  TH1D* histoRapRecDY_RecEff[2][nRecNuisances];
  TH1D* histoRapRecVV_RecEff[2][nEffNuisances];
  TH1D* histoRapRecDY_LepEff[2][nEffNuisances];
  TH1D* histoRapRecVV_LepEff[2][nEffNuisances];
  TH1D* histoRapRecDA_MomRes[2][nMomNuisances];
  TH1D* histoRapRecDY_MomRes[2][nMomNuisances];
  TH1D* histoRapRecEM_MomRes[2][nMomNuisances];
  TH1D* histoRapRecVV_MomRes[2][nMomNuisances];
  TH1D* histoRapRecVV_PDF[2];
  TH1D* histoRapRecDY_PDF[2];
  TH1D* histoRapRecVV_QCD[2];
  TH1D* histoRapRecDY_QCD[2];
  TH1D* histoRapRecVV_QCDPart[2][6];
  TH1D* histoRapRecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoRapRecGen[i]        = new TH2D(Form("histoRapRecGen_%d",i),        Form("histoRapRecGen_%d",i),        nBinRecoRap, xbinsRecoRap, nBinRap, xbinsRap);
   for(int j=0; j<nRecNuisances; j++){
     histoRapRecGen_RecEff[i][j] = new TH2D(Form("histoRapRecGen_RecEff_%d_%d",i,j), Form("histoRapRecGen_RecEff_%d_%d",i,j), nBinRecoRap, xbinsRecoRap, nBinRap, xbinsRap);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoRapRecGen_LepEff[i][j] = new TH2D(Form("histoRapRecGen_LepEff_%d_%d",i,j), Form("histoRapRecGen_LepEff_%d_%d",i,j), nBinRecoRap, xbinsRecoRap, nBinRap, xbinsRap);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoRapRecGen_MomRes[i][j] = new TH2D(Form("histoRapRecGen_MomRes_%d_%d",i,j), Form("histoRapRecGen_MomRes_%d_%d",i,j), nBinRecoRap, xbinsRecoRap, nBinRap, xbinsRap);
   }
   histoRapRecDA[i]  = new TH1D(Form("histoRapRecDA_%d",i),  Form("histoRapRecDA_%d",i),  nBinRecoRap, xbinsRecoRap);
   histoRapRecDY[i]  = new TH1D(Form("histoRapRecDY_%d",i),  Form("histoRapRecDY_%d",i),  nBinRecoRap, xbinsRecoRap);
   histoRapRecEM[i]  = new TH1D(Form("histoRapRecEM_%d",i),  Form("histoRapRecEM_%d",i),  nBinRecoRap, xbinsRecoRap);
   histoRapRecVV[i]  = new TH1D(Form("histoRapRecVV_%d",i),  Form("histoRapRecVV_%d",i),  nBinRecoRap, xbinsRecoRap);

   for(int j=0; j<nRecNuisances; j++){
     histoRapRecDY_RecEff[i][j] = new TH1D(Form("histoRapRecDY_RecEff_%d_%d",i,j), Form("histoRapRecDY_RecEff_%d_%d",i,j), nBinRecoRap, xbinsRecoRap);
     histoRapRecVV_RecEff[i][j] = new TH1D(Form("histoRapRecVV_RecEff_%d_%d",i,j), Form("histoRapRecVV_RecEff_%d_%d",i,j), nBinRecoRap, xbinsRecoRap);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoRapRecDY_LepEff[i][j] = new TH1D(Form("histoRapRecDY_LepEff_%d_%d",i,j), Form("histoRapRecDY_LepEff_%d_%d",i,j), nBinRecoRap, xbinsRecoRap);
     histoRapRecVV_LepEff[i][j] = new TH1D(Form("histoRapRecVV_LepEff_%d_%d",i,j), Form("histoRapRecVV_LepEff_%d_%d",i,j), nBinRecoRap, xbinsRecoRap);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoRapRecDA_MomRes[i][j] = new TH1D(Form("histoRapRecDA_MomRes_%d_%d",i,j), Form("histoRapRecDA_MomRes_%d_%d",i,j), nBinRecoRap, xbinsRecoRap);
     histoRapRecDY_MomRes[i][j] = new TH1D(Form("histoRapRecDY_MomRes_%d_%d",i,j), Form("histoRapRecDY_MomRes_%d_%d",i,j), nBinRecoRap, xbinsRecoRap);
     histoRapRecEM_MomRes[i][j] = new TH1D(Form("histoRapRecEM_MomRes_%d_%d",i,j), Form("histoRapRecEM_MomRes_%d_%d",i,j), nBinRecoRap, xbinsRecoRap);
     histoRapRecVV_MomRes[i][j] = new TH1D(Form("histoRapRecVV_MomRes_%d_%d",i,j), Form("histoRapRecVV_MomRes_%d_%d",i,j), nBinRecoRap, xbinsRecoRap);
   }
   histoRapRecVV_PDF[i] = new TH1D(Form("histoRapRecVV_PDF_%d",i), Form("histoRapRecVV_PDF_%d",i), nBinRecoRap, xbinsRecoRap);
   histoRapRecDY_PDF[i] = new TH1D(Form("histoRapRecDY_PDF_%d",i), Form("histoRapRecDY_PDF_%d",i), nBinRecoRap, xbinsRecoRap);

   histoRapRecVV_QCD[i] = new TH1D(Form("histoRapRecVV_QCD_%d",i), Form("histoRapRecVV_QCD_%d",i), nBinRecoRap, xbinsRecoRap);
   histoRapRecDY_QCD[i] = new TH1D(Form("histoRapRecDY_QCD_%d",i), Form("histoRapRecDY_QCD_%d",i), nBinRecoRap, xbinsRecoRap);

   for(int j=0; j<6; j++){
      histoRapRecVV_QCDPart[i][j] = new TH1D(Form("histoRapRecVV_QCDPart_%d_%d",i,j), Form("histoRapRecVV_QCDPart_%d_%d",i,j), nBinRecoRap, xbinsRecoRap);
      histoRapRecDY_QCDPart[i][j] = new TH1D(Form("histoRapRecDY_QCDPart_%d_%d",i,j), Form("histoRapRecDY_QCDPart_%d_%d",i,j), nBinRecoRap, xbinsRecoRap);
    }
  }

  TH2D* histoPhiStarRecGen[2];
  TH2D* histoPhiStarRecGen_RecEff[2][nRecNuisances];
  TH2D* histoPhiStarRecGen_LepEff[2][nEffNuisances];
  TH2D* histoPhiStarRecGen_MomRes[2][nMomNuisances];
  TH1D* histoPhiStarRecDA[2];
  TH1D* histoPhiStarRecDY[2];
  TH1D* histoPhiStarRecEM[2];
  TH1D* histoPhiStarRecVV[2];
  TH1D* histoPhiStarRecDY_RecEff[2][nRecNuisances];
  TH1D* histoPhiStarRecVV_RecEff[2][nEffNuisances];
  TH1D* histoPhiStarRecDY_LepEff[2][nEffNuisances];
  TH1D* histoPhiStarRecVV_LepEff[2][nEffNuisances];
  TH1D* histoPhiStarRecDA_MomRes[2][nMomNuisances];
  TH1D* histoPhiStarRecDY_MomRes[2][nMomNuisances];
  TH1D* histoPhiStarRecEM_MomRes[2][nMomNuisances];
  TH1D* histoPhiStarRecVV_MomRes[2][nMomNuisances];
  TH1D* histoPhiStarRecVV_PDF[2];
  TH1D* histoPhiStarRecDY_PDF[2];
  TH1D* histoPhiStarRecVV_QCD[2];
  TH1D* histoPhiStarRecDY_QCD[2];
  TH1D* histoPhiStarRecVV_QCDPart[2][6];
  TH1D* histoPhiStarRecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPhiStarRecGen[i]        = new TH2D(Form("histoPhiStarRecGen_%d",i),        Form("histoPhiStarRecGen_%d",i),        nBinRecoPhiStar, xbinsRecoPhiStar, nBinPhiStar, xbinsPhiStar);
   for(int j=0; j<nRecNuisances; j++){
     histoPhiStarRecGen_RecEff[i][j] = new TH2D(Form("histoPhiStarRecGen_RecEff_%d_%d",i,j), Form("histoPhiStarRecGen_RecEff_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar, nBinPhiStar, xbinsPhiStar);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPhiStarRecGen_LepEff[i][j] = new TH2D(Form("histoPhiStarRecGen_LepEff_%d_%d",i,j), Form("histoPhiStarRecGen_LepEff_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar, nBinPhiStar, xbinsPhiStar);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPhiStarRecGen_MomRes[i][j] = new TH2D(Form("histoPhiStarRecGen_MomRes_%d_%d",i,j), Form("histoPhiStarRecGen_MomRes_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar, nBinPhiStar, xbinsPhiStar);
   }
   histoPhiStarRecDA[i]  = new TH1D(Form("histoPhiStarRecDA_%d",i),  Form("histoPhiStarRecDA_%d",i),  nBinRecoPhiStar, xbinsRecoPhiStar);
   histoPhiStarRecDY[i]  = new TH1D(Form("histoPhiStarRecDY_%d",i),  Form("histoPhiStarRecDY_%d",i),  nBinRecoPhiStar, xbinsRecoPhiStar);
   histoPhiStarRecEM[i]  = new TH1D(Form("histoPhiStarRecEM_%d",i),  Form("histoPhiStarRecEM_%d",i),  nBinRecoPhiStar, xbinsRecoPhiStar);
   histoPhiStarRecVV[i]  = new TH1D(Form("histoPhiStarRecVV_%d",i),  Form("histoPhiStarRecVV_%d",i),  nBinRecoPhiStar, xbinsRecoPhiStar);

   for(int j=0; j<nRecNuisances; j++){
     histoPhiStarRecDY_RecEff[i][j] = new TH1D(Form("histoPhiStarRecDY_RecEff_%d_%d",i,j), Form("histoPhiStarRecDY_RecEff_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar);
     histoPhiStarRecVV_RecEff[i][j] = new TH1D(Form("histoPhiStarRecVV_RecEff_%d_%d",i,j), Form("histoPhiStarRecVV_RecEff_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPhiStarRecDY_LepEff[i][j] = new TH1D(Form("histoPhiStarRecDY_LepEff_%d_%d",i,j), Form("histoPhiStarRecDY_LepEff_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar);
     histoPhiStarRecVV_LepEff[i][j] = new TH1D(Form("histoPhiStarRecVV_LepEff_%d_%d",i,j), Form("histoPhiStarRecVV_LepEff_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPhiStarRecDA_MomRes[i][j] = new TH1D(Form("histoPhiStarRecDA_MomRes_%d_%d",i,j), Form("histoPhiStarRecDA_MomRes_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar);
     histoPhiStarRecDY_MomRes[i][j] = new TH1D(Form("histoPhiStarRecDY_MomRes_%d_%d",i,j), Form("histoPhiStarRecDY_MomRes_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar);
     histoPhiStarRecEM_MomRes[i][j] = new TH1D(Form("histoPhiStarRecEM_MomRes_%d_%d",i,j), Form("histoPhiStarRecEM_MomRes_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar);
     histoPhiStarRecVV_MomRes[i][j] = new TH1D(Form("histoPhiStarRecVV_MomRes_%d_%d",i,j), Form("histoPhiStarRecVV_MomRes_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar);
   }
   histoPhiStarRecVV_PDF[i] = new TH1D(Form("histoPhiStarRecVV_PDF_%d",i), Form("histoPhiStarRecVV_PDF_%d",i), nBinRecoPhiStar, xbinsRecoPhiStar);
   histoPhiStarRecDY_PDF[i] = new TH1D(Form("histoPhiStarRecDY_PDF_%d",i), Form("histoPhiStarRecDY_PDF_%d",i), nBinRecoPhiStar, xbinsRecoPhiStar);

   histoPhiStarRecVV_QCD[i] = new TH1D(Form("histoPhiStarRecVV_QCD_%d",i), Form("histoPhiStarRecVV_QCD_%d",i), nBinRecoPhiStar, xbinsRecoPhiStar);
   histoPhiStarRecDY_QCD[i] = new TH1D(Form("histoPhiStarRecDY_QCD_%d",i), Form("histoPhiStarRecDY_QCD_%d",i), nBinRecoPhiStar, xbinsRecoPhiStar);

   for(int j=0; j<6; j++){
      histoPhiStarRecVV_QCDPart[i][j] = new TH1D(Form("histoPhiStarRecVV_QCDPart_%d_%d",i,j), Form("histoPhiStarRecVV_QCDPart_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar);
      histoPhiStarRecDY_QCDPart[i][j] = new TH1D(Form("histoPhiStarRecDY_QCDPart_%d_%d",i,j), Form("histoPhiStarRecDY_QCDPart_%d_%d",i,j), nBinRecoPhiStar, xbinsRecoPhiStar);
    }
  }

  TH2D* histoPtRap0RecGen[2];
  TH2D* histoPtRap0RecGen_RecEff[2][nRecNuisances];
  TH2D* histoPtRap0RecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRap0RecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRap0RecDA[2];
  TH1D* histoPtRap0RecDY[2];
  TH1D* histoPtRap0RecEM[2];
  TH1D* histoPtRap0RecVV[2];
  TH1D* histoPtRap0RecDY_RecEff[2][nRecNuisances];
  TH1D* histoPtRap0RecVV_RecEff[2][nEffNuisances];
  TH1D* histoPtRap0RecDY_LepEff[2][nEffNuisances];
  TH1D* histoPtRap0RecVV_LepEff[2][nEffNuisances];
  TH1D* histoPtRap0RecDA_MomRes[2][nMomNuisances];
  TH1D* histoPtRap0RecDY_MomRes[2][nMomNuisances];
  TH1D* histoPtRap0RecEM_MomRes[2][nMomNuisances];
  TH1D* histoPtRap0RecVV_MomRes[2][nMomNuisances];
  TH1D* histoPtRap0RecVV_PDF[2];
  TH1D* histoPtRap0RecDY_PDF[2];
  TH1D* histoPtRap0RecVV_QCD[2];
  TH1D* histoPtRap0RecDY_QCD[2];
  TH1D* histoPtRap0RecVV_QCDPart[2][6];
  TH1D* histoPtRap0RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRap0RecGen[i]        = new TH2D(Form("histoPtRap0RecGen_%d",i),        Form("histoPtRap0RecGen_%d",i),        nBinRecoPtRap0, xbinsRecoPtRap0, nBinPtRap0, xbinsPtRap0);
   for(int j=0; j<nRecNuisances; j++){
     histoPtRap0RecGen_RecEff[i][j] = new TH2D(Form("histoPtRap0RecGen_RecEff_%d_%d",i,j), Form("histoPtRap0RecGen_RecEff_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0, nBinPtRap0, xbinsPtRap0);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRap0RecGen_LepEff[i][j] = new TH2D(Form("histoPtRap0RecGen_LepEff_%d_%d",i,j), Form("histoPtRap0RecGen_LepEff_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0, nBinPtRap0, xbinsPtRap0);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRap0RecGen_MomRes[i][j] = new TH2D(Form("histoPtRap0RecGen_MomRes_%d_%d",i,j), Form("histoPtRap0RecGen_MomRes_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0, nBinPtRap0, xbinsPtRap0);
   }
   histoPtRap0RecDA[i]  = new TH1D(Form("histoPtRap0RecDA_%d",i),  Form("histoPtRap0RecDA_%d",i),  nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecDY[i]  = new TH1D(Form("histoPtRap0RecDY_%d",i),  Form("histoPtRap0RecDY_%d",i),  nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecEM[i]  = new TH1D(Form("histoPtRap0RecEM_%d",i),  Form("histoPtRap0RecEM_%d",i),  nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecVV[i]  = new TH1D(Form("histoPtRap0RecVV_%d",i),  Form("histoPtRap0RecVV_%d",i),  nBinRecoPtRap0, xbinsRecoPtRap0);

   for(int j=0; j<nRecNuisances; j++){
     histoPtRap0RecDY_RecEff[i][j] = new TH1D(Form("histoPtRap0RecDY_RecEff_%d_%d",i,j), Form("histoPtRap0RecDY_RecEff_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
     histoPtRap0RecVV_RecEff[i][j] = new TH1D(Form("histoPtRap0RecVV_RecEff_%d_%d",i,j), Form("histoPtRap0RecVV_RecEff_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRap0RecDY_LepEff[i][j] = new TH1D(Form("histoPtRap0RecDY_LepEff_%d_%d",i,j), Form("histoPtRap0RecDY_LepEff_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
     histoPtRap0RecVV_LepEff[i][j] = new TH1D(Form("histoPtRap0RecVV_LepEff_%d_%d",i,j), Form("histoPtRap0RecVV_LepEff_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRap0RecDA_MomRes[i][j] = new TH1D(Form("histoPtRap0RecDA_MomRes_%d_%d",i,j), Form("histoPtRap0RecDA_MomRes_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
     histoPtRap0RecDY_MomRes[i][j] = new TH1D(Form("histoPtRap0RecDY_MomRes_%d_%d",i,j), Form("histoPtRap0RecDY_MomRes_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
     histoPtRap0RecEM_MomRes[i][j] = new TH1D(Form("histoPtRap0RecEM_MomRes_%d_%d",i,j), Form("histoPtRap0RecEM_MomRes_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
     histoPtRap0RecVV_MomRes[i][j] = new TH1D(Form("histoPtRap0RecVV_MomRes_%d_%d",i,j), Form("histoPtRap0RecVV_MomRes_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
   }
   histoPtRap0RecVV_PDF[i] = new TH1D(Form("histoPtRap0RecVV_PDF_%d",i), Form("histoPtRap0RecVV_PDF_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecDY_PDF[i] = new TH1D(Form("histoPtRap0RecDY_PDF_%d",i), Form("histoPtRap0RecDY_PDF_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);

   histoPtRap0RecVV_QCD[i] = new TH1D(Form("histoPtRap0RecVV_QCD_%d",i), Form("histoPtRap0RecVV_QCD_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);
   histoPtRap0RecDY_QCD[i] = new TH1D(Form("histoPtRap0RecDY_QCD_%d",i), Form("histoPtRap0RecDY_QCD_%d",i), nBinRecoPtRap0, xbinsRecoPtRap0);

   for(int j=0; j<6; j++){
      histoPtRap0RecVV_QCDPart[i][j] = new TH1D(Form("histoPtRap0RecVV_QCDPart_%d_%d",i,j), Form("histoPtRap0RecVV_QCDPart_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
      histoPtRap0RecDY_QCDPart[i][j] = new TH1D(Form("histoPtRap0RecDY_QCDPart_%d_%d",i,j), Form("histoPtRap0RecDY_QCDPart_%d_%d",i,j), nBinRecoPtRap0, xbinsRecoPtRap0);
    }
  }

  TH2D* histoPtRap1RecGen[2];
  TH2D* histoPtRap1RecGen_RecEff[2][nRecNuisances];
  TH2D* histoPtRap1RecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRap1RecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRap1RecDA[2];
  TH1D* histoPtRap1RecDY[2];
  TH1D* histoPtRap1RecEM[2];
  TH1D* histoPtRap1RecVV[2];
  TH1D* histoPtRap1RecDY_RecEff[2][nRecNuisances];
  TH1D* histoPtRap1RecVV_RecEff[2][nEffNuisances];
  TH1D* histoPtRap1RecDY_LepEff[2][nEffNuisances];
  TH1D* histoPtRap1RecVV_LepEff[2][nEffNuisances];
  TH1D* histoPtRap1RecDA_MomRes[2][nMomNuisances];
  TH1D* histoPtRap1RecDY_MomRes[2][nMomNuisances];
  TH1D* histoPtRap1RecEM_MomRes[2][nMomNuisances];
  TH1D* histoPtRap1RecVV_MomRes[2][nMomNuisances];
  TH1D* histoPtRap1RecVV_PDF[2];
  TH1D* histoPtRap1RecDY_PDF[2];
  TH1D* histoPtRap1RecVV_QCD[2];
  TH1D* histoPtRap1RecDY_QCD[2];
  TH1D* histoPtRap1RecVV_QCDPart[2][6];
  TH1D* histoPtRap1RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRap1RecGen[i]        = new TH2D(Form("histoPtRap1RecGen_%d",i),        Form("histoPtRap1RecGen_%d",i),        nBinRecoPtRap1, xbinsRecoPtRap1, nBinPtRap1, xbinsPtRap1);
   for(int j=0; j<nRecNuisances; j++){
     histoPtRap1RecGen_RecEff[i][j] = new TH2D(Form("histoPtRap1RecGen_RecEff_%d_%d",i,j), Form("histoPtRap1RecGen_RecEff_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1, nBinPtRap1, xbinsPtRap1);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRap1RecGen_LepEff[i][j] = new TH2D(Form("histoPtRap1RecGen_LepEff_%d_%d",i,j), Form("histoPtRap1RecGen_LepEff_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1, nBinPtRap1, xbinsPtRap1);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRap1RecGen_MomRes[i][j] = new TH2D(Form("histoPtRap1RecGen_MomRes_%d_%d",i,j), Form("histoPtRap1RecGen_MomRes_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1, nBinPtRap1, xbinsPtRap1);
   }
   histoPtRap1RecDA[i]  = new TH1D(Form("histoPtRap1RecDA_%d",i),  Form("histoPtRap1RecDA_%d",i),  nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecDY[i]  = new TH1D(Form("histoPtRap1RecDY_%d",i),  Form("histoPtRap1RecDY_%d",i),  nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecEM[i]  = new TH1D(Form("histoPtRap1RecEM_%d",i),  Form("histoPtRap1RecEM_%d",i),  nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecVV[i]  = new TH1D(Form("histoPtRap1RecVV_%d",i),  Form("histoPtRap1RecVV_%d",i),  nBinRecoPtRap1, xbinsRecoPtRap1);

   for(int j=0; j<nRecNuisances; j++){
     histoPtRap1RecDY_RecEff[i][j] = new TH1D(Form("histoPtRap1RecDY_RecEff_%d_%d",i,j), Form("histoPtRap1RecDY_RecEff_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
     histoPtRap1RecVV_RecEff[i][j] = new TH1D(Form("histoPtRap1RecVV_RecEff_%d_%d",i,j), Form("histoPtRap1RecVV_RecEff_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRap1RecDY_LepEff[i][j] = new TH1D(Form("histoPtRap1RecDY_LepEff_%d_%d",i,j), Form("histoPtRap1RecDY_LepEff_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
     histoPtRap1RecVV_LepEff[i][j] = new TH1D(Form("histoPtRap1RecVV_LepEff_%d_%d",i,j), Form("histoPtRap1RecVV_LepEff_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRap1RecDA_MomRes[i][j] = new TH1D(Form("histoPtRap1RecDA_MomRes_%d_%d",i,j), Form("histoPtRap1RecDA_MomRes_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
     histoPtRap1RecDY_MomRes[i][j] = new TH1D(Form("histoPtRap1RecDY_MomRes_%d_%d",i,j), Form("histoPtRap1RecDY_MomRes_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
     histoPtRap1RecEM_MomRes[i][j] = new TH1D(Form("histoPtRap1RecEM_MomRes_%d_%d",i,j), Form("histoPtRap1RecEM_MomRes_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
     histoPtRap1RecVV_MomRes[i][j] = new TH1D(Form("histoPtRap1RecVV_MomRes_%d_%d",i,j), Form("histoPtRap1RecVV_MomRes_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
   }
   histoPtRap1RecVV_PDF[i] = new TH1D(Form("histoPtRap1RecVV_PDF_%d",i), Form("histoPtRap1RecVV_PDF_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecDY_PDF[i] = new TH1D(Form("histoPtRap1RecDY_PDF_%d",i), Form("histoPtRap1RecDY_PDF_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);

   histoPtRap1RecVV_QCD[i] = new TH1D(Form("histoPtRap1RecVV_QCD_%d",i), Form("histoPtRap1RecVV_QCD_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);
   histoPtRap1RecDY_QCD[i] = new TH1D(Form("histoPtRap1RecDY_QCD_%d",i), Form("histoPtRap1RecDY_QCD_%d",i), nBinRecoPtRap1, xbinsRecoPtRap1);

   for(int j=0; j<6; j++){
      histoPtRap1RecVV_QCDPart[i][j] = new TH1D(Form("histoPtRap1RecVV_QCDPart_%d_%d",i,j), Form("histoPtRap1RecVV_QCDPart_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
      histoPtRap1RecDY_QCDPart[i][j] = new TH1D(Form("histoPtRap1RecDY_QCDPart_%d_%d",i,j), Form("histoPtRap1RecDY_QCDPart_%d_%d",i,j), nBinRecoPtRap1, xbinsRecoPtRap1);
    }
  }

  TH2D* histoPtRap2RecGen[2];
  TH2D* histoPtRap2RecGen_RecEff[2][nRecNuisances];
  TH2D* histoPtRap2RecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRap2RecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRap2RecDA[2];
  TH1D* histoPtRap2RecDY[2];
  TH1D* histoPtRap2RecEM[2];
  TH1D* histoPtRap2RecVV[2];
  TH1D* histoPtRap2RecDY_RecEff[2][nRecNuisances];
  TH1D* histoPtRap2RecVV_RecEff[2][nEffNuisances];
  TH1D* histoPtRap2RecDY_LepEff[2][nEffNuisances];
  TH1D* histoPtRap2RecVV_LepEff[2][nEffNuisances];
  TH1D* histoPtRap2RecDA_MomRes[2][nMomNuisances];
  TH1D* histoPtRap2RecDY_MomRes[2][nMomNuisances];
  TH1D* histoPtRap2RecEM_MomRes[2][nMomNuisances];
  TH1D* histoPtRap2RecVV_MomRes[2][nMomNuisances];
  TH1D* histoPtRap2RecVV_PDF[2];
  TH1D* histoPtRap2RecDY_PDF[2];
  TH1D* histoPtRap2RecVV_QCD[2];
  TH1D* histoPtRap2RecDY_QCD[2];
  TH1D* histoPtRap2RecVV_QCDPart[2][6];
  TH1D* histoPtRap2RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRap2RecGen[i]        = new TH2D(Form("histoPtRap2RecGen_%d",i),        Form("histoPtRap2RecGen_%d",i),        nBinRecoPtRap2, xbinsRecoPtRap2, nBinPtRap2, xbinsPtRap2);
   for(int j=0; j<nRecNuisances; j++){
     histoPtRap2RecGen_RecEff[i][j] = new TH2D(Form("histoPtRap2RecGen_RecEff_%d_%d",i,j), Form("histoPtRap2RecGen_RecEff_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2, nBinPtRap2, xbinsPtRap2);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRap2RecGen_LepEff[i][j] = new TH2D(Form("histoPtRap2RecGen_LepEff_%d_%d",i,j), Form("histoPtRap2RecGen_LepEff_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2, nBinPtRap2, xbinsPtRap2);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRap2RecGen_MomRes[i][j] = new TH2D(Form("histoPtRap2RecGen_MomRes_%d_%d",i,j), Form("histoPtRap2RecGen_MomRes_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2, nBinPtRap2, xbinsPtRap2);
   }
   histoPtRap2RecDA[i]  = new TH1D(Form("histoPtRap2RecDA_%d",i),  Form("histoPtRap2RecDA_%d",i),  nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecDY[i]  = new TH1D(Form("histoPtRap2RecDY_%d",i),  Form("histoPtRap2RecDY_%d",i),  nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecEM[i]  = new TH1D(Form("histoPtRap2RecEM_%d",i),  Form("histoPtRap2RecEM_%d",i),  nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecVV[i]  = new TH1D(Form("histoPtRap2RecVV_%d",i),  Form("histoPtRap2RecVV_%d",i),  nBinRecoPtRap2, xbinsRecoPtRap2);

   for(int j=0; j<nRecNuisances; j++){
     histoPtRap2RecDY_RecEff[i][j] = new TH1D(Form("histoPtRap2RecDY_RecEff_%d_%d",i,j), Form("histoPtRap2RecDY_RecEff_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
     histoPtRap2RecVV_RecEff[i][j] = new TH1D(Form("histoPtRap2RecVV_RecEff_%d_%d",i,j), Form("histoPtRap2RecVV_RecEff_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRap2RecDY_LepEff[i][j] = new TH1D(Form("histoPtRap2RecDY_LepEff_%d_%d",i,j), Form("histoPtRap2RecDY_LepEff_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
     histoPtRap2RecVV_LepEff[i][j] = new TH1D(Form("histoPtRap2RecVV_LepEff_%d_%d",i,j), Form("histoPtRap2RecVV_LepEff_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRap2RecDA_MomRes[i][j] = new TH1D(Form("histoPtRap2RecDA_MomRes_%d_%d",i,j), Form("histoPtRap2RecDA_MomRes_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
     histoPtRap2RecDY_MomRes[i][j] = new TH1D(Form("histoPtRap2RecDY_MomRes_%d_%d",i,j), Form("histoPtRap2RecDY_MomRes_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
     histoPtRap2RecEM_MomRes[i][j] = new TH1D(Form("histoPtRap2RecEM_MomRes_%d_%d",i,j), Form("histoPtRap2RecEM_MomRes_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
     histoPtRap2RecVV_MomRes[i][j] = new TH1D(Form("histoPtRap2RecVV_MomRes_%d_%d",i,j), Form("histoPtRap2RecVV_MomRes_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
   }
   histoPtRap2RecVV_PDF[i] = new TH1D(Form("histoPtRap2RecVV_PDF_%d",i), Form("histoPtRap2RecVV_PDF_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecDY_PDF[i] = new TH1D(Form("histoPtRap2RecDY_PDF_%d",i), Form("histoPtRap2RecDY_PDF_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);

   histoPtRap2RecVV_QCD[i] = new TH1D(Form("histoPtRap2RecVV_QCD_%d",i), Form("histoPtRap2RecVV_QCD_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);
   histoPtRap2RecDY_QCD[i] = new TH1D(Form("histoPtRap2RecDY_QCD_%d",i), Form("histoPtRap2RecDY_QCD_%d",i), nBinRecoPtRap2, xbinsRecoPtRap2);

   for(int j=0; j<6; j++){
      histoPtRap2RecVV_QCDPart[i][j] = new TH1D(Form("histoPtRap2RecVV_QCDPart_%d_%d",i,j), Form("histoPtRap2RecVV_QCDPart_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
      histoPtRap2RecDY_QCDPart[i][j] = new TH1D(Form("histoPtRap2RecDY_QCDPart_%d_%d",i,j), Form("histoPtRap2RecDY_QCDPart_%d_%d",i,j), nBinRecoPtRap2, xbinsRecoPtRap2);
    }
  }

  TH2D* histoPtRap3RecGen[2];
  TH2D* histoPtRap3RecGen_RecEff[2][nRecNuisances];
  TH2D* histoPtRap3RecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRap3RecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRap3RecDA[2];
  TH1D* histoPtRap3RecDY[2];
  TH1D* histoPtRap3RecEM[2];
  TH1D* histoPtRap3RecVV[2];
  TH1D* histoPtRap3RecDY_RecEff[2][nRecNuisances];
  TH1D* histoPtRap3RecVV_RecEff[2][nEffNuisances];
  TH1D* histoPtRap3RecDY_LepEff[2][nEffNuisances];
  TH1D* histoPtRap3RecVV_LepEff[2][nEffNuisances];
  TH1D* histoPtRap3RecDA_MomRes[2][nMomNuisances];
  TH1D* histoPtRap3RecDY_MomRes[2][nMomNuisances];
  TH1D* histoPtRap3RecEM_MomRes[2][nMomNuisances];
  TH1D* histoPtRap3RecVV_MomRes[2][nMomNuisances];
  TH1D* histoPtRap3RecVV_PDF[2];
  TH1D* histoPtRap3RecDY_PDF[2];
  TH1D* histoPtRap3RecVV_QCD[2];
  TH1D* histoPtRap3RecDY_QCD[2];
  TH1D* histoPtRap3RecVV_QCDPart[2][6];
  TH1D* histoPtRap3RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRap3RecGen[i]        = new TH2D(Form("histoPtRap3RecGen_%d",i),        Form("histoPtRap3RecGen_%d",i),        nBinRecoPtRap3, xbinsRecoPtRap3, nBinPtRap3, xbinsPtRap3);
   for(int j=0; j<nRecNuisances; j++){
     histoPtRap3RecGen_RecEff[i][j] = new TH2D(Form("histoPtRap3RecGen_RecEff_%d_%d",i,j), Form("histoPtRap3RecGen_RecEff_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3, nBinPtRap3, xbinsPtRap3);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRap3RecGen_LepEff[i][j] = new TH2D(Form("histoPtRap3RecGen_LepEff_%d_%d",i,j), Form("histoPtRap3RecGen_LepEff_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3, nBinPtRap3, xbinsPtRap3);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRap3RecGen_MomRes[i][j] = new TH2D(Form("histoPtRap3RecGen_MomRes_%d_%d",i,j), Form("histoPtRap3RecGen_MomRes_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3, nBinPtRap3, xbinsPtRap3);
   }
   histoPtRap3RecDA[i]  = new TH1D(Form("histoPtRap3RecDA_%d",i),  Form("histoPtRap3RecDA_%d",i),  nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecDY[i]  = new TH1D(Form("histoPtRap3RecDY_%d",i),  Form("histoPtRap3RecDY_%d",i),  nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecEM[i]  = new TH1D(Form("histoPtRap3RecEM_%d",i),  Form("histoPtRap3RecEM_%d",i),  nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecVV[i]  = new TH1D(Form("histoPtRap3RecVV_%d",i),  Form("histoPtRap3RecVV_%d",i),  nBinRecoPtRap3, xbinsRecoPtRap3);

   for(int j=0; j<nRecNuisances; j++){
     histoPtRap3RecDY_RecEff[i][j] = new TH1D(Form("histoPtRap3RecDY_RecEff_%d_%d",i,j), Form("histoPtRap3RecDY_RecEff_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
     histoPtRap3RecVV_RecEff[i][j] = new TH1D(Form("histoPtRap3RecVV_RecEff_%d_%d",i,j), Form("histoPtRap3RecVV_RecEff_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRap3RecDY_LepEff[i][j] = new TH1D(Form("histoPtRap3RecDY_LepEff_%d_%d",i,j), Form("histoPtRap3RecDY_LepEff_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
     histoPtRap3RecVV_LepEff[i][j] = new TH1D(Form("histoPtRap3RecVV_LepEff_%d_%d",i,j), Form("histoPtRap3RecVV_LepEff_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRap3RecDA_MomRes[i][j] = new TH1D(Form("histoPtRap3RecDA_MomRes_%d_%d",i,j), Form("histoPtRap3RecDA_MomRes_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
     histoPtRap3RecDY_MomRes[i][j] = new TH1D(Form("histoPtRap3RecDY_MomRes_%d_%d",i,j), Form("histoPtRap3RecDY_MomRes_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
     histoPtRap3RecEM_MomRes[i][j] = new TH1D(Form("histoPtRap3RecEM_MomRes_%d_%d",i,j), Form("histoPtRap3RecEM_MomRes_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
     histoPtRap3RecVV_MomRes[i][j] = new TH1D(Form("histoPtRap3RecVV_MomRes_%d_%d",i,j), Form("histoPtRap3RecVV_MomRes_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
   }
   histoPtRap3RecVV_PDF[i] = new TH1D(Form("histoPtRap3RecVV_PDF_%d",i), Form("histoPtRap3RecVV_PDF_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecDY_PDF[i] = new TH1D(Form("histoPtRap3RecDY_PDF_%d",i), Form("histoPtRap3RecDY_PDF_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);

   histoPtRap3RecVV_QCD[i] = new TH1D(Form("histoPtRap3RecVV_QCD_%d",i), Form("histoPtRap3RecVV_QCD_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);
   histoPtRap3RecDY_QCD[i] = new TH1D(Form("histoPtRap3RecDY_QCD_%d",i), Form("histoPtRap3RecDY_QCD_%d",i), nBinRecoPtRap3, xbinsRecoPtRap3);

   for(int j=0; j<6; j++){
      histoPtRap3RecVV_QCDPart[i][j] = new TH1D(Form("histoPtRap3RecVV_QCDPart_%d_%d",i,j), Form("histoPtRap3RecVV_QCDPart_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
      histoPtRap3RecDY_QCDPart[i][j] = new TH1D(Form("histoPtRap3RecDY_QCDPart_%d_%d",i,j), Form("histoPtRap3RecDY_QCDPart_%d_%d",i,j), nBinRecoPtRap3, xbinsRecoPtRap3);
    }
  }

  TH2D* histoPtRap4RecGen[2];
  TH2D* histoPtRap4RecGen_RecEff[2][nRecNuisances];
  TH2D* histoPtRap4RecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRap4RecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRap4RecDA[2];
  TH1D* histoPtRap4RecDY[2];
  TH1D* histoPtRap4RecEM[2];
  TH1D* histoPtRap4RecVV[2];
  TH1D* histoPtRap4RecDY_RecEff[2][nRecNuisances];
  TH1D* histoPtRap4RecVV_RecEff[2][nEffNuisances];
  TH1D* histoPtRap4RecDY_LepEff[2][nEffNuisances];
  TH1D* histoPtRap4RecVV_LepEff[2][nEffNuisances];
  TH1D* histoPtRap4RecDA_MomRes[2][nMomNuisances];
  TH1D* histoPtRap4RecDY_MomRes[2][nMomNuisances];
  TH1D* histoPtRap4RecEM_MomRes[2][nMomNuisances];
  TH1D* histoPtRap4RecVV_MomRes[2][nMomNuisances];
  TH1D* histoPtRap4RecVV_PDF[2];
  TH1D* histoPtRap4RecDY_PDF[2];
  TH1D* histoPtRap4RecVV_QCD[2];
  TH1D* histoPtRap4RecDY_QCD[2];
  TH1D* histoPtRap4RecVV_QCDPart[2][6];
  TH1D* histoPtRap4RecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoPtRap4RecGen[i]        = new TH2D(Form("histoPtRap4RecGen_%d",i),        Form("histoPtRap4RecGen_%d",i),        nBinRecoPtRap4, xbinsRecoPtRap4, nBinPtRap4, xbinsPtRap4);
   for(int j=0; j<nRecNuisances; j++){
     histoPtRap4RecGen_RecEff[i][j] = new TH2D(Form("histoPtRap4RecGen_RecEff_%d_%d",i,j), Form("histoPtRap4RecGen_RecEff_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4, nBinPtRap4, xbinsPtRap4);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRap4RecGen_LepEff[i][j] = new TH2D(Form("histoPtRap4RecGen_LepEff_%d_%d",i,j), Form("histoPtRap4RecGen_LepEff_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4, nBinPtRap4, xbinsPtRap4);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRap4RecGen_MomRes[i][j] = new TH2D(Form("histoPtRap4RecGen_MomRes_%d_%d",i,j), Form("histoPtRap4RecGen_MomRes_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4, nBinPtRap4, xbinsPtRap4);
   }
   histoPtRap4RecDA[i]  = new TH1D(Form("histoPtRap4RecDA_%d",i),  Form("histoPtRap4RecDA_%d",i),  nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecDY[i]  = new TH1D(Form("histoPtRap4RecDY_%d",i),  Form("histoPtRap4RecDY_%d",i),  nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecEM[i]  = new TH1D(Form("histoPtRap4RecEM_%d",i),  Form("histoPtRap4RecEM_%d",i),  nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecVV[i]  = new TH1D(Form("histoPtRap4RecVV_%d",i),  Form("histoPtRap4RecVV_%d",i),  nBinRecoPtRap4, xbinsRecoPtRap4);

   for(int j=0; j<nRecNuisances; j++){
     histoPtRap4RecDY_RecEff[i][j] = new TH1D(Form("histoPtRap4RecDY_RecEff_%d_%d",i,j), Form("histoPtRap4RecDY_RecEff_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
     histoPtRap4RecVV_RecEff[i][j] = new TH1D(Form("histoPtRap4RecVV_RecEff_%d_%d",i,j), Form("histoPtRap4RecVV_RecEff_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
   }
   for(int j=0; j<nEffNuisances; j++){
     histoPtRap4RecDY_LepEff[i][j] = new TH1D(Form("histoPtRap4RecDY_LepEff_%d_%d",i,j), Form("histoPtRap4RecDY_LepEff_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
     histoPtRap4RecVV_LepEff[i][j] = new TH1D(Form("histoPtRap4RecVV_LepEff_%d_%d",i,j), Form("histoPtRap4RecVV_LepEff_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoPtRap4RecDA_MomRes[i][j] = new TH1D(Form("histoPtRap4RecDA_MomRes_%d_%d",i,j), Form("histoPtRap4RecDA_MomRes_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
     histoPtRap4RecDY_MomRes[i][j] = new TH1D(Form("histoPtRap4RecDY_MomRes_%d_%d",i,j), Form("histoPtRap4RecDY_MomRes_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
     histoPtRap4RecEM_MomRes[i][j] = new TH1D(Form("histoPtRap4RecEM_MomRes_%d_%d",i,j), Form("histoPtRap4RecEM_MomRes_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
     histoPtRap4RecVV_MomRes[i][j] = new TH1D(Form("histoPtRap4RecVV_MomRes_%d_%d",i,j), Form("histoPtRap4RecVV_MomRes_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
   }
   histoPtRap4RecVV_PDF[i] = new TH1D(Form("histoPtRap4RecVV_PDF_%d",i), Form("histoPtRap4RecVV_PDF_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecDY_PDF[i] = new TH1D(Form("histoPtRap4RecDY_PDF_%d",i), Form("histoPtRap4RecDY_PDF_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);

   histoPtRap4RecVV_QCD[i] = new TH1D(Form("histoPtRap4RecVV_QCD_%d",i), Form("histoPtRap4RecVV_QCD_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);
   histoPtRap4RecDY_QCD[i] = new TH1D(Form("histoPtRap4RecDY_QCD_%d",i), Form("histoPtRap4RecDY_QCD_%d",i), nBinRecoPtRap4, xbinsRecoPtRap4);

   for(int j=0; j<6; j++){
      histoPtRap4RecVV_QCDPart[i][j] = new TH1D(Form("histoPtRap4RecVV_QCDPart_%d_%d",i,j), Form("histoPtRap4RecVV_QCDPart_%d_%d",i,j), nBinRecoPtRap4, xbinsRecoPtRap4);
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
      //if(thePandaFlat.metFilter == 0) continue;
      
      if(thePandaFlat.nLooseLep != 2) continue;

      if(TMath::Abs(thePandaFlat.looseLep1Eta) >= 2.4 || TMath::Abs(thePandaFlat.looseLep2Eta) >= 2.4) continue;

      int theCategory = infileCat_[ifile];

      int lepType = -1;
      if     (abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 0;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 1;
      else if(abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 2;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 2;
      else if(abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 3;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 4;
      else if(abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 5;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 5;
      else {printf("Impossible dilepton combination: %d %d\n",thePandaFlat.looseLep1PdgId,thePandaFlat.looseLep2PdgId); continue;}

      if(lepType >= 3 || (lepType == 2 && theCategory != 0)) continue;
      
      if(lepType == 2 && theCategory == 0) theCategory = 5; // using data e-mu events to estimate non-resonant background

      double thePDGMass[2] = {mass_mu, mass_mu};
      if     (abs(lepType) == 1) {thePDGMass[0] = mass_el; thePDGMass[1] = mass_el;}
      else if(abs(lepType) == 2 && abs(thePandaFlat.looseLep1PdgId)==11) {thePDGMass[0] = mass_el;}
      else if(abs(lepType) == 2 && abs(thePandaFlat.looseLep2PdgId)==11) {thePDGMass[1] = mass_el;}

      TLorentzVector vLoose1,vLoose2;
      vLoose1.SetPtEtaPhiM(thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      vLoose2.SetPtEtaPhiM(thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      bool passLooseSel = ((thePandaFlat.looseLep1SelBit & kLoose) == kLoose) && ((thePandaFlat.looseLep2SelBit & kLoose) == kLoose) && 
                          TMath::Abs((vLoose1+vLoose2).M()-91.1876) < 15 && vLoose1.Pt() > 35 && vLoose2.Pt() > 35;

      double the_eta_sf[2] = {1.0, 1.0}; double the_eta_sf_unc[2] = {0.0, 0.0};
      if(theCategory != 0 && theCategory != 5){
        if(abs(thePandaFlat.looseLep1PdgId)==13){
          double etal = thePandaFlat.looseLep1Eta; if(etal >= 2.4) etal = 2.3999; else if(etal <= -2.4) etal = -2.3999;
          int binEta = scalefactors_Muon_Eta->GetXaxis()->FindFixBin(etal);
          the_eta_sf[0] = scalefactors_Muon_Eta->GetBinContent(binEta);
	  the_eta_sf_unc[0] = scalefactors_Muon_Eta->GetBinError(binEta)/scalefactors_Muon_Eta->GetBinContent(binEta);
        } else {
          double etal = thePandaFlat.looseLep1Eta; if(etal >= 2.5) etal = 2.4999; else if(etal <= -2.5) etal = -2.4999;
          int binEta = scalefactors_Electron_Eta->GetXaxis()->FindFixBin(etal);
          the_eta_sf[0] = scalefactors_Electron_Eta->GetBinContent(binEta);
	  the_eta_sf_unc[0] = scalefactors_Electron_Eta->GetBinError(binEta)/scalefactors_Electron_Eta->GetBinContent(binEta);
	  if(thePandaFlat.looseLep1Pt < 20 || thePandaFlat.looseLep1Pt > 80) the_eta_sf_unc[0] = sqrt(the_eta_sf_unc[0]*the_eta_sf_unc[0] + 0.01*0.01);
        }        
        if(abs(thePandaFlat.looseLep2PdgId)==13){
          double etal = thePandaFlat.looseLep2Eta; if(etal >= 2.4) etal = 2.3999; else if(etal <= -2.4) etal = -2.3999;
          int binEta = scalefactors_Muon_Eta->GetXaxis()->FindFixBin(etal);
          the_eta_sf[1] = scalefactors_Muon_Eta->GetBinContent(binEta);
	  the_eta_sf_unc[1] = scalefactors_Muon_Eta->GetBinError(binEta)/scalefactors_Muon_Eta->GetBinContent(binEta);
        } else {
          double etal = thePandaFlat.looseLep2Eta; if(etal >= 2.5) etal = 2.4999; else if(etal <= -2.5) etal = -2.4999;
          int binEta = scalefactors_Electron_Eta->GetXaxis()->FindFixBin(etal);
          the_eta_sf[1] = scalefactors_Electron_Eta->GetBinContent(binEta);
	  the_eta_sf_unc[1] = scalefactors_Electron_Eta->GetBinError(binEta)/scalefactors_Electron_Eta->GetBinContent(binEta);
	  if(thePandaFlat.looseLep2Pt < 20 || thePandaFlat.looseLep2Pt > 80) the_eta_sf_unc[1] = sqrt(the_eta_sf_unc[1]*the_eta_sf_unc[1] + 0.01*0.01);
        }
      }

      if(passLooseSel) {
        double totalLooseWeight = 1.0;
        if(theCategory != 0 && theCategory != 5){
          totalLooseWeight = thePandaFlat.normalizedWeight * lumi * thePandaFlat.sf_pu *
		             the_eta_sf[0] * the_eta_sf[1] *
		             trigger_sf(thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1PdgId,thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2PdgId);
        }
        if     (theCategory == 0){
	  histo[lepType+22][theCategory]->Fill((vLoose1+vLoose2).Rapidity(),1.0);
	  histo[lepType+24][theCategory]->Fill(TMath::Min(thePandaFlat.looseLep1Eta,thePandaFlat.looseLep2Eta),1.0);
	  histo[lepType+26][theCategory]->Fill(TMath::Max(thePandaFlat.looseLep1Eta,thePandaFlat.looseLep2Eta),1.0);
	  histo[lepType+28][theCategory]->Fill(thePandaFlat.looseLep1Eta,1.0);
	  histo[lepType+28][theCategory]->Fill(thePandaFlat.looseLep2Eta,1.0);
        }
        else if(theCategory != 5){
	  histo[lepType+22][theCategory]->Fill((vLoose1+vLoose2).Rapidity(),totalLooseWeight);
	  histo[lepType+24][theCategory]->Fill(TMath::Min(thePandaFlat.looseLep1Eta,thePandaFlat.looseLep2Eta),totalLooseWeight);
	  histo[lepType+26][theCategory]->Fill(TMath::Max(thePandaFlat.looseLep1Eta,thePandaFlat.looseLep2Eta),totalLooseWeight);
	  histo[lepType+28][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalLooseWeight);
	  histo[lepType+28][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalLooseWeight);
        }
	else {
          for(int theLepType = 0; theLepType<2; theLepType++) {
            double theKeff = k_eff;
            if(theLepType == 1) theKeff = 1.0/k_eff;
	    histo[theLepType+22][theCategory]->Fill((vLoose1+vLoose2).Rapidity(),1.0*theKeff);
	    histo[theLepType+24][theCategory]->Fill(TMath::Min(thePandaFlat.looseLep1Eta,thePandaFlat.looseLep2Eta),1.0*theKeff);
	    histo[theLepType+26][theCategory]->Fill(TMath::Max(thePandaFlat.looseLep1Eta,thePandaFlat.looseLep2Eta),1.0*theKeff);
	    histo[theLepType+28][theCategory]->Fill(thePandaFlat.looseLep1Eta,1.0*theKeff);
	    histo[theLepType+28][theCategory]->Fill(thePandaFlat.looseLep2Eta,1.0*theKeff);
          }
	}
      }

      bool passLepId = ((thePandaFlat.looseLep1SelBit & kMedium) == kMedium) && ((thePandaFlat.looseLep2SelBit & kMedium) == kMedium);
      if(passLepId == false) continue;

      float lepPtSF[2] = {1,1};
      float lepPtSFSyst[2][nMomNuisances] = {1,1,1,1,1,1,1,1,1,1};
      if(theCategory == 0 || theCategory == 5) { // Data
        if(abs(thePandaFlat.looseLep1PdgId)==13) {
	  double rnd[2] = {gRandom->Rndm(), gRandom->Rndm()};
          lepPtSF[0] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 0, 0);
          double muonPtSystA[100],muonPtSystB[1],muonPtSystC[5],muonPtSystD[5];
          double muonPtSystAll[4] = {0,0,0,0};

          for(int i=0; i<100; i++)
            muonPtSystA[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 1, i);
          for(int i=0; i<1; i++)
            muonPtSystB[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 2, i);
          for(int i=0; i<5; i++)
            muonPtSystC[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 4, i);
          for(int i=0; i<5; i++)
            muonPtSystD[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 5, i);

          for(int i=0; i<100; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystA[i]) > muonPtSystAll[0]) muonPtSystAll[0] = TMath::Abs(lepPtSF[0]-muonPtSystA[i]);
          for(int i=0; i<1  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystB[i]) > muonPtSystAll[1]) muonPtSystAll[1] = TMath::Abs(lepPtSF[0]-muonPtSystB[i]);
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystC[i]) > muonPtSystAll[2]) muonPtSystAll[2] = TMath::Abs(lepPtSF[0]-muonPtSystC[i]);
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystD[i]) > muonPtSystAll[3]) muonPtSystAll[3] = TMath::Abs(lepPtSF[0]-muonPtSystD[i]);
          lepPtSFSyst[0][0] = sqrt(muonPtSystAll[0]*muonPtSystAll[0]+muonPtSystAll[1]*muonPtSystAll[1]+
                                   muonPtSystAll[2]*muonPtSystAll[2]+muonPtSystAll[3]*muonPtSystAll[3])/lepPtSF[0] + 1.0;
          lepPtSFSyst[0][1] = muonPtSystAll[0]/lepPtSF[0] + 1.0;
          lepPtSFSyst[0][2] = muonPtSystAll[1]/lepPtSF[0] + 1.0;
          lepPtSFSyst[0][3] = muonPtSystAll[2]/lepPtSF[0] + 1.0;
          lepPtSFSyst[0][4] = muonPtSystAll[3]/lepPtSF[0] + 1.0;
        }
        if(abs(thePandaFlat.looseLep2PdgId)==13) {
	  double rnd[2] = {gRandom->Rndm(), gRandom->Rndm()};
	  lepPtSF[1] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 0, 0);
          double muonPtSystA[100],muonPtSystB[1],muonPtSystC[5],muonPtSystD[5];
          double muonPtSystAll[4] = {0,0,0,0};

          for(int i=0; i<100; i++)
            muonPtSystA[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 1, i);
          for(int i=0; i<1; i++)
            muonPtSystB[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 2, i);
          for(int i=0; i<5; i++)
            muonPtSystC[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 4, i);
          for(int i=0; i<5; i++)
            muonPtSystD[i] = rmcor.kScaleDT(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 5, i);

          for(int i=0; i<100; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystA[i]) > muonPtSystAll[0]) muonPtSystAll[0] = TMath::Abs(lepPtSF[1]-muonPtSystA[i]);
          for(int i=0; i<1  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystB[i]) > muonPtSystAll[1]) muonPtSystAll[1] = TMath::Abs(lepPtSF[1]-muonPtSystB[i]);
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystC[i]) > muonPtSystAll[2]) muonPtSystAll[2] = TMath::Abs(lepPtSF[1]-muonPtSystC[i]);
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystD[i]) > muonPtSystAll[3]) muonPtSystAll[3] = TMath::Abs(lepPtSF[1]-muonPtSystD[i]);
          lepPtSFSyst[1][0] = sqrt(muonPtSystAll[0]*muonPtSystAll[0]+muonPtSystAll[1]*muonPtSystAll[1]+
                                   muonPtSystAll[2]*muonPtSystAll[2]+muonPtSystAll[3]*muonPtSystAll[3])/lepPtSF[1] + 1.0;
          lepPtSFSyst[1][1] = muonPtSystAll[0]/lepPtSF[1] + 1.0;
          lepPtSFSyst[1][2] = muonPtSystAll[1]/lepPtSF[1] + 1.0;
          lepPtSFSyst[1][3] = muonPtSystAll[2]/lepPtSF[1] + 1.0;
          lepPtSFSyst[1][4] = muonPtSystAll[3]/lepPtSF[1] + 1.0;
        }
      } else { // MC
        if(abs(thePandaFlat.looseLep1PdgId)==13) {
	  double rnd[2] = {gRandom->Rndm(), gRandom->Rndm()};
          lepPtSF[0] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, rnd[0], rnd[1], 0, 0);
          double muonPtSyst[11];
          double muonPtSystA[100],muonPtSystB[1],muonPtSystC[5],muonPtSystD[5];
          double muonPtSystAll[4] = {0,0,0,0};

          for(int i=0; i<100; i++)
            muonPtSystA[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, rnd[0], rnd[1], 1, i);

          for(int i=0; i<1; i++)
            muonPtSystB[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, rnd[0], rnd[1], 2, i);
          for(int i=0; i<5; i++)
            muonPtSystC[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, rnd[0], rnd[1], 4, i);
          for(int i=0; i<5; i++)
            muonPtSystD[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep1PdgId)/thePandaFlat.looseLep1PdgId, thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi, 10, rnd[0], rnd[1], 5, i);

          for(int i=0; i<100; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystA[i]) > muonPtSystAll[0]) muonPtSystAll[0] = TMath::Abs(lepPtSF[0]-muonPtSystA[i]);
          for(int i=0; i<1  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystB[i]) > muonPtSystAll[1]) muonPtSystAll[1] = TMath::Abs(lepPtSF[0]-muonPtSystB[i]);
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystC[i]) > muonPtSystAll[2]) muonPtSystAll[2] = TMath::Abs(lepPtSF[0]-muonPtSystC[i]);
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[0]-muonPtSystD[i]) > muonPtSystAll[3]) muonPtSystAll[3] = TMath::Abs(lepPtSF[0]-muonPtSystD[i]);
          lepPtSFSyst[0][0] = sqrt(muonPtSystAll[0]*muonPtSystAll[0]+muonPtSystAll[1]*muonPtSystAll[1]+
                                   muonPtSystAll[2]*muonPtSystAll[2]+muonPtSystAll[3]*muonPtSystAll[3])/lepPtSF[0] + 1.0;
          lepPtSFSyst[0][1] = muonPtSystAll[0]/lepPtSF[0] + 1.0;
          lepPtSFSyst[0][2] = muonPtSystAll[1]/lepPtSF[0] + 1.0;
          lepPtSFSyst[0][3] = muonPtSystAll[2]/lepPtSF[0] + 1.0;
          lepPtSFSyst[0][4] = muonPtSystAll[3]/lepPtSF[0] + 1.0;
        }
        if(abs(thePandaFlat.looseLep2PdgId)==13) {
	  double rnd[2] = {gRandom->Rndm(), gRandom->Rndm()};
	  lepPtSF[1] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 10, rnd[0], rnd[1], 0, 0);
          double muonPtSyst[11];
          double muonPtSystA[100],muonPtSystB[1],muonPtSystC[5],muonPtSystD[5];
          double muonPtSystAll[4] = {0,0,0,0};

          for(int i=0; i<100; i++)
            muonPtSystA[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 10, rnd[0], rnd[1], 1, i);

          for(int i=0; i<1; i++)
            muonPtSystB[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 10, rnd[0], rnd[1], 2, i);
          for(int i=0; i<5; i++)
            muonPtSystC[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 10, rnd[0], rnd[1], 4, i);
          for(int i=0; i<5; i++)
            muonPtSystD[i] = rmcor.kScaleAndSmearMC(-1*abs(thePandaFlat.looseLep2PdgId)/thePandaFlat.looseLep2PdgId, thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi, 10, rnd[0], rnd[1], 5, i);

          for(int i=0; i<100; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystA[i]) > muonPtSystAll[0]) muonPtSystAll[0] = TMath::Abs(lepPtSF[1]-muonPtSystA[i]);
          for(int i=0; i<1  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystB[i]) > muonPtSystAll[1]) muonPtSystAll[1] = TMath::Abs(lepPtSF[1]-muonPtSystB[i]);
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystC[i]) > muonPtSystAll[2]) muonPtSystAll[2] = TMath::Abs(lepPtSF[1]-muonPtSystC[i]);
          for(int i=0; i<5  ; i++) if(TMath::Abs(lepPtSF[1]-muonPtSystD[i]) > muonPtSystAll[3]) muonPtSystAll[3] = TMath::Abs(lepPtSF[1]-muonPtSystD[i]);
          lepPtSFSyst[1][0] = sqrt(muonPtSystAll[0]*muonPtSystAll[0]+muonPtSystAll[1]*muonPtSystAll[1]+
                                   muonPtSystAll[2]*muonPtSystAll[2]+muonPtSystAll[3]*muonPtSystAll[3])/lepPtSF[1] + 1.0;
          lepPtSFSyst[1][1] = muonPtSystAll[0]/lepPtSF[1] + 1.0;
          lepPtSFSyst[1][2] = muonPtSystAll[1]/lepPtSF[1] + 1.0;
          lepPtSFSyst[1][3] = muonPtSystAll[2]/lepPtSF[1] + 1.0;
          lepPtSFSyst[1][4] = muonPtSystAll[3]/lepPtSF[1] + 1.0;
        }
        if(abs(thePandaFlat.looseLep1PdgId)==11) {
          if(TMath::Abs(thePandaFlat.looseLep1Eta) <  1.5) lepPtSF[0] = gRandom->Gaus(1.000,0.013);
          else                                             lepPtSF[0] = gRandom->Gaus(0.993,0.025);
          lepPtSFSyst[0][0] = lepPtSF[0];
          lepPtSFSyst[0][1] = lepPtSF[0];
        }
        if(abs(thePandaFlat.looseLep2PdgId)==11) {
          if(TMath::Abs(thePandaFlat.looseLep2Eta) <  1.5) lepPtSF[1] = gRandom->Gaus(1.000,0.013);
          else                                             lepPtSF[1] = gRandom->Gaus(0.993,0.025);
          lepPtSFSyst[1][0] = lepPtSF[1];
          lepPtSFSyst[1][1] = lepPtSF[1];
        }
      }

      TLorentzVector v1,v2;
      v1.SetPtEtaPhiM(thePandaFlat.looseLep1Pt*lepPtSF[0],thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      v2.SetPtEtaPhiM(thePandaFlat.looseLep2Pt*lepPtSF[1],thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      TLorentzVector vMomRes1[nMomNuisances],vMomRes2[nMomNuisances];
      vMomRes1[0].SetPtEtaPhiM(thePandaFlat.looseLep1Pt*lepPtSF[0]*lepPtSFSyst[0][0],thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      vMomRes2[0].SetPtEtaPhiM(thePandaFlat.looseLep2Pt*lepPtSF[1]*lepPtSFSyst[1][0],thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      vMomRes1[1].SetPtEtaPhiM(thePandaFlat.looseLep1Pt*lepPtSF[0]*lepPtSFSyst[0][1],thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      vMomRes2[1].SetPtEtaPhiM(thePandaFlat.looseLep2Pt*lepPtSF[1]*lepPtSFSyst[1][1],thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      vMomRes1[2].SetPtEtaPhiM(thePandaFlat.looseLep1Pt*lepPtSF[0]*lepPtSFSyst[0][2],thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      vMomRes2[2].SetPtEtaPhiM(thePandaFlat.looseLep2Pt*lepPtSF[1]*lepPtSFSyst[1][2],thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      vMomRes1[3].SetPtEtaPhiM(thePandaFlat.looseLep1Pt*lepPtSF[0]*lepPtSFSyst[0][3],thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      vMomRes2[3].SetPtEtaPhiM(thePandaFlat.looseLep2Pt*lepPtSF[1]*lepPtSFSyst[1][3],thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      vMomRes1[4].SetPtEtaPhiM(thePandaFlat.looseLep1Pt*lepPtSF[0]*lepPtSFSyst[0][4],thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      vMomRes2[4].SetPtEtaPhiM(thePandaFlat.looseLep2Pt*lepPtSF[1]*lepPtSFSyst[1][4],thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      bool passSel = TMath::Abs((v1+v2).M()-91.1876) < 15 && v1.Pt() > 25 && v2.Pt() > 25 && (v1+v2).Pt() > dileptonPtCut;
      bool passSystSel[nMomNuisances] = {TMath::Abs((vMomRes1[0]+vMomRes2[0]).M()-91.1876) < 15 && vMomRes1[0].Pt() > 25 && vMomRes2[0].Pt() > 25 && (vMomRes1[0]+vMomRes2[0]).Pt() > dileptonPtCut,
                                         TMath::Abs((vMomRes1[1]+vMomRes2[1]).M()-91.1876) < 15 && vMomRes1[1].Pt() > 25 && vMomRes2[1].Pt() > 25 && (vMomRes1[0]+vMomRes2[0]).Pt() > dileptonPtCut,
                                         TMath::Abs((vMomRes1[2]+vMomRes2[2]).M()-91.1876) < 15 && vMomRes1[2].Pt() > 25 && vMomRes2[2].Pt() > 25 && (vMomRes1[0]+vMomRes2[0]).Pt() > dileptonPtCut,
                                         TMath::Abs((vMomRes1[3]+vMomRes2[3]).M()-91.1876) < 15 && vMomRes1[3].Pt() > 25 && vMomRes2[3].Pt() > 25 && (vMomRes1[0]+vMomRes2[0]).Pt() > dileptonPtCut,
                                         TMath::Abs((vMomRes1[4]+vMomRes2[4]).M()-91.1876) < 15 && vMomRes1[4].Pt() > 25 && vMomRes2[4].Pt() > 25 && (vMomRes1[0]+vMomRes2[0]).Pt() > dileptonPtCut};

      double ZRecTot  = 0.500;
      double ZRecPt  = TMath::Min((v1+v2).Pt(),1499.999);
      double ZRecPhiStar = TMath::Min(phi_star_eta(v1,v2,thePandaFlat.looseLep1PdgId),49.999);
      if(ZRecPhiStar <= 0.001) ZRecPhiStar = 0.001;
      double ZRecRap = TMath::Abs((v1+v2).Rapidity());

      double ZRecSystTot[nMomNuisances] = {0.500,0.500,0.500,0.500,0.500};
      double ZRecSystPt[nMomNuisances] = {(vMomRes1[0]+vMomRes2[0]).Pt(),
                                          (vMomRes1[1]+vMomRes2[1]).Pt(),
                                          (vMomRes1[2]+vMomRes2[2]).Pt(),
                                          (vMomRes1[3]+vMomRes2[3]).Pt(),
                                          (vMomRes1[4]+vMomRes2[4]).Pt()};

      double ZGenTot = 0; double ZGenPt = 0; double ZGenPhiStar = 0; double ZGenRap = 0; bool passPtFid = false; bool passRapFid = false; bool passPtRapFid[5] = {false, false, false, false, false}; 
      if(thePandaFlat.looseGenLep1PdgId != 0 && thePandaFlat.looseGenLep2PdgId != 0 &&
         thePandaFlat.genLep1Pt > 25 && TMath::Abs(thePandaFlat.genLep1Eta) < 2.4 &&
	 thePandaFlat.genLep2Pt > 25 && TMath::Abs(thePandaFlat.genLep2Eta) < 2.4){
        TLorentzVector vGen1,vGen2;
        vGen1.SetPtEtaPhiM(thePandaFlat.genLep1Pt,thePandaFlat.genLep1Eta,thePandaFlat.genLep1Phi,thePDGMass[0]);
        vGen2.SetPtEtaPhiM(thePandaFlat.genLep2Pt,thePandaFlat.genLep2Eta,thePandaFlat.genLep2Phi,thePDGMass[1]);
	if(TMath::Abs((vGen1+vGen2).M()-91.1876) < 15.0 && (vGen1+vGen2).Pt() >= dileptonPtCut) {
	  passPtFid = true;
	  ZGenTot = 0.500; 
	  ZGenPt = TMath::Min((vGen1+vGen2).Pt(),1499.999);
          ZGenPhiStar = TMath::Min(phi_star_eta(vGen1,vGen2,thePandaFlat.looseGenLep1PdgId),49.999);;
          if(ZGenPhiStar <= 0.001) ZGenPhiStar = 0.001;
          ZGenRap = TMath::Abs((vGen1+vGen2).Rapidity());
          if(ZGenRap < 2.4) { passRapFid = true;}
          if     (ZGenRap < 0.5) passPtRapFid[0] = true;
          else if(ZGenRap < 1.0) passPtRapFid[1] = true;
          else if(ZGenRap < 1.5) passPtRapFid[2] = true;
          //else if(ZGenRap < 2.0) passPtRapFid[3] = true;
          else if(ZGenRap < 2.4) passPtRapFid[3] = true;
          else if(ZGenRap < 2.4) passPtRapFid[4] = true;
	}
      }

      double totalWeight = 1.0; double sfWeightLepEff[2] = {1.0, 1.0};
      double sfSystWeightLepEff[2][nEffNuisances];
      for(int ni=0; ni<2; ni++) for(int nj=0; nj<nEffNuisances; nj++) sfSystWeightLepEff[ni][nj] = 0.0;

      if(theCategory != 0 && theCategory != 5){

        if(abs(thePandaFlat.looseLep1PdgId)==13){
          double etal = thePandaFlat.looseLep1Eta; if(etal >= 2.4) etal = 2.3999; else if(etal <= -2.4) etal = -2.3999;
          int binXT   = scalefactors_Medium_Muon->GetXaxis()->FindFixBin(etal);
          int binYT_c = scalefactors_Medium_Muon->GetYaxis()->FindFixBin(TMath::Min((double)thePandaFlat.looseLep1Pt,getMaxPtForSFs[0]));
          int binYT_s = scalefactors_Medium_Muon->GetYaxis()->FindFixBin(TMath::Min((double)thePandaFlat.looseLep1Pt,getMaxPtForSFs[2]));
	  sfWeightLepEff[0]        = scalefactors_Medium_Muon                   ->GetBinContent(binXT,binYT_c);
          sfSystWeightLepEff[0][0] = scalefactors_Medium_Muon                   ->GetBinError  (binXT,binYT_c);
          sfSystWeightLepEff[0][1] = scalefactors_Medium_Muon_stat_error_hi	->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[0][2] = scalefactors_Medium_Muon_signalFsrTNP	->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[0][3] = scalefactors_Medium_Muon_signalResTNP	->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[0][4] = scalefactors_Medium_Muon_bkgModelTNP	->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[0][5] = scalefactors_Medium_Muon_tagBiasTNP	->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[0][6] = scalefactors_Medium_Muon_generatorChoiceTNP->GetBinContent(binXT,binYT_s);
          //for(int nj=0; nj<nMuSFBins; nj++) sfSystWeightLepEff[0][nj+5] = scalefactors_Medium_Muon_stat_error_hi_bins[nj]->GetBinContent(binXT,binYT_s);
        } else {
          double etal = thePandaFlat.looseLep1Eta; if(etal >= 2.5) etal = 2.4999; else if(etal <= -2.5) etal = -2.4999;
          int binXT   = scalefactors_Medium_Electron->GetXaxis()->FindFixBin(etal);
          int binYT_c = scalefactors_Medium_Electron->GetYaxis()->FindFixBin(TMath::Min((double)thePandaFlat.looseLep1Pt,getMaxPtForSFs[1]));
          int binYT_s = scalefactors_Medium_Electron->GetYaxis()->FindFixBin(TMath::Min((double)thePandaFlat.looseLep1Pt,getMaxPtForSFs[3]));
	  sfWeightLepEff[0]        = scalefactors_Medium_Electron                   ->GetBinContent(binXT,binYT_c);
          sfSystWeightLepEff[0][0] = scalefactors_Medium_Electron                   ->GetBinError  (binXT,binYT_c);
          sfSystWeightLepEff[0][1] = scalefactors_Medium_Electron_stat_error_hi     ->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[0][2] = scalefactors_Medium_Electron_signalFsrTNP      ->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[0][3] = scalefactors_Medium_Electron_signalResTNP      ->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[0][4] = scalefactors_Medium_Electron_bkgModelTNP	    ->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[0][5] = scalefactors_Medium_Electron_tagBiasTNP	    ->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[0][6] = scalefactors_Medium_Electron_generatorChoiceTNP->GetBinContent(binXT,binYT_s);
          //for(int nj=0; nj<nElSFBins; nj++) sfSystWeightLepEff[0][nj+5] = scalefactors_Medium_Electron_stat_error_hi_bins[nj]->GetBinContent(binXT,binYT_s);
        }
        if(abs(thePandaFlat.looseLep2PdgId)==13){
          double etal = thePandaFlat.looseLep2Eta; if(etal >= 2.4) etal = 2.3999; else if(etal <= -2.4) etal = -2.3999;
          int binXT   = scalefactors_Medium_Muon->GetXaxis()->FindFixBin(etal);
          int binYT_c = scalefactors_Medium_Muon->GetYaxis()->FindFixBin(TMath::Min((double)thePandaFlat.looseLep2Pt,getMaxPtForSFs[0]));
          int binYT_s = scalefactors_Medium_Muon->GetYaxis()->FindFixBin(TMath::Min((double)thePandaFlat.looseLep2Pt,getMaxPtForSFs[2]));
	  sfWeightLepEff[1]        = scalefactors_Medium_Muon                   ->GetBinContent(binXT,binYT_c);
          sfSystWeightLepEff[1][0] = scalefactors_Medium_Muon                   ->GetBinError  (binXT,binYT_c);
          sfSystWeightLepEff[1][1] = scalefactors_Medium_Muon_stat_error_hi	->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[1][2] = scalefactors_Medium_Muon_signalFsrTNP	->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[1][3] = scalefactors_Medium_Muon_signalResTNP	->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[1][4] = scalefactors_Medium_Muon_bkgModelTNP	->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[1][5] = scalefactors_Medium_Muon_tagBiasTNP	->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[1][6] = scalefactors_Medium_Muon_generatorChoiceTNP->GetBinContent(binXT,binYT_s);
          //for(int nj=0; nj<nMuSFBins; nj++) sfSystWeightLepEff[1][nj+5] = scalefactors_Medium_Muon_stat_error_hi_bins[nj]->GetBinContent(binXT,binYT_s);
        } else {
          double etal = thePandaFlat.looseLep2Eta; if(etal >= 2.5) etal = 2.4999; else if(etal <= -2.5) etal = -2.4999;
          int binXT   = scalefactors_Medium_Electron->GetXaxis()->FindFixBin(etal);
          int binYT_c = scalefactors_Medium_Electron->GetYaxis()->FindFixBin(TMath::Min((double)thePandaFlat.looseLep2Pt,getMaxPtForSFs[1]));
          int binYT_s = scalefactors_Medium_Electron->GetYaxis()->FindFixBin(TMath::Min((double)thePandaFlat.looseLep2Pt,getMaxPtForSFs[3]));
	  sfWeightLepEff[1]        = scalefactors_Medium_Electron                   ->GetBinContent(binXT,binYT_c);
          sfSystWeightLepEff[1][0] = scalefactors_Medium_Electron                   ->GetBinError  (binXT,binYT_c);
          sfSystWeightLepEff[1][1] = scalefactors_Medium_Electron_stat_error_hi     ->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[1][2] = scalefactors_Medium_Electron_signalFsrTNP      ->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[1][3] = scalefactors_Medium_Electron_signalResTNP      ->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[1][4] = scalefactors_Medium_Electron_bkgModelTNP	    ->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[1][5] = scalefactors_Medium_Electron_tagBiasTNP	    ->GetBinContent(binXT,binYT_s);
          sfSystWeightLepEff[1][6] = scalefactors_Medium_Electron_generatorChoiceTNP->GetBinContent(binXT,binYT_s);
          //for(int nj=0; nj<nElSFBins; nj++) sfSystWeightLepEff[1][nj+5] = scalefactors_Medium_Electron_stat_error_hi_bins[nj]->GetBinContent(binXT,binYT_s);
        }

        /*
        printf("%2d %7.2f %7.2f %7.3f %7.3f %7.4f | %2d %7.2f %7.2f %7.3f %7.3f %7.4f - %7.4f/%7.4f %7.4f/%7.4f\n",abs(thePandaFlat.looseLep1PdgId),thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.sf_medium1,sfWeightLepEff[0],
	100*(thePandaFlat.sf_medium1-sfWeightLepEff[0])/sfWeightLepEff[0],
	                                                                     abs(thePandaFlat.looseLep2PdgId),thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.sf_medium2,sfWeightLepEff[1],
	100*(thePandaFlat.sf_medium2-sfWeightLepEff[1])/sfWeightLepEff[1],
	thePandaFlat.sf_unc1,sfSystWeightLepEff[0][0],thePandaFlat.sf_unc2,sfSystWeightLepEff[1][0]);

        if(ZRecRap>2.3) printf("%2d %7.2f %7.2f %7.4f / %2d %7.2f %7.2f %7.4f\n",abs(thePandaFlat.looseLep1PdgId),thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,sfSystWeightLepEff[0][0],
	                    	                                                 abs(thePandaFlat.looseLep2PdgId),thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,sfSystWeightLepEff[1][0]);
        */

        //totalWeight = thePandaFlat.normalizedWeight * lumi * thePandaFlat.sf_pu *
	//	      thePandaFlat.sf_trk1 * thePandaFlat.sf_medium1 *
	//	      thePandaFlat.sf_trk2 * thePandaFlat.sf_medium2 *
	//	      trigger_sf(thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1PdgId,thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2PdgId);
        totalWeight = thePandaFlat.normalizedWeight * lumi * thePandaFlat.sf_pu *
		      the_eta_sf[0] * sfWeightLepEff[0] *
		      the_eta_sf[1] * sfWeightLepEff[1] *
		      trigger_sf(thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1PdgId,thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2PdgId);
	/*
	double totalSumWeightLepEff[2] = {0,0};
	for(int nl=0; nl<2; nl++) {
  	  for(int ns=1; ns<nEffNuisances; ns++) {
	    totalSumWeightLepEff[nl] = totalSumWeightLepEff[nl] + sfSystWeightLepEff[nl][ns]*sfSystWeightLepEff[nl][ns];
	  }
	  printf("lep(%d) %7.4f-%7.4f=%7.4f ",nl,sfSystWeightLepEff[nl][0],sqrt(totalSumWeightLepEff[nl]),sfSystWeightLepEff[nl][0]-sqrt(totalSumWeightLepEff[nl]));
	}
        printf("\n");
	*/
      }

      if(passSel){
        double the_cos_theta_star = cos_theta_star(v1,v2,(v1+v2));
        if(theCategory != 5){
	  histo[lepType+0][theCategory]->Fill((v1+v2).M(),totalWeight);
	  histo[lepType+2][theCategory]->Fill(ZRecPt,totalWeight);
	  histo[lepType+4][theCategory]->Fill(TMath::Abs(ZRecPt-ZGenPt),totalWeight);
          if     (TMath::Abs(thePandaFlat.looseLep1Eta) <  1.5 && TMath::Abs(thePandaFlat.looseLep2Eta) <  1.5){
	    histo[lepType+ 6][theCategory]->Fill((v1+v2).M(),totalWeight);
          }
          else if(TMath::Abs(thePandaFlat.looseLep1Eta) >= 1.5 && TMath::Abs(thePandaFlat.looseLep2Eta) >= 1.5){
	    histo[lepType+10][theCategory]->Fill((v1+v2).M(),totalWeight);
          }
          else {
	    histo[lepType+ 8][theCategory]->Fill((v1+v2).M(),totalWeight);
          }
	  histo[lepType+12][theCategory]->Fill((v1+v2).Px(),totalWeight);
	  histo[lepType+14][theCategory]->Fill((v1+v2).Py(),totalWeight);
	  histo[lepType+16][theCategory]->Fill((v1+v2).Rapidity(),totalWeight);
	  histo[lepType+18][theCategory]->Fill(the_cos_theta_star,totalWeight);
	  histo[lepType+20][theCategory]->Fill(ZRecPhiStar,totalWeight);

	  histo[lepType+30][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	  if     ((v1+v2).Rapidity() < -2.2) histo[lepType+40][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	  else if((v1+v2).Rapidity() < -2.0) histo[lepType+50][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	  if     (v1.Pt() > 70){
	    histo[lepType+32][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	    if     ((v1+v2).Rapidity() < -2.2) histo[lepType+42][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	    else if((v1+v2).Rapidity() < -2.0) histo[lepType+52][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	  }
	  else if(v1.Pt() > 45){
	    histo[lepType+34][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	    if     ((v1+v2).Rapidity() < -2.2) histo[lepType+44][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	    else if((v1+v2).Rapidity() < -2.0) histo[lepType+54][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	  }
	  else if(v1.Pt() > 35){
	    histo[lepType+36][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	    if     ((v1+v2).Rapidity() < -2.2) histo[lepType+46][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	    else if((v1+v2).Rapidity() < -2.0) histo[lepType+56][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	  }
	  else {
	    histo[lepType+38][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	    if     ((v1+v2).Rapidity() < -2.2) histo[lepType+48][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	    else if((v1+v2).Rapidity() < -2.0) histo[lepType+58][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
	  }

	  histo[lepType+30][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	  if     ((v1+v2).Rapidity() < -2.2) histo[lepType+40][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	  else if((v1+v2).Rapidity() < -2.0) histo[lepType+50][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	  if     (v2.Pt() > 70){
	    histo[lepType+32][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	    if     ((v1+v2).Rapidity() < -2.2) histo[lepType+42][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	    else if((v1+v2).Rapidity() < -2.0) histo[lepType+52][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	  }
	  else if(v2.Pt() > 45){
	    histo[lepType+34][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	    if     ((v1+v2).Rapidity() < -2.2) histo[lepType+44][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	    else if((v1+v2).Rapidity() < -2.0) histo[lepType+54][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	  }
	  else if(v2.Pt() > 35){
	    histo[lepType+36][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	    if     ((v1+v2).Rapidity() < -2.2) histo[lepType+46][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	    else if((v1+v2).Rapidity() < -2.0) histo[lepType+56][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	  }
	  else {
	    histo[lepType+38][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	    if     ((v1+v2).Rapidity() < -2.2) histo[lepType+48][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	    else if((v1+v2).Rapidity() < -2.0) histo[lepType+58][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
	  }
	  histo[lepType+60][theCategory]->Fill(ZRecRap,totalWeight);

        }
        else {
          for(int theLepType = 0; theLepType<2; theLepType++) {
            double theKeff = k_eff;
            if(theLepType == 1) theKeff = 1.0/k_eff;
 	    histo[theLepType+0][theCategory]->Fill((v1+v2).M(),totalWeight*theKeff);
	    histo[theLepType+2][theCategory]->Fill(ZRecPt,totalWeight*theKeff);
	    histo[theLepType+4][theCategory]->Fill(TMath::Abs(ZRecPt-ZGenPt),totalWeight*theKeff);
            if     (TMath::Abs(thePandaFlat.looseLep1Eta) <  1.5 && TMath::Abs(thePandaFlat.looseLep2Eta) <  1.5){
	      histo[theLepType+ 6][theCategory]->Fill((v1+v2).M(),totalWeight*theKeff);
            }
            else if(TMath::Abs(thePandaFlat.looseLep1Eta) >= 1.5 && TMath::Abs(thePandaFlat.looseLep2Eta) >= 1.5){
	      histo[theLepType+10][theCategory]->Fill((v1+v2).M(),totalWeight*theKeff);
            }
            else {
	      histo[theLepType+ 8][theCategory]->Fill((v1+v2).M(),totalWeight*theKeff);
            }
	    histo[theLepType+12][theCategory]->Fill((v1+v2).Px(),totalWeight*theKeff);
	    histo[theLepType+14][theCategory]->Fill((v1+v2).Py(),totalWeight*theKeff);
	    histo[theLepType+16][theCategory]->Fill((v1+v2).Rapidity(),totalWeight*theKeff);
	    histo[theLepType+18][theCategory]->Fill(the_cos_theta_star,totalWeight*theKeff);
	    histo[theLepType+20][theCategory]->Fill(ZRecPhiStar,totalWeight*theKeff);

	    histo[theLepType+30][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	    if     ((v1+v2).Rapidity() < -2.2) histo[theLepType+40][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	    else if((v1+v2).Rapidity() < -2.0) histo[theLepType+50][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	    if     (v1.Pt() > 70){
	      histo[theLepType+32][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	      if     ((v1+v2).Rapidity() < -2.2) histo[theLepType+42][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	      else if((v1+v2).Rapidity() < -2.0) histo[theLepType+52][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	    }
	    else if(v1.Pt() > 45){
	      histo[theLepType+34][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	      if     ((v1+v2).Rapidity() < -2.2) histo[theLepType+44][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	      else if((v1+v2).Rapidity() < -2.0) histo[theLepType+54][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	    }
	    else if(v1.Pt() > 35){
	      histo[theLepType+36][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	      if     ((v1+v2).Rapidity() < -2.2) histo[theLepType+46][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	      else if((v1+v2).Rapidity() < -2.0) histo[theLepType+56][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	    }
	    else {
	      histo[theLepType+38][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	      if     ((v1+v2).Rapidity() < -2.2) histo[theLepType+48][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	      else if((v1+v2).Rapidity() < -2.0) histo[theLepType+58][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight*theKeff);
	    }

	    histo[theLepType+30][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	    if     ((v1+v2).Rapidity() < -2.2) histo[theLepType+40][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	    else if((v1+v2).Rapidity() < -2.0) histo[theLepType+50][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	    if     (v2.Pt() > 70){
	      histo[theLepType+32][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	      if     ((v1+v2).Rapidity() < -2.2) histo[theLepType+42][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	      else if((v1+v2).Rapidity() < -2.0) histo[theLepType+52][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	    }
	    else if(v2.Pt() > 45){
	      histo[theLepType+34][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	      if     ((v1+v2).Rapidity() < -2.2) histo[theLepType+44][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	      else if((v1+v2).Rapidity() < -2.0) histo[theLepType+54][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	    }
	    else if(v2.Pt() > 35){
	      histo[theLepType+36][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	      if     ((v1+v2).Rapidity() < -2.2) histo[theLepType+46][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	      else if((v1+v2).Rapidity() < -2.0) histo[theLepType+56][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	    }
	    else {
	      histo[theLepType+38][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	      if     ((v1+v2).Rapidity() < -2.2) histo[theLepType+48][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	      else if((v1+v2).Rapidity() < -2.0) histo[theLepType+58][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight*theKeff);
	    }
            histo[theLepType+60][theCategory]->Fill(ZRecRap,totalWeight*theKeff);
          }
        }

        double maxQCDscale = 1.0;
        if(theCategory != 0 && theCategory != 5 && thePandaFlat.scale[0] != -1){
          maxQCDscale = (TMath::Abs(1+thePandaFlat.scale[0])+TMath::Abs(1+thePandaFlat.scale[1])+TMath::Abs(1+thePandaFlat.scale[2])+
        		 TMath::Abs(1+thePandaFlat.scale[3])+TMath::Abs(1+thePandaFlat.scale[4])+TMath::Abs(1+thePandaFlat.scale[4]))/6.0;
        }

	if     (theCategory == 2){ // DY
          // Tot
          histoTotRecDY        [lepType]   ->Fill(ZRecTot,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoTotRecDY_RecEff[lepType][nj]->Fill(ZRecTot,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoTotRecDY_LepEff[lepType][nj]->Fill(ZRecTot,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
          histoTotRecDY_PDF    [lepType]   ->Fill(ZRecTot,totalWeight*thePandaFlat.pdfUp);
          histoTotRecDY_QCDPart[lepType][0]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoTotRecDY_QCDPart[lepType][1]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoTotRecDY_QCDPart[lepType][2]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoTotRecDY_QCDPart[lepType][3]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoTotRecDY_QCDPart[lepType][4]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoTotRecDY_QCDPart[lepType][5]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          // Pt
          histoPtRecDY        [lepType]   ->Fill(ZRecPt,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPtRecDY_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRecDY_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
          histoPtRecDY_PDF    [lepType]   ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPtRecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          // PhiStar
          histoPhiStarRecDY        [lepType]   ->Fill(ZRecPhiStar,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPhiStarRecDY_RecEff[lepType][nj]->Fill(ZRecPhiStar,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPhiStarRecDY_LepEff[lepType][nj]->Fill(ZRecPhiStar,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
          histoPhiStarRecDY_PDF    [lepType]   ->Fill(ZRecPhiStar,totalWeight*thePandaFlat.pdfUp);
          histoPhiStarRecDY_QCDPart[lepType][0]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPhiStarRecDY_QCDPart[lepType][1]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPhiStarRecDY_QCDPart[lepType][2]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPhiStarRecDY_QCDPart[lepType][3]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPhiStarRecDY_QCDPart[lepType][4]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPhiStarRecDY_QCDPart[lepType][5]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          // Rap
          if(ZRecRap < 2.4){
            histoRapRecDY        [lepType]   ->Fill(ZRecRap,totalWeight);
            for(int nj=0; nj<nRecNuisances; nj++) histoRapRecDY_RecEff[lepType][nj]->Fill(ZRecRap,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoRapRecDY_LepEff[lepType][nj]->Fill(ZRecRap,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
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
            for(int nj=0; nj<nRecNuisances; nj++) histoPtRap0RecDY_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap0RecDY_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
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
            for(int nj=0; nj<nRecNuisances; nj++) histoPtRap1RecDY_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap1RecDY_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
            histoPtRap1RecDY_PDF    [lepType]   ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
            histoPtRap1RecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRap1RecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRap1RecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRap1RecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRap1RecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRap1RecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
          else if(ZRecRap < 1.5){
            histoPtRap2RecDY	    [lepType]	->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nRecNuisances; nj++) histoPtRap2RecDY_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap2RecDY_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
            histoPtRap2RecDY_PDF    [lepType]	->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
            histoPtRap2RecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRap2RecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRap2RecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRap2RecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRap2RecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRap2RecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
          //else if(ZRecRap < 2.0){
          else if(ZRecRap < 2.4){
            histoPtRap3RecDY	    [lepType]	->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nRecNuisances; nj++) histoPtRap3RecDY_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap3RecDY_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
            histoPtRap3RecDY_PDF    [lepType]	->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
            histoPtRap3RecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRap3RecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRap3RecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRap3RecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRap3RecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRap3RecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
          else if(ZRecRap < 2.4){
            histoPtRap4RecDY	    [lepType]	->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nRecNuisances; nj++) histoPtRap4RecDY_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap4RecDY_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
            histoPtRap4RecDY_PDF    [lepType]	->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
            histoPtRap4RecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRap4RecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRap4RecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRap4RecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRap4RecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRap4RecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
	}
	else if(theCategory == 0){ // Data
          // Tot
          histoTotRecDA[lepType]->Fill(ZRecTot,totalWeight);
          // Pt
          histoPtRecDA[lepType]->Fill(ZRecPt,totalWeight);
          // PhiStar
          histoPhiStarRecDA[lepType]->Fill(ZRecPhiStar,totalWeight);
          // Rap
          if(ZRecRap < 2.4){
          histoRapRecDA[lepType]->Fill(ZRecRap,totalWeight);
          }
          // PtRap0
          if     (ZRecRap < 0.5){
            histoPtRap0RecDA[lepType]->Fill(ZRecPt,totalWeight);
          }
          else if(ZRecRap < 1.0){
            histoPtRap1RecDA[lepType]->Fill(ZRecPt,totalWeight);
          }
          else if(ZRecRap < 1.5){
            histoPtRap2RecDA[lepType]->Fill(ZRecPt,totalWeight);
          }
          //else if(ZRecRap < 2.0){
          else if(ZRecRap < 2.4){
            histoPtRap3RecDA[lepType]->Fill(ZRecPt,totalWeight);
          }
          else if(ZRecRap < 2.4){
            histoPtRap4RecDA[lepType]->Fill(ZRecPt,totalWeight);
          }
	}
	else if(theCategory == 5){ // e-mu Data
          for(int theLepType = 0; theLepType<2; theLepType++) {
            double theKeff = k_eff;
            if(theLepType == 1) theKeff = 1.0/k_eff;
            // Tot
            histoTotRecEM[theLepType]->Fill(ZRecTot,totalWeight*theKeff);
            // Pt
            histoPtRecEM[theLepType]->Fill(ZRecPt,totalWeight*theKeff);
            // PhiStar
            histoPhiStarRecEM[theLepType]->Fill(ZRecPhiStar,totalWeight*theKeff);
            // Rap
            if(ZRecRap < 2.4){
              histoRapRecEM[theLepType]->Fill(ZRecRap,totalWeight*theKeff);
            }
            // PtRap0
            if     (ZRecRap < 0.5){
              histoPtRap0RecEM[theLepType]->Fill(ZRecPt,totalWeight*theKeff);
            }
            else if(ZRecRap < 1.0){
              histoPtRap1RecEM[theLepType]->Fill(ZRecPt,totalWeight*theKeff);
            }
            else if(ZRecRap < 1.5){
              histoPtRap2RecEM[theLepType]->Fill(ZRecPt,totalWeight*theKeff);
            }
            //else if(ZRecRap < 2.0){
            else if(ZRecRap < 2.4){
              histoPtRap3RecEM[theLepType]->Fill(ZRecPt,totalWeight*theKeff);
            }
            else if(ZRecRap < 2.4){
              histoPtRap4RecEM[theLepType]->Fill(ZRecPt,totalWeight*theKeff);
            }
          } // loop over muons and electrons
	}
	else { // VV
          // Tot
          histoTotRecVV        [lepType]->Fill(ZRecTot,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoTotRecVV_RecEff[lepType][nj]->Fill(ZRecTot,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoTotRecVV_LepEff[lepType][nj]->Fill(ZRecTot,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
          histoTotRecVV_PDF    [lepType]->Fill(ZRecTot,totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoTotRecVV_QCDPart[lepType][0]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoTotRecVV_QCDPart[lepType][1]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoTotRecVV_QCDPart[lepType][2]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoTotRecVV_QCDPart[lepType][3]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoTotRecVV_QCDPart[lepType][4]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoTotRecVV_QCDPart[lepType][5]->Fill(ZRecTot,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoTotRecVV_QCDPart[lepType][0]->Fill(ZRecTot,totalWeight);
            histoTotRecVV_QCDPart[lepType][1]->Fill(ZRecTot,totalWeight);
            histoTotRecVV_QCDPart[lepType][2]->Fill(ZRecTot,totalWeight);
            histoTotRecVV_QCDPart[lepType][3]->Fill(ZRecTot,totalWeight);
            histoTotRecVV_QCDPart[lepType][4]->Fill(ZRecTot,totalWeight);
            histoTotRecVV_QCDPart[lepType][5]->Fill(ZRecTot,totalWeight);
          }
          // Pt
          histoPtRecVV        [lepType]->Fill(ZRecPt,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPtRecVV_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRecVV_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
          histoPtRecVV_PDF    [lepType]->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoPtRecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoPtRecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
            histoPtRecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
            histoPtRecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
            histoPtRecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
            histoPtRecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
            histoPtRecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
          }
          // PhiStar
          histoPhiStarRecVV        [lepType]->Fill(ZRecPhiStar,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPhiStarRecVV_RecEff[lepType][nj]->Fill(ZRecPhiStar,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPhiStarRecVV_LepEff[lepType][nj]->Fill(ZRecPhiStar,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
          histoPhiStarRecVV_PDF    [lepType]->Fill(ZRecPhiStar,totalWeight*thePandaFlat.pdfUp);
          if(thePandaFlat.scale[0] != -1){
            histoPhiStarRecVV_QCDPart[lepType][0]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPhiStarRecVV_QCDPart[lepType][1]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPhiStarRecVV_QCDPart[lepType][2]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPhiStarRecVV_QCDPart[lepType][3]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPhiStarRecVV_QCDPart[lepType][4]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPhiStarRecVV_QCDPart[lepType][5]->Fill(ZRecPhiStar,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          } else {
            histoPhiStarRecVV_QCDPart[lepType][0]->Fill(ZRecPhiStar,totalWeight);
            histoPhiStarRecVV_QCDPart[lepType][1]->Fill(ZRecPhiStar,totalWeight);
            histoPhiStarRecVV_QCDPart[lepType][2]->Fill(ZRecPhiStar,totalWeight);
            histoPhiStarRecVV_QCDPart[lepType][3]->Fill(ZRecPhiStar,totalWeight);
            histoPhiStarRecVV_QCDPart[lepType][4]->Fill(ZRecPhiStar,totalWeight);
            histoPhiStarRecVV_QCDPart[lepType][5]->Fill(ZRecPhiStar,totalWeight);
          }
          // Rap
          if(ZRecRap < 2.4){
            histoRapRecVV        [lepType]     ->Fill(ZRecRap,totalWeight);
            for(int nj=0; nj<nRecNuisances; nj++) histoRapRecVV_RecEff[lepType][nj]->Fill(ZRecRap,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoRapRecVV_LepEff[lepType][nj]->Fill(ZRecRap,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
            histoRapRecVV_PDF    [lepType]     ->Fill(ZRecRap,totalWeight*thePandaFlat.pdfUp);
            if(thePandaFlat.scale[0] != -1){
              histoRapRecVV_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
              histoRapRecVV_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
              histoRapRecVV_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
              histoRapRecVV_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
              histoRapRecVV_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
              histoRapRecVV_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
            } else {
              histoRapRecVV_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight);
              histoRapRecVV_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight);
              histoRapRecVV_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight);
              histoRapRecVV_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight);
              histoRapRecVV_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight);
              histoRapRecVV_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight);
            }
          }
          // PtRap0
          if     (ZRecRap < 0.5){
            histoPtRap0RecVV        [lepType]     ->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nRecNuisances; nj++) histoPtRap0RecVV_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap0RecVV_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
            histoPtRap0RecVV_PDF    [lepType]     ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
            if(thePandaFlat.scale[0] != -1){
              histoPtRap0RecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
              histoPtRap0RecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
              histoPtRap0RecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
              histoPtRap0RecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
              histoPtRap0RecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
              histoPtRap0RecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
            } else {
              histoPtRap0RecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
              histoPtRap0RecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
              histoPtRap0RecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
              histoPtRap0RecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
              histoPtRap0RecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
              histoPtRap0RecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
            }
          }
          else if(ZRecRap < 1.0){
            histoPtRap1RecVV        [lepType]     ->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nRecNuisances; nj++) histoPtRap1RecVV_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap1RecVV_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
            histoPtRap1RecVV_PDF    [lepType]     ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
            if(thePandaFlat.scale[0] != -1){
              histoPtRap1RecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
              histoPtRap1RecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
              histoPtRap1RecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
              histoPtRap1RecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
              histoPtRap1RecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
              histoPtRap1RecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
            } else {
              histoPtRap1RecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
              histoPtRap1RecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
              histoPtRap1RecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
              histoPtRap1RecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
              histoPtRap1RecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
              histoPtRap1RecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
            }
          }
          else if(ZRecRap < 1.5){
            histoPtRap2RecVV        [lepType]     ->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nRecNuisances; nj++) histoPtRap2RecVV_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap2RecVV_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
            histoPtRap2RecVV_PDF    [lepType]     ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
            if(thePandaFlat.scale[0] != -1){
              histoPtRap2RecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
              histoPtRap2RecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
              histoPtRap2RecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
              histoPtRap2RecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
              histoPtRap2RecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
              histoPtRap2RecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
            } else {
              histoPtRap2RecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
              histoPtRap2RecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
              histoPtRap2RecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
              histoPtRap2RecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
              histoPtRap2RecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
              histoPtRap2RecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
            }
          }
          //else if(ZRecRap < 2.0){
          else if(ZRecRap < 2.4){
            histoPtRap3RecVV        [lepType]     ->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nRecNuisances; nj++) histoPtRap3RecVV_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap3RecVV_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
            histoPtRap3RecVV_PDF    [lepType]     ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
            if(thePandaFlat.scale[0] != -1){
              histoPtRap3RecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
              histoPtRap3RecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
              histoPtRap3RecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
              histoPtRap3RecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
              histoPtRap3RecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
              histoPtRap3RecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
            } else {
              histoPtRap3RecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
              histoPtRap3RecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
              histoPtRap3RecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
              histoPtRap3RecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
              histoPtRap3RecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
              histoPtRap3RecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
            }
          }
          else if(ZRecRap < 2.4){
            histoPtRap4RecVV        [lepType]     ->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nRecNuisances; nj++) histoPtRap4RecVV_RecEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap4RecVV_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
            histoPtRap4RecVV_PDF    [lepType]     ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
            if(thePandaFlat.scale[0] != -1){
              histoPtRap4RecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
              histoPtRap4RecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
              histoPtRap4RecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
              histoPtRap4RecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
              histoPtRap4RecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
              histoPtRap4RecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
            } else {
              histoPtRap4RecVV_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight);
              histoPtRap4RecVV_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight);
              histoPtRap4RecVV_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight);
              histoPtRap4RecVV_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight);
              histoPtRap4RecVV_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight);
              histoPtRap4RecVV_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight);
            }
          }
	}

	if(theCategory == 2 && passPtFid == true){
          histoTotRecGen       [lepType] ->Fill(ZRecTot,ZGenTot,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoTotRecGen_RecEff [lepType][nj]->Fill(ZRecTot,ZGenTot,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoTotRecGen_LepEff [lepType][nj]->Fill(ZRecTot,ZGenTot,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));

          histoPtRecGen       [lepType] ->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPtRecGen_RecEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));

          histoPhiStarRecGen       [lepType] ->Fill(ZRecPhiStar,ZGenPhiStar,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPhiStarRecGen_RecEff [lepType][nj]->Fill(ZRecPhiStar,ZGenPhiStar,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPhiStarRecGen_LepEff [lepType][nj]->Fill(ZRecPhiStar,ZGenPhiStar,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
	}
	if(theCategory == 2 && passRapFid == true && ZRecRap < 2.4){
          histoRapRecGen       [lepType]->Fill(ZRecRap,ZGenRap,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoRapRecGen_RecEff [lepType][nj]->Fill(ZRecRap,ZGenRap,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoRapRecGen_LepEff [lepType][nj]->Fill(ZRecRap,ZGenRap,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
	}
	if(theCategory == 2 && passPtRapFid[0] == true && ZRecRap >= 0.0 && ZRecRap < 0.5){
          histoPtRap0RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPtRap0RecGen_RecEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRap0RecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
	}
	if(theCategory == 2 && passPtRapFid[1] == true && ZRecRap >= 0.5 && ZRecRap < 1.0){
          histoPtRap1RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPtRap1RecGen_RecEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRap1RecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
	}
	if(theCategory == 2 && passPtRapFid[2] == true && ZRecRap >= 1.0 && ZRecRap < 1.5){
          histoPtRap2RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPtRap2RecGen_RecEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRap2RecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
	}
	//if(theCategory == 2 && passPtRapFid[3] == true && ZRecRap >= 1.5 && ZRecRap < 2.0){
	if(theCategory == 2 && passPtRapFid[3] == true && ZRecRap >= 1.5 && ZRecRap < 2.4){
          histoPtRap3RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPtRap3RecGen_RecEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRap3RecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
	}
	if(theCategory == 2 && passPtRapFid[4] == true && ZRecRap >= 2.0 && ZRecRap < 2.4){
          histoPtRap4RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nRecNuisances; nj++) histoPtRap4RecGen_RecEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+the_eta_sf_unc[0])*(1+the_eta_sf_unc[1]));
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRap4RecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+sfSystWeightLepEff[0][nj])*(1+sfSystWeightLepEff[1][nj]));
	}
      } // end passSel
      
      for(int nc=0; nc<nMomNuisances; nc++){
	if(passSystSel[nc]){
	  if     (theCategory == 2){
            histoTotRecDY_MomRes[lepType][nc] ->Fill(ZRecSystTot[nc],totalWeight);
            if(passPtFid == true) histoTotRecGen_MomRes[lepType][nc]->Fill(ZRecSystTot[nc],ZGenTot,totalWeight);

            histoPtRecDY_MomRes[lepType][nc] ->Fill(ZRecSystPt[nc],totalWeight);
            if(passPtFid == true) histoPtRecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);

            histoPhiStarRecDY_MomRes[lepType][nc] ->Fill(ZRecPhiStar,totalWeight);
            if(passPtFid == true) histoPhiStarRecGen_MomRes[lepType][nc]->Fill(ZRecPhiStar,ZGenPhiStar,totalWeight);

            if(ZRecRap < 2.4) {
              histoRapRecDY_MomRes[lepType][nc] ->Fill(ZRecRap,totalWeight);
              if(passRapFid == true) histoRapRecGen_MomRes[lepType][nc]->Fill(ZRecRap,ZGenRap,totalWeight);
            }
            if     (ZRecRap < 0.5) {if(passPtRapFid[0] == true) histoPtRap0RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap0RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
            else if(ZRecRap < 1.0) {if(passPtRapFid[1] == true) histoPtRap1RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap1RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
            else if(ZRecRap < 1.5) {if(passPtRapFid[2] == true) histoPtRap2RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap2RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
            //else if(ZRecRap < 2.0) {if(passPtRapFid[3] == true) histoPtRap3RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap3RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
            else if(ZRecRap < 2.4) {if(passPtRapFid[3] == true) histoPtRap3RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap3RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
            else if(ZRecRap < 2.4) {if(passPtRapFid[4] == true) histoPtRap4RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap4RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
	  }
	  else if(theCategory == 0){
            histoTotRecDA_MomRes[lepType][nc] ->Fill(ZRecSystTot[nc],totalWeight);
            histoPtRecDA_MomRes[lepType][nc] ->Fill(ZRecSystPt[nc],totalWeight);
            histoPhiStarRecDA_MomRes[lepType][nc] ->Fill(ZRecPhiStar,totalWeight);
            if(ZRecRap < 2.4) {
              histoRapRecDA_MomRes[lepType][nc]->Fill(ZRecRap,totalWeight);
            }
            if     (ZRecRap < 0.5) histoPtRap0RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 1.0) histoPtRap1RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 1.5) histoPtRap2RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            //else if(ZRecRap < 2.0) histoPtRap3RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 2.4) histoPtRap3RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 2.4) histoPtRap4RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
	  }
	  else if(theCategory == 5){
            for(int theLepType = 0; theLepType<2; theLepType++) {
              double theKeff = k_eff;
              if(theLepType == 1) theKeff = 1.0/k_eff;
              histoTotRecEM_MomRes[theLepType][nc] ->Fill(ZRecSystTot[nc],totalWeight*theKeff);
              histoPtRecEM_MomRes[theLepType][nc] ->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              histoPhiStarRecEM_MomRes[theLepType][nc] ->Fill(ZRecPhiStar,totalWeight*theKeff);
              if(ZRecRap < 2.4) {
        	histoRapRecEM_MomRes[theLepType][nc]->Fill(ZRecRap,totalWeight*theKeff);
              }
              if     (ZRecRap < 0.5) histoPtRap0RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              else if(ZRecRap < 1.0) histoPtRap1RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              else if(ZRecRap < 1.5) histoPtRap2RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              //else if(ZRecRap < 2.0) histoPtRap3RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              else if(ZRecRap < 2.4) histoPtRap3RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              else if(ZRecRap < 2.4) histoPtRap4RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
            }
	  }
	  else {
            histoTotRecVV_MomRes[lepType][nc] ->Fill(ZRecSystTot[nc],totalWeight);
            histoPtRecVV_MomRes[lepType][nc] ->Fill(ZRecSystPt[nc],totalWeight);
            histoPhiStarRecVV_MomRes[lepType][nc] ->Fill(ZRecPhiStar,totalWeight);
            if(ZRecRap < 2.4) {
              histoRapRecVV_MomRes[lepType][nc]->Fill(ZRecRap,totalWeight);
            }
            if     (ZRecRap < 0.5) histoPtRap0RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 1.0) histoPtRap1RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 1.5) histoPtRap2RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            //else if(ZRecRap < 2.0) histoPtRap3RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 2.4) histoPtRap3RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 2.4) histoPtRap4RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
	  }
	} // end passSystSel[nc]
      } // end loop over passSystSel
    } // end event loop
  } // end samples loop

  // Tot
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDTot(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoTotRecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoTotRecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoTotRecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoTotRecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoTotRecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoTotRecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoTotRecVV[ntype]->GetSumOfWeights(),
           histoTotRecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoTotRecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoTotRecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoTotRecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoTotRecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoTotRecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoTotRecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoTot+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoTotRecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoTotRecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoTotRecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoTotRecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoTotRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoTotRecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoTotRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoTotRecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoTotRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoTotRecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoTotRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoTotRecDY[ntype]->GetBinContent(nb));
      }
      if(histoTotRecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoTotRecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoTotRecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoTotRecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoTotRecVV_QCD[ntype]->SetBinContent(nb, histoTotRecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoTotRecDY_QCD[ntype]->SetBinContent(nb, histoTotRecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // Pt
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPt(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRecVV[ntype]->GetSumOfWeights(),
           histoPtRecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPt+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRecVV_QCD[ntype]->SetBinContent(nb, histoPtRecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRecDY_QCD[ntype]->SetBinContent(nb, histoPtRecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PhiStar
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPhiStar(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPhiStarRecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoPhiStarRecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoPhiStarRecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPhiStarRecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoPhiStarRecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoPhiStarRecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoPhiStarRecVV[ntype]->GetSumOfWeights(),
           histoPhiStarRecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPhiStarRecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPhiStarRecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPhiStarRecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPhiStarRecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPhiStarRecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPhiStarRecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPhiStar+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPhiStarRecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoPhiStarRecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPhiStarRecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPhiStarRecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPhiStarRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPhiStarRecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPhiStarRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPhiStarRecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPhiStarRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPhiStarRecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPhiStarRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPhiStarRecDY[ntype]->GetBinContent(nb));
      }
      if(histoPhiStarRecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPhiStarRecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPhiStarRecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPhiStarRecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPhiStarRecVV_QCD[ntype]->SetBinContent(nb, histoPhiStarRecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPhiStarRecDY_QCD[ntype]->SetBinContent(nb, histoPhiStarRecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // Rap
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDRap(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoRapRecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoRapRecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoRapRecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoRapRecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoRapRecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoRapRecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoRapRecVV[ntype]->GetSumOfWeights(),
           histoRapRecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoRapRecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoRapRecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoRapRecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoRapRecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoRapRecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoRapRecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRap+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoRapRecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoRapRecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoRapRecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoRapRecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoRapRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapRecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoRapRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapRecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoRapRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapRecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoRapRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapRecDY[ntype]->GetBinContent(nb));
      }
      if(histoRapRecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoRapRecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoRapRecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoRapRecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoRapRecVV_QCD[ntype]->SetBinContent(nb, histoRapRecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoRapRecDY_QCD[ntype]->SetBinContent(nb, histoRapRecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PtRap0
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPtRap0(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRap0RecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap0RecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap0RecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap0RecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap0RecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap0RecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap0RecVV[ntype]->GetSumOfWeights(),
           histoPtRap0RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap0RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap0RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap0RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap0RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap0RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap0RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPtRap0+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRap0RecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap0RecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRap0RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap0RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRap0RecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap0RecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRap0RecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap0RecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRap0RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap0RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRap0RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap0RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRap0RecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRap0RecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRap0RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRap0RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRap0RecVV_QCD[ntype]->SetBinContent(nb, histoPtRap0RecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRap0RecDY_QCD[ntype]->SetBinContent(nb, histoPtRap0RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PtRap1
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPtRap1(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRap1RecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap1RecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap1RecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap1RecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap1RecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap1RecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap1RecVV[ntype]->GetSumOfWeights(),
           histoPtRap1RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap1RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap1RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap1RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap1RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap1RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap1RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPtRap1+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRap1RecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap1RecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRap1RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap1RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRap1RecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap1RecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRap1RecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap1RecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRap1RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap1RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRap1RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap1RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRap1RecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRap1RecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRap1RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRap1RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRap1RecVV_QCD[ntype]->SetBinContent(nb, histoPtRap1RecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRap1RecDY_QCD[ntype]->SetBinContent(nb, histoPtRap1RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PtRap2
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPtRap2(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRap2RecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap2RecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap2RecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap2RecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap2RecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap2RecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap2RecVV[ntype]->GetSumOfWeights(),
           histoPtRap2RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap2RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap2RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap2RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap2RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap2RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap2RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPtRap2+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRap2RecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap2RecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRap2RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap2RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRap2RecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap2RecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRap2RecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap2RecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRap2RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap2RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRap2RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap2RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRap2RecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRap2RecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRap2RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRap2RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRap2RecVV_QCD[ntype]->SetBinContent(nb, histoPtRap2RecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRap2RecDY_QCD[ntype]->SetBinContent(nb, histoPtRap2RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PtRap3
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPtRap3(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRap3RecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap3RecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap3RecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap3RecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap3RecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap3RecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap3RecVV[ntype]->GetSumOfWeights(),
           histoPtRap3RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap3RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap3RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap3RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap3RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap3RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap3RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPtRap3+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRap3RecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap3RecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRap3RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap3RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRap3RecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap3RecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRap3RecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap3RecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRap3RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap3RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRap3RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap3RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRap3RecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRap3RecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRap3RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRap3RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRap3RecVV_QCD[ntype]->SetBinContent(nb, histoPtRap3RecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRap3RecDY_QCD[ntype]->SetBinContent(nb, histoPtRap3RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // PtRap4
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDPtRap4(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoPtRap4RecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap4RecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap4RecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap4RecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap4RecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap4RecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap4RecVV[ntype]->GetSumOfWeights(),
           histoPtRap4RecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoPtRap4RecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoPtRap4RecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoPtRap4RecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoPtRap4RecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoPtRap4RecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoPtRap4RecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRecoPtRap4+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoPtRap4RecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap4RecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoPtRap4RecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoPtRap4RecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoPtRap4RecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap4RecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoPtRap4RecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap4RecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoPtRap4RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap4RecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoPtRap4RecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoPtRap4RecDY[ntype]->GetBinContent(nb));
      }
      if(histoPtRap4RecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoPtRap4RecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoPtRap4RecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoPtRap4RecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoPtRap4RecVV_QCD[ntype]->SetBinContent(nb, histoPtRap4RecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoPtRap4RecDY_QCD[ntype]->SetBinContent(nb, histoPtRap4RecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  char output[200];
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    sprintf(output,"histoDY%dzll_%d%s.root",whichDY,thePlot,isNoDYName.Data());	
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
  sprintf(output,"histoDY%dzllTotRecGen%s.root",whichDY,isNoDYName.Data()); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoTotRecGen[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoTotRecGen_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoTotRecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoTotRecGen_MomRes[i][nj]->Write();
    histoTotRecDA[i]->Write();
    histoTotRecDY[i]->Write();
    histoTotRecEM[i]->Write();
    histoTotRecVV[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoTotRecDY_RecEff[i][nj]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoTotRecVV_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoTotRecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoTotRecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoTotRecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoTotRecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoTotRecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoTotRecVV_MomRes[i][nj]->Write();
    histoTotRecVV_PDF[i]->Write();
    histoTotRecDY_PDF[i]->Write();
    histoTotRecVV_QCD[i]->Write();
    histoTotRecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRecGen%s.root",whichDY,isNoDYName.Data()); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRecGen[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRecGen_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRecGen_MomRes[i][nj]->Write();
    histoPtRecDA[i]->Write();
    histoPtRecDY[i]->Write();
    histoPtRecEM[i]->Write();
    histoPtRecVV[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRecDY_RecEff[i][nj]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRecVV_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRecVV_MomRes[i][nj]->Write();
    histoPtRecVV_PDF[i]->Write();
    histoPtRecDY_PDF[i]->Write();
    histoPtRecVV_QCD[i]->Write();
    histoPtRecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllRapRecGen%s.root",whichDY,isNoDYName.Data());
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoRapRecGen[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoRapRecGen_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoRapRecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapRecGen_MomRes[i][nj]->Write();
    histoRapRecDA[i]->Write();
    histoRapRecDY[i]->Write();
    histoRapRecEM[i]->Write();
    histoRapRecVV[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoRapRecDY_RecEff[i][nj]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoRapRecVV_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoRapRecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoRapRecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapRecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapRecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapRecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapRecVV_MomRes[i][nj]->Write();
    histoRapRecVV_PDF[i]->Write();
    histoRapRecDY_PDF[i]->Write();
    histoRapRecVV_QCD[i]->Write();
    histoRapRecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPhiStarRecGen%s.root",whichDY,isNoDYName.Data()); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPhiStarRecGen[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPhiStarRecGen_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPhiStarRecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPhiStarRecGen_MomRes[i][nj]->Write();
    histoPhiStarRecDA[i]->Write();
    histoPhiStarRecDY[i]->Write();
    histoPhiStarRecEM[i]->Write();
    histoPhiStarRecVV[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPhiStarRecDY_RecEff[i][nj]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPhiStarRecVV_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPhiStarRecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPhiStarRecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPhiStarRecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPhiStarRecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPhiStarRecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPhiStarRecVV_MomRes[i][nj]->Write();
    histoPhiStarRecVV_PDF[i]->Write();
    histoPhiStarRecDY_PDF[i]->Write();
    histoPhiStarRecVV_QCD[i]->Write();
    histoPhiStarRecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap0RecGen%s.root",whichDY,isNoDYName.Data()); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRap0RecGen[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap0RecGen_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap0RecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap0RecGen_MomRes[i][nj]->Write();
    histoPtRap0RecDA[i]->Write();
    histoPtRap0RecDY[i]->Write();
    histoPtRap0RecEM[i]->Write();
    histoPtRap0RecVV[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap0RecDY_RecEff[i][nj]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap0RecVV_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap0RecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap0RecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap0RecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap0RecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap0RecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap0RecVV_MomRes[i][nj]->Write();
    histoPtRap0RecVV_PDF[i]->Write();
    histoPtRap0RecDY_PDF[i]->Write();
    histoPtRap0RecVV_QCD[i]->Write();
    histoPtRap0RecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap1RecGen%s.root",whichDY,isNoDYName.Data()); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRap1RecGen[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap1RecGen_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap1RecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap1RecGen_MomRes[i][nj]->Write();
    histoPtRap1RecDA[i]->Write();
    histoPtRap1RecDY[i]->Write();
    histoPtRap1RecEM[i]->Write();
    histoPtRap1RecVV[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap1RecDY_RecEff[i][nj]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap1RecVV_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap1RecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap1RecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap1RecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap1RecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap1RecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap1RecVV_MomRes[i][nj]->Write();
    histoPtRap1RecVV_PDF[i]->Write();
    histoPtRap1RecDY_PDF[i]->Write();
    histoPtRap1RecVV_QCD[i]->Write();
    histoPtRap1RecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap2RecGen%s.root",whichDY,isNoDYName.Data()); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRap2RecGen[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap2RecGen_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap2RecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap2RecGen_MomRes[i][nj]->Write();
    histoPtRap2RecDA[i]->Write();
    histoPtRap2RecDY[i]->Write();
    histoPtRap2RecEM[i]->Write();
    histoPtRap2RecVV[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap2RecDY_RecEff[i][nj]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap2RecVV_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap2RecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap2RecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap2RecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap2RecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap2RecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap2RecVV_MomRes[i][nj]->Write();
    histoPtRap2RecVV_PDF[i]->Write();
    histoPtRap2RecDY_PDF[i]->Write();
    histoPtRap2RecVV_QCD[i]->Write();
    histoPtRap2RecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap3RecGen%s.root",whichDY,isNoDYName.Data()); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRap3RecGen[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap3RecGen_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap3RecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap3RecGen_MomRes[i][nj]->Write();
    histoPtRap3RecDA[i]->Write();
    histoPtRap3RecDY[i]->Write();
    histoPtRap3RecEM[i]->Write();
    histoPtRap3RecVV[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap3RecDY_RecEff[i][nj]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap3RecVV_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap3RecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap3RecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap3RecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap3RecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap3RecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap3RecVV_MomRes[i][nj]->Write();
    histoPtRap3RecVV_PDF[i]->Write();
    histoPtRap3RecDY_PDF[i]->Write();
    histoPtRap3RecVV_QCD[i]->Write();
    histoPtRap3RecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap4RecGen%s.root",whichDY,isNoDYName.Data()); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRap4RecGen[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap4RecGen_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap4RecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap4RecGen_MomRes[i][nj]->Write();
    histoPtRap4RecDA[i]->Write();
    histoPtRap4RecDY[i]->Write();
    histoPtRap4RecEM[i]->Write();
    histoPtRap4RecVV[i]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap4RecDY_RecEff[i][nj]->Write();
    for(int nj=0; nj<nRecNuisances; nj++) histoPtRap4RecVV_RecEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap4RecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap4RecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap4RecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap4RecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap4RecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap4RecVV_MomRes[i][nj]->Write();
    histoPtRap4RecVV_PDF[i]->Write();
    histoPtRap4RecDY_PDF[i]->Write();
    histoPtRap4RecVV_QCD[i]->Write();
    histoPtRap4RecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

}
