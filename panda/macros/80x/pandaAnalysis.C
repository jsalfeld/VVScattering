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
  double k_eff = 0.5 * sqrt(20285930./12446486.);
  //*******************************************************
  //Input Files
  //*******************************************************
  TString filesPath    = "/data/t3home000/ceballos/panda/v_004_0/";
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
  else if(whichDY == 3)
 {infileName_.push_back(Form("%sDYJetsToLL_Pt0To50.root",filesPath.Data()));   infileCat_.push_back(2);
  infileName_.push_back(Form("%sDYJetsToLL_Pt50To100.root",filesPath.Data())); infileCat_.push_back(2);
  infileName_.push_back(Form("%sDYJetsToLL_Pt100To250.root",filesPath.Data()));infileCat_.push_back(2);
  infileName_.push_back(Form("%sDYJetsToLL_Pt250To400.root",filesPath.Data()));infileCat_.push_back(2);
  infileName_.push_back(Form("%sDYJetsToLL_Pt400To650.root",filesPath.Data()));infileCat_.push_back(2);
  infileName_.push_back(Form("%sDYJetsToLL_Pt650ToInf.root",filesPath.Data()));infileCat_.push_back(2);
  }
  infileName_.push_back(Form("%sDYJetsToLL_M-10to50.root" ,filesPath.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sTT2L.root" ,filesPath.Data()));                infileCat_.push_back(3);
  //infileName_.push_back(Form("%sTW.root" ,filesPath.Data()));                  infileCat_.push_back(3);

  infileName_.push_back(Form("%sqqZZ.root" ,filesPath.Data()));                infileCat_.push_back(4);
  infileName_.push_back(Form("%sggZZ.root" ,filesPath.Data()));                infileCat_.push_back(4);
  infileName_.push_back(Form("%sWZ.root" ,filesPath.Data()));                  infileCat_.push_back(4);
  infileName_.push_back(Form("%sVVV.root" ,filesPath.Data()));                 infileCat_.push_back(4);
  infileName_.push_back(Form("%sTTV.root" ,filesPath.Data()));                 infileCat_.push_back(4);
  //infileName_.push_back(Form("%sWGstar.root" ,filesPath.Data()));              infileCat_.push_back(4);
  //infileName_.push_back(Form("%sVG.root" ,filesPath.Data()));                  infileCat_.push_back(6);
  //infileName_.push_back(Form("%sH125.root" ,filesPath.Data()));                infileCat_.push_back(7);
/*
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
*/
  const int nBinPt     = 37; Float_t xbinsPt[nBinPt+1]         = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,350,400,450,500,1000};

  const int nBinPtRap0 = 37; Float_t xbinsPtRap0[nBinPtRap0+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,350,400,450,500,1000};
  const int nBinPtRap1 = 37; Float_t xbinsPtRap1[nBinPtRap1+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,350,400,450,500,1000};
  const int nBinPtRap2 = 37; Float_t xbinsPtRap2[nBinPtRap2+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,350,400,450,500,1000};
  const int nBinPtRap3 = 37; Float_t xbinsPtRap3[nBinPtRap3+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,350,400,450,500,1000};
  const int nBinPtRap4 = 37; Float_t xbinsPtRap4[nBinPtRap4+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,350,400,450,500,1000};

  const int nBinRecoPt     = 74; Float_t xbinsRecoPt[nBinRecoPt+1]         = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300,325,350,375,400,425,450,475,500, 750,1000};

  const int nBinRecoPtRap0 = 74; Float_t xbinsRecoPtRap0[nBinRecoPtRap0+1] = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300,325,350,375,400,425,450,475,500, 750,1000};

  const int nBinRecoPtRap1 = 74; Float_t xbinsRecoPtRap1[nBinRecoPtRap1+1] = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300,325,350,375,400,425,450,475,500, 750,1000};

  const int nBinRecoPtRap2 = 74; Float_t xbinsRecoPtRap2[nBinRecoPtRap2+1] = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300,325,350,375,400,425,450,475,500, 750,1000};

  const int nBinRecoPtRap3 = 74; Float_t xbinsRecoPtRap3[nBinRecoPtRap3+1] = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300,325,350,375,400,425,450,475,500, 750,1000};

  const int nBinRecoPtRap4 = 74; Float_t xbinsRecoPtRap4[nBinRecoPtRap4+1] = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
                                                                              10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
                                                                              28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
                                                                             220, 235,250, 275,300,325,350,375,400,425,450,475,500, 750,1000};

  const int nBinRap = 24; Float_t xbinsRap[nBinRap+1]; for(int i=0; i<=nBinRap;i++) xbinsRap[i] = i * 0.1;

  TFile *fLepton_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_dylan_37ifb.root"));
  TH2D* scalefactors_Medium_Muon                       	   = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon");                          scalefactors_Medium_Muon			     ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_syst_error_combined   	   = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_syst_error_combined");      scalefactors_Medium_Muon_syst_error_combined       ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_stat_error_hi         	   = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_stat_error_hi");	          scalefactors_Medium_Muon_stat_error_hi	     ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_syst_error_alt_bkg    	   = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_syst_error_alt_bkg");	  scalefactors_Medium_Muon_syst_error_alt_bkg        ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_syst_error_generator  	   = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_syst_error_generator");     scalefactors_Medium_Muon_syst_error_generator      ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_syst_error_alt_signal 	   = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_syst_error_alt_signal");    scalefactors_Medium_Muon_syst_error_alt_signal     ->SetDirectory(0);
  TH2D* scalefactors_Medium_Muon_syst_error_alt_tag    	   = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Muon_syst_error_alt_tag");	  scalefactors_Medium_Muon_syst_error_alt_tag        ->SetDirectory(0);
  const int nMuSFBins = 72;
  TH2D* scalefactors_Medium_Muon_stat_error_hi_bins[nMuSFBins];
  for(int nj=0; nj<nMuSFBins; nj++) {scalefactors_Medium_Muon_stat_error_hi_bins[nj] =(TH2D*)fLepton_SF->Get(Form("scalefactors_Medium_Muon_stat_error_hi_bins_%d",nj));scalefactors_Medium_Muon_stat_error_hi_bins[nj]->SetDirectory(0);}
  TH2D* scalefactors_Medium_Electron                       = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron");		          scalefactors_Medium_Electron		             ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_syst_error_combined   = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_syst_error_combined");  scalefactors_Medium_Electron_syst_error_combined   ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_stat_error_hi         = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_stat_error_hi");	  scalefactors_Medium_Electron_stat_error_hi         ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_syst_error_alt_bkg    = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_syst_error_alt_bkg");   scalefactors_Medium_Electron_syst_error_alt_bkg    ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_syst_error_generator  = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_syst_error_generator"); scalefactors_Medium_Electron_syst_error_generator  ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_syst_error_alt_signal = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_syst_error_alt_signal");scalefactors_Medium_Electron_syst_error_alt_signal ->SetDirectory(0);
  TH2D* scalefactors_Medium_Electron_syst_error_alt_tag    = (TH2D*)fLepton_SF->Get("scalefactors_Medium_Electron_syst_error_alt_tag");   scalefactors_Medium_Electron_syst_error_alt_tag    ->SetDirectory(0);
  const int nElSFBins = 140;
  TH2D* scalefactors_Medium_Electron_stat_error_hi_bins[nElSFBins];
  for(int nj=0; nj<nElSFBins; nj++) {scalefactors_Medium_Electron_stat_error_hi_bins[nj] = (TH2D*)fLepton_SF->Get(Form("scalefactors_Medium_Electron_stat_error_hi_bins_%d",nj)); scalefactors_Medium_Electron_stat_error_hi_bins[nj]->SetDirectory(0);}

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

  const int nEffNuisances = 145;
  const int nMomNuisances = 5;

  TH2D* histoPtRecGen[2];
  TH2D* histoPtRecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRecDA[2];
  TH1D* histoPtRecDY[2];
  TH1D* histoPtRecEM[2];
  TH1D* histoPtRecVV[2];
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
  TH2D* histoRapRecGen_LepEff[2][nEffNuisances];
  TH2D* histoRapRecGen_MomRes[2][nMomNuisances];
  TH1D* histoRapRecDA[2];
  TH1D* histoRapRecDY[2];
  TH1D* histoRapRecEM[2];
  TH1D* histoRapRecVV[2];
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
   histoRapRecGen[i]        = new TH2D(Form("histoRapRecGen_%d",i),        Form("histoRapRecGen_%d",i),        nBinRap, xbinsRap, nBinRap, xbinsRap);
   for(int j=0; j<nEffNuisances; j++){
     histoRapRecGen_LepEff[i][j] = new TH2D(Form("histoRapRecGen_LepEff_%d_%d",i,j), Form("histoRapRecGen_LepEff_%d_%d",i,j), nBinRap, xbinsRap, nBinRap, xbinsRap);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoRapRecGen_MomRes[i][j] = new TH2D(Form("histoRapRecGen_MomRes_%d_%d",i,j), Form("histoRapRecGen_MomRes_%d_%d",i,j), nBinRap, xbinsRap, nBinRap, xbinsRap);
   }
   histoRapRecDA[i]  = new TH1D(Form("histoRapRecDA_%d",i),  Form("histoRapRecDA_%d",i),  nBinRap, xbinsRap);
   histoRapRecDY[i]  = new TH1D(Form("histoRapRecDY_%d",i),  Form("histoRapRecDY_%d",i),  nBinRap, xbinsRap);
   histoRapRecEM[i]  = new TH1D(Form("histoRapRecEM_%d",i),  Form("histoRapRecEM_%d",i),  nBinRap, xbinsRap);
   histoRapRecVV[i]  = new TH1D(Form("histoRapRecVV_%d",i),  Form("histoRapRecVV_%d",i),  nBinRap, xbinsRap);

   for(int j=0; j<nEffNuisances; j++){
     histoRapRecDY_LepEff[i][j] = new TH1D(Form("histoRapRecDY_LepEff_%d_%d",i,j), Form("histoRapRecDY_LepEff_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapRecVV_LepEff[i][j] = new TH1D(Form("histoRapRecVV_LepEff_%d_%d",i,j), Form("histoRapRecVV_LepEff_%d_%d",i,j), nBinRap, xbinsRap);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoRapRecDA_MomRes[i][j] = new TH1D(Form("histoRapRecDA_MomRes_%d_%d",i,j), Form("histoRapRecDA_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapRecDY_MomRes[i][j] = new TH1D(Form("histoRapRecDY_MomRes_%d_%d",i,j), Form("histoRapRecDY_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapRecEM_MomRes[i][j] = new TH1D(Form("histoRapRecEM_MomRes_%d_%d",i,j), Form("histoRapRecEM_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapRecVV_MomRes[i][j] = new TH1D(Form("histoRapRecVV_MomRes_%d_%d",i,j), Form("histoRapRecVV_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
   }
   histoRapRecVV_PDF[i] = new TH1D(Form("histoRapRecVV_PDF_%d",i), Form("histoRapRecVV_PDF_%d",i), nBinRap, xbinsRap);
   histoRapRecDY_PDF[i] = new TH1D(Form("histoRapRecDY_PDF_%d",i), Form("histoRapRecDY_PDF_%d",i), nBinRap, xbinsRap);

   histoRapRecVV_QCD[i] = new TH1D(Form("histoRapRecVV_QCD_%d",i), Form("histoRapRecVV_QCD_%d",i), nBinRap, xbinsRap);
   histoRapRecDY_QCD[i] = new TH1D(Form("histoRapRecDY_QCD_%d",i), Form("histoRapRecDY_QCD_%d",i), nBinRap, xbinsRap);

   for(int j=0; j<6; j++){
      histoRapRecVV_QCDPart[i][j] = new TH1D(Form("histoRapRecVV_QCDPart_%d_%d",i,j), Form("histoRapRecVV_QCDPart_%d_%d",i,j), nBinRap, xbinsRap);
      histoRapRecDY_QCDPart[i][j] = new TH1D(Form("histoRapRecDY_QCDPart_%d_%d",i,j), Form("histoRapRecDY_QCDPart_%d_%d",i,j), nBinRap, xbinsRap);
    }
  }

  TH2D* histoRapPRecGen[2];
  TH2D* histoRapPRecGen_LepEff[2][nEffNuisances];
  TH2D* histoRapPRecGen_MomRes[2][nMomNuisances];
  TH1D* histoRapPRecDA[2];
  TH1D* histoRapPRecDY[2];
  TH1D* histoRapPRecEM[2];
  TH1D* histoRapPRecVV[2];
  TH1D* histoRapPRecDY_LepEff[2][nEffNuisances];
  TH1D* histoRapPRecVV_LepEff[2][nEffNuisances];
  TH1D* histoRapPRecDA_MomRes[2][nMomNuisances];
  TH1D* histoRapPRecDY_MomRes[2][nMomNuisances];
  TH1D* histoRapPRecEM_MomRes[2][nMomNuisances];
  TH1D* histoRapPRecVV_MomRes[2][nMomNuisances];
  TH1D* histoRapPRecVV_PDF[2];
  TH1D* histoRapPRecDY_PDF[2];
  TH1D* histoRapPRecVV_QCD[2];
  TH1D* histoRapPRecDY_QCD[2];
  TH1D* histoRapPRecVV_QCDPart[2][6];
  TH1D* histoRapPRecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoRapPRecGen[i]        = new TH2D(Form("histoRapPRecGen_%d",i),        Form("histoRapPRecGen_%d",i),        nBinRap, xbinsRap, nBinRap, xbinsRap);
   for(int j=0; j<nEffNuisances; j++){
     histoRapPRecGen_LepEff[i][j] = new TH2D(Form("histoRapPRecGen_LepEff_%d_%d",i,j), Form("histoRapPRecGen_LepEff_%d_%d",i,j), nBinRap, xbinsRap, nBinRap, xbinsRap);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoRapPRecGen_MomRes[i][j] = new TH2D(Form("histoRapPRecGen_MomRes_%d_%d",i,j), Form("histoRapPRecGen_MomRes_%d_%d",i,j), nBinRap, xbinsRap, nBinRap, xbinsRap);
   }
   histoRapPRecDA[i]  = new TH1D(Form("histoRapPRecDA_%d",i),  Form("histoRapPRecDA_%d",i),  nBinRap, xbinsRap);
   histoRapPRecDY[i]  = new TH1D(Form("histoRapPRecDY_%d",i),  Form("histoRapPRecDY_%d",i),  nBinRap, xbinsRap);
   histoRapPRecEM[i]  = new TH1D(Form("histoRapPRecEM_%d",i),  Form("histoRapPRecEM_%d",i),  nBinRap, xbinsRap);
   histoRapPRecVV[i]  = new TH1D(Form("histoRapPRecVV_%d",i),  Form("histoRapPRecVV_%d",i),  nBinRap, xbinsRap);

   for(int j=0; j<nEffNuisances; j++){
     histoRapPRecDY_LepEff[i][j] = new TH1D(Form("histoRapPRecDY_LepEff_%d_%d",i,j), Form("histoRapPRecDY_LepEff_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapPRecVV_LepEff[i][j] = new TH1D(Form("histoRapPRecVV_LepEff_%d_%d",i,j), Form("histoRapPRecVV_LepEff_%d_%d",i,j), nBinRap, xbinsRap);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoRapPRecDA_MomRes[i][j] = new TH1D(Form("histoRapPRecDA_MomRes_%d_%d",i,j), Form("histoRapPRecDA_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapPRecDY_MomRes[i][j] = new TH1D(Form("histoRapPRecDY_MomRes_%d_%d",i,j), Form("histoRapPRecDY_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapPRecEM_MomRes[i][j] = new TH1D(Form("histoRapPRecEM_MomRes_%d_%d",i,j), Form("histoRapPRecEM_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapPRecVV_MomRes[i][j] = new TH1D(Form("histoRapPRecVV_MomRes_%d_%d",i,j), Form("histoRapPRecVV_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
   }
   histoRapPRecVV_PDF[i] = new TH1D(Form("histoRapPRecVV_PDF_%d",i), Form("histoRapPRecVV_PDF_%d",i), nBinRap, xbinsRap);
   histoRapPRecDY_PDF[i] = new TH1D(Form("histoRapPRecDY_PDF_%d",i), Form("histoRapPRecDY_PDF_%d",i), nBinRap, xbinsRap);

   histoRapPRecVV_QCD[i] = new TH1D(Form("histoRapPRecVV_QCD_%d",i), Form("histoRapPRecVV_QCD_%d",i), nBinRap, xbinsRap);
   histoRapPRecDY_QCD[i] = new TH1D(Form("histoRapPRecDY_QCD_%d",i), Form("histoRapPRecDY_QCD_%d",i), nBinRap, xbinsRap);

   for(int j=0; j<6; j++){
      histoRapPRecVV_QCDPart[i][j] = new TH1D(Form("histoRapPRecVV_QCDPart_%d_%d",i,j), Form("histoRapPRecVV_QCDPart_%d_%d",i,j), nBinRap, xbinsRap);
      histoRapPRecDY_QCDPart[i][j] = new TH1D(Form("histoRapPRecDY_QCDPart_%d_%d",i,j), Form("histoRapPRecDY_QCDPart_%d_%d",i,j), nBinRap, xbinsRap);
    }
  }

  TH2D* histoRapMRecGen[2];
  TH2D* histoRapMRecGen_LepEff[2][nEffNuisances];
  TH2D* histoRapMRecGen_MomRes[2][nMomNuisances];
  TH1D* histoRapMRecDA[2];
  TH1D* histoRapMRecDY[2];
  TH1D* histoRapMRecEM[2];
  TH1D* histoRapMRecVV[2];
  TH1D* histoRapMRecDY_LepEff[2][nEffNuisances];
  TH1D* histoRapMRecVV_LepEff[2][nEffNuisances];
  TH1D* histoRapMRecDA_MomRes[2][nMomNuisances];
  TH1D* histoRapMRecDY_MomRes[2][nMomNuisances];
  TH1D* histoRapMRecEM_MomRes[2][nMomNuisances];
  TH1D* histoRapMRecVV_MomRes[2][nMomNuisances];
  TH1D* histoRapMRecVV_PDF[2];
  TH1D* histoRapMRecDY_PDF[2];
  TH1D* histoRapMRecVV_QCD[2];
  TH1D* histoRapMRecDY_QCD[2];
  TH1D* histoRapMRecVV_QCDPart[2][6];
  TH1D* histoRapMRecDY_QCDPart[2][6];
  for(int i=0; i<2; i++){
   histoRapMRecGen[i]        = new TH2D(Form("histoRapMRecGen_%d",i),        Form("histoRapMRecGen_%d",i),        nBinRap, xbinsRap, nBinRap, xbinsRap);
   for(int j=0; j<nEffNuisances; j++){
     histoRapMRecGen_LepEff[i][j] = new TH2D(Form("histoRapMRecGen_LepEff_%d_%d",i,j), Form("histoRapMRecGen_LepEff_%d_%d",i,j), nBinRap, xbinsRap, nBinRap, xbinsRap);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoRapMRecGen_MomRes[i][j] = new TH2D(Form("histoRapMRecGen_MomRes_%d_%d",i,j), Form("histoRapMRecGen_MomRes_%d_%d",i,j), nBinRap, xbinsRap, nBinRap, xbinsRap);
   }
   histoRapMRecDA[i]  = new TH1D(Form("histoRapMRecDA_%d",i),  Form("histoRapMRecDA_%d",i),  nBinRap, xbinsRap);
   histoRapMRecDY[i]  = new TH1D(Form("histoRapMRecDY_%d",i),  Form("histoRapMRecDY_%d",i),  nBinRap, xbinsRap);
   histoRapMRecEM[i]  = new TH1D(Form("histoRapMRecEM_%d",i),  Form("histoRapMRecEM_%d",i),  nBinRap, xbinsRap);
   histoRapMRecVV[i]  = new TH1D(Form("histoRapMRecVV_%d",i),  Form("histoRapMRecVV_%d",i),  nBinRap, xbinsRap);

   for(int j=0; j<nEffNuisances; j++){
     histoRapMRecDY_LepEff[i][j] = new TH1D(Form("histoRapMRecDY_LepEff_%d_%d",i,j), Form("histoRapMRecDY_LepEff_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapMRecVV_LepEff[i][j] = new TH1D(Form("histoRapMRecVV_LepEff_%d_%d",i,j), Form("histoRapMRecVV_LepEff_%d_%d",i,j), nBinRap, xbinsRap);
   }
   for(int j=0; j<nMomNuisances; j++){
     histoRapMRecDA_MomRes[i][j] = new TH1D(Form("histoRapMRecDA_MomRes_%d_%d",i,j), Form("histoRapMRecDA_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapMRecDY_MomRes[i][j] = new TH1D(Form("histoRapMRecDY_MomRes_%d_%d",i,j), Form("histoRapMRecDY_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapMRecEM_MomRes[i][j] = new TH1D(Form("histoRapMRecEM_MomRes_%d_%d",i,j), Form("histoRapMRecEM_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
     histoRapMRecVV_MomRes[i][j] = new TH1D(Form("histoRapMRecVV_MomRes_%d_%d",i,j), Form("histoRapMRecVV_MomRes_%d_%d",i,j), nBinRap, xbinsRap);
   }
   histoRapMRecVV_PDF[i] = new TH1D(Form("histoRapMRecVV_PDF_%d",i), Form("histoRapMRecVV_PDF_%d",i), nBinRap, xbinsRap);
   histoRapMRecDY_PDF[i] = new TH1D(Form("histoRapMRecDY_PDF_%d",i), Form("histoRapMRecDY_PDF_%d",i), nBinRap, xbinsRap);

   histoRapMRecVV_QCD[i] = new TH1D(Form("histoRapMRecVV_QCD_%d",i), Form("histoRapMRecVV_QCD_%d",i), nBinRap, xbinsRap);
   histoRapMRecDY_QCD[i] = new TH1D(Form("histoRapMRecDY_QCD_%d",i), Form("histoRapMRecDY_QCD_%d",i), nBinRap, xbinsRap);

   for(int j=0; j<6; j++){
      histoRapMRecVV_QCDPart[i][j] = new TH1D(Form("histoRapMRecVV_QCDPart_%d_%d",i,j), Form("histoRapMRecVV_QCDPart_%d_%d",i,j), nBinRap, xbinsRap);
      histoRapMRecDY_QCDPart[i][j] = new TH1D(Form("histoRapMRecDY_QCDPart_%d_%d",i,j), Form("histoRapMRecDY_QCDPart_%d_%d",i,j), nBinRap, xbinsRap);
    }
  }

  TH2D* histoPtRap0RecGen[2];
  TH2D* histoPtRap0RecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRap0RecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRap0RecDA[2];
  TH1D* histoPtRap0RecDY[2];
  TH1D* histoPtRap0RecEM[2];
  TH1D* histoPtRap0RecVV[2];
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
  TH2D* histoPtRap1RecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRap1RecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRap1RecDA[2];
  TH1D* histoPtRap1RecDY[2];
  TH1D* histoPtRap1RecEM[2];
  TH1D* histoPtRap1RecVV[2];
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
  TH2D* histoPtRap2RecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRap2RecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRap2RecDA[2];
  TH1D* histoPtRap2RecDY[2];
  TH1D* histoPtRap2RecEM[2];
  TH1D* histoPtRap2RecVV[2];
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
  TH2D* histoPtRap3RecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRap3RecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRap3RecDA[2];
  TH1D* histoPtRap3RecDY[2];
  TH1D* histoPtRap3RecEM[2];
  TH1D* histoPtRap3RecVV[2];
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
  TH2D* histoPtRap4RecGen_LepEff[2][nEffNuisances];
  TH2D* histoPtRap4RecGen_MomRes[2][nMomNuisances];
  TH1D* histoPtRap4RecDA[2];
  TH1D* histoPtRap4RecDY[2];
  TH1D* histoPtRap4RecEM[2];
  TH1D* histoPtRap4RecVV[2];
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
      if(thePandaFlat.metFilter == 0) continue;
      
      if(thePandaFlat.nLooseLep != 2) continue;

      int theCategory = infileCat_[ifile];
      float lepPtSF[2] = {1,1};
      float lepPtSFSyst[2][nMomNuisances] = {1,1,1,1,1,1,1,1,1,1};
      if(theCategory == 0) { // Data
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
      bool passSel = TMath::Abs((v1+v2).M()-91.1876) < 15 && v1.Pt() > 25 && v2.Pt() > 25;
      bool passSystSel[nMomNuisances] = {TMath::Abs((vMomRes1[0]+vMomRes2[0]).M()-91.1876) < 15 && vMomRes1[0].Pt() > 25 && vMomRes2[0].Pt() > 25,
                                         TMath::Abs((vMomRes1[1]+vMomRes2[1]).M()-91.1876) < 15 && vMomRes1[1].Pt() > 25 && vMomRes2[1].Pt() > 25,
                                         TMath::Abs((vMomRes1[2]+vMomRes2[2]).M()-91.1876) < 15 && vMomRes1[2].Pt() > 25 && vMomRes2[2].Pt() > 25,
                                         TMath::Abs((vMomRes1[3]+vMomRes2[3]).M()-91.1876) < 15 && vMomRes1[3].Pt() > 25 && vMomRes2[3].Pt() > 25,
                                         TMath::Abs((vMomRes1[4]+vMomRes2[4]).M()-91.1876) < 15 && vMomRes1[4].Pt() > 25 && vMomRes2[4].Pt() > 25};

      double ZRecPt  = (v1+v2).Pt();
      double ZRecRap = TMath::Abs((v1+v2).Rapidity());

      double ZRecSystPt[nMomNuisances] = {(vMomRes1[0]+vMomRes2[0]).Pt(),
                                          (vMomRes1[1]+vMomRes2[1]).Pt(),
                                          (vMomRes1[2]+vMomRes2[2]).Pt(),
                                          (vMomRes1[3]+vMomRes2[3]).Pt(),
                                          (vMomRes1[4]+vMomRes2[4]).Pt()};
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
      double totalWeightLepEff[2][nEffNuisances];
      for(int ni=0; ni<2; ni++) for(int nj=0; nj<nEffNuisances; nj++) totalWeightLepEff[ni][nj] = 0.0;

      if(theCategory != 0 && theCategory != 5){
        totalWeight = thePandaFlat.normalizedWeight * lumi * thePandaFlat.sf_pu *
	              thePandaFlat.sf_trk1 * thePandaFlat.sf_medium1 *
		      thePandaFlat.sf_trk2 * thePandaFlat.sf_medium2;
        if(abs(thePandaFlat.looseLep1PdgId)==13){
          double etal = thePandaFlat.looseLep1Eta; if(etal >= 2.4) etal = 2.3999; else if(etal <= -2.4) etal = -2.3999;
          double ptl = thePandaFlat.looseLep1Pt; if(ptl>=1000) ptl = 999.999;
          int binXT = scalefactors_Medium_Muon->GetXaxis()->FindFixBin(etal);
          int binYT = scalefactors_Medium_Muon->GetYaxis()->FindFixBin(ptl);
          totalWeightLepEff[0][0] = scalefactors_Medium_Muon_syst_error_combined  ->GetBinContent(binXT,binYT);
          totalWeightLepEff[0][1] = scalefactors_Medium_Muon_syst_error_alt_bkg   ->GetBinContent(binXT,binYT);
          totalWeightLepEff[0][2] = scalefactors_Medium_Muon_syst_error_generator ->GetBinContent(binXT,binYT);
          totalWeightLepEff[0][3] = scalefactors_Medium_Muon_syst_error_alt_signal->GetBinContent(binXT,binYT);
          totalWeightLepEff[0][4] = scalefactors_Medium_Muon_syst_error_alt_tag   ->GetBinContent(binXT,binYT);
          for(int nj=0; nj<nMuSFBins; nj++) totalWeightLepEff[0][nj+5] = scalefactors_Medium_Muon_stat_error_hi_bins[nj]->GetBinContent(binXT,binYT);
        } else {
          double etal = thePandaFlat.looseLep1Eta; if(etal >= 2.4) etal = 2.3999; else if(etal <= -2.4) etal = -2.3999;
          double ptl = thePandaFlat.looseLep1Pt; if(ptl>=1000) ptl = 999.999;
          int binXT = scalefactors_Medium_Electron->GetXaxis()->FindFixBin(etal);
          int binYT = scalefactors_Medium_Electron->GetYaxis()->FindFixBin(ptl);
          totalWeightLepEff[0][0] = scalefactors_Medium_Electron_syst_error_combined  ->GetBinContent(binXT,binYT);
          totalWeightLepEff[0][1] = scalefactors_Medium_Electron_syst_error_alt_bkg   ->GetBinContent(binXT,binYT);
          totalWeightLepEff[0][2] = scalefactors_Medium_Electron_syst_error_generator ->GetBinContent(binXT,binYT);
          totalWeightLepEff[0][3] = scalefactors_Medium_Electron_syst_error_alt_signal->GetBinContent(binXT,binYT);
          totalWeightLepEff[0][4] = scalefactors_Medium_Electron_syst_error_alt_tag   ->GetBinContent(binXT,binYT);
          for(int nj=0; nj<nElSFBins; nj++) totalWeightLepEff[0][nj+5] = scalefactors_Medium_Electron_stat_error_hi_bins[nj]->GetBinContent(binXT,binYT);
        }
        if(abs(thePandaFlat.looseLep2PdgId)==13){
          double etal = thePandaFlat.looseLep2Eta; if(etal >= 2.4) etal = 2.3999; else if(etal <= -2.4) etal = -2.3999;
          double ptl = thePandaFlat.looseLep2Pt; if(ptl>=1000) ptl = 999.999;
          int binXT = scalefactors_Medium_Muon->GetXaxis()->FindFixBin(etal);
          int binYT = scalefactors_Medium_Muon->GetYaxis()->FindFixBin(ptl);
          totalWeightLepEff[1][0] = scalefactors_Medium_Muon_syst_error_combined  ->GetBinContent(binXT,binYT);
          totalWeightLepEff[1][1] = scalefactors_Medium_Muon_syst_error_alt_bkg   ->GetBinContent(binXT,binYT);
          totalWeightLepEff[1][2] = scalefactors_Medium_Muon_syst_error_generator ->GetBinContent(binXT,binYT);
          totalWeightLepEff[1][3] = scalefactors_Medium_Muon_syst_error_alt_signal->GetBinContent(binXT,binYT);
          totalWeightLepEff[1][4] = scalefactors_Medium_Muon_syst_error_alt_tag   ->GetBinContent(binXT,binYT);
          for(int nj=0; nj<nMuSFBins; nj++) totalWeightLepEff[1][nj+5] = scalefactors_Medium_Muon_stat_error_hi_bins[nj]->GetBinContent(binXT,binYT);
        } else {
          double etal = thePandaFlat.looseLep2Eta; if(etal >= 2.4) etal = 2.3999; else if(etal <= -2.4) etal = -2.3999;
          double ptl = thePandaFlat.looseLep2Pt; if(ptl>=1000) ptl = 999.999;
          int binXT = scalefactors_Medium_Electron->GetXaxis()->FindFixBin(etal);
          int binYT = scalefactors_Medium_Electron->GetYaxis()->FindFixBin(ptl);
          totalWeightLepEff[1][0] = scalefactors_Medium_Electron_syst_error_combined  ->GetBinContent(binXT,binYT);
          totalWeightLepEff[1][1] = scalefactors_Medium_Electron_syst_error_alt_bkg   ->GetBinContent(binXT,binYT);
          totalWeightLepEff[1][2] = scalefactors_Medium_Electron_syst_error_generator ->GetBinContent(binXT,binYT);
          totalWeightLepEff[1][3] = scalefactors_Medium_Electron_syst_error_alt_signal->GetBinContent(binXT,binYT);
          totalWeightLepEff[1][4] = scalefactors_Medium_Electron_syst_error_alt_tag   ->GetBinContent(binXT,binYT);
          for(int nj=0; nj<nElSFBins; nj++) totalWeightLepEff[1][nj+5] = scalefactors_Medium_Electron_stat_error_hi_bins[nj]->GetBinContent(binXT,binYT);
        }
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
            double theKeff = k_eff;
            if(theLepType == 1) theKeff = 1.0/k_eff;
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
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRecDY_LepEff[lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
          histoPtRecDY_PDF    [lepType]   ->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
          histoPtRecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
          histoPtRecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          // Rap
          if(ZRecRap < 2.4){
            histoRapRecDY        [lepType]   ->Fill(ZRecRap,totalWeight);
            for(int nj=0; nj<nEffNuisances; nj++) histoRapRecDY_LepEff[lepType][nj]->Fill(ZRecRap,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
            histoRapRecDY_PDF    [lepType]   ->Fill(ZRecRap,totalWeight*thePandaFlat.pdfUp);
            histoRapRecDY_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoRapRecDY_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoRapRecDY_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoRapRecDY_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoRapRecDY_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoRapRecDY_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
            if(thePandaFlat.looseLep1PdgId < 0){
              histoRapPRecDY        [lepType]   ->Fill(ZRecRap,totalWeight);
              for(int nj=0; nj<nEffNuisances; nj++) histoRapPRecDY_LepEff[lepType][nj]->Fill(ZRecRap,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
              histoRapPRecDY_PDF    [lepType]   ->Fill(ZRecRap,totalWeight*thePandaFlat.pdfUp);
              histoRapPRecDY_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
              histoRapPRecDY_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
              histoRapPRecDY_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
              histoRapPRecDY_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
              histoRapPRecDY_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
              histoRapPRecDY_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
            } else {
              histoRapMRecDY        [lepType]   ->Fill(ZRecRap,totalWeight);
              for(int nj=0; nj<nEffNuisances; nj++) histoRapMRecDY_LepEff[lepType][nj]->Fill(ZRecRap,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
              histoRapMRecDY_PDF    [lepType]   ->Fill(ZRecRap,totalWeight*thePandaFlat.pdfUp);
              histoRapMRecDY_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
              histoRapMRecDY_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
              histoRapMRecDY_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
              histoRapMRecDY_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
              histoRapMRecDY_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
              histoRapMRecDY_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
            }
          }
          // PtRap0
          if      (ZRecRap < 0.5){
            histoPtRap0RecDY        [lepType]   ->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap0RecDY_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap1RecDY_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap2RecDY_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
            histoPtRap2RecDY_PDF    [lepType]	->Fill(ZRecPt,totalWeight*thePandaFlat.pdfUp);
            histoPtRap2RecDY_QCDPart[lepType][0]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
            histoPtRap2RecDY_QCDPart[lepType][1]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
            histoPtRap2RecDY_QCDPart[lepType][2]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
            histoPtRap2RecDY_QCDPart[lepType][3]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
            histoPtRap2RecDY_QCDPart[lepType][4]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
            histoPtRap2RecDY_QCDPart[lepType][5]->Fill(ZRecPt,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
          }
          else if(ZRecRap < 2.0){
            histoPtRap3RecDY	    [lepType]	->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap3RecDY_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap4RecDY_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
          // Pt
          histoPtRecDA[lepType]->Fill(ZRecPt,totalWeight);
          // Rap
          if(ZRecRap < 2.4){
          histoRapRecDA[lepType]->Fill(ZRecRap,totalWeight);
            if(thePandaFlat.looseLep1PdgId < 0){
              histoRapPRecDA	  [lepType]->Fill(ZRecRap,totalWeight);
            } else {
              histoRapMRecDA	  [lepType]->Fill(ZRecRap,totalWeight);
            }
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
          else if(ZRecRap < 2.0){
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
            // Pt
            histoPtRecEM[theLepType]->Fill(ZRecPt,totalWeight*theKeff);
            // Rap
            if(ZRecRap < 2.4){
              histoRapRecEM[theLepType]->Fill(ZRecRap,totalWeight*theKeff);
            if(thePandaFlat.looseLep1PdgId < 0){
              histoRapPRecEM[theLepType]->Fill(ZRecRap,totalWeight*theKeff);
            } else {
              histoRapMRecEM[theLepType]->Fill(ZRecRap,totalWeight*theKeff);
            }
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
            else if(ZRecRap < 2.0){
              histoPtRap3RecEM[theLepType]->Fill(ZRecPt,totalWeight*theKeff);
            }
            else if(ZRecRap < 2.4){
              histoPtRap4RecEM[theLepType]->Fill(ZRecPt,totalWeight*theKeff);
            }
          } // loop over muons and electrons
	}
	else { // VV
          // Pt
          histoPtRecVV        [lepType]->Fill(ZRecPt,totalWeight);
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRecVV_LepEff [lepType][nj]   ->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
          // Rap
          if(ZRecRap < 2.4){
            histoRapRecVV        [lepType]     ->Fill(ZRecRap,totalWeight);
            for(int nj=0; nj<nEffNuisances; nj++) histoRapRecVV_LepEff [lepType][nj]->Fill(ZRecRap,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
            if(thePandaFlat.looseLep1PdgId < 0){
              histoRapPRecVV        [lepType]     ->Fill(ZRecRap,totalWeight);
              for(int nj=0; nj<nEffNuisances; nj++) histoRapPRecVV_LepEff [lepType][nj]->Fill(ZRecRap,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
              histoRapPRecVV_PDF    [lepType]     ->Fill(ZRecRap,totalWeight*thePandaFlat.pdfUp);
              if(thePandaFlat.scale[0] != -1){
        	histoRapPRecVV_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
        	histoRapPRecVV_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
        	histoRapPRecVV_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
        	histoRapPRecVV_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
        	histoRapPRecVV_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
        	histoRapPRecVV_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
              } else {
        	histoRapPRecVV_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight);
        	histoRapPRecVV_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight);
        	histoRapPRecVV_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight);
        	histoRapPRecVV_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight);
        	histoRapPRecVV_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight);
        	histoRapPRecVV_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight);
              }
            } else {
              histoRapMRecVV        [lepType]     ->Fill(ZRecRap,totalWeight);
              for(int nj=0; nj<nEffNuisances; nj++) histoRapMRecVV_LepEff [lepType][nj]->Fill(ZRecRap,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
              histoRapMRecVV_PDF    [lepType]     ->Fill(ZRecRap,totalWeight*thePandaFlat.pdfUp);
              if(thePandaFlat.scale[0] != -1){
        	histoRapMRecVV_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[0])/maxQCDscale);
        	histoRapMRecVV_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[1])/maxQCDscale);
        	histoRapMRecVV_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[2])/maxQCDscale);
        	histoRapMRecVV_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[3])/maxQCDscale);
        	histoRapMRecVV_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[4])/maxQCDscale);
        	histoRapMRecVV_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight*TMath::Abs(1+thePandaFlat.scale[5])/maxQCDscale);
              } else {
        	histoRapMRecVV_QCDPart[lepType][0]->Fill(ZRecRap,totalWeight);
        	histoRapMRecVV_QCDPart[lepType][1]->Fill(ZRecRap,totalWeight);
        	histoRapMRecVV_QCDPart[lepType][2]->Fill(ZRecRap,totalWeight);
        	histoRapMRecVV_QCDPart[lepType][3]->Fill(ZRecRap,totalWeight);
        	histoRapMRecVV_QCDPart[lepType][4]->Fill(ZRecRap,totalWeight);
        	histoRapMRecVV_QCDPart[lepType][5]->Fill(ZRecRap,totalWeight);
              }
            }
          }
          // PtRap0
          if     (ZRecRap < 0.5){
            histoPtRap0RecVV        [lepType]     ->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap0RecVV_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap1RecVV_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap2RecVV_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
          else if(ZRecRap < 2.0){
            histoPtRap3RecVV        [lepType]     ->Fill(ZRecPt,totalWeight);
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap3RecVV_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
            for(int nj=0; nj<nEffNuisances; nj++) histoPtRap4RecVV_LepEff [lepType][nj]->Fill(ZRecPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
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
          histoPtRecGen       [lepType] ->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
	}
	if(theCategory == 2 && passRapFid == true && ZRecRap < 2.4){
          histoRapRecGen       [lepType]->Fill(ZRecRap,ZGenRap,totalWeight);
          for(int nj=0; nj<nEffNuisances; nj++) histoRapRecGen_LepEff [lepType][nj]->Fill(ZRecRap,ZGenRap,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
          if(thePandaFlat.looseLep1PdgId < 0){
            histoRapPRecGen       [lepType]->Fill(ZRecRap,ZGenRap,totalWeight);
            for(int nj=0; nj<nEffNuisances; nj++) histoRapPRecGen_LepEff [lepType][nj]->Fill(ZRecRap,ZGenRap,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
          } else {
            histoRapMRecGen       [lepType]->Fill(ZRecRap,ZGenRap,totalWeight);
            for(int nj=0; nj<nEffNuisances; nj++) histoRapMRecGen_LepEff [lepType][nj]->Fill(ZRecRap,ZGenRap,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
          }
	}
	if(theCategory == 2 && passPtRapFid[0] == true && ZRecRap >= 0.0 && ZRecRap < 0.5){
          histoPtRap0RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRap0RecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
	}
	if(theCategory == 2 && passPtRapFid[1] == true && ZRecRap >= 0.5 && ZRecRap < 1.0){
          histoPtRap1RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRap1RecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
	}
	if(theCategory == 2 && passPtRapFid[2] == true && ZRecRap >= 1.0 && ZRecRap < 1.5){
          histoPtRap2RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRap2RecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
	}
	if(theCategory == 2 && passPtRapFid[3] == true && ZRecRap >= 1.5 && ZRecRap < 2.0){
          histoPtRap3RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRap3RecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
	}
	if(theCategory == 2 && passPtRapFid[4] == true && ZRecRap >= 2.0 && ZRecRap < 2.4){
          histoPtRap4RecGen       [lepType]->Fill(ZRecPt,ZGenPt,totalWeight);
          for(int nj=0; nj<nEffNuisances; nj++) histoPtRap4RecGen_LepEff [lepType][nj]->Fill(ZRecPt,ZGenPt,totalWeight*(1+totalWeightLepEff[0][nj])*(1+totalWeightLepEff[1][nj]));
	}
      } // end passSel
      
      for(int nc=0; nc<nMomNuisances; nc++){
	if(passSystSel[nc]){
	  if     (theCategory == 2){
            histoPtRecDY_MomRes[lepType][nc] ->Fill(ZRecSystPt[nc],totalWeight);
            histoPtRecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);
            if(ZRecRap < 2.4) {
              histoRapRecDY_MomRes[lepType][nc] ->Fill(ZRecRap,totalWeight);
              histoRapRecGen_MomRes[lepType][nc]->Fill(ZRecRap,ZGenRap,totalWeight);
              if(thePandaFlat.looseLep1PdgId < 0){
        	histoRapPRecDY_MomRes[lepType][nc] ->Fill(ZRecRap,totalWeight);
        	histoRapPRecGen_MomRes[lepType][nc]->Fill(ZRecRap,ZGenRap,totalWeight);
              } else {
        	histoRapMRecDY_MomRes[lepType][nc] ->Fill(ZRecRap,totalWeight);
        	histoRapMRecGen_MomRes[lepType][nc]->Fill(ZRecRap,ZGenRap,totalWeight);
              }
            }
            if     (ZRecRap < 0.5) {histoPtRap0RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap0RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
            else if(ZRecRap < 1.0) {histoPtRap1RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap1RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
            else if(ZRecRap < 1.5) {histoPtRap2RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap2RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
            else if(ZRecRap < 2.0) {histoPtRap3RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap3RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
            else if(ZRecRap < 2.4) {histoPtRap4RecGen_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],ZGenPt,totalWeight);histoPtRap4RecDY_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);}
	  }
	  else if(theCategory == 0){
            histoPtRecDA_MomRes[lepType][nc] ->Fill(ZRecSystPt[nc],totalWeight);
            if(ZRecRap < 2.4) {
              histoRapRecDA_MomRes[lepType][nc]->Fill(ZRecRap,totalWeight);
              if(thePandaFlat.looseLep1PdgId < 0){
        	histoRapPRecDA_MomRes[lepType][nc]->Fill(ZRecRap,totalWeight);
              } else {
        	histoRapMRecDA_MomRes[lepType][nc]->Fill(ZRecRap,totalWeight);
              }
            }
            if     (ZRecRap < 0.5) histoPtRap0RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 1.0) histoPtRap1RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 1.5) histoPtRap2RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 2.0) histoPtRap3RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 2.4) histoPtRap4RecDA_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
	  }
	  else if(theCategory == 5){
            for(int theLepType = 0; theLepType<2; theLepType++) {
              double theKeff = k_eff;
              if(theLepType == 1) theKeff = 1.0/k_eff;
              histoPtRecEM_MomRes[theLepType][nc] ->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              if(ZRecRap < 2.4) {
        	histoRapRecEM_MomRes[theLepType][nc]->Fill(ZRecRap,totalWeight*theKeff);
        	if(thePandaFlat.looseLep1PdgId < 0){
        	  histoRapPRecEM_MomRes[theLepType][nc]->Fill(ZRecRap,totalWeight*theKeff);
        	} else {
        	  histoRapMRecEM_MomRes[theLepType][nc]->Fill(ZRecRap,totalWeight*theKeff);
        	}
              }
              if     (ZRecRap < 0.5) histoPtRap0RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              else if(ZRecRap < 1.0) histoPtRap1RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              else if(ZRecRap < 1.5) histoPtRap2RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              else if(ZRecRap < 2.0) histoPtRap3RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
              else if(ZRecRap < 2.4) histoPtRap4RecEM_MomRes[theLepType][nc]->Fill(ZRecSystPt[nc],totalWeight*theKeff);
            }
	  }
	  else {
            histoPtRecVV_MomRes[lepType][nc] ->Fill(ZRecSystPt[nc],totalWeight);
            if(ZRecRap < 2.4) {
              histoRapRecVV_MomRes[lepType][nc]->Fill(ZRecRap,totalWeight);
              if(thePandaFlat.looseLep1PdgId < 0){
        	histoRapPRecVV_MomRes[lepType][nc]->Fill(ZRecRap,totalWeight);
              } else {
        	histoRapMRecVV_MomRes[lepType][nc]->Fill(ZRecRap,totalWeight);
              }
            }
            if     (ZRecRap < 0.5) histoPtRap0RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 1.0) histoPtRap1RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 1.5) histoPtRap2RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 2.0) histoPtRap3RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
            else if(ZRecRap < 2.4) histoPtRap4RecVV_MomRes[lepType][nc]->Fill(ZRecSystPt[nc],totalWeight);
	  }
	} // end passSystSel[nc]
      } // end loop over passSystSel
    } // end event loop
  } // end samples loop

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

  // RapP
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDRapP(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoRapPRecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoRapPRecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoRapPRecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoRapPRecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoRapPRecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoRapPRecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoRapPRecVV[ntype]->GetSumOfWeights(),
           histoRapPRecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoRapPRecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoRapPRecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoRapPRecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoRapPRecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoRapPRecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoRapPRecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRap+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoRapPRecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoRapPRecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoRapPRecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoRapPRecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoRapPRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapPRecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoRapPRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapPRecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoRapPRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapPRecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoRapPRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapPRecDY[ntype]->GetBinContent(nb));
      }
      if(histoRapPRecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoRapPRecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoRapPRecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoRapPRecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoRapPRecVV_QCD[ntype]->SetBinContent(nb, histoRapPRecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoRapPRecDY_QCD[ntype]->SetBinContent(nb, histoRapPRecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
    }
  }

  // RapM
  for(int ntype=0; ntype<2; ntype++){
    printf("QCDRapM(%d): (%f/%f/%f/%f/%f/%f->%f) (%f/%f/%f/%f/%f/%f->%f)\n",ntype,
           histoRapMRecVV_QCDPart[ntype][0]->GetSumOfWeights(),histoRapMRecVV_QCDPart[ntype][1]->GetSumOfWeights(),histoRapMRecVV_QCDPart[ntype][2]->GetSumOfWeights(),
           histoRapMRecVV_QCDPart[ntype][3]->GetSumOfWeights(),histoRapMRecVV_QCDPart[ntype][4]->GetSumOfWeights(),histoRapMRecVV_QCDPart[ntype][5]->GetSumOfWeights(),histoRapMRecVV[ntype]->GetSumOfWeights(),
           histoRapMRecDY_QCDPart[ntype][0]->GetSumOfWeights(),histoRapMRecDY_QCDPart[ntype][1]->GetSumOfWeights(),histoRapMRecDY_QCDPart[ntype][2]->GetSumOfWeights(),
           histoRapMRecDY_QCDPart[ntype][3]->GetSumOfWeights(),histoRapMRecDY_QCDPart[ntype][4]->GetSumOfWeights(),histoRapMRecDY_QCDPart[ntype][5]->GetSumOfWeights(),histoRapMRecDY[ntype]->GetSumOfWeights());
  }
  for(int ntype=0; ntype<2; ntype++){
    for(int nb=1; nb<=nBinRap+1; nb++){
      // QCD study
      double systQCDScale[2] = {TMath::Abs(histoRapMRecVV_QCDPart[ntype][0]->GetBinContent(nb)-histoRapMRecVV[ntype]->GetBinContent(nb)),
  			        TMath::Abs(histoRapMRecDY_QCDPart[ntype][0]->GetBinContent(nb)-histoRapMRecDY[ntype]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histoRapMRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapMRecVV[ntype]->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histoRapMRecVV_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapMRecVV[ntype]->GetBinContent(nb));
        if(TMath::Abs(histoRapMRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapMRecDY[ntype]->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histoRapMRecDY_QCDPart[ntype][nqcd]->GetBinContent(nb)-histoRapMRecDY[ntype]->GetBinContent(nb));
      }
      if(histoRapMRecVV[ntype]->GetBinContent(nb) > 0) systQCDScale[0] = 1.0+systQCDScale[0]/histoRapMRecVV[ntype]->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histoRapMRecDY[ntype]->GetBinContent(nb) > 0) systQCDScale[1] = 1.0+systQCDScale[1]/histoRapMRecDY[ntype]->GetBinContent(nb); else systQCDScale[1] = 1;

      histoRapMRecVV_QCD[ntype]->SetBinContent(nb, histoRapMRecVV[ntype]->GetBinContent(nb)*systQCDScale[0]);
      histoRapMRecDY_QCD[ntype]->SetBinContent(nb, histoRapMRecDY[ntype]->GetBinContent(nb)*systQCDScale[1]);
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
  for(int i=0; i<2; i++){
    histoPtRecGen[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRecGen_MomRes[i][nj]->Write();
    histoPtRecDA[i]->Write();
    histoPtRecDY[i]->Write();
    histoPtRecEM[i]->Write();
    histoPtRecVV[i]->Write();
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
  sprintf(output,"histoDY%dzllRapRecGen.root",whichDY);
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoRapRecGen[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoRapRecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapRecGen_MomRes[i][nj]->Write();
    histoRapRecDA[i]->Write();
    histoRapRecDY[i]->Write();
    histoRapRecEM[i]->Write();
    histoRapRecVV[i]->Write();
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
  sprintf(output,"histoDY%dzllRapPRecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoRapPRecGen[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoRapPRecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapPRecGen_MomRes[i][nj]->Write();
    histoRapPRecDA[i]->Write();
    histoRapPRecDY[i]->Write();
    histoRapPRecEM[i]->Write();
    histoRapPRecVV[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoRapPRecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoRapPRecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapPRecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapPRecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapPRecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapPRecVV_MomRes[i][nj]->Write();
    histoRapPRecVV_PDF[i]->Write();
    histoRapPRecDY_PDF[i]->Write();
    histoRapPRecVV_QCD[i]->Write();
    histoRapPRecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllRapMRecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoRapMRecGen[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoRapMRecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapMRecGen_MomRes[i][nj]->Write();
    histoRapMRecDA[i]->Write();
    histoRapMRecDY[i]->Write();
    histoRapMRecEM[i]->Write();
    histoRapMRecVV[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoRapMRecDY_LepEff[i][nj]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoRapMRecVV_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapMRecDA_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapMRecDY_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapMRecEM_MomRes[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoRapMRecVV_MomRes[i][nj]->Write();
    histoRapMRecVV_PDF[i]->Write();
    histoRapMRecDY_PDF[i]->Write();
    histoRapMRecVV_QCD[i]->Write();
    histoRapMRecDY_QCD[i]->Write();
  }
  outFilePlots->Close();
  }

  {
  sprintf(output,"histoDY%dzllPtRap0RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRap0RecGen[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap0RecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap0RecGen_MomRes[i][nj]->Write();
    histoPtRap0RecDA[i]->Write();
    histoPtRap0RecDY[i]->Write();
    histoPtRap0RecEM[i]->Write();
    histoPtRap0RecVV[i]->Write();
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
  sprintf(output,"histoDY%dzllPtRap1RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRap1RecGen[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap1RecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap1RecGen_MomRes[i][nj]->Write();
    histoPtRap1RecDA[i]->Write();
    histoPtRap1RecDY[i]->Write();
    histoPtRap1RecEM[i]->Write();
    histoPtRap1RecVV[i]->Write();
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
  sprintf(output,"histoDY%dzllPtRap2RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRap2RecGen[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap2RecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap2RecGen_MomRes[i][nj]->Write();
    histoPtRap2RecDA[i]->Write();
    histoPtRap2RecDY[i]->Write();
    histoPtRap2RecEM[i]->Write();
    histoPtRap2RecVV[i]->Write();
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
  sprintf(output,"histoDY%dzllPtRap3RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRap3RecGen[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap3RecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap3RecGen_MomRes[i][nj]->Write();
    histoPtRap3RecDA[i]->Write();
    histoPtRap3RecDY[i]->Write();
    histoPtRap3RecEM[i]->Write();
    histoPtRap3RecVV[i]->Write();
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
  sprintf(output,"histoDY%dzllPtRap4RecGen.root",whichDY); 
  TFile* outFilePlots = new TFile(output,"recreate");
  outFilePlots->cd();
  for(int i=0; i<2; i++){
    histoPtRap4RecGen[i]->Write();
    for(int nj=0; nj<nEffNuisances; nj++) histoPtRap4RecGen_LepEff[i][nj]->Write();
    for(int nj=0; nj<nMomNuisances; nj++) histoPtRap4RecGen_MomRes[i][nj]->Write();
    histoPtRap4RecDA[i]->Write();
    histoPtRap4RecDY[i]->Write();
    histoPtRap4RecEM[i]->Write();
    histoPtRap4RecVV[i]->Write();
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
