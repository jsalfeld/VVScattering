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

#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"

#include "MitAnalysisRunII/macros/factors.h"

enum selType                     { SIGSEL, nSelTypes};
TString selTypeName[nSelTypes]= { "SIGSEL"};

enum systType                     {METUP=0, METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"METUP","METDOWN"};

double mcPrescale = 1.0;
//const TString typeLepSel = "medium";
//const TString type3rdLepSel = "medium";
bool usePureMC = false;

void wzAnalysis(
 double minMass =  76.1876,
 double maxMass = 106.1876,
 bool applyBtagging = true,
 TString typeLepSel = "medium",
 TString type3rdLepSel = "default"
 ){

  Int_t period = 1;
  TString filesPath  = "/scratch5/ceballos/ntuples_weights/";
  Double_t lumi = 0.0715;
  if(period == 1) lumi = 2.2;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  if      (period==1){
  puPath = "/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src/MitAnalysisRunII/data/puWeights_13TeV_25ns.root";

  infilenamev.push_back(Form("%sdata_AOD_Run2015C1_25ns.root",filesPath.Data()));												  infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D3_25ns.root",filesPath.Data()));												  infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D4_25ns.root",filesPath.Data()));												  infilecatv.push_back(0);
  if(usePureMC == true){
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));					          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		  		          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));         infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  		  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  	          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); 		          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(1);
  /////infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); 		  infilecatv.push_back(1);
  /////infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(1);
  }
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(2);

  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));                          infilecatv.push_back(3);

  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				          infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToZZTo2e2mu_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToZZTo2e2tau_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToZZTo2mu2tau_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToZZTo4e_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToZZTo4mu_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToZZTo4tau_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));        infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));         infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(5);

  }
  else {assert(0);}
  
  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
 
  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 40;
  const int histBins = 6;
  const int allStates = 5;
  TH1D* histo[allStates][allPlots][histBins];
  TString processName[histBins] = {"..Data", "....EM", "Zgamma", "....WZ", "...ZZ", "...VVV"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  2) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  3 && thePlot <=  3) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 150.0;}
    else if(thePlot >=  4 && thePlot <=  4) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  5 && thePlot <=  5) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >=  6 && thePlot <=  6) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot =   4.0;}
    else if(thePlot >=  8 && thePlot <= 11) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 12 && thePlot <= 12) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 13 && thePlot <= 13) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot =  40.0;}
    else if(thePlot >= 14 && thePlot <= 14) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot >= 15 && thePlot <= 15) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >= 16 && thePlot <= 17) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 18 && thePlot <= 18) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 20 && thePlot <= 20) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot >= 21 && thePlot <= 22) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >= 23 && thePlot <= 23) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >= 24 && thePlot <= 24) {nBinPlot =  50; xminPlot =-0.5; xmaxPlot =  49.5;}
    else if(thePlot >= 25 && thePlot <= 31) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 32 && thePlot <= 32) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot >= 33 && thePlot <= 36) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) {
      for(int j=0; j<allStates; j++) histo[j][thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    }
    histos->Reset();histos->Clear();
  }

  TString ECMsb  = "13TeV2015";
  const int nBinMVA = 4; Float_t xbins[nBinMVA+1] = {0, 1, 2, 3, 4};
  //const int nBinMVA = 1; Float_t xbins[nBinMVA+1] = {0, 1};
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();

  TH1D *histo_Data   = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_Zg     = (TH1D*) histoMVA->Clone("histo_Zg"); 
  TH1D *histo_VVV    = (TH1D*) histoMVA->Clone("histo_VVV");	 
  TH1D *histo_WZ     = (TH1D*) histoMVA->Clone("histo_WZ");	 
  TH1D *histo_ZZ     = (TH1D*) histoMVA->Clone("histo_ZZ");	 
  TH1D *histo_FakeM  = (TH1D*) histoMVA->Clone("histo_FakeM");	 
  TH1D *histo_FakeE  = (TH1D*) histoMVA->Clone("histo_FakeE");	 

  char finalStateName[2],effMName[10],effEName[10],momMName[10],momEName[10];
  sprintf(effMName,"CMS_eff_m");sprintf(momMName,"CMS_scale_m");
  sprintf(effEName,"CMS_eff_e");sprintf(momEName,"CMS_scale_e");
  sprintf(finalStateName,"3l");

  TH1D* histo_Zg_CMS_MVAZHStatBoundingUp           = new TH1D( Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVAZHStatBoundingDown         = new TH1D( Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp         = new TH1D( Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown       = new TH1D( Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingUp           = new TH1D( Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingDown         = new TH1D( Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingUp           = new TH1D( Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingDown         = new TH1D( Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingDown->Sumw2();
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingUp     = new TH1D( Form("histo_FakeM_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_FakeM_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeM_CMS_MVAFakeMStatBoundingUp  ->Sumw2();
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingDown   = new TH1D( Form("histo_FakeM_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_FakeM_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeM_CMS_MVAFakeMStatBoundingDown->Sumw2();
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingUp     = new TH1D( Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->Sumw2();
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingDown   = new TH1D( Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingDown->Sumw2();

  TH1D* histo_Diff = new TH1D("dummy", "dummy",1000,-1,1); histo_Diff->Sumw2();

  TH1D* histo_Zg_CMS_QCDScaleBounding[6];
  TH1D* histo_VVV_CMS_QCDScaleBounding[6];
  TH1D* histo_WZ_CMS_QCDScaleBounding[6];
  TH1D* histo_ZZ_CMS_QCDScaleBounding[6];
  for(int nb=0; nb<6; nb++){
    histo_Zg_CMS_QCDScaleBounding[nb]      = new TH1D(Form("histo_Zg_QCDScale_f%d",nb),      Form("histo_Zg_QCDScale_f%d",nb),nBinMVA, xbins);      histo_Zg_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_VVV_CMS_QCDScaleBounding[nb]     = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),     Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WZ_CMS_QCDScaleBounding[nb]	   = new TH1D(Form("histo_WZ_QCDScale_f%d",nb),      Form("histo_WZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_WZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ZZ_CMS_QCDScaleBounding[nb]	   = new TH1D(Form("histo_ZZ_QCDScale_f%d",nb),      Form("histo_ZZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_ZZ_CMS_QCDScaleBounding[nb]->Sumw2();
  }
  TH1D* histo_Zg_CMS_PDFBounding[102];
  TH1D* histo_VVV_CMS_PDFBounding[102];
  TH1D* histo_WZ_CMS_PDFBounding[102];
  TH1D* histo_ZZ_CMS_PDFBounding[102];
  for(int nb=0; nb<102; nb++){
    histo_Zg_CMS_PDFBounding[nb]      = new TH1D(Form("histo_Zg_PDF_f%d",nb),      Form("histo_Zg_PDF_f%d",nb),     nBinMVA, xbins); histo_Zg_CMS_PDFBounding[nb]->Sumw2();
    histo_VVV_CMS_PDFBounding[nb]     = new TH1D(Form("histo_VVV_PDF_f%d",nb),     Form("histo_VVV_PDF_f%d",nb),    nBinMVA, xbins); histo_VVV_CMS_PDFBounding[nb]->Sumw2();
    histo_WZ_CMS_PDFBounding[nb]      = new TH1D(Form("histo_WZ_PDF_f%d",nb),      Form("histo_WZ_PDF_f%d",nb),     nBinMVA, xbins); histo_WZ_CMS_PDFBounding[nb]->Sumw2();
    histo_ZZ_CMS_PDFBounding[nb]      = new TH1D(Form("histo_ZZ_PDF_f%d",nb),      Form("histo_ZZ_PDF_f%d",nb),     nBinMVA, xbins); histo_ZZ_CMS_PDFBounding[nb]->Sumw2();
  }

  TH1D* histo_Zg_CMS_MVAZHStatBoundingBinUp[nBinMVA];
  TH1D* histo_Zg_CMS_MVAZHStatBoundingBinDown[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinUp[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinDown[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinUp[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinDown[nBinMVA];
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[nBinMVA];
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[nBinMVA];
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nBinMVA];
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++){
    histo_Zg_CMS_MVAZHStatBoundingBinUp[nb]         = new TH1D(Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingBinUp[nb]   ->Sumw2();
    histo_Zg_CMS_MVAZHStatBoundingBinDown[nb]       = new TH1D(Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingBinDown[nb] ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	    = new TH1D(Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	    = new TH1D(Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	    = new TH1D(Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	    = new TH1D(Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	    = new TH1D(Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]      ->Sumw2();
    histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[nb]   = new TH1D(Form("histo_FakeM_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_FakeM_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[nb]       ->Sumw2();
    histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[nb] = new TH1D(Form("histo_FakeM_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_FakeM_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[nb]	 ->Sumw2();
    histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nb]   = new TH1D(Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nb]       ->Sumw2();
    histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nb] = new TH1D(Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nb]	 ->Sumw2();
  }

  TH1D* histo_Zg_CMS_MVALepEffMBoundingUp        = new TH1D( Form("histo_Zg_%sUp",effMName)  , Form("histo_Zg_%sUp",effMName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVALepEffMBoundingDown      = new TH1D( Form("histo_Zg_%sDown",effMName), Form("histo_Zg_%sDown",effMName), nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effMName)  , Form("histo_VVV_%sUp",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effMName), Form("histo_VVV_%sDown",effMName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_WZ_%sUp",effMName)  , Form("histo_WZ_%sUp",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_WZ_%sDown",effMName), Form("histo_WZ_%sDown",effMName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_ZZ_%sUp",effMName)  , Form("histo_ZZ_%sUp",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_ZZ_%sDown",effMName), Form("histo_ZZ_%sDown",effMName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_MVALepEffMBoundingAvg       = new TH1D( Form("histo_Zg_%sAvg",effMName)  , Form("histo_Zg_%sAvg",effMName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effMName)  , Form("histo_VVV_%sAvg",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_WZ_%sAvg",effMName)  , Form("histo_WZ_%sAvg",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_ZZ_%sAvg",effMName)  , Form("histo_ZZ_%sAvg",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();

  TH1D* histo_Zg_CMS_MVALepEffEBoundingUp        = new TH1D( Form("histo_Zg_%sUp",effEName)  , Form("histo_Zg_%sUp",effEName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVALepEffEBoundingDown      = new TH1D( Form("histo_Zg_%sDown",effEName), Form("histo_Zg_%sDown",effEName), nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effEName)  , Form("histo_VVV_%sUp",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effEName), Form("histo_VVV_%sDown",effEName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_WZ_%sUp",effEName)  , Form("histo_WZ_%sUp",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_WZ_%sDown",effEName), Form("histo_WZ_%sDown",effEName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_ZZ_%sUp",effEName)  , Form("histo_ZZ_%sUp",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_ZZ_%sDown",effEName), Form("histo_ZZ_%sDown",effEName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_MVALepEffEBoundingAvg       = new TH1D( Form("histo_Zg_%sAvg",effEName)  , Form("histo_Zg_%sAvg",effEName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effEName)  , Form("histo_VVV_%sAvg",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_WZ_%sAvg",effEName)  , Form("histo_WZ_%sAvg",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_ZZ_%sAvg",effEName)  , Form("histo_ZZ_%sAvg",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();

  TH1D* histo_Zg_CMS_MVAMETBoundingUp           = new TH1D( Form("histo_Zg_CMS_scale_metUp")  , Form("histo_Zg_CMS_scale_metUp")  , nBinMVA, xbins); histo_Zg_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVAMETBoundingDown         = new TH1D( Form("histo_Zg_CMS_scale_metDown"), Form("histo_Zg_CMS_scale_metDown"), nBinMVA, xbins); histo_Zg_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingUp   	= new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingDown 	= new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown"), nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_WZ_CMS_scale_metUp")  , Form("histo_WZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_WZ_CMS_scale_metDown"), Form("histo_WZ_CMS_scale_metDown"), nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_scale_metUp")  , Form("histo_ZZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_scale_metDown"), Form("histo_ZZ_CMS_scale_metDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingDown->Sumw2();

  unsigned int numberOfLeptons = 3;

  double totalEventsProcess[50];
  std::vector<double> sumEventsProcess(infilenamev.size(), 0.0);

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");
    TTree *the_PDF_tree   = (TTree*)the_input_file.FindObjectAny("pdfReweight");

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareJets eventJets;
    eventJets.setBranchAddresses(the_input_tree);

    BareLeptons eventLeptons;
    eventLeptons.setBranchAddresses(the_input_tree);

    BareTaus eventTaus;
    eventTaus.setBranchAddresses(the_input_tree);

    BareMet eventMet;
    eventMet.SetExtend();
    eventMet.setBranchAddresses(the_input_tree);

    BareTrigger eventTrigger;
    eventTrigger.setBranchAddresses(the_input_tree);

    BareVertex eventVertex;
    eventVertex.setBranchAddresses(the_input_tree);

    BareMonteCarlo eventMonteCarlo;
    eventMonteCarlo.setBranchAddresses(the_input_tree);

    TNamed *triggerNames = (TNamed*)the_input_file.FindObjectAny("triggerNames");
    char **tokens;
    size_t numtokens;
    tokens = strsplit(triggerNames->GetTitle(), ",", &numtokens);
    if(infilecatv[ifile] == 0){
      for (int i = 0; i < (int)numtokens; i++) {
        printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
      }
    }
    else {
      printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());
    }

    char weightDef[256];
    int initPDFTag = -1;
    if(the_PDF_tree) {
      the_PDF_tree->SetBranchAddress("weightDef", &weightDef);
      for (int i=0; i<the_PDF_tree->GetEntries(); ++i) {
        the_PDF_tree->GetEntry(i);
        char **tokensPDF;
        size_t numtokensPDF;
        tokensPDF = strsplit(weightDef, " = ", &numtokensPDF);
        for (int k = 0; k < (int)numtokensPDF; k++) if(strcmp(tokensPDF[k],"292201") == 0||
	                                               strcmp(tokensPDF[k],"292001") == 0||
						       strcmp(tokensPDF[k],"260001") == 0) {initPDFTag = i; break;}
	if(initPDFTag != -1) break;
      }
    }
    if(infilecatv[ifile] != 0 && initPDFTag == -1) {
      printf("PDFTAG PROBLEM\n");
      if(the_PDF_tree) {
        printf("PDFTree Entries: %d\n",(int)the_PDF_tree->GetEntries());
        for (int i=0; i<the_PDF_tree->GetEntries(); ++i) {
          the_PDF_tree->GetEntry(i);
	  printf("PDF(%d): %s\n",i,weightDef);
        }
      }
      else {
        printf("PDFTree not available\n");
      }
      //return;
    }

    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_input_tree->GetEntry(i);

      if(i%1000000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      Bool_t passFilter[11] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
      for (int nt = 0; nt <(int)numtokens; nt++) {
        if((*eventTrigger.triggerFired)[nt] == 0) continue;
        if((strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
           (strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
           (strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*") 	    == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")	    == 0) ||
           (strcmp(tokens[nt],"HLT_IsoMu27_v*") 				    == 0) ||
           (strcmp(tokens[nt],"HLT_IsoMu20_v*") 				    == 0) ||
           (strcmp(tokens[nt],"HLT_IsoTkMu20_v*")				    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")	    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele23_WPLoose_Gsf_v*")			    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele22_eta2p1_WP75_Gsf_v*")			    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele27_WPLoose_Gsf_v*")			    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele27_WP85_Gsf_v*")  			    == 0)
           ) passFilter[1] = kTRUE;
      }

      if(passFilter[0] == kFALSE) continue;
      if(passFilter[1] == kFALSE) continue;
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep); if(numberOfLeptons == 4) goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}
      }
      if(idLep.size()!=idTight.size()) {assert(1); return;}
      if(idLep.size()==numberOfLeptons) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;

      if(goodIsTight == idTight.size()) passFilter[3] = kTRUE;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 10 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() <= 10) continue;

      double dPhiLepMETMin = 999.;
      int signQ = 0;
      double systTotLep[2] = {1.0, 1.0}; // m/e
      for(unsigned nl=0; nl<idLep.size(); nl++){
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) systTotLep[0] = systTotLep[0] * 1.02;
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 11) systTotLep[1] = systTotLep[1] * 1.02;
        signQ = signQ + (int)(*eventLeptons.pdgId)[idLep[nl]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]);
        if(dPhiLepMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0]))))
           dPhiLepMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));      
      }
      double minMET  = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      double minPMET = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      if(dPhiLepMETMin < TMath::Pi()/2) minPMET = minPMET * sin(dPhiLepMETMin);

      passFilter[4] = TMath::Abs(signQ) == 1;
      if(passFilter[4] == kFALSE) continue;

      double minMassll = 999.0;
      double minMassZ = 999.0;
      double mass3l = 0.0;
      double deltaRllMin = 999.0;
      int type3l = 0;
      int tagZ[3] = {-1,-1,-1};
      for(unsigned nl0=0; nl0<idLep.size()-1; nl0++){
        for(unsigned nl1=nl0+1; nl1<idLep.size(); nl1++){
          double deltaRllAux = ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl0]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl1]]));
          if(deltaRllAux < deltaRllMin) deltaRllMin = deltaRllAux;

	  if((int)(*eventLeptons.pdgId)[idLep[nl0]] * (int)(*eventLeptons.pdgId)[idLep[nl1]] > 0) continue;
          TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl1])) ) ));

	  if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl1]]) &&
	     TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassZ-91.1876)) {
	     minMassZ = dilepAux.M();tagZ[0]=nl0;tagZ[1]=nl1;
	     if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]]) == 13) type3l = 0;
	     else                                                         type3l = 1;
	  }

	  if(minMassll > dilepAux.M()) minMassll = dilepAux.M();
        }
      }

      vector<int> idJet;
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      double theHT = 0;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;        

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(dPhiJetMET   == -1) dPhiJetMET   = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 10 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) continue;

        theHT = theHT + ((TLorentzVector*)(*eventJets.p4)[nj])->Pt();
	idJet.push_back(nj);
      }

      for(unsigned nl0=0; nl0<idLep.size(); nl0++){
        if((int)nl0==tagZ[0]||(int)nl0==tagZ[1]) continue;
        tagZ[2] = nl0;
        break;
      }
      
      if(tagZ[0] == -1 || tagZ[1] == -1 || tagZ[2] == -1) continue;
      
      bool tight3rdLepId = true;
      if(idTight[tagZ[2]] == 1){
        tight3rdLepId = selectIdIsoCut(type3rdLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Eta()),(double)(*eventLeptons.iso)[idLep[tagZ[2]]],(int)(*eventLeptons.selBits)[idLep[tagZ[2]]]);
      }

      if(tight3rdLepId == false) continue;

      if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]) == 13) type3l += 0;
      else							       type3l += 2;
      TLorentzVector trilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + 
        		      ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) + 
        		      ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[2])) ) ));
      mass3l = trilep.M();

      passFilter[ 5] = minMassll > 4; //deltaRllMin > 0.1;
      passFilter[ 6] = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 30;
      passFilter[ 7] = minMassZ > minMass && minMassZ < maxMass;
      passFilter[ 8] = ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt() > 20;
      passFilter[ 9] = mass3l > 100;
      passFilter[10] = true;
      if(applyBtagging) passFilter[10] = bDiscrMax < 0.97;

      if(passFilter[7]==kTRUE && (tagZ[0] == tagZ[1] || tagZ[0] == tagZ[2] || tagZ[1] == tagZ[2])) {printf("ZPROBLEM!\n");assert(0);return;}

      bool passNMinusOne[6] = {                 passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10],
                               passFilter[5] &&                  passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10],
                               passFilter[5] && passFilter[6] &&                  passFilter[8] && passFilter[9] && passFilter[10],
                               passFilter[5] && passFilter[6] && passFilter[7] &&                  passFilter[9] && passFilter[10],
                               passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8]                  && passFilter[10],
			       passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9]};

      bool passAllCuts[1] = {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};

      bool controlSel[3] = {passFilter[5] && passFilter[6] && !passFilter[7] && passFilter[9] && bDiscrMax > 0.605 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20,
                            passFilter[5] &&                  !passFilter[7] && passFilter[9] && bDiscrMax > 0.605 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20,
			    passFilter[5] && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() > 20};

      bool passVBF = passAllCuts[0] && idJet.size() >= 2 &&
                     (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M() > 500;

      double deltaPhiLeptonMet = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtLN = TMath::Sqrt(2.0*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiLeptonMet)));

      double deltaPhiTriLeptonMet = TMath::Abs(trilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtEvent = TMath::Sqrt(2.0*trilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiTriLeptonMet)));

     bool passSystCuts[nSystTypes] = {
          passFilter[5] && (double)(*eventMet.ptJESUP)[0]   > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10],
	  passFilter[5] && (double)(*eventMet.ptJESDOWN)[0] > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]
     };

      // begin event weighting
      vector<bool> isGenDupl;
      int numberQuarks[2] = {0,0};
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 4 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[0]++;
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 5 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[1]++;
        isGenDupl.push_back(0);
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) isGenDupl[ngen0] = 1;
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) continue;
        for(int ngen1=ngen0+1; ngen1<eventMonteCarlo.p4->GetEntriesFast(); ngen1++) {
          if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen1])) < 0.02) {
	    isGenDupl[ngen0] = 1;
	    break;
	  }
        }
      }
      vector<int> isGenLep; unsigned int goodIsGenLep = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        bool isGenLepton = false;
        for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	  if(isGenDupl[ngen] == 1) continue;
          if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.1) {
	    isGenLepton = true;
	    break;
	  }
	}
	if(isGenLepton == true) {isGenLep.push_back(1); goodIsGenLep++;}
	else                    {isGenLep.push_back(0);}
      }

      // luminosity
      double theLumi  = 1.0; if(infilecatv[ifile] != 0) theLumi  = lumi;
      // pile-up
      double puWeight = 1.0; if(infilecatv[ifile] != 0) puWeight = nPUScaleFactor(fhDPU, (double)eventVertex.npv);
      // lepton efficiency
      double effSF = 1.0;
      if(infilecatv[ifile] != 0){
        for(unsigned int nl=0; nl<idLep.size(); nl++){
	  if(tagZ[2] != (int)nl)
          effSF = effSF * effScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
          else
	  effSF = effSF * effScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,type3rdLepSel.Data());
        }
      }

      // fake rate
      unsigned int typeFakeLepton[2] = {0,0};
      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false){
        if     ((infilecatv[ifile] == 0 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add Z+jets from data
          for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    if(tagZ[2] != (int)nl)
	    fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    else
	    fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,type3rdLepSel.Data());
	    theCategory = 1;
	    if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) typeFakeLepton[0]++;
	    else                                                        typeFakeLepton[1]++;
	  }
          if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-3) fakeSF = -1.0 * fakeSF; // triple fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) fakeSF = +1.0 * fakeSF; // double fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) fakeSF = -1.0 * fakeSF; // single fake, MC
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-3) fakeSF = +1.0 * fakeSF; // triple fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) fakeSF = -1.0 * fakeSF; // double fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) fakeSF = +1.0 * fakeSF; // single fake, data
	  if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) {
	    if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 13){
	      if     (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479 && 
	             ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() < 15) fakeSF = 2.6 * fakeSF;
	      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479 && 
	             ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() < 20) fakeSF = 1.6 * fakeSF;
	      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479 && 
	             ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() < 15) fakeSF = 3.7 * fakeSF;
	      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479 && 
	             ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() < 20) fakeSF = 2.0 * fakeSF;
            }
          }
	}
        else if(infilecatv[ifile] != 0 && infilecatv[ifile] != 2 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] == 2 && goodIsTight != idTight.size()){ // remove Z+gamma, fakeable objects
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
	  fakeSF = 1.0;
        }
        else if(infilecatv[ifile] == 0 || infilecatv[ifile] == 2){ // data or Z+gamma with all good leptons
	  fakeSF = 1.0;
        }
	else {
	  printf("PROBLEM: %d %d %d %d %d\n",infilecatv[ifile],goodIsGenLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
	  assert(0);
	}
      }
      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      //double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale;
      if(totalWeight == 0) continue;
      // end event weighting
      if((infilecatv[ifile] != 0 || theCategory == 0) && passAllCuts[0]) sumEventsProcess[ifile] += totalWeight;

      for(int thePlot=0; thePlot<allPlots; thePlot++){
	double theVar = 0.0;
	bool makePlot = false;
 	if     (thePlot ==  0 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(mtLN,199.999);}
	else if(thePlot ==  1 && passNMinusOne[0])                       {makePlot = true;theVar = TMath::Min(minMassll,199.999);}
	else if(thePlot ==  2 && passNMinusOne[1])                       {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot ==  3 && passNMinusOne[2])                       {makePlot = true;theVar = TMath::Max(TMath::Min(minMassZ,149.999),50.001);}
	else if(thePlot ==  4 && passNMinusOne[3])                       {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt(),199.999);}
	else if(thePlot ==  5 && passNMinusOne[4])                       {makePlot = true;theVar = TMath::Min(mass3l,399.999);}
	else if(thePlot ==  6 && passNMinusOne[5])                       {makePlot = true;theVar = TMath::Min(bDiscrMax,0.999);}
	else if(thePlot ==  7 && passNMinusOne[0])                       {makePlot = true;theVar = TMath::Min(deltaRllMin,3.999);}
	else if(thePlot ==  8 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(minPMET,199.999);}
	else if(thePlot ==  9 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Pt(),199.999);}
	else if(thePlot == 10 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Pt(),199.999);}
	else if(thePlot == 11 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(minMET,199.999);}
	else if(thePlot == 12 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 13 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min((double)eventEvent.rho,39.999);}
	else if(thePlot == 14 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	else if(thePlot == 15 && passAllCuts[0])                         {makePlot = true;theVar = (double)type3l;}
	else if(thePlot == 16 && passAllCuts[0])                         {makePlot = true;theVar = dPhiJetMET*180/TMath::Pi();}
	else if(thePlot == 17 && passAllCuts[0])                         {makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	else if(thePlot == 18 && controlSel[0])                          {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 19 && controlSel[1])                          {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 20 && controlSel[2])                          {makePlot = true;theVar = TMath::Max(TMath::Min(mass3l,249.999),50.001);}
	else if(thePlot == 21 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(mtEvent,499.999);}
	else if(thePlot == 22 && passVBF)                                {makePlot = true;theVar = TMath::Min(mtEvent,499.999);}
	else if(thePlot == 23 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(theHT,399.999);}
	else if(thePlot == 24 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min((double)(numberQuarks[0]+10*numberQuarks[1]),49.499);}
	else if(thePlot == 25 && passNMinusOne[1])                       {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 26 && passNMinusOne[3])	                 {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt(),199.999);}
	else if(thePlot == 27 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(mtLN,199.999);}
	else if(thePlot == 28 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 11 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 29 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 13 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 30 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 11 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 31 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 13 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 32 && controlSel[2])                          {makePlot = true;theVar = TMath::Max(TMath::Min(mass3l,249.999),50.001);}
	else if(thePlot == 33 && passNMinusOne[5] && numberQuarks[1] == 0                  ) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 34 && passNMinusOne[5] && numberQuarks[1] == 0 && passFilter[10]) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 35 && passNMinusOne[5] && numberQuarks[1]  > 0                  ) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 36 && passNMinusOne[5] && numberQuarks[1]  > 0 && passFilter[10]) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}

	if(makePlot) histo[     4][thePlot][theCategory]->Fill(theVar,totalWeight);
	if(makePlot) histo[type3l][thePlot][theCategory]->Fill(theVar,totalWeight);
      }
      
      if(1) {
	double MVAVar = (double)type3l;
	//double MVAVar = 0;
        if     (theCategory == 0){
	  if(passAllCuts[SIGSEL]) histo_Data->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 1){
	  if(passAllCuts[SIGSEL]) {
	    if(typeFakeLepton[0]+typeFakeLepton[1] != goodIsTight == idTight.size()) {printf("PROBLEMFake %d %d %d %d\n",typeFakeLepton[0],typeFakeLepton[1],goodIsTight,(int)idTight.size()); return;}
	    if     (typeFakeLepton[0] > typeFakeLepton[1]) histo_FakeM->Fill(MVAVar,totalWeight);
	    else if(typeFakeLepton[0] < typeFakeLepton[1]) histo_FakeE->Fill(MVAVar,totalWeight);
	    else {
	      histo_FakeM->Fill(MVAVar,totalWeight/2.);
	      histo_FakeE->Fill(MVAVar,totalWeight/2.);
	    }
	  }
        }
        else if(theCategory == 2){
	  if(passAllCuts[SIGSEL]) {
	     histo_Zg->Fill(MVAVar,totalWeight);

	     histo_Zg_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_Zg_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_Zg_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_Zg_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_Zg_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_Zg_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_Zg_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Zg_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Zg_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Zg_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Zg_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_Zg_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
	  }
          if(passSystCuts[METUP])  histo_Zg_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_Zg_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 3){
	  if(passAllCuts[SIGSEL]) {
	     histo_WZ->Fill(MVAVar,totalWeight);
	     histo_WZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_WZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_WZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_WZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_WZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_WZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else
	     for(int npdf=0; npdf<102; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_WZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_WZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_WZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_WZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
	  }
          if(passSystCuts[METUP])  histo_WZ_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_WZ_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
	}
        else if(theCategory == 4){
	  if(passAllCuts[SIGSEL]) {
	     histo_ZZ->Fill(MVAVar,totalWeight);
	     histo_ZZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_ZZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_ZZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_ZZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_ZZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_ZZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
          }
          if(passSystCuts[METUP])  histo_ZZ_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_ZZ_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 5){
	  if(passAllCuts[SIGSEL]) {
	     histo_VVV->Fill(MVAVar,totalWeight);
	     histo_VVV_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_VVV_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_VVV_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_VVV_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_VVV_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_VVV_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
          }
          if(passSystCuts[METUP])  histo_VVV_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_VVV_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
	else {
	  printf("CATEGORY PROBLEM!\n"); return;
	}
      } // making data cards
    }
    printf("eff_cuts: %f\n",sumEventsProcess[ifile]);

  } // end of chain

  double sumEvents = 0;
  for(int np=1; np<histBins; np++) sumEvents += histo[4][0][np]->GetSumOfWeights();
  printf("yields: %f |",histo[4][0][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[4][0][np]->GetSumOfWeights());
  printf(" = %.3f\n",sumEvents);
  double sumEventsType[5] = {0,0,0,0,0};
  double sumEventsTypeE[5] = {0,0,0,0,0};
  printf("                  all                 mmm                 eem                 mme                 eee\n");
  printf("-----------------------------------------------------------------------------------------------------------\n");
  for(int np=0; np<histBins; np++) {
    double sumEvents = 0.; double sumEventsE = 0.;
    for(int i=1; i<=4; i++) {sumEvents = sumEvents + histo[4][15][np]->GetBinContent(i); sumEventsE = sumEventsE + histo[4][15][np]->GetBinError(i)*histo[4][15][np]->GetBinError(i);}
    if(np!=0){
    sumEventsType[0] = sumEventsType[0] + sumEvents;                          sumEventsTypeE[0] = sumEventsTypeE[0] + sumEventsE;
    sumEventsType[1] = sumEventsType[1] + histo[4][15][np]->GetBinContent(1); sumEventsTypeE[1] = sumEventsTypeE[1] + histo[4][15][np]->GetBinError(1) * histo[4][15][np]->GetBinError(1);
    sumEventsType[2] = sumEventsType[2] + histo[4][15][np]->GetBinContent(2); sumEventsTypeE[2] = sumEventsTypeE[2] + histo[4][15][np]->GetBinError(2) * histo[4][15][np]->GetBinError(2);
    sumEventsType[3] = sumEventsType[3] + histo[4][15][np]->GetBinContent(3); sumEventsTypeE[3] = sumEventsTypeE[3] + histo[4][15][np]->GetBinError(3) * histo[4][15][np]->GetBinError(3);
    sumEventsType[4] = sumEventsType[4] + histo[4][15][np]->GetBinContent(4); sumEventsTypeE[4] = sumEventsTypeE[4] + histo[4][15][np]->GetBinError(4) * histo[4][15][np]->GetBinError(4);
    }
    sumEventsE = sqrt(sumEventsE);
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    processName[np].Data(),sumEvents,sumEventsE,
    histo[4][15][np]->GetBinContent(1),histo[4][15][np]->GetBinError(1),histo[4][15][np]->GetBinContent(2),histo[4][15][np]->GetBinError(2),
    histo[4][15][np]->GetBinContent(3),histo[4][15][np]->GetBinError(3),histo[4][15][np]->GetBinContent(4),histo[4][15][np]->GetBinError(4));
    if(np==0)
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }
    printf("-----------------------------------------------------------------------------------------------------------\n");
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    "   all",sumEventsType[0],sqrt(sumEventsTypeE[0]),
    sumEventsType[1],sqrt(sumEventsTypeE[1]),sumEventsType[2],sqrt(sumEventsTypeE[2]),
    sumEventsType[3],sqrt(sumEventsTypeE[3]),sumEventsType[4],sqrt(sumEventsTypeE[4]));
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    for(int j=0; j<allStates; j++){
      char output[200];
      sprintf(output,"histowz_nice_%d_%d.root",j,thePlot);	  
      TFile* outFilePlotsNote = new TFile(output,"recreate");
      outFilePlotsNote->cd();
      for(int np=0; np<histBins; np++) histo[j][thePlot][np]->Write();
      outFilePlotsNote->Close();
    }
  }
  printf("QCD Corr: WZ(%f:%f/%f/%f/%f/%f/%f) ZZ(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) ZH(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_WZ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_Zg->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[5]->GetSumOfWeights());

  for(int i=1; i<=histo_Zg->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_Zg_CMS_MVAZHStatBoundingUp        ->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorUp  *histo_Zg->GetBinError(i),0.000001));
    histo_Zg_CMS_MVAZHStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorDown*histo_Zg->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp      ->SetBinContent(i,TMath::Max(histo_VVV    ->GetBinContent(i)+factorUp  *histo_VVV    ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_VVV    ->GetBinContent(i)+factorDown*histo_VVV    ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_WZ     ->GetBinContent(i)+factorUp  *histo_WZ     ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_WZ     ->GetBinContent(i)+factorDown*histo_WZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_ZZ     ->GetBinContent(i)+factorUp  *histo_ZZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_ZZ     ->GetBinContent(i)+factorDown*histo_ZZ     ->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_FakeM  ->GetBinContent(i)+factorUp  *histo_FakeM  ->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingDown->SetBinContent(i,TMath::Max(histo_FakeM  ->GetBinContent(i)+factorDown*histo_FakeM  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorUp  *histo_FakeE  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingDown->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorDown*histo_FakeE  ->GetBinError(i),0.000001));

    histo_Zg_CMS_MVAZHStatBoundingBinUp[i-1]        ->Add(histo_Zg     ); histo_Zg_CMS_MVAZHStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorUp  *histo_Zg->GetBinError(i),0.000001));
    histo_Zg_CMS_MVAZHStatBoundingBinDown[i-1]      ->Add(histo_Zg     ); histo_Zg_CMS_MVAZHStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorDown*histo_Zg->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	    ->Add(histo_VVV    ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]     ->SetBinContent(i,TMath::Max(histo_VVV	->GetBinContent(i)+factorUp  *histo_VVV    ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]    ->Add(histo_VVV    ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]   ->SetBinContent(i,TMath::Max(histo_VVV	->GetBinContent(i)+factorDown*histo_VVV    ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]	    ->Add(histo_WZ     ); histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]       ->SetBinContent(i,TMath::Max(histo_WZ	->GetBinContent(i)+factorUp  *histo_WZ     ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]	    ->Add(histo_WZ     ); histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]     ->SetBinContent(i,TMath::Max(histo_WZ	->GetBinContent(i)+factorDown*histo_WZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]	    ->Add(histo_ZZ     ); histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]       ->SetBinContent(i,TMath::Max(histo_ZZ	->GetBinContent(i)+factorUp  *histo_ZZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]	    ->Add(histo_ZZ     ); histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]     ->SetBinContent(i,TMath::Max(histo_ZZ	->GetBinContent(i)+factorDown*histo_ZZ     ->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[i-1]  ->Add(histo_FakeM  ); histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[i-1]    ->SetBinContent(i,TMath::Max(histo_FakeM  ->GetBinContent(i)+factorUp  *histo_FakeM  ->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[i-1]->Add(histo_FakeM  ); histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[i-1]  ->SetBinContent(i,TMath::Max(histo_FakeM  ->GetBinContent(i)+factorDown*histo_FakeM  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[i-1]  ->Add(histo_FakeE  ); histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[i-1]    ->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorUp  *histo_FakeE  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[i-1]->Add(histo_FakeE  ); histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[i-1]  ->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorDown*histo_FakeE  ->GetBinError(i),0.000001));
  }
  char outputLimits[200];
  sprintf(outputLimits,"wz3l%2s_input_%4s.root",finalStateName,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data   ->Write();
  histo_Zg     ->Write();
  histo_VVV    ->Write();
  histo_WZ     ->Write();
  histo_ZZ     ->Write();
  histo_FakeM  ->Write();
  histo_FakeE  ->Write();
  cout << histo_Data   ->GetSumOfWeights() << " ";
  cout << histo_Zg     ->GetSumOfWeights() << " ";
  cout << histo_VVV    ->GetSumOfWeights() << " ";
  cout << histo_WZ     ->GetSumOfWeights() << " ";
  cout << histo_ZZ     ->GetSumOfWeights() << " ";
  cout << histo_FakeM  ->GetSumOfWeights() << " ";
  cout << histo_FakeE  ->GetSumOfWeights() << " ";
  cout << endl;
  printf("uncertainties Stat\n");
  histo_Zg_CMS_MVAZHStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg	->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAZHStatBoundingUp   ->GetBinContent(i)/histo_Zg   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVAZHStatBoundingDown      ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg	->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAZHStatBoundingDown ->GetBinContent(i)/histo_Zg	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingDown    ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingUp	     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingDown      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingUp	     ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingDown      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_FakeM_CMS_MVAFakeMStatBoundingUp  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_FakeM->GetBinContent(i)>0)printf("%5.1f ",histo_FakeM_CMS_MVAFakeMStatBoundingUp	     ->GetBinContent(i)/histo_FakeM   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_FakeM_CMS_MVAFakeMStatBoundingDown->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_FakeM->GetBinContent(i)>0)printf("%5.1f ",histo_FakeM_CMS_MVAFakeMStatBoundingDown      ->GetBinContent(i)/histo_FakeM   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_FakeE->GetBinContent(i)>0)printf("%5.1f ",histo_FakeE_CMS_MVAFakeEStatBoundingUp	     ->GetBinContent(i)/histo_FakeE   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_FakeE_CMS_MVAFakeEStatBoundingDown->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_FakeE->GetBinContent(i)>0)printf("%5.1f ",histo_FakeE_CMS_MVAFakeEStatBoundingDown      ->GetBinContent(i)/histo_FakeE   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffM\n");
  histo_Zg_CMS_MVALepEffMBoundingUp       ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVALepEffMBoundingDown     ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffE\n");
  histo_Zg_CMS_MVALepEffEBoundingUp       ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVALepEffEBoundingDown     ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties MET\n");
  histo_Zg_CMS_MVAMETBoundingUp           ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_Zg      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVAMETBoundingDown         ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAMETBoundingDown   ->GetBinContent(i)/histo_Zg   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  outFileLimits->Close();

  double lumiE = 1.046;
  double systLepResE[4] = {1.01,1.01,1.01,1.01};
  double systLepResM[4] = {1.01,1.01,1.01,1.01};
  double syst_btag = 1.02;
  for(int nb=1; nb<=nBinMVA; nb++){
     // QCD study
    double systQCDScale[4] = {TMath::Abs(histo_Zg_CMS_QCDScaleBounding[0] ->GetBinContent(nb)-histo_Zg ->GetBinContent(nb)),
                              TMath::Abs(histo_VVV_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_VVV->GetBinContent(nb)),
                              TMath::Abs(histo_WZ_CMS_QCDScaleBounding[0] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb)),
                              TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[0] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb))};
    for(int nqcd=1; nqcd<6; nqcd++) {
      if(TMath::Abs(histo_Zg_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb) -histo_Zg ->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_Zg_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_Zg ->GetBinContent(nb));
      if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb));
      if(TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb) -histo_WZ ->GetBinContent(nb)) > systQCDScale[2]) systQCDScale[2] = TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb));
      if(TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb) -histo_ZZ ->GetBinContent(nb)) > systQCDScale[3]) systQCDScale[3] = TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb));
    }                 
    if(histo_Zg ->GetBinContent(nb) != 0) systQCDScale[0] = 1 + systQCDScale[0]/histo_Zg ->GetBinContent(nb); else systQCDScale[0] = 1;
    if(histo_VVV->GetBinContent(nb) != 0) systQCDScale[1] = 1 + systQCDScale[1]/histo_VVV->GetBinContent(nb); else systQCDScale[1] = 1;
    if(histo_WZ ->GetBinContent(nb) != 0) systQCDScale[2] = 1 + systQCDScale[2]/histo_WZ ->GetBinContent(nb); else systQCDScale[2] = 1;
    if(histo_ZZ ->GetBinContent(nb) != 0) systQCDScale[3] = 1 + systQCDScale[3]/histo_ZZ ->GetBinContent(nb); else systQCDScale[3] = 1;
    printf("QCDScale(%d): %f %f %f %f\n",nb,systQCDScale[0],systQCDScale[1],systQCDScale[2],systQCDScale[3]);
    
    // PDF study
    double systPDF[4];
    histo_Diff->Reset();
    if(histo_Zg ->GetBinContent(nb) > 0) for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_Zg_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_Zg ->GetBinContent(nb))/histo_Zg ->GetBinContent(nb));
    systPDF[0] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    if(histo_VVV->GetBinContent(nb) > 0) for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_VVV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))/histo_VVV->GetBinContent(nb));
    systPDF[1] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    if(histo_WZ ->GetBinContent(nb) > 0) for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_WZ_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb))/histo_WZ ->GetBinContent(nb));
    systPDF[2] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    if(histo_ZZ ->GetBinContent(nb) > 0) for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_ZZ_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb))/histo_ZZ ->GetBinContent(nb));
    systPDF[3] = 1.0+histo_Diff->GetRMS();
    printf("PDF(%d): %f %f %f %f\n",nb,systPDF[0],systPDF[1],systPDF[2],systPDF[3]);

    double systLepEffM[4] = {1.0,1.0,1.0,1.0};
    if     (histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_Zg_CMS_MVALepEffMBoundingUp   ->GetBinContent(nb) > 0) systLepEffM[0] = histo_Zg_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_Zg_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[0] = histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_Zg_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffMBoundingDown     ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_WZ_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_WZ_CMS_MVALepEffMBoundingDown	   ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_ZZ_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_ZZ_CMS_MVALepEffMBoundingDown	   ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);

    double systLepEffE[4] = {1.0,1.0,1.0,1.0};
    if     (histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_Zg_CMS_MVALepEffEBoundingUp   ->GetBinContent(nb) > 0) systLepEffE[0] = histo_Zg_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_Zg_CMS_MVALepEffEBoundingDown ->GetBinContent(nb) > 0) systLepEffE[0] = histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_Zg_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[1] = histo_VVV_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffEBoundingDown     ->GetBinContent(nb) > 0) systLepEffE[1] = histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_WZ_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[2] = histo_WZ_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_WZ_CMS_MVALepEffEBoundingDown	   ->GetBinContent(nb) > 0) systLepEffE[2] = histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_WZ_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_ZZ_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[3] = histo_ZZ_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_ZZ_CMS_MVALepEffEBoundingDown	   ->GetBinContent(nb) > 0) systLepEffE[3] = histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingDown->GetBinContent(nb);

    double systMet[4] = {1.0,1.0,1.0,1.0};
    if     (histo_Zg->GetBinContent(nb) > 0 && histo_Zg_CMS_MVAMETBoundingUp   ->GetBinContent(nb) > 0) systMet[0] = histo_Zg_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    else if(histo_Zg->GetBinContent(nb) > 0 && histo_Zg_CMS_MVAMETBoundingDown ->GetBinContent(nb) > 0) systMet[0] = histo_Zg->GetBinContent(nb)/histo_Zg_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[1] = histo_VVV_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAMETBoundingDown     ->GetBinContent(nb) > 0) systMet[1] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_WZ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[2] = histo_WZ_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVAMETBoundingDown      ->GetBinContent(nb) > 0) systMet[2] = histo_WZ->GetBinContent(nb)/histo_WZ_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[3] = histo_ZZ_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVAMETBoundingDown      ->GetBinContent(nb) > 0) systMet[3] = histo_ZZ->GetBinContent(nb)/histo_ZZ_CMS_MVAMETBoundingDown->GetBinContent(nb);

    char outputLimitsShape[200];
    sprintf(outputLimitsShape,"histo_limits_wz3l%2s_%4s_bin%d.txt",finalStateName,ECMsb.Data(),nb-1);
    ofstream newcardShape;
    newcardShape.open(outputLimitsShape);
    newcardShape << Form("imax 1 number of channels\n");
    newcardShape << Form("jmax * number of background\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");
    newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
    newcardShape << Form("bin wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
    newcardShape << Form("process Zg VVV WZ ZZ FakeM FakeE\n");
    newcardShape << Form("process 1 2 0 3 4 5\n");
    newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",histo_Zg->GetBinContent(nb),histo_VVV->GetBinContent(nb),histo_WZ->GetBinContent(nb),histo_ZZ->GetBinContent(nb),histo_FakeM->GetBinContent(nb),histo_FakeE->GetBinContent(nb));
    newcardShape << Form("lumi_%4s                               lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE);			
    newcardShape << Form("%s                                     lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",effMName,systLepEffM[0],systLepEffM[1],systLepEffM[2],systLepEffM[3]);
    newcardShape << Form("%s                                     lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",effEName,systLepEffE[0],systLepEffE[1],systLepEffE[2],systLepEffE[3]);
    newcardShape << Form("%s                                     lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",momMName,systLepResM[0],systLepResM[1],systLepResM[2],systLepResM[3]);
    newcardShape << Form("%s                                     lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",momEName,systLepResE[0],systLepResE[1],systLepResE[2],systLepResE[3]);
    newcardShape << Form("CMS_scale_met                          lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",systMet[0],systMet[1],systMet[2],systMet[3]);
    newcardShape << Form("CMS_eff_b                              lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",syst_btag,syst_btag,syst_btag,syst_btag);
    newcardShape << Form("pdf_qqbar                              lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",systPDF[0],systPDF[1],systPDF[2],systPDF[3]);
    newcardShape << Form("QCDscale_VVV		                 lnN    -     %7.5f   -     -	  -    -  \n",systQCDScale[1]); 	   
    newcardShape << Form("QCDscale_VV		                 lnN  %7.5f     -   %7.5f %7.5f   -    -  \n",systQCDScale[0],systQCDScale[2],systQCDScale[3]); 	   
    newcardShape << Form("CMS_wz3l_FakeMSyst_%4s                 lnN -       -     -	 -   %7.5f  -  \n",ECMsb.Data(),1.36);  	
    newcardShape << Form("CMS_wz3l_FakeESyst_%4s                 lnN -       -     -	 -     -  %7.5f\n",ECMsb.Data(),1.36);  	
    if(histo_Zg->GetBinContent(nb)        > 0) newcardShape << Form("CMS_wz3l%s_MVAZgStatBounding_%s_Bin%d	  lnN    %7.5f   -    -    -	-    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_Zg    ->GetBinError(nb)/histo_Zg    ->GetBinContent(nb));
    if(histo_VVV->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%d       lnN      -  %7.5f   -    -	-    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_VVV   ->GetBinError(nb)/histo_VVV   ->GetBinContent(nb));
    if(histo_WZ->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAWZStatBounding_%s_Bin%d	  lnN      -     -  %7.5f  -	-    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_WZ    ->GetBinError(nb)/histo_WZ    ->GetBinContent(nb));
    if(histo_ZZ->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAZZStatBounding_%s_Bin%d	  lnN      -     -    -  %7.5f  -    -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_ZZ    ->GetBinError(nb)/histo_ZZ    ->GetBinContent(nb));
    if(histo_FakeM->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%d	  lnN      -     -    -    -  %7.5f  -  \n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_FakeM ->GetBinError(nb)/histo_FakeM ->GetBinContent(nb));
    if(histo_FakeE->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_wz3l%s_MVAFameEStatBounding_%s_Bin%d	  lnN      -     -    -    -	-  %7.5f\n",finalStateName,ECMsb.Data(),nb-1,1.0+histo_FakeE ->GetBinError(nb)/histo_FakeE ->GetBinContent(nb));
    newcardShape.close();
  }
}
