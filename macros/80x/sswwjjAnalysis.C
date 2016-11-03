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

#include "MitAnalysisRunII/macros/80x/factors.h"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"

enum selType                     {SIGSEL=0, TOPSEL,   WZSEL,   DILSEL,   WW2J,   ZSEL, nSelTypes};
TString selTypeName[nSelTypes]= {"SIGSEL", "TOPSEL", "WZSEL", "DILSEL", "WW2J", "ZSEL"};

enum systType                     {JESUP=0, JESDOWN,  METUP,  METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN"};

bool isMINIAOD = false;
int whichSkim = 6;
double mcPrescale = 1.0;
bool usePureMC = false;
int period = 1;
const TString typeLepSel = "default";
const bool usePUPPI = false;

void sswwjjAnalysis(
 ){

  TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_80x/met_";
  if(isMINIAOD) filesPathDA = "/scratch5/dhsu/ntuples_goodrun_80x/met_";
  TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_80x/met_";
  Double_t lumi = 12.9;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  TString triggerSuffix = "*";
  if(isMINIAOD) triggerSuffix = "";
  if      (period==1){
  puPath = "MitAnalysisRunII/data/80x/puWeights_80x_13p0ifb_62_64_66.root";

  //data samples
  if(isMINIAOD) {
    infilenamev.push_back(Form("%sdata_Run2016B_skim.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016C_skim.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016D_skim.root",filesPathDA.Data())); infilecatv.push_back(0);
  } else {
    infilenamev.push_back(Form("%sdata_Run2016B.root",filesPathDA.Data()));   infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016C.root",filesPathDA.Data()));   infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016D.root",filesPathDA.Data()));   infilecatv.push_back(0);
  }

  //MC samples
  //signal: EWK + QCD
  infilenamev.push_back(Form("%sWpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sWpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                    infilecatv.push_back(1);

  //QCD to be subtracted from signal
  infilenamev.push_back(Form("%sWpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                    infilecatv.push_back(-1);

  //QCD to be added to background
  infilenamev.push_back(Form("%sWpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                    infilecatv.push_back(2);

  //VV
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                        infilecatv.push_back(3); 
  infilenamev.push_back(Form("%sWLLJJToLNu_M-60_EWK_13TeV-madgraph-pythia8+RunIISpring16DR80-premix_withHLT_80X_mcRun2_asymptotic_v14-v1+AODSIM.root",filesPathMC.Data()));                        infilecatv.push_back(3);

  ////infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                                   infilecatv.push_back(3);
  ////infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(3);
  ////infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(3);
  ////infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(3);
  ////infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(3);
  ////infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));  		   infilecatv.push_back(3);
  ////infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data())); 		   infilecatv.push_back(3);
  ////infilenamev.push_back(Form("%sZZJJTo4L_EWK_13TeV-madgraph-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(3);

  //VVV
  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  ////infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(4);
  ////infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(4);
  ////infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  ////infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  ////infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(4);
  ////infilenamev.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));  		   infilecatv.push_back(4);

  //Wrong sign
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                                            infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                                        infilecatv.push_back(5);

  ////infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                       infilecatv.push_back(5);

  ////infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                                    infilecatv.push_back(5);
  ////infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                       infilecatv.push_back(5);
  ////infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(5);
  ////infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                   infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sTT_TuneCUETP8M1_13TeV-powheg-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext3-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(5);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(5);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));        infilecatv.push_back(5);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1+AODSIM.root",filesPathMC.Data()));       infilecatv.push_back(5);

  //Wgamma
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                   infilecatv.push_back(6);
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                   infilecatv.push_back(6);

  //DPS
  infilenamev.push_back(Form("%sWW_DoubleScattering_13TeV-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                                      infilecatv.push_back(7);

  //Non-prompt leptons
  //infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                infilecatv.push_back(8)
  }
  else {assert(0); return;}

  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}
  
  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);

  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
  delete fPUFile;

  TFile *fTrackElectronReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_13p0ifb.root"));
  TH2D *fhDeltrksf= (TH2D*)(fTrackElectronReco_SF->Get("scalefactors_Reco_Electron")); assert(fhDeltrksf); fhDeltrksf->SetDirectory(0);
  delete fTrackElectronReco_SF;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_13p0ifb.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("scalefactors_Medium_Electron"));
  if(typeLepSel == "medium_mva") fhDElMediumSF = (TH2D*)(fElSF->Get("scalefactors_MediumMVA_Electron"));
  TH2D *fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_Tight_Electron"));
  if(typeLepSel == "default_mva") fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_TightMVA_Electron"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF->SetDirectory(0);
  delete fElSF;

  TFile *fTrackMuonReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/trackMuReco_SF.root"));
  TH1D *fhDmutrksfptg10 = (TH1D*)(fTrackMuonReco_SF->Get("mutrksfptg10")); assert(fhDmutrksfptg10); fhDmutrksfptg10->SetDirectory(0);
  //TH1D *fhDmutrksfptl10 = (TH1D*)(fTrackMuonReco_SF->Get("mutrksfptl10")); assert(fhDmutrksfptl10); fhDmutrksfptl10->SetDirectory(0);
  delete fTrackMuonReco_SF;

  //TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x.root"));
  //TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("scalefactors_Tight_Muon")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonID_Z_RunBCD_prompt80X_7p65.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  delete fMuSF;

  TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonIso_Z_RunBCD_prompt80X_7p65.root"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
  delete fMuIsoSF;

  TFile *fWWPtRatio = TFile::Open(Form("MitAnalysisRunII/data/74x/MyRatioWWpTHistogramAll.root"));
  TH1D *fhDWWPtRatio           = (TH1D*)(fWWPtRatio->Get("wwpt"));
  assert(fhDWWPtRatio	       );
  fhDWWPtRatio  	->SetDirectory(0);
  delete fWWPtRatio;

  TString ECMsb  = "13TeV2016";
  const int nBinMVA = 6; Float_t xbins[nBinMVA+1] = {500, 700, 900, 1000, 1500, 2000, 3000};
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();

  TH1D *histo_Data  = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_EWK   = (TH1D*) histoMVA->Clone("histo_EWK"); 
  TH1D *histo_QCD   = (TH1D*) histoMVA->Clone("histo_QCD");    
  TH1D *histo_VV    = (TH1D*) histoMVA->Clone("histo_VV");
  TH1D *histo_VVV   = (TH1D*) histoMVA->Clone("histo_VVV");
  TH1D *histo_WS    = (TH1D*) histoMVA->Clone("histo_WS");
  TH1D *histo_WG    = (TH1D*) histoMVA->Clone("histo_WG");
  TH1D *histo_DPS   = (TH1D*) histoMVA->Clone("histo_DPS");
  TH1D *histo_FakeM = (TH1D*) histoMVA->Clone("histo_FakeM");  
  TH1D *histo_FakeE = (TH1D*) histoMVA->Clone("histo_FakeE");  

  double totalFakeDataCount[4][4];
  for(int i=0; i<4; i++) for(int j=0; j<4; j++) totalFakeDataCount[i][j] = 0;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 20;
  const int histBins = 10;
  TH1D* histo[7][allPlots][histBins];
  TString processName[histBins] = {".Data", "EWKWW", "QCDWW", "...VV", "..VVV", "...WS", "...WG", "..DPS", "FakeM", "FakeE"};

  for(int nModel=0; nModel<7; nModel++){
    for(int thePlot=0; thePlot<allPlots; thePlot++){
      if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 2000;}
      else if(thePlot >=  1 && thePlot <=  1) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot = 8;}
      else if(thePlot >=  2 && thePlot <=  2) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot = 6.5;}
      else if(thePlot >=  3 && thePlot <=  4) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400;}
      else if(thePlot >=  5 && thePlot <=  5) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot = 39.5;}
      else if(thePlot >=  6 && thePlot <=  6) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot = 3.5;}
      else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot = 40;}
      else if(thePlot >=  8 && thePlot <=  8) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot = 50;}
      else if(thePlot >=  9 && thePlot <=  9) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 1.0;}
      else if(thePlot >= 10 && thePlot <= 10) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot = 3.5;}
      else if(thePlot >= 11 && thePlot <= 12) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 2000;}
      else if(thePlot >= 13 && thePlot <= 14) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot = 8;}
      else if(thePlot >= 15 && thePlot <= 16) {nBinPlot = 100; xminPlot =-5.0; xmaxPlot = 5.0;}
      else if(thePlot >= 17 && thePlot <= 17) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200;}
      TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
      histos->Sumw2();
      for(int i=0; i<histBins; i++) histo[nModel][thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
      histos->Reset();histos->Clear();
    }
  }

  double bgdDecay[6][nSelTypes][histBins],weiDecay[6][nSelTypes][histBins];
  for(int nModel=0; nModel<6; nModel++) { for(unsigned int i=0; i<nSelTypes; i++) { for(int j=0; j<histBins; j++) {       
    bgdDecay[nModel][i][j] = 0.0; weiDecay[nModel][i][j] = 0.0; 
  }}}

 const int numberCuts = 11;
  TH1D* hDWWLL[7];
  hDWWLL[0] = new TH1D("hDWWLL0", "hDWWLL0", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[1] = new TH1D("hDWWLL1", "hDWWLL1", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[2] = new TH1D("hDWWLL2", "hDWWLL2", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[3] = new TH1D("hDWWLL3", "hDWWLL3", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[4] = new TH1D("hDWWLL4", "hDWWLL4", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[5] = new TH1D("hDWWLL5", "hDWWLL5", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[6] = new TH1D("hDWWLL6", "hDWWLL6", numberCuts, -0.5, numberCuts-0.5);
  TString cutName[numberCuts] = {"#lep, sign, flavor","        lep pT cut","  #jets with pT>30","         btag-veto","          tau veto",
                                 "           mll cut","      loose Z veto","            Z veto","           met cut","           mjj cut",
				 "    deltaetajj cut"};
  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");
    TTree *the_PDF_tree   = (TTree*)the_input_file.FindObjectAny("pdfReweight");
    TTree *the_SelBit_tree= (TTree*)the_input_file.FindObjectAny("SelBit_tree");

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareJets eventJets;
    if(usePUPPI) eventJets.setBranchAddresses(the_input_tree,"puppi");
    else         eventJets.setBranchAddresses(the_input_tree);

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

    unsigned int selBit_= 0;
    the_SelBit_tree->SetBranchAddress("selBit", &selBit_);
    hDWWLL[0]->Scale(0.0);hDWWLL[1]->Scale(0.0);hDWWLL[2]->Scale(0.0);hDWWLL[3]->Scale(0.0);hDWWLL[4]->Scale(0.0);hDWWLL[5]->Scale(0.0);hDWWLL[6]->Scale(0.0);
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_SelBit_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      if((selBit_ & 0x1<<whichSkim) == 0) continue;

      the_input_tree->GetEntry(i);
 
      int sumEvol[6] = {-1, -1, -1, -1, -1, -1};
      Bool_t passFilterSig[numberCuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      Bool_t passFilterCR1[numberCuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      Bool_t passFilterCR2[numberCuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      Bool_t passFilterCR3[numberCuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      Bool_t passPresel = kFALSE;
      Bool_t passSignalRegion = kFALSE;
      Bool_t passControlRegionTop  = kFALSE;
      Bool_t passControlRegionWZ   = kFALSE;
      Bool_t passControlRegionDi   = kFALSE;
      Bool_t passControlRegionZ    = kFALSE;
      Bool_t passControlRegionSS2j = kFALSE;

      if(infilecatv[ifile] == 0) {
        for (int nt = 0; nt <(int)numtokens; nt++) {
          if((*eventTrigger.triggerFired)[nt] == 0) continue;
          if((strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
             (strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
             (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data())) == 0) ||
             (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data())) == 0) ||
             (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
             (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v%s",triggerSuffix.Data())) 	          == 0) ||
             (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v%s",triggerSuffix.Data()))	          == 0) ||
             (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v%s",triggerSuffix.Data())) 	          == 0) ||
             (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v%s",triggerSuffix.Data()))	          == 0) ||
             (strcmp(tokens[nt],Form("HLT_IsoMu20_v%s",triggerSuffix.Data())) 				          == 0) ||
             (strcmp(tokens[nt],Form("HLT_IsoTkMu20_v%s",triggerSuffix.Data())) 				  == 0) ||
             (strcmp(tokens[nt],Form("HLT_IsoMu22_v%s",triggerSuffix.Data())) 				          == 0) ||
             (strcmp(tokens[nt],Form("HLT_IsoTkMu22_v%s",triggerSuffix.Data())) 				  == 0) ||
             (strcmp(tokens[nt],Form("HLT_IsoMu24_v%s",triggerSuffix.Data()))				          == 0) ||
             (strcmp(tokens[nt],Form("HLT_IsoTkMu24_v%s",triggerSuffix.Data()))				          == 0) ||
             (strcmp(tokens[nt],Form("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))	  == 0) ||
             (strcmp(tokens[nt],Form("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))	  == 0) ||
             (strcmp(tokens[nt],Form("HLT_Ele25_eta2p1_WPTight_Gsf_v%s",triggerSuffix.Data()))                    == 0) ||
             (strcmp(tokens[nt],Form("HLT_Ele27_eta2p1_WPLoose_Gsf_v%s",triggerSuffix.Data()))                    == 0) ||
             (strcmp(tokens[nt],Form("HLT_Ele27_WPTight_Gsf_v%s",triggerSuffix.Data()))			          == 0) ||
             (strcmp(tokens[nt],Form("HLT_Ele35_WPLoose_Gsf_v%s",triggerSuffix.Data()))			          == 0)
           ) passPresel = kTRUE;
        }
      } else { passPresel = kTRUE;}

      //acess lepton info
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0; vector<int> idLepLoose;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        //if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep]) == 11 && ((int)(*eventLeptons.selBits)[nlep] & BareLeptons::EleTripleCharge)  != BareLeptons::EleTripleCharge) continue;
        //if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep]) == 11 && ((int)(*eventLeptons.selBits)[nlep] & BareLeptons::EleNoMissingHits)  != BareLeptons::EleNoMissingHits) continue;
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                                   {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake )     {idTight.push_back(0); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP)    {idSoft.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepBaseline)== BareLeptons::LepBaseline){idLepLoose.push_back(nlep);}
      }

      if(idLep.size()!=idTight.size()) {assert(1); return;}

      if(idLep.size() == 2 || idLep.size() == 3) passPresel = passPresel && kTRUE;
      else                                       passPresel = passPresel && kFALSE;

      if(passPresel == kFALSE) continue; // two or three leptons in the event

      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 25 && 
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 25) passPresel = passPresel && kTRUE;
      else                                                           passPresel = passPresel && kFALSE;

      if(passPresel == kFALSE) continue; // ptl1/l2 > 25

      //typeSel = 0(m+m+), 1(e+e+), 2(e+m+/m+e+) 3(m-m-), 4(e-e-), 5(e-m-/m-e-)
      int typeSel = -1;
      if(idTight.size() >= 2){
        if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) typeSel = 0;
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) typeSel = 1;
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])!=TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]))         typeSel = 2;
        else {assert(1); return;}
	if((int)(*eventLeptons.pdgId)[idLep[0]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) > 0) typeSel = typeSel + 3;
      }
      else {assert(1); printf("Not possible %d %d %d\n",(int)idLep.size(),(int)idTight.size(),goodIsTight); return;}                                                                                                     ;

      //signQ = 0(opposite sign), +/-2(same sign)
      int signQ = -1;
      signQ = (int)(*eventLeptons.pdgId)[idLep[0]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) + (int)(*eventLeptons.pdgId)[idLep[1]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]);

      //#lep, same sign or opposite sign, flavor combination choice
      if(idTight.size() == 2 && signQ != 0){
	passFilterSig[0] = kTRUE;
	passFilterCR1[0] = kTRUE;
      }
      if(idTight.size() == 2 && signQ == 0 && (typeSel == 2 || typeSel == 5))  passFilterCR3[0] = kTRUE;
      if(idTight.size() == 3 )  					       passFilterCR2[0] = kTRUE; 

      //lepton pT cut
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 25 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 25){
	passFilterSig[1] = kTRUE;
	passFilterCR1[1] = kTRUE;
	passFilterCR3[1] = kTRUE;
	if(idTight.size() == 3 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() > 10) passFilterCR2[1] = kTRUE;
      }

      //access jet info
      vector<int> idJet,idJetUp,idJetDown;
      double bDiscrMax = 0.0;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 15) continue;

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.4) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 20) { 
	   if ((float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];
        }
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()      > 30) {idJet.push_back(nj);}
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*1.03 > 30) idJetUp.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*0.97 > 30) idJetDown.push_back(nj);
      }

      //#Jet with pT > 30 GeV
      if(idJet.size() >= 2){
	passFilterSig[2] = kTRUE;
        passFilterCR1[2] = kTRUE;
        passFilterCR2[2] = kTRUE;
        passFilterCR3[2] = kTRUE;
      }

      int numberGoodTaus = 0;
      for(int ntau=0; ntau<eventTaus.p4->GetEntriesFast(); ntau++) {
        if(((TLorentzVector*)(*eventTaus.p4)[ntau])->Pt() <= 18.0 ||
           TMath::Abs(((TLorentzVector*)(*eventTaus.p4)[ntau])->Eta()) >= 2.3) continue;
        bool isElMu = false;
        for(unsigned nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventTaus.p4)[ntau])) < 0.1) {
            isElMu = true;
            break;
          }
        }
        if(isElMu == false &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::TauDecayModeFinding	 ) == BareTaus::TauDecayModeFinding &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::TauDecayModeFindingNewDMs) == BareTaus::TauDecayModeFindingNewDMs &&
           (double)(*eventTaus.iso)[ntau] < 5.0){
          numberGoodTaus++;
        }
      }

      //b-tag selection
      if(bDiscrMax < 0.700 && idSoft.size() == 0){
	passFilterSig[3] = kTRUE;
        passFilterCR2[3] = kTRUE;
        passFilterCR3[3] = kTRUE;
      } else{
        passFilterCR1[3] = kTRUE;
      }

      //tau veto
      if(numberGoodTaus == 0){
	passFilterSig[4] = kTRUE;
        passFilterCR1[4] = kTRUE;
        passFilterCR2[4] = kTRUE;
        passFilterCR3[4] = kTRUE;
      }

      //Mll cut
      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) ));
      if(dilep.M() > 50){
	passFilterSig[5] = kTRUE;
	passFilterCR1[5] = kTRUE;
	passFilterCR3[5] = kTRUE;
      }
      passFilterCR2[5] = kTRUE;


      double minMassZ = 999.0;
      for(unsigned nl0=0; nl0<idLep.size(); nl0++){
        for(unsigned nl1=0; nl1<idLepLoose.size(); nl1++){

	  if((int)(*eventLeptons.pdgId)[idLep[nl0]] * (int)(*eventLeptons.pdgId)[idLepLoose[nl1]] > 0) continue;
          TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLepLoose[nl1])) ) ));

	  if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLepLoose[nl1]]) &&
	     TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassZ-91.1876)) {
	     minMassZ = dilepAux.M();
	  }
        }
      }

      // Loose Z veto
      if(TMath::Abs(minMassZ-91.1876) > 15){
	passFilterSig[6] = kTRUE;
        passFilterCR1[6] = kTRUE;
        passFilterCR3[6] = kTRUE;     
      }
      passFilterCR2[6] = kTRUE;
      
      //Z veto
      if(((typeSel == 1 || typeSel == 4) && (TMath::Abs(dilep.M()-91.1876) > 15.0)) || !(typeSel == 1 || typeSel == 4)){
	passFilterSig[7] = kTRUE;
	passFilterCR1[7] = kTRUE;
      }
      if(idTight.size() == 3){
	for(int i1=0; i1<3; i1++) {
	  for(int i2=i1+1; i2<3; i2++) {
	    if((int)(*eventLeptons.pdgId)[idLep[i1]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[i1]]) + (int)(*eventLeptons.pdgId)[idLep[i2]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[i2]])==0){
	      TLorentzVector dilepCR2(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[i1])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[i2])) ) ));
	      if((TMath::Abs(dilepCR2.M()-91.1876) < 15.0)) passFilterCR2[7] = kTRUE;
	    }
	  }
	}
      }
      passFilterCR3[7] = kTRUE;

      //acess met info
      TLorentzVector theMET;
      if(usePUPPI){
        theMET.SetPx((double)eventMet.metPuppi->Px());
        theMET.SetPy((double)eventMet.metPuppi->Py());
      } else {
        theMET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px());
        theMET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py());
      }

      //met cut
      if(theMET.Pt() > 40){
        passFilterSig[8] = kTRUE;
        passFilterCR1[8] = kTRUE;
        passFilterCR2[8] = kTRUE;
        passFilterCR3[8] = kTRUE;
      }

      TLorentzVector dijet; double deltaEtaJJ = 0.0;
      if(idJet.size() >= 2){
	//Mjj cut
	dijet = ( *(TLorentzVector*)(*eventJets.p4)[idJet[0]] ) + ( *(TLorentzVector*)(*eventJets.p4)[idJet[1]] );
	if(dijet.M() > 500){
          passFilterSig[9] = kTRUE;
          passFilterCR1[9] = kTRUE;
          passFilterCR2[9] = kTRUE;
	}
	else if(dijet.M() > 100){
          passFilterCR3[9] = kTRUE;
	}

	//deltaEta jj cut
	deltaEtaJJ = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta());
	if(deltaEtaJJ > 2.5){
          passFilterSig[10] = kTRUE;
          passFilterCR1[10] = kTRUE;
          passFilterCR2[10] = kTRUE;
	}
	passFilterCR3[10] = kTRUE;
      }

      //                      #lep, sign, flavor  lep pT cut          #jets with pT>30    btag-veto           tau veto            mll cut             loose Z veto        Z veto              met cut             mjj cut             deltaetajj cut
      passSignalRegion      = passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] && passFilterSig[9] && passFilterSig[10];
      passControlRegionTop  = passFilterCR1[0] && passFilterCR1[1] && passFilterCR1[2] && passFilterCR1[3] && passFilterCR1[4] && passFilterCR1[5] && passFilterCR1[6] && passFilterCR1[7] && passFilterCR1[8] && passFilterCR1[9] && passFilterCR1[10];
      passControlRegionWZ   = passFilterCR2[0] && passFilterCR2[1] && passFilterCR2[2] && passFilterCR2[3] && passFilterCR2[4] && passFilterCR2[5] && passFilterCR2[6] && passFilterCR2[7] && passFilterCR2[8] && passFilterCR2[9] && passFilterCR2[10];
      passControlRegionDi   = passFilterCR3[0] && passFilterCR3[1] && passFilterCR3[2] && passFilterCR3[3] && passFilterCR3[4] && passFilterCR3[5] && passFilterCR3[6] && passFilterCR3[7] && passFilterCR3[8] && passFilterCR3[9] && passFilterCR3[10];
      passControlRegionSS2j = passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8];
      passControlRegionZ    = passFilterSig[0] && passFilterSig[1] &&                     passFilterSig[3] && passFilterSig[4] &&                     passFilterSig[6];

      bool passNMinusOne[7] = {
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] &&                     passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] && passFilterSig[9] && passFilterSig[10], // btag veto
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] &&                     passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] && passFilterSig[9] && passFilterSig[10], // tau veto
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] &&                     passFilterSig[6] &&                     passFilterSig[8] && passFilterSig[9] && passFilterSig[10], // mll cut && Z veto
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] &&                     passFilterSig[7] && passFilterSig[8] && passFilterSig[9] && passFilterSig[10], // Loose Z veto
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] &&                     passFilterSig[9] && passFilterSig[10], // met cut
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] &&                     passFilterSig[10], // mjj cut
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] && passFilterSig[9]                       // detajj cut
                               };

      bool totalSel = kTRUE;
      for(int isel=0; isel<numberCuts; isel++) {
        totalSel = totalSel && passFilterSig[isel];
	if(totalSel == kTRUE) sumEvol[typeSel]++;
      }

      bool passAllCuts[nSelTypes] = {passSignalRegion,passControlRegionTop,passControlRegionWZ,passControlRegionDi,passControlRegionSS2j,passControlRegionZ};                 

      // begin event weighting
      vector<bool> isGenDupl;
      vector<int>wBoson;
      vector<int>zBoson;
      int numberQuarks[2] = {0,0};
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 4 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[0]++;
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 5 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[1]++;
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 24) {
	  wBoson.push_back(ngen0);
	}
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23) {
	  zBoson.push_back(ngen0);
	}
        isGenDupl.push_back(0);
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) isGenDupl[ngen0] = 1;
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) continue;
        for(int ngen1=ngen0+1; ngen1<eventMonteCarlo.p4->GetEntriesFast(); ngen1++) {
	  if((int)(*eventMonteCarlo.pdgId)[ngen0] != (int)(*eventMonteCarlo.pdgId)[ngen1]) continue;
          if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen1])) < 0.02) {
	    isGenDupl[ngen0] = 1;
	    break;
	  }
        }
      }

      int genLep = 0;
      int numberGoodGenLep[2] = {0,0};
      for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	if(isGenDupl[ngen] == 1) continue;
        genLep++;
	numberGoodGenLep[0]++;
	if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt() <= 10 ||
	   TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()) >= 2.5) continue;
	numberGoodGenLep[1]++;
      }

      vector<int> isGenLep; unsigned int goodIsGenRSLep = 0; unsigned int goodIsGenWSLep = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        bool isGenRSLepton = false;
        bool isGenWSLepton = false;
        for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	  if(isGenDupl[ngen] == 1) continue;
          if((int)(*eventLeptons.pdgId)[idLep[nl]] == (int)(*eventMonteCarlo.pdgId)[ngen] &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.3) {
	    isGenRSLepton = true;
	  }
          else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.3) {
	    isGenWSLepton = true;
	  }
	}
	if     (isGenRSLepton == true) {isGenLep.push_back(1); goodIsGenRSLep++;}
	else if(isGenWSLepton == true) {isGenLep.push_back(2); goodIsGenWSLep++;}
	else                           {isGenLep.push_back(0);}
      }

      // trigger efficiency
      double trigEff = 1.0;
      //if(infilecatv[ifile] != 0) {
      //  trigEff = trigLookup.GetExpectedTriggerEfficiency(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),
      //  						  ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),
      //  						 TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]));
      //}
      // luminosity
      double theLumi  = 1.0; if(infilecatv[ifile] != 0) theLumi  = lumi;
      // pile-up
      double puWeight     = 1.0; if(infilecatv[ifile] != 0) puWeight     = nPUScaleFactor(fhDPU    , (double)eventMonteCarlo.puTrueInt);
      double puWeightUp   = 1.0; if(infilecatv[ifile] != 0) puWeightUp   = nPUScaleFactor(fhDPUUp  , (double)eventMonteCarlo.puTrueInt);
      double puWeightDown = 1.0; if(infilecatv[ifile] != 0) puWeightDown = nPUScaleFactor(fhDPUDown, (double)eventMonteCarlo.puTrueInt);
      // lepton efficiency
      double effSF = 1.0;
      if(infilecatv[ifile] != 0){
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
	  	typeLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,true,fhDMuIsoSF);
        }
      }

      int theCategory = infilecatv[ifile];

      // wrong sign
      if((theCategory == 3 || theCategory == 4) && goodIsGenWSLep > 0) theCategory  = 5;

      // fake rate
      double fakeSF = 1.0;
      if(usePureMC == false){
	if((infilecatv[ifile] == 0 || infilecatv[ifile] == 6 || (goodIsGenRSLep+goodIsGenWSLep) == isGenLep.size()) && goodIsTight != idTight.size()){
            unsigned int typeFakeLepton[2] = {0,0};
            for(unsigned int nl=0; nl<idLep.size(); nl++){
              if(idTight[nl] == 1) continue;
              fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
              theCategory = 8;
	      if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) typeFakeLepton[0]++;
	      else                                                        typeFakeLepton[1]++;
            }
            if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) fakeSF =  1.0 * fakeSF; // double fake, MC
            else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) fakeSF = -1.0 * fakeSF; // single fake, MC
            else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) fakeSF = -1.0 * fakeSF; // double fake, data
            else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) fakeSF =  1.0 * fakeSF; // single fake, data
            if(typeFakeLepton[0] < typeFakeLepton[1]) theCategory = 9;
	}
	else if(infilecatv[ifile] != 0 && infilecatv[ifile] != 6 && (goodIsGenRSLep+goodIsGenWSLep) != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
	}
	else if(infilecatv[ifile] != 0 && (goodIsGenRSLep+goodIsGenWSLep) == isGenLep.size()){ // MC with all good leptons
          fakeSF = 1.0;
	}
	else if(infilecatv[ifile] == 0 || infilecatv[ifile] == 6){ // data or W+gamma with all good leptons
          fakeSF = 1.0;
	}
	else {
          printf("PROBLEMFAKES: %d %d %d %d %d %d\n",infilecatv[ifile],goodIsGenRSLep,goodIsGenWSLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
          assert(0);
	}
      }

      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale*trigEff;

      if(theCategory == -1) {theCategory = 1; totalWeight = -1.0 * totalWeight;}

      if(totalWeight == 0) continue;
      // end event weighting

      for(int nl=0; nl <=sumEvol[typeSel]; nl++) if(fakeSF == 1) {hDWWLL[typeSel]->Fill((double)nl,totalWeight);hDWWLL[6]->Fill((double)nl,totalWeight);}

      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i]) {
          bgdDecay[typeSel][i][theCategory] += totalWeight;
          weiDecay[typeSel][i][theCategory] += totalWeight*totalWeight;
        }
      }

      for(int thePlot=0; thePlot<allPlots; thePlot++){
        double theVar = 0.0;
        bool makePlot = false;
        if     (thePlot ==  0 && passNMinusOne[5])     {makePlot = true;theVar = TMath::Min(dijet.M(),1999.999);}
        else if(thePlot ==  1 && passNMinusOne[6])     {makePlot = true;theVar = TMath::Min(deltaEtaJJ,7.999);}
        else if(thePlot ==  2 && passSignalRegion)     {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
        else if(thePlot ==  3 && passNMinusOne[4])     {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),399.999);}
        else if(thePlot ==  4 && passNMinusOne[2])     {makePlot = true;theVar = TMath::Min(dilep.M(),399.999);}
        else if(thePlot ==  5 && passSignalRegion)     {makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	else if(thePlot ==  6 && passNMinusOne[1])     {makePlot = true;theVar = TMath::Min((double)numberGoodTaus,3.499);}
	else if(thePlot ==  7 && passNMinusOne[3])     {makePlot = true;theVar = TMath::Min(TMath::Abs(minMassZ-91.1876),39.999);}
	else if(thePlot ==  8 && passSignalRegion)     {makePlot = true;theVar = (double)(numberGoodGenLep[1]+10*numberGoodGenLep[0]);}
	else if(thePlot ==  9 && passNMinusOne[0])     {makePlot = true;theVar = TMath::Min(TMath::Max(bDiscrMax,0.001),0.999);}
	else if(thePlot == 10 && passNMinusOne[0])     {makePlot = true;theVar = TMath::Min((double)idSoft.size(),3.499);}
        else if(thePlot == 11 && passControlRegionTop) {makePlot = true;theVar = TMath::Min(dijet.M(),1999.999);}
        else if(thePlot == 12 && passControlRegionWZ)  {makePlot = true;theVar = TMath::Min(dijet.M(),1999.999);}
        else if(thePlot == 13 && passControlRegionTop) {makePlot = true;theVar = TMath::Min(deltaEtaJJ,7.999);}
        else if(thePlot == 14 && passControlRegionWZ)  {makePlot = true;theVar = TMath::Min(deltaEtaJJ,7.999);}
        else if(thePlot == 15 && passControlRegionSS2j){makePlot = true;theVar = ((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta();}
        else if(thePlot == 16 && passControlRegionSS2j){makePlot = true;theVar = ((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta();}
        else if(thePlot == 17 && passControlRegionZ)   {makePlot = true;theVar = TMath::Min(dilep.M(),199.999);}
        if(makePlot) histo[typeSel][thePlot][theCategory]->Fill(theVar,totalWeight);
        if(makePlot) histo[6][thePlot][theCategory]->Fill(theVar,totalWeight);
      }

      if(passSignalRegion==kTRUE) {
	double MVAVar = 1;//(double)dijet.M();
	//double MVAVar = 0;
        if     (theCategory == 0){
	  histo_Data->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 1){
	     histo_EWK->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 2){
	     histo_QCD->Fill(MVAVar,-totalWeight);
        }
        else if(theCategory == 3){
	     histo_VV->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 4){
	     histo_VVV->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 5){
	     histo_WS->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 6){
	     histo_WG->Fill(MVAVar,totalWeight);
        }
	else if(theCategory == 7){
	     histo_DPS->Fill(MVAVar,totalWeight);
	}
        else if(theCategory == 8){
	     histo_FakeM->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 9){
	     histo_FakeE->Fill(MVAVar,totalWeight);
        }
	else {
	  printf("CATEGORY PROBLEM!\n"); return;
	}
      } // making data cards
    }
    printf("                          mm+          ee+          em+          mm-          ee-          em-          all\n");
    for(int nc=0; nc<numberCuts; nc++){
      printf("(%15s): %10.2f   %10.2f   %10.2f   %10.2f   %10.2f   %10.2f   %10.2f\n",cutName[nc].Data(),hDWWLL[0]->GetBinContent(nc+1),hDWWLL[1]->GetBinContent(nc+1),hDWWLL[2]->GetBinContent(nc+1),
                                                                                                hDWWLL[3]->GetBinContent(nc+1),hDWWLL[4]->GetBinContent(nc+1),hDWWLL[5]->GetBinContent(nc+1),
												hDWWLL[6]->GetBinContent(nc+1));
    }
  } // end of chain

  for(int nModel=0; nModel<7; nModel++){
    for(int thePlot=0; thePlot<allPlots; thePlot++){
      char output[200];
      sprintf(output,"histowwss_%d_%d.root",nModel,thePlot);	  
      TFile* outFilePlotsNote = new TFile(output,"recreate");
      outFilePlotsNote->cd();
      for(int np=0; np<histBins; np++) histo[nModel][thePlot][np]->Write();
      outFilePlotsNote->Close();
    }
  }

  for(int ns=0; ns<nSelTypes; ns++) {
    printf("Selection: %s\n",selTypeName[ns].Data());
    double sumEventsType[7] = {0,0,0,0,0,0,0}; double sumEventsTypeE[7] = {0,0,0,0,0,0,0};
    for(int np=0; np<histBins; np++) {
      double tot = 0; double totE = 0;
      for(int nModel=0; nModel<6; nModel++) {
        tot = tot + bgdDecay[nModel][ns][np]; totE = totE + weiDecay[nModel][ns][np];
        if(np!=0){
          sumEventsType[nModel] = sumEventsType[nModel] + bgdDecay[nModel][ns][np]; sumEventsTypeE[nModel] = sumEventsTypeE[nModel] + weiDecay[nModel][ns][np];
	  sumEventsType[     6] = sumEventsType[     6] + bgdDecay[nModel][ns][np]; sumEventsTypeE[     6] = sumEventsTypeE[     6] + weiDecay[nModel][ns][np];
        }
      }
      printf("(%5s): %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f => %8.2f +/- %6.2f\n",
      processName[np].Data(),bgdDecay[0][ns][np],sqrt(weiDecay[0][ns][np]),bgdDecay[1][ns][np],sqrt(weiDecay[1][ns][np]),
        		     bgdDecay[2][ns][np],sqrt(weiDecay[2][ns][np]),bgdDecay[3][ns][np],sqrt(weiDecay[3][ns][np]),
        		     bgdDecay[4][ns][np],sqrt(weiDecay[4][ns][np]),bgdDecay[5][ns][np],sqrt(weiDecay[5][ns][np]),
			     tot,sqrt(totE));
    }
    printf("(..bkg): %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f => %8.2f +/- %6.2f\n",
    	   sumEventsType[0],sqrt(sumEventsTypeE[0]),sumEventsType[1],sqrt(sumEventsTypeE[1]),
           sumEventsType[2],sqrt(sumEventsTypeE[2]),sumEventsType[3],sqrt(sumEventsTypeE[3]),
           sumEventsType[4],sqrt(sumEventsTypeE[4]),sumEventsType[5],sqrt(sumEventsTypeE[5]),
	   sumEventsType[6],sqrt(sumEventsTypeE[6]));
  }
}
