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

#include "MitAnalysisRunII/macros/76x/factors.h"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"

#include "WWAnalysis/resummation/WWpTreweight.h"

enum selType                     { SIGSEL, nSelTypes};
TString selTypeName[nSelTypes]= { "SIGSEL"};

enum systType                     {JESUP=0, JESDOWN,  METUP,  METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN"};

double mcPrescale = 1.0;
bool usePureMC = false; 
//const bool useDYMVA = false;
//const TString typeLepSel = "default";
//const bool usePUPPI = true;

double topNorm[3]  = {0.85,1.02,0.99};
double topNormE[3] = {0.10,0.08,0.01};

void wwAnalysis(
 unsigned int nJetsType = 0,
 bool useDYMVA = false,
 bool usePUPPI = false,
 TString typeLepSel = "default"
 ){

  Int_t period = 1;
  TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_76x/met_";
  TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_76x/met_";
  Double_t lumi = 2.318;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  if      (period==1){
  puPath = "MitAnalysisRunII/data/76x/puWeights_76x.root";

  infilenamev.push_back(Form("%sdata_AOD_Run2015C_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);

  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                                            infilecatv.push_back(1);

  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(2);

  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));						   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(3);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(3);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(4);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(4);

  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  				   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	 	   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(5);

  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3+AODSIM.root",filesPathMC.Data())); 			   infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 			   infilecatv.push_back(6);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                         infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(6);

  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(7);
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                   infilecatv.push_back(7);

  //infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(9);

  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(11);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(11);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(11);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                            infilecatv.push_back(11);
  //infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 		   infilecatv.push_back(11);
  //infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(11);
  }
  else {assert(0); return;}

  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}
  
  //infilenamev.clear();infilecatv.clear();
  //infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg-herwigpp+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(11);
  //infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_herwigpp+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(11);

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
  
  WWpTreweight theWWpTreweight(
              "WWAnalysis/resummation/central.dat",
              "WWAnalysis/resummation/resum_up.dat",
              "WWAnalysis/resummation/resum_down.dat",
              "WWAnalysis/resummation/scale_up.dat",
              "WWAnalysis/resummation/scale_down.dat",
              "WWAnalysis/resummation/nnlo_central.dat",
              "WWAnalysis/resummation/powheg_2l2nu_nlo.dat",
              "WWAnalysis/resummation/powheg_2l2nu_qup_nlo.dat",
              "WWAnalysis/resummation/powheg_2l2nu_qdown_nlo.dat",
              "WWAnalysis/resummation/powheg_2l2nu_sup_nlo.dat",
              "WWAnalysis/resummation/powheg_2l2nu_sdown_nlo.dat",
              "WWAnalysis/resummation/powheg_2l2nu_nnlo.dat");

  float dymva_= -999.;
  unsigned int nlep_= -1;
  unsigned int njets_= -1;

  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Medium_ele"));
  TH2D *fhDElTightSF  = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Tight_ele"));
  TH2D *fhDElMediumMVASF = (TH2D*)(fElSF->Get("unfactorized_scalefactors_MediumMVA_ele"));
  TH2D *fhDElTightMVASF  = (TH2D*)(fElSF->Get("unfactorized_scalefactors_TightMVA_ele"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  assert(fhDElMediumMVASF);
  assert(fhDElTightMVASF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF ->SetDirectory(0);
  fhDElMediumMVASF->SetDirectory(0);
  fhDElTightMVASF ->SetDirectory(0);
  delete fElSF;

  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Medium_mu"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Iso_mu"));
  assert(fhDMuMediumSF);
  assert(fhDMuIsoSF);
  fhDMuMediumSF->SetDirectory(0);
  fhDMuIsoSF->SetDirectory(0);
  delete fMuSF;

  TFile *fWWPtRatio = TFile::Open(Form("/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src/MitAnalysisRunII/data/74x/MyRatioWWpTHistogramAll.root"));
  TH1D *fhDWWPtRatio           = (TH1D*)(fWWPtRatio->Get("wwpt"));
  TH1D *fhDWWPtRatio_scaleup   = (TH1D*)(fWWPtRatio->Get("wwpt_scaleup"));
  TH1D *fhDWWPtRatio_scaledown = (TH1D*)(fWWPtRatio->Get("wwpt_scaledown"));
  TH1D *fhDWWPtRatio_resumup   = (TH1D*)(fWWPtRatio->Get("wwpt_resumup"));
  TH1D *fhDWWPtRatio_resumdown = (TH1D*)(fWWPtRatio->Get("wwpt_resumdown"));
  assert(fhDWWPtRatio	       );
  assert(fhDWWPtRatio_scaleup  );
  assert(fhDWWPtRatio_scaledown);
  assert(fhDWWPtRatio_resumup  );
  assert(fhDWWPtRatio_resumdown);
  fhDWWPtRatio  	->SetDirectory(0);
  fhDWWPtRatio_scaleup  ->SetDirectory(0);
  fhDWWPtRatio_scaledown->SetDirectory(0);
  fhDWWPtRatio_resumup  ->SetDirectory(0);
  fhDWWPtRatio_resumdown->SetDirectory(0);
  delete fWWPtRatio;

  TString ECMsb  = "13TeV2015";
  const int nBinMVA = 2; Float_t xbins[nBinMVA+1] = {0, 1, 2};
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();

  TH1D *histo_Data   = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_qqWW   = (TH1D*) histoMVA->Clone("histo_qqWW"); 
  TH1D *histo_ggWW   = (TH1D*) histoMVA->Clone("histo_ggWW");
  TH1D *histo_Top    = (TH1D*) histoMVA->Clone("histo_Top");	 
  TH1D *histo_DY     = (TH1D*) histoMVA->Clone("histo_DY");
  TH1D *histo_VV     = (TH1D*) histoMVA->Clone("histo_VV");
  TH1D *histo_VVV    = (TH1D*) histoMVA->Clone("histo_VVV");
  TH1D *histo_WG     = (TH1D*) histoMVA->Clone("histo_WG");
  TH1D *histo_WGS    = (TH1D*) histoMVA->Clone("histo_WGS");
  TH1D *histo_WjetsM = (TH1D*) histoMVA->Clone("histo_WjetsM");	 
  TH1D *histo_WjetsE = (TH1D*) histoMVA->Clone("histo_WjetsE");	 
  TH1D *histo_Higgs  = (TH1D*) histoMVA->Clone("histo_Higgs");	 

  double totalFakeDataCount[4][4];
  for(int i=0; i<4; i++) for(int j=0; j<4; j++) totalFakeDataCount[i][j] = 0;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 10;
  const int histBins = 12;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {".Data", ".qqWW", ".ggWW", "..Top", "...DY", "...VV", "..VVV", "...WG", "..WGS", "WjetsM", "WjetsE", "Higgs"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 100; xminPlot =200.0;xmaxPlot = 600.0;}
    else if(thePlot >=  1 && thePlot <=  1) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >=  2 && thePlot <=  2) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Reset();histos->Clear();
  }

  char finalStateName[2],effMName[10],effEName[10],momMName[10],momEName[10];
  sprintf(effMName,"CMS_eff_m");sprintf(momMName,"CMS_scale_m");
  sprintf(effEName,"CMS_eff_e");sprintf(momEName,"CMS_scale_e");
  sprintf(finalStateName,"lnln");

  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingUp       = new TH1D( Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingDown     = new TH1D( Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingUp       = new TH1D( Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingDown     = new TH1D( Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVATopStatBoundingUp         = new TH1D( Form("histo_Top_CMS_ww%s_MVATopStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Top_CMS_ww%s_MVATopStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVATopStatBoundingDown       = new TH1D( Form("histo_Top_CMS_ww%s_MVATopStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Top_CMS_ww%s_MVATopStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVADYStatBoundingUp           = new TH1D( Form("histo_DY_CMS_ww%s_MVADYStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_DY_CMS_ww%s_MVADYStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_DY_CMS_MVADYStatBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVADYStatBoundingDown         = new TH1D( Form("histo_DY_CMS_ww%s_MVADYStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_DY_CMS_ww%s_MVADYStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_DY_CMS_MVADYStatBoundingDown->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingUp     = new TH1D( Form("histo_Zjets_CMS_ww%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_ww%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingDown   = new TH1D( Form("histo_Zjets_CMS_ww%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_ww%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVAVVStatBoundingUp           = new TH1D( Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVAVVStatBoundingDown         = new TH1D( Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp         = new TH1D( Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown       = new TH1D( Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAWGStatBoundingUp           = new TH1D( Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WG_CMS_MVAWGStatBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAWGStatBoundingDown         = new TH1D( Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WG_CMS_MVAWGStatBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVAWGSStatBoundingUp         = new TH1D( Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WGS_CMS_MVAWGSStatBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVAWGSStatBoundingDown       = new TH1D( Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WGS_CMS_MVAWGSStatBoundingDown->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingUp   = new TH1D( Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingUp  ->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingDown = new TH1D( Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingDown->Sumw2();
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingUp   = new TH1D( Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingUp  ->Sumw2();
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingDown = new TH1D( Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingUp     = new TH1D( Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingDown   = new TH1D( Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVAWW_nlo       = new TH1D( Form("histo_qqWW_CMS_MVAWW_nlo"), Form("histo_qqWW_CMS_MVAWW_nlo"), nBinMVA, xbins); histo_qqWW_CMS_MVAWW_nlo->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWW_qup_nlo   = new TH1D( Form("histo_qqWW_CMS_MVAWW_qup_nlo"), Form("histo_qqWW_CMS_MVAWW_qup_nlo"), nBinMVA, xbins); histo_qqWW_CMS_MVAWW_qup_nlo->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWW_qdown_nlo = new TH1D( Form("histo_qqWW_CMS_MVAWW_qdown_nlo"), Form("histo_qqWW_CMS_MVAWW_qdown_nlo"), nBinMVA, xbins); histo_qqWW_CMS_MVAWW_qdown_nlo->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWW_sup_nlo   = new TH1D( Form("histo_qqWW_CMS_MVAWW_sup_nlo"), Form("histo_qqWW_CMS_MVAWW_sup_nlo"), nBinMVA, xbins); histo_qqWW_CMS_MVAWW_sup_nlo->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWW_sdown_nlo = new TH1D( Form("histo_qqWW_CMS_MVAWW_sdown_nlo"), Form("histo_qqWW_CMS_MVAWW_sdown_nlo"), nBinMVA, xbins); histo_qqWW_CMS_MVAWW_sdown_nlo->Sumw2();

  TH1D* histo_Diff = new TH1D("dummy", "dummy",1000,-1,1); histo_Diff->Sumw2();

  TH1D* histo_qqWW_CMS_QCDScaleBounding[6];
  TH1D* histo_ggWW_CMS_QCDScaleBounding[6];
  TH1D* histo_Top_CMS_QCDScaleBounding[6];
  TH1D* histo_DY_CMS_QCDScaleBounding[6];
  TH1D* histo_VV_CMS_QCDScaleBounding[6];
  TH1D* histo_VVV_CMS_QCDScaleBounding[6];
  TH1D* histo_WG_CMS_QCDScaleBounding[6];
  TH1D* histo_WGS_CMS_QCDScaleBounding[6];
  TH1D* histo_Higgs_CMS_QCDScaleBounding[6];
  for(int nb=0; nb<6; nb++){
    histo_qqWW_CMS_QCDScaleBounding[nb]  = new TH1D(Form("histo_qqWW_QCDScale_f%d",nb),     Form("histo_qqWW_QCDScale_f%d",nb),nBinMVA, xbins);     histo_qqWW_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ggWW_CMS_QCDScaleBounding[nb]  = new TH1D(Form("histo_ggWW_QCDScale_f%d",nb),     Form("histo_ggWW_QCDScale_f%d",nb),nBinMVA, xbins);     histo_ggWW_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_Top_CMS_QCDScaleBounding[nb]   = new TH1D(Form("histo_Top_QCDScale_f%d",nb),     Form("histo_Top_QCDScale_f%d",nb),nBinMVA, xbins);     histo_Top_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_DY_CMS_QCDScaleBounding[nb]    = new TH1D(Form("histo_DY_QCDScale_f%d",nb),     Form("histo_DY_QCDScale_f%d",nb),nBinMVA, xbins);     histo_DY_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_VV_CMS_QCDScaleBounding[nb]    = new TH1D(Form("histo_VV_QCDScale_f%d",nb),     Form("histo_VV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_VVV_CMS_QCDScaleBounding[nb]   = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),     Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WG_CMS_QCDScaleBounding[nb]    = new TH1D(Form("histo_WG_QCDScale_f%d",nb),     Form("histo_WG_QCDScale_f%d",nb),nBinMVA, xbins);     histo_WG_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WGS_CMS_QCDScaleBounding[nb]   = new TH1D(Form("histo_WGS_QCDScale_f%d",nb),     Form("histo_WGS_QCDScale_f%d",nb),nBinMVA, xbins);     histo_WGS_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_Higgs_CMS_QCDScaleBounding[nb] = new TH1D(Form("histo_Higgs_QCDScale_f%d",nb),      Form("histo_Higgs_QCDScale_f%d",nb),nBinMVA, xbins);      histo_Higgs_CMS_QCDScaleBounding[nb]->Sumw2();
  }
  TH1D* histo_qqWW_CMS_PDFBounding[102];
  TH1D* histo_ggWW_CMS_PDFBounding[102];
  TH1D* histo_Top_CMS_PDFBounding[102];
  TH1D* histo_DY_CMS_PDFBounding[102];
  TH1D* histo_VV_CMS_PDFBounding[102];
  TH1D* histo_VVV_CMS_PDFBounding[102];
  TH1D* histo_WG_CMS_PDFBounding[102];
  TH1D* histo_WGS_CMS_PDFBounding[102];
  TH1D* histo_Higgs_CMS_PDFBounding[102];
  for(int nb=0; nb<102; nb++){
    histo_qqWW_CMS_PDFBounding[nb]   = new TH1D(Form("histo_qqWW_PDF_f%d",nb),     Form("histo_qqWW_PDF_f%d",nb),    nBinMVA, xbins); histo_qqWW_CMS_PDFBounding[nb]->Sumw2();
    histo_ggWW_CMS_PDFBounding[nb]   = new TH1D(Form("histo_ggWW_PDF_f%d",nb),     Form("histo_ggWW_PDF_f%d",nb),    nBinMVA, xbins); histo_ggWW_CMS_PDFBounding[nb]->Sumw2();
    histo_Top_CMS_PDFBounding[nb]    = new TH1D(Form("histo_Top_PDF_f%d",nb),     Form("histo_Top_PDF_f%d",nb),    nBinMVA, xbins); histo_Top_CMS_PDFBounding[nb]->Sumw2();
    histo_DY_CMS_PDFBounding[nb]     = new TH1D(Form("histo_DY_PDF_f%d",nb),     Form("histo_DY_PDF_f%d",nb),    nBinMVA, xbins); histo_DY_CMS_PDFBounding[nb]->Sumw2();
    histo_VV_CMS_PDFBounding[nb]     = new TH1D(Form("histo_VV_PDF_f%d",nb),     Form("histo_VV_PDF_f%d",nb),    nBinMVA, xbins); histo_VV_CMS_PDFBounding[nb]->Sumw2();
    histo_VVV_CMS_PDFBounding[nb]    = new TH1D(Form("histo_VVV_PDF_f%d",nb),     Form("histo_VVV_PDF_f%d",nb),    nBinMVA, xbins); histo_VVV_CMS_PDFBounding[nb]->Sumw2();
    histo_WG_CMS_PDFBounding[nb]     = new TH1D(Form("histo_WG_PDF_f%d",nb),     Form("histo_WG_PDF_f%d",nb),    nBinMVA, xbins); histo_WG_CMS_PDFBounding[nb]->Sumw2();
    histo_WGS_CMS_PDFBounding[nb]    = new TH1D(Form("histo_WGS_PDF_f%d",nb),     Form("histo_WGS_PDF_f%d",nb),    nBinMVA, xbins); histo_WGS_CMS_PDFBounding[nb]->Sumw2();
    histo_Higgs_CMS_PDFBounding[nb]  = new TH1D(Form("histo_Higgs_PDF_f%d",nb),     Form("histo_Higgs_PDF_f%d",nb),    nBinMVA, xbins); histo_Higgs_CMS_PDFBounding[nb]->Sumw2();
  }

  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[nBinMVA];
  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[nBinMVA];
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingBinUp[nBinMVA];
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingBinDown[nBinMVA];
  TH1D* histo_Top_CMS_MVATopStatBoundingBinUp[nBinMVA];
  TH1D* histo_Top_CMS_MVATopStatBoundingBinDown[nBinMVA];
  TH1D* histo_DY_CMS_MVADYStatBoundingBinUp[nBinMVA];
  TH1D* histo_DY_CMS_MVADYStatBoundingBinDown[nBinMVA];
  TH1D* histo_VV_CMS_MVAVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VV_CMS_MVAVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_WG_CMS_MVAWGStatBoundingBinUp[nBinMVA];
  TH1D* histo_WG_CMS_MVAWGStatBoundingBinDown[nBinMVA];
  TH1D* histo_WGS_CMS_MVAWGSStatBoundingBinUp[nBinMVA];
  TH1D* histo_WGS_CMS_MVAWGSStatBoundingBinDown[nBinMVA];
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[nBinMVA];
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[nBinMVA];
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[nBinMVA];
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++){
    histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[nb]	    = new TH1D(Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[nb]      ->Sumw2();
    histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[nb]   = new TH1D(Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[nb]    ->Sumw2();
    histo_ggWW_CMS_MVAggWWStatBoundingBinUp[nb]	    = new TH1D(Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingBinUp[nb]      ->Sumw2();
    histo_ggWW_CMS_MVAggWWStatBoundingBinDown[nb]   = new TH1D(Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingBinDown[nb]    ->Sumw2();
    histo_Top_CMS_MVATopStatBoundingBinUp[nb]	    = new TH1D(Form("histo_Top_CMS_ww%s_MVATopStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Top_CMS_ww%s_MVATopStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingBinUp[nb]      ->Sumw2();
    histo_Top_CMS_MVATopStatBoundingBinDown[nb]	    = new TH1D(Form("histo_Top_CMS_ww%s_MVATopStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_Top_CMS_ww%s_MVATopStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingBinDown[nb]    ->Sumw2();
    histo_DY_CMS_MVADYStatBoundingBinUp[nb]	    = new TH1D(Form("histo_DY_CMS_ww%s_MVADYStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_DY_CMS_ww%s_MVADYStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_DY_CMS_MVADYStatBoundingBinUp[nb]      ->Sumw2();
    histo_DY_CMS_MVADYStatBoundingBinDown[nb]	    = new TH1D(Form("histo_DY_CMS_ww%s_MVADYStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_DY_CMS_ww%s_MVADYStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_DY_CMS_MVADYStatBoundingBinDown[nb]    ->Sumw2();
    histo_VV_CMS_MVAVVStatBoundingBinUp[nb]	    = new TH1D(Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VV_CMS_MVAVVStatBoundingBinDown[nb]	    = new TH1D(Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	    = new TH1D(Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	    = new TH1D(Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WG_CMS_MVAWGStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WG_CMS_MVAWGStatBoundingBinUp[nb]      ->Sumw2();
    histo_WG_CMS_MVAWGStatBoundingBinDown[nb]	    = new TH1D(Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WG_CMS_MVAWGStatBoundingBinDown[nb]    ->Sumw2();
    histo_WGS_CMS_MVAWGSStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WGS_CMS_MVAWGSStatBoundingBinUp[nb]      ->Sumw2();
    histo_WGS_CMS_MVAWGSStatBoundingBinDown[nb]	    = new TH1D(Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WGS_CMS_MVAWGSStatBoundingBinDown[nb]    ->Sumw2();
    histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[nb]   = new TH1D(Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[nb]      ->Sumw2();
    histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[nb] = new TH1D(Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[nb]    ->Sumw2();
    histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[nb]   = new TH1D(Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[nb]      ->Sumw2();
    histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[nb] = new TH1D(Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[nb]    ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]   = new TH1D(Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]      ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb] = new TH1D(Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb]    ->Sumw2();
  }

  TH1D* histo_qqWW_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_qqWW_%sUp",effMName)  , Form("histo_qqWW_%sUp",effMName)  , nBinMVA, xbins); histo_qqWW_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_qqWW_%sDown",effMName), Form("histo_qqWW_%sDown",effMName), nBinMVA, xbins); histo_qqWW_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_ggWW_%sUp",effMName)  , Form("histo_ggWW_%sUp",effMName)  , nBinMVA, xbins); histo_ggWW_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_ggWW_%sDown",effMName), Form("histo_ggWW_%sDown",effMName), nBinMVA, xbins); histo_ggWW_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_Top_%sUp",effMName)  , Form("histo_Top_%sUp",effMName)  , nBinMVA, xbins); histo_Top_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_Top_%sDown",effMName), Form("histo_Top_%sDown",effMName), nBinMVA, xbins); histo_Top_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_DY_%sUp",effMName)  , Form("histo_DY_%sUp",effMName)  , nBinMVA, xbins); histo_DY_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_DY_%sDown",effMName), Form("histo_DY_%sDown",effMName), nBinMVA, xbins); histo_DY_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_VV_%sUp",effMName)  , Form("histo_VV_%sUp",effMName)  , nBinMVA, xbins); histo_VV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_VV_%sDown",effMName), Form("histo_VV_%sDown",effMName), nBinMVA, xbins); histo_VV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effMName)  , Form("histo_VVV_%sUp",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effMName), Form("histo_VVV_%sDown",effMName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_WG_%sUp",effMName)  , Form("histo_WG_%sUp",effMName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_WG_%sDown",effMName), Form("histo_WG_%sDown",effMName), nBinMVA, xbins); histo_WG_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_WGS_%sUp",effMName)  , Form("histo_WGS_%sUp",effMName)  , nBinMVA, xbins); histo_WGS_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_WGS_%sDown",effMName), Form("histo_WGS_%sDown",effMName), nBinMVA, xbins); histo_WGS_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_Higgs_%sUp",effMName)  , Form("histo_Higgs_%sUp",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_Higgs_%sDown",effMName), Form("histo_Higgs_%sDown",effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_qqWW_%sAvg",effMName)  , Form("histo_qqWW_%sAvg",effMName)  , nBinMVA, xbins); histo_qqWW_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_ggWW_%sAvg",effMName)  , Form("histo_ggWW_%sAvg",effMName)  , nBinMVA, xbins); histo_ggWW_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_Top_%sAvg",effMName)  , Form("histo_Top_%sAvg",effMName)  , nBinMVA, xbins); histo_Top_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_DY_%sAvg",effMName)  , Form("histo_DY_%sAvg",effMName)  , nBinMVA, xbins); histo_DY_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_VV_%sAvg",effMName)  , Form("histo_VV_%sAvg",effMName)  , nBinMVA, xbins); histo_VV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effMName)  , Form("histo_VVV_%sAvg",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_WG_%sAvg",effMName)  , Form("histo_WG_%sAvg",effMName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_WGS_%sAvg",effMName)  , Form("histo_WGS_%sAvg",effMName)  , nBinMVA, xbins); histo_WGS_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingAvg    = new TH1D( Form("histo_Higgs_%sAvg",effMName)  , Form("histo_Higgs_%sAvg",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingAvg  ->Sumw2();

  TH1D* histo_qqWW_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_qqWW_%sUp",effEName)  , Form("histo_qqWW_%sUp",effEName)  , nBinMVA, xbins); histo_qqWW_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_qqWW_%sDown",effEName), Form("histo_qqWW_%sDown",effEName), nBinMVA, xbins); histo_qqWW_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_ggWW_%sUp",effEName)  , Form("histo_ggWW_%sUp",effEName)  , nBinMVA, xbins); histo_ggWW_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_ggWW_%sDown",effEName), Form("histo_ggWW_%sDown",effEName), nBinMVA, xbins); histo_ggWW_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_Top_%sUp",effEName)  , Form("histo_Top_%sUp",effEName)  , nBinMVA, xbins); histo_Top_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_Top_%sDown",effEName), Form("histo_Top_%sDown",effEName), nBinMVA, xbins); histo_Top_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_DY_%sUp",effEName)  , Form("histo_DY_%sUp",effEName)  , nBinMVA, xbins); histo_DY_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_DY_%sDown",effEName), Form("histo_DY_%sDown",effEName), nBinMVA, xbins); histo_DY_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_VV_%sUp",effEName)  , Form("histo_VV_%sUp",effEName)  , nBinMVA, xbins); histo_VV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_VV_%sDown",effEName), Form("histo_VV_%sDown",effEName), nBinMVA, xbins); histo_VV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effEName)  , Form("histo_VVV_%sUp",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effEName), Form("histo_VVV_%sDown",effEName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_WG_%sUp",effEName)  , Form("histo_WG_%sUp",effEName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_WG_%sDown",effEName), Form("histo_WG_%sDown",effEName), nBinMVA, xbins); histo_WG_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_WGS_%sUp",effEName)  , Form("histo_WGS_%sUp",effEName)  , nBinMVA, xbins); histo_WGS_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_WGS_%sDown",effEName), Form("histo_WGS_%sDown",effEName), nBinMVA, xbins); histo_WGS_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_Higgs_%sUp",effEName)  , Form("histo_Higgs_%sUp",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_Higgs_%sDown",effEName), Form("histo_Higgs_%sDown",effEName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_qqWW_%sAvg",effEName)  , Form("histo_qqWW_%sAvg",effEName)  , nBinMVA, xbins); histo_qqWW_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_ggWW_%sAvg",effEName)  , Form("histo_ggWW_%sAvg",effEName)  , nBinMVA, xbins); histo_ggWW_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_Top_%sAvg",effEName)  , Form("histo_Top_%sAvg",effEName)  , nBinMVA, xbins); histo_Top_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_DY_%sAvg",effEName)  , Form("histo_DY_%sAvg",effEName)  , nBinMVA, xbins); histo_DY_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_VV_%sAvg",effEName)  , Form("histo_VV_%sAvg",effEName)  , nBinMVA, xbins); histo_VV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effEName)  , Form("histo_VVV_%sAvg",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_WG_%sAvg",effEName)  , Form("histo_WG_%sAvg",effEName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_WGS_%sAvg",effEName)  , Form("histo_WGS_%sAvg",effEName)  , nBinMVA, xbins); histo_WGS_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingAvg    = new TH1D( Form("histo_Higgs_%sAvg",effEName)  , Form("histo_Higgs_%sAvg",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingAvg  ->Sumw2();

  TH1D* histo_qqWW_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_qqWW_CMS_scale_metUp")  , Form("histo_qqWW_CMS_scale_metUp")  , nBinMVA, xbins); histo_qqWW_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_qqWW_CMS_scale_metDown"), Form("histo_qqWW_CMS_scale_metDown"), nBinMVA, xbins); histo_qqWW_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_ggWW_CMS_scale_metUp")  , Form("histo_ggWW_CMS_scale_metUp")  , nBinMVA, xbins); histo_ggWW_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_ggWW_CMS_scale_metDown"), Form("histo_ggWW_CMS_scale_metDown"), nBinMVA, xbins); histo_ggWW_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_Top_CMS_scale_metUp")  , Form("histo_Top_CMS_scale_metUp")  , nBinMVA, xbins); histo_Top_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_Top_CMS_scale_metDown"), Form("histo_Top_CMS_scale_metDown"), nBinMVA, xbins); histo_Top_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_DY_CMS_scale_metUp")  , Form("histo_DY_CMS_scale_metUp")  , nBinMVA, xbins); histo_DY_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_DY_CMS_scale_metDown"), Form("histo_DY_CMS_scale_metDown"), nBinMVA, xbins); histo_DY_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_VV_CMS_scale_metUp")  , Form("histo_VV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_VV_CMS_scale_metDown"), Form("histo_VV_CMS_scale_metDown"), nBinMVA, xbins); histo_VV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown"), nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_WG_CMS_scale_metUp")  , Form("histo_WG_CMS_scale_metUp")  , nBinMVA, xbins); histo_WG_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_WG_CMS_scale_metDown"), Form("histo_WG_CMS_scale_metDown"), nBinMVA, xbins); histo_WG_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_WGS_CMS_scale_metUp")  , Form("histo_WGS_CMS_scale_metUp")  , nBinMVA, xbins); histo_WGS_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_WGS_CMS_scale_metDown"), Form("histo_WGS_CMS_scale_metDown"), nBinMVA, xbins); histo_WGS_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_Higgs_CMS_scale_metUp")  , Form("histo_Higgs_CMS_scale_metUp")  , nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_Higgs_CMS_scale_metDown"), Form("histo_Higgs_CMS_scale_metDown"), nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_qqWW_CMS_scale_jUp")  , Form("histo_qqWW_CMS_scale_jUp")  , nBinMVA, xbins); histo_qqWW_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_qqWW_CMS_scale_jDown"), Form("histo_qqWW_CMS_scale_jDown"), nBinMVA, xbins); histo_qqWW_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_ggWW_CMS_scale_jUp")  , Form("histo_ggWW_CMS_scale_jUp")  , nBinMVA, xbins); histo_ggWW_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_ggWW_CMS_scale_jDown"), Form("histo_ggWW_CMS_scale_jDown"), nBinMVA, xbins); histo_ggWW_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_Top_CMS_scale_jUp")  , Form("histo_Top_CMS_scale_jUp")  , nBinMVA, xbins); histo_Top_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_Top_CMS_scale_jDown"), Form("histo_Top_CMS_scale_jDown"), nBinMVA, xbins); histo_Top_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_DY_CMS_scale_jUp")  , Form("histo_DY_CMS_scale_jUp")  , nBinMVA, xbins); histo_DY_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_DY_CMS_scale_jDown"), Form("histo_DY_CMS_scale_jDown"), nBinMVA, xbins); histo_DY_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_VV_CMS_scale_jUp")  , Form("histo_VV_CMS_scale_jUp")  , nBinMVA, xbins); histo_VV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_VV_CMS_scale_jDown"), Form("histo_VV_CMS_scale_jDown"), nBinMVA, xbins); histo_VV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_VVV_CMS_scale_jUp")  , Form("histo_VVV_CMS_scale_jUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_VVV_CMS_scale_jDown"), Form("histo_VVV_CMS_scale_jDown"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_WG_CMS_scale_jUp")  , Form("histo_WG_CMS_scale_jUp")  , nBinMVA, xbins); histo_WG_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_WG_CMS_scale_jDown"), Form("histo_WG_CMS_scale_jDown"), nBinMVA, xbins); histo_WG_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_WGS_CMS_scale_jUp")  , Form("histo_WGS_CMS_scale_jUp")  , nBinMVA, xbins); histo_WGS_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_WGS_CMS_scale_jDown"), Form("histo_WGS_CMS_scale_jDown"), nBinMVA, xbins); histo_WGS_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_Higgs_CMS_scale_jUp")  , Form("histo_Higgs_CMS_scale_jUp")  , nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_Higgs_CMS_scale_jDown"), Form("histo_Higgs_CMS_scale_jDown"), nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_PUBoundingUp   	         = new TH1D( Form("histo_qqWW_CMS_puUp")  , Form("histo_qqWW_CMS_puUp")  , nBinMVA, xbins); histo_qqWW_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_PUBoundingDown 	         = new TH1D( Form("histo_qqWW_CMS_puDown"), Form("histo_qqWW_CMS_puDown"), nBinMVA, xbins); histo_qqWW_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_PUBoundingUp   	         = new TH1D( Form("histo_ggWW_CMS_puUp")  , Form("histo_ggWW_CMS_puUp")  , nBinMVA, xbins); histo_ggWW_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_PUBoundingDown 	         = new TH1D( Form("histo_ggWW_CMS_puDown"), Form("histo_ggWW_CMS_puDown"), nBinMVA, xbins); histo_ggWW_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_PUBoundingUp   	         = new TH1D( Form("histo_Top_CMS_puUp")  , Form("histo_Top_CMS_puUp")  , nBinMVA, xbins); histo_Top_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_PUBoundingDown 	         = new TH1D( Form("histo_Top_CMS_puDown"), Form("histo_Top_CMS_puDown"), nBinMVA, xbins); histo_Top_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_PUBoundingUp   	         = new TH1D( Form("histo_DY_CMS_puUp")  , Form("histo_DY_CMS_puUp")  , nBinMVA, xbins); histo_DY_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_PUBoundingDown 	         = new TH1D( Form("histo_DY_CMS_puDown"), Form("histo_DY_CMS_puDown"), nBinMVA, xbins); histo_DY_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_PUBoundingUp   	         = new TH1D( Form("histo_VV_CMS_puUp")  , Form("histo_VV_CMS_puUp")  , nBinMVA, xbins); histo_VV_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_PUBoundingDown 	         = new TH1D( Form("histo_VV_CMS_puDown"), Form("histo_VV_CMS_puDown"), nBinMVA, xbins); histo_VV_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingUp   	         = new TH1D( Form("histo_VVV_CMS_puUp")  , Form("histo_VVV_CMS_puUp")  , nBinMVA, xbins); histo_VVV_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingDown 	         = new TH1D( Form("histo_VVV_CMS_puDown"), Form("histo_VVV_CMS_puDown"), nBinMVA, xbins); histo_VVV_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_PUBoundingUp   	         = new TH1D( Form("histo_WG_CMS_puUp")  , Form("histo_WG_CMS_puUp")  , nBinMVA, xbins); histo_WG_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_PUBoundingDown 	         = new TH1D( Form("histo_WG_CMS_puDown"), Form("histo_WG_CMS_puDown"), nBinMVA, xbins); histo_WG_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_PUBoundingUp   	         = new TH1D( Form("histo_WGS_CMS_puUp")  , Form("histo_WGS_CMS_puUp")  , nBinMVA, xbins); histo_WGS_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_PUBoundingDown 	         = new TH1D( Form("histo_WGS_CMS_puDown"), Form("histo_WGS_CMS_puDown"), nBinMVA, xbins); histo_WGS_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingUp   	         = new TH1D( Form("histo_Higgs_CMS_puUp")  , Form("histo_Higgs_CMS_puUp")  , nBinMVA, xbins); histo_Higgs_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingDown 	         = new TH1D( Form("histo_Higgs_CMS_puDown"), Form("histo_Higgs_CMS_puDown"), nBinMVA, xbins); histo_Higgs_CMS_PUBoundingDown->Sumw2();

  double bgdDecay[nSelTypes*4][histBins],weiDecay[nSelTypes*2][histBins];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(int j=0; j<histBins; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }

  const int numberCuts = 10;
  TH1D* hDWWLL[3];
  hDWWLL[0] = new TH1D("hDWWLL0", "hDWWLL0", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[1] = new TH1D("hDWWLL1", "hDWWLL1", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[2] = new TH1D("hDWWLL2", "hDWWLL2", numberCuts, -0.5, numberCuts-0.5);
  TString cutName[numberCuts] = {"ptl>20/20", "3rd lepton veto", "mll>12", "MET>20", "Z veto", "minPMET>20/45", "ptll>30/45", "btag-veto", "Nsoftmu=0", "Njets=0"};
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
 
    if(useDYMVA == true){
      the_input_tree->SetBranchAddress("dymva", &dymva_);
      the_input_tree->SetBranchAddress("nlep", &nlep_ );
      the_input_tree->SetBranchAddress("njets", &njets_);
    }

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
    if(infilecatv[ifile] != 0 && initPDFTag == -1 && infilenamev[ifile].Contains("powheg") == false) {
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

    hDWWLL[0]->Scale(0.0);hDWWLL[1]->Scale(0.0);hDWWLL[2]->Scale(0.0);
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_input_tree->GetEntry(i);

      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      int sumEvol[3] = {-1, -1, -1};
      Bool_t passFilter[numberCuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      Bool_t passPresel = kFALSE;
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
           ) passPresel = kTRUE;
      }
      //if(infilecatv[ifile] != 0) passPresel = kTRUE; // do not apply trigger filters to MC

      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}
      }

      if(idLep.size()!=idTight.size()) {assert(1); return;}

      if(idLep.size() >= 2) passPresel = passPresel && kTRUE;
      else                  passPresel = passPresel && kFALSE;

      if(passPresel == kFALSE) continue;

      int typeSel = -1;
      if(idTight.size() >= 2){
        if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) typeSel = 0;
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) typeSel = 1;
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])!=TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]))         typeSel = 2;
        else {assert(1); return;}                                                                                                     ;
      }
      else {assert(1); printf("Not possible %d %d %d\n",(int)idLep.size(),(int)idTight.size(),goodIsTight); return;}                                                                                                     ;

      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 20 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20) passFilter[0] = kTRUE;

      if((goodIsTight == idTight.size() || usePureMC == false) && idTight.size() == 2) passFilter[1] = kTRUE;

      TLorentzVector theMET;
      if(usePUPPI){
        theMET.SetPx((double)eventMet.metPuppi->Px());
        theMET.SetPy((double)eventMet.metPuppi->Py());
      } else {
        theMET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px());
        theMET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py());
      }
      double dPhiLepMETMin = 999.;
      int signQ = 0;
      double systTotLep[2] = {1.0, 1.0}; // m/e
      for(unsigned nl=0; nl<idLep.size(); nl++){
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) systTotLep[0] = systTotLep[0] * 1.02;
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 11) systTotLep[1] = systTotLep[1] * 1.02;
        signQ = signQ + (int)(*eventLeptons.pdgId)[idLep[nl]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]);
        if(dPhiLepMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(theMET)))
           dPhiLepMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(theMET));      
      }
      double minMET  = TMath::Min(theMET.Pt(),(double)eventMet.trackMet->Pt());
      double minPMET[3] = {TMath::Min(theMET.Pt(),(double)eventMet.trackMet->Pt()),
                           TMath::Min((double)(*eventMet.ptJESUP)[0],(double)eventMet.trackMet->Pt()),
                           TMath::Min((double)(*eventMet.ptJESDOWN)[0],(double)eventMet.trackMet->Pt())};
      if(usePUPPI){ // no changes so far!
        minPMET[1] = TMath::Min((double)(*eventMet.ptJESUP)[0],(double)eventMet.trackMet->Pt());
	minPMET[2] = TMath::Min((double)(*eventMet.ptJESDOWN)[0],(double)eventMet.trackMet->Pt());
      }
      if(dPhiLepMETMin < TMath::Pi()/2) {minPMET[0] = minPMET[0] * sin(dPhiLepMETMin);minPMET[1] = minPMET[1] * sin(dPhiLepMETMin);minPMET[2] = minPMET[2] * sin(dPhiLepMETMin);}

      passFilter[0] = passFilter[0] * (signQ == 0);

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 
      if(dilep.M() > 12) passFilter[2] = kTRUE;

      vector<int> idJet,idJetUp,idJetDown;
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;        

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(dPhiJetMET   == -1) dPhiJetMET   = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(theMET))*180./TMath::Pi();
        if(dPhiJetDiLep == -1) dPhiJetDiLep = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventJets.p4)[nj])))*180./TMath::Pi();

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 15 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()      > 30) idJet.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*1.05 > 30) idJetUp.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*0.95 > 30) idJetDown.push_back(nj);
      }

      if(useDYMVA == true){
        if(nlep_ != idLep.size()) {printf("PROBLEM nlep %d != %d\n",(int)nlep_,(int)idLep.size()); assert(1); return;}
        if(njets_ != idJet.size()) {printf("PROBLEM njet %d != %d\n",(int)njets_,(int)idJet.size()); assert(1); return;}
      }

      int typePair = 0;
      if     (idLep.size() == 2 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])) typePair = 1;

      if(theMET.Pt() > 20) passFilter[3] = kTRUE;
      passFilter[4] = kTRUE;
      if(minPMET[0] > 20) passFilter[5] = kTRUE;
      if(dilep.Pt() > 30) passFilter[6] = kTRUE;
      if(typePair == 1){
        if(TMath::Abs(dilep.M()-91.1876) <= 15.0) passFilter[4] = kFALSE;
        if     (useDYMVA == false && minPMET[0] <= 45) passFilter[5] = kFALSE;
        else if(useDYMVA == true && dymva_ < 0.3) passFilter[5] = kFALSE;
        if(dilep.Pt() <= 45) passFilter[6] = kFALSE;
      }
      if(bDiscrMax < 0.560) passFilter[7] = kTRUE;
      if(idSoft.size() == 0) passFilter[8] = kTRUE;
      if(idJet.size() == nJetsType) passFilter[9] = kTRUE;

      bool totalSel = kTRUE;
      for(int isel=0; isel<numberCuts; isel++) {
        totalSel = totalSel && passFilter[isel];
	if(totalSel == kTRUE) sumEvol[typeSel]++;
      }
      bool passAllCuts[1] = {totalSel};

      bool passSystCuts[nSystTypes] = {
        passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5]                                           && passFilter[6] && passFilter[7] && passFilter[8] && idJetUp.size()   == nJetsType,
        passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5]                                           && passFilter[6] && passFilter[7] && passFilter[8] && idJetDown.size() == nJetsType,
        passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && (minPMET[1] > 20 && (minPMET[1] > 45 || typePair == 0)) && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9],
        passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && (minPMET[2] > 20 && (minPMET[2] > 45 || typePair == 0)) && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9]};
	
      // begin event weighting
      vector<bool> isGenDupl;
      vector<int>wBoson;
      vector<int>zBoson;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
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
      double thePtwwWeight[5] = {1.0,1.0,1.0,1.0,1.0};
      if(infilecatv[ifile] == 1 && wBoson.size() == 2){
        TLorentzVector wwSystem(( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(wBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(wBoson[1])) ) )); 
        Int_t nptwwbin[5] = {fhDWWPtRatio	   ->GetXaxis()->FindBin(TMath::Min(wwSystem.Pt(),499.999)),
	                     fhDWWPtRatio_scaleup  ->GetXaxis()->FindBin(TMath::Min(wwSystem.Pt(),499.999)),
	                     fhDWWPtRatio_scaledown->GetXaxis()->FindBin(TMath::Min(wwSystem.Pt(),499.999)),
	                     fhDWWPtRatio_resumup  ->GetXaxis()->FindBin(TMath::Min(wwSystem.Pt(),499.999)),
	                     fhDWWPtRatio_resumdown->GetXaxis()->FindBin(TMath::Min(wwSystem.Pt(),499.999))};
        thePtwwWeight[0] = fhDWWPtRatio 	 ->GetBinContent(nptwwbin[0]);
	thePtwwWeight[1] = fhDWWPtRatio_scaleup  ->GetBinContent(nptwwbin[1]);
	thePtwwWeight[2] = fhDWWPtRatio_scaledown->GetBinContent(nptwwbin[2]);
	thePtwwWeight[3] = fhDWWPtRatio_resumup  ->GetBinContent(nptwwbin[3]);
	thePtwwWeight[4] = fhDWWPtRatio_resumdown->GetBinContent(nptwwbin[4]);
      }
      vector<int> isGenLep; unsigned int goodIsGenLep = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        bool isGenLepton = false;
        for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	  if(isGenDupl[ngen] == 1) continue;
          if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.3) {
	    isGenLepton = true;
	    break;
	  }
	}
	if(isGenLepton == true) {isGenLep.push_back(1); goodIsGenLep++;}
	else                    {isGenLep.push_back(0);}
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
		typeLepSel.Data(),fhDMuMediumSF,fhDMuIsoSF,fhDElMediumSF,fhDElTightSF,fhDElMediumMVASF,fhDElTightMVASF);
        }
      }

      // fake rate
      int type2l = -1;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == 13 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]) == 13) type2l = 0;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == 13 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]) == 11) type2l = 1;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == 11 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]) == 13) type2l = 2;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == 11 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]) == 11) type2l = 3;
      else {printf("IMPOSSIBLE TYPE2L\n");}
      int nFakeCount = 0;
      unsigned int typeFakeLepton[2] = {0,0};
      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false){
        if     (theCategory == 9){ // remove W+jets from MC
          fakeSF = 0.0;
        }
        else if(theCategory == 4 && goodIsTight != idTight.size()){ // remove Z+jets from MC as fakeable objects
          fakeSF = 0.0;
        }
        else if((infilecatv[ifile] == 0 || infilecatv[ifile] == 7 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add W+jets from data or W+gamma from MC
	  for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) nFakeCount = nFakeCount + 1;
	    else                                                        nFakeCount = nFakeCount + 2;

	    fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    theCategory = 9;
	    if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) typeFakeLepton[0]++;
	    else                                                        typeFakeLepton[1]++;
          }
          if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) fakeSF =  1.0 * fakeSF; // double fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) fakeSF = -1.0 * fakeSF; // single fake, MC
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) fakeSF = -1.0 * fakeSF; // double fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) fakeSF =  1.0 * fakeSF; // single fake, data
        }
        else if(infilecatv[ifile] != 0 && infilecatv[ifile] != 7 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
	  fakeSF = 1.0;
        }
        else if(infilecatv[ifile] == 0 || infilecatv[ifile] == 7){ // data or W+gamma with all good leptons
	  fakeSF = 1.0;
        }
	else {
	  printf("PROBLEM: %d %d %d %d %d\n",infilecatv[ifile],goodIsGenLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
	  assert(0);
	}
      }
      if(typeFakeLepton[0] < typeFakeLepton[1]) theCategory = 10;
      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale*trigEff*thePtwwWeight[0];
      //double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale*trigEff;

      // top-quark estimation
      if(infilecatv[ifile] == 3) totalWeight = totalWeight * topNorm[TMath::Min((int)nJetsType,2)];

      // z pt correction
      if(infilecatv[ifile] == 4 && zBoson.size() == 1) totalWeight = totalWeight * zpt_correction(((TLorentzVector*)(*eventMonteCarlo.p4)[zBoson[0]])->Pt(), 0);

      if(totalWeight == 0) continue;
      // end event weighting

      if(passAllCuts[0] && infilecatv[ifile] == 0) totalFakeDataCount[type2l][nFakeCount] = totalFakeDataCount[type2l][nFakeCount] + 1;

      for(int nl=0; nl <=sumEvol[typeSel]; nl++) if(fakeSF == 1) hDWWLL[typeSel]->Fill((double)nl,totalWeight);

      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i]) {
          bgdDecay[i+typePair*nSelTypes][theCategory] += totalWeight;
          weiDecay[i+typePair*nSelTypes][theCategory] += totalWeight*totalWeight;
        }
      }

      if(1) {
	double MVAVar = (double)typePair;
	//double MVAVar = 0;
        if     (theCategory == 0){
	  if(passAllCuts[SIGSEL]) histo_Data->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 1){ // qqWW
	  if(passAllCuts[SIGSEL]) {
	     histo_qqWW->Fill(MVAVar,totalWeight);

             histo_qqWW_CMS_MVAWW_nlo      ->Fill(MVAVar,totalWeight);
             histo_qqWW_CMS_MVAWW_qup_nlo  ->Fill(MVAVar,totalWeight*thePtwwWeight[1]);
             histo_qqWW_CMS_MVAWW_qdown_nlo->Fill(MVAVar,totalWeight*thePtwwWeight[2]);
             histo_qqWW_CMS_MVAWW_sup_nlo  ->Fill(MVAVar,totalWeight*thePtwwWeight[3]);
             histo_qqWW_CMS_MVAWW_sdown_nlo->Fill(MVAVar,totalWeight*thePtwwWeight[4]);

	     histo_qqWW_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_qqWW_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_qqWW_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_qqWW_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_qqWW_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_qqWW_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_qqWW_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_qqWW_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_qqWW_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_qqWW_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_qqWW_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_qqWW_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_qqWW_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_qqWW_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_qqWW_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_qqWW_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_qqWW_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_qqWW_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_qqWW_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_qqWW_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_qqWW_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 2){ // ggWW
	  if(passAllCuts[SIGSEL]) {
	     histo_ggWW->Fill(MVAVar,totalWeight);

	     histo_ggWW_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_ggWW_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_ggWW_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_ggWW_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_ggWW_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_ggWW_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_ggWW_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_ggWW_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_ggWW_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_ggWW_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ggWW_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ggWW_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_ggWW_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_ggWW_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_ggWW_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_ggWW_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ggWW_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_ggWW_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ggWW_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ggWW_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_ggWW_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 3){ // Top
	  if(passAllCuts[SIGSEL]) {
	     histo_Top->Fill(MVAVar,totalWeight);

	     histo_Top_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_Top_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_Top_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_Top_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_Top_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_Top_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_Top_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true && (int)(*eventMonteCarlo.pdfRwgt).size() > 0)
	     for(int npdf=0; npdf<102; npdf++) histo_Top_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_Top_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_Top_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Top_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Top_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Top_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Top_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_Top_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_Top_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_Top_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_Top_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_Top_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_Top_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_Top_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 4){ // DY
	  if(passAllCuts[SIGSEL]) {
	     histo_DY->Fill(MVAVar,totalWeight);

	     histo_DY_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_DY_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_DY_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_DY_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_DY_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_DY_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_DY_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_DY_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_DY_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_DY_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_DY_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_DY_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_DY_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_DY_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_DY_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_DY_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_DY_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_DY_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_DY_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_DY_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_DY_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 5){ // VV
	  if(passAllCuts[SIGSEL]) {
	     histo_VV->Fill(MVAVar,totalWeight);

	     histo_VV_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_VV_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_VV_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_VV_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_VV_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_VV_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_VV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_VV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_VV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_VV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_VV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_VV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_VV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_VV_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_VV_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_VV_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_VV_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_VV_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_VV_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 6){ // VVV
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
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_VVV_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_VVV_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_VVV_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_VVV_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_VVV_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_VVV_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 7){ // WG
	  if(passAllCuts[SIGSEL]) {
	     histo_WG->Fill(MVAVar,totalWeight);

	     histo_WG_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_WG_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_WG_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_WG_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_WG_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_WG_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_WG_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_WG_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_WG_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_WG_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WG_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WG_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_WG_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_WG_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_WG_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_WG_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_WG_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_WG_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_WG_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_WG_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_WG_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 8){ // WGS
	  if(passAllCuts[SIGSEL]) {
	     histo_WGS->Fill(MVAVar,totalWeight);

	     histo_WGS_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_WGS_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_WGS_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_WGS_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_WGS_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_WGS_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_WGS_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_WGS_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_WGS_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_WGS_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WGS_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WGS_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_WGS_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_WGS_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_WGS_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_WGS_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_WGS_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_WGS_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_WGS_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_WGS_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_WGS_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 9){ // WjetsM
	  if(passAllCuts[SIGSEL]) {
	     histo_WjetsM->Fill(MVAVar,totalWeight);
	  }
        }
        else if(theCategory == 10){ // WjetsE
	  if(passAllCuts[SIGSEL]) {
	     histo_WjetsE->Fill(MVAVar,totalWeight);
	  }
        }
        else if(theCategory == 11){ // Higgs
	  if(passAllCuts[SIGSEL]) {
	     histo_Higgs->Fill(MVAVar,totalWeight);

	     histo_Higgs_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_Higgs_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_Higgs_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_Higgs_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_Higgs_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_Higgs_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_Higgs_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Higgs_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Higgs_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Higgs_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Higgs_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_Higgs_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_Higgs_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_Higgs_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[JESUP])  histo_Higgs_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_Higgs_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_Higgs_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_Higgs_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
	else {
	  printf("CATEGORY PROBLEM!\n"); return;
	}
      } // making data cards

    }
    printf("                   mm           ee           em\n");
    for(int nc=0; nc<numberCuts; nc++){
      printf("(%15s): %10.2f   %10.2f   %10.2f\n",cutName[nc].Data(),hDWWLL[0]->GetBinContent(nc+1),hDWWLL[1]->GetBinContent(nc+1),hDWWLL[2]->GetBinContent(nc+1));
    }
  } // end of chain

  printf("----------------------totalFakeDataCount--------------------------------\n");
  for(int ni=0; ni<4; ni++) {
    printf("(%d): ",ni);
    for(int nj=0; nj<4; nj++) printf("%6.1f ",totalFakeDataCount[ni][nj]);
    printf("\n");
  }
  printf("                    em                     ee/mm                     ll\n");
  printf("----------------------------------------------------------------------------------\n");
  for(int ns=0; ns<nSelTypes; ns++) {
    printf("Selection: %s\n",selTypeName[ns].Data());
    double sumEventsBckType[3] = {0,0,0}; double sumEventsBckTypeE[3] = {0,0,0};
    double sumEventsSigType[3] = {0,0,0}; double sumEventsSigTypeE[3] = {0,0,0};
    for(int np=0; np<histBins; np++) {       
       bgdDecay[ns+nSelTypes*2][np] = bgdDecay[ns+nSelTypes*0][np] + bgdDecay[ns+nSelTypes*1][np];
       weiDecay[ns+nSelTypes*2][np] = weiDecay[ns+nSelTypes*0][np] + weiDecay[ns+nSelTypes*1][np];
       printf("(%6s): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
       processName[np].Data(),bgdDecay[ns+nSelTypes*0][np],sqrt(weiDecay[ns+nSelTypes*0][np]),bgdDecay[ns+nSelTypes*1][np],sqrt(weiDecay[ns+nSelTypes*1][np]),
                              bgdDecay[ns+nSelTypes*2][np],sqrt(weiDecay[ns+nSelTypes*2][np]));
       if(np!=0 && np!=1 && np !=2){
         sumEventsBckType[0] = sumEventsBckType[0] + bgdDecay[ns+nSelTypes*0][np]; sumEventsBckTypeE[0] = sumEventsBckTypeE[0] + weiDecay[ns+nSelTypes*0][np];
         sumEventsBckType[1] = sumEventsBckType[1] + bgdDecay[ns+nSelTypes*1][np]; sumEventsBckTypeE[1] = sumEventsBckTypeE[1] + weiDecay[ns+nSelTypes*1][np];
         sumEventsBckType[2] = sumEventsBckType[2] + bgdDecay[ns+nSelTypes*2][np]; sumEventsBckTypeE[2] = sumEventsBckTypeE[2] + weiDecay[ns+nSelTypes*2][np];
       }
       if(np==1 || np==2){
         sumEventsSigType[0] = sumEventsSigType[0] + bgdDecay[ns+nSelTypes*0][np]; sumEventsSigTypeE[0] = sumEventsSigTypeE[0] + weiDecay[ns+nSelTypes*0][np];
         sumEventsSigType[1] = sumEventsSigType[1] + bgdDecay[ns+nSelTypes*1][np]; sumEventsSigTypeE[1] = sumEventsSigTypeE[1] + weiDecay[ns+nSelTypes*1][np];
         sumEventsSigType[2] = sumEventsSigType[2] + bgdDecay[ns+nSelTypes*2][np]; sumEventsSigTypeE[2] = sumEventsSigTypeE[2] + weiDecay[ns+nSelTypes*2][np];
       }
    }
    printf("(...bkg): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
           sumEventsBckType[0],sqrt(sumEventsBckTypeE[0]),sumEventsBckType[1],sqrt(sumEventsBckTypeE[1]),sumEventsBckType[2],sqrt(sumEventsBckTypeE[2]));
    printf("(....WW): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
           sumEventsSigType[0],sqrt(sumEventsSigTypeE[0]),sumEventsSigType[1],sqrt(sumEventsSigTypeE[1]),sumEventsSigType[2],sqrt(sumEventsSigTypeE[2]));
    printf("(.total): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
           sumEventsBckType[0]+sumEventsSigType[0],sqrt(sumEventsBckTypeE[0]+sumEventsSigTypeE[0]),sumEventsBckType[1]+sumEventsSigType[1],sqrt(sumEventsBckTypeE[1]+sumEventsSigTypeE[1]),sumEventsBckType[2]+sumEventsSigType[2],sqrt(sumEventsBckTypeE[2]+sumEventsSigTypeE[2]));
    printf("--------------------------------------------------------------------------------\n");
  }
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    char output[200];
    sprintf(output,"histoww_nice%d_%d.root",nJetsType,thePlot);	  
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }
  printf("QCD Corr: qqWW(%f:%f/%f/%f/%f/%f/%f) ggWW(%f:%f/%f/%f/%f/%f/%f) Top(%f:%f/%f/%f/%f/%f/%f) DY(%f:%f/%f/%f/%f/%f/%f) VV(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) WG(%f:%f/%f/%f/%f/%f/%f) WGS(%f:%f/%f/%f/%f/%f/%f) Higgs(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_qqWW->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ggWW->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_Top->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_DY->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VV->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_WG->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_WGS->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_Higgs->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[5]->GetSumOfWeights());

  printf("ptWW  raw: nlo(%f)up/down(%f/%f/%f/%f)nnlo(%f)\n",histo_qqWW_CMS_MVAWW_nlo->GetSumOfWeights(),
                                     histo_qqWW_CMS_MVAWW_qup_nlo->GetSumOfWeights(),histo_qqWW_CMS_MVAWW_qdown_nlo->GetSumOfWeights(),
                                     histo_qqWW_CMS_MVAWW_sup_nlo->GetSumOfWeights(),histo_qqWW_CMS_MVAWW_sdown_nlo->GetSumOfWeights(),
				     histo_qqWW->GetSumOfWeights());

  double ratio_nnlo_nlo = histo_qqWW->GetSumOfWeights()/histo_qqWW_CMS_MVAWW_nlo->GetSumOfWeights();
  histo_qqWW_CMS_MVAWW_nlo      ->Scale(ratio_nnlo_nlo);
  histo_qqWW_CMS_MVAWW_qup_nlo  ->Scale(ratio_nnlo_nlo);
  histo_qqWW_CMS_MVAWW_qdown_nlo->Scale(ratio_nnlo_nlo);
  histo_qqWW_CMS_MVAWW_sup_nlo  ->Scale(ratio_nnlo_nlo);
  histo_qqWW_CMS_MVAWW_sdown_nlo->Scale(ratio_nnlo_nlo);

  printf("ptWW norm: nlo(%f)up/down(%f/%f/%f/%f)nnlo(%f)\n",histo_qqWW_CMS_MVAWW_nlo->GetSumOfWeights(),
                                     histo_qqWW_CMS_MVAWW_qup_nlo->GetSumOfWeights(),histo_qqWW_CMS_MVAWW_qdown_nlo->GetSumOfWeights(),
                                     histo_qqWW_CMS_MVAWW_sup_nlo->GetSumOfWeights(),histo_qqWW_CMS_MVAWW_sdown_nlo->GetSumOfWeights(),
				     histo_qqWW->GetSumOfWeights());

  for(int i=1; i<=histo_qqWW->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_qqWW_CMS_MVAqqWWStatBoundingUp       ->SetBinContent(i,TMath::Max(histo_qqWW    ->GetBinContent(i)+factorUp  *histo_qqWW    ->GetBinError(i),0.000001));
    histo_qqWW_CMS_MVAqqWWStatBoundingDown     ->SetBinContent(i,TMath::Max(histo_qqWW    ->GetBinContent(i)+factorDown*histo_qqWW    ->GetBinError(i),0.000001));
    histo_ggWW_CMS_MVAggWWStatBoundingUp       ->SetBinContent(i,TMath::Max(histo_ggWW    ->GetBinContent(i)+factorUp  *histo_ggWW    ->GetBinError(i),0.000001));
    histo_ggWW_CMS_MVAggWWStatBoundingDown     ->SetBinContent(i,TMath::Max(histo_ggWW    ->GetBinContent(i)+factorDown*histo_ggWW    ->GetBinError(i),0.000001));
    histo_Top_CMS_MVATopStatBoundingUp         ->SetBinContent(i,TMath::Max(histo_Top	 ->GetBinContent(i)+factorUp  *histo_Top    ->GetBinError(i),0.000001));
    histo_Top_CMS_MVATopStatBoundingDown       ->SetBinContent(i,TMath::Max(histo_Top	 ->GetBinContent(i)+factorDown*histo_Top    ->GetBinError(i),0.000001));
    histo_DY_CMS_MVADYStatBoundingUp           ->SetBinContent(i,TMath::Max(histo_DY    ->GetBinContent(i)+factorUp  *histo_DY    ->GetBinError(i),0.000001));
    histo_DY_CMS_MVADYStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_DY    ->GetBinContent(i)+factorDown*histo_DY    ->GetBinError(i),0.000001));
    histo_VV_CMS_MVAVVStatBoundingUp           ->SetBinContent(i,TMath::Max(histo_VV    ->GetBinContent(i)+factorUp  *histo_VV    ->GetBinError(i),0.000001));
    histo_VV_CMS_MVAVVStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_VV    ->GetBinContent(i)+factorDown*histo_VV    ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp         ->SetBinContent(i,TMath::Max(histo_VVV    ->GetBinContent(i)+factorUp  *histo_VVV    ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown       ->SetBinContent(i,TMath::Max(histo_VVV    ->GetBinContent(i)+factorDown*histo_VVV    ->GetBinError(i),0.000001));
    histo_WG_CMS_MVAWGStatBoundingUp           ->SetBinContent(i,TMath::Max(histo_WG    ->GetBinContent(i)+factorUp  *histo_WG    ->GetBinError(i),0.000001));
    histo_WG_CMS_MVAWGStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_WG    ->GetBinContent(i)+factorDown*histo_WG    ->GetBinError(i),0.000001));
    histo_WGS_CMS_MVAWGSStatBoundingUp         ->SetBinContent(i,TMath::Max(histo_WGS    ->GetBinContent(i)+factorUp  *histo_WGS    ->GetBinError(i),0.000001));
    histo_WGS_CMS_MVAWGSStatBoundingDown       ->SetBinContent(i,TMath::Max(histo_WGS    ->GetBinContent(i)+factorDown*histo_WGS    ->GetBinError(i),0.000001));
    histo_WjetsM_CMS_MVAWjetsMStatBoundingUp   ->SetBinContent(i,TMath::Max(histo_WjetsM    ->GetBinContent(i)+factorUp  *histo_WjetsM    ->GetBinError(i),0.000001));
    histo_WjetsM_CMS_MVAWjetsMStatBoundingDown ->SetBinContent(i,TMath::Max(histo_WjetsM    ->GetBinContent(i)+factorDown*histo_WjetsM    ->GetBinError(i),0.000001));
    histo_WjetsE_CMS_MVAWjetsEStatBoundingUp   ->SetBinContent(i,TMath::Max(histo_WjetsE    ->GetBinContent(i)+factorUp  *histo_WjetsE    ->GetBinError(i),0.000001));
    histo_WjetsE_CMS_MVAWjetsEStatBoundingDown ->SetBinContent(i,TMath::Max(histo_WjetsE    ->GetBinContent(i)+factorDown*histo_WjetsE    ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingUp     ->SetBinContent(i,TMath::Max(histo_Higgs    ->GetBinContent(i)+factorUp  *histo_Higgs    ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingDown   ->SetBinContent(i,TMath::Max(histo_Higgs    ->GetBinContent(i)+factorDown*histo_Higgs    ->GetBinError(i),0.000001));
  }
  char outputLimits[200];
  sprintf(outputLimits,"ww%4s_%dj_input_%s.root",finalStateName,nJetsType,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data ->Write();
  histo_qqWW ->Write();
  histo_ggWW ->Write();
  histo_Top  ->Write();
  histo_DY   ->Write();
  histo_VV   ->Write();
  histo_VVV  ->Write();
  histo_WG   ->Write();
  histo_WGS  ->Write();
  histo_WjetsM->Write();
  histo_WjetsE->Write();
  histo_Higgs->Write();
  cout << histo_Data  ->GetSumOfWeights() << " ";
  cout << histo_qqWW  ->GetSumOfWeights() << " ";
  cout << histo_ggWW  ->GetSumOfWeights() << " ";
  cout << histo_Top   ->GetSumOfWeights() << " ";
  cout << histo_DY    ->GetSumOfWeights() << " ";
  cout << histo_VV    ->GetSumOfWeights() << " ";
  cout << histo_VVV   ->GetSumOfWeights() << " ";
  cout << histo_WG    ->GetSumOfWeights() << " ";
  cout << histo_WGS   ->GetSumOfWeights() << " ";
  cout << histo_WjetsM->GetSumOfWeights() << " ";
  cout << histo_WjetsE->GetSumOfWeights() << " ";
  cout << histo_Higgs ->GetSumOfWeights() << " ";
  cout << endl;
  printf("uncertainties Stat\n");
  histo_qqWW_CMS_MVAqqWWStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW   ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAqqWWStatBoundingUp	   ->GetBinContent(i)/histo_qqWW  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAqqWWStatBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW   ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAqqWWStatBoundingDown    ->GetBinContent(i)/histo_qqWW  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAggWWStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW   ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAggWWStatBoundingUp	   ->GetBinContent(i)/histo_ggWW  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAggWWStatBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW   ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAggWWStatBoundingDown    ->GetBinContent(i)/histo_ggWW  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVATopStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVATopStatBoundingUp	->GetBinContent(i)/histo_Top  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVATopStatBoundingDown    ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVATopStatBoundingDown	->GetBinContent(i)/histo_Top  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVADYStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY     ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVADYStatBoundingUp      ->GetBinContent(i)/histo_DY  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVADYStatBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY     ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVADYStatBoundingDown    ->GetBinContent(i)/histo_DY  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAVVStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV     ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAVVStatBoundingUp      ->GetBinContent(i)/histo_VV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAVVStatBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV     ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAVVStatBoundingDown    ->GetBinContent(i)/histo_VV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingDown    ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAWGStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG     ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAWGStatBoundingUp      ->GetBinContent(i)/histo_WG  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAWGStatBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG     ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAWGStatBoundingDown    ->GetBinContent(i)/histo_WG  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAWGSStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAWGSStatBoundingUp	->GetBinContent(i)/histo_WGS  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAWGSStatBoundingDown    ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAWGSStatBoundingDown	->GetBinContent(i)/histo_WGS  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsM_CMS_MVAWjetsMStatBoundingUp  ->Write();for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsM ->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsM_CMS_MVAWjetsMStatBoundingUp      ->GetBinContent(i)/histo_WjetsM  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsM_CMS_MVAWjetsMStatBoundingDown->Write();for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsM ->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsM_CMS_MVAWjetsMStatBoundingDown    ->GetBinContent(i)/histo_WjetsM  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsE_CMS_MVAWjetsEStatBoundingUp  ->Write();for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsE ->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsE_CMS_MVAWjetsEStatBoundingUp      ->GetBinContent(i)/histo_WjetsE  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsE_CMS_MVAWjetsEStatBoundingDown->Write();for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsE ->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsE_CMS_MVAWjetsEStatBoundingDown    ->GetBinContent(i)/histo_WjetsE  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAHiggsStatBoundingUp      ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAHiggsStatBoundingDown->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAHiggsStatBoundingDown    ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffM\n");
  histo_qqWW_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_qqWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_qqWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_ggWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_ggWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_Top_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_Top_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_DY_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_DY_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_WG_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_WG_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_WGS_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_WGS_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffMBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffE\n");
  histo_qqWW_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_qqWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_qqWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_ggWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_ggWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_Top_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_Top_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_DY_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_DY_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_WG_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_WG_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_WGS_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_WGS_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffEBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties MET\n");
  histo_qqWW_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JES\n");
  histo_qqWW_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVAJESBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAJESBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAJESBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties PU\n");
  histo_qqWW_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_PUBoundingUp	   ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_PUBoundingDown    ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_PUBoundingUp	   ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_PUBoundingDown    ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_PUBoundingUp	 ->GetBinContent(i)/histo_Top	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_PUBoundingDown	 ->GetBinContent(i)/histo_Top	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_PUBoundingUp      ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_PUBoundingDown    ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_PUBoundingUp      ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_PUBoundingDown    ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingUp	 ->GetBinContent(i)/histo_VVV	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingDown	 ->GetBinContent(i)/histo_VVV	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_PUBoundingUp      ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_PUBoundingDown    ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_PUBoundingUp	 ->GetBinContent(i)/histo_WGS	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_PUBoundingDown	 ->GetBinContent(i)/histo_WGS	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_PUBoundingUp      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_PUBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_PUBoundingDown    ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties WW NNLO\n");
  histo_qqWW_CMS_MVAWW_qup_nlo	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWW_qup_nlo         ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAWW_qdown_nlo	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWW_qdown_nlo       ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAWW_sup_nlo	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWW_sup_nlo         ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAWW_sdown_nlo	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWW_sdown_nlo       ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  outFileLimits->Close();
  double lumiE = 1.027;
  double systLepResE[9] = {1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01};
  double systLepResM[9] = {1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01};
  double syst_btag = 1.02;
  for(int nb=1; nb<=nBinMVA; nb++){
     // QCD study
    double systQCDScale[9] = {TMath::Abs(histo_qqWW_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_qqWW->GetBinContent(nb)),
                              TMath::Abs(histo_ggWW_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_ggWW->GetBinContent(nb)),
                              TMath::Abs(histo_Top_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_Top->GetBinContent(nb)),
                              TMath::Abs(histo_DY_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_DY->GetBinContent(nb)),
                              TMath::Abs(histo_VV_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_VV->GetBinContent(nb)),
                              TMath::Abs(histo_VVV_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_VVV->GetBinContent(nb)),
                              TMath::Abs(histo_WG_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_WG->GetBinContent(nb)),
                              TMath::Abs(histo_WGS_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_WGS->GetBinContent(nb)),
                              TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb))};
    for(int nqcd=1; nqcd<6; nqcd++) {
      if(TMath::Abs(histo_qqWW_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_qqWW->GetBinContent(nb)) >   systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_qqWW_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_qqWW->GetBinContent(nb));
      if(TMath::Abs(histo_ggWW_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_ggWW->GetBinContent(nb)) >   systQCDScale[1]) systQCDScale[1] = TMath::Abs(histo_ggWW_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_ggWW->GetBinContent(nb));
      if(TMath::Abs(histo_Top_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Top->GetBinContent(nb)) >     systQCDScale[2]) systQCDScale[2] = TMath::Abs(histo_Top_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Top->GetBinContent(nb));
      if(TMath::Abs(histo_DY_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_DY->GetBinContent(nb)) >       systQCDScale[3]) systQCDScale[3] = TMath::Abs(histo_DY_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_DY->GetBinContent(nb));
      if(TMath::Abs(histo_VV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VV->GetBinContent(nb)) >       systQCDScale[4]) systQCDScale[4] = TMath::Abs(histo_VV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VV->GetBinContent(nb));
      if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb)) >     systQCDScale[5]) systQCDScale[5] = TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb));
      if(TMath::Abs(histo_WG_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_WG->GetBinContent(nb)) >       systQCDScale[6]) systQCDScale[6] = TMath::Abs(histo_WG_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_WG->GetBinContent(nb));
      if(TMath::Abs(histo_WGS_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_WGS->GetBinContent(nb)) >     systQCDScale[7]) systQCDScale[7] = TMath::Abs(histo_WGS_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_WGS->GetBinContent(nb));
      if(TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb)) > systQCDScale[8]) systQCDScale[8] = TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb));
    }
    if(histo_qqWW ->GetBinContent(nb) != 0) systQCDScale[0] = 1 + systQCDScale[0]/histo_qqWW ->GetBinContent(nb); else systQCDScale[0] = 1;
    if(histo_ggWW ->GetBinContent(nb) != 0) systQCDScale[1] = 1 + systQCDScale[1]/histo_ggWW ->GetBinContent(nb); else systQCDScale[1] = 1;
    if(histo_Top  ->GetBinContent(nb) != 0) systQCDScale[2] = 1 + systQCDScale[2]/histo_Top  ->GetBinContent(nb); else systQCDScale[2] = 1;
    if(histo_DY   ->GetBinContent(nb) != 0) systQCDScale[3] = 1 + systQCDScale[3]/histo_DY   ->GetBinContent(nb); else systQCDScale[3] = 1;
    if(histo_VV   ->GetBinContent(nb) != 0) systQCDScale[4] = 1 + systQCDScale[4]/histo_VV   ->GetBinContent(nb); else systQCDScale[4] = 1;
    if(histo_VVV  ->GetBinContent(nb) != 0) systQCDScale[5] = 1 + systQCDScale[5]/histo_VVV  ->GetBinContent(nb); else systQCDScale[5] = 1;
    if(histo_WG   ->GetBinContent(nb) != 0) systQCDScale[6] = 1 + systQCDScale[6]/histo_WG   ->GetBinContent(nb); else systQCDScale[6] = 1;
    if(histo_WGS  ->GetBinContent(nb) != 0) systQCDScale[7] = 1 + systQCDScale[7]/histo_WGS  ->GetBinContent(nb); else systQCDScale[7] = 1;
    if(histo_Higgs->GetBinContent(nb) != 0) systQCDScale[8] = 1 + systQCDScale[8]/histo_Higgs->GetBinContent(nb); else systQCDScale[8] = 1;
    printf("QCDScale(%d): %f %f %f %f %f %f %f %f %f\n",nb,systQCDScale[0],systQCDScale[1],systQCDScale[2],systQCDScale[3],systQCDScale[4],systQCDScale[5],systQCDScale[6],systQCDScale[7],systQCDScale[8]);
    
    // TOP and DY QCDScale change
    systQCDScale[2] = 1.0 + topNormE[TMath::Min((int)nJetsType,2)]/topNorm[TMath::Min((int)nJetsType,2)];
    systQCDScale[3] = 1.10;

    // PDF study
    double systPDF[9];
    histo_Diff->Reset();
    if(histo_qqWW->GetBinContent(nb) > 0) for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_qqWW_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_qqWW->GetBinContent(nb))/histo_qqWW->GetBinContent(nb));
    systPDF[0] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_ggWW->GetBinContent(nb) > 0)for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_ggWW_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_ggWW->GetBinContent(nb))/histo_ggWW->GetBinContent(nb));
    systPDF[1] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_Top->GetBinContent(nb) > 0)for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_Top_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_Top->GetBinContent(nb))/histo_Top->GetBinContent(nb));
    systPDF[2] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_DY->GetBinContent(nb) > 0)for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_DY_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_DY->GetBinContent(nb))/histo_DY->GetBinContent(nb));
    systPDF[3] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_VV->GetBinContent(nb) > 0)for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_VV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VV->GetBinContent(nb))/histo_VV->GetBinContent(nb));
    systPDF[4] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_VVV->GetBinContent(nb) > 0)for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_VVV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))/histo_VVV->GetBinContent(nb));
    systPDF[5] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_WG->GetBinContent(nb) > 0)for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_WG_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_WG->GetBinContent(nb))/histo_WG->GetBinContent(nb));
    systPDF[6] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_WGS->GetBinContent(nb) > 0)for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_WGS_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_WGS->GetBinContent(nb))/histo_WGS->GetBinContent(nb));
    systPDF[7] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_Higgs->GetBinContent(nb) > 0)for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_Higgs_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb))/histo_Higgs->GetBinContent(nb));
    systPDF[8] = 1.0+histo_Diff->GetRMS();

    printf("PDF(%d): %f %f %f %f %f %f %f %f %f\n",nb,systPDF[0],systPDF[1],systPDF[2],systPDF[3],systPDF[4],systPDF[5],systPDF[6],systPDF[7],systPDF[8]);

    double systLepEffM[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_qqWW_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[0] = histo_qqWW_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_qqWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_qqWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_qqWW_CMS_MVALepEffMBoundingDown   ->GetBinContent(nb) > 0) systLepEffM[0] = histo_qqWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_qqWW_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_ggWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_ggWW_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[1] = histo_ggWW_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ggWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_ggWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_ggWW_CMS_MVALepEffMBoundingDown   ->GetBinContent(nb) > 0) systLepEffM[1] = histo_ggWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ggWW_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_Top_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_Top_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[2] = histo_Top_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_Top_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_Top_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_Top_CMS_MVALepEffMBoundingDown     ->GetBinContent(nb) > 0) systLepEffM[2] = histo_Top_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_Top_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_DY_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_DY_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[3] = histo_DY_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_DY_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_DY_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_DY_CMS_MVALepEffMBoundingDown       ->GetBinContent(nb) > 0) systLepEffM[3] = histo_DY_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_DY_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_VV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VV_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[4] = histo_VV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_VV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VV_CMS_MVALepEffMBoundingDown       ->GetBinContent(nb) > 0) systLepEffM[4] = histo_VV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[5] = histo_VVV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffMBoundingDown     ->GetBinContent(nb) > 0) systLepEffM[5] = histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_WG_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[6] = histo_WG_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_WG_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_WG_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffMBoundingDown       ->GetBinContent(nb) > 0) systLepEffM[6] = histo_WG_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_WG_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_WGS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_WGS_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[7] = histo_WGS_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_WGS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_WGS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_WGS_CMS_MVALepEffMBoundingDown     ->GetBinContent(nb) > 0) systLepEffM[7] = histo_WGS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_WGS_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffMBoundingUp   ->GetBinContent(nb) > 0) systLepEffM[8] = histo_Higgs_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[8] = histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffMBoundingDown->GetBinContent(nb);

    double systLepEffE[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_qqWW_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[0] = histo_qqWW_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_qqWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_qqWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_qqWW_CMS_MVALepEffEBoundingDown   ->GetBinContent(nb) > 0) systLepEffE[0] = histo_qqWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_qqWW_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_ggWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_ggWW_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[1] = histo_ggWW_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_ggWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_ggWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_ggWW_CMS_MVALepEffEBoundingDown   ->GetBinContent(nb) > 0) systLepEffE[1] = histo_ggWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_ggWW_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_Top_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_Top_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[2] = histo_Top_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_Top_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_Top_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_Top_CMS_MVALepEffEBoundingDown     ->GetBinContent(nb) > 0) systLepEffE[2] = histo_Top_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_Top_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_DY_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_DY_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[3] = histo_DY_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_DY_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_DY_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_DY_CMS_MVALepEffEBoundingDown       ->GetBinContent(nb) > 0) systLepEffE[3] = histo_DY_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_DY_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_VV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VV_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[4] = histo_VV_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_VV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_VV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VV_CMS_MVALepEffEBoundingDown       ->GetBinContent(nb) > 0) systLepEffE[4] = histo_VV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_VV_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[5] = histo_VVV_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffEBoundingDown     ->GetBinContent(nb) > 0) systLepEffE[5] = histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_WG_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[6] = histo_WG_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_WG_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_WG_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffEBoundingDown       ->GetBinContent(nb) > 0) systLepEffE[6] = histo_WG_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_WG_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_WGS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_WGS_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[7] = histo_WGS_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_WGS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_WGS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_WGS_CMS_MVALepEffEBoundingDown     ->GetBinContent(nb) > 0) systLepEffE[7] = histo_WGS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_WGS_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffEBoundingUp   ->GetBinContent(nb) > 0) systLepEffE[8] = histo_Higgs_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffEBoundingDown ->GetBinContent(nb) > 0) systLepEffE[8] = histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffEBoundingDown->GetBinContent(nb);

    double systMet[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[0] = histo_qqWW_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAMETBoundingDown   ->GetBinContent(nb) > 0) systMet[0] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[1] = histo_ggWW_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_ggWW->GetBinContent(nb);
    else if(histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVAMETBoundingDown   ->GetBinContent(nb) > 0) systMet[1] = histo_ggWW->GetBinContent(nb)/histo_ggWW_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[2] = histo_Top_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_Top->GetBinContent(nb);
    else if(histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVAMETBoundingDown     ->GetBinContent(nb) > 0) systMet[2] = histo_Top->GetBinContent(nb)/histo_Top_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[3] = histo_DY_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_DY->GetBinContent(nb);
    else if(histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVAMETBoundingDown       ->GetBinContent(nb) > 0) systMet[3] = histo_DY->GetBinContent(nb)/histo_DY_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[4] = histo_VV_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_VV->GetBinContent(nb);
    else if(histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVAMETBoundingDown       ->GetBinContent(nb) > 0) systMet[4] = histo_VV->GetBinContent(nb)/histo_VV_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[5] = histo_VVV_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAMETBoundingDown     ->GetBinContent(nb) > 0) systMet[5] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[6] = histo_WG_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_WG->GetBinContent(nb);
    else if(histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVAMETBoundingDown       ->GetBinContent(nb) > 0) systMet[6] = histo_WG->GetBinContent(nb)/histo_WG_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[7] = histo_WGS_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_WGS->GetBinContent(nb);
    else if(histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVAMETBoundingDown     ->GetBinContent(nb) > 0) systMet[7] = histo_WGS->GetBinContent(nb)/histo_WGS_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAMETBoundingUp   ->GetBinContent(nb) > 0) systMet[8] = histo_Higgs_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAMETBoundingDown ->GetBinContent(nb) > 0) systMet[8] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_MVAMETBoundingDown->GetBinContent(nb);

    double systJes[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[0] = histo_qqWW_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAJESBoundingDown   ->GetBinContent(nb) > 0) systJes[0] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[1] = histo_ggWW_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_ggWW->GetBinContent(nb);
    else if(histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVAJESBoundingDown   ->GetBinContent(nb) > 0) systJes[1] = histo_ggWW->GetBinContent(nb)/histo_ggWW_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[2] = histo_Top_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_Top->GetBinContent(nb);
    else if(histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVAJESBoundingDown     ->GetBinContent(nb) > 0) systJes[2] = histo_Top->GetBinContent(nb)/histo_Top_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[3] = histo_DY_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_DY->GetBinContent(nb);
    else if(histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVAJESBoundingDown       ->GetBinContent(nb) > 0) systJes[3] = histo_DY->GetBinContent(nb)/histo_DY_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[4] = histo_VV_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_VV->GetBinContent(nb);
    else if(histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVAJESBoundingDown       ->GetBinContent(nb) > 0) systJes[4] = histo_VV->GetBinContent(nb)/histo_VV_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[5] = histo_VVV_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAJESBoundingDown     ->GetBinContent(nb) > 0) systJes[5] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[6] = histo_WG_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_WG->GetBinContent(nb);
    else if(histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVAJESBoundingDown       ->GetBinContent(nb) > 0) systJes[6] = histo_WG->GetBinContent(nb)/histo_WG_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[7] = histo_WGS_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_WGS->GetBinContent(nb);
    else if(histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVAJESBoundingDown     ->GetBinContent(nb) > 0) systJes[7] = histo_WGS->GetBinContent(nb)/histo_WGS_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAJESBoundingUp   ->GetBinContent(nb) > 0) systJes[8] = histo_Higgs_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAJESBoundingDown ->GetBinContent(nb) > 0) systJes[8] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_MVAJESBoundingDown->GetBinContent(nb);

    double systPU[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_PUBoundingUp	    ->GetBinContent(nb) > 0) systPU[0] = histo_qqWW_CMS_PUBoundingUp->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_PUBoundingDown       ->GetBinContent(nb) > 0) systPU[0] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_PUBoundingUp	    ->GetBinContent(nb) > 0) systPU[1] = histo_ggWW_CMS_PUBoundingUp->GetBinContent(nb)/histo_ggWW->GetBinContent(nb);
    else if(histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_PUBoundingDown       ->GetBinContent(nb) > 0) systPU[1] = histo_ggWW->GetBinContent(nb)/histo_ggWW_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_PUBoundingUp	    ->GetBinContent(nb) > 0) systPU[2] = histo_Top_CMS_PUBoundingUp->GetBinContent(nb)/histo_Top->GetBinContent(nb);
    else if(histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_PUBoundingDown         ->GetBinContent(nb) > 0) systPU[2] = histo_Top->GetBinContent(nb)/histo_Top_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_PUBoundingUp	            ->GetBinContent(nb) > 0) systPU[3] = histo_DY_CMS_PUBoundingUp->GetBinContent(nb)/histo_DY->GetBinContent(nb);
    else if(histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_PUBoundingDown           ->GetBinContent(nb) > 0) systPU[3] = histo_DY->GetBinContent(nb)/histo_DY_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_PUBoundingUp	            ->GetBinContent(nb) > 0) systPU[4] = histo_VV_CMS_PUBoundingUp->GetBinContent(nb)/histo_VV->GetBinContent(nb);
    else if(histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_PUBoundingDown           ->GetBinContent(nb) > 0) systPU[4] = histo_VV->GetBinContent(nb)/histo_VV_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_PUBoundingUp	    ->GetBinContent(nb) > 0) systPU[5] = histo_VVV_CMS_PUBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_PUBoundingDown         ->GetBinContent(nb) > 0) systPU[5] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_PUBoundingUp	            ->GetBinContent(nb) > 0) systPU[6] = histo_WG_CMS_PUBoundingUp->GetBinContent(nb)/histo_WG->GetBinContent(nb);
    else if(histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_PUBoundingDown           ->GetBinContent(nb) > 0) systPU[6] = histo_WG->GetBinContent(nb)/histo_WG_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_PUBoundingUp	    ->GetBinContent(nb) > 0) systPU[7] = histo_WGS_CMS_PUBoundingUp->GetBinContent(nb)/histo_WGS->GetBinContent(nb);
    else if(histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_PUBoundingDown         ->GetBinContent(nb) > 0) systPU[7] = histo_WGS->GetBinContent(nb)/histo_WGS_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_PUBoundingUp       ->GetBinContent(nb) > 0) systPU[8] = histo_Higgs_CMS_PUBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_PUBoundingDown     ->GetBinContent(nb) > 0) systPU[8] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_PUBoundingDown->GetBinContent(nb);
    for(int npu=0; npu<9; npu++) if(systPU[npu] > 1.02) systPU[npu] = 1.02;
    for(int npu=0; npu<9; npu++) if(systPU[npu] < 0.98) systPU[npu] = 0.98;

    double systWWNNLO[2] = {1.0,1.0};
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAWW_qup_nlo	    ->GetBinContent(nb) > 0) systWWNNLO[0] = histo_qqWW_CMS_MVAWW_qup_nlo->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAWW_qdown_nlo      ->GetBinContent(nb) > 0) systWWNNLO[0] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_MVAWW_qdown_nlo->GetBinContent(nb);
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAWW_sup_nlo	    ->GetBinContent(nb) > 0) systWWNNLO[1] = histo_qqWW_CMS_MVAWW_sup_nlo->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAWW_sdown_nlo      ->GetBinContent(nb) > 0) systWWNNLO[1] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_MVAWW_sdown_nlo->GetBinContent(nb);


    char outputLimitsShape[200];
    sprintf(outputLimitsShape,"histo_limits_ww%s_%dj_%s_bin%d.txt",finalStateName,nJetsType,ECMsb.Data(),nb-1);
    ofstream newcardShape;
    newcardShape.open(outputLimitsShape);
    newcardShape << Form("imax 1 number of channels\n");
    newcardShape << Form("jmax * number of background\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");
    newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
    newcardShape << Form("bin ww%2s%4s%d ww%2s%4s%d ww%2s%4s%d ww%2s%4s%d ww%2s%4s%d ww%2s%4s%d ww%2s%4s%d ww%2s%4s%d ww%2s%4s%d ww%2s%4s%d ww%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
    //                            0    1    2   3  4  5   6  7   8     9      10
    newcardShape << Form("process qqWW ggWW Top DY VV VVV WG WGS Higgs WjetsM WjetsE\n");
    newcardShape << Form("process -1 0 1 2 3 4 5 6 7 8 9\n");
    newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f %8.5f %8.5f %8.5f  %8.5f  %8.5f %8.5f %8.5f\n",histo_qqWW->GetBinContent(nb),histo_ggWW->GetBinContent(nb),histo_Top->GetBinContent(nb),histo_DY->GetBinContent(nb),histo_VV->GetBinContent(nb),histo_VVV->GetBinContent(nb),histo_WG->GetBinContent(nb),histo_WGS->GetBinContent(nb),histo_Higgs->GetBinContent(nb),histo_WjetsM->GetBinContent(nb),histo_WjetsE->GetBinContent(nb));
    newcardShape << Form("lumi_%4s                               lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE);		     
    newcardShape << Form("%s                                     lnN  %7.5f %7.5f   -   %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",effMName,systLepEffM[0],systLepEffM[1],systLepEffM[3],systLepEffM[4],systLepEffM[5],systLepEffM[6],systLepEffM[7],systLepEffM[8]);
    newcardShape << Form("%s                                     lnN  %7.5f %7.5f   -   %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",effEName,systLepEffE[0],systLepEffE[1],systLepEffE[3],systLepEffE[4],systLepEffE[5],systLepEffE[6],systLepEffE[7],systLepEffE[8]);
    newcardShape << Form("%s                                     lnN  %7.5f %7.5f   -   %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",momMName,systLepResM[0],systLepResM[1],systLepResM[3],systLepResM[4],systLepResM[5],systLepResM[6],systLepResM[7],systLepResM[8]);
    newcardShape << Form("%s                                     lnN  %7.5f %7.5f   -   %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",momEName,systLepResE[0],systLepResE[1],systLepResE[3],systLepResE[4],systLepResE[5],systLepResE[6],systLepResE[7],systLepResE[8]);
    newcardShape << Form("CMS_pu                                 lnN  %7.5f %7.5f   -   %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",systPU[0],systPU[1],systPU[3],systPU[4],systPU[5],systPU[6],systPU[7],systPU[8]);
    newcardShape << Form("CMS_scale_met                          lnN  %7.5f %7.5f   -   %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",systMet[0],systMet[1],systMet[3],systMet[4],systMet[5],systMet[6],systMet[7],systMet[8]);
    newcardShape << Form("CMS_scale_j                            lnN  %7.5f %7.5f   -   %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",systJes[0],systJes[1],systJes[3],systJes[4],systJes[5],systJes[6],systJes[7],systJes[8]);
    newcardShape << Form("CMS_eff_b                              lnN  %7.5f %7.5f   -   %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",syst_btag,syst_btag,syst_btag,syst_btag,syst_btag,syst_btag,syst_btag,syst_btag);
    newcardShape << Form("pdf_qqbar                              lnN  %7.5f  -      -   %7.5f %7.5f %7.5f %7.5f %7.5f   -    -    -  \n",systPDF[0],systPDF[3],systPDF[4],systPDF[5],systPDF[6],systPDF[7]);
    newcardShape << Form("pdf_gg                                 lnN    -   %7.5f   -     -	-     -     -	  -   %7.5f  -    -  \n",1.01,systPDF[8]);
    newcardShape << Form("QCDscale_VVV		                 lnN    -     -     -     -     -   %7.5f   -	  -	-    -    -  \n",systQCDScale[5]);	    
    newcardShape << Form("QCDscale_ggVV		                 lnN    -   %7.5f   -     -     -     -     -     -     -    -    -  \n",1.15);	    
    newcardShape << Form("QCDscale_qqVV		                 lnN  %7.5f   -     -     -   %7.5f   -   %7.5f	%7.5f	-    -    -  \n",systQCDScale[0],systQCDScale[4],systQCDScale[6],systQCDScale[7]);	    
    newcardShape << Form("QCDscale_ggH		                 lnN    -     -     -     -     -     -     -     -   %7.5f  -    -  \n",systQCDScale[8]);	    
    newcardShape << Form("norm_Top_%dj		                 lnN    -     -   %7.5f   -     -     -     -     -     -    -    -  \n",nJetsType,systQCDScale[2]);	    
    newcardShape << Form("norm_DY_%dj		                 lnN    -     -     -   %7.5f   -     -     -     -     -    -    -  \n",nJetsType,systQCDScale[3]);	    
    newcardShape << Form("norm_WjetsM		                 lnN    -     -     -     -     -     -     -     -     -   %7.5f -  \n",1.30);	    
    newcardShape << Form("norm_WjetsE		                 lnN    -     -     -     -     -     -     -     -     -    -  %7.5f\n",1.30);	    
    newcardShape << Form("WWNNLO_resum		                 lnN  %7.5f   -     -     -     -     -     -     -	-    -    -  \n",systWWNNLO[0]);	    
    newcardShape << Form("WWNNLO_scale		                 lnN  %7.5f   -     -     -     -     -     -     -	-    -    -  \n",systWWNNLO[1]);	    
    newcardShape << Form("UEPS                                   lnN  %7.5f %7.5f   -     -     -     -     -     -     -    -    -  \n",1.03,1.03);	    
    if(histo_qqWW->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAqqWWStatBounding_%s_Bin%d    lnN  %7.5f   -     -     -     -     -     -     -     -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_qqWW->GetBinError(nb)/histo_qqWW->GetBinContent(nb));
    if(histo_ggWW->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAggWWStatBounding_%s_Bin%d    lnN    -   %7.5f   -     -     -     -     -     -	 -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_ggWW->GetBinError(nb)/histo_ggWW->GetBinContent(nb));
    if(histo_Top->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVATopStatBounding_%s_Bin%d     lnN    -     -   %7.5f   -     -     -     -     -	 -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_Top->GetBinError(nb)/histo_Top->GetBinContent(nb));
    if(histo_DY->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVADYStatBounding_%s_Bin%d      lnN    -     -     -   %7.5f   -     -     -     -	 -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_DY->GetBinError(nb)/histo_DY->GetBinContent(nb));
    if(histo_VV->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAVVStatBounding_%s_Bin%d      lnN    -     -     -     -   %7.5f   -     -     -	 -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_VV->GetBinError(nb)/histo_VV->GetBinContent(nb));
    if(histo_VVV->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAVVVStatBounding_%s_Bin%d     lnN    -     -     -     -     -   %7.5f   -     -     -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_VVV->GetBinError(nb)/histo_VVV->GetBinContent(nb));
    if(histo_WG->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAWGStatBounding_%s_Bin%d      lnN    -     -     -     -     -     -   %7.5f   -     -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_WG->GetBinError(nb)/histo_WG->GetBinContent(nb));
    if(histo_WGS->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAWGSStatBounding_%s_Bin%d     lnN    -     -     -     -     -     -     -   %7.5f   -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_WGS->GetBinError(nb)/histo_WGS->GetBinContent(nb));
    if(histo_Higgs->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAHiggsStatBounding_%s_Bin%d   lnN    -     -     -     -     -     -     -     -   %7.5f   -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_Higgs->GetBinError(nb)/histo_Higgs->GetBinContent(nb));
    if(histo_WjetsM ->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAWjetsMStatBounding_%s_Bin%d  lnN    -     -     -     -     -     -     -     -	 -   %7.5f -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_WjetsM ->GetBinError(nb)/histo_WjetsM ->GetBinContent(nb));
    if(histo_WjetsE ->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAWjetsEStatBounding_%s_Bin%d  lnN    -     -     -     -     -     -     -     -	 -     - %7.5f\n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+histo_WjetsE ->GetBinError(nb)/histo_WjetsE ->GetBinContent(nb));
    newcardShape.close();
  }
}
