#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

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

bool usePureMC = false; 
double mcPrescale = 5.0;
const bool useDYMVA = false;
const bool doTriggerStudy = true;
const TString typeLepSel = "default";

void baseAnalysis(
 Int_t nsel = 4,
 Int_t typeSel = 4,
 Int_t period = 1
 ){

  TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_76x/met_";
  TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_76x/met_";
  Double_t lumi = 2.318;

  if(nsel == 2) usePureMC = true;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  if(period==1){
  puPath = "MitAnalysisRunII/data/76x/puWeights_76x.root";
  infilenamev.push_back(Form("%sdata_AOD_Run2015C_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);

  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                                            infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(1);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(2);

  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));						   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(3);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(3);

  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  				   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3+AODSIM.root",filesPathMC.Data())); 			   infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                         infilecatv.push_back(4);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(5);

  if(nsel!=11){
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(6);
  } else {
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                   infilecatv.push_back(6);
  infilenamev.push_back(Form("%sWGstarToLNuMuMu_012Jets_13TeV-madgraph+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 		           infilecatv.push_back(6);
  infilenamev.push_back(Form("%sWGstarToLNuEE_012Jets_13TeV-madgraph+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 		                   infilecatv.push_back(6);
  }

  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(7);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(7);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(7);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                            infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 		   infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(7);

  }
  else {assert(0);}
  
  //infilenamev.clear();infilecatv.clear();
  //infilenamev.push_back(Form("/scratch5/ceballos/ntuples_weights_74x/WWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root")); infilecatv.push_back(1);

  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
 
  float dymva_= -999.;
  unsigned int nlep_= -1;
  unsigned int njets_= -1;

  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
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

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 67;
  const int histBins = 8;
  TH1D* histo[allPlots][histBins];

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  5) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot ==  6 || thePlot == 11) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot =  40.0;}
    else if(thePlot >=  8 && thePlot <=  8) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot >=  9 && thePlot <=  9) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >= 10 && thePlot <= 10) {nBinPlot =  50; xminPlot =-0.5; xmaxPlot =  49.5;}
    else if(thePlot >= 12 && thePlot <= 14) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 15 && thePlot <= 18) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >= 20 && thePlot <= 20) {nBinPlot = 200; xminPlot =-1.0; xmaxPlot =   1.0;}
    else if(thePlot >= 21 && thePlot <= 23) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 24 && thePlot <= 24) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   5.0;}
    else if(thePlot >= 25 && thePlot <= 27) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 28 && thePlot <= 28) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 29 && thePlot <= 31) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 32 && thePlot <= 32) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 33 && thePlot <= 35) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 36 && thePlot <= 36) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 37 && thePlot <= 40) {nBinPlot = 400; xminPlot = -50; xmaxPlot =  50.0;}
    else if(thePlot >= 41 && thePlot <= 42) {nBinPlot = 200; xminPlot = -TMath::Pi(); xmaxPlot = TMath::Pi();}
    else if(thePlot >= 43 && thePlot <= 46) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 47 && thePlot <= 50) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >= 51 && thePlot <= 52) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >= 53 && thePlot <= 53) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot =   4.0;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Clear();
  }
  TH1D* histo_rinout_met[2][3][2];
  TH1D* histo_rinout_dym[2][3][2];
  for(int i=0; i<2; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<2; k++){
        histo_rinout_met[i][j][k] = new TH1D(Form("histo_rinout_met_%d_%d_%d",i,j,k), Form("histo_rinout_met_%d_%d_%d",i,j,k), 4, -0.5, 3.5);
        histo_rinout_dym[i][j][k] = new TH1D(Form("histo_rinout_dym_%d_%d_%d",i,j,k), Form("histo_rinout_dym_%d_%d_%d",i,j,k), 4, -0.5, 3.5);
      }
    }
  }

  unsigned int numberOfLeptons = 2;
  double ptl2nd = 20;
  if     (nsel == 2)  {numberOfLeptons = 2; ptl2nd = 10;}
  else if(nsel == 4)  {numberOfLeptons = 3; ptl2nd = 10;}
  else if(nsel == 7)  {numberOfLeptons = 4; ptl2nd = 10;}
  else if(nsel == 11) {numberOfLeptons = 3; ptl2nd = 20;}

  double totalEventsProcess[50];
  std::vector<double> sumEventsProcess(infilenamev.size(), 0.0);

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");

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

    Int_t nPassTrigger[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    Int_t nPassCuts[10] = {0,0,0,0,0,0,0,0,0,0};
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
    //for (int i=0; i<10000; ++i) {
      the_input_tree->GetEntry(i);

      if(i%1000000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      Bool_t passFilter[10] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
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
        if(doTriggerStudy){
          if(strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*")  == 0) nPassTrigger[ 0]++;
          if(strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*")  == 0) nPassTrigger[ 1]++;
          if(strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) nPassTrigger[ 2]++;
          if(strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) nPassTrigger[ 3]++;
          if(strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*")  	     == 0) nPassTrigger[ 4]++;
          if(strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")	     == 0) nPassTrigger[ 5]++;
          if(strcmp(tokens[nt],"HLT_IsoMu27_v*")  				     == 0) nPassTrigger[ 6]++;
          if(strcmp(tokens[nt],"HLT_IsoMu20_v*")  				     == 0) nPassTrigger[ 7]++;
          if(strcmp(tokens[nt],"HLT_IsoTkMu20_v*")				     == 0) nPassTrigger[ 8]++;
          if(strcmp(tokens[nt],"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")	     == 0) nPassTrigger[ 9]++;
          if(strcmp(tokens[nt],"HLT_Ele23_WPLoose_Gsf_v*")			     == 0) nPassTrigger[10]++;
          if(strcmp(tokens[nt],"HLT_Ele22_eta2p1_WP75_Gsf_v*")			     == 0) nPassTrigger[11]++;
          if(strcmp(tokens[nt],"HLT_Ele27_WPLoose_Gsf_v*")			     == 0) nPassTrigger[12]++;
          if(strcmp(tokens[nt],"HLT_Ele27_WP85_Gsf_v*")				     == 0) nPassTrigger[13]++;
        }
      }
      //if(infilecatv[ifile] != 0) passFilter[1] = kTRUE; // do not apply trigger filters to MC
      if(passFilter[0] == kTRUE) nPassCuts[0]++;
      if(passFilter[0] == kFALSE) continue;
      if(passFilter[1] == kTRUE) nPassCuts[1]++;
      if(passFilter[1] == kFALSE) continue;
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}
      }
      if(idLep.size()!=idTight.size()) {assert(1); return;}
      if(idLep.size()==numberOfLeptons) passFilter[2] = kTRUE;
      if(passFilter[2] == kTRUE) nPassCuts[2]++;
      if(passFilter[2] == kFALSE) continue;

      if(goodIsTight == idTight.size()) passFilter[3] = kTRUE;
      if(usePureMC == false || passFilter[3] ==  kTRUE) nPassCuts[3]++;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= ptl2nd) continue;

      if     (typeSel == 0 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) passFilter[4] = kTRUE;
      else if(typeSel == 3 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) passFilter[4] = kTRUE;
      else if(typeSel == 5 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]))         passFilter[4] = kTRUE;
      else if(typeSel == 6 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])!=TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]))         passFilter[4] = kTRUE;
      else if(typeSel == 4)                                                                                                               passFilter[4] = kTRUE;

      double dPhiLepMETMin = 999.;
      int signQ = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        signQ = signQ + (int)(*eventLeptons.pdgId)[idLep[nl]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]);
        if(dPhiLepMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0]))))
           dPhiLepMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));      
      }
      double minMET  = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      double minPMET = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      if(dPhiLepMETMin < TMath::Pi()/2) minPMET = minPMET * sin(dPhiLepMETMin);

      if     (numberOfLeptons == 2 && nsel == 6)            passFilter[4] = passFilter[4] * (signQ != 0);
      else if(numberOfLeptons == 2 || numberOfLeptons == 4) passFilter[4] = passFilter[4] * (signQ == 0);
      else if(numberOfLeptons == 3)                         passFilter[4] = passFilter[4] * (TMath::Abs(signQ) == 1);
      else {printf("PROBLEM numberOfLeptons\n");assert(0);return;}

      if(passFilter[4] == kTRUE) nPassCuts[4]++;
      if(passFilter[4] == kFALSE) continue;

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 
      double minMassll = 999.0;
      double minMassZ = 999.0;
      double mass3l = 0.0;
      double deltaRllMin = 999.0;
      double minMassSF = 999.0;
      double deltaRllMinWGS = 999.0;
      int type3l = 0;int type3lWGS = 0;
      int tagZ[3] = {-1,-1,0};int tagWGS[3] = {-1,-1,0};
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
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl1]]) &&
	     dilepAux.M() < minMassSF) {
	     minMassSF = dilepAux.M();tagWGS[0]=nl0;tagWGS[1]=nl1;
	     if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]]) == 13) type3lWGS = 0;
	     else                                                         type3lWGS = 1;
	  }

	  if(minMassll > dilepAux.M()) minMassll = dilepAux.M();
        }
      }
      if(minMassll > 12 || (numberOfLeptons == 4 && minMassll > 4) || numberOfLeptons == 3) passFilter[5] = kTRUE;

      vector<int> idJet;
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

        if(dPhiJetMET   == -1) dPhiJetMET   = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])))*180./TMath::Pi();
        if(dPhiJetDiLep == -1) dPhiJetDiLep = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventJets.p4)[nj])))*180./TMath::Pi();

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 10 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) continue;
        
	idJet.push_back(nj);
      }
      double dPhiDiLepMET = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])))*180./TMath::Pi();
      double ptFrac = TMath::Abs(dilep.Pt()-((TLorentzVector*)(*eventMet.p4)[0])->Pt())/dilep.Pt();

      if(useDYMVA == true){
        if(nlep_ != idLep.size()) {printf("PROBLEM nlep %d != %d\n",(int)nlep_,(int)idLep.size()); assert(1); return;}
        if(njets_ != idJet.size()) {printf("PROBLEM njet %d != %d\n",(int)njets_,(int)idJet.size()); assert(1); return;}
      }

      if(nsel == 0 || nsel == 6){ // WW selection (opposite-sign / same-sign)
        if(idJet.size() == 0) passFilter[6] = kTRUE;
	passFilter[6] = kTRUE;
	if(bDiscrMax < 0.605 && idSoft.size() == 0) 
	passFilter[7] = kTRUE;
	if(minPMET > 20) passFilter[8] = kTRUE;
	if(dilep.Pt() > 30) 
	passFilter[9] = kTRUE;
	if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])){
	  if(TMath::Abs(dilep.M()-91.1876) <= 15.0) passFilter[5] = kFALSE;
	  if(minPMET <= 45) passFilter[8] = kFALSE;
	  if(dilep.Pt() <= 45) passFilter[9] = kFALSE;
	}
        //passFilter[5] = kTRUE;
        //passFilter[6] = kTRUE;
        //passFilter[7] = kTRUE;
        //passFilter[8] = kTRUE;
        //passFilter[9] = kTRUE;
      }
      else if(nsel == 1 || nsel == 9){ // ttbar selection
        if     (idJet.size() == 0 && nsel == 1) passFilter[6] = kTRUE;
        else if(idJet.size() == 1 && nsel == 9) passFilter[6] = kTRUE;
	if(!(bDiscrMax < 0.605 && idSoft.size() == 0)) passFilter[7] = kTRUE;
	if(minPMET > 20) passFilter[8] = kTRUE;
	if(dilep.Pt() > 30) passFilter[9] = kTRUE;
	if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])){
	  if(TMath::Abs(dilep.M()-91.1876) <= 15.0) passFilter[5] = kFALSE;
	  if(minPMET <= 45) passFilter[8] = kFALSE;
	  if(dilep.Pt() <= 45) passFilter[9] = kFALSE;
	}
      }
      else if(nsel == 2 || nsel == 10){ // Z selection in mZ mass window, anti b-tagging against top events for e-mu
        passFilter[6] = kTRUE;
	if((bDiscrMax < 0.605 && idSoft.size() == 0) || TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])) passFilter[7] = kTRUE;
	passFilter[8] = kTRUE;
	passFilter[9] = kTRUE;
	if   (TMath::Abs(dilep.M()-91.1876)>15.0 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])) passFilter[5] = kFALSE;
	else                                                                                                                                            passFilter[5] = kTRUE;              
      }
      else if(nsel == 3) { // Z selection
        passFilter[6] = kTRUE;
        passFilter[7] = kTRUE;
        passFilter[8] = kTRUE;
        passFilter[9] = kTRUE;
      }
      else if(nsel == 4) { // WZ selection
        for(unsigned nl0=0; nl0<idLep.size(); nl0++){
	  if((int)nl0==tagZ[0]||(int)nl0==tagZ[1]) continue;
	  tagZ[2] = nl0;
	  break;
	}
	if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]) == 13) type3l += 0;
	else							         type3l += 2;
        TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + 
	                          ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) + 
	                          ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[2])) ) ));
        mass3l = dilepAux.M();
	passFilter[6] = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 30;
        passFilter[7] = TMath::Abs(minMassZ-91.1876) < 30.0;
        passFilter[8] = ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt() > 20;
	passFilter[9] = mass3l > 100 && deltaRllMin > 0.1;
	if(passFilter[7]==kTRUE && (tagZ[0] == tagZ[1] || tagZ[0] == tagZ[2] || tagZ[1] == tagZ[2])) {printf("ZPROBLEM!\n");assert(0);return;}
      }
      else if(nsel == 5){ // Z(ll)H(inv)
	if   (TMath::Abs(dilep.M()-91.1876)>15.0) passFilter[5] = kFALSE;
	else                                      passFilter[5] = kTRUE;              
        passFilter[6] = idJet.size() == 0;
	passFilter[7] = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 100.;
	passFilter[8] = ptFrac < 0.5 && dPhiDiLepMET > 150.;
	passFilter[9] = bDiscrMax < 0.605 && idSoft.size() == 0;
      }
      else if(nsel == 7) { // ZZ4l selection
        TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + 
	                          ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) + 
	                          ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[2])) ) + 
	                          ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[3])) ) ));
        mass3l = dilepAux.M();
	int typeL[2] = {0,0};
        for(unsigned nl0=0; nl0<idLep.size(); nl0++){
          if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]]) == 11) typeL[0]++;
          else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]]) == 13) typeL[1]++;
          else {printf("ZZ4lPROBLEM!\n");assert(0);return;}
        }
	passFilter[6] = typeL[0] == 4 || typeL[1] == 4 || (typeL[0] == 2 && typeL[1] == 2);
        passFilter[7] = kTRUE;
        passFilter[8] = kTRUE;
        passFilter[9] = kTRUE;
	if     (typeL[0] == 4) type3l = 0;
	else if(typeL[1] == 4) type3l = 1;
	else if(typeL[0] == 2) type3l = 2;
	else                   type3l = 3;
      }
      else if(nsel == 8){ // WW same-flavor like selection for all final states with loose MET
        passFilter[6] = dilep.M() > 20;
	if(bDiscrMax < 0.605 && idSoft.size() == 0) passFilter[7] = kTRUE;
	if(minPMET > 20) passFilter[8] = kTRUE;
	if(dilep.Pt() > 45) passFilter[9] = kTRUE;
      }
      else if(nsel == 11){ // wg* selection
        for(unsigned nl0=0; nl0<idLep.size(); nl0++){
	  if((int)nl0==tagWGS[0]||(int)nl0==tagWGS[1]) continue;
	  tagWGS[2] = nl0;
	  break;
	}
	if(tagWGS[0]!=-1){
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagWGS[2]]]) == 13) type3lWGS += 0;
	  else							           type3lWGS += 2;
          TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + 
	                            ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) + 
	                            ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[2])) ) ));
          mass3l = dilepAux.M();
	  deltaRllMinWGS = ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagWGS[0]]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagWGS[2]]]));
	  if(deltaRllMinWGS > ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagWGS[1]]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagWGS[2]]]))) 
	     deltaRllMinWGS = ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagWGS[1]]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagWGS[2]]]));
	  passFilter[6] = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 30;
          passFilter[7] = bDiscrMax < 0.935;
          passFilter[8] = ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagWGS[2]]])->Pt() > 20;
          passFilter[9] = minMassSF < 200.0;
	 if(tagWGS[0] == tagWGS[1] || tagWGS[0] == tagWGS[2] || tagWGS[1] == tagWGS[2]) {printf("WGSPROBLEM!\n");assert(0);return;}
	}
      }
      else {assert(0); return;}
      //printf("5 %d %f\n",(int)eventEvent.eventNum,minMassll);
      if(passFilter[5] == kTRUE) nPassCuts[5]++;
      if(passFilter[5] == kFALSE) continue;
      //printf("6 %d %f\n",(int)eventEvent.eventNum,((TLorentzVector*)(*eventMet.p4)[0])->Pt());
      if(passFilter[6] == kTRUE) nPassCuts[6]++;
      if(passFilter[6] == kFALSE) continue;
      //printf("7 %d %f\n",(int)eventEvent.eventNum,TMath::Abs(minMassZ-91.1876));
      if(passFilter[7] == kTRUE) nPassCuts[7]++;
      if(passFilter[7] == kFALSE) continue;
      //printf("8 %d %f\n",(int)eventEvent.eventNum,((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt());
      if(passFilter[8] == kTRUE) nPassCuts[8]++;
      if(passFilter[8] == kFALSE) continue;
      //printf("9 %d %f %f\n",(int)eventEvent.eventNum,mass3l,deltaRllMin);
      if(passFilter[9] == kTRUE) nPassCuts[9]++;
      if(passFilter[9] == kFALSE) continue;
      double deltaPhiDileptonMet = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtW = TMath::Sqrt(2.0*dilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiDileptonMet)));

      double deltaPhiLeptonMet = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtLN = TMath::Sqrt(2.0*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagWGS[2]]])->Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiLeptonMet)));

      vector<int> idGenJet;
      //printf("----------------------------\n");
      for(int ngenj=0; ngenj<eventMonteCarlo.jetP4->GetEntriesFast(); ngenj++) {
        //printf("ngenj(%d) %f %f %f\n",ngenj,((TLorentzVector*)(*eventMonteCarlo.jetP4)[ngenj])->Pt(),
	//                                    ((TLorentzVector*)(*eventMonteCarlo.jetP4)[ngenj])->Eta(),
	//				    ((TLorentzVector*)(*eventMonteCarlo.jetP4)[ngenj])->Phi());
        Bool_t isGenLepton = kFALSE;
        for(int ngenl=0; ngenl<eventMonteCarlo.p4->GetEntriesFast(); ngenl++) {
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngenl]) != 12 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngenl]) != 14 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngenl]) != 16 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngenl]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngenl]) != 13) continue;
          if(((TLorentzVector*)(*eventMonteCarlo.jetP4)[ngenj])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngenl])) < 0.3) isGenLepton = kTRUE;
        }
        if(isGenLepton == kTRUE) continue;
        if(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[ngenj])->Eta()) >= 5) continue;
        idGenJet.push_back(ngenj);
        //printf("iso genjet\n");
      }
      //for(int ngenl=0; ngenl<eventMonteCarlo.p4->GetEntriesFast(); ngenl++) {
      //  printf("ngenl(%d) %f %f %f %d\n",ngenl,((TLorentzVector*)(*eventMonteCarlo.p4)[ngenl])->Pt(),
     //	                                    ((TLorentzVector*)(*eventMonteCarlo.p4)[ngenl])->Eta(),
	//				    ((TLorentzVector*)(*eventMonteCarlo.p4)[ngenl])->Phi(),
	//				    (int)(*eventMonteCarlo.pdgId)[ngenl]);
      //}
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
      double puWeight = 1.0; if(infilecatv[ifile] != 0) puWeight = nPUScaleFactor(fhDPU, (double)eventMonteCarlo.puTrueInt);
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
      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false){
        if     (theCategory == 5){ // remove W+jets from MC
          fakeSF = 0.0;
        }
        else if(theCategory == 2 && goodIsTight != idTight.size()){ // remove Z+jets from MC as fakeable objects
          fakeSF = 0.0;
        }
        else if((infilecatv[ifile] == 0 || infilecatv[ifile] == 6 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add W+jets from data
          for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    effSF = effSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    theCategory = 5;
          }
          if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) effSF =  1.0 * effSF; // double fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) effSF = -1.0 * effSF; // single fake, MC
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) effSF = -1.0 * effSF; // double fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) effSF =  1.0 * effSF; // single fake, data
        }
        else if(infilecatv[ifile] != 0 && infilecatv[ifile] != 6 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
	  fakeSF = 1.0;
        }
        else if(infilecatv[ifile] == 0 || infilecatv[ifile] == 6){ // data or W+gamma with all good leptons
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
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale*trigEff;
      //printf("AAA %f %f %f %f %f %f %f\n",eventMonteCarlo.mcWeight,theLumi,puWeight,effSF,fakeSF,theMCPrescale,trigEff);
      if(totalWeight == 0) continue;
      // end event weighting
      if(infilecatv[ifile] != 0 || theCategory == 0) sumEventsProcess[ifile] += totalWeight;

      TVector2 dilv(dilep.Px(), dilep.Py());

      TVector2 metv(((TLorentzVector*)(*eventMet.p4)[0])->Px(), ((TLorentzVector*)(*eventMet.p4)[0])->Py());
      TVector2 utv_met = -1.*(metv+dilv);
      double phiv_met = utv_met.DeltaPhi(dilv);
      double the_upara_met = utv_met.Mod()*TMath::Cos(phiv_met);
      double the_uperp_met = utv_met.Mod()*TMath::Sin(phiv_met);

      TVector2 trkmetv((double)eventMet.trackMet->Px(), (double)eventMet.trackMet->Py());
      TVector2 utv_trkmet = -1.*(metv+dilv);
      double phiv_trkmet = utv_trkmet.DeltaPhi(dilv);
      double the_upara_trkmet = utv_trkmet.Mod()*TMath::Cos(phiv_trkmet);
      double the_uperp_trkmet = utv_trkmet.Mod()*TMath::Sin(phiv_trkmet);
      if(nsel == 8){
	double theRinoutBin[2] = {-1, -1}; int nJetsBin = -1;
	if     (minPMET > 20 && minPMET <= 25    ) theRinoutBin[0] = 0;
	else if(minPMET > 25 && minPMET <= 35    ) theRinoutBin[0] = 1;
	else if(minPMET > 35 && minPMET <= 45    ) theRinoutBin[0] = 2;
	else if(minPMET > 45		         ) theRinoutBin[0] = 3;
	if     (dymva_ > -0.20 && dymva_ <=  0.10 && idJet.size() == 0) theRinoutBin[1] = 0;
	else if(dymva_ >  0.10 && dymva_ <=  0.20 && idJet.size() == 0) theRinoutBin[1] = 1;
	else if(dymva_ >  0.20 && dymva_ <=  0.30 && idJet.size() == 0) theRinoutBin[1] = 2;
	else if(dymva_ >  0.30		          && idJet.size() == 0) theRinoutBin[1] = 3;
	else if(dymva_ > -0.20 && dymva_ <= -0.05 && idJet.size() == 1) theRinoutBin[1] = 0;
	else if(dymva_ > -0.05 && dymva_ <=  0.05 && idJet.size() == 1) theRinoutBin[1] = 1;
	else if(dymva_ >  0.05 && dymva_ <=  0.20 && idJet.size() == 1) theRinoutBin[1] = 2;
	else if(dymva_ >  0.20		          && idJet.size() == 1) theRinoutBin[1] = 3;
	else if(dymva_ > -0.20 && dymva_ <= -0.05 && idJet.size() == 2) theRinoutBin[1] = 0;
	else if(dymva_ > -0.05 && dymva_ <=  0.05 && idJet.size() == 2) theRinoutBin[1] = 1;
	else if(dymva_ >  0.05 && dymva_ <=  0.15 && idJet.size() == 2) theRinoutBin[1] = 2;
	else if(dymva_ >  0.15		          && idJet.size() == 2) theRinoutBin[1] = 3;
	if     (idJet.size() == 0) nJetsBin = 0;
	else if(idJet.size() == 1) nJetsBin = 1;
	else if(idJet.size() == 2) nJetsBin = 2;
	int theInOut = 1; if(TMath::Abs(dilep.M()-91.1876) <= 15.0) theInOut = 0;
	double leptonPair = 1;if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])!=TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])) leptonPair = -1;
	if(theCategory == 4) leptonPair = -1 * leptonPair;
	int plotType = -1;
	if     (theCategory == 0 || theCategory == 4) plotType = 0;
	else if(theCategory == 2) plotType = 1;
        // no emu substraction for Z MC
	if(theCategory == 2 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])!=TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])) nJetsBin = -1;
	
	if(theRinoutBin[0] >= 0 && plotType >= 0 && nJetsBin >= 0) histo_rinout_met[plotType][nJetsBin][theInOut]->Fill(theRinoutBin[0],totalWeight*leptonPair);
	if(theRinoutBin[1] >= 0 && plotType >= 0 && nJetsBin >= 0) histo_rinout_dym[plotType][nJetsBin][theInOut]->Fill(theRinoutBin[1],totalWeight*leptonPair);
      }

      for(int thePlot=0; thePlot<allPlots; thePlot++){
	double theVar = 0.0;
        bool makePlot = false;
	if     (thePlot ==  0) {makePlot = true;theVar = TMath::Min(minMassll,199.999);}
	else if(thePlot ==  1) {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot ==  2) {makePlot = true;theVar = TMath::Min((double)eventMet.trackMet->Pt(),199.999);}
	else if(thePlot ==  3) {makePlot = true;theVar = TMath::Min(minPMET,199.999);}
	else if(thePlot ==  4) {makePlot = true;theVar = TMath::Min(dilep.Pt(),199.999);}
	else if(thePlot ==  5) {makePlot = true;theVar = TMath::Min(mtW,199.999);}
	else if(thePlot ==  6) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot ==  7) {makePlot = true;theVar = TMath::Min((double)eventEvent.rho,39.999);}
	else if(thePlot ==  8) {makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	else if(thePlot ==  9) {makePlot = true;theVar = TMath::Min(bDiscrMax,0.999);}
	else if(thePlot == 10) {makePlot = true;theVar = TMath::Min((double)(numberQuarks[0]+10*numberQuarks[1]),49.499);}
	else if(thePlot == 11 && bDiscrMax < 0.97 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 12) {makePlot = true;theVar = TMath::Min(mass3l,199.999);}
	else if(thePlot == 13) {makePlot = true;theVar = TMath::Min(mtLN,199.999);}
	else if(thePlot == 14) {makePlot = true;theVar = TMath::Min(minMET,199.999);}
	else if(thePlot == 15) {makePlot = true;theVar = dPhiJetMET;}
	else if(thePlot == 16) {makePlot = true;theVar = dPhiJetDiLep;}
	else if(thePlot == 17) {makePlot = true;theVar = dPhiDiLepMET;}
	else if(thePlot == 18) {makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	else if(thePlot == 19) {makePlot = true;theVar = TMath::Min(ptFrac,0.999);}
	else if(thePlot == 20) {makePlot = true;theVar = TMath::Max(TMath::Min((double)dymva_,0.999),-0.999);}
	else if(thePlot == 21) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),199.999);}
	else if(thePlot == 22) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),199.999);}
	else if(thePlot == 23 && idGenJet.size() >= 1) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[0]])->Pt(),199.999);}
	else if(thePlot == 24 && idGenJet.size() >= 1 && ((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[0]])->Pt() > 30) {makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[0]])->Eta());}
	else if(thePlot == 25 && idJet.size() == 0 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 26 && idJet.size() == 0 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min((double)eventMet.trackMet->Pt(),199.999);}
	else if(thePlot == 27 && idJet.size() == 0 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min(minMET,199.999);}
	else if(thePlot == 28 && idJet.size() == 0 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	else if(thePlot == 29 && idJet.size() == 1 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 30 && idJet.size() == 1 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min((double)eventMet.trackMet->Pt(),199.999);}
	else if(thePlot == 31 && idJet.size() == 1 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min(minMET,199.999);}
	else if(thePlot == 32 && idJet.size() == 1 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	else if(thePlot == 33 && idJet.size() == 2 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 34 && idJet.size() == 2 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min((double)eventMet.trackMet->Pt(),199.999);}
	else if(thePlot == 35 && idJet.size() == 2 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min(minMET,199.999);}
	else if(thePlot == 36 && idJet.size() == 2 && minPMET > 20 && dilep.Pt() > 45 && bDiscrMax < 0.605 && idSoft.size() == 0) {makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	else if(thePlot == 37 && idJet.size() == 0) {makePlot = true;theVar = the_upara_met;}
	else if(thePlot == 38 && idJet.size() == 0) {makePlot = true;theVar = the_uperp_met;}
	else if(thePlot == 39 && idJet.size() == 0) {makePlot = true;theVar = the_upara_trkmet;}
	else if(thePlot == 40 && idJet.size() == 0) {makePlot = true;theVar = the_uperp_trkmet;}
	else if(thePlot == 41) {makePlot = true;theVar = ((TLorentzVector*)(*eventMet.p4)[0])->Phi();}
	else if(thePlot == 42) {makePlot = true;theVar = (double)eventMet.trackMet->Py();}
	else if(thePlot == 43 && numberQuarks[1] == 0                                          ) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 44 && numberQuarks[1] == 0 && bDiscrMax < 0.97 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 45 && numberQuarks[1]  > 0                                          ) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 46 && numberQuarks[1]  > 0 && bDiscrMax < 0.97 && idSoft.size() == 0) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 47 && type3lWGS == 0) {makePlot = true;theVar = TMath::Min(minMassSF,99.999);}
	else if(thePlot == 48 && type3lWGS == 1) {makePlot = true;theVar = TMath::Min(minMassSF,99.999);}
	else if(thePlot == 49 && type3lWGS == 2) {makePlot = true;theVar = TMath::Min(minMassSF,99.999);}
	else if(thePlot == 50 && type3lWGS == 3) {makePlot = true;theVar = TMath::Min(minMassSF,99.999);}
	else if(thePlot == 51) {makePlot = true;theVar = type3l;}
	else if(thePlot == 52) {makePlot = true;theVar = type3lWGS;}
	else if(thePlot == 53) {makePlot = true;theVar = TMath::Min(deltaRllMinWGS,3.999);}
	if(makePlot == true) histo[thePlot][theCategory]->Fill(theVar,totalWeight);
      }
    }
    printf("eff_cuts(%f): ",sumEventsProcess[ifile]);
    for(int nc=0; nc<10; nc++){
      double nminusone = the_input_tree->GetEntries();
      if(nc>0) nminusone = nPassCuts[nc-1];
      printf("(%d): %6.3f(%4.1f) | ",nc,100*(double)nPassCuts[nc]/the_input_all->GetEntries()*theMCPrescale,100*(double)nPassCuts[nc]/nminusone);
    }
    printf("\n");

    if(doTriggerStudy){
      printf("trigger_cuts: ");
      for(unsigned int nc=0; nc<sizeof(nPassTrigger)/sizeof(*nPassTrigger); nc++){
        printf("%4.1f ",100*(double)nPassTrigger[nc]/the_input_tree->GetEntries());
      }
      printf("\n");
    }
  } // end of chain

  if(nsel == 8){
    for(int i=1; i<=histo_rinout_met[0][0][0]->GetNbinsX(); i++){
      printf("MET MC: %7.1f/%7.1f=%.3f %7.1f/%7.1f=%.3f %7.1f/%7.1f=%.3f | DA: %7.1f/%7.1f=%.3f %7.1f/%7.1f=%.3f %7.1f/%7.1f=%.3f\n",
      histo_rinout_met[1][0][1]->GetBinContent(i),histo_rinout_met[1][0][0]->GetBinContent(i),histo_rinout_met[1][0][1]->GetBinContent(i)/histo_rinout_met[1][0][0]->GetBinContent(i),
      histo_rinout_met[1][1][1]->GetBinContent(i),histo_rinout_met[1][1][0]->GetBinContent(i),histo_rinout_met[1][1][1]->GetBinContent(i)/histo_rinout_met[1][1][0]->GetBinContent(i),
      histo_rinout_met[1][2][1]->GetBinContent(i),histo_rinout_met[1][2][0]->GetBinContent(i),histo_rinout_met[1][2][1]->GetBinContent(i)/histo_rinout_met[1][2][0]->GetBinContent(i),
      histo_rinout_met[0][0][1]->GetBinContent(i),histo_rinout_met[0][0][0]->GetBinContent(i),histo_rinout_met[0][0][1]->GetBinContent(i)/histo_rinout_met[0][0][0]->GetBinContent(i),
      histo_rinout_met[0][1][1]->GetBinContent(i),histo_rinout_met[0][1][0]->GetBinContent(i),histo_rinout_met[0][1][1]->GetBinContent(i)/histo_rinout_met[0][1][0]->GetBinContent(i),
      histo_rinout_met[0][2][1]->GetBinContent(i),histo_rinout_met[0][2][0]->GetBinContent(i),histo_rinout_met[0][2][1]->GetBinContent(i)/histo_rinout_met[0][2][0]->GetBinContent(i)
      );
    }
    for(int i=1; i<=histo_rinout_met[0][0][0]->GetNbinsX(); i++){
      printf("dym MC: %7.1f/%7.1f=%.3f %7.1f/%7.1f=%.3f %7.1f/%7.1f=%.3f | DA: %7.1f/%7.1f=%.3f %7.1f/%7.1f=%.3f %7.1f/%7.1f=%.3f\n",
      histo_rinout_dym[1][0][1]->GetBinContent(i),histo_rinout_dym[1][0][0]->GetBinContent(i),histo_rinout_dym[1][0][1]->GetBinContent(i)/histo_rinout_dym[1][0][0]->GetBinContent(i),
      histo_rinout_dym[1][1][1]->GetBinContent(i),histo_rinout_dym[1][1][0]->GetBinContent(i),histo_rinout_dym[1][1][1]->GetBinContent(i)/histo_rinout_dym[1][1][0]->GetBinContent(i),
      histo_rinout_dym[1][2][1]->GetBinContent(i),histo_rinout_dym[1][2][0]->GetBinContent(i),histo_rinout_dym[1][2][1]->GetBinContent(i)/histo_rinout_dym[1][2][0]->GetBinContent(i),
      histo_rinout_dym[0][0][1]->GetBinContent(i),histo_rinout_dym[0][0][0]->GetBinContent(i),histo_rinout_dym[0][0][1]->GetBinContent(i)/histo_rinout_dym[0][0][0]->GetBinContent(i),
      histo_rinout_dym[0][1][1]->GetBinContent(i),histo_rinout_dym[0][1][0]->GetBinContent(i),histo_rinout_dym[0][1][1]->GetBinContent(i)/histo_rinout_dym[0][1][0]->GetBinContent(i),
      histo_rinout_dym[0][2][1]->GetBinContent(i),histo_rinout_dym[0][2][0]->GetBinContent(i),histo_rinout_dym[0][2][1]->GetBinContent(i)/histo_rinout_dym[0][2][0]->GetBinContent(i)
      );
    }
  }


  double sumEvents = 0;
  for(int np=1; np<histBins; np++) sumEvents += histo[0][np]->GetSumOfWeights();
  printf("yields: %f |",histo[0][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[0][np]->GetSumOfWeights());
  printf(" = %.3f\n",sumEvents);
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    char output[200];
    sprintf(output,"histo_nice_%d_%d_%d.root",thePlot,nsel,typeSel);	  
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }
}
