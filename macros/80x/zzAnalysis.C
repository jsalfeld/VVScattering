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

#include "MitAnalysisRunII/macros/80x/factors.h"

int whichSkim = 2;
bool usePureMC = false; 
double mcPrescale = 1.0;
const TString typeLepSel = "medium";

void zzAnalysis(
 ){

  Int_t period = 1;
  TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_80x/met_";
  TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_80x/met_";
  Double_t lumi = 12.9;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";

  puPath = "MitAnalysisRunII/data/80x/puWeights_80x.root";
  infilenamev.push_back(Form("%sdata_Run2016B.root",filesPathDA.Data()));											      infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_Run2016C.root",filesPathDA.Data()));											      infilecatv.push_back(0);

  if(usePureMC == true){
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                                            infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));						   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1+RAWAODSIM.root",filesPathMC.Data()));   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1+RAWAODSIM.root",filesPathMC.Data()));      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1+RAWAODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1+RAWAODSIM.root",filesPathMC.Data()));                infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));  				   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data())); 			   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(1);
  //infilenamev.push_back(Form("%stZq_nunu_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data())); 		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data())); 		           infilecatv.push_back(1);
  }

  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));	 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));	 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(5);

  infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data())); 		   infilecatv.push_back(6);
  infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix+RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3_ext1-v1+RAWAODSIM.root",filesPathMC.Data()));  infilecatv.push_back(6);
  infilenamev.push_back(Form("%sGluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8+RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1+RAWAODSIM.root",filesPathMC.Data()));  infilecatv.push_back(6);
  infilenamev.push_back(Form("%sVBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8+RunIISpring16DR80-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1+RAWAODSIM.root",filesPathMC.Data()));  infilecatv.push_back(6);
  
  //infilenamev.clear();infilecatv.clear();
  //infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	 infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sTT_TuneCUETP8M1_13TeV-powheg-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));     infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sTTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	 infilecatv.push_back(4);

  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
 
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
  TH2D *fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_Tight_Electron"));
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

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 14;
  const int histBins = 7;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {"..Data", "....EM", "Zgamma", "....WZ", "...ZZ", "...VVV", ".Higgs"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >=  1 && thePlot <=  1) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >=  2 && thePlot <=  2) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >=  3 && thePlot <=  4) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  5 && thePlot <=  5) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot =   4.0;}
    else if(thePlot >=  6 && thePlot <=  6) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >=  8 && thePlot <=  9) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 10 && thePlot <= 10) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >= 11 && thePlot <= 13) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Clear();
  }

  unsigned int numberOfLeptons = 4;

  double totalEventsProcess[50];
  std::vector<double> sumEventsProcess(infilenamev.size(), 0.0);

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    //TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");
    TTree *the_SelBit_tree= (TTree*)the_input_file.FindObjectAny("SelBit_tree");

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareJets eventJets;
    eventJets.setBranchAddresses(the_input_tree);

    BareLeptons eventLeptons;
    eventLeptons.setBranchAddresses(the_input_tree);

    //BareTaus eventTaus;
    //eventTaus.setBranchAddresses(the_input_tree);

    BareMet eventMet;
    eventMet.SetExtend();
    eventMet.setBranchAddresses(the_input_tree);

    BareTrigger eventTrigger;
    eventTrigger.setBranchAddresses(the_input_tree);

    //BareVertex eventVertex;
    //eventVertex.setBranchAddresses(the_input_tree);

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
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_SelBit_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      if((selBit_ & 0x1<<whichSkim) == 0) continue;

      the_input_tree->GetEntry(i);

      Bool_t passFilter[10] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
      if(infilecatv[ifile] == 0) {
	for (int nt = 0; nt <(int)numtokens; nt++) {
          if((*eventTrigger.triggerFired)[nt] == 0) continue;
          if((strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
             (strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
             (strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
             (strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*") 	              == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*")	      == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*") 	      == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")	      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoMu20_v*") 				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoTkMu20_v*") 				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoMu22_v*") 				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoTkMu22_v*") 				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoMu24_v*")				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoTkMu24_v*")				      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")	      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")	      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele25_eta2p1_WPTight_Gsf_v*")                    == 0) ||
             (strcmp(tokens[nt],"HLT_Ele27_eta2p1_WPLoose_Gsf_v*")                    == 0) ||
             (strcmp(tokens[nt],"HLT_Ele27_WPTight_Gsf_v*")			      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele35_WPLoose_Gsf_v*")			      == 0)
             ) passFilter[1] = kTRUE;
	}
      } else { passFilter[1] = kTRUE;}

      if(passFilter[0] == kFALSE) continue;
      if(passFilter[1] == kFALSE) continue;
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep); }
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}
      }
      if(idLep.size()!=idTight.size()) {assert(1); return;}
      if(idLep.size()==numberOfLeptons) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;

      if(goodIsTight == idTight.size()) passFilter[3] = kTRUE;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() <= 10 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[3]])->Pt() <= 10) continue;

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

      passFilter[4] = TMath::Abs(signQ) == 0;
      if(passFilter[4] == kFALSE) continue;

      double minMassll = 999.0;
      double minMassZ[2] = {999.0, 999.0};
      double deltaRllMin = 999.0;
      int tagZ[4] = {-1,-1,-1,-1};
      for(unsigned nl0=0; nl0<idLep.size()-1; nl0++){
        for(unsigned nl1=nl0+1; nl1<idLep.size(); nl1++){
          double deltaRllAux = ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl0]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl1]]));
          if(deltaRllAux < deltaRllMin) deltaRllMin = deltaRllAux;

	  if((int)(*eventLeptons.pdgId)[idLep[nl0]] * (int)(*eventLeptons.pdgId)[idLep[nl1]] > 0) continue;
          TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl1])) ) ));

	  if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl1]])){
	    if(TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassZ[0]-91.1876)) {
	      minMassZ[1] = minMassZ[0];tagZ[2]=tagZ[0];tagZ[3]=tagZ[1];
	      minMassZ[0] = dilepAux.M();tagZ[0]=nl0;tagZ[1]=nl1;
	    }
	    else if(TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassZ[1]-91.1876)) {
	      minMassZ[1] = dilepAux.M();tagZ[2]=nl0;tagZ[3]=nl1;
	    }
	  }

	  if(minMassll > dilepAux.M()) minMassll = dilepAux.M();
        }
      }

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

        if(dPhiJetMET   == -1) dPhiJetMET   = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 10 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) continue;
        
	idJet.push_back(nj);
      }

      TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + 
        			( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) + 
        			( *(TLorentzVector*)(eventLeptons.p4->At(idLep[2])) ) + 
        			( *(TLorentzVector*)(eventLeptons.p4->At(idLep[3])) ) ));
      double mass4l = dilepAux.M();
      int typeL[2] = {0,0};
      for(unsigned nl0=0; nl0<idLep.size(); nl0++){
        if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]]) == 11) typeL[0]++;
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]]) == 13) typeL[1]++;
        else {printf("ZZ4lPROBLEM!\n");assert(0);return;}
      }
      passFilter[5] = typeL[0] == 4 || typeL[1] == 4 || (typeL[0] == 2 && typeL[1] == 2);
      passFilter[6] = minMassll > 4.0;
      passFilter[7] = minMassZ[0] > 60 && minMassZ[0] < 120 && minMassZ[1] > 60 && minMassZ[1] < 120;
      
      int type4l = 0;
      if     (typeL[0] == 4) type4l = 0;
      else if(typeL[1] == 4) type4l = 1;
      else if(typeL[0] == 2) type4l = 2;
      else		     type4l = 3;

      bool passNMinusOne[2] = {                 passFilter[6],
                               passFilter[5]                 };
      bool passAllCuts = passFilter[5] && passFilter[6];
      bool passZZSel = passFilter[5] && passFilter[6] && passFilter[7];

      TLorentzVector theFakeMET[2];
      if(passAllCuts){
        theFakeMET[0].SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Px()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Px());
        theFakeMET[0].SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Py()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Py());
        theFakeMET[1].SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Px()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Px());
        theFakeMET[1].SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Py()-((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[3]]])->Py());
      }

      // begin event weighting
      vector<int>zzBoson;
      vector<bool> isGenDupl; double bosonPtMin = 1000000000; bool isBosonFound = false;
      int numberQuarks[2] = {0,0};
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23) zzBoson.push_back(ngen0);
        if((TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23||TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 24) &&
	   ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() < bosonPtMin) {bosonPtMin = ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt(); isBosonFound = true;}
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 4 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[0]++;
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 5 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[1]++;
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

      // fake rate
      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false){
        if     ((infilecatv[ifile] == 0 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add Z+jets from data
	  for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    theCategory = 1;
          }
          if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-3) fakeSF = -1.0 * fakeSF; // triple fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) fakeSF = +1.0 * fakeSF; // double fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) fakeSF = -1.0 * fakeSF; // single fake, MC
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-3) fakeSF = +1.0 * fakeSF; // triple fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) fakeSF = -1.0 * fakeSF; // double fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) fakeSF = +1.0 * fakeSF; // single fake, data
        }
        else if(infilecatv[ifile] != 0 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
	  fakeSF = 1.0;
        }
        else if(infilecatv[ifile] == 0){ // data or Z+gamma with all good leptons
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
      double totalWeight = mcWeight*theLumi*puWeight*fakeSF*theMCPrescale;
      if(totalWeight == 0) continue;

      double theZZCorr[2] {1,1};
      if(theCategory == 4 && infilenamev[ifile].Contains("GluGlu") == kFALSE) {
        float GENmZZ = 0.0;
	if(zzBoson.size() >= 2) GENmZZ = ( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zzBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(zzBoson[1])) ) ).M();
	theZZCorr[0] = weightEWKCorr(bosonPtMin,1);
	theZZCorr[1] = kfactor_qqZZ_qcd_M(GENmZZ);
        totalWeight = totalWeight * (theZZCorr[0]*theZZCorr[1]);
      }

      // end event weighting
      if((infilecatv[ifile] != 0 || theCategory == 0) && passAllCuts) sumEventsProcess[ifile] += totalWeight;

      for(int thePlot=0; thePlot<allPlots; thePlot++){
	double theVar = 0.0;
	bool makePlot = false;
	if     (thePlot ==  0 && passAllCuts)      {makePlot = true;theVar = TMath::Min((double)mass4l,399.999);}
	else if(thePlot ==  1 && passNMinusOne[0]) {makePlot = true;theVar = type4l;}
	else if(thePlot ==  2 && passNMinusOne[1]) {makePlot = true;theVar = TMath::Min(minMassll,99.999);}
	else if(thePlot ==  3 && passAllCuts)      {makePlot = true;theVar = TMath::Min(minMassZ[0],199.999);}
	else if(thePlot ==  4 && passAllCuts)      {makePlot = true;theVar = TMath::Min(minMassZ[1],199.999);}
	else if(thePlot ==  5 && passAllCuts)      {makePlot = true;theVar = TMath::Min(deltaRllMin,3.999);}
	else if(thePlot ==  6 && passAllCuts)      {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot ==  7 && passAllCuts)      {makePlot = true;theVar = TMath::Min(bDiscrMax,0.999);}
	else if(thePlot ==  8 && passAllCuts)      {makePlot = true;theVar = dPhiJetMET*180/TMath::Pi();}
	else if(thePlot ==  9 && passAllCuts)      {makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	else if(thePlot == 10 && passZZSel)        {makePlot = true;theVar = TMath::Min((double)mass4l,399.999);}
	else if(thePlot == 11 && passZZSel)        {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 12 && passZZSel)	   {makePlot = true;theVar = TMath::Min(TMath::Max(theFakeMET[0].Pt(),theFakeMET[1].Pt()),199.999);}
	else if(thePlot == 13 && passZZSel)	   {makePlot = true;theVar = TMath::Min(TMath::Min(theFakeMET[0].Pt(),theFakeMET[1].Pt()),199.999);}

	if(makePlot) histo[thePlot][theCategory]->Fill(theVar,totalWeight);
      }
    }
    printf("eff_cuts: %f\n",sumEventsProcess[ifile]);

  } // end of chain

  double sumEvents = 0;
  for(int np=1; np<histBins; np++) sumEvents += histo[0][np]->GetSumOfWeights();
  printf("yields: %f |",histo[0][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[0][np]->GetSumOfWeights());
  printf(" = %.3f\n",sumEvents);
  double sumEventsType[5] = {0,0,0,0,0};
  double sumEventsTypeE[5] = {0,0,0,0,0};
  printf("                  all                 4e                  4m                  2e2m\n");
  printf("---------------------------------------------------------------------------------------\n");
  for(int np=0; np<histBins; np++) {
    double sumEvents = 0.; double sumEventsE = 0.;
    for(int i=1; i<=3; i++) {sumEvents = sumEvents + histo[1][np]->GetBinContent(i); sumEventsE = sumEventsE + histo[1][np]->GetBinError(i)*histo[1][np]->GetBinError(i);}
    if(np!=0){
    sumEventsType[0] = sumEventsType[0] + sumEvents;                      sumEventsTypeE[0] = sumEventsTypeE[0] + sumEventsE * sumEventsE;
    sumEventsType[1] = sumEventsType[1] + histo[1][np]->GetBinContent(1); sumEventsTypeE[1] = sumEventsTypeE[1] + histo[1][np]->GetBinError(1) * histo[1][np]->GetBinError(1);
    sumEventsType[2] = sumEventsType[2] + histo[1][np]->GetBinContent(2); sumEventsTypeE[2] = sumEventsTypeE[2] + histo[1][np]->GetBinError(2) * histo[1][np]->GetBinError(2);
    sumEventsType[3] = sumEventsType[3] + histo[1][np]->GetBinContent(3); sumEventsTypeE[3] = sumEventsTypeE[3] + histo[1][np]->GetBinError(3) * histo[1][np]->GetBinError(3);
    sumEventsType[4] = sumEventsType[4] + histo[1][np]->GetBinContent(4); sumEventsTypeE[4] = sumEventsTypeE[4] + histo[1][np]->GetBinError(4) * histo[1][np]->GetBinError(4);
    }
    sumEventsE = sqrt(sumEventsE);
    if(sumEvents>0)
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    processName[np].Data(),sumEvents,sumEventsE,
    histo[1][np]->GetBinContent(1),histo[1][np]->GetBinError(1),
    histo[1][np]->GetBinContent(2),histo[1][np]->GetBinError(2),
    histo[1][np]->GetBinContent(3),histo[1][np]->GetBinError(3));
    if(np==0)
    printf("---------------------------------------------------------------------------------------\n");
  }
    printf("---------------------------------------------------------------------------------------\n");
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    "   all",
    sumEventsType[0],sqrt(sumEventsTypeE[0]),
    sumEventsType[1],sqrt(sumEventsTypeE[1]),
    sumEventsType[2],sqrt(sumEventsTypeE[2]),
    sumEventsType[3],sqrt(sumEventsTypeE[3]));
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    char output[200];
    sprintf(output,"histozz_nice_%d.root",thePlot);	  
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }
}
