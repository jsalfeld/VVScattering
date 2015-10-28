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

#include "MitAnalysisRunII/macros/factors.h"

bool usePureMC = true; 
double mcPrescale = 1.0;
enum selType                     {ZSEL=0,    SIGSEL,     WWSEL,      WWLOOSESEL,      BTAGSEL, nSelTypes};
TString selTypeName[nSelTypes]= {"ZSEL",    "SIGSEL",   "WWSEL",    "WWLOOSESEL",    "BTAGSEL"};

void zhAnalysis(
 Int_t period = 1,
 Int_t typeSel = 3
 ){

  TString filesPath  = "/scratch5/ceballos/ntuples_weights/";
  Double_t lumi = 0.0715;
  if(period == 1) lumi = 0.5947;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  if      (period==0){
  puPath = "/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src/MitAnalysisRunII/data/puWeights_13TeV_50ns.root";
  infilenamev.push_back(Form("%sdata_AOD_50ns.root",filesPath.Data()));														  infilecatv.push_back(0);
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1+AODSIM.root",filesPath.Data()));	  infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));  	  infilecatv.push_back(2);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(3);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data())); infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZ_TuneCUETP8M1_13TeV-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));					  infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZ_TuneCUETP8M1_13TeV-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));					  infilecatv.push_back(4);
  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1+AODSIM.root",filesPath.Data()));		  infilecatv.push_back(5);
  assert(0);
  }
  else if(period==1){
  puPath = "/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src/MitAnalysisRunII/data/puWeights_13TeV_25ns.root";
  infilenamev.push_back(Form("%sdata_AOD_25ns.root",filesPath.Data())); 													  infilecatv.push_back(0);

  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));					          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));      infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  		  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_amcatnloFXFX_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  	          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); 		          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(1);
  /////infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); 		  infilecatv.push_back(1);
  /////infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));        infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));         infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(1);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));         infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(2);
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(2);

  infilenamev.push_back(Form("%sWZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));                          infilecatv.push_back(3);

  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				          infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				          infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(5);

  infilenamev.push_back(Form("%sZH_ZToMM_HToInvisible_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(6);
  infilenamev.push_back(Form("%sZH_ZToEE_HToInvisible_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(6);

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
  const int allPlots = 15;
  const int histBins = 7;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {"..Data", "....EM", "...DY", "...WZ", "....ZZ", "...VVV", ".Higgs"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >=  1 && thePlot <=  1) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >=  2 && thePlot <=  2) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >=  3 && thePlot <=  3) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  4 && thePlot <=  4) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   2.0;}
    else if(thePlot >=  5 && thePlot <=  5) {nBinPlot = 180; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >=  6 && thePlot <=  6) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >=  7 && thePlot <=  9) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 10 && thePlot <= 10) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot =  40.0;}
    else if(thePlot >= 11 && thePlot <= 11) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot >= 12 && thePlot <= 14) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Clear();
  }

  double nSigCut[nSelTypes*4],nSigECut[nSelTypes*4];
  double bgdDecay[nSelTypes*4][histBins],weiDecay[nSelTypes*4][histBins];
  for(unsigned int i=0; i<nSelTypes*4; i++) {
    nSigCut[i] = 0.0; nSigECut[i] = 0.0;
    for(int j=0; j<histBins; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }

  unsigned int numberOfLeptons = 2;

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

    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_input_tree->GetEntry(i);

      if(i%1000000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      Bool_t passFilter[5] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
      for (int nt = 0; nt <(int)numtokens; nt++) {
	if((strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*")  == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*")  == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*")             == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")           == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_IsoMu27_v*")                                     == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_IsoMu20_v*")				            == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_IsoTkMu20_v*")				    == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")       == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_Ele27_eta2p1_WP75_Gsf_v*")                       == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_Ele27_eta2p1_WPLoose_Gsf_v*")                    == 0 && (*eventTrigger.triggerFired)[nt] == 1)
           ) passFilter[1] = kTRUE;
      }

      if(passFilter[0] == kFALSE) continue;
      if(passFilter[1] == kFALSE) continue;
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut("medium",TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
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
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 20) continue;

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

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 

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

      int typePair = 0;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) typePair = 1;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) typePair = 2;

      double dPhiDiLepMET = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double ptFrac = TMath::Abs(dilep.Pt()-((TLorentzVector*)(*eventMet.p4)[0])->Pt())/dilep.Pt();
      double deltaPhiDileptonMet = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtW = TMath::Sqrt(2.0*dilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiDileptonMet)));

      bool passZMass     = TMath::Abs(dilep.M()-91.1876) < 15.0;
      bool passNjets     = idJet.size() == 0;
      bool passMET       = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 120.;
      bool passPTFrac    = ptFrac < 0.5;
      bool passDPhiZMET  = dPhiDiLepMET*180./TMath::Pi() > 150.;
      bool passBtagVeto  = bDiscrMax < 0.605 && idSoft.size() == 0;
      bool passPTLL      = dilep.Pt() > 60;
      bool passMETMin    = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 70.;

      bool passZMassLarge = TMath::Abs(dilep.M()-91.1876) < 30.0;
      bool passZMassSB    = (dilep.M() > 110.0 && dilep.M() < 200.0);

      bool passNMinusOne[7] = {             passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL,
                               passZMass &&              passMET && passPTFrac && passDPhiZMET && passBtagVeto && passPTLL,
                               passZMass && passNjets &&            passPTFrac && passDPhiZMET && passBtagVeto && passPTLL,
                               passZMass && passNjets && passMET &&               passDPhiZMET && passBtagVeto && passPTLL,
                               passZMass && passNjets && passMET && passPTFrac                 && passBtagVeto && passPTLL,
                               passZMass && passNjets && passMET && passPTFrac && passDPhiZMET                 && passPTLL,
                               passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && passBtagVeto};
      bool passAllCuts[nSelTypes] = {                   passZMass && passNjets,
                                                        passZMass && passNjets && passMET && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL,
                                     passZMassLarge && !passZMass && passNjets && passMET && passPTFrac && passDPhiZMET &&  passBtagVeto && passPTLL,
                                     passZMassSB    && !passZMass && passNjets && passMETMin                            && !passBtagVeto            ,
                                                        passZMass && passNjets && passMET && passPTFrac && passDPhiZMET && !passBtagVeto && passPTLL};

      // begin event weighting
      vector<bool> isGenDupl;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
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
          effSF = effSF * effScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period);
        }
      }

      // fake rate
      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false){
        printf("NEED TO WORK ON IT IF WE WANT TO USE IT\n");return;
        if     (theCategory == 5){ // remove W+jets from MC
          fakeSF = 0.0;
        }
        else if(theCategory == 2 && goodIsTight != idTight.size()){ // remove Z+jets from MC as fakeable objects
          fakeSF = 0.0;
        }
        else if((infilecatv[ifile] == 0 || infilecatv[ifile] == 6 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add W+jets from data
          for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    effSF = effSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period);
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
      //double totalWeight = eventMonteCarlo.mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale;
      double totalWeight = eventMonteCarlo.mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale;
      if(totalWeight == 0) continue;
      // end event weighting
      if((infilecatv[ifile] != 0 || theCategory == 0) && passAllCuts[SIGSEL]) sumEventsProcess[ifile] += totalWeight;

      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i]) {
          bgdDecay[i+typePair*nSelTypes][theCategory] += totalWeight;
          weiDecay[i+typePair*nSelTypes][theCategory] += totalWeight*totalWeight;
        }
      }

      if((typeSel == typePair) || (typeSel == 3 && (typePair == 1 || typePair == 2))) {
	for(int thePlot=0; thePlot<allPlots; thePlot++){
	  double theVar = 0.0;
	  bool makePlot = false;
	  if     (thePlot ==  0 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min(mtW,499.999);}
	  else if(thePlot ==  1 && passNMinusOne[0])   {makePlot = true;theVar = TMath::Min(TMath::Abs(dilep.M()-91.1876),99.999);}
	  else if(thePlot ==  2 && passNMinusOne[1])   {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	  else if(thePlot ==  3 && passNMinusOne[2])   {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	  else if(thePlot ==  4 && passNMinusOne[3])   {makePlot = true;theVar = TMath::Min(ptFrac,1.999);}
	  else if(thePlot ==  5 && passNMinusOne[4])   {makePlot = true;theVar = dPhiDiLepMET*180/TMath::Pi();}
	  else if(thePlot ==  6 && passNMinusOne[5])   {makePlot = true;theVar = TMath::Max(TMath::Min(bDiscrMax,0.999),0.001);}
	  else if(thePlot ==  7 && passNMinusOne[6])   {makePlot = true;theVar = TMath::Min(dilep.Pt(),199.999);}
	  else if(thePlot ==  8 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),199.999);}
	  else if(thePlot ==  9 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),199.999);}
	  else if(thePlot == 10 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min((double)eventEvent.rho,39.999);}
	  else if(thePlot == 11 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	  else if(thePlot == 12 && passAllCuts[SIGSEL]){makePlot = true;theVar = dPhiJetMET*180/TMath::Pi();}
	  else if(thePlot == 13 && passAllCuts[SIGSEL]){makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	  else if(thePlot == 14 && passAllCuts[SIGSEL]){makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]]))*180/TMath::Pi();}

	  if(makePlot) histo[thePlot][theCategory]->Fill(theVar,totalWeight);
	}
      }
    }
    printf("eff_cuts: %f\n",sumEventsProcess[ifile]);

  } // end of chain

  // "-1" to remove the Higgs contribution
  double sumEvents = 0;
  for(int np=1; np<histBins-1; np++) sumEvents += histo[0][np]->GetSumOfWeights();
  printf("yields: %f |",histo[0][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[0][np]->GetSumOfWeights());
  printf(" = %.3f\n",sumEvents);

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    char output[200];
    sprintf(output,"histozh_nice_%d.root",thePlot);	  
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }

  printf("                    em                      mm                      ee                      ll\n");
  printf("-----------------------------------------------------------------------------------------------------------\n");
  for(int ns=0; ns<nSelTypes; ns++) {
    printf("Selection: %s\n",selTypeName[ns].Data());
    double sumEventsType[4] = {0,0,0,0}; double sumEventsTypeE[4] = {0,0,0,0};
    for(int np=0; np<histBins; np++) {       
       bgdDecay[ns+nSelTypes*3][np] = bgdDecay[ns+nSelTypes*1][np] + bgdDecay[ns+nSelTypes*2][np];
       weiDecay[ns+nSelTypes*3][np] = weiDecay[ns+nSelTypes*1][np] + weiDecay[ns+nSelTypes*2][np];
       printf("(%6s): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
       processName[np].Data(),bgdDecay[ns+nSelTypes*0][np],sqrt(weiDecay[ns+nSelTypes*0][np]),bgdDecay[ns+nSelTypes*1][np],sqrt(weiDecay[ns+nSelTypes*1][np]),
                              bgdDecay[ns+nSelTypes*2][np],sqrt(weiDecay[ns+nSelTypes*2][np]),bgdDecay[ns+nSelTypes*3][np],sqrt(weiDecay[ns+nSelTypes*3][np]));
       if(np!=0 && np!=6){
         sumEventsType[0] = sumEventsType[0] + bgdDecay[ns+nSelTypes*0][np]; sumEventsTypeE[0] = sumEventsTypeE[0] + weiDecay[ns+nSelTypes*0][np];
         sumEventsType[1] = sumEventsType[1] + bgdDecay[ns+nSelTypes*1][np]; sumEventsTypeE[1] = sumEventsTypeE[1] + weiDecay[ns+nSelTypes*1][np];
         sumEventsType[2] = sumEventsType[2] + bgdDecay[ns+nSelTypes*2][np]; sumEventsTypeE[2] = sumEventsTypeE[2] + weiDecay[ns+nSelTypes*2][np];
         sumEventsType[3] = sumEventsType[3] + bgdDecay[ns+nSelTypes*3][np]; sumEventsTypeE[3] = sumEventsTypeE[3] + weiDecay[ns+nSelTypes*3][np];
       }
    }
    printf("(...bkg): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
           sumEventsType[0],sqrt(sumEventsTypeE[0]),sumEventsType[1],sqrt(sumEventsTypeE[1]),
	   sumEventsType[2],sqrt(sumEventsTypeE[2]),sumEventsType[3],sqrt(sumEventsTypeE[3]));
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }

  // This is the alternative method
  double kEoverM  = sqrt(bgdDecay[ZSEL+nSelTypes*2][0]/bgdDecay[ZSEL+nSelTypes*1][0]);
  double NemFact[3] = {(1.0/kEoverM)*0.5, (kEoverM)*0.5, (kEoverM+1.0/kEoverM)*0.5};

  // This is the default method
  double NemFact_FromMLLSB[3] = {NemFact[0], NemFact[1], NemFact[2]}; double NemFact_FromMLLSBE[3] = {1.0, 1.0, 1.0};
  if(bgdDecay[WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[WWLOOSESEL+nSelTypes*1][0]) > 0) {
    NemFact_FromMLLSB[0]  = (bgdDecay[WWLOOSESEL+nSelTypes*1][0])/bgdDecay[WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[0] = sqrt(1./bgdDecay[WWLOOSESEL+nSelTypes*1][0]+1./bgdDecay[WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[0];
  }
  if(bgdDecay[WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[WWLOOSESEL+nSelTypes*2][0]) > 0) {
    NemFact_FromMLLSB[1]  = (bgdDecay[WWLOOSESEL+nSelTypes*2][0])/bgdDecay[WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[1] = sqrt(1./bgdDecay[WWLOOSESEL+nSelTypes*2][0]+1./bgdDecay[WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[1];
  }
  if(bgdDecay[WWLOOSESEL+nSelTypes*0][0] > 0 && (bgdDecay[WWLOOSESEL+nSelTypes*3][0]) > 0) {
    NemFact_FromMLLSB[2]  = (bgdDecay[WWLOOSESEL+nSelTypes*3][0])/bgdDecay[WWLOOSESEL+nSelTypes*0][0];
    NemFact_FromMLLSBE[2] = sqrt(1./bgdDecay[WWLOOSESEL+nSelTypes*3][0]+1./bgdDecay[WWLOOSESEL+nSelTypes*0][0])*NemFact_FromMLLSB[2];
  }

  printf("(mm) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[0],NemFact_FromMLLSB[0],NemFact_FromMLLSBE[0]);
  printf("(ee) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[1],NemFact_FromMLLSB[1],NemFact_FromMLLSBE[1]);
  printf("(ll) kEoverM: %f ---> NemFact: %f | NemFact_FromMLLSB = %f +/- %f\n", kEoverM,NemFact[2],NemFact_FromMLLSB[2],NemFact_FromMLLSBE[2]);

  // There uncertainties: closure test, different between default and alternative method, and data statistics
  double EMSystTotal[3] = {1.0,1.0,1.0}; double EMSyst[2][3] = {0.0,0.0,0.0,0.0,0.0,0.0};

  EMSyst[0][0] = bgdDecay[SIGSEL+nSelTypes*1][1]/bgdDecay[SIGSEL+nSelTypes*0][1]/NemFact_FromMLLSB[0];
  if(EMSyst[0][0] < 1.0) EMSyst[0][0] = 1/EMSyst[0][0]; EMSyst[0][0] = EMSyst[0][0] - 1.0;

  EMSyst[0][1] = bgdDecay[SIGSEL+nSelTypes*2][1]/bgdDecay[SIGSEL+nSelTypes*0][1]/NemFact_FromMLLSB[1];
  if(EMSyst[0][1] < 1.0) EMSyst[0][1] = 1/EMSyst[0][1]; EMSyst[0][1] = EMSyst[0][1] - 1.0;

  EMSyst[0][2] = bgdDecay[SIGSEL+nSelTypes*3][1]/bgdDecay[SIGSEL+nSelTypes*0][1]/NemFact_FromMLLSB[2];
  if(EMSyst[0][2] < 1.0) EMSyst[0][2] = 1/EMSyst[0][2]; EMSyst[0][2] = EMSyst[0][2] - 1.0;

  EMSyst[1][0] = NemFact_FromMLLSB[0]/NemFact[0];
  if(EMSyst[1][0] < 1.0) EMSyst[1][0] = 1/EMSyst[1][0]; EMSyst[1][0] = EMSyst[1][0] - 1.0;

  EMSyst[1][1] = NemFact_FromMLLSB[1]/NemFact[1];
  if(EMSyst[1][1] < 1.0) EMSyst[1][1] = 1/EMSyst[1][1]; EMSyst[1][1] = EMSyst[1][1] - 1.0;

  EMSyst[1][2] = NemFact_FromMLLSB[2]/NemFact[2];
  if(EMSyst[1][2] < 1.0) EMSyst[1][2] = 1/EMSyst[1][2]; EMSyst[1][2] = EMSyst[1][2] - 1.0;

  if(bgdDecay[SIGSEL][0] > 0) EMSystTotal[0] = sqrt(EMSyst[0][0]*EMSyst[0][0] + EMSyst[1][0]*EMSyst[1][0] + 1/bgdDecay[SIGSEL][0]);
  else                        EMSystTotal[0] = sqrt(EMSyst[0][0]*EMSyst[0][0] + EMSyst[1][0]*EMSyst[1][0] + 1.0);

  if(bgdDecay[SIGSEL][0] > 0) EMSystTotal[1] = sqrt(EMSyst[0][1]*EMSyst[0][1] + EMSyst[1][1]*EMSyst[1][1] + 1/bgdDecay[SIGSEL][0]);
  else                        EMSystTotal[1] = sqrt(EMSyst[0][1]*EMSyst[0][1] + EMSyst[1][1]*EMSyst[1][1] + 1.0);

  if(bgdDecay[SIGSEL][0] > 0) EMSystTotal[2] = sqrt(EMSyst[0][2]*EMSyst[0][2] + EMSyst[1][2]*EMSyst[1][2] + 1/bgdDecay[SIGSEL][0]);
  else                        EMSystTotal[2] = sqrt(EMSyst[0][2]*EMSyst[0][2] + EMSyst[1][2]*EMSyst[1][2] + 1.0);

  double EMbkg = bgdDecay[SIGSEL][2]+bgdDecay[SIGSEL][3]+bgdDecay[SIGSEL][4]+bgdDecay[SIGSEL][5];

  printf("(mm) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[SIGSEL+nSelTypes*1][1],sqrt(weiDecay[SIGSEL+nSelTypes*1][4]),
	 bgdDecay[SIGSEL][1]*NemFact_FromMLLSB[0]  ,sqrt(weiDecay[SIGSEL][4])*NemFact_FromMLLSB[0],bgdDecay[SIGSEL][0],EMbkg,EMSystTotal[0],EMSyst[0][0],EMSyst[1][0],sqrt(EMSystTotal[0]*EMSystTotal[0] - EMSyst[0][0]*EMSyst[0][0] - EMSyst[1][0]*EMSyst[1][0]));

  printf("(ee) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[SIGSEL+nSelTypes*2][1],sqrt(weiDecay[SIGSEL+nSelTypes*2][4]),
	 bgdDecay[SIGSEL][1]*NemFact_FromMLLSB[1]  ,sqrt(weiDecay[SIGSEL][4])*NemFact_FromMLLSB[1],bgdDecay[SIGSEL][0],EMbkg,EMSystTotal[1],EMSyst[0][1],EMSyst[1][1],sqrt(EMSystTotal[1]*EMSystTotal[1] - EMSyst[0][1]*EMSyst[0][1] - EMSyst[1][1]*EMSyst[1][1]));

  printf("(ll) EM MC: %8.3f +/- %5.3f --> EM Prediction: %8.3f +/- %5.3f, EM data/bkg: %f/%f --> syst: %f (%f,%f,%f)\n",
         bgdDecay[SIGSEL+nSelTypes*3][1],sqrt(weiDecay[SIGSEL+nSelTypes*3][4]),
	 bgdDecay[SIGSEL][1]*NemFact_FromMLLSB[2]  ,sqrt(weiDecay[SIGSEL][4])*NemFact_FromMLLSB[2],bgdDecay[SIGSEL][0],EMbkg,EMSystTotal[2],EMSyst[0][2],EMSyst[1][2],sqrt(EMSystTotal[2]*EMSystTotal[2] - EMSyst[0][2]*EMSyst[0][2] - EMSyst[1][2]*EMSyst[1][2]));
}
