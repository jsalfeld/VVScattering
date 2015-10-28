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

double mcPrescale = 1.0;
const bool useDYMVA = false;

void wwAnalysis(
 Int_t period = 1
 ){

  TString filesPath  = "/scratch5/ceballos/ntuples_weights/";
  Double_t lumi = 1;

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
  //infilenamev.push_back(Form("%sdata_AOD_25ns.root",filesPath.Data()));												          infilecatv.push_back(0);
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	  infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3+AODSIM.root",filesPath.Data()));		  infilecatv.push_back(2);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	  infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));					  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));					  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		  infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		  infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_amcatnloFXFX_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		  infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  		  infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); 			  infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				  infilecatv.push_back(7);
  /////infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); 		  infilecatv.push_back(7);
  /////infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); infilecatv.push_back(7);

  //infilenamev.push_back(Form("%sHWminusJ_HToWW_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  			  infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sHWplusJ_HToWW_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				  infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sHZJ_HToWW_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				  infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sGluGluZH_HToWW_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  			  infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sWplusHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  			  infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sWminusHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); 			  infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sZHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				  infilecatv.push_back(7);
  }
  else {assert(0);}

  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
 
  float dymva_= -999.;
  unsigned int nlep_= -1;
  unsigned int njets_= -1;

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

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

    const int numberCuts = 10;
    TH1D* hDWWLL[3];
    hDWWLL[0] = new TH1D("hDWWLL0", "hDWWLL0", numberCuts, -0.5, numberCuts-0.5);
    hDWWLL[1] = new TH1D("hDWWLL1", "hDWWLL1", numberCuts, -0.5, numberCuts-0.5);
    hDWWLL[2] = new TH1D("hDWWLL2", "hDWWLL2", numberCuts, -0.5, numberCuts-0.5);
 
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

    hDWWLL[0]->Clear();hDWWLL[1]->Clear();hDWWLL[2]->Clear();
    TString cutName[numberCuts] = {"ptl>20/20", "3rd lepton veto", "mll>12", "MET>20", "Z veto", "minPMET>20/45", "ptll>30/45", "btag-veto", "Nsoftmu=0", "Njets=0"};
    for (int i=0; i<int(the_input_tree->GetEntries()/1); ++i) {
      the_input_tree->GetEntry(i);

      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      int sumEvol[3] = {-1, -1, -1};
      Bool_t passFilter[numberCuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      Bool_t passPresel = kFALSE;
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
           ) passPresel = kTRUE;
      }

      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut("default",TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}
      }

      if(idLep.size()!=idTight.size()) {assert(1); return;}

      if(goodIsTight >= 2) passPresel = passPresel && kTRUE;
      else                 passPresel = passPresel && kFALSE;

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

      if(goodIsTight == idTight.size() && idTight.size() == 2) passFilter[1] = kTRUE;

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

       passFilter[0] = passFilter[0] * (signQ == 0);

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 
      if(dilep.M() > 12) passFilter[2] = kTRUE;

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

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 15 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) continue;
        
	idJet.push_back(nj);
      }

      if(useDYMVA == true){
        if(nlep_ != idLep.size()) {printf("PROBLEM nlep %d != %d\n",(int)nlep_,(int)idLep.size()); assert(1); return;}
        if(njets_ != idJet.size()) {printf("PROBLEM njet %d != %d\n",(int)njets_,(int)idJet.size()); assert(1); return;}
      }

      if(((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 20) passFilter[3] = kTRUE;
      passFilter[4] = kTRUE;
      if(minPMET > 20) passFilter[5] = kTRUE;
      if(dilep.Pt() > 30) passFilter[6] = kTRUE;
      if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])){
        if(TMath::Abs(dilep.M()-91.1876) <= 15.0) passFilter[4] = kFALSE;
        if(minPMET <= 45) passFilter[5] = kFALSE;
        if(dilep.Pt() <= 45) passFilter[6] = kFALSE;
      }
      if(bDiscrMax < 0.605) passFilter[7] = kTRUE;
      if(idSoft.size() == 0) passFilter[8] = kTRUE;
      if(idJet.size() == 0) passFilter[9] = kTRUE;

      bool totalSel = kTRUE;
      for(int isel=0; isel<numberCuts; isel++) {
        totalSel = totalSel && passFilter[isel];
	if(totalSel == kTRUE) sumEvol[typeSel]++;
      }

      // luminosity
      double theLumi  = 1.0; if(infilecatv[ifile] != 0) theLumi  = lumi;
      // pile-up
      double puWeight = 1.0; if(infilecatv[ifile] != 0) puWeight = nPUScaleFactor(fhDPU, (double)eventVertex.npv);

      double totalWeight = eventMonteCarlo.mcWeight*theLumi;

      for(int nl=0; nl <=sumEvol[typeSel]; nl++) hDWWLL[typeSel]->Fill((double)nl,totalWeight);
    }
    printf("                   mm        ee        em\n");
    for(int nc=0; nc<numberCuts; nc++){
      printf("(%15s): %7.2f   %7.2f   %7.2f\n",cutName[nc].Data(),hDWWLL[0]->GetBinContent(nc+1),hDWWLL[1]->GetBinContent(nc+1),hDWWLL[2]->GetBinContent(nc+1));
    }

  } // end of chain
}
