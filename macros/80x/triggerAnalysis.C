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
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"

#include "MitAnalysisRunII/macros/80x/factors.h"

void triggerAnalysis(
bool applyFullLepSel = false,
TString typeLepSel = "medium",
int typeAna = 0
){

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;

  if(typeAna == 0){
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/MET_Run2016B.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/MET_Run2016C.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/MET_Run2016D.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/MET_Run2016E.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/MET_Run2016F.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/MET_Run2016G.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/MET_Run2016H.root");

  } else {
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/merging_80x/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_0000.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/merging_80x/WWTo2L2Nu_13TeV-powheg_0000.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/merging_80x/GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8_0000.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/merging_80x/ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_0000.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/onlymet/DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/onlymet/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/onlymet/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/onlymet/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root");
    infilenamev.push_back("root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/merging_80x/onlymet/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_0000.root");
  }

  // Initializations
  double nPassTrigger[3][25];
  for(int i=0; i<3; i++){
    for(int j=0; j<25; j++){
      nPassTrigger[i][j] = 0.0;
    }
  }
  double nPassSel[4][2];
  for(int i=0; i<4; i++){
    for(int j=0; j<2; j++){
      nPassSel[i][j] = 0.0;
    }
  }

  TString triggerSuffix = "";

  double xminPlot = 0.0;
  double xmaxPlot = 1.0;
  int nBinPlot = 200;
  const int histBins = 3;
  const int allPlots = 6;
  TH1D* histo[histBins][allPlots];

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  5) {nBinPlot = 72; xminPlot =-0.5; xmaxPlot = 71.5;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[i][thePlot] = (TH1D*) histos->Clone(Form("histo_%d",i));
    histos->Clear();
  }

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());
    TFile *the_input_file = TFile::Open(infilenamev[ifile].Data());
    if(!the_input_file) continue;
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareJets eventJets;
    eventJets.setBranchAddresses(the_input_tree);

    BareLeptons eventLeptons;
    eventLeptons.setBranchAddresses(the_input_tree);

    BareMet eventMet;
    eventMet.SetExtend();
    eventMet.setBranchAddresses(the_input_tree);

    BareTrigger eventTrigger;
    eventTrigger.setBranchAddresses(the_input_tree);

    BareVertex eventVertex;
    eventVertex.setBranchAddresses(the_input_tree);

    BareMonteCarlo eventMonteCarlo;
    eventMonteCarlo.SetExtend();
    eventMonteCarlo.setBranchAddresses(the_input_tree);

    TNamed *triggerNames = (TNamed*)the_input_file->FindObjectAny("triggerNames");
    char **tokens;
    size_t numtokens;
    tokens = strsplit(triggerNames->GetTitle(), ",", &numtokens);

    double theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      the_input_tree->GetEntry(i);

      vector<int> idLep; vector<int> idTight; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
      }

      if(idLep.size() < 2) continue;

      if(applyFullLepSel == true && idLep.size() != goodIsTight) continue;

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) ));

      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 10 ||
         dilep.M() <= 12) continue;

      double totalWeight = 1.0;

      int typePair = -1;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) {typePair = 0;}
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) {typePair = 1;}
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])!=TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]))         {typePair = 2;}
      else {assert(1); return;}

      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 25.0) nPassSel[typePair][0] = nPassSel[typePair][0] + totalWeight;
      if(idLep.size() >= 3 && 
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 25.0) nPassSel[       3][0] = nPassSel[       3][0] + totalWeight;

      int iPt[2] = {-1, -1};
      if     (((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 25) iPt[0] = 0;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 30) iPt[0] = 1;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 50) iPt[0] = 2;
      else								  iPt[0] = 3;
      if     (((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() < 15) iPt[1] = 0;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() < 20) iPt[1] = 1;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() < 25) iPt[1] = 2;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() < 30) iPt[1] = 3;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() < 50) iPt[1] = 4;
      else								  iPt[1] = 5;

      int iEta[2] = {-1, -1};
      if    (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) < 1.5) iEta[0] = 0;
      else									       iEta[0] = 1;
      if    (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()) < 1.5) iEta[1] = 0;
      else									       iEta[1] = 1;

      double theBin = -1;
      if     (iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin = 0;
      else if(iPt[0] == 0 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin = 1;
      else if(iPt[0] == 0 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin = 2;
      else if(iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin = 3;
      else if(iPt[0] == 0 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin = 4;
      else if(iPt[0] == 0 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin = 5;
      else if(iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin = 6;
      else if(iPt[0] == 0 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin = 7;
      else if(iPt[0] == 0 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin = 8;
      else if(iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin = 9;
      else if(iPt[0] == 0 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =10;
      else if(iPt[0] == 0 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =11;

      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin =12;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin =13;
      else if(iPt[0] == 1 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin =14;
      else if(iPt[0] == 1 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 0) theBin =15;
      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin =16;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin =17;
      else if(iPt[0] == 1 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin =18;
      else if(iPt[0] == 1 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 1) theBin =19;
      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin =20;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin =21;
      else if(iPt[0] == 1 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin =22;
      else if(iPt[0] == 1 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 0) theBin =23;
      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin =24;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =25;
      else if(iPt[0] == 1 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =26;
      else if(iPt[0] == 1 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 1) theBin =27;

      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin =28;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin =29;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin =30;
      else if(iPt[0] == 2 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 0) theBin =31;
      else if(iPt[0] == 2 && iPt[1] == 4 && iEta[0] == 0 && iEta[1] == 0) theBin =32;
      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin =33;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin =34;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin =35;
      else if(iPt[0] == 2 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 1) theBin =36;
      else if(iPt[0] == 2 && iPt[1] == 4 && iEta[0] == 0 && iEta[1] == 1) theBin =37;
      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin =38;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin =39;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin =40;
      else if(iPt[0] == 2 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 0) theBin =41;
      else if(iPt[0] == 2 && iPt[1] == 4 && iEta[0] == 1 && iEta[1] == 0) theBin =42;
      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin =43;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =44;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =45;
      else if(iPt[0] == 2 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 1) theBin =46;
      else if(iPt[0] == 2 && iPt[1] == 4 && iEta[0] == 1 && iEta[1] == 1) theBin =47;

      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin =48;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin =49;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin =50;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 0) theBin =51;
      else if(iPt[0] == 3 && iPt[1] == 4 && iEta[0] == 0 && iEta[1] == 0) theBin =52;
      else if(iPt[0] == 3 && iPt[1] == 5 && iEta[0] == 0 && iEta[1] == 0) theBin =53;
      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin =54;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin =55;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin =56;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 1) theBin =57;
      else if(iPt[0] == 3 && iPt[1] == 4 && iEta[0] == 0 && iEta[1] == 1) theBin =58;
      else if(iPt[0] == 3 && iPt[1] == 5 && iEta[0] == 0 && iEta[1] == 1) theBin =59;
      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin =60;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin =61;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin =62;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 0) theBin =63;
      else if(iPt[0] == 3 && iPt[1] == 4 && iEta[0] == 1 && iEta[1] == 0) theBin =64;
      else if(iPt[0] == 3 && iPt[1] == 5 && iEta[0] == 1 && iEta[1] == 0) theBin =65;
      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin =66;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =67;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =68;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 1) theBin =69;
      else if(iPt[0] == 3 && iPt[1] == 4 && iEta[0] == 1 && iEta[1] == 1) theBin =70;
      else if(iPt[0] == 3 && iPt[1] == 5 && iEta[0] == 1 && iEta[1] == 1) theBin =71;
 
      else {printf("IMPOSSIBLE\n");}

      Bool_t passFilter = kFALSE;
      Bool_t passMETFilter = kFALSE;
      for (int nt = 0; nt <(int)numtokens; nt++) {
        if((*eventTrigger.triggerFired)[nt] == 0) continue;
        if((strcmp(tokens[nt],Form("HLT_Ele25_eta2p1_WPTight_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Ele27_eta2p1_WPLoose_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_IsoMu24_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_IsoTkMu24_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Ele27_WPTight_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Ele30_WPTight_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Ele35_WPLoose_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_IsoMu22_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_IsoTkMu22_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu45_eta2p1_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu50_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0)
           ) passFilter = kTRUE;

        if((strcmp(tokens[nt],Form("HLT_PFMET170_NoiseCleaned%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMETNoMu110_NoiseCleaned_PFMHTNoMu110_IDTight%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMET170_HBHECleaned%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMET170_JetIdCleaned%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMET170_NotCleaned%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMET170_HBHE_BeamHaloCleaned%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140_v%s",triggerSuffix.Data()))  == 0) ||
           (strcmp(tokens[nt],Form("HLT_PFMET120_PFMHT120_IDTight%s",triggerSuffix.Data()))  == 0)
           ) passMETFilter = kTRUE;

        if(strcmp(tokens[nt],Form("HLT_Ele25_eta2p1_WPTight_Gsf_v%s",triggerSuffix.Data()))			 == 0) nPassTrigger[typePair][ 0]++;
        if(strcmp(tokens[nt],Form("HLT_Ele27_eta2p1_WPLoose_Gsf_v%s",triggerSuffix.Data()))			 == 0) nPassTrigger[typePair][ 1]++;
        if(strcmp(tokens[nt],Form("HLT_IsoMu24_v%s",triggerSuffix.Data()))					 == 0) nPassTrigger[typePair][ 2]++;
        if(strcmp(tokens[nt],Form("HLT_IsoTkMu24_v%s",triggerSuffix.Data()))					 == 0) nPassTrigger[typePair][ 3]++;
        if(strcmp(tokens[nt],Form("HLT_Ele27_WPTight_Gsf_v%s",triggerSuffix.Data()))				 == 0) nPassTrigger[typePair][ 4]++;
        if(strcmp(tokens[nt],Form("HLT_Ele30_WPTight_Gsf_v%s",triggerSuffix.Data()))				 == 0) nPassTrigger[typePair][ 5]++;
        if(strcmp(tokens[nt],Form("HLT_Ele35_WPLoose_Gsf_v%s",triggerSuffix.Data()))				 == 0) nPassTrigger[typePair][ 6]++;
        if(strcmp(tokens[nt],Form("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))	 == 0) nPassTrigger[typePair][ 7]++;
        if(strcmp(tokens[nt],Form("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v%s",triggerSuffix.Data()))		 == 0) nPassTrigger[typePair][ 8]++;
        if(strcmp(tokens[nt],Form("HLT_IsoMu22_v%s",triggerSuffix.Data()))					 == 0) nPassTrigger[typePair][ 9]++;
        if(strcmp(tokens[nt],Form("HLT_IsoTkMu22_v%s",triggerSuffix.Data()))					 == 0) nPassTrigger[typePair][10]++;
        if(strcmp(tokens[nt],Form("HLT_Mu45_eta2p1_v%s",triggerSuffix.Data()))  				 == 0) nPassTrigger[typePair][11]++;
        if(strcmp(tokens[nt],Form("HLT_Mu50_v%s",triggerSuffix.Data())) 					 == 0) nPassTrigger[typePair][12]++;
        if(strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v%s",triggerSuffix.Data())) 		 == 0) nPassTrigger[typePair][13]++;
        if(strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v%s",triggerSuffix.Data()))		 == 0) nPassTrigger[typePair][14]++;
        if(strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v%s",triggerSuffix.Data()))		 == 0) nPassTrigger[typePair][15]++;
        if(strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v%s",triggerSuffix.Data()))		 == 0) nPassTrigger[typePair][16]++;
        if(strcmp(tokens[nt],Form("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))== 0) nPassTrigger[typePair][17]++;
        if(strcmp(tokens[nt],Form("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))   == 0) nPassTrigger[typePair][18]++;
        if(strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))== 0) nPassTrigger[typePair][19]++;
        if(strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))   == 0) nPassTrigger[typePair][20]++;
        if(strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data())) == 0) nPassTrigger[typePair][21]++;
        if(strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))	 == 0) nPassTrigger[typePair][22]++;
        if(strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data())) == 0) nPassTrigger[typePair][23]++;
        if(strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))	 == 0) nPassTrigger[typePair][24]++;

      } // end trigger paths loop

      histo[typePair][0]->Fill(theBin,totalWeight);

      if(passMETFilter == true) histo[typePair][2]->Fill(theBin,totalWeight);
      else                      histo[typePair][4]->Fill(theBin,totalWeight);

      if(passFilter){
        if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 25.0) nPassSel[typePair][1] = nPassSel[typePair][1] + totalWeight;
        if(idLep.size() >= 3 && 
           ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 25.0) nPassSel[       3][1] = nPassSel[    3][1] + totalWeight;
        histo[typePair][1]->Fill(theBin,totalWeight);
        if(passMETFilter == true) histo[typePair][3]->Fill(theBin,totalWeight);
        else                      histo[typePair][5]->Fill(theBin,totalWeight);
      }
    } // end events loop
  } // end chain loop

  for(int i=0; i<4; i++){
    if(nPassSel[i][0] > 0){
      printf("trigger_eff(%d): %f +/- %f\n",i,nPassSel[i][1]/nPassSel[i][0],sqrt(nPassSel[i][1]/nPassSel[i][0]*(1-nPassSel[i][1]/nPassSel[i][0])/nPassSel[i][0]));
      if(i != 3){
        printf("trigger_cuts: ");
        for(unsigned int nc=0; nc<25; nc++){
          printf("%4.1f ",100*(double)nPassTrigger[i][nc]/nPassSel[i][0]);
        }
        printf("\n");
      }
    }
  } // end loop over typePair (4)

  for(int theType=0; theType<histBins; theType++){
    histo[theType][1]->Divide(histo[theType][0]);
    histo[theType][3]->Divide(histo[theType][2]);
    histo[theType][5]->Divide(histo[theType][4]);
  }
  char output[200];
  sprintf(output,"histo_trigger_study_eff_%d_%s_%d.root",(int)applyFullLepSel,typeLepSel.Data(),typeAna);	
  TFile* outFilePlotsEff = new TFile(output,"recreate");
  outFilePlotsEff->cd();
  for(int theType=0; theType<histBins; theType++){
    histo[theType][1]->Write();
  }
  outFilePlotsEff->Close();

  sprintf(output,"histo_trigger_study_evt_%d_%s_%d.root",(int)applyFullLepSel,typeLepSel.Data(),typeAna);
  TFile* outFilePlotsEvt = new TFile(output,"recreate");
  outFilePlotsEvt->cd();
  for(int theType=0; theType<histBins; theType++){
    histo[theType][0]->Write();
  }
  outFilePlotsEvt->Close();

  sprintf(output,"histo_trigger_study_eff_%d_%s_%d_met.root",(int)applyFullLepSel,typeLepSel.Data(),typeAna);	
  TFile* outFilePlotsEffMET = new TFile(output,"recreate");
  outFilePlotsEffMET->cd();
  for(int theType=0; theType<histBins; theType++){
    histo[theType][3]->Write();
  }
  outFilePlotsEffMET->Close();

  sprintf(output,"histo_trigger_study_evt_%d_%s_%d_met.root",(int)applyFullLepSel,typeLepSel.Data(),typeAna);
  TFile* outFilePlotsEvtMET = new TFile(output,"recreate");
  outFilePlotsEvtMET->cd();
  for(int theType=0; theType<histBins; theType++){
    histo[theType][2]->Write();
  }
  outFilePlotsEvtMET->Close();

  sprintf(output,"histo_trigger_study_eff_%d_%s_%d_nomet.root",(int)applyFullLepSel,typeLepSel.Data(),typeAna);	
  TFile* outFilePlotsEffNOMET = new TFile(output,"recreate");
  outFilePlotsEffNOMET->cd();
  for(int theType=0; theType<histBins; theType++){
    histo[theType][5]->Write();
  }
  outFilePlotsEffNOMET->Close();

  sprintf(output,"histo_trigger_study_evt_%d_%s_%d_nomet.root",(int)applyFullLepSel,typeLepSel.Data(),typeAna);
  TFile* outFilePlotsEvtNOMET = new TFile(output,"recreate");
  outFilePlotsEvtNOMET->cd();
  for(int theType=0; theType<histBins; theType++){
    histo[theType][4]->Write();
  }
  outFilePlotsEvtNOMET->Close();
}
