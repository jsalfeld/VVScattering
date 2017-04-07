#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"
#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"

#include "MitAnalysisRunII/macros/80x/factors.h"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"

int whichSkim = 0;
double mcPrescale = 1.0;
bool isMINIAOD = true;

void ZllAnalysis(
 TString typeLepSel = "medium"
){

  TString filesPathDA = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/";
  TString filesPathMC  = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/";
  Double_t lumi = 35.9;

  const int etaBins = 5;
  const int ptBins = 8;

  double denFRDAMU[etaBins][ptBins] = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};
  double numFRDAMU[etaBins][ptBins] = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};
  double denFRBGMU[etaBins][ptBins] = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};
  double numFRBGMU[etaBins][ptBins] = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};
  double denFRDAEL[etaBins][ptBins] = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};
  double numFRDAEL[etaBins][ptBins] = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};
  double denFRBGEL[etaBins][ptBins] = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};
  double numFRBGEL[etaBins][ptBins] = {0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString triggerSuffix = "*";
  if(isMINIAOD) triggerSuffix = "";
  TString puPath = "MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root";
  if(isMINIAOD) {
    infilenamev.push_back(Form("%sdata_Run2016B.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016C.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016D.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016E.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016F.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016G.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016H.root",filesPathDA.Data())); infilecatv.push_back(0);
  } else {
  }

  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data())); infilecatv.push_back(1);

  if(infilenamev.size() != infilecatv.size()) assert(0);

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
  delete fPUFile;

  TFile *fTrackElectronReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root"));
  TH2D *fhDeltrksf= (TH2D*)(fTrackElectronReco_SF->Get("scalefactors_Reco_Electron")); assert(fhDeltrksf); fhDeltrksf->SetDirectory(0);
  delete fTrackElectronReco_SF;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("scalefactors_Medium_Electron"));
  TH2D *fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_Tight_Electron"));
  if(typeLepSel == "medium_mva") fhDElMediumSF = (TH2D*)(fElSF->Get("scalefactors_MediumMVA_Electron"));
  if(typeLepSel == "default_mva") fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_TightMVA_Electron"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF->SetDirectory(0);
  delete fElSF;

  TFile *fElVeryTightSF = TFile::Open(Form("MitAnalysisRunII/data/80x/veryTightSF_37ifb.root"));
  TH1D *fhDVeryTightSF = (TH1D*)(fElVeryTightSF->Get("veryTightSF"));
  assert(fhDVeryTightSF);
  fhDVeryTightSF->SetDirectory(0);
  delete fElVeryTightSF;

  TFile *fTrackMuonReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/Tracking_EfficienciesAndSF_BCDEFGH.root"));
  TH1D *fhDmutrksfptg10 = (TH1D*)(fTrackMuonReco_SF->Get("ratio_eff_eta3_dr030e030_corr")); assert(fhDmutrksfptg10); fhDmutrksfptg10->SetDirectory(0);
  delete fTrackMuonReco_SF;

  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("scalefactors_TightId_Muon")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  //TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonID_Z_RunBCD_prompt80X_7p65.root"));
  //TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  delete fMuSF;

  TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("scalefactors_Iso_MuonTightId")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
  //TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonIso_Z_RunBCD_prompt80X_7p65.root"));
  //TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
  delete fMuIsoSF;

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());

    TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    //TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");
    TTree *the_SelBit_tree= (TTree*)the_input_file.FindObjectAny("SelBit_tree");

    BareMonteCarlo eventMonteCarlo;
    eventMonteCarlo.setBranchAddresses(the_input_tree);

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareLeptons eventLeptons;
    eventLeptons.setBranchAddresses(the_input_tree);

    BareTrigger eventTrigger;
    eventTrigger.setBranchAddresses(the_input_tree);

    BareVertex eventVertex;
    eventVertex.setBranchAddresses(the_input_tree);

    TNamed *triggerNames = (TNamed*)the_input_file.FindObjectAny("triggerNames");
    char **tokens;
    size_t numtokens;
    tokens = strsplit(triggerNames->GetTitle(), ",", &numtokens);
    if(ifile == 0){
      for (int i = 0; i < (int)numtokens; i++) {
        printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
      }
    }
    else {
    }

    unsigned int selBit_= 0;
    the_SelBit_tree->SetBranchAddress("selBit", &selBit_);
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    if(the_input_tree->GetEntries() != the_SelBit_tree->GetEntries()) {printf("BIG SKIMMING FAILURE\n"); return;}
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_SelBit_tree->GetEntry(i);
      if(i%1000000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      if((selBit_ & 0x1<<whichSkim) == 0) continue;

      the_input_tree->GetEntry(i);

      Bool_t passFilter[5] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
      if(infilecatv[ifile] != 999) {
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
             ) passFilter[1] = kTRUE;
	}
      } else { passFilter[1] = kTRUE;}

      if(passFilter[0] == kFALSE) continue;
      if(passFilter[1] == kFALSE) continue;

      vector<int> idLep; vector<int> idTight;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
      }
      if(idLep.size()==2) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;

      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 25 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 12) continue;

      passFilter[3] = (int)(*eventLeptons.pdgId)[idLep[0]]*(int)(*eventLeptons.pdgId)[idLep[1]] < 0;
      if(passFilter[3] == kFALSE) continue;

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 

      if(TMath::Abs(dilep.M()-91.1876)<15.0) passFilter[4] = kTRUE;  	    
      if(passFilter[4] == kFALSE) continue;

      int typeSel = -1;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) {typeSel = 0;}
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) {typeSel = 1;}
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])!=TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]))	  {typeSel = 2;}
      else {printf("IMPOSSIBLE TYPESEL!\n"); return;}

      int iPt[2] = {-1, -1};
      if     (((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <  15) iPt[0] = 0;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <  20) iPt[0] = 1;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <  25) iPt[0] = 2;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <  30) iPt[0] = 3;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <  40) iPt[0] = 4;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <  60) iPt[0] = 5;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 100) iPt[0] = 6;
      else                                                                 iPt[0] = 7;

      if     (((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <  15) iPt[1] = 0;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <  20) iPt[1] = 1;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <  25) iPt[1] = 2;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <  30) iPt[1] = 3;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <  40) iPt[1] = 4;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <  60) iPt[1] = 5;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() < 100) iPt[1] = 6;
      else                                                                 iPt[1] = 7;

      int iEta[2] = {-1, -1};
      if     (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) < 0.5) iEta[0] = 0;
      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) < 1.0) iEta[0] = 1;
      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) < 1.5) iEta[0] = 2;
      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) < 2.0) iEta[0] = 3;
      else                                                                              iEta[0] = 4;

      if     (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()) < 0.5) iEta[1] = 0;
      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()) < 1.0) iEta[1] = 1;
      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()) < 1.5) iEta[1] = 2;
      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()) < 2.0) iEta[1] = 3;
      else                                                                              iEta[1] = 4;

      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
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
		typeLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,true,fhDMuIsoSF,fhDVeryTightSF,true);
        }
      }

      double totalWeight = mcWeight*theLumi*puWeight*effSF*theMCPrescale;
      if(typeSel == 2) totalWeight = -0.5 * totalWeight;

      if(typeSel == 0 || typeSel == 2){ // mm or em
	if(infilecatv[ifile] == 0) {
          denFRDAMU[iEta[0]][iPt[0]] = denFRDAMU[iEta[0]][iPt[0]] + totalWeight;
          denFRDAMU[iEta[1]][iPt[1]] = denFRDAMU[iEta[1]][iPt[1]] + totalWeight;
          if(idTight[0]) numFRDAMU[iEta[0]][iPt[0]] = numFRDAMU[iEta[0]][iPt[0]] + totalWeight;
          if(idTight[1]) numFRDAMU[iEta[1]][iPt[1]] = numFRDAMU[iEta[1]][iPt[1]] + totalWeight;
	}
	else {
          denFRBGMU[iEta[0]][iPt[0]] = denFRBGMU[iEta[0]][iPt[0]] + totalWeight;
          denFRBGMU[iEta[1]][iPt[1]] = denFRBGMU[iEta[1]][iPt[1]] + totalWeight;
          if(idTight[0]) numFRBGMU[iEta[0]][iPt[0]] = numFRBGMU[iEta[0]][iPt[0]] + totalWeight;
          if(idTight[1]) numFRBGMU[iEta[1]][iPt[1]] = numFRBGMU[iEta[1]][iPt[1]] + totalWeight;
	}
      }

      if(typeSel == 1 || typeSel == 2){ // ee or em
	if(infilecatv[ifile] == 0) {
          denFRDAEL[iEta[0]][iPt[0]] = denFRDAEL[iEta[0]][iPt[0]] + totalWeight;
          denFRDAEL[iEta[1]][iPt[1]] = denFRDAEL[iEta[1]][iPt[1]] + totalWeight;
          if(idTight[0]) numFRDAEL[iEta[0]][iPt[0]] = numFRDAEL[iEta[0]][iPt[0]] + totalWeight;
          if(idTight[1]) numFRDAEL[iEta[1]][iPt[1]] = numFRDAEL[iEta[1]][iPt[1]] + totalWeight;
	}
	else {
          denFRBGEL[iEta[0]][iPt[0]] = denFRBGEL[iEta[0]][iPt[0]] + totalWeight;
          denFRBGEL[iEta[1]][iPt[1]] = denFRBGEL[iEta[1]][iPt[1]] + totalWeight;
          if(idTight[0]) numFRBGEL[iEta[0]][iPt[0]] = numFRBGEL[iEta[0]][iPt[0]] + totalWeight;
          if(idTight[1]) numFRBGEL[iEta[1]][iPt[1]] = numFRBGEL[iEta[1]][iPt[1]] + totalWeight;
	}
      }
    }
  } // end of chain

  double sumTot[2] = {0.,0.};
  printf("*********Muons*********\n");
  sumTot[0] = 0.; sumTot[1] = 0.;
  for(int iEta=0; iEta<etaBins; iEta++){
    for(int iPt=0; iPt<ptBins; iPt++){
      sumTot[0] = sumTot[0] + numFRDAMU[iEta][iPt];
      sumTot[1] = sumTot[1] + denFRDAMU[iEta][iPt];
      printf("(%d,%d): %9.1f/%9.1f=%4.3f | ",iPt,iEta,numFRDAMU[iEta][iPt],denFRDAMU[iEta][iPt],numFRDAMU[iEta][iPt]/denFRDAMU[iEta][iPt]);
      if(iPt==ptBins-1) printf("\n");
    }
  }
  printf("sumTotDA = %f / %f = %f\n",sumTot[0],sumTot[1],sumTot[0]/sumTot[1]);

  sumTot[0] = 0.; sumTot[1] = 0.;
  for(int iEta=0; iEta<etaBins; iEta++){
    for(int iPt=0; iPt<ptBins; iPt++){
      sumTot[0] = sumTot[0] + numFRBGMU[iEta][iPt];
      sumTot[1] = sumTot[1] + denFRBGMU[iEta][iPt];
      printf("(%d,%d): %9.1f/%9.1f=%4.3f | ",iPt,iEta,numFRBGMU[iEta][iPt],denFRBGMU[iEta][iPt],numFRBGMU[iEta][iPt]/denFRBGMU[iEta][iPt]);
      if(iPt==ptBins-1) printf("\n");
    }
  }
  printf("sumTotBG = %f / %f = %f\n",sumTot[0],sumTot[1],sumTot[0]/sumTot[1]);

  printf("double prompt_rate_DAMU[%d][%d] = {\n",etaBins,ptBins);
  for(int iEta=0; iEta<etaBins; iEta++){
    for(int iPt=0; iPt<ptBins; iPt++){
      printf("%4.3f",numFRDAMU[iEta][iPt]/denFRDAMU[iEta][iPt]);
      if(iPt!=ptBins-1||iEta!=etaBins-1) printf(",");
      if(iPt==ptBins-1) printf("\n");
    }
  }
  printf("};\n");

  printf("double prompt_rate_BGMU[%d][%d] = {\n",etaBins,ptBins);
  for(int iEta=0; iEta<etaBins; iEta++){
    for(int iPt=0; iPt<ptBins; iPt++){
      printf("%4.3f",numFRBGMU[iEta][iPt]/denFRBGMU[iEta][iPt]);
      if(iPt!=ptBins-1||iEta!=etaBins-1) printf(",");
      if(iPt==ptBins-1) printf("\n");
    }
  }
  printf("};\n");

  printf("*********Electrons*********\n");
  sumTot[0] = 0.; sumTot[1] = 0.;
  for(int iEta=0; iEta<etaBins; iEta++){
    for(int iPt=0; iPt<ptBins; iPt++){
      sumTot[0] = sumTot[0] + numFRDAEL[iEta][iPt];
      sumTot[1] = sumTot[1] + denFRDAEL[iEta][iPt];
      printf("(%d,%d): %9.1f/%9.1f=%4.3f | ",iPt,iEta,numFRDAEL[iEta][iPt],denFRDAEL[iEta][iPt],numFRDAEL[iEta][iPt]/denFRDAEL[iEta][iPt]);
      if(iPt==ptBins-1) printf("\n");
    }
  }
  printf("sumTotDA = %f / %f = %f\n",sumTot[0],sumTot[1],sumTot[0]/sumTot[1]);

  sumTot[0] = 0.; sumTot[1] = 0.;
  for(int iEta=0; iEta<etaBins; iEta++){
    for(int iPt=0; iPt<ptBins; iPt++){
      sumTot[0] = sumTot[0] + numFRBGEL[iEta][iPt];
      sumTot[1] = sumTot[1] + denFRBGEL[iEta][iPt];
      printf("(%d,%d): %9.1f/%9.1f=%4.3f | ",iPt,iEta,numFRBGEL[iEta][iPt],denFRBGEL[iEta][iPt],numFRBGEL[iEta][iPt]/denFRBGEL[iEta][iPt]);
      if(iPt==ptBins-1) printf("\n");
    }
  }
  printf("sumTotBG = %f / %f = %f\n",sumTot[0],sumTot[1],sumTot[0]/sumTot[1]);

  printf("double prompt_rate_DAEL[%d][%d] = {\n",etaBins,ptBins);
  for(int iEta=0; iEta<etaBins; iEta++){
    for(int iPt=0; iPt<ptBins; iPt++){
      printf("%4.3f",numFRDAEL[iEta][iPt]/denFRDAEL[iEta][iPt]);
      if(iPt!=ptBins-1||iEta!=etaBins-1) printf(",");
      if(iPt==ptBins-1) printf("\n");
    }
  }
  printf("};\n");

  printf("double prompt_rate_BGEL[%d][%d] = {\n",etaBins,ptBins);
  for(int iEta=0; iEta<etaBins; iEta++){
    for(int iPt=0; iPt<ptBins; iPt++){
      printf("%4.3f",numFRBGEL[iEta][iPt]/denFRBGEL[iEta][iPt]);
      if(iPt!=ptBins-1||iEta!=etaBins-1) printf(",");
      if(iPt==ptBins-1) printf("\n");
    }
  }
  printf("};\n");
}
