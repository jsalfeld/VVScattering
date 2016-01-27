#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "DYMVA_WWxsec_2015/dymva_0jet_BDT.class.C"
#include "DYMVA_WWxsec_2015/dymva_1jet_BDT.class.C"
#include "DYMVA_WWxsec_2015/dymva_jets_BDT.class.C"

#include "NeroProducer/Core/interface/BareAll.hpp"
#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BarePhotons.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"
#include "MitAna/Utils/interface/SimpleTable.h"

#include "MitAnalysisRunII/macros/factors.h"
#include "WWAnalysis/resummation/WWpTreweight.h"

// filterType = -1 ==> no filter
//            - 0  ==> ptl1>20 && ptl2>10
//            - 1  ==> QCD enriched sample <--- no PDF info
//            - 2  ==> ptl1>20 && ptl2>10 && min(met,trackMet) > 20, adding DYMVA
//            - 3  ==> ptl1>30, met>30 <--- no PDF info
//            - 4  ==> one photon, ptg>60 <--- no PDF info

void makeSkimSample(
 TString input_file     = "nero_old.root",
 TString outputFileName = "nero_new.root",
 TString processName    = "data",
 Int_t   filterType     = -1,
 Int_t   useFakeMET     = 0
 ){

  bool fillPDFInfo = true;
  if(filterType == 1 || filterType == 3 || filterType == 4) fillPDFInfo = false;

  TFile *the_input_file = TFile::Open(input_file.Data());
  TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
  TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");
  TTree *the_PDF_tree   = (TTree*)the_input_file->FindObjectAny("pdfReweight");
  TH1D* the_pdfReweightSums = (TH1D*)the_input_file->FindObjectAny("pdfReweightSums");

  std::vector<std::string> theInputVars_0j;
  std::vector<std::string> theInputVars_1j;
  std::vector<std::string> theInputVars_2j;
  const char* inputVars_0j[] = { "met", "metsig", "uperp", "upara", "nGoodVertices", "dilep_pt", "min_mt", "max_mt", "min_lep_met_dphi", "max_lep_met_dphi" };
  const char* inputVars_1j[] = { "met", "metsig", "jet1_met_dphi", "upara", "uperp", "nGoodVertices", "dilep_pt", "min_mt", "max_mt", "min_lep_met_dphi", "max_lep_met_dphi", "jet1_pt", "dilep_jet1_dphi" };
  const char* inputVars_2j[] = { "met", "metsig", "min_jet_met_dphi", "max_jet_met_dphi", "upara", "uperp", "nGoodVertices", "dilep_pt", "min_mt", "max_mt", "min_lep_met_dphi", "max_lep_met_dphi", "jet1_pt", "jet2_pt", "dilep_jet1_dphi", "dilep_jet2_dphi", "jet1_jet2_dphi" };
  for (unsigned int i=0;i<10;++i) theInputVars_0j.push_back(inputVars_0j[i]);
  for (unsigned int i=0;i<13;++i) theInputVars_1j.push_back(inputVars_1j[i]);
  for (unsigned int i=0;i<17;++i) theInputVars_2j.push_back(inputVars_2j[i]);
  dymva_0j::ReadBDT rbdtDy_0j = dymva_0j::ReadBDT(theInputVars_0j);
  dymva_1j::ReadBDT rbdtDy_1j = dymva_1j::ReadBDT(theInputVars_1j);
  dymva_2j::ReadBDT rbdtDy_2j = dymva_2j::ReadBDT(theInputVars_2j);

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

  mithep::SimpleTable xstab("$CMSSW_BASE/src/MitAnalysisRunII/data/xs.dat");
  Double_t crossSection = xstab.Get(processName.Data()) * 1000.0; // cross section in fb
  printf("crossSection(%s): %f fb\n",processName.Data(),crossSection);

  BareAll eventAll0;
  eventAll0.setBranchAddresses(the_input_all);

  BareMonteCarlo eventMonteCarlo;
  eventMonteCarlo.setBranchAddresses(the_input_tree);

  BareLeptons eventLeptons;
  eventLeptons.setBranchAddresses(the_input_tree);

  BarePhotons eventPhotons;
  eventPhotons.SetExtend();
  eventPhotons.setBranchAddresses(the_input_tree);

  BareMet eventMet;
  eventMet.SetExtend();
  eventMet.setBranchAddresses(the_input_tree);

  BareTrigger eventTrigger;
  eventTrigger.setBranchAddresses(the_input_tree);

  BareJets eventJets;
  eventJets.setBranchAddresses(the_input_tree);

  BareVertex eventVertex;
  eventVertex.setBranchAddresses(the_input_tree);

  float dymva_= -999.;
  unsigned int nlep_= 0;
  unsigned int njets_= 0;
  if(filterType == 2){
    the_input_tree->Branch( "dymva", &dymva_);     
    the_input_tree->Branch( "nlep", &nlep_);     
    the_input_tree->Branch( "njets", &njets_);     
  }

  TFile *outputFile = new TFile(outputFileName.Data(), "RECREATE");
  outputFile->cd();
  TTree *normalizedTree0   = the_input_all ->CloneTree(0);
  TTree *normalizedTree1   = the_input_tree->CloneTree(0);
  TTree *theClone_PDF_tree;

  if(the_PDF_tree && fillPDFInfo) {
    theClone_PDF_tree = (TTree*)the_PDF_tree->CloneTree(0);
    for (int i=0; i<the_PDF_tree->GetEntries(); ++i) {
      the_PDF_tree->GetEntry(i);
      theClone_PDF_tree->Fill();
    }
  }

  UInt_t N_all  = the_input_all->GetEntries();
  UInt_t N_pass = 0;
  Double_t sumAllEvents = 0;
  Double_t sumPassEvents = 0;
  Double_t sumPassEventsNoNNLO = 0;
  for (int i=0; i<the_input_all->GetEntries(); ++i) {
    the_input_all->GetEntry(i);
    if(i%100000==0) printf("eventsAll %d out of %d\n",i,(int)the_input_all->GetEntries());
    if(processName.CompareTo("data") != 0 && eventAll0.mcWeight == 0) {assert(0); printf("PROBLEM, event weight == 0\n"); return;}
    if(processName.CompareTo("data") != 0) sumAllEvents = sumAllEvents + eventAll0.mcWeight;
    else                                   sumAllEvents = sumAllEvents + 1.0;
    normalizedTree0->Fill(); 
  }

  TNamed *triggerNames = (TNamed*)the_input_file->FindObjectAny("triggerNames");
  printf("triggerNames: %s %s\n",triggerNames->GetName(),triggerNames->GetTitle());
  char **tokens;
  size_t numtokens;
  tokens = strsplit(triggerNames->GetTitle(), ",", &numtokens);
  for (int i = 0; i < (int)numtokens; i++) {
    printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
  }

  for (int i=0; i<the_input_tree->GetEntries(); ++i) {
    the_input_tree->GetEntry(i);
    if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

    vector<int> idOnlyLep;
    for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
      if(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt() <= 10) continue;
      
      if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)     == BareLeptons::LepFake     ||
         ((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepMedium)   == BareLeptons::LepMedium   ||
         ((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepTight)    == BareLeptons::LepTight    ||
         ((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepMediumIP) == BareLeptons::LepMediumIP ||
         ((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepTightIP)  == BareLeptons::LepTightIP) {idOnlyLep.push_back(nlep);}
    }

    vector<int> idLep; vector<int> idTight;
    for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
      if(selectIdIsoCut("default",TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
         TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep]))
        											     {idTight.push_back(1); idLep.push_back(nlep);}
      else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
    }

    Bool_t passFilter = kFALSE;
    if     (filterType == -1) passFilter = kTRUE;
    else if(filterType == 0){
      if(idOnlyLep.size() >= 2 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idOnlyLep[0]])->Pt() > 20 && 
         ((TLorentzVector*)(*eventLeptons.p4)[idOnlyLep[1]])->Pt() > 10) passFilter = kTRUE;
    }
    else if(filterType == 1){
      bool passTrigger = kFALSE;
      for (int nt = 0; nt < (int)numtokens; nt++) {
	if((strcmp(tokens[nt],"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
            strcmp(tokens[nt],"HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
            strcmp(tokens[nt],"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
            strcmp(tokens[nt],"HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
	    strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_v*")  == 0 ||              
	    strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_v*") == 0 ||              
	    strcmp(tokens[nt],"HLT_Mu24_TrkIsoVVL_v*") == 0 ||              
	    strcmp(tokens[nt],"HLT_Mu34_TrkIsoVVL_v*") == 0) &&
	    (*eventTrigger.triggerFired)[nt] == 1) passTrigger = kTRUE;
      }
      if(passTrigger == kTRUE &&
         idOnlyLep.size() >= 1 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idOnlyLep[0]])->Pt() > 10 && 
	 ((TLorentzVector*)(*eventMet.p4)[0])->Pt() < 30.0) passFilter = kTRUE;
    }
    else if(filterType == 2){
      dymva_= -999.;
      nlep_= 0;
      njets_= 0;
      if(idLep.size() >= 2 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 20 && 
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 10 && 
	 TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt()) > 20.0){
	
	passFilter = kTRUE;

        // Begin DYMVA implementation
	TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 
        double the_met = ((TLorentzVector*)(*eventMet.p4)[0])->Pt();
	double the_metsig = ((TLorentzVector*)(*eventMet.p4)[0])->Pt()/sqrt((double)eventMet.sumEtRaw);
	TVector2 metv(((TLorentzVector*)(*eventMet.p4)[0])->Px(), ((TLorentzVector*)(*eventMet.p4)[0])->Py());
        TVector2 dilv(dilep.Px(), dilep.Py());
        TVector2 utv = -1.*(metv+dilv);
        double phiv = utv.DeltaPhi(dilv);
        double the_upara = utv.Mod()*TMath::Cos(phiv);
        double the_uperp = utv.Mod()*TMath::Sin(phiv);
	double the_nGoodVertices =  (double)eventVertex.npv;
	double the_dilep_pt = dilep.Pt();
        double the_min_mt = 999999.;
        double the_max_mt = -1.;
        double the_min_lep_met_dphi = 999.;
        double the_max_lep_met_dphi = -1.;
        for(unsigned nl=0; nl<idLep.size(); nl++){
          double deltaPhiLeptonMet = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
          if(the_min_lep_met_dphi > deltaPhiLeptonMet) the_min_lep_met_dphi = deltaPhiLeptonMet;      
          if(the_max_lep_met_dphi < deltaPhiLeptonMet) the_max_lep_met_dphi = deltaPhiLeptonMet;      

          double mtW = TMath::Sqrt(2.0*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiLeptonMet)));
	  if(the_min_mt > mtW) the_min_mt = mtW;
	  if(the_max_mt < mtW) the_max_mt = mtW;
        }

        vector<int> idJet;
        for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
          bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
          //if(passId == false) continue;

          Bool_t isLepton = kFALSE;
          for(unsigned int nl=0; nl<idLep.size(); nl++){
            if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	  }
	  if(isLepton == kTRUE) continue;

          if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) continue;

          idJet.push_back(nj);
        }

	double the_jet1_met_dphi = -1;
	double the_jet1_pt = -1;
	double the_dilep_jet1_dphi = -1;
	if(idJet.size() >= 1){
	  the_jet1_met_dphi = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
	  the_jet1_pt = ((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Pt();
	  the_dilep_jet1_dphi = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->DeltaPhi(dilep));
	}
	double the_min_jet_met_dphi = 999.;
	double the_max_jet_met_dphi = -1.;
	double the_jet2_pt = -1;
	double the_dilep_jet2_dphi = -1;
	double the_jet1_jet2_dphi = -1.;
	if(idJet.size() >= 2){
	  double dPhiJetMET = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
	  if(the_min_jet_met_dphi > dPhiJetMET) the_min_jet_met_dphi = dPhiJetMET;
	  if(the_max_jet_met_dphi < dPhiJetMET) the_max_jet_met_dphi = dPhiJetMET;
          dPhiJetMET = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[1]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
	  if(the_min_jet_met_dphi > dPhiJetMET) the_min_jet_met_dphi = dPhiJetMET;
	  if(the_max_jet_met_dphi < dPhiJetMET) the_max_jet_met_dphi = dPhiJetMET;
	  the_jet2_pt = ((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Pt();
	  the_dilep_jet2_dphi = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[1]])->DeltaPhi(dilep));
	  the_jet1_jet2_dphi =  TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->DeltaPhi(*((TLorentzVector*)(*eventJets.p4)[idJet[1]])));
        }

        if     (idJet.size() == 0) {
          std::vector<double> theInputVals;
          const double inputVals[] = {the_met, the_metsig, the_uperp, the_upara, the_nGoodVertices, the_dilep_pt, the_min_mt, the_max_mt, the_min_lep_met_dphi, the_max_lep_met_dphi};
          for (int i=0;i<10;++i) theInputVals.push_back(inputVals[i]);
          dymva_ = rbdtDy_0j.GetMvaValue(theInputVals);
        }
	else if(idJet.size() == 1) {
          std::vector<double> theInputVals;
          const double inputVals[] = {the_met, the_metsig, the_jet1_met_dphi, the_upara, the_uperp, the_nGoodVertices, the_dilep_pt, the_min_mt, the_max_mt, the_min_lep_met_dphi, the_max_lep_met_dphi, the_jet1_pt, the_dilep_jet1_dphi};
          for (int i=0;i<13;++i) theInputVals.push_back(inputVals[i]);
          dymva_ = rbdtDy_1j.GetMvaValue(theInputVals);
        }
	else {
          std::vector<double> theInputVals;
          const double inputVals[] = {the_met, the_metsig, the_min_jet_met_dphi, the_max_jet_met_dphi, the_upara, the_uperp, the_nGoodVertices, the_dilep_pt, the_min_mt, the_max_mt, the_min_lep_met_dphi, the_max_lep_met_dphi, the_jet1_pt, the_jet2_pt, the_dilep_jet1_dphi, the_dilep_jet2_dphi, the_jet1_jet2_dphi};
          for (int i=0;i<17;++i) theInputVals.push_back(inputVals[i]);
          dymva_ = rbdtDy_2j.GetMvaValue(theInputVals);
        }
	nlep_ = idLep.size();
	njets_ = idJet.size();
        // End DYMVA implementation

      } // pass filter 2
    } // end filter 2
    else if(filterType == 3){
      if(idLep.size() == 1 &&
        ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 30 && 
	((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 30.0) passFilter = kTRUE;
    } // end filter 3
    else if(filterType == 4){ // begin filter 4
      vector<int> idPho;
      for(int npho=0; npho<eventPhotons.p4->GetEntriesFast(); npho++) {
        if(TMath::Abs(((TLorentzVector*)(*eventPhotons.p4)[npho])->Eta()) >= 2.4) continue;
        if(((TLorentzVector*)(*eventPhotons.p4)[npho])->Pt() <= 60) continue;
        if((double)(*eventPhotons.r9)[npho] > 0.9 &&
	   ((int)(*eventPhotons.selBits)[npho] & BarePhotons::PhoTight)== BarePhotons::PhoTight){idPho.push_back(npho);}
      }
      bool passTrigger = kFALSE;
      for (int nt = 0; nt < (int)numtokens; nt++) {
	if((strcmp(tokens[nt],"HLT_Photon50_R9Id90_HE10_IsoM_v*") == 0 ||
            strcmp(tokens[nt],"HLT_Photon75_R9Id90_HE10_IsoM_v*") == 0 ||
            strcmp(tokens[nt],"HLT_Photon90_R9Id90_HE10_IsoM_v*") == 0 ||
            strcmp(tokens[nt],"HLT_Photon120_R9Id90_HE10_IsoM_v*") == 0) &&
	    (*eventTrigger.triggerFired)[nt] == 1) passTrigger = kTRUE;
      }
      if(passTrigger == kTRUE &&
         idPho.size() >= 1) passFilter = kTRUE;
    } // pass filter 4

    if(passFilter == kFALSE) continue;
    N_pass++;

    if(useFakeMET == 1){
      vector<int>idBoson;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 24) {
	  idBoson.push_back(ngen0);
	}
      }
      if(idBoson.size() == 2){
        TLorentzVector wwSystem(( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(idBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(idBoson[1])) ) )); 
        ((TLorentzVector*)(*eventMet.p4)[0])->SetXYZT(wwSystem.Px(), wwSystem.Py(), wwSystem.Pz(), wwSystem.E());
      }
    }

    if(eventMonteCarlo.mcWeight == 0 && processName.CompareTo("data") != 0) {printf("PROBLEM WIH WEIGTHS\n");return;}

    double thePtwwWeight = 1.0;
    //if(processName.CompareTo("wwlnln") == 0) {
    if(processName.CompareTo("dummy") == 0) {
      vector<int>idBoson;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23||TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 24) {
	  idBoson.push_back(ngen0);
	}
      }
      if(idBoson.size() == 2){
        TLorentzVector wwSystem(( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(idBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(idBoson[1])) ) )); 
        thePtwwWeight = theWWpTreweight.reweight(wwSystem.Pt(), 5); // NNLO wwwights
      }
      if(eventMonteCarlo.mcWeight != 0) sumPassEventsNoNNLO = sumPassEventsNoNNLO + eventMonteCarlo.mcWeight;
      else                              sumPassEventsNoNNLO++;
    }
    
    if(eventMonteCarlo.mcWeight != 0) sumPassEvents = sumPassEvents + eventMonteCarlo.mcWeight * thePtwwWeight;
    else                              sumPassEvents++;

    // weight per event in fb
    if(processName.CompareTo("data") != 0){
      eventMonteCarlo.mcWeight = eventMonteCarlo.mcWeight * crossSection / sumAllEvents;
    }
    else {
      eventMonteCarlo.mcWeight = 1.0;
    }

    normalizedTree1->Fill(); 
  }
  printf("N pass/all = %d / %d = %f | Sum pass/all = %f / %f = %f\n",N_pass,N_all,(double)N_pass/N_all,sumPassEvents,sumAllEvents,sumPassEvents/sumAllEvents);
  if(processName.CompareTo("wwlnln") == 0) printf("sumPassEventsNoNNLO: %f\n",sumPassEventsNoNNLO);
  normalizedTree0->Write();
  normalizedTree1->Write();
  if(the_PDF_tree && fillPDFInfo) theClone_PDF_tree->Write();
  if(the_pdfReweightSums && fillPDFInfo) the_pdfReweightSums->Write();
  triggerNames->Write();
  outputFile->Close();
}
