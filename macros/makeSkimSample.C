#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "NeroProducer/Core/interface/BareAll.hpp"
#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "MitAna/Utils/interface/SimpleTable.h"

#include "MitAnalysisRunII/macros/factors.h"

// filterType = -1 ==> no filter
//            - 0  ==> ptl1>20 && ptl2>10
//            - 1  ==> QCD enriched sample
//            - 2  ==> ptl1>20 && ptl2>10 && min(met,trackMet) > 20

void makeSkimSample(
 TString input_file     = "nero_old.root",
 TString outputFileName = "nero_new.root",
 TString processName    = "data",
 Int_t   filterType     = -1
 ){

  TFile *the_input_file = TFile::Open(input_file.Data());
  TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
  TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");

  mithep::SimpleTable xstab("$CMSSW_BASE/src/MitAnalysisRunII/data/xs.dat");
  Double_t crossSection = xstab.Get(processName.Data()) * 1000.0; // cross section in fb
  printf("crossSection(%s): %f fb\n",processName.Data(),crossSection);

  BareAll eventAll0;
  eventAll0.setBranchAddresses(the_input_all);

  BareMonteCarlo eventMonteCarlo;
  eventMonteCarlo.setBranchAddresses(the_input_tree);

  BareLeptons eventLeptons;
  eventLeptons.setBranchAddresses(the_input_tree);

  BareMet eventMet;
  eventMet.SetExtend();
  eventMet.setBranchAddresses(the_input_tree);

  BareTrigger eventTrigger;
  eventTrigger.setBranchAddresses(the_input_tree);

  TFile *outputFile = new TFile(outputFileName.Data(), "RECREATE");
  outputFile->cd();
  TTree *normalizedTree0 = the_input_all ->CloneTree(0);
  TTree *normalizedTree1 = the_input_tree->CloneTree(0);

  UInt_t N_all  = the_input_all->GetEntries();
  UInt_t N_pass = 0;
  Double_t sumAllEvents = 0;
  Double_t sumPassEvents = 0;
  for (int i=0; i<the_input_all->GetEntries(); ++i) {
    the_input_all->GetEntry(i);
    if(i%100000==0) printf("eventsAll %d out of %d\n",i,(int)the_input_all->GetEntries());
    if(eventAll0.mcWeight == 0) assert(0);
    if(processName.CompareTo("data") != 0) sumAllEvents = sumAllEvents + eventAll0.mcWeight/TMath::Abs(eventAll0.mcWeight);
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

    Bool_t passFilter = kFALSE;
    if     (filterType == -1) passFilter = kTRUE;
    else if(filterType == 0){
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
         ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
         ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter = kTRUE;
    }
    else if(filterType == 1){
      bool passTrigger = kFALSE;
      for (int nt = 0; nt < (int)numtokens; nt++) {
	if((strcmp(tokens[nt],"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
            strcmp(tokens[nt],"HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
            strcmp(tokens[nt],"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
            strcmp(tokens[nt],"HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
            strcmp(tokens[nt],"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0         ||
            strcmp(tokens[nt],"HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*") == 0         ||
	    strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_v*")  == 0 ||              
	    strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_v*") == 0 ||              
	    strcmp(tokens[nt],"HLT_Mu24_TrkIsoVVL_v*") == 0 ||              
	    strcmp(tokens[nt],"HLT_Mu34_TrkIsoVVL_v*") == 0) &&
	    (*eventTrigger.triggerFired)[nt] == 1) passTrigger = kTRUE;
      }
      if(passTrigger == kTRUE &&
         eventLeptons.p4->GetEntriesFast() >= 1 &&
         ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 10 && 
	 (double)eventMet.pfMet_e3p0->Pt() < 30.0) passFilter = kTRUE;
    }
    else if(filterType == 2){
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
         ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
         ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10 && 
	 TMath::Min((double)eventMet.pfMet_e3p0->Pt(),(double)eventMet.trackMet->Pt()) > 20.0) passFilter = kTRUE;
    }

    if(passFilter == kFALSE) continue;

    N_pass++;
    sumPassEvents = sumPassEvents + eventMonteCarlo.mcWeight/TMath::Abs(eventMonteCarlo.mcWeight);

    // weight per event in fb
    if(processName.CompareTo("data") != 0){
      eventMonteCarlo.mcWeight = eventMonteCarlo.mcWeight/TMath::Abs(eventMonteCarlo.mcWeight) * crossSection / sumAllEvents;
    }

    normalizedTree1->Fill(); 
  }
  printf("N pass/all = %d / %d = %f | Sum pass/all = %f / %f = %f\n",N_pass,N_all,(double)N_pass/N_all,sumPassEvents,sumAllEvents,sumPassEvents/sumAllEvents);
  normalizedTree0->Write();
  normalizedTree1->Write();
  triggerNames->Write();
  outputFile->Close();
}
