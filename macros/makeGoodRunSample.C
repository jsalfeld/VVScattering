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
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

void makeGoodRunSample(
 TString input_file     = "nero_old.root",
 TString outputFileName = "nero_new.root",
 string jsonFile        = "MitAnalysisRunII/json/Cert_246908-251642_13TeV_PromptReco_Collisions15_JSON.txt"
 ){

  TFile *the_input_file = TFile::Open(input_file.Data());
  TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
  TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");

  BareAll eventAll;
  eventAll.setBranchAddresses(the_input_tree);

  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile(jsonFile.c_str()); 

  UInt_t N_all  = the_input_tree->GetEntries();
  UInt_t N_good = 0;

  TNamed *triggerNames = (TNamed*)the_input_file->FindObjectAny("triggerNames");
  printf("triggerNames: %s %s\n",triggerNames->GetName(),triggerNames->GetTitle());

  TFile *outputFile = new TFile(outputFileName.Data(), "RECREATE");
  outputFile->cd();
  TTree *normalizedTree0 = the_input_all ->CloneTree(0);
  TTree *normalizedTree1 = the_input_tree->CloneTree(0);

  //dubplicate check
  std::map<Int_t, std::set<Int_t> > DoubleChecker;
  UInt_t doubleCount = 0;
  
  for (int i=0; i<the_input_tree->GetEntries(); ++i) {
    the_input_tree->GetEntry(i);
    if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

    //------------------------------
    // check if double counting //Remove for MC
    //------------------------------
    Bool_t DuplicateEvent = kFALSE;
    std::map<Int_t, std::set<Int_t> >::iterator runner = DoubleChecker.find(eventAll.runNum);
    if (runner == DoubleChecker.end()){
      std::set<Int_t> evtTemp;
      evtTemp.insert(eventAll.runNum);
      DoubleChecker.insert( make_pair(eventAll.runNum, evtTemp));
    }
    else{
      std::set<Int_t>::iterator evter = (*runner).second.find(eventAll.eventNum);
      if (evter == (*runner).second.end()){
        (*runner).second.insert(eventAll.eventNum);
      }
      else { DuplicateEvent = kTRUE;
      }
    }
    if(DuplicateEvent) doubleCount++;
    if(DuplicateEvent) continue;

    mithep::RunLumiRangeMap::RunLumiPairType rl(eventAll.runNum, eventAll.lumiNum);      

    if(!rlrm.HasRunLumi(rl)) continue;

    N_good++;
    normalizedTree0->Fill(); 
    normalizedTree1->Fill(); 
  }
  printf("N good/all = %d / %d = %f | duplicates: %d\n",N_good,N_all-doubleCount,(double)N_good/(N_all-doubleCount),doubleCount);
  normalizedTree0->Write();
  normalizedTree1->Write();
  triggerNames->Write();
  outputFile->Close();
}
