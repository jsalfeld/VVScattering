#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <set>
#include <fstream>
#include <stdexcept>
#include <string>
#include <cstdlib>
#include <iostream>

#include "NeroProducer/Core/interface/BareAll.hpp"
#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

typedef std::map<Int_t, std::set<Int_t> > EventList;

EventList
readEventList(char const* _fileName)
{
  EventList list;
  ifstream listFile(_fileName);
  if (!listFile.is_open())
    throw std::runtime_error(_fileName);

  unsigned iL(0);
  std::string line;
  while (true) {
    std::getline(listFile, line);
    if (!listFile.good())
      break;
    
    if (line.find(":") == std::string::npos || line.find(":") == line.rfind(":"))
      continue;

    unsigned run(std::atoi(line.substr(0, line.find(":")).c_str()));
    unsigned event(std::atoi(line.substr(line.rfind(":") + 1).c_str()));

    list[run].insert(event);

    ++iL;
  }

  std::cout << "Loaded " << iL << " events" << std::endl;

  return list;
}

void makeGoodRunSample(
 TString input_file     = "nero_old.root",
 TString outputFileName = "nero_new.root",
 string jsonFile        = "VVScattering/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt",
 string badFile         = "json/csc2015_Dec01.txt"
 ){

  EventList theBadEventList = readEventList(badFile.c_str());

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
  UInt_t badCount = 0;
  
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

    Bool_t BadEvent = kFALSE;
    std::map<Int_t, std::set<Int_t> >::iterator rbadItr = theBadEventList.find(eventAll.runNum);
    if (rbadItr != theBadEventList.end()) {
      std::set<Int_t>::iterator evbadter = (*rbadItr).second.find(eventAll.eventNum);
      if (evbadter != (*rbadItr).second.end()) BadEvent = kTRUE;
    }
    if(BadEvent) badCount++;
    if(BadEvent) continue;

    N_good++;
    normalizedTree0->Fill(); 
    normalizedTree1->Fill(); 
  }
  printf("N good/all = %d / %d = %f | duplicates: %d | bad: %d\n",N_good,N_all-doubleCount,(double)N_good/(N_all-doubleCount),doubleCount,badCount);
  normalizedTree0->Write();
  normalizedTree1->Write();
  triggerNames->Write();
  outputFile->Close();
}
