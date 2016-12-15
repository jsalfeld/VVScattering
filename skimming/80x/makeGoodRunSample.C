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

typedef std::map<ULong64_t, std::set<ULong64_t> > EventList;

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
 string jsonFile        = "MitAnalysisRunII/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt"
 ){

  TFile *the_input_file = TFile::Open(input_file.Data());
  TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
  TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");
  TTree *the_SelBit_tree  = (TTree*)the_input_file->FindObjectAny("SelBit_tree");

  BareAll eventAll;
  eventAll.setBranchAddresses(the_input_tree);

  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile(jsonFile.c_str()); 

  ULong64_t N_all  = the_input_tree->GetEntries();
  ULong64_t N_good = 0;

  TNamed *triggerNames = (TNamed*)the_input_file->FindObjectAny("triggerNames");
  printf("triggerNames: %s %s\n",triggerNames->GetName(),triggerNames->GetTitle());

  TFile *outputFile = new TFile(outputFileName.Data(), "RECREATE");
  outputFile->cd();
  TTree *normalizedTree0 = the_input_all ->CloneTree(0);
  TTree *normalizedTree1 = the_input_tree->CloneTree(0);
  TTree *selBitTree;
  if(the_SelBit_tree) selBitTree = the_SelBit_tree->CloneTree(0);
  else                printf("No selBitTree exists\n");

  //dubplicate check
  std::map<ULong64_t, std::set<ULong64_t> > DoubleChecker;
  ULong64_t doubleCount = 0;
  
  for (int i=0; i<the_input_tree->GetEntries(); ++i) {
    the_input_tree->GetEntry(i);
    the_SelBit_tree->GetEntry(i);
    if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

    //------------------------------
    // check if double counting //Remove for MC
    //------------------------------
    Bool_t DuplicateEvent = kFALSE;
    std::map<ULong64_t, std::set<ULong64_t> >::iterator runner = DoubleChecker.find(eventAll.runNum);
    if (runner == DoubleChecker.end()){
      std::set<ULong64_t> evtTemp;
      evtTemp.insert(eventAll.runNum);
      DoubleChecker.insert( make_pair(eventAll.runNum, evtTemp));
    }
    else{
      std::set<ULong64_t>::iterator evter = (*runner).second.find(eventAll.eventNum);
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
    if(the_SelBit_tree) selBitTree->Fill();
  }
  printf("N good/all = %llu / %llu = %f | duplicates: %llu\n",N_good,N_all-doubleCount,(double)N_good/(N_all-doubleCount),doubleCount);
  normalizedTree0->Write();
  normalizedTree1->Write();
  if(the_SelBit_tree) selBitTree->Write();
  triggerNames->Write();
  outputFile->Close();
}
