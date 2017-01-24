#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <iostream>
#include <fstream>

void readRootFiles(){

  TString filesPath = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/skim_80x/DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/170120_103436/0000/";

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;
  for(int i=0; i<231; i++){
    infilenamev.push_back(Form("%sNeroNtuples_%d.root",filesPath.Data(),i));
  }

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());
    TFile *the_input_file = TFile::Open(infilenamev[ifile].Data());
    if(the_input_file){
      the_input_file->Map();
    }
  }
}
