#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <iostream>
#include <fstream>

void makePUHist(){

  TString filesPath0 = "root://eoscms.cern.ch//eos/cms/store/user/ceballos/Nero/v2.1/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/170112_100003/0000/";
  TString filesPath1 = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/v2.1/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/170112_143437/0000/";
  TString filesPath2 = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/v2.1/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/TTToSemilepton_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/170112_143641/0000/";

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;
  for(int i=0; i<500; i++){
    infilenamev.push_back(Form("%sNeroNtuples_%d.root",filesPath0.Data(),i));
  }
  for(int i=0; i<500; i++){
    infilenamev.push_back(Form("%sNeroNtuples_%d.root",filesPath1.Data(),i));
  }
  for(int i=0; i<500; i++){
    infilenamev.push_back(Form("%sNeroNtuples_%d.root",filesPath2.Data(),i));
  }

  int nBinPlot      = 100;
  double xminPlot   = 0.0;
  double xmaxPlot   = 100.0;
  TH1D* histos = new TH1D("pileup", "pileup", nBinPlot, xminPlot, xmaxPlot);

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());
    TFile *the_input_file = TFile::Open(infilenamev[ifile].Data());
    if(the_input_file){
      TH1D *histoAux = (TH1D*)the_input_file->FindObjectAny("hDPileup");
      histos->Add(histoAux);
    }
  }

  char output[200];
  sprintf(output,"MC_NPUTrue_80x_Moriond17.root");     
  TFile* outFilePlotsNote = new TFile(output,"recreate");
  outFilePlotsNote->cd();
  histos->Write();
  outFilePlotsNote->Close();
}
