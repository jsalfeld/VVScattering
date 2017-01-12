#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <iostream>
#include <fstream>

void makePUHist(){

  TString filesPath = "root://eoscms.cern.ch//eos/cms/store/user/ceballos/Nero/v2.1/WGstarToLNuMuMu_012Jets_13TeV-madgraph/WGstarToLNuMuMu_012Jets_13TeV-madgraph/170112_093713/0000/";

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;
  for(int i=10; i<15; i++){
    infilenamev.push_back(Form("%sNeroNtuples_%d.root",filesPath.Data(),i));
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
