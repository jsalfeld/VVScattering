#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1D.h>

void getNPUTrue(TString sourceFileName = "/scratch5/ceballos/ntuples_noweights_76x/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1+AODSIM.root",
                TString outputFileName = "MC_NPUTrue_76x.root") {
  int reBin = 10;

  TFile *denumerator = TFile::Open(sourceFileName.Data());
  TH1D * hDenumerator = (TH1D*)denumerator->Get("hNPUTrue");
  hDenumerator->Sumw2();
  hDenumerator->Rebin(reBin);
  hDenumerator->Scale(1./hDenumerator->Integral());

  TFile output(outputFileName.Data(),"RECREATE");
  hDenumerator->Write();
}
