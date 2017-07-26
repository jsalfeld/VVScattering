void mergeZptHist() {
TString inputFolder = "/data/t3home000/ceballos/panda/v_004_0";

TFile *_file[6];
_file[0] = TFile::Open(Form("%s/DYJetsToLL_Pt0To50.root"  ,inputFolder.Data()));
_file[1] = TFile::Open(Form("%s/DYJetsToLL_Pt50To100.root",inputFolder.Data()));
_file[2] = TFile::Open(Form("%s/DYJetsToLL_Pt100To250.root",inputFolder.Data()));
_file[3] = TFile::Open(Form("%s/DYJetsToLL_Pt250To400.root",inputFolder.Data()));
_file[4] = TFile::Open(Form("%s/DYJetsToLL_Pt400To650.root",inputFolder.Data()));
_file[5] = TFile::Open(Form("%s/DYJetsToLL_Pt650ToInf.root",inputFolder.Data()));

double xs[6] = {5387.789,357.200,81.2000,3.06200,0.39050,0.03613};

TH1D *hDITotalMCWeight[6];

TH1D *hDIDilPtMM[6];
TH1D *hDIDilPtEE[6];
TH1D *hDIDilPt2MM[6];
TH1D *hDIDilPt2EE[6];
TH1D *hDIDilRapMM[6];
TH1D *hDIDilRapEE[6];
TH1D *hDIDilPtRap0MM[6];
TH1D *hDIDilPtRap0EE[6];
TH1D *hDIDilPtRap1MM[6];
TH1D *hDIDilPtRap1EE[6];
TH1D *hDIDilPtRap2MM[6];
TH1D *hDIDilPtRap2EE[6];
TH1D *hDIDilPtRap3MM[6];
TH1D *hDIDilPtRap3EE[6];
TH1D *hDIDilPtRap4MM[6];
TH1D *hDIDilPtRap4EE[6];

for(int i=0; i<6; i++){
  hDITotalMCWeight[i] = (TH1D*)_file[i]->Get("hDTotalMCWeight");	 
  hDIDilPtMM[i]       = (TH1D*)_file[i]->Get("hDDilPtMM");     hDIDilPtMM[i]    ->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());	 
  hDIDilPtEE[i]       = (TH1D*)_file[i]->Get("hDDilPtEE");     hDIDilPtEE[i]    ->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());	 
  hDIDilPt2MM[i]      = (TH1D*)_file[i]->Get("hDDilPt2MM");    hDIDilPt2MM[i]   ->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());	
  hDIDilPt2EE[i]      = (TH1D*)_file[i]->Get("hDDilPt2EE");    hDIDilPt2EE[i]   ->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());	
  hDIDilRapMM[i]      = (TH1D*)_file[i]->Get("hDDilRapMM");    hDIDilRapMM[i]   ->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());	
  hDIDilRapEE[i]      = (TH1D*)_file[i]->Get("hDDilRapEE");    hDIDilRapEE[i]   ->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());	
  hDIDilPtRap0MM[i]   = (TH1D*)_file[i]->Get("hDDilPtRap0MM"); hDIDilPtRap0MM[i]->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());
  hDIDilPtRap0EE[i]   = (TH1D*)_file[i]->Get("hDDilPtRap0EE"); hDIDilPtRap0EE[i]->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());
  hDIDilPtRap1MM[i]   = (TH1D*)_file[i]->Get("hDDilPtRap1MM"); hDIDilPtRap1MM[i]->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());
  hDIDilPtRap1EE[i]   = (TH1D*)_file[i]->Get("hDDilPtRap1EE"); hDIDilPtRap1EE[i]->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());
  hDIDilPtRap2MM[i]   = (TH1D*)_file[i]->Get("hDDilPtRap2MM"); hDIDilPtRap2MM[i]->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());
  hDIDilPtRap2EE[i]   = (TH1D*)_file[i]->Get("hDDilPtRap2EE"); hDIDilPtRap2EE[i]->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());
  hDIDilPtRap3MM[i]   = (TH1D*)_file[i]->Get("hDDilPtRap3MM"); hDIDilPtRap3MM[i]->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());
  hDIDilPtRap3EE[i]   = (TH1D*)_file[i]->Get("hDDilPtRap3EE"); hDIDilPtRap3EE[i]->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());
  hDIDilPtRap4MM[i]   = (TH1D*)_file[i]->Get("hDDilPtRap4MM"); hDIDilPtRap4MM[i]->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());
  hDIDilPtRap4EE[i]   = (TH1D*)_file[i]->Get("hDDilPtRap4EE"); hDIDilPtRap4EE[i]->Scale(xs[i]/hDITotalMCWeight[i]->GetSumOfWeights());
}

TH1D *hDDilPtMM     = (TH1D*)hDIDilPtMM    [0]->Clone();
TH1D *hDDilPtEE     = (TH1D*)hDIDilPtEE    [0]->Clone();
TH1D *hDDilPt2MM    = (TH1D*)hDIDilPt2MM   [0]->Clone();
TH1D *hDDilPt2EE    = (TH1D*)hDIDilPt2EE   [0]->Clone();
TH1D *hDDilRapMM    = (TH1D*)hDIDilRapMM   [0]->Clone();
TH1D *hDDilRapEE    = (TH1D*)hDIDilRapEE   [0]->Clone();
TH1D *hDDilPtRap0MM = (TH1D*)hDIDilPtRap0MM[0]->Clone();
TH1D *hDDilPtRap0EE = (TH1D*)hDIDilPtRap0EE[0]->Clone();
TH1D *hDDilPtRap1MM = (TH1D*)hDIDilPtRap1MM[0]->Clone();
TH1D *hDDilPtRap1EE = (TH1D*)hDIDilPtRap1EE[0]->Clone();
TH1D *hDDilPtRap2MM = (TH1D*)hDIDilPtRap2MM[0]->Clone();
TH1D *hDDilPtRap2EE = (TH1D*)hDIDilPtRap2EE[0]->Clone();
TH1D *hDDilPtRap3MM = (TH1D*)hDIDilPtRap3MM[0]->Clone();
TH1D *hDDilPtRap3EE = (TH1D*)hDIDilPtRap3EE[0]->Clone();
TH1D *hDDilPtRap4MM = (TH1D*)hDIDilPtRap4MM[0]->Clone();
TH1D *hDDilPtRap4EE = (TH1D*)hDIDilPtRap4EE[0]->Clone();

for(int i=1; i<6; i++){
  hDDilPtMM	->Add(hDIDilPtMM    [i]);   
  hDDilPtEE	->Add(hDIDilPtEE    [i]);   
  hDDilPt2MM	->Add(hDIDilPt2MM   [i]);  
  hDDilPt2EE	->Add(hDIDilPt2EE   [i]);  
  hDDilRapMM	->Add(hDIDilRapMM   [i]);  
  hDDilRapEE	->Add(hDIDilRapEE   [i]);  
  hDDilPtRap0MM ->Add(hDIDilPtRap0MM[i]);
  hDDilPtRap0EE ->Add(hDIDilPtRap0EE[i]);
  hDDilPtRap1MM ->Add(hDIDilPtRap1MM[i]);
  hDDilPtRap1EE ->Add(hDIDilPtRap1EE[i]);
  hDDilPtRap2MM ->Add(hDIDilPtRap2MM[i]);
  hDDilPtRap2EE ->Add(hDIDilPtRap2EE[i]);
  hDDilPtRap3MM ->Add(hDIDilPtRap3MM[i]);
  hDDilPtRap3EE ->Add(hDIDilPtRap3EE[i]);
  hDDilPtRap4MM ->Add(hDIDilPtRap4MM[i]);
  hDDilPtRap4EE ->Add(hDIDilPtRap4EE[i]);
}

TFile myOutputFile("genZpt.root","RECREATE");
  hDDilPtMM	->Write(); 
  hDDilPtEE	->Write(); 
  hDDilPt2MM	->Write();
  hDDilPt2EE	->Write();
  hDDilRapMM	->Write();
  hDDilRapEE	->Write();
  hDDilPtRap0MM ->Write();
  hDDilPtRap0EE ->Write();
  hDDilPtRap1MM ->Write();
  hDDilPtRap1EE ->Write();
  hDDilPtRap2MM ->Write();
  hDDilPtRap2EE ->Write();
  hDDilPtRap3MM ->Write();
  hDDilPtRap3EE ->Write();
  hDDilPtRap4MM ->Write();
  hDDilPtRap4EE ->Write();
myOutputFile.Close();

}
