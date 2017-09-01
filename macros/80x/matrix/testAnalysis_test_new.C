#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
void helper_function(int me=0,int DY=1){//27
  Double_t lumiE=1.025;
  // int gen_bin=5;//11
  int gen_bin=37;
  int reco_bin=2*gen_bin;;
  int bck=3;
  TH1D* px;
  int bin1=1;
  int lep=144;
  int mom=4;
  double prot=0.000001;
  const int nBinRecoPt     = 74; Float_t xbinsRecoPt[nBinRecoPt+1]         = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
									       10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
									       28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
									       220, 235,250, 275,300,325,350,375,400,425,450,475,500, 750,1000};

  /****** input data *******/
  TString path="/afs/cern.ch/user/c/christiw/public/panda/Unfolding/macros/inputs/";
  TFile* file = new TFile(Form("%shistoDY%dzllPtRecGen.root",path.Data(),DY), "read");  if(!file) {printf("File does not exist\n"); return;}

  TH2D* h1=(TH2D*)file->Get(Form("histoPtRecGen_%d",me));
  TH1D* histo_Data=(TH1D*) file->Get(Form("histoPtRecDA_%d",me));
  TH1D* histo_dy=(TH1D*) file->Get(Form("histoPtRecDY_%d",me));

  TH1D* histo_dy_LepEff[lep];
  TH1D* histo_dy_MomRes[mom];
  for(int i=1;i<=mom;i++){
    histo_dy_MomRes[i-1]=(TH1D*) file->Get(Form("histoPtRecDY_MomRes_%d_%d",me,i));
  }
  for(int i=1;i<=lep;i++){
    histo_dy_LepEff[i-1]=(TH1D*) file->Get(Form("histoPtRecDY_LepEff_%d_%d",me,i));
  }


  TH1D* resonant=(TH1D*) file->Get(Form("histoPtRecVV_%d",me));
  TH1D* nonres=(TH1D*) file->Get(Form("histoPtRecEM_%d",me));
 
  TH1D* temp=new TH1D(Form("histo"),  Form("histo"),  nBinRecoPt, xbinsRecoPt);


  TH1D* nonfiducial=(TH1D*) file->Get(Form("histoPtRecDY_%d",me));
  nonfiducial->SetName("histo_nondifucial");


  for(int i=bin1*2-1; i<=reco_bin+2*bin1-2; i++) { // RECO
    double allRecBin = histo_dy->GetBinContent(i);
    double fidRecBin = 0;
    for(int j=bin1; j<=gen_bin+bin1-1; j++) { // GEN
      fidRecBin = fidRecBin + h1->GetBinContent(i,j);
   }
    double gendiff = allRecBin - fidRecBin;
    if(gendiff <= 0) printf("Bin %d totally efficient\n",i);
    nonfiducial->SetBinContent(i-bin1*2+2,TMath::Max(gendiff,prot));
  }

  TH2D* histo_sig_LepEff[lep];                                                                                  
  TH2D* histo_sig_MomRes[mom];
  for(int i=1;i<=mom;i++){
    histo_sig_MomRes[i-1]=(TH2D*)file->Get(Form("histoPtRecGen_MomRes_%d_%d",me,i));;
  }
  for(int i=1;i<=lep;i++){
    histo_sig_LepEff[i-1]=(TH2D*)file->Get(Form("histoPtRecGen_LepEff_%d_%d",me,i));;
  }
                                                                         
  TH1D* histo_resonant_LepEff[lep];
  TH1D* histo_resonant_MomRes[mom];                              
  for(int i=1;i<=mom;i++){
    histo_resonant_MomRes[i-1]=(TH1D*)file->Get(Form("histoPtRecVV_MomRes_%d_%d",me,i));;
  }
  for(int i=1;i<=lep;i++){
    histo_resonant_LepEff[i-1]=(TH1D*)file->Get(Form("histoPtRecVV_LepEff_%d_%d",me,i));;
  }


  TH1D* histo_resonant_PDF=(TH1D*) file->Get(Form("histoPtRecVV_PDF_%d",me));
  TH1D* histo_resonant_QCD=(TH1D*) file->Get(Form("histoPtRecVV_QCD_%d",me));



  /*************new root file*********/


  char outputLimits[200];
  sprintf(outputLimits,"input/data_test_%d%d.root",me,DY);
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  int nBinX=h1->GetNbinsX();
  int nBinY=h1->GetNbinsY();
  Double_t integrals[nBinY];

  /* 37 signals */
  for(int i=bin1; i<=gen_bin+bin1-1; i++){
    px=h1->ProjectionX(Form("histo_s%d",i),i,i);
    for(int j=2*bin1-1; j<=reco_bin+2*bin1-2; j++){
      temp->SetBinContent(j-2*bin1+2,TMath::Max(px->GetBinContent(j),prot));
    }
    px->Reset();
    px->Clear();
    px=(TH1D*)temp->Clone(Form("histo_s%d",i));
    integrals[i-1]=px->Integral();
    px->Write();
  }
  histo_dy->Write();



  /* data */
  histo_Data->SetName("histo_Data");
  histo_Data->Write();
  double integral_data=histo_Data->Integral();
  resonant->SetName("histo_resonant");
  nonres->SetName("histo_nonres");
  nonfiducial->SetName("histo_nonfiducial");
  resonant->Write();
  nonres->Write();
  nonfiducial->Write();

  /****** uncertainties **********/

  //build sig momres and lepeff
  TH1D* MomRes;
  TH1D* LepEff;
  for(int k=1;k<=lep;k++){
    for(int i=bin1; i<=gen_bin+bin1-1; i++){
      LepEff=histo_sig_LepEff[k-1]->ProjectionX(Form("histo_s%d_LepEff%dUp",i,k),i,i);
      LepEff->SetName(Form("histo_s%d_LepEff%dUp",i,k));
      LepEff->Write();
    }
  }  
  for(int k=1;k<=mom;k++){
    for(int i=bin1; i<=gen_bin+bin1-1; i++){
      MomRes=histo_sig_MomRes[k-1]->ProjectionX(Form("histo_s%d_MomRes%dUp",i,k),i,i);
      MomRes->SetName(Form("histo_s%d_MomRes%dUp",i,k));
      MomRes->Write();
    }
  }

  //build nonfiducial momres and lepeff                                                                                                                        
  for(int k=1;k<=lep;k++){
  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++) {
    double error1=TMath::Abs(histo_dy_LepEff[k-1]->GetBinContent(i)-histo_dy->GetBinContent(i));
    double fidRecBin = 0;
    double fidRecAlt=0;
    for(int j=bin1; j<=gen_bin+bin1-1; j++) { // Gen                                                                                                                                     
      fidRecAlt = fidRecBin + histo_sig_LepEff[k-1]->GetBinContent(i,j);
      fidRecBin=fidRecAlt+h1->GetBinContent(i,j);

    }
    double error2=TMath::Abs(fidRecAlt-fidRecBin);
    double error=TMath::Sqrt(error1*error1+error2*error2);
    LepEff->SetBinContent(i,TMath::Max(error+nonfiducial->GetBinContent(i),prot));
  }
  LepEff->SetName(Form("histo_nonfiducial_LepEff%dUp",k));
  LepEff->Write();
  }
  for(int k=1;k<=mom;k++){

  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++) {
    double error1=TMath::Abs(histo_dy_MomRes[k-1]->GetBinContent(i)-histo_dy->GetBinContent(i));
    double fidRecBin = 0;
    double fidRecAlt=0;
    for(int j=bin1; j<=gen_bin+bin1-1; j++) { // Gen
      fidRecAlt = fidRecBin + histo_sig_MomRes[k-1]->GetBinContent(i,j);
      fidRecBin=fidRecAlt+h1->GetBinContent(i,j);
    }
    double error2=TMath::Abs(fidRecAlt-fidRecBin);
    double error=TMath::Sqrt(error1*error1+error2*error2);
    MomRes->SetBinContent(i,TMath::Max(error+nonfiducial->GetBinContent(i),prot));
  }

  MomRes->SetName(Form("histo_nonfiducial_MomRes%dUp",k));
  MomRes->Write();
  }

  //uncertainty for background

  for(int k=1;k<=mom;k++){
    histo_resonant_MomRes[k-1]->SetName(Form("histo_resonant_MomRes%dUp",k));  
    histo_resonant_MomRes[k-1]->Write();
 }
  for(int k=1;k<=lep;k++){
    histo_resonant_LepEff[k-1]->SetName(Form("histo_resonant_LepEff%dUp",k));
    histo_resonant_LepEff[k-1]->Write();
  }

  histo_resonant_PDF->SetName("histo_resonant_PDFUp");
  histo_resonant_QCD->SetName("histo_resonant_QCDUp");
  histo_resonant_PDF->Write();
  histo_resonant_QCD->Write();


  /************ statistical uncertainty ****************/

  TH1D* stat_up=(TH1D*)histo_resonant_QCD->Clone("stat_up");
  TH1D* stat_down=(TH1D*)histo_resonant_QCD->Clone("stat_down");
  TH1D* stat_temp;

  for(int k=bin1; k<=gen_bin+bin1-1; k++){
    stat_temp=(TH1D*) outFileLimits->Get(Form("histo_s%d",k));
    for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
      for(int j=2*bin1-1; j<=reco_bin+2*bin1-2; j++){
	stat_up->SetBinContent(j-bin1*2+2,TMath::Max(stat_temp->GetBinContent(j-bin1*2+2),prot));
	stat_down->SetBinContent(j-bin1*2+2,TMath::Max(stat_temp->GetBinContent(j-bin1*2+2),prot));
      if (i==j){
	stat_up->SetBinContent(j-bin1*2+2,TMath::Max(stat_temp->GetBinContent(j-bin1*2+2)+TMath::Abs(stat_temp->GetBinError(j-bin1*2+2)),prot));
	stat_down->SetBinContent(j-bin1*2+2,TMath::Max(stat_temp->GetBinContent(j-bin1*2+2)-TMath::Abs(stat_temp->GetBinError(j-bin1*2+2)),prot));

      }
    }
    stat_down->SetName(Form("histo_s%d_stat_reco%d_gen%dDown",k,i,k));
    stat_up->SetName(Form("histo_s%d_stat_reco%d_gen%dUp",k,i,k));
    stat_down->Write();
    stat_up->Write();


  }
  }



  outFileLimits->Close();








}
void testAnalysis_test_new(){
  for(int i=0;i<=1;i++){
    for(int j=0;j<=3;j++){
      helper_function(i,j);
    }
  }

}
