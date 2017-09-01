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

void testAnalysis_test200(int me=0,int DY=3,int alt=0,int bin1=27,int gen_bin=11){//27
  Double_t lumiE=1.025;
  // int gen_bin=5;//11
  int reco_bin=2*gen_bin;;
  int bck=3;
  TH1D* temp=new TH1D("histo","histo",reco_bin,0,gen_bin);
  TH1D* px;
  int lep=144;
  int mom=4;
  /****** input data *******/
  TString path="/afs/cern.ch/user/c/christiw/public/panda/Unfolding/macros/inputs/";
  TFile* file = new TFile(Form("%shistoDY%dzllPtRecGen.root",path.Data(),DY), "read");  if(!file) {printf("File does not exist\n"); return;}
  TFile* file_alt = new TFile(Form("%shistoDY%dzllPtRecGen.root",path.Data(),alt), "read");  if(!file) {printf("File Alternative does not exist\n"); return;}
  TH2D* h1_alt=(TH2D*)file_alt->Get(Form("histoPtRecGen_%d",me));

  TH2D* h1=(TH2D*)file->Get(Form("histoPtRecGen_%d",me));
  TH1D* histo_Data=(TH1D*) file->Get(Form("histoPtRecDA_%d",me));
  TH1D* histo_dy=(TH1D*) file->Get(Form("histoPtRecDY_%d",me));
  TH1D* histo_dy_alt=(TH1D*) file_alt->Get(Form("histoPtRecDY_%d",me));

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
  TH1D* resonant_alt=(TH1D*) file_alt->Get(Form("histoPtRecVV_%d",me));
  TH1D* nonres_alt=(TH1D*) file_alt->Get(Form("histoPtRecEM_%d",me));

  TH1D* nonfiducial=new TH1D("histo_nonfiducial","histo_nonfiducial",reco_bin,0,gen_bin);
  for(int i=bin1*2-1; i<=reco_bin+2*bin1-2; i++) { // RECO
    double allRecBin = histo_dy->GetBinContent(i);
    double fidRecBin = 0;
    for(int j=bin1; j<=gen_bin+bin1-1; j++) { // GEN
      fidRecBin = fidRecBin + h1->GetBinContent(i,j);
   }
    double gendiff = allRecBin - fidRecBin;
    if(gendiff <= 0) printf("Bin %d totally efficient\n",i);
    nonfiducial->SetBinContent(i-bin1*2+2,TMath::Max(gendiff,0.0000001));
  }

  TH1D* nonfiducial_alt=new TH1D("histo_nonfiducial_mcDown","histo_nonfiducial_mcDown",reco_bin,0,gen_bin);
  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++) { // RECO                                                                                                             
    double allRecBin = histo_dy_alt->GetBinContent(i);
    double fidRecBin = 0;
    for(int j=bin1; j<=gen_bin+bin1-1; j++) { // GEN                                                                                         
      fidRecBin = fidRecBin + h1_alt->GetBinContent(i,j);
    }
    double gendiff = allRecBin - fidRecBin;
    if(gendiff <= 0) printf("Bin %d totally efficient\n",i);
    nonfiducial_alt->SetBinContent(i-2*bin1+2,TMath::Max(gendiff,0.0000001));
  }
  TH2D* histo_sig_LepEff[lep];
  TH2D* histo_sig_MomRes[mom];
  for(int i=1;i<=mom;i++){
    histo_sig_MomRes[i-1]=(TH2D*)file->Get(Form("histoPtRecGen_MomRes_%d_%d",me,i));;
  }
  for(int i=1;i<=lep;i++){
    histo_sig_LepEff[i-1]=(TH2D*)file->Get(Form("histoPtRecGen_LepEff_%d_%d",me,i));;
  }  TH1D* histo_nonres_MomRes=(TH1D*) file->Get(Form("histoPtRecEM_MomRes_%d",me));
  TH1D* histo_resonant_PDF=(TH1D*) file->Get(Form("histoPtRecVV_PDF_%d",me));
  TH1D* histo_resonant_QCD=(TH1D*) file->Get(Form("histoPtRecVV_QCD_%d",me));

  TH1D* histo_resonant_LepEff[lep];
  TH1D* histo_resonant_MomRes[mom];
  for(int i=1;i<=mom;i++){
    histo_resonant_MomRes[i-1]=(TH1D*)file->Get(Form("histoPtRecVV_MomRes_%d_%d",me,i));;
  }
  for(int i=1;i<=lep;i++){
    histo_resonant_LepEff[i-1]=(TH1D*)file->Get(Form("histoPtRecVV_LepEff_%d_%d",me,i));;
  }


  /********** cut to reco # of bins ************/

  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    temp->SetBinContent(i-2*bin1+2,TMath::Max(histo_Data->GetBinContent(i),0.0000001));
  }
  histo_Data->Reset();
  histo_Data->Clear();
  histo_Data=(TH1D*)temp->Clone("histo_Data");

  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    temp->SetBinContent(i-2*bin1+2,TMath::Max(nonres->GetBinContent(i),0.0000001));
  }
  nonres->Reset();
  nonres->Clear();
  nonres=(TH1D*)temp->Clone("histo_nonres");

  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    temp->SetBinContent(i-2*bin1+2,TMath::Max(resonant->GetBinContent(i),0.000001));
  }
  resonant->Reset();
  resonant->Clear();
  resonant=(TH1D*)temp->Clone("histo_resonant");

  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    temp->SetBinContent(i-2*bin1+2,TMath::Max(histo_resonant_PDF->GetBinContent(i),0.0000001));
  }
  histo_resonant_PDF->Reset();
  histo_resonant_PDF->Clear();
  histo_resonant_PDF=(TH1D*)temp->Clone("histo_resonant_PDF");

  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    temp->SetBinContent(i-bin1*2+2,TMath::Max(histo_resonant_QCD->GetBinContent(i),0.0000001));
  }
  histo_resonant_QCD->Reset();
  histo_resonant_QCD->Clear();
  histo_resonant_QCD=(TH1D*)temp->Clone("histo_resonant_QCD");


  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    temp->SetBinContent(i-bin1*2+2,TMath::Max(nonres_alt->GetBinContent(i),0.0000001));
  }
  nonres_alt->Reset();
  nonres_alt->Clear();
  nonres_alt=(TH1D*)temp->Clone("histo_nonres_mcUp");

  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    temp->SetBinContent(i-bin1*2+2,TMath::Max(resonant_alt->GetBinContent(i),0.00000001));
  }
  resonant_alt->Reset();
  resonant_alt->Clear();
  resonant_alt=(TH1D*)temp->Clone("histo_resonant_mcUp");
  for(int j=0;j<mom;j++){
  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    temp->SetBinContent(i-bin1*2+2,TMath::Max(histo_resonant_MomRes[j]->GetBinContent(i),0.0000001));
  }
  histo_resonant_MomRes[j]->Reset();
  histo_resonant_MomRes[j]->Clear();
  histo_resonant_MomRes[j]=(TH1D*)temp->Clone(Form("histo_resonant_MomRes%dUp",j+1));
  }
 for(int j=0;j<lep;j++){
  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    temp->SetBinContent(i-bin1*2+2,TMath::Max(histo_resonant_LepEff[j]->GetBinContent(i),0.0000001));
  }
  histo_resonant_LepEff[j]->Reset();
  histo_resonant_LepEff[j]->Clear();
  histo_resonant_LepEff[j]=(TH1D*)temp->Clone(Form("histo_resonant_LepEff%dUp",j+1));
  }
 
  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    temp->SetBinContent(i-bin1*2+2,TMath::Max(histo_dy->GetBinContent(i),0.0000001));
  }
  histo_dy->Reset();
  histo_dy->Clear();
  histo_dy=(TH1D*)temp->Clone("histo_dy");



  /*************new root file*********/


  char outputLimits[200];
  sprintf(outputLimits,"data_test_%d%d%d.root",me,DY,alt);
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  int nBinX=h1->GetNbinsX();
  int nBinY=h1->GetNbinsY();
  Double_t integrals[nBinY];

  /* 37 signals */
  for(int i=bin1; i<=gen_bin+bin1-1; i++){
    px=h1->ProjectionX(Form("histo_s%d",i),i,i);
    for(int j=2*bin1-1; j<=reco_bin+2*bin1-2; j++){
      temp->SetBinContent(j-2*bin1+2,TMath::Max(px->GetBinContent(j),0.0000001));
    }
    px->Reset();
    px->Clear();
    px=(TH1D*)temp->Clone(Form("histo_s%d",i));
    integrals[i-1]=px->Integral();
    px->Write();
  }
  histo_dy->Write();

  /*********** MC uncertainty **************/
  for(int i=bin1; i<=gen_bin+bin1-1; i++){
    px=h1_alt->ProjectionX(Form("histo_s%d_mcDown",i),i,i);
    for(int j=2*bin1-1; j<=reco_bin+2*bin1-2; j++){
      temp->SetBinContent(j-2*bin1+2,TMath::Max(px->GetBinContent(j),0.0000001));
    }
    px->Reset();
    px->Clear();
    px=(TH1D*)temp->Clone(Form("histo_s%d_mcDown",i));
    px->Write();
  }


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
  nonfiducial_alt->Write();
  /****** uncertainties **********/

  //build sig momres and lepeff
  TH1D* MomRes;
  TH1D* LepEff;
  for(int k=1;k<=lep;k++){
    for(int i=bin1; i<=gen_bin+bin1-1; i++){
      LepEff=histo_sig_LepEff[k-1]->ProjectionX(Form("histo_s%d_LepEff%dUp",i,k),i,i);
      for(int j=2*bin1-1; j<=reco_bin+2*bin1-2; j++){
	temp->SetBinContent(j-2*bin1+2,TMath::Max(LepEff->GetBinContent(j),0.0000001));
      }
      LepEff->Reset();LepEff->Clear();
      LepEff=(TH1D*)temp->Clone(Form("histo_s%d_LepEff%dUp",i,k));
      LepEff->Write();
    }
  }
  for(int k=1;k<=mom;k++){
    for(int i=bin1; i<=gen_bin+bin1-1; i++){
      MomRes=histo_sig_MomRes[k-1]->ProjectionX(Form("histo_s%d_MomRes%dUp",i,k),i,i);
      for(int j=2*bin1-1; j<=reco_bin+2*bin1-2; j++){
        temp->SetBinContent(j-2*bin1+2,TMath::Max(MomRes->GetBinContent(j),0.0000001));
      } 
      MomRes->Reset();MomRes->Clear();
      MomRes=(TH1D*)temp->Clone(Form("histo_s%d_MomRes%dUp",i,k));
      MomRes->Write();
    }
  }  //build nonfiducial momres and lepeff                        
  //uncertainty for bck
  histo_resonant_PDF->SetName("histo_resonant_PDFUp");
  histo_resonant_QCD->SetName("histo_resonant_QCDUp");
  for(int k=1;k<=mom;k++){
    histo_resonant_MomRes[k-1]->SetName(Form("histo_resonant_MomRes%dUp",k));
    histo_resonant_MomRes[k-1]->Write();
  }
  for(int k=1;k<=144;k++){
    histo_resonant_LepEff[k-1]->SetName(Form("histo_resonant_LepEff%dUp",k));
    histo_resonant_LepEff[k-1]->Write();
  }

  histo_resonant_PDF->Write();
  histo_resonant_QCD->Write();


  /******** produce uncertainty down histo ***********/
  double mean;
  double up;
  double diff;
  double down;
  TH1D* histo_resonant_QCD_down=(TH1D*)histo_resonant_QCD->Clone("histo_resonant_QCDDown");
  TH1D* histo_resonant_PDF_down=(TH1D*)histo_resonant_PDF->Clone("histo_resonant_PDFDown");
  //  histo_resonant_PDF_down->Reset();histo_resonant_PDF_down->Clear();
  //histo_resonant_QCD_down->Reset();histo_resonant_QCD_down->Clear();

  
  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    mean = resonant->GetBinContent(i-2*bin1+2);
    up=histo_resonant_PDF->GetBinContent(i-2*bin1+2);
    diff=up-mean;
    histo_resonant_PDF_down->SetBinContent(i-2*bin1+2,TMath::Max(mean-diff,0.000001));
    up=histo_resonant_QCD->GetBinContent(i-2*bin1+2);
    diff=up-mean;
    histo_resonant_QCD_down->SetBinContent(i-bin1*2+2,TMath::Max(mean-diff,0.000001));

  }
  histo_resonant_PDF_down->Write();
  histo_resonant_QCD_down->Write();

  TH1D* histo_resonant_LepEff_down=(TH1D*)histo_resonant_LepEff[0]->Clone("histo_resonant_LepEffDown");
  TH1D* histo_resonant_MomRes_down=(TH1D*)histo_resonant_MomRes[0]->Clone("histo_resonant_MomResDown");  
  for(int j=0;j<lep;j++){
    for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
      mean = resonant->GetBinContent(i-2*bin1+2);
      up=histo_resonant_LepEff[j]->GetBinContent(i-2*bin1+2);
      diff=up-mean;
      histo_resonant_LepEff_down->SetBinContent(i-bin1*2+2,TMath::Max(mean-diff,0.000001));
    }
    histo_resonant_LepEff_down->SetName(Form("histo_resonant_LepEff%dDown",j+1));
    histo_resonant_LepEff_down->Write();

  }
  for(int j=0;j<mom;j++){
    for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
      mean = resonant->GetBinContent(i-2*bin1+2);
      up=histo_resonant_MomRes[j]->GetBinContent(i-2*bin1+2);
      diff=up-mean;
      histo_resonant_MomRes_down->SetBinContent(i-bin1*2+2,TMath::Max(mean-diff,0.000001));
     
    }
    histo_resonant_MomRes_down->SetName(Form("histo_resonant_MomRes%dDown",j+1));
    histo_resonant_MomRes_down->Write();

  }


  /********* nondifucial LepEff and MomRes ********/
  /*  TH1D* histo_nonfiducial_LepEff_down=(TH1D*)histo_resonant_LepEff->Clone("histo_nonfiducial_LepEffDown");
  TH1D* histo_nonfiducial_MomRes_down=(TH1D*)histo_resonant_MomRes->Clone("histo_nonfiducial_MomResDown");
  TH1D* histo_nonfiducial_LepEff_up=(TH1D*) outFileLimits->Get("histo_nonfiducial_LepEffUp");
  TH1D* histo_nonfiducial_MomRes_up=(TH1D*) outFileLimits->Get("histo_nonfiducial_MomResUp");


  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    mean = nonfiducial->GetBinContent(i-2*bin1+2);
    up=histo_nonfiducial_MomRes_up->GetBinContent(i-2*bin1+2);
    diff=up-mean;
    histo_nonfiducial_MomRes_down->SetBinContent(i-bin1*2+2,TMath::Max(mean-diff,0.000001));
    up=histo_nonfiducial_LepEff_up->GetBinContent(i-2*bin1+2);
    diff=up-mean;
    histo_nonfiducial_LepEff_down->SetBinContent(i-bin1*2+2,TMath::Max(mean-diff,0.000001));

  }
  histo_nonfiducial_MomRes_down->Write();
  histo_nonfiducial_LepEff_down->Write();*/
  /************ MC up uncertainty ***************/

  TH1D* histo_mean;
  TH1D* histo_mc_down;
  TH1D* histo_mc_up;
  for(int i=bin1; i<=gen_bin+bin1-1; i++){
    
    histo_mean=(TH1D*) outFileLimits->Get(Form("histo_s%d",i));
    histo_mc_down=(TH1D*) outFileLimits->Get(Form("histo_s%d_mcDown",i));
    histo_mc_up=(TH1D*)histo_mean->Clone(Form("histo_s%d_mcUp",i));
    for(int j=2*bin1-1; j<=reco_bin+2*bin1-2; j++){
      mean=histo_mean->GetBinContent(j-2*bin1+2);
      down=histo_mc_down->GetBinContent(j-2*bin1+2);
      diff=mean-down;
      histo_mc_up->SetBinContent(j-bin1*2+2,TMath::Max(mean+diff,0.000001));
    }
    histo_mc_up->Write();


  }

  /******** MC nonfiducial up uncertainty ***********/

  for(int j=2*bin1-1; j<=reco_bin+2*bin1-2; j++){
    mean=nonfiducial->GetBinContent(j-2*bin1+2);
    down=nonfiducial_alt->GetBinContent(j-2*bin1+2);
    diff=mean-down;
    histo_mc_up->SetBinContent(j-2*bin1+2,TMath::Max(mean+diff,0.000001));
  }
  histo_mc_up->SetName("histo_nonfiducial_mcUp");
  histo_mc_up->Write();

  /****************** LepEff & MomRes Down ***************************/

  TH1D* histo_LepEff_up;
  TH1D* histo_MomRes_up;
  TH1D* histo_LepEff_down=(TH1D*)temp->Clone("histo_LepEff_Down");

  TH1D* histo_MomRes_down=(TH1D*)temp->Clone("histo_MomRes_down");


  for(int k=1;k<=lep;k++){
    for(int i=bin1; i<=gen_bin+bin1-1; i++){
      histo_mean=(TH1D*) outFileLimits->Get(Form("histo_s%d",i));
      histo_LepEff_up=(TH1D*) outFileLimits->Get(Form("histo_s%d_LepEff%dUp",i,k));
      for(int j=2*bin1-1; j<=reco_bin+2*bin1-2; j++){
	mean=histo_mean->GetBinContent(j-2*bin1+2);
	up=histo_LepEff_up->GetBinContent(j-2*bin1+2);
	diff=up-mean;
	histo_LepEff_down->SetBinContent(j-bin1*2+2,TMath::Max(mean-diff,0.000001));
      }
      histo_LepEff_down->SetName(Form("histo_s%d_LepEff%dDown",i,k));
      histo_LepEff_down->Write();
      
    }
  }
  for(int k=1;k<=mom;k++){

    for(int i=bin1; i<=gen_bin+bin1-1; i++){
      histo_mean=(TH1D*) outFileLimits->Get(Form("histo_s%d",i));
      histo_MomRes_up=(TH1D*) outFileLimits->Get(Form("histo_s%d_MomRes%dUp",i,k));
      histo_MomRes_down=(TH1D*) outFileLimits->Get(Form("histo_s%d",i));
      histo_MomRes_down->Scale(2);
      histo_MomRes_down->Add(histo_MomRes_up,-1);
      histo_MomRes_down->SetName(Form("histo_s%d_MomRes%dDown",i,k));
      histo_MomRes_down->Write();
    }
  }

  /************ statistical uncertainty ****************/

  TH1D* stat_up=new TH1D("stat_up","stat_up",reco_bin,0,gen_bin);
  TH1D* stat_down=new TH1D("stat_down","stat_down",reco_bin,0,gen_bin);
  TH1D* stat_temp;

  for(int k=bin1; k<=gen_bin+bin1-1; k++){
    stat_temp=(TH1D*) outFileLimits->Get(Form("histo_s%d",k));
    for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
      for(int j=2*bin1-1; j<=reco_bin+2*bin1-2; j++){
	stat_up->SetBinContent(j-bin1*2+2,TMath::Max(stat_temp->GetBinContent(j-bin1*2+2),0.000001));
	stat_down->SetBinContent(j-bin1*2+2,TMath::Max(stat_temp->GetBinContent(j-bin1*2+2),0.000001));
      if (i==j){
	stat_up->SetBinContent(j-bin1*2+2,TMath::Max(stat_temp->GetBinContent(j-bin1*2+2)+TMath::Abs(stat_temp->GetBinError(j-bin1*2+2)),0.000001));
	stat_down->SetBinContent(j-bin1*2+2,TMath::Max(stat_temp->GetBinContent(j-bin1*2+2)-TMath::Abs(stat_temp->GetBinError(j-bin1*2+2)),0.000001));

      }
    }
    stat_down->SetName(Form("histo_s%d_stat_reco%d_gen%dDown",k,i,k));
    stat_up->SetName(Form("histo_s%d_stat_reco%d_gen%dUp",k,i,k));
    stat_down->Write();
    stat_up->Write();


  }
  }






  outFileLimits->Close();








  /******make data cards**********/
   char outputLimitsShape[200];                                            
   sprintf(outputLimitsShape,"test_200_%d%d%d.txt",me,DY,0);
  ofstream newcardShape;
  newcardShape.open(outputLimitsShape);
  newcardShape << Form("imax * number of channels\n");
  newcardShape << Form("jmax * number of background minus 1\n");
  newcardShape << Form("kmax * number of nuisance parameters\n");

  newcardShape << Form("shapes    *   *   %s  histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n",outputLimits);
  newcardShape << Form("shapes data_obs * %s  histo_Data \n",outputLimits);

  newcardShape << Form("Observation -1\n");

  newcardShape << Form("bin   ");
  for (int i=bin1; i<=gen_bin+bin1-1+bck; i++){
    newcardShape << Form("ch1  ");
  }
  newcardShape << Form("\n");
 
  newcardShape << Form("process  ");
  for (int i=bin1; i<=gen_bin+bin1-1; i++){
    newcardShape << Form("s%d  ",i);
  }
  newcardShape << Form("nonres resonant nonfiducial\n ");

  newcardShape << Form("process  ");
  for(int i=0;i>=-1*gen_bin+1;i--){
    newcardShape << Form("%d  ",i);
  }
  newcardShape << Form("1  2  3\n");
  

  newcardShape << Form("rate  ");
  for(int i=1;i<=gen_bin;i++){
    newcardShape << Form("-1   ");
  }
  newcardShape << Form("-1   -1   -1\n");

  newcardShape << Form("lumi    lnN     ");
  for (int i=1;i<=gen_bin;i++){
    newcardShape << Form("%6.3f  ",lumiE);
  }
  newcardShape << Form("-- %6.3f %6.3f \n",lumiE,lumiE);

  newcardShape << Form("QCD    shape     ");
  for(int i=1;i<=gen_bin;i++){
    newcardShape << Form("--  ");
  }
  newcardShape << Form("-- 1.0 --\n");
    

    newcardShape << Form("PDF    shape     ");
  for(int i=1;i<=gen_bin;i++){
    newcardShape << Form("--  ");
  }
  newcardShape << Form("-- 1.0  --\n");
  for(int j=1;j<=mom;j++){
    newcardShape << Form("MomRes%d    shape     ",j);
  for(int i=1;i<=gen_bin;i++){
    newcardShape << Form("1.0  ");
  }
  newcardShape << Form("-- 1.0  --\n");
  }
  for(int j=1;j<=lep;j++){
    newcardShape << Form("LepEff%d    shape     ",j);
  for(int i=1;i<=gen_bin;i++){
    newcardShape << Form("1.0  ");
  }
  newcardShape << Form("-- 1.0  --\n");
  }
  /* newcardShape << Form("mc    shape     ");
  for(int i=1;i<=gen_bin;i++){
    newcardShape << Form("1.0  ");
  }
  newcardShape << Form("-- -- 1.0 \n");*/
  
  newcardShape << Form ("nonres lnN  ");
  for(int i=1;i<=gen_bin;i++){
    newcardShape << Form("--  ");
  }
  newcardShape << Form("1.05  --  --\n");

  for(int i=2*bin1-1; i<=reco_bin+2*bin1-2; i++){
    for(int j=bin1; j<=gen_bin+bin1-1; j++){
      newcardShape << Form("stat_reco%d_gen%d     shape    ",i,j);

      for(int k=1;k<j-bin1+1;k++){
	newcardShape << Form("--  ");
      }
      newcardShape << Form("1.0  ");
      for(int k=j-bin1+2;k<=gen_bin+bck;k++){
	newcardShape << Form("--  ");
      }
      newcardShape << Form ("\n");

    }
  }





}
