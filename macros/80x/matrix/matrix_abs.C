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
void helper_function(int me=0,int DY=3){
  Double_t lumiE=1.025;
  int gen_bin=37;
  int reco_bin=74;
  int bck=3;
  TH1D* px;
  int lep=144;
  int mom=4;
  double min=0.0000001;
  const int nBinRecoPt     = 74; Float_t xbinsRecoPt[nBinRecoPt+1]         = { 0, 0.5,  1, 1.5,  2, 2.5, 3, 3.5,  4,4.5,  5,5.5,  6,6.5,  7,7.5,  8, 8.5,   9, 9.5,
									       10,10.5, 11,11.5, 12,12.5,13,13.5, 14, 15, 16, 17, 18, 19, 20, 21, 22,23.5,  25,26.5,
									       28,  30, 32,  35, 37,  40,43,  48, 52, 59, 65, 75, 86,100,120,140,160, 175, 190, 205,
									       220, 235,250, 275,300,325,350,375,400,425,450,475,500, 750,1000};
  /****** input data *******/
  TString path="/afs/cern.ch/user/c/christiw/public/panda/Unfolding/macros/inputs/";
  TFile* file_alt = new TFile(Form("data_test_%d%d.root",me,DY), "read");  if(!file_alt) {printf("File does not exist\n"); return;} //argument 0,0,1
  TFile* file = new TFile(Form("%shistoDY%dzllPtRecGen.root",path.Data(),DY), "read");  if(!file) {printf("File does not exist\n"); return;}
  TH1D* resonant=(TH1D*) file->Get(Form("histoPtRecVV_%d",me));
  TH1D* nonres=(TH1D*) file->Get(Form("histoPtRecEM_%d",me));
  TH2D* covariance_qcd = new TH2D(Form("covariance_qcd_%d",me), Form("covariance_qcd_%d",me), nBinRecoPt, xbinsRecoPt, nBinRecoPt, xbinsRecoPt);
  TH2D* covariance_pdf = new TH2D(Form("covariance_pdf_%d",me), Form("covariance_pdf_%d",me), nBinRecoPt, xbinsRecoPt, nBinRecoPt, xbinsRecoPt);
  TH2D* covariance_lumi = new TH2D(Form("covariance_lumi_%d",me), Form("covariance_lumi_%d",me), nBinRecoPt, xbinsRecoPt, nBinRecoPt, xbinsRecoPt);
  TH2D* covariance_stat[gen_bin];
  TH2D* covariance_stat_sum=new TH2D(Form("covariance_statsum_%d",me), Form("covariance_statsum_%d",me), nBinRecoPt, xbinsRecoPt, nBinRecoPt, xbinsRecoPt);
  TH2D* covariance_momres_sum=new TH2D(Form("covariance_momressum_%d",me), Form("covariance_momressum_%d",me), nBinRecoPt, xbinsRecoPt, nBinRecoPt, xbinsRecoPt);
  TH2D* covariance_lepeff_sum=new TH2D(Form("covariance_lepeffsum_%d",me), Form("covariance_lepeffsum_%d",me), nBinRecoPt, xbinsRecoPt, nBinRecoPt, xbinsRecoPt);
  

  TH2D* covariance_lepeff[lep];
  TH2D* covariance_momres[mom];
 
  TH2D* covariance_all=new TH2D(Form("covariance_all_%d",me), Form("covariance_all_%d",me), nBinRecoPt, xbinsRecoPt, nBinRecoPt, xbinsRecoPt);
  for(int i=1;i<=gen_bin;i++){
    covariance_stat[i-1] = new TH2D(Form("covariance_stat_%d%d",me,i), Form("covariance_stat_%d_%d",me,i), nBinRecoPt, xbinsRecoPt, nBinRecoPt, xbinsRecoPt);
  }
  for(int i=1;i<=lep;i++){
    covariance_lepeff[i-1] = new TH2D(Form("covariance_lepeff_%d%d",me,i), Form("covariance_lepeff_%d_%d",me,i), nBinRecoPt, xbinsRecoPt, nBinRecoPt, xbinsRecoPt);
  }
  for(int i=1;i<=mom;i++){
    covariance_momres[i-1] = new TH2D(Form("covariance_momres_%d%d",me,i), Form("covariance_momres_%d_%d",me,i), nBinRecoPt, xbinsRecoPt, nBinRecoPt, xbinsRecoPt);
  }

			 

  TH1D* histo_dy=(TH1D*) file->Get(Form("histoPtRecDY_%d",me));
  double sigma_qcd[reco_bin];
  double sigma_pdf[reco_bin];
  double sigma_lumi[reco_bin];
  double sigma_stat[reco_bin];
  double sigma_momres[reco_bin];
  double sigma_lepeff[reco_bin];
			 
  /******************* QCD & PDF ******************/
  TH1D* pdf=(TH1D*) file_alt->Get(Form("histo_resonant_PDFUp"));
  TH1D* qcd=(TH1D*) file_alt->Get(Form("histo_resonant_QCDUp"));

  for (int i=1;i<=reco_bin;i++){
    if(resonant->GetBinError(i)/resonant->GetBinContent(i)<0.2){
      sigma_pdf[i-1]=1.0*(pdf->GetBinContent(i)-resonant->GetBinContent(i));
      sigma_qcd[i-1]=1.0*(qcd->GetBinContent(i)-resonant->GetBinContent(i));
    }
    else{
      //sigma_pdf[i-1]=0.0;
      //sigma_qcd[i-1]=0.0;
    }
    
  }
  for (int i=1;i<=reco_bin;i++){
    for (int j=1;j<=reco_bin;j++){
      covariance_pdf->SetBinContent(i,j,sigma_pdf[i-1]*sigma_pdf[j-1]);
      covariance_qcd->SetBinContent(i,j,sigma_qcd[i-1]*sigma_qcd[j-1]);
    }
  }
  /*************** luminosity ********************/
  TH1D* temp;
  TH1D* normal;
  double tempbin;
  for (int i=1;i<=reco_bin;i++){
    for (int j=1;j<=reco_bin;j++){
      covariance_lumi->SetBinContent(i,j,0);
    }
  }
  for (int i=1;i<=gen_bin;i++){
    normal=(TH1D*) file_alt->Get(Form("histo_s%d",i));

      for(int k=1;k<=reco_bin;k++){
	if(normal->GetBinError(k)/normal->GetBinContent(k)<0.2){
	  sigma_lumi[k-1]=0.025*normal->GetBinContent(k);
	}
	else{
	  sigma_lumi[k-1]=0.025*normal->GetBinContent(k);
	}
	//if (sigma_lumi[k-1]<=min){sigma_lumi[k-1]=0.0;}

      }
      for(int i=1;i<=reco_bin;i++){
	for(int k=1;k<=reco_bin;k++){
	  tempbin=covariance_lumi->GetBinContent(k,i);
	  covariance_lumi->SetBinContent(k,i,tempbin+sigma_lumi[k-1]*sigma_lumi[i-1]);
	}
      }
  }
  
  normal=(TH1D*) file_alt->Get("histo_nonfiducial");
  for(int k=1;k<=reco_bin;k++){
    sigma_lumi[k-1]=0.025*normal->GetBinContent(k);
  }
  for(int i=1;i<=reco_bin;i++){
    for(int k=1;k<=reco_bin;k++){
      tempbin=covariance_lumi->GetBinContent(k,i);
      covariance_lumi->SetBinContent(k,i,tempbin+sigma_lumi[k-1]*sigma_lumi[i-1]);
    }
  }

  normal=(TH1D*) file_alt->Get("histo_resonant");
  for(int k=1;k<=reco_bin;k++){
    sigma_lumi[k-1]=0.025*normal->GetBinContent(k);
  }
  for(int i=1;i<=reco_bin;i++){
    for(int k=1;k<=reco_bin;k++){
      tempbin=covariance_lumi->GetBinContent(k,i);
      covariance_lumi->SetBinContent(k,i,tempbin+sigma_lumi[k-1]*sigma_lumi[i-1]);
    }
    }

  /************************* LepEff ******************/

  for(int k=1;k<=lep;k++){
    for (int i=1;i<=reco_bin;i++){
      for (int j=1;j<=reco_bin;j++){
	covariance_lepeff[k-1]->SetBinContent(i,j,0.0);
      }
    }
  }
  for (int i=1;i<=gen_bin;i++){
    normal=(TH1D*) file_alt->Get(Form("histo_s%d",i));

    for (int j=1;j<=lep;j++){
      temp=(TH1D*) file_alt->Get(Form("histo_s%d_LepEff%dUp",i,j));
      for(int k=1;k<=reco_bin;k++){
        sigma_lepeff[k-1]=1.0*temp->GetBinContent(k)-normal->GetBinContent(k);
        //sigma_lepeff[k-1]=1.0*sigma_lepeff[k-1]/normal->GetBinContent(k);
        //if(normal->GetBinError(k)/normal->GetBinContent(k)>0.2){sigma_lepeff[k-1]=0.0;}
      }
      for(int k=1;k<=reco_bin;k++){
        for(int q=1;q<=reco_bin;q++){
	  tempbin=covariance_lepeff[j-1]->GetBinContent(k,q);
	  covariance_lepeff[j-1]->SetBinContent(k,q,tempbin+sigma_lepeff[k-1]*sigma_lepeff[q-1]);
        }
      }
    }
  }
  normal=(TH1D*) file_alt->Get("histo_resonant");
  for (int i=1;i<=lep;i++){
    temp=(TH1D*) file_alt->Get(Form("histo_resonant_LepEff%dUp",i));
    for(int k=1;k<=reco_bin;k++){
      sigma_lepeff[k-1]=1.0*temp->GetBinContent(k)-normal->GetBinContent(k);
      //if(normal->GetBinError(k)/normal->GetBinContent(k)>0.2){sigma_lepeff[k-1]=0.0;}
    }
    for(int k=1;k<=reco_bin;k++){
      for(int q=1;q<=reco_bin;q++){
	tempbin=covariance_lepeff[i-1]->GetBinContent(k,q);
	covariance_lepeff[i-1]->SetBinContent(k,q,tempbin+sigma_lepeff[k-1]*sigma_lepeff[q-1]);
      }
    }
  }
  /************************* MomRes ******************/

  for(int k=1;k<=mom;k++){
    for (int i=1;i<=reco_bin;i++){
      for (int j=1;j<=reco_bin;j++){
	covariance_momres[k-1]->SetBinContent(i,j,0.0);
      }
    }
  }
  for (int i=1;i<=gen_bin;i++){
    normal=(TH1D*) file_alt->Get(Form("histo_s%d",i));

    for (int j=1;j<=mom;j++){
      temp=(TH1D*) file_alt->Get(Form("histo_s%d_MomRes%dUp",i,j));
      for(int k=1;k<=reco_bin;k++){
        sigma_momres[k-1]=1.0*temp->GetBinContent(k)-normal->GetBinContent(k);
        //if(normal->GetBinError(k)/normal->GetBinContent(k)>0.2){sigma_momres[k-1]=0.0;}
      }
      for(int k=1;k<=reco_bin;k++){
        for(int q=1;q<=reco_bin;q++){
	  tempbin=covariance_momres[j-1]->GetBinContent(k,q);
	  covariance_momres[j-1]->SetBinContent(k,q,tempbin+sigma_momres[k-1]*sigma_momres[q-1]);
        }
      }
    }
  }
  normal=(TH1D*) file_alt->Get("histo_resonant");
  for (int i=1;i<=mom;i++){
    temp=(TH1D*) file_alt->Get(Form("histo_resonant_MomRes%dUp",i));
    for(int k=1;k<=reco_bin;k++){
      sigma_momres[k-1]=1.0*temp->GetBinContent(k)-normal->GetBinContent(k);
      //if(normal->GetBinError(k)/normal->GetBinContent(k)>0.2){sigma_momres[k-1]=0.0;}
    }
    for(int k=1;k<=reco_bin;k++){
      for(int q=1;q<=reco_bin;q++){
	tempbin=covariance_momres[i-1]->GetBinContent(k,q);
	covariance_momres[i-1]->SetBinContent(k,q,tempbin+sigma_momres[k-1]*sigma_momres[q-1]);
      }
    }
  }
			  



  /******************* stat  ******************/
  for(int k=1;k<=gen_bin;k++){
  for (int i=1;i<=reco_bin;i++){
    for (int j=1;j<=reco_bin;j++){
      covariance_stat[k-1]->SetBinContent(i,j,0.0);
    }
  }}
  for (int i=1;i<=gen_bin;i++){
    normal=(TH1D*) file_alt->Get(Form("histo_s%d",i));

    for (int j=1;j<=reco_bin;j++){
      temp=(TH1D*) file_alt->Get(Form("histo_s%d_stat_reco%d_gen%dUp",i,j,i));
      for(int k=1;k<=reco_bin;k++){
	sigma_stat[k-1]=1.0*temp->GetBinContent(k)-normal->GetBinContent(k);
	//if(normal->GetBinError(k)/normal->GetBinContent(k)>0.2){sigma_stat[k-1]=0.0;}
      }
      for(int k=1;k<=reco_bin;k++){
	for(int q=1;q<=reco_bin;q++){
	tempbin=covariance_stat[i-1]->GetBinContent(k,q);
	covariance_stat[i-1]->SetBinContent(k,q,tempbin+sigma_stat[k-1]*sigma_stat[q-1]);
	}
      }
    }
  }

  /********** sum all stat ***********/
  for (int i=1;i<=reco_bin;i++){
    for (int j=1;j<=reco_bin;j++){
      covariance_stat_sum->SetBinContent(i,j,0.0);
    }}
  for (int i=1;i<=gen_bin;i++){
      covariance_stat_sum->Add(covariance_stat[i-1],1);

  }
  /********** sum all momres ***********/
  for (int i=1;i<=reco_bin;i++){
    for (int j=1;j<=reco_bin;j++){
      covariance_momres_sum->SetBinContent(i,j,0.0);
    }}
  for (int i=1;i<=mom;i++){
      covariance_momres_sum->Add(covariance_momres[i-1],1);

  }
  /********** sum all lepeff ***********/
  for (int i=1;i<=reco_bin;i++){
    for (int j=1;j<=reco_bin;j++){
      covariance_lepeff_sum->SetBinContent(i,j,0.0);
    }}
  for (int i=1;i<=lep;i++){
      covariance_lepeff_sum->Add(covariance_lepeff[i-1],1);

  }







  /********** ALL  ***********/
  for (int i=1;i<=reco_bin;i++){
    for (int j=1;j<=reco_bin;j++){
      covariance_all->SetBinContent(i,j,0.0);
    }
  }
  covariance_all->Add(covariance_stat_sum,1);
  covariance_all->Add(covariance_qcd,1);
  covariance_all->Add(covariance_pdf,1);
  covariance_all->Add(covariance_lumi,1);
  covariance_all->Add(covariance_lepeff_sum,1);
  covariance_all->Add(covariance_momres_sum,1);



  /*************new root file*********/


  char outputLimits[200];
  sprintf(outputLimits,"matrix_abs_%d%d.root",me,DY);
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  covariance_lumi->Write();  
  covariance_qcd->Write();
  covariance_pdf->Write();
  covariance_all->Write();
  covariance_stat_sum->Write();
  covariance_momres_sum->Write();		       
  covariance_lepeff_sum->Write();
  outFileLimits->Close();
  
}
void matrix_abs(){
  for(int i=0;i<=1;i++){
    for(int j=0;j<=3;j++){
      helper_function(i,j);
      printf("%d%d,",i,j);
    }
  }

}
