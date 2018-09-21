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
#include "TRandom2.h"



void showYields(){

  TFile* fileIncl = TFile::Open("histowz_nice_4_2.root");
  //	TFile* fileIncl = TFile::Open("sample_23_histowz_nice_4_15.root");
    	TH1F*  Incl[8] ;    
    	for(int i=0; i<8; i++){
	  Incl[i]= (TH1F*)fileIncl->Get(Form("histo%d",i));
	}

	
	
		

   const int histBins = 8;
  double sumEventsType[5] = {0,0,0,0,0};
  double sumEventsTypeE[5] = {0,0,0,0,0};
TString processName[histBins] = {"..Data", ".Fakes", "Zgamma", "....WZ", "....ZZ", "...VVV","..EWWZ" ,".Higgs"};



sumEventsType[0] = 0.;sumEventsType[1] = 0.;sumEventsType[2] = 0.;sumEventsType[3] = 0.;sumEventsType[4] = 0.;
  sumEventsTypeE[0] = 0.;sumEventsTypeE[1] = 0.;sumEventsTypeE[2] = 0.;sumEventsTypeE[3] = 0.;sumEventsTypeE[4] = 0.;
  

  cout<<endl<<"-------   INCLUSIVE YIELDS    ---------"<<endl;
  printf("                  all      &           mmm      &           eem        &         mme    &             eee \\\ \n");
  printf("-----------------------------------------------------------------------------------------------------------\n");
for(int np=0; np<histBins; np++) {
    double sumEvents = 0.; double sumEventsE = 0.;
    for(int i=1; i<=4; i++) {sumEvents = sumEvents + Incl[np]->GetBinContent(i); sumEventsE = sumEventsE + Incl[np]->GetBinError(i)*Incl[np]->GetBinError(i);}
    if(np!=0){
    sumEventsType[0] = sumEventsType[0] + sumEvents;                          sumEventsTypeE[0] = sumEventsTypeE[0] + sumEventsE;
    sumEventsType[1] = sumEventsType[1] + Incl[np]->GetBinContent(1); sumEventsTypeE[1] = sumEventsTypeE[1] + Incl[np]->GetBinError(1) * Incl[np]->GetBinError(1);
    sumEventsType[2] = sumEventsType[2] + Incl[np]->GetBinContent(2); sumEventsTypeE[2] = sumEventsTypeE[2] + Incl[np]->GetBinError(2) * Incl[np]->GetBinError(2);
    sumEventsType[3] = sumEventsType[3] + Incl[np]->GetBinContent(3); sumEventsTypeE[3] = sumEventsTypeE[3] + Incl[np]->GetBinError(3) * Incl[np]->GetBinError(3);
    sumEventsType[4] = sumEventsType[4] + Incl[np]->GetBinContent(4); sumEventsTypeE[4] = sumEventsTypeE[4] + Incl[np]->GetBinError(4) * Incl[np]->GetBinError(4);
    }
    sumEventsE = sqrt(sumEventsE);
    printf("(%6s)& %7.2f $\\pm$ %5.2f & %7.2f $\\pm$ %5.2f & %7.2f $\\pm$ %5.2f & %7.2f $\\pm$ %5.2f & %7.2f $\\pm$ %5.2f \\\ \n",
    processName[np].Data(),sumEvents,sumEventsE,
    Incl[np]->GetBinContent(1),Incl[np]->GetBinError(1),Incl[np]->GetBinContent(2),Incl[np]->GetBinError(2),
    Incl[np]->GetBinContent(3),Incl[np]->GetBinError(3),Incl[np]->GetBinContent(4),Incl[np]->GetBinError(4));
    if(np==0)
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }
    printf("-----------------------------------------------------------------------------------------------------------\n");
    printf("(%6s)& %7.2f $\\pm$ %5.2f & %7.2f $\\pm$ %5.2f & %7.2f $\\pm$ %5.2f & %7.2f $\\pm$ %5.2f & %7.2f $\\pm$  %5.2f \\\ \n",
    "   all",sumEventsType[0],sqrt(sumEventsTypeE[0]),
    sumEventsType[1],sqrt(sumEventsTypeE[1]),sumEventsType[2],sqrt(sumEventsTypeE[2]),
    sumEventsType[3],sqrt(sumEventsTypeE[3]),sumEventsType[4],sqrt(sumEventsTypeE[4]));

	sumEventsType[0] = 0.;sumEventsType[1] = 0.;sumEventsType[2] = 0.;sumEventsType[3] = 0.;sumEventsType[4] = 0.;
  sumEventsTypeE[0] = 0.;sumEventsTypeE[1] = 0.;sumEventsTypeE[2] = 0.;sumEventsTypeE[3] = 0.;sumEventsTypeE[4] = 0.;
	
}
