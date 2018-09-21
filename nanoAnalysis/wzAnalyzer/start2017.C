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




//#include "MitAnalysisRunII/macros/80x/factors.h"

//#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"


double superVar(double mjjIn, double jetDetaJJ, TLorentzVector Zp4, TLorentzVector Wp4, double Zep3l, bool isHiggsX=false ){

double bin = -9.;

 if(!isHiggsX){
   if((jetDetaJJ <= 2.5     || mjjIn<=500 || Zep3l > 2.5))      bin=0;

   if(jetDetaJJ <= 3.5 	&& jetDetaJJ > 2.5 && Zep3l < 2.5){
     if(mjjIn>500. 	&& mjjIn<=1000.)       	bin=1.;
     if(mjjIn>1000. 	&& mjjIn<=1350.) 	bin=2.;
     if(mjjIn>1350. 	&& mjjIn<=1750.) 	bin=3.;
     if(mjjIn>1750.)             		bin=4.;
   }

   if(jetDetaJJ <= 5. 	&& jetDetaJJ > 3.5 && Zep3l < 2.5){
     if(mjjIn>500. 	&& mjjIn<=1000.) 	bin=5.;
     if(mjjIn>1000. 	&& mjjIn<=1350.) 	bin=6.;
     if(mjjIn>1350. 	&& mjjIn<=1750.) 	bin=7.;
     if(mjjIn>1750.)             		bin=8.;
   }

   if( jetDetaJJ > 5. && Zep3l < 2.5){
     if(mjjIn>500. 	&& mjjIn<=1000.) 	bin=9.;
     if(mjjIn>1000. 	&& mjjIn<=1350.) 	bin=10.;
     if(mjjIn>1350. 	&& mjjIn<=1750.) 	bin=11.;
     if(mjjIn>1750.)             		bin=12.;
  }
 }

if(isHiggsX){
   if(jetDetaJJ <= 2.5     || mjjIn<=500 )      bin=-1;
   if(jetDetaJJ  > 2.5     && mjjIn >500 ){
    bin = TMath::Sqrt( TMath::Power(Zp4.Et() + Wp4.Et(),2) - TMath::Power(Zp4.Px() + Wp4.Px(),2) - TMath::Power(Zp4.Py() + Wp4.Py(),2));
   }
}

if( bin == -9 ) cout<<"superVar DANGER: "<<"mJJ: "<<mjjIn<<",deltaJJ: "<<jetDetaJJ<<",zep3l:  "<<Zep3l<<endl;

return bin;
}





void wzScattering2017(
 
 ){

  
    infilenamev.push_back(Form("%sdata_Run2016B.root",filesPathDA.Data())); infilecatv.push_back(0);
       TH1D* histo_Higgs_CMS_MVALepEffMperMassBoundingUp[11];
  TH1D* histo_Higgs_CMS_MVALepEffMperMassBoundingAvg[11];
  TH1D* histo_Higgs_CMS_MVALepEffEperMassBoundingDown[11];
  TH1D* histo_Higgs_CMS_MVALepEffEperMassBoundingUp[11];
  TH1D* histo_Higgs_CMS_MVALepEffEperMassBoundingAvg[11];

  for(int i=0; i<11; i++)
    { 
	histo_Higgs_CMS_MVALepEffMperMassBoundingDown[i]  	 = new TH1D( Form("histo_Higgs_"+SignalSuffix[i]+"_%sDown",effMName), Form("histo_Higgs_"+SignalSuffix[i]+"_%sDown",effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMperMassBoundingDown[i]->Sumw2();
	histo_Higgs_CMS_MVALepEffMperMassBoundingUp[i]  	 = new TH1D( Form("histo_Higgs_"+SignalSuffix[i]+"_%sUp",effMName), Form("histo_Higgs_"+SignalSuffix[i]+"_%sUp",effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMperMassBoundingUp[i]->Sumw2();
	histo_Higgs_CMS_MVALepEffMperMassBoundingAvg[i]  	 = new TH1D( Form("histo_Higgs_"+SignalSuffix[i]+"_%sAvg",effMName), Form("histo_Higgs_"+SignalSuffix[i]+"_%sAvg",effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMperMassBoundingAvg[i]->Sumw2();

	histo_Higgs_CMS_MVALepEffEperMassBoundingDown[i]  	 = new TH1D( Form("histo_Higgs_"+SignalSuffix[i]+"_%sDown",effEName), Form("histo_Higgs_"+SignalSuffix[i]+"_%sDown",effEName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEperMassBoundingDown[i]->Sumw2();
	histo_Higgs_CMS_MVALepEffEperMassBoundingUp[i]  	 = new TH1D( Form("histo_Higgs_"+SignalSuffix[i]+"_%sUp",effEName), Form("histo_Higgs_"+SignalSuffix[i]+"_%sUp",effEName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEperMassBoundingUp[i]->Sumw2();
	histo_Higgs_CMS_MVALepEffEperMassBoundingAvg[i]  	 = new TH1D( Form("histo_Higgs_"+SignalSuffix[i]+"_%sAvg",effEName), Form("histo_Higgs_"+SignalSuffix[i]+"_%sAvg",effEName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEperMassBoundingAvg[i]->Sumw2();
	}
  

  TH1D* histo_Zg_CMS_MVALepEffEBoundingUp        = new TH1D( Form("histo_Zg_%sUp",effEName)  , Form("histo_Zg_%sUp",effEName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVALepEffEBoundingDown      = new TH1D( Form("histo_Zg_%sDown",effEName), Form("histo_Zg_%sDown",effEName), nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effEName)  , Form("histo_VVV_%sUp",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effEName), Form("histo_VVV_%sDown",effEName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_WZ_%sUp",effEName)  , Form("histo_WZ_%sUp",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_WZ_%sDown",effEName), Form("histo_WZ_%sDown",effEName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_EWWZ_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_EWWZ_%sUp",effEName)  , Form("histo_EWWZ_%sUp",effEName)  , nBinMVA, xbins); histo_EWWZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_EWWZ_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_EWWZ_%sDown",effEName), Form("histo_EWWZ_%sDown",effEName), nBinMVA, xbins); histo_EWWZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_ZZ_%sUp",effEName)  , Form("histo_ZZ_%sUp",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_ZZ_%sDown",effEName), Form("histo_ZZ_%sDown",effEName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_Higgs_%sUp",effEName)  , Form("histo_Higgs_%sUp",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_Higgs_%sDown",effEName), Form("histo_Higgs_%sDown",effEName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_Zg_CMS_MVALepEffEBoundingAvg       = new TH1D( Form("histo_Zg_%sAvg",effEName)  , Form("histo_Zg_%sAvg",effEName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effEName)  , Form("histo_VVV_%sAvg",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_WZ_%sAvg",effEName)  , Form("histo_WZ_%sAvg",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();

 TH1D* histo_EWWZ_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_EWWZ_%sAvg",effEName)  , Form("histo_EWWZ_%sAvg",effEName)  , nBinMVA, xbins); histo_EWWZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();

  TH1D* histo_ZZ_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_ZZ_%sAvg",effEName)  , Form("histo_ZZ_%sAvg",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_Higgs_%sAvg",effEName)  , Form("histo_Higgs_%sAvg",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingAvg  ->Sumw2();

  TH1D* histo_Zg_CMS_MVAMETBoundingUp           = new TH1D( Form("histo_Zg_CMS_scale_metUp")  , Form("histo_Zg_CMS_scale_metUp")  	, nBinMVA, xbins); histo_Zg_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVAMETBoundingDown         = new TH1D( Form("histo_Zg_CMS_scale_metDown"), Form("histo_Zg_CMS_scale_metDown")	, nBinMVA, xbins); histo_Zg_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingUp   	= new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  	, nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingDown 	= new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown")	, nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_WZ_CMS_scale_metUp")  , Form("histo_WZ_CMS_scale_metUp")  	, nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_WZ_CMS_scale_metDown"), Form("histo_WZ_CMS_scale_metDown")	, nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingDown->Sumw2();

TH1D* histo_EWWZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_EWWZ_CMS_scale_metUp")  , Form("histo_EWWZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_EWWZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_EWWZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_EWWZ_CMS_scale_metDown"), Form("histo_EWWZ_CMS_scale_metDown"), nBinMVA, xbins); histo_EWWZ_CMS_MVAMETBoundingDown->Sumw2();

  TH1D* histo_ZZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_scale_metUp")  , Form("histo_ZZ_CMS_scale_metUp")  , 	nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_scale_metDown"), Form("histo_ZZ_CMS_scale_metDown"), 	nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_Higgs_CMS_scale_metUp")  , Form("histo_Higgs_CMS_scale_metUp"),	nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_Higgs_CMS_scale_metDown"), Form("histo_Higgs_CMS_scale_metDown"),nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_puUp               = new TH1D( Form("histo_Zg_CMS_puUp")  		, Form("histo_Zg_CMS_puUp")  , 	nBinMVA, xbins); histo_Zg_CMS_puUp  ->Sumw2();
  TH1D* histo_Zg_CMS_puDown             = new TH1D( Form("histo_Zg_CMS_puDown")		, Form("histo_Zg_CMS_puDown"), 	nBinMVA, xbins); histo_Zg_CMS_puDown->Sumw2();
  TH1D* histo_VVV_CMS_puUp           	= new TH1D( Form("histo_VVV_CMS_puUp")  	, Form("histo_VVV_CMS_puUp")  , nBinMVA, xbins); histo_VVV_CMS_puUp  ->Sumw2();
  TH1D* histo_VVV_CMS_puDown         	= new TH1D( Form("histo_VVV_CMS_puDown")	, Form("histo_VVV_CMS_puDown"), nBinMVA, xbins); histo_VVV_CMS_puDown->Sumw2();
  TH1D* histo_WZ_CMS_puUp            	= new TH1D( Form("histo_WZ_CMS_puUp")  		, Form("histo_WZ_CMS_puUp")  , 	nBinMVA, xbins); histo_WZ_CMS_puUp  ->Sumw2();
  TH1D* histo_WZ_CMS_puDown  	        = new TH1D( Form("histo_WZ_CMS_puDown")		, Form("histo_WZ_CMS_puDown"), 	nBinMVA, xbins); histo_WZ_CMS_puDown->Sumw2();

  TH1D* histo_EWWZ_CMS_puUp            	= new TH1D( Form("histo_EWWZ_CMS_puUp")  	, Form("histo_EWWZ_CMS_puUp")  , nBinMVA, xbins); histo_EWWZ_CMS_puUp  ->Sumw2();
  TH1D* histo_EWWZ_CMS_puDown  	        = new TH1D( Form("histo_EWWZ_CMS_puDown")	, Form("histo_EWWZ_CMS_puDown"), nBinMVA, xbins); histo_EWWZ_CMS_puDown->Sumw2();

  TH1D* histo_ZZ_CMS_puUp    	        = new TH1D( Form("histo_ZZ_CMS_puUp")  		, Form("histo_ZZ_CMS_puUp")  , 	nBinMVA, xbins); histo_ZZ_CMS_puUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_puDown  	        = new TH1D( Form("histo_ZZ_CMS_puDown")		, Form("histo_ZZ_CMS_puDown"), 	nBinMVA, xbins); histo_ZZ_CMS_puDown->Sumw2();
  TH1D* histo_Higgs_CMS_puUp[11];
  TH1D* histo_Higgs_CMS_puDown[11]; 
  
  for(unsigned int k=0; k<11; k++){
	histo_Higgs_CMS_puUp[k]    	        = new TH1D( "histo_Higgs_"+SignalSuffix[k]+"_CMS_puUp"  	, "histo_Higgs_"+SignalSuffix[k]+"_CMS_puUp"      	, nBinMVA, xbins); histo_Higgs_CMS_puUp[k]  ->Sumw2();
	histo_Higgs_CMS_puDown[k]    	        = new TH1D( "histo_Higgs_"+SignalSuffix[k]+"_CMS_puDown"  	, "histo_Higgs_"+SignalSuffix[k]+"_CMS_puDown"  	, nBinMVA, xbins); histo_Higgs_CMS_puDown[k]  ->Sumw2();

  }

  unsigned int numberOfLeptons = 3;

  double totalEventsProcess[50];
  std::vector<double> sumEventsProcess(infilenamev.size(), 0.0);
  

  TH1F* cutFlow = new TH1F("cutFlow"+TString(sampleID),"cutFlow"+TString(sampleID), 10, 0.,10.);
  TH1F* dRZ = new TH1F("dRZ"+TString(sampleID),"dRZ"+TString(sampleID), 100, 0.,3.14);
  cout<<"number of samples: "	<<infilenamev.size()<<endl;
//cout<<"sample ID: "		<<ifile<<"  " <<sampleID<<endl;
//  cout<<"sample name: "	<<infilenamev[ifile]<<endl;

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    cout<<"number of samples: "	<<infilenamev.size()<<endl;
    cout<<"sample ID: "		<<ifile<<"  " <<sampleID<<endl;
    cout<<"sample name: "	<<infilenamev[ifile]<<endl;
    if(ifile != sampleID && sampleID > -1) continue;
  
    
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());

    TFile *the_input_file=TFile::Open(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("Events");
    TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("Runs");
    //TTree *the_SelBit_tree= (TTree*)the_input_file->FindObjectAny("SelBit_tree");
    //TTree *the_PDF_tree   = (TTree*)the_input_file->FindObjectAny("pdfReweight");

    
    //cout<<"here"<<endl;


    TClonesArray Electron_pt;
    the_input_tree->SetBranchAddress("Electron_pt",Electron_pt);    

    
    
    
    
    double theMCPrescale = mcPrescale;
    int numberOfEvents = int(the_input_tree->GetEntries()/theMCPrescale);
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    
    for (int i=0; i < numberOfEvents; ++i) {
    
       


      the_input_tree->GetEntry(i);

      
      int initPDFTag = 0;
      if((*eventMonteCarlo.pdfRwgt).size() == 0) {initPDFTag = -1;}
      
      //else  {cout<<"GOOD PDF WEIGTHS !!"<<endl;}

      // WZ QCD initial
      if(infilecatv[ifile] == 3) {
	histo_WZ_CMS_QCDScaleInitial[0] = histo_WZ_CMS_QCDScaleInitial[0] + TMath::Abs((double)eventMonteCarlo.r1f2);
	histo_WZ_CMS_QCDScaleInitial[1] = histo_WZ_CMS_QCDScaleInitial[1] + TMath::Abs((double)eventMonteCarlo.r1f5);
	histo_WZ_CMS_QCDScaleInitial[2] = histo_WZ_CMS_QCDScaleInitial[2] + TMath::Abs((double)eventMonteCarlo.r2f1);
	histo_WZ_CMS_QCDScaleInitial[3] = histo_WZ_CMS_QCDScaleInitial[3] + TMath::Abs((double)eventMonteCarlo.r2f2);
	histo_WZ_CMS_QCDScaleInitial[4] = histo_WZ_CMS_QCDScaleInitial[4] + TMath::Abs((double)eventMonteCarlo.r5f1);
	histo_WZ_CMS_QCDScaleInitial[5] = histo_WZ_CMS_QCDScaleInitial[5] + TMath::Abs((double)eventMonteCarlo.r5f5);
        histo_WZ_CMS_QCDScaleInitial[6] = histo_WZ_CMS_QCDScaleInitial[6] + 1.0;
      }

      Bool_t passFilter[11] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
            for (int nt = 0; nt <(int)numtokens; nt++) {
	     if((*eventTrigger.triggerFired)[nt] == 0) continue;
             if(
	      (strcmp(tokens[nt],Form("HLT_Ele25_eta2p1_WPTight_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Ele27_eta2p1_WPLoose_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_IsoMu24_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_IsoTkMu24_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Ele27_WPTight_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Ele30_WPTight_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Ele35_WPLoose_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_IsoMu22_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_IsoTkMu22_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu45_eta2p1_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu50_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0)
		)passFilter[1]=kTRUE;
	      
	      
	   }

      

      	if(passFilter[0] == kFALSE) continue;
	cutFlow->Fill(2.);    
	
	if(passFilter[1] == kFALSE) continue;
	cutFlow->Fill(3.);    
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {

        if(selectIdIsoCut(typeLepSel.Data(),
		TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   	TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),
		(double)(*eventLeptons.iso)[nlep],
		(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))                                                                    
			{idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}

        else if( (((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) ){idTight.push_back(0); idLep.push_back(nlep); if(numberOfLeptons == 4) goodIsTight++;}
        else if( (((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP )){idSoft.push_back(nlep);}
      }
      // cout<<idLep.size()<<endl;
      if(idLep.size()!=idTight.size()) {assert(1); return;}
      if(idLep.size()==numberOfLeptons) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;
      cutFlow->Fill(4.);    
      

      bool passLeptonPtCut = false;
      //cout<<"pT 1: "<<((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt()<<endl;
      //cout<<"pT 2: "<<((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt()<<endl;
      //cout<<"pT 3: "<<((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt()<<endl;

      if(goodIsTight == idTight.size()) passFilter[3] = kTRUE;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 25 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() <= 15) continue;
      cutFlow->Fill(5.);    
      

      double dPhiLepMETMin = 999.;
      int signQ = 0;
      double systTotLep[2] = {1.0, 1.0}; // m/e
      for(unsigned nl=0; nl<idLep.size(); nl++){
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) systTotLep[0] = systTotLep[0] * 1.015;
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 11) systTotLep[1] = systTotLep[1] * 1.02;
        signQ = signQ + (int)(*eventLeptons.pdgId)[idLep[nl]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]);
        if(dPhiLepMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0]))))
           dPhiLepMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));      
      }
     
      double minMET  = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      double minPMET = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      if(dPhiLepMETMin < TMath::Pi()/2) minPMET = minPMET * sin(dPhiLepMETMin);
     
      passFilter[4] = TMath::Abs(signQ) == 1;
      if(passFilter[4] == kFALSE) continue;
      cutFlow->Fill(6.);    

      int nFakeCount = 0;
      double minMassll = 999.0;
      double minMassZ = 999.0;
      double mass3l = 0.0;
      double deltaRllMin = 999.0;
      int type3l = 0;
      int tagZ[3] = {-1,-1,-1};
      for(unsigned nl0=0; nl0<idLep.size()-1; nl0++){
        for(unsigned nl1=nl0+1; nl1<idLep.size(); nl1++){
          double deltaRllAux = ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl0]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl1]]));
          if(deltaRllAux < deltaRllMin) deltaRllMin = deltaRllAux;

	  if((int)(*eventLeptons.pdgId)[idLep[nl0]] * (int)(*eventLeptons.pdgId)[idLep[nl1]] > 0) continue;
          TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl1])) ) ));
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl1]]) &&
	     TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassZ-91.1876)) {
	     minMassZ = dilepAux.M();tagZ[0]=nl0;tagZ[1]=nl1;
	     if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]]) == 13) type3l = 0;
	     else                                                         type3l = 1;
	  }

	  if(minMassll > dilepAux.M()) minMassll = dilepAux.M();
        }
      }

      dRZ->Fill((double)((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->DeltaR( *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])));
 
      vector<int> idJet;vector<int> idJet50; vector<int> idJetJESUP; vector<int> idJetJESDOWN; vector<int> idJetJERUP; vector<int> idJetJERDOWN; double btagjet[2] = {0., 0.};
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      double theHT = 0;
      
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
                

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.4) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(dPhiJetMET   == -1) dPhiJetMET   = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 30 && (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() >= jetIDpTcut) {	
	  if     ((float)(*eventJets.bDiscr)[nj] > btagjet[0]) {btagjet[1] = btagjet[0]; btagjet[0] = (float)(*eventJets.bDiscr)[nj];}
	  else if((float)(*eventJets.bDiscr)[nj] > btagjet[1]) {btagjet[1] = (float)(*eventJets.bDiscr)[nj];}

	  theHT = theHT + ((TLorentzVector*)(*eventJets.p4)[nj])->Pt();
	  idJet.push_back(nj);
	  if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()>50.) idJet50.push_back(nj);
	}

	//JES and JER Uncertainties	
	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*(1+(float)(*eventJets.unc)[nj]) > jetIDpTcut) idJetJESUP.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*(1-(float)(*eventJets.unc)[nj]) > jetIDpTcut) idJetJESDOWN.push_back(nj);

	double sf_jer[2] = {1,1};
	if(infilecatv[ifile] != 0){
	  if((float)(*eventJets.ptResUncCentral)[nj] > 0 ){
	    sf_jer[0] = (float)(*eventJets.ptResUncUp)[nj]   / (float)(*eventJets.ptResUncCentral)[nj];
	    sf_jer[1] = (float)(*eventJets.ptResUncDown)[nj] / (float)(*eventJets.ptResUncCentral)[nj];
          }
	}
          if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*sf_jer[0] > jetIDpTcut) idJetJERUP.push_back(nj);
          if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*sf_jer[1] > jetIDpTcut) idJetJERDOWN.push_back(nj);

      }


      for(unsigned nl0=0; nl0<idLep.size(); nl0++){
        if((int)nl0==tagZ[0]||(int)nl0==tagZ[1]) continue;
        tagZ[2] = nl0;
        break;
      }
      
      if(tagZ[0] == -1 || tagZ[1] == -1 || tagZ[2] == -1) continue;
      
      bool tight3rdLepId = true;
      if(idTight[tagZ[2]] == 1){
        tight3rdLepId = selectIdIsoCut(type3rdLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt()),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Eta()),(double)(*eventLeptons.iso)[idLep[tagZ[2]]],(int)(*eventLeptons.selBits)[idLep[tagZ[2]]],(double)(*eventLeptons.mva)[idLep[tagZ[2]]]);
      }
      
      if(tight3rdLepId == false) continue;

      if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]) == 13) type3l += 0;
      else							       type3l += 2;
      TLorentzVector trilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + 
        		      ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) + 
        		      ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[2])) ) ));

      mass3l = trilep.M();
      double eta3l = trilep.Eta();
      
      
      passFilter[ 5] = minMassll > 4;
      passFilter[ 6] = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 30;
      passFilter[ 7] = minMassZ > minMass && minMassZ < maxMass;
      passFilter[ 8] = ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt() > 20;
      passFilter[ 9] = mass3l > 100;
      passFilter[10] = true;
      if(applyBtagging) passFilter[10] = bDiscrMax < 0.9535;

      if(passFilter[7]==kTRUE && (tagZ[0] == tagZ[1] || tagZ[0] == tagZ[2] || tagZ[1] == tagZ[2])) {printf("ZPROBLEM!\n");assert(0);return;}

      bool passNMinusOne[6] = {                 passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10],
                               passFilter[5] &&                  passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10],
                               passFilter[5] && passFilter[6] &&                  passFilter[8] && passFilter[9] && passFilter[10],
                               passFilter[5] && passFilter[6] && passFilter[7] &&                  passFilter[9] && passFilter[10],
                               passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8]                  && passFilter[10],
			       passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9]};

      bool passAllCuts[1] = {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};
     

      bool controlSel[4] = {passFilter[5] && passFilter[6] && !passFilter[7] && passFilter[9] && bDiscrMax > 0.605 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20,
                            passFilter[5] &&                  !passFilter[7] && passFilter[9] && bDiscrMax > 0.605 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20,
			    passFilter[5] && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() > 20,
			    passFilter[5] && !passFilter[6] && passFilter[7] && passFilter[9] &&                       ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20};

      bool passTTZSel[2] = {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && idJet.size() >= 4 && btagjet[0] > 0.560 && btagjet[1] > 0.560,
                            passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && idJet.size() == 3 && btagjet[0] > 0.800 && btagjet[1] > 0.560};


      double deltaPhiLeptonMet    = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtLN                 = TMath::Sqrt(2.0*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiLeptonMet)));
      double deltaPhiTriLeptonMet = TMath::Abs(trilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtEvent              = TMath::Sqrt(2.0*trilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiTriLeptonMet)));

     bool passSystCuts[nSystTypes] = {
          passFilter[5] && (double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt()   > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10],
	  passFilter[5] && (double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt()   > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]
     };


      // begin event weighting
      vector<bool> isGenDupl;
      int numberQuarks[2] = {0,0};
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 4 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[0]++;
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 5 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[1]++;
        isGenDupl.push_back(0);
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) isGenDupl[ngen0] = 1;
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) continue;
        for(int ngen1=ngen0+1; ngen1<eventMonteCarlo.p4->GetEntriesFast(); ngen1++) {
	  if((int)(*eventMonteCarlo.pdgId)[ngen0] != (int)(*eventMonteCarlo.pdgId)[ngen1]) continue;
          if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen1])) < 0.02) {
	    isGenDupl[ngen0] = 1;
	    break;
	  }
        }
      }
      vector<int> isGenLep; unsigned int goodIsGenLep = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        bool isGenLepton = false;
        for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	  if(isGenDupl[ngen] == 1) continue;
          if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.3) {
	    isGenLepton = true;
	    break;
	  }
	}
	if(isGenLepton == true) {isGenLep.push_back(1); goodIsGenLep++;}
	else                    {isGenLep.push_back(0);}
      }

      

      double pT_0 = ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt();
      double px_0 = ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Px();
      double py_0 = ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Px();

      double pT_1 = ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt();
      double px_1 = ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Px();
      double py_1 = ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Px();

      double pT_2 = ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt();
      double px_2 = ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Px();
      double py_2 = ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Px();

      double met_t = ((TLorentzVector*)(*eventMet.p4)[0])->Pt();
      double met_x = ((TLorentzVector*)(*eventMet.p4)[0])->Px();
      double met_y = ((TLorentzVector*)(*eventMet.p4)[0])->Py();
     
	
     TLorentzVector dijetjesUp,dijetjesDown;
      double deltaEtaJJjesUp = 0; double deltaEtaJJjesDown = 0;
      if(idJetJESUP.size() >= 2 ) {	
	dijetjesUp.SetPx(((TLorentzVector*)(*eventJets.p4)[idJetJESUP[0]])->Px()*(1+(float)(*eventJets.unc)[idJetJESUP[0]])+((TLorentzVector*)(*eventJets.p4)[idJetJESUP[1]])->Px()*(1+(float)(*eventJets.unc)[idJetJESUP[1]]));
	dijetjesUp.SetPy(((TLorentzVector*)(*eventJets.p4)[idJetJESUP[0]])->Py()*(1+(float)(*eventJets.unc)[idJetJESUP[0]])+((TLorentzVector*)(*eventJets.p4)[idJetJESUP[1]])->Py()*(1+(float)(*eventJets.unc)[idJetJESUP[1]]));
	dijetjesUp.SetPz(((TLorentzVector*)(*eventJets.p4)[idJetJESUP[0]])->Pz()*(1+(float)(*eventJets.unc)[idJetJESUP[0]])+((TLorentzVector*)(*eventJets.p4)[idJetJESUP[1]])->Pz()*(1+(float)(*eventJets.unc)[idJetJESUP[1]]));
	dijetjesUp.SetE (((TLorentzVector*)(*eventJets.p4)[idJetJESUP[0]])-> E()*(1+(float)(*eventJets.unc)[idJetJESUP[0]])+((TLorentzVector*)(*eventJets.p4)[idJetJESUP[1]])-> E()*(1+(float)(*eventJets.unc)[idJetJESUP[1]]));
	deltaEtaJJjesUp = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJetJESUP[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJetJESUP[1]])->Eta());
      }
      if(idJetJESDOWN.size() >= 2) {
	dijetjesDown.SetPx(((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[0]])->Px()*(1-(float)(*eventJets.unc)[idJetJESDOWN[0]])+((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[1]])->Px()*(1-(float)(*eventJets.unc)[idJetJESDOWN[1]]));
	dijetjesDown.SetPy(((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[0]])->Py()*(1-(float)(*eventJets.unc)[idJetJESDOWN[0]])+((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[1]])->Py()*(1-(float)(*eventJets.unc)[idJetJESDOWN[1]]));
	dijetjesDown.SetPz(((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[0]])->Pz()*(1-(float)(*eventJets.unc)[idJetJESDOWN[0]])+((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[1]])->Pz()*(1-(float)(*eventJets.unc)[idJetJESDOWN[1]]));
	dijetjesDown.SetE (((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[0]])-> E()*(1-(float)(*eventJets.unc)[idJetJESDOWN[0]])+((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[1]])-> E()*(1-(float)(*eventJets.unc)[idJetJESDOWN[1]]));
        deltaEtaJJjesDown = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[1]])->Eta());
      }

      TLorentzVector dijetjerUp,dijetjerDown;
      double deltaEtaJJjerUp = 0; double deltaEtaJJjerDown = 0;
      if(idJetJERUP.size() >= 2 && infilecatv[ifile]!=0) {
        double sfjets_jer[2] = {1,1};
        if((float)(*eventJets.ptResUncCentral)[idJetJERUP[0]] > 0) sfjets_jer[0] = (float)(*eventJets.ptResUncUp)[idJetJERUP[0]] / (float)(*eventJets.ptResUncCentral)[idJetJERUP[0]];
        if((float)(*eventJets.ptResUncCentral)[idJetJERUP[1]] > 0) sfjets_jer[1] = (float)(*eventJets.ptResUncUp)[idJetJERUP[1]] / (float)(*eventJets.ptResUncCentral)[idJetJERUP[1]];
	dijetjerUp.SetPx(((TLorentzVector*)(*eventJets.p4)[idJetJERUP[0]])->Px()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJetJERUP[1]])->Px()*sfjets_jer[1]);
	dijetjerUp.SetPy(((TLorentzVector*)(*eventJets.p4)[idJetJERUP[0]])->Py()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJetJERUP[1]])->Py()*sfjets_jer[1]);
	dijetjerUp.SetPz(((TLorentzVector*)(*eventJets.p4)[idJetJERUP[0]])->Pz()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJetJERUP[1]])->Pz()*sfjets_jer[1]);
	dijetjerUp.SetE (((TLorentzVector*)(*eventJets.p4)[idJetJERUP[0]])-> E()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJetJERUP[1]])-> E()*sfjets_jer[1]);
	deltaEtaJJjerUp = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJetJERUP[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJetJERUP[1]])->Eta());
      }
      if(idJetJERDOWN.size() >= 2 && infilecatv[ifile]!=0) {
        double sfjets_jer[2] = {1,1};
        if((float)(*eventJets.ptResUncCentral)[idJetJERDOWN[0]] > 0) sfjets_jer[0] = (float)(*eventJets.ptResUncDown)[idJetJERDOWN[0]] / (float)(*eventJets.ptResUncCentral)[idJetJERDOWN[0]];
        if((float)(*eventJets.ptResUncCentral)[idJetJERDOWN[1]] > 0) sfjets_jer[1] = (float)(*eventJets.ptResUncDown)[idJetJERDOWN[1]] / (float)(*eventJets.ptResUncCentral)[idJetJERDOWN[1]];
	dijetjerDown.SetPx(((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[0]])->Px()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[1]])->Px()*sfjets_jer[1]);
	dijetjerDown.SetPy(((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[0]])->Py()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[1]])->Py()*sfjets_jer[1]);
	dijetjerDown.SetPz(((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[0]])->Pz()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[1]])->Pz()*sfjets_jer[1]);
	dijetjerDown.SetE (((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[0]])-> E()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[1]])-> E()*sfjets_jer[1]);
	deltaEtaJJjerDown = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[1]])->Eta());
	}

    
     double transverseMass = TMath::Sqrt( TMath::Power(pT_0+pT_1+pT_2+met_t,2) - TMath::Power(px_0+px_1+px_2+met_x,2) - TMath::Power(py_0+py_1+py_2+met_y,2));

     TLorentzVector wp4 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + *((TLorentzVector*)(*eventMet.p4)[0]);
     TLorentzVector zp4 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]]) + *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]) ;
     double transverseMass2 = TMath::Sqrt( TMath::Power(zp4.Et() + wp4.Et(),2) - TMath::Power(zp4.Px() + wp4.Px(),2) - TMath::Power(zp4.Py() + wp4.Py(),2));

     //JES UNCERTAINTY
     TLorentzVector  wp4JESUP 	= *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + *((TLorentzVector*)(*eventMet.p4)[0])*((double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt()/met_t);
     TLorentzVector  wp4JESDOWN = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + *((TLorentzVector*)(*eventMet.p4)[0])*((double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])  ->Pt()/met_t);
     TLorentzVector  wp4JERUP 	= *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + *((TLorentzVector*)(*eventMet.p4)[0])*((double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JerUp])  ->Pt()/met_t);
     TLorentzVector  wp4JERDOWN = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + *((TLorentzVector*)(*eventMet.p4)[0])*((double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JerDown])  ->Pt()/met_t);
     TLorentzVector wp4LEPScaleUP = (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]))*1.01 + *((TLorentzVector*)(*eventMet.p4)[0]);
     TLorentzVector zp4LEPScaleUP = (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]]))*1.01 + (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]))*1.01 ;  
     TLorentzVector wp4LEPScaleDOWN = (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]))*0.99 + *((TLorentzVector*)(*eventMet.p4)[0]);
     TLorentzVector zp4LEPScaleDOWN = (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]]))*0.99 + (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]))*0.99 ;

     //JER UNCERTAINTY
     double metJESUP 		= ((double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt());//r.Gaus(met_t,0.1*met_t);
     double metJESDOWN 		= ((double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt());//r.Gaus(met_t,0.1*met_t);
     double metJERUP		= ((double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JerUp])  ->Pt());//r.Gaus(met_t,0.1*met_t);
     double metJERDOWN	 	= ((double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JerDown])->Pt());//r.Gaus(met_t,0.1*met_t);
     
     bool passAllCutsJESUP 	= {passFilter[5] &&  metJESUP   > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};
     bool passAllCutsJESDOWN    = {passFilter[5] &&  metJESDOWN > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};
     bool passAllCutsJERUP 	= {passFilter[5] &&  metJERUP   > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};
     bool passAllCutsJERDOWN    = {passFilter[5] &&  metJERDOWN > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};

     bool passAllCutsZEnriched  = {passFilter[5] &&     met_t < 30   && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};

     TLorentzVector dijet;
     double dijetDeltaEta	=0;
     bool passVBFLoose 		= false;
     bool passVBFLooseJESUP 	= false;
     bool passVBFLooseJESDOWN	= false;
     bool passVBFLooseJERUP 	= false;
     bool passVBFLooseJERDOWN 	= false;
     bool passVBFLooseLEPScaleUP	= false;
     bool passVBFLooseLEPScaleDOWN	= false;

     bool passVBF 		= false;
     bool passVBFzep 		= false;
     bool passVBFzep2Jet50 	= false;
     bool passVBFJESUP 		= false;
     bool passVBFJESDOWN	= false;
     bool passVBFJERUP 		= false;
     bool passVBFJERDOWN 	= false;	
     bool passVBFLEPScaleUP  	= false;
     bool passVBFLEPScaleDOWN 	= false;
     bool passVBFLooseMJJLess500 = false;
     bool passVBFLooseDEtaLess25 = false;
     bool pass2JetControl=false;
     double zep3l = 999;

     if(idJet.size()>=2){ 
	dijet 		= 	(( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) )   + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) ));
	dijetDeltaEta 	= 	TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()- ((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta());
	zep3l = TMath::Abs(eta3l - 0.5*(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta() + ((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()));

	passVBFLoose 		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>100.;     
      	passVBFLooseJESUP 	= passAllCutsJESUP 		&& idJetJESUP.size()   >= 2 	&& dijetjesUp.M()	>100.;
      	passVBFLooseJESDOWN 	= passAllCutsJESDOWN	 	&& idJetJESDOWN.size() >= 2 	&& dijetjesDown.M()	>100.;
      	passVBFLooseJERUP 	= passAllCutsJERUP 		&& idJetJERUP.size()   >= 2     && dijetjerUp.M()	>100.;
      	passVBFLooseJERDOWN 	= passAllCutsJERDOWN		&& idJetJERDOWN.size() >= 2 	&& dijetjerDown.M()	>100.;
      	passVBFLooseLEPScaleUP   = passVBFLoose;
      	passVBFLooseLEPScaleDOWN = passVBFLoose;

      	passVBF 		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>500.	&& dijetDeltaEta > 2.5;     
        passVBFzep 		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>500.	&& dijetDeltaEta > 2.5	&& zep3l<2.5;     
	passVBFzep2Jet50 	= passAllCuts[0] 		&& idJet50.size() >= 2 		&&  	(( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M()	>500.	       	&& TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()- ((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()) > 2.5				&& zep3l<2.5; 
    
      	passVBFJESUP 		= passAllCutsJESUP 		&& idJetJESUP.size()   >= 2 	&& dijetjesUp.M()	>500.	&& dijetDeltaEta > 2.5;
      	passVBFJESDOWN 		= passAllCutsJESDOWN	 	&& idJetJESDOWN.size() >= 2 	&& dijetjesDown.M()	>500.	&& dijetDeltaEta > 2.5;
      	passVBFJERUP 		= passAllCutsJERUP 		&& idJetJERUP.size()   >= 2     && dijetjerUp.M()	>500.	&& dijetDeltaEta > 2.5;
      	passVBFJERDOWN 		= passAllCutsJERDOWN		&& idJetJERDOWN.size() >= 2 	&& dijetjerDown.M()	>500.	&& dijetDeltaEta > 2.5;
      	passVBFLEPScaleUP  	= passVBF;
      	passVBFLEPScaleDOWN 	= passVBF;
	
	passVBFLooseMJJLess500  = passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		<500.;
	passVBFLooseDEtaLess25  = passAllCuts[0] 		&& idJet.size() >= 2 						&& dijetDeltaEta < 2.5;
	pass2JetControl		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>100.   &&(dijetDeltaEta < 2.5 || dijet.M()<500.) ;

	
	//cout<<zep3l<<endl;

      }
 
     
      bool pass1Jet = passAllCuts[0] && idJet.size() >= 1;
      
      // luminosity
      double theLumi  = 1.0; if(infilecatv[ifile] != 0) theLumi  = lumi;

      // pile-up
      double puWeight     = 1.0; if(infilecatv[ifile] != 0) puWeight     = nPUScaleFactor(fhDPU    , (double)eventMonteCarlo.puTrueInt);
      double puWeightUp   = 1.0; if(infilecatv[ifile] != 0) puWeightUp   = nPUScaleFactor(fhDPUUp  , (double)eventMonteCarlo.puTrueInt);
      double puWeightDown = 1.0; if(infilecatv[ifile] != 0) puWeightDown = nPUScaleFactor(fhDPUDown, (double)eventMonteCarlo.puTrueInt);

      // lepton efficiency
      double effSF = 1.0;
        if(infilecatv[ifile] != 0){
         for(unsigned int nl=0; nl<idLep.size(); nl++){
	   if(tagZ[2] != (int)nl)
	     effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
					     ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
					     typeLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,true,
					      fhDMuIsoSF,fhDVeryTightSF,true);

	   else
	     effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
					     ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
					     type3rdLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,
					      true,fhDMuIsoSF,fhDVeryTightSF,true);
        }
      }

      double trigEff = 1.0;


      // fake rate
      nFakeCount = 0;
      unsigned int typeFakeLepton[2] = {0,0};
     int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false && infilecatv[ifile]<7){
        if     ((infilecatv[ifile] == 0 || goodIsGenLep == isGenLep.size() || (goodIsGenLep < isGenLep.size() && infilecatv[ifile] == 2)) && goodIsTight != idTight.size() ){ // add Z+jets from data
	  for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    if(tagZ[0] == (int)nl) nFakeCount = nFakeCount + 1;
	    if(tagZ[1] == (int)nl) nFakeCount = nFakeCount + 2;
	    if(tagZ[2] == (int)nl) nFakeCount = nFakeCount + 4;

	    if(tagZ[2] != (int)nl) fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    else                   fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,type3rdLepSel.Data());

	    theCategory = 1;
	    if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) typeFakeLepton[0]++;
	    else                                                        typeFakeLepton[1]++;
	  }
          if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-3) fakeSF = -1.0 * fakeSF; // triple fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) fakeSF = +1.0 * fakeSF; // double fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) fakeSF = -1.0 * fakeSF; // single fake, MC
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-3) fakeSF = +1.0 * fakeSF; // triple fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) fakeSF = -1.0 * fakeSF; // double fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) fakeSF = +1.0 * fakeSF; // single fake, data	  
	}
        else if(infilecatv[ifile] != 0 && infilecatv[ifile] != 2 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;// make 0. 0 is correct, 1.0 is just for test
        }
        else if(infilecatv[ifile] == 2 && goodIsTight != idTight.size()){ // remove Z+gamma, fakeable objects
          fakeSF = 0.0;// make 0. 0 is correct, 1.0 is just for test
	 }
        else if(infilecatv[ifile] != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
	  fakeSF = 1.0;
        }
        else if(infilecatv[ifile] == 0 || infilecatv[ifile] == 2){ // data or Z+gamma with all good leptons
	  fakeSF = 1.0;
        }
	else {
	  printf("PROBLEM: %d %d %d %d %d\n",infilecatv[ifile],goodIsGenLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
	  assert(0);
	}
      }
 

      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;

      
      double totalWeight 			= trigEff*mcWeight*theLumi*effSF*puWeight*fakeSF*theMCPrescale;   
      if(infilecatv[ifile] == 7) totalWeight	= trigEff*mcWeight*theLumi*effSF*puWeight*fakeSF*theMCPrescale/1000.;
      if(totalWeight == 0) continue;
      if(theCategory == 0 ) totalWeight = 1.0;
      // end event weighting
      if((infilecatv[ifile] != 0 || theCategory == 0) && passAllCuts[0]) sumEventsProcess[ifile] += totalWeight;
      if(passAllCuts[0] && infilecatv[ifile] == 0) totalFakeDataCount[type3l][nFakeCount]++;

      
      double testVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[0])) )->Pt() + (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[1])) )->Pt() + (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[2])) )->Pt();
      double mvaLoose=-999.;
   if(theCategory == 0 && passAllCuts[0]){	  
	    eventNumbers[type3l].push_back(eventEvent.eventNum);
	    runNumbers[type3l].push_back(eventEvent.runNum);
	    lsNumbers[type3l].push_back(eventEvent.lumiNum);
	}
 

      for(int thePlot=0; thePlot<allPlots; thePlot++){
	//cout<<thePlot<<endl;
	double theVar = 0.0;
	bool makePlot = false;
 	if     (thePlot ==  0 && passAllCuts[0])                          		{makePlot = true;theVar = TMath::Min(transverseMass2,1199.999);}
	else if(thePlot ==  1 && passNMinusOne[0])                      		 {makePlot = true;theVar = TMath::Min(minMassll,199.999);}
	else if(thePlot ==  2 && pass2JetControl)                        		{makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot ==  3 && passNMinusOne[2])                       		{makePlot = true;theVar = TMath::Max(TMath::Min(minMassZ,149.999),50.001);}
	else if(thePlot ==  4 && passNMinusOne[3])                       		{makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt(),199.999);}
	else if(thePlot ==  5 && passNMinusOne[4])                       		{makePlot = true;theVar = TMath::Min(mass3l,399.999);}
	else if(thePlot ==  6 && passNMinusOne[5])                       		{makePlot = true;theVar = TMath::Min(bDiscrMax,0.999);}
	else if(thePlot ==  7 && passNMinusOne[0])                       		{makePlot = true;theVar = TMath::Min(deltaRllMin,3.999);}
	else if(thePlot ==  8 && passAllCuts[0])                         		{makePlot = true;theVar = TMath::Min(minPMET,199.999);}
	else if(thePlot ==  9 && passAllCuts[0])                         		{makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Pt(),199.999);}
	else if(thePlot == 10 && passAllCuts[0])                         		{makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Pt(),199.999);}
	else if(thePlot == 11 && passAllCuts[0])                         		{makePlot = true;theVar = TMath::Min(minMET,199.999);}
	else if(thePlot == 12 && controlSel[3] )                        		{makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 13 && passAllCuts[0])                         		{makePlot = true;theVar = TMath::Min((double)eventEvent.rho,39.999);}
	else if(thePlot == 14 && passAllCuts[0])                         		{makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	else if(thePlot == 15 && passAllCuts[0])                         		{makePlot = true;theVar = (double)type3l;}
	//else if(thePlot == 16 && passAllCuts[0])                         {makePlot = true;theVar = dPhiJetMET*180/TMath::Pi();}
	else if(thePlot == 16 && passVBFLoose)                           		{makePlot = true;theVar = testVar*180/TMath::Pi();}
	else if(thePlot == 17 && passAllCuts[0])                         		{makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	else if(thePlot == 18 && controlSel[0])                          		{makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 19 && controlSel[1])                          		{makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 20 && controlSel[2])                          		{makePlot = true;theVar = TMath::Max(TMath::Min(mass3l,249.999),50.001);}
	else if(thePlot == 22 && passAllCuts[0])                         		{makePlot = true;theVar = TMath::Min(mtEvent,999.999);}
	else if(thePlot == 21 && passVBFLoose)                           		{makePlot = true;theVar = TMath::Min(zep3l,7.);}
	else if(thePlot == 23 && passTTZSel[0])                          		{makePlot = true;theVar = TMath::Min(theHT,999.999);}
	else if(thePlot == 24 && passTTZSel[1])                          		{makePlot = true;theVar = TMath::Min(theHT,999.999);}
	else if(thePlot == 25 && passAllCuts[0])                         		{makePlot = true;theVar = TMath::Min((double)(numberQuarks[0]+10*numberQuarks[1]),49.499);}
	else if(thePlot == 26 && passNMinusOne[3])	                 		{makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt(),199.999);}
	else if(thePlot == 27 && passAllCuts[0])                         		{makePlot = true;theVar = TMath::Min(mtLN,199.999);}
	else if(thePlot == 28 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 11 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 29 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 13 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 30 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 11 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 31 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 13 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 32 && controlSel[2])                          {makePlot = true;theVar = TMath::Max(TMath::Min(mass3l,249.999),50.001);}
	else if(thePlot == 33 && passNMinusOne[5] && numberQuarks[1] == 0                  ) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 34 && passNMinusOne[5] && numberQuarks[1] == 0 && passFilter[10]) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}

	else if(thePlot == 35 && passAllCutsZEnriched)			 		{makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 36 && passNMinusOne[5] && numberQuarks[1]  > 0 && passFilter[10]) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 37 && pass2JetControl)                         		{makePlot = true;theVar = TMath::Min(transverseMass2,1199.999);}
	else if(thePlot == 38 && passVBFLoose)                            		{makePlot = true;theVar = TMath::Min((( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M(),1999.999);}
	else if(thePlot == 39 && passVBFLoose)                            		{makePlot = true;theVar = TMath::Min(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()),7.999);}
	//else if(thePlot == 40 && passVBFLoose)                            {makePlot = true;theVar = (double)type3l;}
	else if(thePlot == 40 && pass2JetControl)                            		{makePlot = true;theVar = (double)type3l;}
	else if(thePlot == 41 && pass2JetControl)                        		 {makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[0])) )->Pt();}
        else if(thePlot == 42 && pass2JetControl)                         		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[1])) )->Pt();}
        else if(thePlot == 43 && pass2JetControl)                        		 {makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[0])) )->Eta();}
        else if(thePlot == 44 && pass2JetControl)                         		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[1])) )->Eta();}
        else if(thePlot == 45 && pass2JetControl)                         		{makePlot = true;theVar = (double) TMath::Min((( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M(),1000.999);}
        else if(thePlot == 46 && pass2JetControl)                          		{makePlot = true;theVar = (double) TMath::Min(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()),5.49999);}
	else if(thePlot == 47 && pass2JetControl )                         		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[0])) )->Pt();}
	else if(thePlot == 48 && pass2JetControl )                         		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[1])) )->Pt();}
	else if(thePlot == 49 && pass2JetControl )                         		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[2])) )->Pt();}
        else if(thePlot == 50 && pass2JetControl )                        		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[0])) )->Eta();}
	else if(thePlot == 51 &&  pass2JetControl )                         		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[1])) )->Eta();}
	else if(thePlot == 52 &&  pass2JetControl)                         		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[2])) )->Eta();}
        else if(thePlot == 53 && passVBF && theCategory !=0 )               		{makePlot = true;theVar = (double)  type3l;}
	else if(thePlot == 54 && passVBF )                                  		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(tagZ[1])) )->Pt();}
	else if(thePlot == 55 && pass1Jet)                                 		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[0])) )->Pt();}
	else if(thePlot == 56 && passVBF )                          			{makePlot = true;theVar = (double)transverseMass2;}
	else if(thePlot == 57 && passVBFLooseMJJLess500 )                            	{makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta());}
	else if(thePlot == 58 && passVBFLooseDEtaLess25 )                            	{makePlot = true;theVar = (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M();}
        else if(thePlot == 59 && passVBFLooseMJJLess500 )                            	{makePlot = true;theVar = (double)TMath::Min(superVar(dijet.M(),	dijetDeltaEta,	zp4,	wp4, zep3l,	isHiggs),13.);}
	else if(thePlot == 60 && passVBFzep && theCategory != 0)                            	{makePlot = true;theVar = (double)type3l;}
	else if(thePlot == 61 && passVBFzep2Jet50 && theCategory != 0)                            	{makePlot = true;theVar = (double)type3l;}
	else if(thePlot == 62 && passVBFLooseDEtaLess25 )                            	{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[0])) )->Pt();}
	else if(thePlot == 63 && passVBFLooseDEtaLess25 )                            	{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[1])) )->Pt();}
	else if(thePlot == 64 && passVBFLooseDEtaLess25 )                            	{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[2])) )->Pt();}
	else if(thePlot == 65 && passVBFLooseDEtaLess25 )                            	{makePlot = true;theVar = (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M();}
	else if(thePlot == 66 && passVBFLooseMJJLess500 )                            	{makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta());}
	if( (makePlot && theCategory<7) || (makePlot && theCategory == 7 && infilenamev[ifile].Contains("M900")) ) {histo[     4][thePlot][theCategory]->Fill(theVar,totalWeight);}
      //if( (makePlot && theCategory<6) || (makePlot && theCategory == 6 && infilenamev[ifile].Contains("M900")) ) histo[     4][thePlot][theCategory]->Fill(theVar,totalWeight);
	if( (makePlot && theCategory<7) || (makePlot && theCategory == 7 && infilenamev[ifile].Contains("M900")) ) histo[type3l][thePlot][theCategory]->Fill(theVar,totalWeight);
      }

	// Avoid QCD scale and PDF weights that are anomalous high
      double maxQCDscale = (TMath::Abs((double)eventMonteCarlo.r1f2)+TMath::Abs((double)eventMonteCarlo.r1f5)+
	TMath::Abs((double)eventMonteCarlo.r2f1)+ TMath::Abs((double)eventMonteCarlo.r2f2)+TMath::Abs((double)eventMonteCarlo.r5f1)+
	TMath::Abs((double)eventMonteCarlo.r5f5))/6.0;

	double PDFAvg = 0.0;
      	if(infilecatv[ifile] != 0){
        	if(initPDFTag != -1)
        	for(int npdf=0; npdf<100; npdf++) PDFAvg = PDFAvg + TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]);
        	PDFAvg = PDFAvg/100.0;
      	}

      
      if(1 && (passVBFLoose || passVBFLooseJESDOWN || passVBFLooseJESUP || passVBFLooseJERUP || passVBFLooseJERDOWN || passVBFLooseLEPScaleUP || passVBFLooseLEPScaleDOWN) 
	 && ( /*theCategory!=0*/ 1 	/*|| (theCategory==0 	&& !passVBF)*/)) {
	//	if(1 && passVBFLoose && theCategory!=0) {


	double MVAVar 			= -999.;
	double MVAVarLEPScaleUP 	= -999.;
	double MVAVarLEPScaleDOWN 	= -999.;
	double MVAVarJESDOWN 		= -999.;
	double MVAVarJESUP   		= -999.;
	double MVAVarJERUP 		= -999.;
	double MVAVarJERDOWN 		= -999.;
	double bound      	        = 12.;
	if(isHiggs) bound 		= 1999.9;
      
	if(passVBFLoose){
	 
	  MVAVar 		= (double)TMath::Min(superVar(dijet.M(),	dijetDeltaEta,	zp4,	wp4 , zep3l,	isHiggs ),	bound);
	  //MVAVarLEPScaleUP 	= TMath::Min(MVAVar, bound);
	  //MVAVarLEPScaleDOWN 	= TMath::Min(MVAVar, bound);
	}          	
	if(passVBFLooseJESDOWN)	 	MVAVarJESDOWN 	= (double)TMath::Min(superVar(dijetjesDown.M()	,deltaEtaJJjesDown,	zp4,	wp4JESDOWN	,zep3l, isHiggs)	,bound);
	if(passVBFLooseJESUP)		MVAVarJESUP    	= (double)TMath::Min(superVar(dijetjesUp.M()	,deltaEtaJJjesUp,	zp4,	wp4JESUP	,zep3l,	isHiggs)	,bound);
	if(passVBFLooseJERDOWN) 	MVAVarJERDOWN 	= (double)TMath::Min(superVar(dijetjerDown.M()	,deltaEtaJJjerDown,	zp4,	wp4JERDOWN	,zep3l,	isHiggs)	,bound);
	if(passVBFLooseJERUP)		MVAVarJERUP    	= (double)TMath::Min(superVar(dijetjerUp.M()	,deltaEtaJJjerUp,	zp4,	wp4JESUP	,zep3l,	isHiggs)	,bound);
	if(passVBFLooseLEPScaleDOWN) 	MVAVarLEPScaleDOWN 	= (double)TMath::Min(superVar(dijet.M(),	dijetDeltaEta,	zp4LEPScaleDOWN,	wp4LEPScaleDOWN	, zep3l,	isHiggs)	,bound);
	if(passVBFLooseLEPScaleUP)	MVAVarLEPScaleUP    	= (double)TMath::Min(superVar(dijet.M(),	dijetDeltaEta,	zp4LEPScaleUP,		wp4LEPScaleUP	, zep3l,	isHiggs)	,bound);
	

	/////-------------------------        

        if     (theCategory == 0){
	  if(passVBFLoose) {histo_Data->Fill(MVAVar,totalWeight);}
	}

        else if(theCategory == 1){
	  if(passVBFLoose) {
		if(typeFakeLepton[0]+typeFakeLepton[1] != goodIsTight == idTight.size()) 
			{printf("PROBLEMFake %d %d %d %d\n",typeFakeLepton[0],typeFakeLepton[1],goodIsTight,(int)idTight.size()); return;}
	     	histo_Fake->Fill(MVAVar,totalWeight);
	  }
        }

        else if(theCategory == 2 ){
	  if(passVBFLoose) {
	   	  
	     histo_Zg->Fill(MVAVar,totalWeight);

	     histo_Zg_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_Zg_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_Zg_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_Zg_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_Zg_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_Zg_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<100; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_Zg_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Zg_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Zg_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Zg_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Zg_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_Zg_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_Zg_CMS_puUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_Zg_CMS_puDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
	  if(passVBFLooseJESUP) 	histo_Zg_CMS_scale_jUp		->Fill(MVAVarJESUP,totalWeight);
	  if(passVBFLooseJESDOWN)	histo_Zg_CMS_scale_jDown	->Fill(MVAVarJESDOWN,totalWeight);
	  if(passVBFLooseJERUP)		histo_Zg_CMS_res_jUp		->Fill(MVAVarJERUP,totalWeight);
	  if(passVBFLooseJERDOWN) 	histo_Zg_CMS_res_jDown		->Fill(MVAVarJERDOWN,totalWeight);
	  if(passVBFLooseLEPScaleUP) 	histo_Zg_CMS_scale_lUp		->Fill(MVAVarLEPScaleUP,totalWeight);
	  if(passVBFLooseLEPScaleDOWN) 	histo_Zg_CMS_scale_lDown	->Fill(MVAVarLEPScaleDOWN,totalWeight);
	}

        else if(theCategory == 3){
	  if(passVBFLoose) {

	    if(MVAVar == -999.) cout<<"THIS CANNOT BE"<<endl;
	     histo_WZ->Fill(MVAVar,totalWeight);
	     
	     histo_WZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<100; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
             else
	     for(int npdf=0; npdf<100; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_WZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_WZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_WZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_WZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_WZ_CMS_puUp  	->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_WZ_CMS_puDown	->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
	  if(passVBFLooseJESDOWN) 	histo_WZ_CMS_scale_jDown->Fill(MVAVarJESDOWN,totalWeight);
	  if(passVBFLooseJESUP) 	histo_WZ_CMS_scale_jUp->Fill(MVAVarJESUP,totalWeight);
	  if(passVBFLooseJERDOWN) 	histo_WZ_CMS_res_jDown->Fill(MVAVarJERDOWN,totalWeight);
	  if(passVBFLooseJERUP) 	histo_WZ_CMS_res_jUp->Fill(MVAVarJERUP,totalWeight);
	  if(passVBFLooseLEPScaleUP) 	histo_WZ_CMS_scale_lUp->Fill(MVAVarLEPScaleUP,totalWeight);
	  if(passVBFLooseLEPScaleDOWN) 	histo_WZ_CMS_scale_lDown->Fill(MVAVarLEPScaleDOWN,totalWeight);
          
	}
        else if(theCategory == 4){
	  if(passVBFLoose) {
	     histo_ZZ->Fill(MVAVar,totalWeight);
	     histo_ZZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_ZZ_CMS_puUp  	->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZZ_CMS_puDown	->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
          }
	  if(passVBFLooseJESDOWN) 	histo_ZZ_CMS_scale_jDown	->Fill(MVAVarJESDOWN,totalWeight);
          if(passVBFLooseJESUP) 	histo_ZZ_CMS_scale_jUp		->Fill(MVAVarJESUP,totalWeight);
	  if(passVBFLooseJERDOWN) 	histo_ZZ_CMS_res_jDown		->Fill(MVAVarJERDOWN,totalWeight);
          if(passVBFLooseJERUP) 	histo_ZZ_CMS_res_jUp		->Fill(MVAVarJERUP,totalWeight);
          if(passVBFLooseLEPScaleUP) 	histo_ZZ_CMS_scale_lUp		->Fill(MVAVarLEPScaleUP,totalWeight);
	  if(passVBFLooseLEPScaleDOWN) 	histo_ZZ_CMS_scale_lDown	->Fill(MVAVarLEPScaleDOWN,totalWeight);
          
        }
        else if(theCategory == 5){
	  if(passVBFLoose) {
	     histo_VVV->Fill(MVAVar,totalWeight);
	     histo_VVV_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	       for(int npdf=0; npdf<100; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<100; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_VVV_CMS_puUp  		->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_VVV_CMS_puDown		->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
          }
	  if(passVBFLooseJESDOWN) 	histo_VVV_CMS_scale_jDown	->Fill(MVAVarJESDOWN,totalWeight);
	  if(passVBFLooseJESUP) 	histo_VVV_CMS_scale_jUp		->Fill(MVAVarJESUP,totalWeight);
	  if(passVBFLooseJERDOWN) 	histo_VVV_CMS_res_jDown		->Fill(MVAVarJERDOWN,totalWeight);
	  if(passVBFLooseJERUP) 	histo_VVV_CMS_res_jUp		->Fill(MVAVarJERUP,totalWeight);
	  if(passVBFLooseLEPScaleUP) 	histo_VVV_CMS_scale_lUp		->Fill(MVAVarLEPScaleUP,totalWeight);
	  if(passVBFLooseLEPScaleDOWN) 	histo_VVV_CMS_scale_lDown	->Fill(MVAVarLEPScaleDOWN,totalWeight); 
	}

	else if(theCategory == 6){
	  if(passVBFLoose) {
	    if(MVAVar == -999.) cout<<"THIS CANNOT BE"<<endl;
	     histo_EWWZ->Fill(MVAVar,totalWeight);
	     histo_EWWZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_EWWZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_EWWZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_EWWZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_EWWZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_EWWZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	       for(int npdf=0; npdf<100; npdf++) histo_EWWZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
             else if(infilenamev[ifile].Contains("powheg") == true)
	       for(int npdf=0; npdf<100; npdf++) histo_EWWZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
             else
	       for(int npdf=0; npdf<100; npdf++) histo_EWWZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_EWWZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_EWWZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_EWWZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_EWWZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_EWWZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_EWWZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_EWWZ_CMS_puUp  	->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_EWWZ_CMS_puDown	->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  	}
	 if(passVBFLooseJESDOWN) 	histo_EWWZ_CMS_scale_jDown	->Fill(MVAVarJESDOWN,totalWeight);
       	 if(passVBFLooseJESUP) 		histo_EWWZ_CMS_scale_jUp	->Fill(MVAVarJESUP,totalWeight);
	 if(passVBFLooseJERDOWN) 	histo_EWWZ_CMS_res_jDown	->Fill(MVAVarJERDOWN,totalWeight);
	 if(passVBFLooseJERUP) 		histo_EWWZ_CMS_res_jUp		->Fill(MVAVarJERUP,totalWeight);
	 if(passVBFLooseLEPScaleUP) 	histo_EWWZ_CMS_scale_lUp	->Fill(MVAVarLEPScaleUP,totalWeight);
	 if(passVBFLooseLEPScaleDOWN) 	histo_EWWZ_CMS_scale_lDown	->Fill(MVAVarLEPScaleDOWN,totalWeight);		
       }

       else if( theCategory == 7){
	for(int i=0; i<11; i++) { 
	  if(infilenamev[ifile].Contains("_"+SignalSuffix[i]+"_")){	 
	 	 if(passVBFLoose) {	  
		   //cout<<initPDFTag<<endl;
		   //cout<<"totalWeight: "<<totalWeight<<" PDFweight*totalWeight: "<<totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[5+initPDFTag])<<endl; 
		   histo_Higgs_[i]->Fill(MVAVar,totalWeight);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
		   if(initPDFTag != -1)
		     for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBoundingPerMass[i][npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
		   else if(infilenamev[ifile].Contains("powheg") == true)
		     for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBoundingPerMass[i][npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
		   else
		     for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBoundingPerMass[i][npdf]->Fill(MVAVar,totalWeight);
		   histo_Higgs_CMS_MVALepEffMperMassBoundingAvg[i] ->Fill(MVAVar,totalWeight*1.00);
		   histo_Higgs_CMS_MVALepEffEperMassBoundingAvg[i] ->Fill(MVAVar,totalWeight*1.00);
		   histo_Higgs_CMS_MVALepEffMperMassBoundingUp[i]  ->Fill(MVAVar,totalWeight*systTotLep[0]);
		   histo_Higgs_CMS_MVALepEffEperMassBoundingUp[i]  ->Fill(MVAVar,totalWeight*systTotLep[1]);
		   histo_Higgs_CMS_MVALepEffMperMassBoundingDown[i]->Fill(MVAVar,totalWeight/systTotLep[0]);
		   histo_Higgs_CMS_MVALepEffEperMassBoundingDown[i]->Fill(MVAVar,totalWeight/systTotLep[1]);
		   histo_Higgs_CMS_puUp[i]  			->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
		   histo_Higgs_CMS_puDown[i]			->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	
	      	}
	       if(passVBFLooseJESDOWN) 		histo_Higgs_CMS_scale_jDown_[i]	->Fill(MVAVarJESDOWN,totalWeight);
	       if(passVBFLooseJESUP) 		histo_Higgs_CMS_scale_jUp_[i]	->Fill(MVAVarJESUP,totalWeight);
	       if(passVBFLooseJERDOWN) 		histo_Higgs_CMS_res_jDown_[i]	->Fill(MVAVarJERDOWN,totalWeight);
	       if(passVBFLooseJERUP) 		histo_Higgs_CMS_res_jUp_[i]	->Fill(MVAVarJERUP,totalWeight);
	       if(passVBFLooseLEPScaleUP) 	histo_Higgs_CMS_scale_lUp_[i]	->Fill(MVAVarLEPScaleUP,totalWeight);
	       if(passVBFLooseLEPScaleDOWN) 	histo_Higgs_CMS_scale_lDown_[i]	->Fill(MVAVarLEPScaleDOWN,totalWeight);
		   		
	      }

          }
	}
	else {
	  cout<<"the category: "<<theCategory<<endl;
	  
	  printf("CATEGORY PROBLEM!\n"); return;
	}
      } 
	// making data cards
    }
    printf("eff_cuts: %f\n",sumEventsProcess[ifile]);

  } // end of chain

  printf("----------------------totalFakeDataCount--------------------------------\n");
  printf("      TTT    FTT    TFT    FFT    TTF    FTF    TFF    FFF\n");
  for(int ni=0; ni<4; ni++) {
    printf("(%d): ",ni);
    for(int nj=0; nj<9; nj++) printf("%6.1f ",totalFakeDataCount[ni][nj]);
    printf("\n");
  }

  double sumEvents = 0;
  for(int np=1; np<histBins; np++) sumEvents += histo[4][0][np]->GetSumOfWeights();
  printf("yields: %f |",histo[4][0][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[4][0][np]->GetSumOfWeights());
  printf(" = %.3f\n",sumEvents);
  double sumEventsType[5] = {0,0,0,0,0};
  double sumEventsTypeE[5] = {0,0,0,0,0};

  cout<<endl<<"-------   TWO JET CONTROL YIELDS    ---------"<<endl;

  printf("                  all                 mmm                 eem                 mme                 eee\n");
  printf("-----------------------------------------------------------------------------------------------------------\n");
  for(int np=0; np<histBins; np++) {
    double sumEvents = 0.; double sumEventsE = 0.;
    for(int i=1; i<=4; i++) {sumEvents = sumEvents + histo[4][40][np]->GetBinContent(i); sumEventsE = sumEventsE + histo[4][40][np]->GetBinError(i)*histo[4][40][np]->GetBinError(i);}
    if(np!=0){
    sumEventsType[0] = sumEventsType[0] + sumEvents;                          sumEventsTypeE[0] = sumEventsTypeE[0] + sumEventsE;
    sumEventsType[1] = sumEventsType[1] + histo[4][40][np]->GetBinContent(1); sumEventsTypeE[1] = sumEventsTypeE[1] + histo[4][40][np]->GetBinError(1) * histo[4][40][np]->GetBinError(1);
    sumEventsType[2] = sumEventsType[2] + histo[4][40][np]->GetBinContent(2); sumEventsTypeE[2] = sumEventsTypeE[2] + histo[4][40][np]->GetBinError(2) * histo[4][40][np]->GetBinError(2);
    sumEventsType[3] = sumEventsType[3] + histo[4][40][np]->GetBinContent(3); sumEventsTypeE[3] = sumEventsTypeE[3] + histo[4][40][np]->GetBinError(3) * histo[4][40][np]->GetBinError(3);
    sumEventsType[4] = sumEventsType[4] + histo[4][40][np]->GetBinContent(4); sumEventsTypeE[4] = sumEventsTypeE[4] + histo[4][40][np]->GetBinError(4) * histo[4][40][np]->GetBinError(4);
    }
    sumEventsE = sqrt(sumEventsE);
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    processName[np].Data(),sumEvents,sumEventsE,
    histo[4][40][np]->GetBinContent(1),histo[4][40][np]->GetBinError(1),histo[4][40][np]->GetBinContent(2),histo[4][40][np]->GetBinError(2),
    histo[4][40][np]->GetBinContent(3),histo[4][40][np]->GetBinError(3),histo[4][40][np]->GetBinContent(4),histo[4][40][np]->GetBinError(4));
    if(np==0)
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }
    printf("-----------------------------------------------------------------------------------------------------------\n");
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    "   all",sumEventsType[0],sqrt(sumEventsTypeE[0]),
    sumEventsType[1],sqrt(sumEventsTypeE[1]),sumEventsType[2],sqrt(sumEventsTypeE[2]),
    sumEventsType[3],sqrt(sumEventsTypeE[3]),sumEventsType[4],sqrt(sumEventsTypeE[4]));

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    for(int j=4; j<allStates; j++){
      char output[200];
      sprintf(output,outputDirectory+"sample_%d_histowz_nice_%d_%d.root",sampleID,j,thePlot);	  
      TFile* outFilePlotsNote = new TFile(output,"recreate");
      outFilePlotsNote->cd();
      for(int np=0; np<histBins; np++) histo[j][thePlot][np]->Write();
      outFilePlotsNote->Close();
    }
  }


  sumEventsType[0] = 0.;sumEventsType[1] = 0.;sumEventsType[2] = 0.;sumEventsType[3] = 0.;sumEventsType[4] = 0.;
  sumEventsTypeE[0] = 0.;sumEventsTypeE[1] = 0.;sumEventsTypeE[2] = 0.;sumEventsTypeE[3] = 0.;sumEventsTypeE[4] = 0.;

cout<<endl<<"-------   TIGHT VBF YIELDS    ---------"<<endl;
printf("                  all                 mmm                 eem                 mme                 eee\n");
printf("-----------------------------------------------------------------------------------------------------------\n");
for(int np=0; np<histBins; np++) {
    double sumEvents = 0.; double sumEventsE = 0.;
    for(int i=1; i<=4; i++) {sumEvents = sumEvents + histo[4][53][np]->GetBinContent(i); sumEventsE = sumEventsE + histo[4][53][np]->GetBinError(i)*histo[4][53][np]->GetBinError(i);}
    if(np!=0){
    sumEventsType[0] = sumEventsType[0] + sumEvents;                          sumEventsTypeE[0] = sumEventsTypeE[0] + sumEventsE;
    sumEventsType[1] = sumEventsType[1] + histo[4][53][np]->GetBinContent(1); sumEventsTypeE[1] = sumEventsTypeE[1] + histo[4][53][np]->GetBinError(1) * histo[4][53][np]->GetBinError(1);
    sumEventsType[2] = sumEventsType[2] + histo[4][53][np]->GetBinContent(2); sumEventsTypeE[2] = sumEventsTypeE[2] + histo[4][53][np]->GetBinError(2) * histo[4][53][np]->GetBinError(2);
    sumEventsType[3] = sumEventsType[3] + histo[4][53][np]->GetBinContent(3); sumEventsTypeE[3] = sumEventsTypeE[3] + histo[4][53][np]->GetBinError(3) * histo[4][53][np]->GetBinError(3);
    sumEventsType[4] = sumEventsType[4] + histo[4][53][np]->GetBinContent(4); sumEventsTypeE[4] = sumEventsTypeE[4] + histo[4][53][np]->GetBinError(4) * histo[4][53][np]->GetBinError(4);
    }
    sumEventsE = sqrt(sumEventsE);
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    processName[np].Data(),sumEvents,sumEventsE,
    histo[4][53][np]->GetBinContent(1),histo[4][53][np]->GetBinError(1),histo[4][53][np]->GetBinContent(2),histo[4][53][np]->GetBinError(2),
    histo[4][53][np]->GetBinContent(3),histo[4][53][np]->GetBinError(3),histo[4][53][np]->GetBinContent(4),histo[4][53][np]->GetBinError(4));
    if(np==0)
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }
    printf("-----------------------------------------------------------------------------------------------------------\n");
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    "   all",sumEventsType[0],sqrt(sumEventsTypeE[0]),
    sumEventsType[1],sqrt(sumEventsTypeE[1]),sumEventsType[2],sqrt(sumEventsTypeE[2]),
    sumEventsType[3],sqrt(sumEventsTypeE[3]),sumEventsType[4],sqrt(sumEventsTypeE[4]));


  sumEventsType[0] = 0.;sumEventsType[1] = 0.;sumEventsType[2] = 0.;sumEventsType[3] = 0.;sumEventsType[4] = 0.;
  sumEventsTypeE[0] = 0.;sumEventsTypeE[1] = 0.;sumEventsTypeE[2] = 0.;sumEventsTypeE[3] = 0.;sumEventsTypeE[4] = 0.;


  cout<<endl<<"-------   INCLUSIVE YIELDS    ---------"<<endl;
  printf("                  all                 mmm                 eem                 mme                 eee\n");
  printf("-----------------------------------------------------------------------------------------------------------\n");
for(int np=0; np<histBins; np++) {
    double sumEvents = 0.; double sumEventsE = 0.;
    for(int i=1; i<=4; i++) {sumEvents = sumEvents + histo[4][15][np]->GetBinContent(i); sumEventsE = sumEventsE + histo[4][15][np]->GetBinError(i)*histo[4][15][np]->GetBinError(i);}
    if(np!=0){
    sumEventsType[0] = sumEventsType[0] + sumEvents;                          sumEventsTypeE[0] = sumEventsTypeE[0] + sumEventsE;
    sumEventsType[1] = sumEventsType[1] + histo[4][15][np]->GetBinContent(1); sumEventsTypeE[1] = sumEventsTypeE[1] + histo[4][15][np]->GetBinError(1) * histo[4][15][np]->GetBinError(1);
    sumEventsType[2] = sumEventsType[2] + histo[4][15][np]->GetBinContent(2); sumEventsTypeE[2] = sumEventsTypeE[2] + histo[4][15][np]->GetBinError(2) * histo[4][15][np]->GetBinError(2);
    sumEventsType[3] = sumEventsType[3] + histo[4][15][np]->GetBinContent(3); sumEventsTypeE[3] = sumEventsTypeE[3] + histo[4][15][np]->GetBinError(3) * histo[4][15][np]->GetBinError(3);
    sumEventsType[4] = sumEventsType[4] + histo[4][15][np]->GetBinContent(4); sumEventsTypeE[4] = sumEventsTypeE[4] + histo[4][15][np]->GetBinError(4) * histo[4][15][np]->GetBinError(4);
    }
    sumEventsE = sqrt(sumEventsE);
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    processName[np].Data(),sumEvents,sumEventsE,
    histo[4][15][np]->GetBinContent(1),histo[4][15][np]->GetBinError(1),histo[4][15][np]->GetBinContent(2),histo[4][15][np]->GetBinError(2),
    histo[4][15][np]->GetBinContent(3),histo[4][15][np]->GetBinError(3),histo[4][15][np]->GetBinContent(4),histo[4][15][np]->GetBinError(4));
    if(np==0)
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }
    printf("-----------------------------------------------------------------------------------------------------------\n");
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    "   all",sumEventsType[0],sqrt(sumEventsTypeE[0]),
    sumEventsType[1],sqrt(sumEventsTypeE[1]),sumEventsType[2],sqrt(sumEventsTypeE[2]),
    sumEventsType[3],sqrt(sumEventsTypeE[3]),sumEventsType[4],sqrt(sumEventsTypeE[4]));
  
    /*printf("QCD Init: WZ(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_WZ_CMS_QCDScaleInitial[6],histo_WZ_CMS_QCDScaleInitial[0],histo_WZ_CMS_QCDScaleInitial[1],histo_WZ_CMS_QCDScaleInitial[2],histo_WZ_CMS_QCDScaleInitial[3],histo_WZ_CMS_QCDScaleInitial[4],histo_WZ_CMS_QCDScaleInitial[5]);
  for(int nj=0; nj<6; nj++) histo_WZ_CMS_QCDScaleInitial[nj] = histo_WZ_CMS_QCDScaleInitial[nj] / histo_WZ_CMS_QCDScaleInitial[6];
  printf("QCD RateInit: WZ(%f/%f/%f/%f/%f/%f)\n",
  histo_WZ_CMS_QCDScaleInitial[0],histo_WZ_CMS_QCDScaleInitial[1],histo_WZ_CMS_QCDScaleInitial[2],histo_WZ_CMS_QCDScaleInitial[3],histo_WZ_CMS_QCDScaleInitial[4],histo_WZ_CMS_QCDScaleInitial[5]);*/

  // correcting by initial normalization
  /*for(int nj=0; nj<6; nj++) histo_WZ_CMS_QCDScaleBounding[nj]->Scale(1.0/histo_WZ_CMS_QCDScaleInitial[nj]);*/

  /*printf("QCD Corr: WZ(%f:%f/%f/%f/%f/%f/%f) ZZ(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) ZH(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_WZ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_Zg->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[5]->GetSumOfWeights());*/

  

  char outputLimits[200]; 
  sprintf(outputLimits,outputDirectory+"sample_%d_wz3l%2s_input_%4s.root",sampleID,finalStateName,ECMsb.Data());
 
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data   ->Write();
  /*for (int i=1;i<10;i++){
    if(histo_Zg->GetBinContent(i)<=0.){histo_Zg->SetBinContent(i,0.000001);histo_Zg->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_scale_jUp->GetBinContent(i)<=0.){histo_Zg_CMS_scale_jUp->SetBinContent(i,0.000001);histo_Zg_CMS_scale_jUp->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_scale_jDown->GetBinContent(i)<=0.){histo_Zg_CMS_scale_jDown->SetBinContent(i,0.000001);histo_Zg_CMS_scale_jDown->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_res_jUp->GetBinContent(i)<=0.){histo_Zg_CMS_res_jUp->SetBinContent(i,0.000001);histo_Zg_CMS_res_jUp->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_res_jDown->GetBinContent(i)<=0.){histo_Zg_CMS_res_jDown->SetBinContent(i,0.000001);histo_Zg_CMS_res_jDown->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_scale_lUp->GetBinContent(i)<=0.){histo_Zg_CMS_scale_lUp->SetBinContent(i,0.000001);histo_Zg_CMS_scale_lUp->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_scale_lDown->GetBinContent(i)<=0.){histo_Zg_CMS_scale_lDown->SetBinContent(i,0.000001); histo_Zg_CMS_scale_lDown->SetBinError(i,0.2*1.6);}
    if(histo_ZZ->GetBinContent(i)<=0.){histo_ZZ->SetBinContent(i,0.000001);histo_ZZ->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_scale_jUp->GetBinContent(i)<=0.){histo_ZZ_CMS_scale_jUp->SetBinContent(i,0.000001);histo_ZZ_CMS_scale_jUp->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_scale_jDown->GetBinContent(i)<=0.){histo_ZZ_CMS_scale_jDown->SetBinContent(i,0.000001);histo_ZZ_CMS_scale_jDown->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_res_jUp->GetBinContent(i)<=0.){histo_ZZ_CMS_res_jUp->SetBinContent(i,0.000001);histo_ZZ_CMS_res_jUp->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_res_jDown->GetBinContent(i)<=0.){histo_ZZ_CMS_res_jDown->SetBinContent(i,0.000001);histo_ZZ_CMS_res_jDown->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_scale_lUp->GetBinContent(i)<=0.){histo_ZZ_CMS_scale_lUp->SetBinContent(i,0.000001);histo_ZZ_CMS_scale_lUp->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_scale_lDown->GetBinContent(i)<=0.){histo_ZZ_CMS_scale_lDown->SetBinContent(i,0.000001);histo_ZZ_CMS_scale_lDown->SetBinError(i,0.2*1.6);}
    if(histo_WZ->GetBinContent(i)<=0.){histo_WZ->SetBinContent(i,0.000001);histo_WZ->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_scale_jUp->GetBinContent(i)<=0.){histo_WZ_CMS_scale_jUp->SetBinContent(i,0.000001);histo_WZ_CMS_scale_jUp->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_scale_jDown->GetBinContent(i)<=0.){histo_WZ_CMS_scale_jDown->SetBinContent(i,0.000001);histo_WZ_CMS_scale_jDown->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_res_jUp->GetBinContent(i)<=0.){histo_WZ_CMS_res_jUp->SetBinContent(i,0.000001);histo_WZ_CMS_res_jUp->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_res_jDown->GetBinContent(i)<=0.){histo_WZ_CMS_res_jDown->SetBinContent(i,0.000001);histo_WZ_CMS_res_jDown->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_scale_lUp->GetBinContent(i)<=0.){histo_WZ_CMS_scale_lUp->SetBinContent(i,0.000001);histo_WZ_CMS_scale_lUp->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_scale_lDown->GetBinContent(i)<=0.){histo_WZ_CMS_scale_lDown->SetBinContent(i,0.000001);histo_WZ_CMS_scale_lDown->SetBinError(i,0.2*1.6);}
    if(histo_EWWZ->GetBinContent(i)<=0.){histo_EWWZ->SetBinContent(i,0.000001);histo_EWWZ->SetBinError(i,0.2*1.6);}
    if(histo_EWWZ_CMS_scale_jUp->GetBinContent(i)<=0.){histo_EWWZ_CMS_scale_jUp->SetBinContent(i,0.000001);histo_EWWZ_CMS_scale_jUp->SetBinError(i,0.2*1.6);}
    if(histo_EWWZ_CMS_scale_jDown->GetBinContent(i)<=0.){histo_EWWZ_CMS_scale_jDown->SetBinContent(i,0.000001);histo_EWWZ_CMS_scale_jDown->SetBinError(i,0.2*1.6);}
    if(histo_EWWZ_CMS_res_jUp->GetBinContent(i)<=0.){histo_EWWZ_CMS_res_jUp->SetBinContent(i,0.000001);histo_EWWZ_CMS_res_jUp->SetBinError(i,0.2*1.6);}
    if(histo_EWWZ_CMS_res_jDown->GetBinContent(i)<=0.){histo_EWWZ_CMS_res_jDown->SetBinContent(i,0.000001);histo_EWWZ_CMS_res_jDown->SetBinError(i,0.2*1.6);}
    if(histo_EWWZ_CMS_scale_lUp->GetBinContent(i)<=0.){histo_EWWZ_CMS_scale_lUp->SetBinContent(i,0.000001);histo_EWWZ_CMS_scale_lUp->SetBinError(i,0.2*1.6);}
    if(histo_EWWZ_CMS_scale_lDown->GetBinContent(i)<=0.){histo_EWWZ_CMS_scale_lDown->SetBinContent(i,0.000001);histo_EWWZ_CMS_scale_lDown->SetBinError(i,0.2*1.6);}
    if(histo_VVV->GetBinContent(i)<=0.){histo_VVV->SetBinContent(i,0.000001);histo_VVV->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_scale_jUp->GetBinContent(i)<=0.){histo_VVV_CMS_scale_jUp->SetBinContent(i,0.000001);histo_VVV_CMS_scale_jUp->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_scale_jDown->GetBinContent(i)<=0.){histo_VVV_CMS_scale_jDown->SetBinContent(i,0.000001);histo_VVV_CMS_scale_jDown->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_res_jUp->GetBinContent(i)<=0.){histo_VVV_CMS_res_jUp->SetBinContent(i,0.000001);histo_VVV_CMS_res_jUp->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_res_jDown->GetBinContent(i)<=0.){histo_VVV_CMS_res_jDown->SetBinContent(i,0.000001);histo_VVV_CMS_res_jDown->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_scale_lUp->GetBinContent(i)<=0.){histo_VVV_CMS_scale_lUp->SetBinContent(i,0.000001);histo_VVV_CMS_scale_lUp->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_scale_lDown->GetBinContent(i)<=0.){histo_VVV_CMS_scale_lDown->SetBinContent(i,0.000001);histo_VVV_CMS_scale_lDown->SetBinError(i,0.2*1.6);}
    }*/



  for(int nb=1; nb<=nBinMVA; nb++){
    //QCD scale study
    	histo_Zg_CMS_QCDScaleBoundingUp		->SetBinContent	(nb, histo_Zg->GetBinContent(nb));
    	histo_Zg_CMS_QCDScaleBoundingDown	->SetBinContent	(nb, histo_Zg->GetBinContent(nb));
    	histo_VVV_CMS_QCDScaleBoundingUp	->SetBinContent	(nb, histo_VVV->GetBinContent(nb));
    	histo_VVV_CMS_QCDScaleBoundingDown	->SetBinContent	(nb, histo_VVV->GetBinContent(nb));
    	histo_WZ_CMS_QCDScaleBoundingUp		->SetBinContent	(nb, histo_WZ->GetBinContent(nb));
    	histo_WZ_CMS_QCDScaleBoundingDown	->SetBinContent	(nb, histo_WZ->GetBinContent(nb));
	histo_EWWZ_CMS_QCDScaleBoundingUp	->SetBinContent	(nb, histo_EWWZ->GetBinContent(nb));
    	histo_EWWZ_CMS_QCDScaleBoundingDown	->SetBinContent	(nb, histo_EWWZ->GetBinContent(nb));
    	histo_ZZ_CMS_QCDScaleBoundingUp		->SetBinContent	(nb, histo_ZZ->GetBinContent(nb));
    	histo_ZZ_CMS_QCDScaleBoundingDown	->SetBinContent	(nb, histo_ZZ->GetBinContent(nb));
    	histo_Higgs_CMS_QCDScaleBoundingUp	->SetBinContent	(nb, histo_Higgs->GetBinContent(nb));
    	histo_Higgs_CMS_QCDScaleBoundingDown	->SetBinContent	(nb, histo_Higgs->GetBinContent(nb));

    	for(int k=0;k<11; k++) {
      		histo_Higgs_CMS_QCDScaleBoundingUp_[k] 	-> SetBinContent(nb, histo_Higgs_[k]->GetBinContent(nb)); 
      		histo_Higgs_CMS_QCDScaleBoundingDown_[k] 	-> SetBinContent(nb, histo_Higgs_[k]->GetBinContent(nb)); 
      	}

	for(int nqcd=1; nqcd<6; nqcd++) {
	  	if(TMath::Abs(histo_Zg_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb))  > histo_Zg_CMS_QCDScaleBoundingUp 	->GetBinContent(nb))           	histo_Zg_CMS_QCDScaleBoundingUp->SetBinContent(nb, histo_Zg_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb));
	  	if(TMath::Abs(histo_Zg_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb))  < histo_Zg_CMS_QCDScaleBoundingDown 	->GetBinContent(nb))           	histo_Zg_CMS_QCDScaleBoundingDown->SetBinContent(nb, histo_Zg_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb));
	  	if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)) > histo_VVV_CMS_QCDScaleBoundingUp 	->GetBinContent(nb))           	histo_VVV_CMS_QCDScaleBoundingUp->SetBinContent(nb, histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb));
	  	if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)) < histo_VVV_CMS_QCDScaleBoundingDown 	->GetBinContent(nb))           	histo_VVV_CMS_QCDScaleBoundingDown->SetBinContent(nb, histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb));
	  	if(TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb))  > histo_WZ_CMS_QCDScaleBoundingUp 	->GetBinContent(nb))           	histo_WZ_CMS_QCDScaleBoundingUp->SetBinContent(nb, histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb));
	  	if(TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb))  < histo_WZ_CMS_QCDScaleBoundingDown 	->GetBinContent(nb))           	histo_WZ_CMS_QCDScaleBoundingDown->SetBinContent(nb, histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb));
		if(TMath::Abs(histo_EWWZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb))  > histo_EWWZ_CMS_QCDScaleBoundingUp 	->GetBinContent(nb))           	histo_EWWZ_CMS_QCDScaleBoundingUp->SetBinContent(nb, histo_EWWZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb));
	  	if(TMath::Abs(histo_EWWZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb))  < histo_EWWZ_CMS_QCDScaleBoundingDown 	->GetBinContent(nb))           	histo_EWWZ_CMS_QCDScaleBoundingDown->SetBinContent(nb, histo_EWWZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb));
	  	if(TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb))  > histo_ZZ_CMS_QCDScaleBoundingUp 	->GetBinContent(nb))           	histo_ZZ_CMS_QCDScaleBoundingUp->SetBinContent(nb, histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb));
	  	if(TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb))  < histo_ZZ_CMS_QCDScaleBoundingDown 	->GetBinContent(nb))           	histo_ZZ_CMS_QCDScaleBoundingDown->SetBinContent(nb, histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb));
         
	  for(int i=0;i<11; i++) {
	    	if(TMath::Abs(histo_Higgs_CMS_QCDScaleBoundingPerMass[i][nqcd]->GetBinContent(nb)) > histo_Higgs_CMS_QCDScaleBoundingUp_[i]   ->GetBinContent(nb))           	histo_Higgs_CMS_QCDScaleBoundingUp_[i]->SetBinContent(nb, histo_Higgs_CMS_QCDScaleBoundingPerMass[i][nqcd]->GetBinContent(nb));
	    	if(TMath::Abs(histo_Higgs_CMS_QCDScaleBoundingPerMass[i][nqcd]->GetBinContent(nb)) < histo_Higgs_CMS_QCDScaleBoundingDown_[i] ->GetBinContent(nb))             	histo_Higgs_CMS_QCDScaleBoundingDown_[i]->SetBinContent(nb, histo_Higgs_CMS_QCDScaleBoundingPerMass[i][nqcd]->GetBinContent(nb));
	
	  }

       } 

	// PDF study
    

    histo_Diff->Reset();
    if(histo_Zg ->GetBinContent(nb) > 0) {
      for(int npdf=1; npdf<100; npdf++) histo_Diff->Fill((histo_Zg_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_Zg ->GetBinContent(nb))/histo_Zg ->GetBinContent(nb));}
    histo_Zg_CMS_pdf_qqbarUp	->SetBinContent(nb, histo_Zg ->GetBinContent(nb) 	+ histo_Diff->GetRMS()*histo_Zg ->GetBinContent(nb));
    histo_Zg_CMS_pdf_qqbarDown	->SetBinContent(nb, histo_Zg ->GetBinContent(nb) 	- histo_Diff->GetRMS()*histo_Zg ->GetBinContent(nb));
    
    histo_Diff->Reset();
    if(histo_VVV->GetBinContent(nb) > 0) {
      for(int npdf=1; npdf<100; npdf++) histo_Diff->Fill((histo_VVV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))/histo_VVV->GetBinContent(nb));}
    histo_VVV_CMS_pdf_qqbarUp	->SetBinContent(nb, histo_VVV ->GetBinContent(nb) + histo_Diff->GetRMS()*histo_VVV ->GetBinContent(nb));
    histo_VVV_CMS_pdf_qqbarDown	->SetBinContent(nb, histo_VVV ->GetBinContent(nb) - histo_Diff->GetRMS()*histo_VVV ->GetBinContent(nb));
    
    histo_Diff->Reset();
    if(histo_WZ ->GetBinContent(nb) > 0) {
      for(int npdf=1; npdf<100; npdf++) histo_Diff->Fill((histo_WZ_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb))/histo_WZ ->GetBinContent(nb));}
    histo_WZ_CMS_pdf_qqbarUp	->SetBinContent(nb, histo_WZ ->GetBinContent(nb) + histo_Diff->GetRMS()*histo_WZ ->GetBinContent(nb));
    histo_WZ_CMS_pdf_qqbarDown	->SetBinContent(nb, histo_WZ ->GetBinContent(nb) - histo_Diff->GetRMS()*histo_WZ ->GetBinContent(nb));

     histo_Diff->Reset();
    if(histo_EWWZ ->GetBinContent(nb) > 0) {
      for(int npdf=1; npdf<100; npdf++) histo_Diff->Fill((histo_EWWZ_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_EWWZ ->GetBinContent(nb))/histo_EWWZ ->GetBinContent(nb));}
    histo_EWWZ_CMS_pdf_qqbarUp	->SetBinContent(nb, histo_EWWZ ->GetBinContent(nb) + histo_Diff->GetRMS()*histo_EWWZ ->GetBinContent(nb));
    histo_EWWZ_CMS_pdf_qqbarDown->SetBinContent(nb, histo_EWWZ ->GetBinContent(nb) - histo_Diff->GetRMS()*histo_EWWZ ->GetBinContent(nb));
    
    histo_Diff->Reset();
    if(histo_ZZ ->GetBinContent(nb) > 0) {
      for(int npdf=1; npdf<100; npdf++) histo_Diff->Fill((histo_ZZ_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb))/histo_ZZ ->GetBinContent(nb));}
    histo_ZZ_CMS_pdf_qqbarUp	->SetBinContent(nb, histo_ZZ ->GetBinContent(nb) + histo_Diff->GetRMS()*histo_ZZ ->GetBinContent(nb));
    histo_ZZ_CMS_pdf_qqbarDown	->SetBinContent(nb, histo_ZZ ->GetBinContent(nb) - histo_Diff->GetRMS()*histo_ZZ ->GetBinContent(nb));
    
    histo_Diff->Reset();		
	
    for (int i=0; i<11; i++){
      for(int npdf=1; npdf<100; npdf++) histo_Diff->Fill((histo_Higgs_CMS_PDFBoundingPerMass[i][npdf] ->GetBinContent(nb)-histo_Higgs_[i] ->GetBinContent(nb))/histo_Higgs_[i] ->GetBinContent(nb));
    histo_Higgs_CMS_pdf_qqbarUp_[i]		->SetBinContent(nb, histo_Higgs_[i] ->GetBinContent(nb) + histo_Diff->GetRMS()*histo_Higgs_[i] ->GetBinContent(nb));
    histo_Higgs_CMS_pdf_qqbarDown_[i]		->SetBinContent(nb, histo_Higgs_[i] ->GetBinContent(nb) - histo_Diff->GetRMS()*histo_Higgs_[i] ->GetBinContent(nb));
    histo_Diff->Reset();
    }

  		

  }

  histo_Zg     				->Write();
  histo_Zg_CMS_scale_jUp     		->Write();
  histo_Zg_CMS_scale_jDown     		->Write();
  histo_Zg_CMS_res_jUp     		->Write();
  histo_Zg_CMS_res_jDown     		->Write();
  histo_Zg_CMS_scale_lUp     		->Write();
  histo_Zg_CMS_scale_lDown     		->Write();
  histo_Zg_CMS_QCDScaleBoundingDown	->Write();
  histo_Zg_CMS_QCDScaleBoundingUp	->Write();
  histo_Zg_CMS_pdf_qqbarDown		->Write();
  histo_Zg_CMS_pdf_qqbarUp		->Write();
  histo_Zg_CMS_MVALepEffMBoundingDown	->Write();
  histo_Zg_CMS_MVALepEffMBoundingUp	->Write();
  histo_Zg_CMS_MVALepEffEBoundingDown	->Write();
  histo_Zg_CMS_MVALepEffEBoundingUp	->Write();
  histo_Zg_CMS_puUp     		->Write();
  histo_Zg_CMS_puDown     		->Write();

  histo_VVV    				->Write();
  histo_VVV_CMS_scale_jUp     		->Write();
  histo_VVV_CMS_scale_jDown     	->Write();
  histo_VVV_CMS_res_jUp     		->Write();
  histo_VVV_CMS_res_jDown     	 ->Write();
  histo_VVV_CMS_scale_lUp     		->Write();
  histo_VVV_CMS_scale_lDown     	->Write();
  histo_VVV_CMS_QCDScaleBoundingDown	->Write();
  histo_VVV_CMS_QCDScaleBoundingUp	->Write();
  histo_VVV_CMS_pdf_qqbarDown		->Write();
  histo_VVV_CMS_pdf_qqbarUp		->Write();
  histo_VVV_CMS_MVALepEffMBoundingDown	->Write();
  histo_VVV_CMS_MVALepEffMBoundingUp	->Write();
  histo_VVV_CMS_MVALepEffEBoundingDown	->Write();
  histo_VVV_CMS_MVALepEffEBoundingUp	->Write();
  histo_VVV_CMS_puUp     		->Write();
  histo_VVV_CMS_puDown     		->Write();


  histo_WZ     				->Write();
  histo_WZ_CMS_scale_jUp    		 ->Write();
  histo_WZ_CMS_scale_jDown     		->Write();
  histo_WZ_CMS_res_jUp     		->Write();
  histo_WZ_CMS_res_jDown     		->Write();
  histo_WZ_CMS_scale_lUp     		->Write();
  histo_WZ_CMS_scale_lDown     		->Write();
  histo_WZ_CMS_QCDunc_scaleDown		->Write();
  histo_WZ_CMS_QCDunc_scaleUp		->Write();
  histo_WZ_CMS_ewkscaleDown		->Write();
  histo_WZ_CMS_ewkscaleUp		->Write();
  histo_WZ_CMS_QCDScaleBoundingDown	->Write();
  histo_WZ_CMS_QCDScaleBoundingUp	->Write();
  histo_WZ_CMS_pdf_qqbarDown		->Write();
  histo_WZ_CMS_pdf_qqbarUp		->Write();
  histo_WZ_CMS_MVALepEffMBoundingDown	->Write();
  histo_WZ_CMS_MVALepEffMBoundingUp	->Write();
  histo_WZ_CMS_MVALepEffEBoundingDown	->Write();
  histo_WZ_CMS_MVALepEffEBoundingUp	->Write();
  histo_WZ_CMS_puUp     		->Write();
  histo_WZ_CMS_puDown     		->Write();

  histo_EWWZ    			->Write();
  histo_EWWZ_CMS_scale_jUp     		->Write();
  histo_EWWZ_CMS_scale_jDown     	->Write();
  histo_EWWZ_CMS_res_jUp     		->Write();
  histo_EWWZ_CMS_res_jDown     	        ->Write();
  histo_EWWZ_CMS_scale_lUp     		->Write();
  histo_EWWZ_CMS_scale_lDown     	->Write();
  histo_EWWZ_CMS_QCDScaleBoundingDown	->Write();
  histo_EWWZ_CMS_QCDScaleBoundingUp	->Write();
  histo_EWWZ_CMS_pdf_qqbarDown		->Write();
  histo_EWWZ_CMS_pdf_qqbarUp		->Write();
  histo_EWWZ_CMS_MVALepEffMBoundingDown	->Write();
  histo_EWWZ_CMS_MVALepEffMBoundingUp	->Write();
  histo_EWWZ_CMS_MVALepEffEBoundingDown	->Write();
  histo_EWWZ_CMS_MVALepEffEBoundingUp	->Write();
  histo_EWWZ_CMS_puUp     		->Write();
  histo_EWWZ_CMS_puDown     		->Write();


  histo_ZZ     				->Write();
  histo_ZZ_CMS_scale_jUp     		->Write();
  histo_ZZ_CMS_scale_jDown     		->Write();
  histo_ZZ_CMS_res_jUp     		->Write();
  histo_ZZ_CMS_res_jDown     		->Write();
  histo_ZZ_CMS_scale_lUp     		->Write();
  histo_ZZ_CMS_scale_lDown     		->Write();
  histo_ZZ_CMS_QCDScaleBoundingDown	->Write();
  histo_ZZ_CMS_QCDScaleBoundingUp	->Write();
  histo_ZZ_CMS_pdf_qqbarDown		->Write();
  histo_ZZ_CMS_pdf_qqbarUp		->Write();
  histo_ZZ_CMS_MVALepEffMBoundingDown	->Write();
  histo_ZZ_CMS_MVALepEffMBoundingUp	->Write();
  histo_ZZ_CMS_MVALepEffEBoundingDown	->Write();
  histo_ZZ_CMS_MVALepEffEBoundingUp	->Write();
  histo_ZZ_CMS_puUp     		->Write();
  histo_ZZ_CMS_puDown     		->Write();

  histo_Fake  				->Write();
  histo_FakeE  				->Write();

  for(int i=0; i<11; i++) {
	histo_Higgs_[i]  ->Write();
	histo_Higgs_CMS_scale_jUp_[i]  	  	->Write();
	histo_Higgs_CMS_scale_jDown_[i]  	->Write();
	histo_Higgs_CMS_res_jUp_[i]  		->Write();
	histo_Higgs_CMS_res_jDown_[i]  	->Write();
	histo_Higgs_CMS_scale_lUp_[i]  		->Write();
	histo_Higgs_CMS_scale_lDown_[i]  	->Write();
	histo_Higgs_CMS_QCDScaleBoundingUp_[i] 	->Write();
	histo_Higgs_CMS_QCDScaleBoundingDown_[i]->Write();
	histo_Higgs_CMS_pdf_qqbarDown_[i]	->Write();
  	histo_Higgs_CMS_pdf_qqbarUp_[i]	        ->Write();
	histo_Higgs_CMS_MVALepEffMperMassBoundingDown[i]->Write();
	histo_Higgs_CMS_MVALepEffMperMassBoundingUp[i]	->Write();
        histo_Higgs_CMS_MVALepEffEperMassBoundingDown[i]->Write();
	histo_Higgs_CMS_MVALepEffEperMassBoundingUp[i]	->Write();
	histo_Higgs_CMS_puDown[i]		->Write();
  	histo_Higgs_CMS_puUp[i]	        	->Write();
	}

  histo_Higgs_CMS_QCDScaleBoundingPerMass[1][0]->Write();
  histo_Higgs_CMS_QCDScaleBoundingPerMass[1][1]->Write();
  histo_Higgs_CMS_QCDScaleBoundingPerMass[1][2]->Write();
  histo_Higgs_CMS_QCDScaleBoundingPerMass[1][3]->Write();
  histo_Higgs_CMS_QCDScaleBoundingPerMass[1][4]->Write();
  histo_Higgs_CMS_QCDScaleBoundingPerMass[1][5]->Write();
  /*histo_EWWZ_CMS_QCDScaleBounding[0]->Write();
  histo_EWWZ_CMS_QCDScaleBounding[1]->Write();
  histo_EWWZ_CMS_QCDScaleBounding[2]->Write();
  histo_EWWZ_CMS_QCDScaleBounding[3]->Write();
  histo_EWWZ_CMS_QCDScaleBounding[4]->Write();
  histo_EWWZ_CMS_QCDScaleBounding[5]->Write();*/


  /*histo_Higgs_CMS_QCDScaleBoundingUp->Write();
    histo_Higgs_CMS_QCDScaleBoundingDown->Write(); */

  cout << histo_Data   		->GetSumOfWeights() << " ";
  cout << histo_Zg     		->GetSumOfWeights() << " ";
  cout << histo_VVV    		->GetSumOfWeights() << " ";
  cout << histo_WZ     		->GetSumOfWeights() << " ";
  cout << histo_EWWZ   		->GetSumOfWeights() << " ";
  cout << histo_ZZ     		->GetSumOfWeights() << " ";
  cout << histo_Fake   		->GetSumOfWeights() << " ";
  cout << histo_FakeE  		->GetSumOfWeights() << " ";
  cout << histo_Higgs_[1]  	->GetSumOfWeights() << " ";
  cout << endl;
  
  
  outFileLimits->Close();

  

  double lumiE = 1.025;
  double lumiEstab = 1.058;
  double systLepResE[4] = {1.01,1.01,1.01,1.01};
  double systLepResM[4] = {1.01,1.01,1.01,1.01};
for(int nb=1; nb<=2; nb++){
    if(nb<2){
      char outputLimitsShape[200];
      sprintf(outputLimitsShape,outputDirectory+"histo_limits_Higgswz3l%2s_%4s.txt",finalStateName,ECMsb.Data());
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("shapes * * wz3l3lHig_input_13TeV2016.root histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n");
      newcardShape << Form("shapes data_obs * wz3l3lHig_input_13TeV2016.root histo_Data\n");
      newcardShape << Form("shapes Higgs * wz3l3lHig_input_13TeV2016.root histo_Higgs_M$MASS histo_Higgs_M$MASS_$SYSTEMATIC\n");
      newcardShape << Form("Observation %d\n", -1);
      newcardShape << Form("bin wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
      newcardShape << Form("process Zg VVV WZ ZZ Fake EWWZ Higgs\n");
      newcardShape << Form("process 1 2 5 3 4 6 0 \n");
      newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f %8.5f  \n",-1.,-1.,-1.,-1.,-1.,-1.,-1.) ;
      newcardShape << Form("%s                               lnN  	%7.5f   %7.5f   %7.5f   %7.5f    -    %7.5f    %7.5f   \n","lumi_13TeV"         	,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE);
      newcardShape << Form("%s		                     lnN    	  -      -      %7.5f     -      -     -         -     \n","norm_WZ_13TeV2016"		,1.05);	
      newcardShape << Form("%s				     lnN  	%7.5f   %7.5f     -     %7.5f    -    %7.5f    %7.5f   \n","CMS_eff_b_mistag_13TeV2016" ,0.98,0.98,0.98,0.98,0.98);
      newcardShape << Form("%s		                     lnN   	  -       -       -       -    %7.5f    -        -     \n","CMS_wz3l_FakeSys" 		,1.30);	
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n",effMName			,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n",effEName			,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n","CMS_pu" 			,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n","CMS_scale_j"		,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n","CMS_res_j"		,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n","CMS_scale_l"		,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("%s	                             shape   	%7.5f     -       -      -      -       -        - 	\n","CMS_QCDScale_Zg"         ,1.);
      newcardShape << Form("%s	                             shape   	  -     %7.5f     -      -      -       -        - 	\n","CMS_QCDScale_VVV"        ,1.);
      newcardShape << Form("%s	                             shape        -       -      %7.5f   -      -       -        - 	\n","CMS_QCDScale_WZ"         ,1.);
      newcardShape << Form("%s	                             shape   	  -       -        -    %7.5f   -       -        - 	\n","CMS_QCDScale_ZZ"         ,1.);
      newcardShape << Form("%s	                             shape        -       -        -      -     -     %7.5f      - 	\n","CMS_QCDScale_EWWZ"       ,1.);
      newcardShape << Form("%s	                             shape        -       -        -      -     -       -      %7.5f 	\n","CMS_QCDScale_Higgs"      ,1.);
      //newcardShape << Form("%s                               lnN      	  -       -      -        -     -        -      %7.5f   \n","ewk_qqbar"			,1.20);
      newcardShape << Form("%s	                             shape   	%7.5f   %7.5f    %7.5f   %7.5f  -     %7.5f    %7.5f   \n","CMS_pdf_qqbar"                 ,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("* autoMCStats 0 \n");
      newcardShape.close();
    }
    else{
      char outputLimitsShape[200];
      sprintf(outputLimitsShape,outputDirectory+"histo_limits_EWWZwz3l%2s_%4s.txt",finalStateName,ECMsb.Data());
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("shapes * * wz3l3lEW_input_13TeV2016.root histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n");
      newcardShape << Form("shapes data_obs * wz3l3lEW_input_13TeV2016.root histo_Data\n");
      newcardShape << Form("Observation %d\n", -1);
      newcardShape << Form("bin wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
      newcardShape << Form("process Zg VVV WZ ZZ Fake EWWZ\n");
      newcardShape << Form("process 1 2 5 3 4 0 \n");
      newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",-1.,-1.,-1.,-1.,-1.,-1.) ;
      newcardShape << Form("%s                               lnN  	%7.5f   %7.5f   %7.5f   %7.5f    -    %7.5f       \n","lumi_13TeV"         	,lumiE,lumiE,lumiE,lumiE,lumiE);
      newcardShape << Form("%s		                   lnN    	  -      -      %7.5f     -      -     -          \n","norm_WZ_13TeV2016"		,1.05);	
      newcardShape << Form("%s				   lnN  	%7.5f   %7.5f     -     %7.5f    -    %7.5f       \n","CMS_eff_b_mistag_13TeV2016",0.98,0.98,0.98,0.98);
      newcardShape << Form("%s		                   lnN   	  -       -       -       -    %7.5f    -         \n","CMS_wz3l_FakeSys" 		,1.30);	
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n",effMName		        ,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n",effEName			,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n","CMS_pu" 			,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n","CMS_scale_j"	        ,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n","CMS_res_j"		,1.,1.,1.,1.,1.);
      //newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n","CMS_scale_l"		,1.,1.,1.,1.,1.);
      newcardShape << Form("%s	                             shape   	%7.5f     -       -      -      -       -         \n","CMS_QCDScale_Zg"         ,1.);
      newcardShape << Form("%s	                             shape   	  -     %7.5f     -      -      -       -         \n","CMS_QCDScale_VVV"        ,1.);
      newcardShape << Form("%s	                             shape        -       -      %7.5f   -      -       -         \n","CMS_QCDScale_WZ"         ,1.);
      newcardShape << Form("%s	                             shape   	  -       -        -    %7.5f   -       -         \n","CMS_QCDScale_ZZ"         ,1.);
      newcardShape << Form("%s	                             shape        -       -        -      -     -     %7.5f       \n","CMS_QCDScale_EWWZ"       ,1.);
      newcardShape << Form("%s	                             shape   	%7.5f   %7.5f    %7.5f  %7.5f   -     %7.5f       \n","CMS_pdf_qqbar"           ,1.,1.,1.,1.,1.);
      newcardShape << Form("* autoMCStats 0 \n");
      newcardShape.close();
    }
  


  }

  cout<<"histo_WZ "<<float(histo_WZ->Integral()/histo_WZ_CMS_puUp->Integral())<<endl;
  cout<<"histo_EWWZ "<<float(histo_EWWZ->Integral()/histo_EWWZ_CMS_puUp->Integral())<<endl;
  cout<<"histo_ZZ "<<float(histo_ZZ->Integral()/histo_ZZ_CMS_puUp->Integral())<<endl;
  cout<<"histo_VVV "<<float(histo_VVV->Integral()/histo_VVV_CMS_puUp->Integral())<<endl;
  cout<<"histo_Zg "<<float(histo_Zg->Integral()/histo_Zg_CMS_puUp->Integral())<<endl;
    
  /*for(int j=0; j<4;j++){
  char outputEventPerCat[200];
  sprintf(outputEventPerCat,outputDirectory+"outEvents_Type%d.txt",j);
  ofstream outEvents;
  outEvents.open(outputEventPerCat);
  for(unsigned int i=0; i<eventNumbers[j].size(); i++){
    outEvents<<Form("%ld %ld %ld\n",eventNumbers[j][i],lsNumbers[j][i],runNumbers[j][i]);
  }  
  outEvents.close();
  }*/
}


