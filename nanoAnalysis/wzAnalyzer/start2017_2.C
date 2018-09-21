#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include "TRandom2.h"
#include "TLorentzVector.h"
#include "wzAnalyzer/factors.h"



double superVar(double mjjIn, double jetDetaJJ, TLorentzVector Zp4, TLorentzVector Wp4, double Zep3l, bool isHiggsX=false ){

double bin = -9.;
 cout<<mjjIn<<" "<<jetDetaJJ<<" "<<Zep3l<<" "<<isHiggsX<<endl;
 if(!isHiggsX){
   if((jetDetaJJ <= 2.5 || mjjIn<=500 || Zep3l > 2.5))      bin=0;

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






void start2017_2(){


  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString>	infilenamev;
  vector<Int_t> 	infilecatv;
  vector<double> 	infilexsecv;
  float                 lumi=41290;//in pb

  double minMass =  76;
  double maxMass = 106;
  bool isHiggs=false;
  TString outputDirectory = "./";
  int sampleID=0;
  TString SignalSuffix[11] = {"M200","M300","M400","M500","M600","M700","M800","M900","M1000","M1500","M2000"};

  //xSecs in pb!
	
  //infilenamev.push_back("/eos/cms/store/user/jsalfeld/mergedDATAnanoALLSingleEleSingleMuFiltered.root"); infilecatv.push_back(0); infilexsecv.push_back(1);
  //infilenamev.push_back("/eos/cms/store/user/jsalfeld/DoubleMuon/doubleMuonbbbB.root"); infilecatv.push_back(0); infilexsecv.push_back(1);
  //infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/SingleMuon/SingleMuon70630e8a/merged.root"); infilecatv.push_back(0); infilexsecv.push_back(1);
  // infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/SingleElectron/SingleElectron70630e8a/merged.root"); infilecatv.push_back(0); infilexsecv.push_back(1);
 
  //infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/mergedSingleElectronSingleMuon2filtered.root"); infilecatv.push_back(0); infilexsecv.push_back(1);
infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8/NanoTestPost6/180429_134358/0000/tree_1.root"); infilecatv.push_back(6); infilexsecv.push_back(0.0175);

infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/mergedSingleElectronSingleMuon2filteredAndDoubleEleMuonfiltered.root"); infilecatv.push_back(0); infilexsecv.push_back(1);

  infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/merged.root"); infilecatv.push_back(3); infilexsecv.push_back(4.715*1.109);


  infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8/merged.root"); infilecatv.push_back(4); infilexsecv.push_back(0.001720*2.3);
  infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8/NanoTestPost6/180327_094001/0000/tree_2.root"); infilecatv.push_back(4); infilexsecv.push_back(0.001720*2.3);
infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8/merged.root"); infilecatv.push_back(4); infilexsecv.push_back(0.001720*2.3);
  infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8/merged.root"); infilecatv.push_back(4); infilexsecv.push_back(0.001720*2.3);
infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/ZZTo4L_13TeV_powheg_pythia8/NanoTestPost6/180327_094722/0000/tree_21.root"); infilecatv.push_back(4); infilexsecv.push_back(1.256);

infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/merged.root"); infilecatv.push_back(5); infilexsecv.push_back(0.2086);
infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/merged.root"); infilecatv.push_back(5); infilexsecv.push_back(0.05565);
infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/NanoTestPost6/180327_093750/0000/tree_1.root"); infilecatv.push_back(5); infilexsecv.push_back(0.01398);
infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/NanoTestPost6/180327_093833/0000/tree_2.root"); infilecatv.push_back(5); infilexsecv.push_back(0.0758);
infilenamev.push_back("/eos/cms/store/user/jsalfeld/DATAnano2/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/merged.root"); infilecatv.push_back(5); infilexsecv.push_back(0.2043);


 int dataevents=0;

//Lepton Scale Factors
 TFile* eleRECOSFfile = TFile::Open("egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root");
 TFile* eleIDSFfile = TFile::Open("egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root");
 TH2F* histEleRECOSF=(TH2F*) eleRECOSFfile->Get("EGamma_SF2D");
 TH2F* histEleIDSF=(TH2F*) eleIDSFfile->Get("EGamma_SF2D");

 TFile* muIDSFfile = TFile::Open("histEffMuID.root");
 TH2F* histMuIDSF=(TH2F*) muIDSFfile->Get("histEffMuID");

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 12;
  const int histBins = 8;
  const int allStates = 5;

  TH1D* checkTheID = new TH1D("checkTheID","checkTheID", 10, 0, 10);

  TH1D* histo[allStates][allPlots][histBins];
  TString processName[histBins] = {"..Data", ".Fakes", "Zgamma", "....WZ", "....ZZ", "...VVV","..EWWZ" ,".Higgs"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 50; 	xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >=  1 && thePlot <=  1) {nBinPlot = 40; 	xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  2 && thePlot <=  2) {nBinPlot = 4; 	xminPlot =-0.5; xmaxPlot = 3.5;}
    else if(thePlot >=  3 && thePlot <=  3) {nBinPlot = 7; 	xminPlot =-0.5; xmaxPlot = 6.5;}
    else if(thePlot >=  4 && thePlot <=  4) {nBinPlot = 50; 	xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  5 && thePlot <=  5) {nBinPlot = 50; 	xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  6 && thePlot <=  6) {nBinPlot = 7; 	xminPlot = -0.5; xmaxPlot = 6.5;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot = 50; 	xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >=  8 && thePlot <=  8) {nBinPlot = 20; 	xminPlot = -5.0; xmaxPlot = 5.0;}
    else if(thePlot >=  9 && thePlot <=  9) {nBinPlot = 20; 	xminPlot = -3.0; xmaxPlot = 3.0;}
    else if(thePlot >=  10 && thePlot <=  10) {nBinPlot = 20; 	xminPlot = -3.0; xmaxPlot = 3.0;}
    else if(thePlot >=  11 && thePlot <=  11) {nBinPlot = 60; 	xminPlot = 0.0; xmaxPlot = 61.0;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) {
      for(int j=0; j<allStates; j++) histo[j][thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    }
    histos->Reset();histos->Clear();
  }


  TString ECMsb  = "13TeV2016";
  const int nBinMVA = 13; Float_t xbins[nBinMVA+1] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13. }; 
  Float_t xbinsX[nBinMVA+1] = {-1.,0.,50.,100.,150.,200.,250.,300.,400.,500.,700.,1000.,1500.,2000.};
  if( isHiggs)  {for(int k=0; k<nBinMVA+1; k++) xbins[k] = xbinsX[k];}
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();

  TH1D *histo_Data   = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_Zg     = (TH1D*) histoMVA->Clone("histo_Zg"); 
  TH1D *histo_VVV    = (TH1D*) histoMVA->Clone("histo_VVV");	 
  TH1D *histo_WZ     = (TH1D*) histoMVA->Clone("histo_WZ");	 
  TH1D *histo_ZZ     = (TH1D*) histoMVA->Clone("histo_ZZ");	 
  TH1D *histo_Fake   = (TH1D*) histoMVA->Clone("histo_Fake");	 
  TH1D *histo_FakeE  = (TH1D*) histoMVA->Clone("histo_FakeE");	 
  TH1D *histo_Higgs  = (TH1D*) histoMVA->Clone("histo_Higgs");	
  TH1D *histo_EWWZ   = (TH1D*) histoMVA->Clone("histo_EWWZ");	
  TH1D *histo_Higgs_[11];

  for(int i = 0; i<11 ; i++)
	{
    		histo_Higgs_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]);	 
	}

  char finalStateName[2],effMName[10],effEName[10],momMName[10],momEName[10];
  sprintf(effMName,"CMS_eff_m");sprintf(momMName,"CMS_scale_m");
  sprintf(effEName,"CMS_eff_e");sprintf(momEName,"CMS_scale_e");
  sprintf(finalStateName,"3l");


  TH1D* histo_Zg_CMS_scale_jUp            = 	new TH1D( "histo_Zg_CMS_scale_jUp"  	, "histo_Zg_CMS_scale_jUp"  	, nBinMVA, xbins);	histo_Zg_CMS_scale_jUp  ->Sumw2();
  TH1D* histo_Zg_CMS_scale_jDown          = 	new TH1D( "histo_Zg_CMS_scale_jDown"  	, "histo_Zg_CMS_scale_jDown"  	, nBinMVA, xbins); 	histo_Zg_CMS_scale_jDown  ->Sumw2();
  TH1D* histo_VVV_CMS_scale_jUp           = 	new TH1D( "histo_VVV_CMS_scale_jUp"  	, "histo_VVV_CMS_scale_jUp"  	, nBinMVA, xbins); 	histo_VVV_CMS_scale_jUp  ->Sumw2();
  TH1D* histo_VVV_CMS_scale_jDown         = 	new TH1D( "histo_VVV_CMS_scale_jDown"  	, "histo_VVV_CMS_scale_jDown"  	, nBinMVA, xbins); 	histo_VVV_CMS_scale_jDown  ->Sumw2();
  TH1D* histo_WZ_CMS_scale_jUp            = 	new TH1D( "histo_WZ_CMS_scale_jUp"  	, "histo_WZ_CMS_scale_jUp"  	, nBinMVA, xbins); 	histo_WZ_CMS_scale_jUp  ->Sumw2();
  TH1D* histo_WZ_CMS_scale_jDown          = 	new TH1D( "histo_WZ_CMS_scale_jDown"  	, "histo_WZ_CMS_scale_jDown"  	, nBinMVA, xbins); 	histo_WZ_CMS_scale_jDown  ->Sumw2();
   TH1D* histo_EWWZ_CMS_scale_jUp            = 	new TH1D( "histo_EWWZ_CMS_scale_jUp"  	, "histo_EWWZ_CMS_scale_jUp"  	, nBinMVA, xbins); 	histo_EWWZ_CMS_scale_jUp  ->Sumw2();
  TH1D* histo_EWWZ_CMS_scale_jDown          = 	new TH1D( "histo_EWWZ_CMS_scale_jDown"  , "histo_EWWZ_CMS_scale_jDown"  , nBinMVA, xbins); 	histo_EWWZ_CMS_scale_jDown  ->Sumw2();
  TH1D* histo_ZZ_CMS_scale_jUp            = 	new TH1D( "histo_ZZ_CMS_scale_jUp"  	, "histo_ZZ_CMS_scale_jUp"  	, nBinMVA, xbins); 	histo_ZZ_CMS_scale_jUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_scale_jDown          = 	new TH1D( "histo_ZZ_CMS_scale_jDown"  	, "histo_ZZ_CMS_scale_jDown"  	, nBinMVA, xbins); 	histo_ZZ_CMS_scale_jDown  ->Sumw2();
  TH1D *histo_Higgs_CMS_scale_jUp_[11];
  TH1D *histo_Higgs_CMS_scale_jDown_[11];
  for(int i = 0; i<11 ; i++)
	{
    		histo_Higgs_CMS_scale_jUp_[i] 	= (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_scale_jUp");	 
		histo_Higgs_CMS_scale_jDown_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_scale_jDown");	 
	}

  TH1D* histo_Zg_CMS_res_jUp            = 	new TH1D( "histo_Zg_CMS_res_jUp"  	, "histo_Zg_CMS_res_jUp"  	, nBinMVA, xbins);	histo_Zg_CMS_res_jUp  ->Sumw2();
  TH1D* histo_Zg_CMS_res_jDown          = 	new TH1D( "histo_Zg_CMS_res_jDown"  	, "histo_Zg_CMS_res_jDown"  	, nBinMVA, xbins); 	histo_Zg_CMS_res_jDown  ->Sumw2();
  TH1D* histo_VVV_CMS_res_jUp           = 	new TH1D( "histo_VVV_CMS_res_jUp"  	, "histo_VVV_CMS_res_jUp"  	, nBinMVA, xbins); 	histo_VVV_CMS_res_jUp  ->Sumw2();
  TH1D* histo_VVV_CMS_res_jDown         = 	new TH1D( "histo_VVV_CMS_res_jDown"  	, "histo_VVV_CMS_res_jDown"  	, nBinMVA, xbins); 	histo_VVV_CMS_res_jDown  ->Sumw2();
  TH1D* histo_WZ_CMS_res_jUp            = 	new TH1D( "histo_WZ_CMS_res_jUp"  	, "histo_WZ_CMS_res_jUp"  	, nBinMVA, xbins);     	histo_WZ_CMS_res_jUp  ->Sumw2();
  TH1D* histo_WZ_CMS_res_jDown          = 	new TH1D( "histo_WZ_CMS_res_jDown"  	, "histo_WZ_CMS_res_jDown"  	, nBinMVA, xbins); 	histo_WZ_CMS_res_jDown  ->Sumw2();
  TH1D* histo_EWWZ_CMS_res_jUp            = 	new TH1D( "histo_EWWZ_CMS_res_jUp"  	, "histo_EWWZ_CMS_res_jUp"  	, nBinMVA, xbins);     	histo_EWWZ_CMS_res_jUp  ->Sumw2();
  TH1D* histo_EWWZ_CMS_res_jDown          = 	new TH1D( "histo_EWWZ_CMS_res_jDown"    , "histo_EWWZ_CMS_res_jDown"    , nBinMVA, xbins); 	histo_EWWZ_CMS_res_jDown  ->Sumw2();
  TH1D* histo_ZZ_CMS_res_jUp            = 	new TH1D( "histo_ZZ_CMS_res_jUp"  	, "histo_ZZ_CMS_res_jUp"  	, nBinMVA, xbins);     	histo_ZZ_CMS_res_jUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_res_jDown          = 	new TH1D( "histo_ZZ_CMS_res_jDown"  	, "histo_ZZ_CMS_res_jDown"  	, nBinMVA, xbins); 	histo_ZZ_CMS_res_jDown  ->Sumw2();
  TH1D *histo_Higgs_CMS_res_jUp_[11];
  TH1D *histo_Higgs_CMS_res_jDown_[11];
  for(int i = 0; i<11 ; i++)
	{
    		histo_Higgs_CMS_res_jUp_[i] 	= (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_res_jUp");	 
		histo_Higgs_CMS_res_jDown_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_res_jDown");	 
	}

  TH1D* histo_Zg_CMS_scale_lUp            = 	new TH1D( "histo_Zg_CMS_scale_lUp" 	, "histo_Zg_CMS_scale_lUp"  	, nBinMVA, xbins);	histo_Zg_CMS_scale_lUp  ->Sumw2();
  TH1D* histo_Zg_CMS_scale_lDown          = 	new TH1D( "histo_Zg_CMS_scale_lDown"  	, "histo_Zg_CMS_scale_lDown"  	, nBinMVA, xbins); 	histo_Zg_CMS_scale_lDown  ->Sumw2();
  TH1D* histo_VVV_CMS_scale_lUp           = 	new TH1D( "histo_VVV_CMS_scale_lUp"  	, "histo_VVV_CMS_scale_lUp"  	, nBinMVA, xbins); 	histo_VVV_CMS_scale_lUp  ->Sumw2();
  TH1D* histo_VVV_CMS_scale_lDown         = 	new TH1D( "histo_VVV_CMS_scale_lDown"  	, "histo_VVV_CMS_scale_lDown"  	, nBinMVA, xbins);	histo_VVV_CMS_scale_lDown  ->Sumw2();
  TH1D* histo_WZ_CMS_scale_lUp            = 	new TH1D( "histo_WZ_CMS_scale_lUp"  	, "histo_WZ_CMS_scale_lUp"  	, nBinMVA, xbins); 	histo_WZ_CMS_scale_lUp  ->Sumw2();
  TH1D* histo_WZ_CMS_scale_lDown          = 	new TH1D( "histo_WZ_CMS_scale_lDown"  	, "histo_WZ_CMS_scale_lDown"  	, nBinMVA, xbins); 	histo_WZ_CMS_scale_lDown  ->Sumw2();
  TH1D* histo_EWWZ_CMS_scale_lUp            = 	new TH1D( "histo_EWWZ_CMS_scale_lUp"  	, "histo_EWWZ_CMS_scale_lUp"  	, nBinMVA, xbins);   	histo_EWWZ_CMS_scale_lUp  ->Sumw2();
  TH1D* histo_EWWZ_CMS_scale_lDown          = 	new TH1D( "histo_EWWZ_CMS_scale_lDown"  , "histo_EWWZ_CMS_scale_lDown"  , nBinMVA, xbins); 	histo_EWWZ_CMS_scale_lDown  ->Sumw2();
  TH1D* histo_ZZ_CMS_scale_lUp            = 	new TH1D( "histo_ZZ_CMS_scale_lUp"  	, "histo_ZZ_CMS_scale_lUp"  	, nBinMVA, xbins);    	histo_ZZ_CMS_scale_lUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_scale_lDown          = 	new TH1D( "histo_ZZ_CMS_scale_lDown"  	, "histo_ZZ_CMS_scale_lDown"  	, nBinMVA, xbins); 	histo_ZZ_CMS_scale_lDown  ->Sumw2();
  TH1D *histo_Higgs_CMS_scale_lUp_[11];
  TH1D *histo_Higgs_CMS_scale_lDown_[11];

  for(int i = 0; i<11 ; i++)
	{
    		histo_Higgs_CMS_scale_lUp_[i] 	= (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_scale_lUp");	 
		histo_Higgs_CMS_scale_lDown_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_scale_lDown");	 
	}
  TH1D* histo_WZ_CMS_ewkscaleUp           = 		new TH1D( "histo_WZ_CMS_ewkscaleUp"  		, "histo_WZ_CMS_ewkscaleUp"  		, nBinMVA, xbins); 	histo_WZ_CMS_ewkscaleUp  ->Sumw2();
  TH1D* histo_WZ_CMS_ewkscaleDown         = 		new TH1D( "histo_WZ_CMS_ewkscaleDown"		, "histo_WZ_CMS_ewkscaleDown"  		, nBinMVA, xbins); 	histo_WZ_CMS_ewkscaleDown  ->Sumw2();
  TH1D* histo_WZ_CMS_QCDunc_scaleUp       = 		new TH1D( "histo_WZ_CMS_QCDunc_scaleUp"  	, "histo_WZ_CMS_QCDunc_scaleUp"  	, nBinMVA, xbins); 	histo_WZ_CMS_QCDunc_scaleUp  ->Sumw2();
  TH1D* histo_WZ_CMS_QCDunc_scaleDown     = 		new TH1D( "histo_WZ_CMS_QCDunc_scaleDown"  	, "histo_WZ_CMS_QCDunc_scaleDown"  	, nBinMVA, xbins);	histo_WZ_CMS_QCDunc_scaleDown  ->Sumw2();



  

  TH1D* histo_Diff = new TH1D("dummy", "dummy",1000,-1,1); histo_Diff->Sumw2();

  double histo_WZ_CMS_QCDScaleInitial[7] = {0,0,0,0,0,0,0};
  TH1D* histo_Zg_CMS_QCDScaleBounding[6];
  TH1D* histo_VVV_CMS_QCDScaleBounding[6];
  TH1D* histo_WZ_CMS_QCDScaleBounding[6];
  TH1D* histo_EWWZ_CMS_QCDScaleBounding[6];
  TH1D* histo_ZZ_CMS_QCDScaleBounding[6];
  TH1D* histo_Higgs_CMS_QCDScaleBounding[6];
  TH1D* histo_Higgs_CMS_QCDScaleBoundingPerMass[11][6];
  

  for(int nb=0; nb<6; nb++){
    histo_Zg_CMS_QCDScaleBounding[nb]      = new TH1D(Form("histo_Zg_QCDScale_f%d",nb),      Form("histo_Zg_QCDScale_f%d",nb),nBinMVA, xbins);      histo_Zg_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_VVV_CMS_QCDScaleBounding[nb]     = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),     Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WZ_CMS_QCDScaleBounding[nb]	   = new TH1D(Form("histo_WZ_QCDScale_f%d",nb),      Form("histo_WZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_WZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_EWWZ_CMS_QCDScaleBounding[nb]	   = new TH1D(Form("histo_EWWZ_QCDScale_f%d",nb),    Form("histo_EWWZ_QCDScale_f%d",nb),nBinMVA, xbins);    histo_EWWZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ZZ_CMS_QCDScaleBounding[nb]	   = new TH1D(Form("histo_ZZ_QCDScale_f%d",nb),      Form("histo_ZZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_ZZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_Higgs_CMS_QCDScaleBounding[nb]   = new TH1D(Form("histo_Higgs_QCDScale_f%d",nb),   Form("histo_Higgs_QCDScale_f%d",nb),nBinMVA, xbins);   histo_Higgs_CMS_QCDScaleBounding[nb]->Sumw2();

    for(int i = 0; i<11 ; i++)
	{
	  histo_Higgs_CMS_QCDScaleBoundingPerMass[i][nb] =  new TH1D(Form("histo_Higgs_"+SignalSuffix[i]+"_QCDScale_f%d",nb),      Form("histo_Higgs_"+SignalSuffix[i]+"_QCDScale_f%d",nb),nBinMVA, xbins);      histo_Higgs_CMS_QCDScaleBoundingPerMass[i][nb]->Sumw2();	 
	}
  }

  TH1D* histo_Zg_CMS_QCDScaleBoundingUp            = 	new TH1D( "histo_Zg_CMS_QCDScale_ZgUp"  	, "histo_Zg_CMS_QCDScale_ZgUp"  	, nBinMVA, xbins);	histo_Zg_CMS_QCDScaleBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_QCDScaleBoundingDown          = 	new TH1D( "histo_Zg_CMS_QCDScale_ZgDown"  	, "histo_Zg_CMS_QCDScale_ZgDown"  	, nBinMVA, xbins); 	histo_Zg_CMS_QCDScaleBoundingDown  ->Sumw2();

  TH1D* histo_VVV_CMS_QCDScaleBoundingUp           = 	new TH1D( "histo_VVV_CMS_QCDScale_VVVUp"  	, "histo_VVV_CMS_QCDScale_VVVUp"  	, nBinMVA, xbins);      histo_VVV_CMS_QCDScaleBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_QCDScaleBoundingDown         = 	new TH1D( "histo_VVV_CMS_QCDScale_VVVDown"  	, "histo_VVV_CMS_QCDScale_VVVDown"  	, nBinMVA, xbins); 	histo_VVV_CMS_QCDScaleBoundingDown  ->Sumw2();

  TH1D* histo_WZ_CMS_QCDScaleBoundingUp            = 	new TH1D( "histo_WZ_CMS_QCDScale_WZUp"  	, "histo_WZ_CMS_QCDScale_WZUp"  	, nBinMVA, xbins); 	histo_WZ_CMS_QCDScaleBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_QCDScaleBoundingDown          = 	new TH1D( "histo_WZ_CMS_QCDScale_WZDown"  	, "histo_WZ_CMS_QCDScale_WZDown"  	, nBinMVA, xbins); 	histo_WZ_CMS_QCDScaleBoundingDown  ->Sumw2();

  TH1D* histo_EWWZ_CMS_QCDScaleBoundingUp          =   	new TH1D( "histo_EWWZ_CMS_QCDScale_EWWZUp"  	, "histo_EWWZ_CMS_QCDScale_EWWZUp"  	, nBinMVA, xbins); 	histo_EWWZ_CMS_QCDScaleBoundingUp  ->Sumw2();
  TH1D* histo_EWWZ_CMS_QCDScaleBoundingDown        = 	new TH1D( "histo_EWWZ_CMS_QCDScale_EWWZDown"    , "histo_EWWZ_CMS_QCDScale_EWWZDown"  	, nBinMVA, xbins); 	histo_EWWZ_CMS_QCDScaleBoundingDown  ->Sumw2();

  TH1D* histo_ZZ_CMS_QCDScaleBoundingUp            =   	new TH1D( "histo_ZZ_CMS_QCDScale_ZZUp"  	, "histo_ZZ_CMS_QCDScale_ZZUp"  	, nBinMVA, xbins); 	histo_ZZ_CMS_QCDScaleBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_QCDScaleBoundingDown          = 	new TH1D( "histo_ZZ_CMS_QCDScale_ZZDown"  	, "histo_ZZ_CMS_QCDScale_ZZDown"  	, nBinMVA, xbins); 	histo_ZZ_CMS_QCDScaleBoundingDown  ->Sumw2();

  TH1D* histo_Higgs_CMS_QCDScaleBoundingUp         =   	new TH1D( "histo_Higgs_CMS_QCDScaleBoundingUp"  , "histo_Higgs_CMS_QCDScaleBoundingUp"  , nBinMVA, xbins); 	histo_Higgs_CMS_QCDScaleBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_QCDScaleBoundingDown       = 	new TH1D( "histo_Higgs_CMS_QCDScaleBoundingDown"  , "histo_Higgs_CMS_QCDScaleBoundingDown"  , nBinMVA, xbins); 	histo_Higgs_CMS_QCDScaleBoundingDown  ->Sumw2();

  TH1D *histo_Higgs_CMS_QCDScaleBoundingUp_[11];
  TH1D *histo_Higgs_CMS_QCDScaleBoundingDown_[11];

  for(int i = 0; i<11 ; i++)
	{
	  histo_Higgs_CMS_QCDScaleBoundingUp_[i]   = new TH1D( "histo_Higgs_"+SignalSuffix[i]+"_CMS_QCDScale_HiggsUp"  		, "histo_Higgs_"+SignalSuffix[i]+"_CMS_QCDScale_HiggsUp"  , nBinMVA, xbins); 		histo_Higgs_CMS_QCDScaleBoundingUp_[i]->Sumw2();	 
	  histo_Higgs_CMS_QCDScaleBoundingDown_[i] = new TH1D( "histo_Higgs_"+SignalSuffix[i]+"_CMS_QCDScale_HiggsDown"  	, "histo_Higgs_"+SignalSuffix[i]+"_CMS_QCDScale_HiggsDown", nBinMVA, xbins); 		histo_Higgs_CMS_QCDScaleBoundingDown_[i]->Sumw2();
	}


  TH1D* histo_Zg_CMS_PDFBounding[102];
  TH1D* histo_VVV_CMS_PDFBounding[102];
  TH1D* histo_WZ_CMS_PDFBounding[102];
  TH1D* histo_EWWZ_CMS_PDFBounding[102];
  TH1D* histo_ZZ_CMS_PDFBounding[102];
  TH1D* histo_Higgs_CMS_PDFBoundingPerMass[11][102];
  TH1D* histo_Higgs_CMS_PDFBounding[102];

  for(int nb=0; nb<102; nb++){
    histo_Zg_CMS_PDFBounding[nb]      = new TH1D(Form("histo_Zg_PDF_f%d",nb) ,      Form("histo_Zg_PDF_f%d",nb),     	nBinMVA, xbins); histo_Zg_CMS_PDFBounding[nb]->Sumw2();
    histo_VVV_CMS_PDFBounding[nb]     = new TH1D(Form("histo_VVV_PDF_f%d",nb),      Form("histo_VVV_PDF_f%d",nb),    	nBinMVA, xbins); histo_VVV_CMS_PDFBounding[nb]->Sumw2();
    histo_WZ_CMS_PDFBounding[nb]      = new TH1D(Form("histo_WZ_PDF_f%d",nb) ,      Form("histo_WZ_PDF_f%d",nb),     	nBinMVA, xbins); histo_WZ_CMS_PDFBounding[nb]->Sumw2();
    histo_EWWZ_CMS_PDFBounding[nb]    = new TH1D(Form("histo_EWWZ_PDF_f%d",nb) ,    Form("histo_EWWZ_PDF_f%d",nb),     	nBinMVA, xbins); histo_EWWZ_CMS_PDFBounding[nb]->Sumw2();
    histo_ZZ_CMS_PDFBounding[nb]      = new TH1D(Form("histo_ZZ_PDF_f%d",nb) ,      Form("histo_ZZ_PDF_f%d",nb),     	nBinMVA, xbins); histo_ZZ_CMS_PDFBounding[nb]->Sumw2();
   
    histo_Higgs_CMS_PDFBounding[nb]   = new TH1D(Form("histo_Higgs_PDF_f%d",nb),    Form("histo_Higgs_PDF_f%d",nb),     nBinMVA, xbins); histo_Higgs_CMS_PDFBounding[nb]->Sumw2();

    for(int i = 0; i<11 ; i++)
	{
	  histo_Higgs_CMS_PDFBoundingPerMass[i][nb] =  new TH1D(Form("histo_Higgs_"+SignalSuffix[i]+"_PDF_f%d",nb),      Form("histo_Higgs_"+SignalSuffix[i]+"_PDF_f%d",nb),nBinMVA, xbins);       histo_Higgs_CMS_PDFBoundingPerMass[i][nb]->Sumw2();	 
	}
  }

  TH1D* histo_Zg_CMS_pdf_qqbarUp            = 	new TH1D( "histo_Zg_CMS_pdf_qqbarUp"  	, "histo_Zg_CMS_pdf_qqbarUp"  		, nBinMVA, xbins);	histo_Zg_CMS_pdf_qqbarUp  ->Sumw2();
  TH1D* histo_Zg_CMS_pdf_qqbarDown          = 	new TH1D( "histo_Zg_CMS_pdf_qqbarDown"  , "histo_Zg_CMS_pdf_qqbarDown"  	, nBinMVA, xbins); 	histo_Zg_CMS_pdf_qqbarDown  ->Sumw2();
  TH1D* histo_VVV_CMS_pdf_qqbarUp           = 	new TH1D( "histo_VVV_CMS_pdf_qqbarUp"  	, "histo_VVV_CMS_pdf_qqbarUp" 	 	, nBinMVA, xbins); 	histo_VVV_CMS_pdf_qqbarUp  ->Sumw2();
  TH1D* histo_VVV_CMS_pdf_qqbarDown         = 	new TH1D( "histo_VVV_CMS_pdf_qqbarDown"  , "histo_VVV_CMS_pdf_qqbarDown" 	, nBinMVA, xbins); 	histo_VVV_CMS_pdf_qqbarDown  ->Sumw2();
  TH1D* histo_WZ_CMS_pdf_qqbarUp            = 	new TH1D( "histo_WZ_CMS_pdf_qqbarUp"  	, "histo_WZ_CMS_pdf_qqbarUp"  		, nBinMVA, xbins); 	histo_WZ_CMS_pdf_qqbarUp  ->Sumw2();
  TH1D* histo_WZ_CMS_pdf_qqbarDown          = 	new TH1D( "histo_WZ_CMS_pdf_qqbarDown"  , "histo_WZ_CMS_pdf_qqbarDown"  	, nBinMVA, xbins); 	histo_WZ_CMS_pdf_qqbarDown  ->Sumw2();
  TH1D* histo_EWWZ_CMS_pdf_qqbarUp          = 	new TH1D( "histo_EWWZ_CMS_pdf_qqbarUp"  , "histo_EWWZ_CMS_pdf_qqbarUp"  	, nBinMVA, xbins); 	histo_EWWZ_CMS_pdf_qqbarUp  ->Sumw2();
  TH1D* histo_EWWZ_CMS_pdf_qqbarDown        =   new TH1D( "histo_EWWZ_CMS_pdf_qqbarDown"  , "histo_EWWZ_CMS_pdf_qqbarDown"  	, nBinMVA, xbins); 	histo_EWWZ_CMS_pdf_qqbarDown  ->Sumw2();

  TH1D* histo_ZZ_CMS_pdf_qqbarUp           =   	new TH1D( "histo_ZZ_CMS_pdf_qqbarUp"  	, "histo_ZZ_CMS_pdf_qqbarUp"  		, nBinMVA, xbins); 	histo_ZZ_CMS_pdf_qqbarUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_pdf_qqbarDown         = 	new TH1D( "histo_ZZ_CMS_pdf_qqbarDown"  , "histo_ZZ_CMS_pdf_qqbarDown"  	, nBinMVA, xbins); 	histo_ZZ_CMS_pdf_qqbarDown  ->Sumw2();

  TH1D* histo_Higgs_CMS_pdf_qqbarUp       = 	new TH1D( "histo_Higgs_CMS_pdf_qqbarUp"  	, "histo_Higgs_CMS_pdf_qqbarUp"  	, nBinMVA, xbins);      histo_Higgs_CMS_pdf_qqbarUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_pdf_qqbarDown     = 	new TH1D( "histo_Higgs_CMS_pdf_qqbarDown"  	, "histo_Higgs_CMS_pdf_qqbarDown"  	, nBinMVA, xbins);	 histo_Higgs_CMS_pdf_qqbarDown  ->Sumw2();

  TH1D *histo_Higgs_CMS_pdf_qqbarUp_[11];
  TH1D *histo_Higgs_CMS_pdf_qqbarDown_[11];

  for(int i = 0; i<11 ; i++)
	{
	  histo_Higgs_CMS_pdf_qqbarUp_[i]   = new TH1D( "histo_Higgs_"+SignalSuffix[i]+"_CMS_pdf_qqbarUp"  	, "histo_Higgs_"+SignalSuffix[i]+"_CMS_pdf_qqbarUp"  , nBinMVA, xbins); 		histo_Higgs_CMS_pdf_qqbarUp_[i]->Sumw2();	 
	  histo_Higgs_CMS_pdf_qqbarDown_[i] = new TH1D( "histo_Higgs_"+SignalSuffix[i]+"_CMS_pdf_qqbarDown"  	, "histo_Higgs_"+SignalSuffix[i]+"_CMS_pdf_qqbarDown"  , nBinMVA, xbins); 		histo_Higgs_CMS_pdf_qqbarDown_[i]->Sumw2();
	}

  
  TH1D* histo_Zg_CMS_MVALepEffMBoundingUp        = new TH1D( Form("histo_Zg_%sUp",effMName)  , Form("histo_Zg_%sUp",effMName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVALepEffMBoundingDown      = new TH1D( Form("histo_Zg_%sDown",effMName), Form("histo_Zg_%sDown",effMName), nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effMName)  , Form("histo_VVV_%sUp",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effMName), Form("histo_VVV_%sDown",effMName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_WZ_%sUp",effMName)  , Form("histo_WZ_%sUp",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_WZ_%sDown",effMName), Form("histo_WZ_%sDown",effMName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingDown->Sumw2();

  TH1D* histo_EWWZ_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_EWWZ_%sUp",effMName)  , Form("histo_EWWZ_%sUp",effMName)  , nBinMVA, xbins); histo_EWWZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_EWWZ_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_EWWZ_%sDown",effMName), Form("histo_EWWZ_%sDown",effMName), nBinMVA, xbins); histo_EWWZ_CMS_MVALepEffMBoundingDown->Sumw2();

  TH1D* histo_ZZ_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_ZZ_%sUp",effMName)  , Form("histo_ZZ_%sUp",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_ZZ_%sDown",effMName), Form("histo_ZZ_%sDown",effMName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingDown->Sumw2();
 TH1D* histo_Higgs_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_Higgs_%sUp",effMName)  , Form("histo_Higgs_%sUp",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_Higgs_%sDown",effMName), Form("histo_Higgs_%sDown",effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_Zg_CMS_MVALepEffMBoundingAvg       = new TH1D( Form("histo_Zg_%sAvg",effMName)  , Form("histo_Zg_%sAvg",effMName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effMName)  , Form("histo_VVV_%sAvg",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_WZ_%sAvg",effMName)  , Form("histo_WZ_%sAvg",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_EWWZ_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_EWWZ_%sAvg",effMName)  , Form("histo_EWWZ_%sAvg",effMName)  , nBinMVA, xbins); histo_EWWZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_ZZ_%sAvg",effMName)  , Form("histo_ZZ_%sAvg",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_Higgs_%sAvg",effMName)  , Form("histo_Higgs_%sAvg",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingAvg  ->Sumw2();


  TH1D* histo_Higgs_CMS_MVALepEffMperMassBoundingDown[11];
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
  
  //for loop over samples
  for(unsigned int ifile=0; ifile<1./*infilenamev.size()*/; ifile++){
  cout<<infilenamev[ifile]<<endl;
  TFile *the_input_file = TFile::Open((TString)infilenamev[ifile]);
  TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("Events");
  TTree *the_input_runs = (TTree*)the_input_file->FindObjectAny("Runs");


  //cross section weight
  float xSecLumiWeight= 1.;
  Double_t genEventCountn;
  Double_t genEventCount=0;
  if(infilecatv[ifile] != 0){
    the_input_runs->SetBranchAddress("genEventSumw",&genEventCountn);
    unsigned int numberOfSubSamples = int(the_input_runs->GetEntries());
    cout<<"numberOfSubSamples: "<<numberOfSubSamples<<endl;
    for(unsigned int n=0; n<numberOfSubSamples; n++){
      the_input_runs->GetEntry(n);
      genEventCount = genEventCount + genEventCountn;
    }
    xSecLumiWeight=infilexsecv[ifile]*lumi/genEventCount;
  }

 int numberOfEvents = the_input_tree->GetEntriesFast();
 cout<<"number of Events: "<<numberOfEvents<<endl;
 cout<<"xSecLumiWeight: "<<xSecLumiWeight<<endl;

 UInt_t  nElectron;
 Float_t Electron_pt[50];
 Float_t Electron_eta[50];
 Float_t Electron_phi[50]; 
 Float_t Electron_mass[50]; 
 Int_t   Electron_pdgId[50]; 
 Int_t   Electron_cutBased[50]; 
 Int_t   Electron_cutBasedHLTPreSel[50];
 UChar_t Electron_genPartFlav[50];
 Float_t Electron_dxy[50]; 
 Float_t Electron_dz[50]; 
 the_input_tree->SetBranchAddress("nElectron",&nElectron); 
 the_input_tree->SetBranchAddress("Electron_pt",Electron_pt); 
 the_input_tree->SetBranchAddress("Electron_eta",Electron_eta); 
 the_input_tree->SetBranchAddress("Electron_phi",Electron_phi); 
 the_input_tree->SetBranchAddress("Electron_mass",Electron_mass); 
 the_input_tree->SetBranchAddress("Electron_pdgId",Electron_pdgId); 
 the_input_tree->SetBranchAddress("Electron_cutBased",Electron_cutBased); 
 the_input_tree->SetBranchAddress("Electron_cutBased_HLTPreSel",Electron_cutBasedHLTPreSel); 
 the_input_tree->SetBranchAddress("Electron_dxy",Electron_dxy); 
 the_input_tree->SetBranchAddress("Electron_dz",Electron_dz); 

 UInt_t  nMuon;
 Float_t Muon_pt[50];
 Float_t Muon_eta[50];
 Float_t Muon_phi[50]; 
 Float_t Muon_mass[50];
 Int_t   Muon_pdgId[50]; 
 Bool_t  Muon_tightId[50];  
 Float_t Muon_pfRelIso04_all[50]; 
 UChar_t Muon_genPartFlav[50];
 Float_t Muon_dxy[50]; 
 Float_t Muon_dz[50]; 
 the_input_tree->SetBranchAddress("nMuon",&nMuon); 
 the_input_tree->SetBranchAddress("Muon_pt",Muon_pt); 
 the_input_tree->SetBranchAddress("Muon_eta",Muon_eta); 
 the_input_tree->SetBranchAddress("Muon_phi",Muon_phi); 
 the_input_tree->SetBranchAddress("Muon_mass",Muon_mass); 
 the_input_tree->SetBranchAddress("Muon_pdgId",Muon_pdgId); 
 the_input_tree->SetBranchAddress("Muon_tightId",Muon_tightId); 
 the_input_tree->SetBranchAddress("Muon_pfRelIso04_all",Muon_pfRelIso04_all); 
 the_input_tree->SetBranchAddress("Muon_dxy",Muon_dxy); 
 the_input_tree->SetBranchAddress("Muon_dz",Muon_dz); 
 

 UInt_t  nJet;
 Float_t Jet_pt[200];
 Float_t Jet_pt_jesUp[200];
 Float_t Jet_pt_jesDown[200];
 Float_t Jet_pt_jerUp[200];
 Float_t Jet_pt_jerDown[200];
 Float_t Jet_eta[200];
 Float_t Jet_phi[200];
 Float_t Jet_mass[200];
 Float_t Jet_mass_jesUp[200];
 Float_t Jet_mass_jesDown[200];
 Float_t Jet_mass_jerUp[200];
 Float_t Jet_mass_jerDown[200];
 Float_t Jet_btagCSVV2[200];
 the_input_tree->SetBranchAddress("nJet"		,&nJet); 
 the_input_tree->SetBranchAddress("Jet_pt"		,Jet_pt);
 the_input_tree->SetBranchAddress("Jet_pt_jesTotalUp"	,Jet_pt_jesUp);
 the_input_tree->SetBranchAddress("Jet_pt_jesTotalDown"	,Jet_pt_jesDown);
 the_input_tree->SetBranchAddress("Jet_pt_jerUp"	,Jet_pt_jerUp);
 the_input_tree->SetBranchAddress("Jet_pt_jerDown"	,Jet_pt_jerDown);
 the_input_tree->SetBranchAddress("Jet_eta"		,Jet_eta); 
 the_input_tree->SetBranchAddress("Jet_phi"		,Jet_phi); 
 the_input_tree->SetBranchAddress("Jet_mass"		,Jet_mass); 
 the_input_tree->SetBranchAddress("Jet_mass_jesTotalUp"	,Jet_mass_jesUp); 
 the_input_tree->SetBranchAddress("Jet_mass_jesTotalDown",Jet_mass_jesDown); 
 the_input_tree->SetBranchAddress("Jet_mass_jerUp"	,Jet_mass_jerUp); 
 the_input_tree->SetBranchAddress("Jet_mass_jerDown"	,Jet_mass_jerDown); 
 the_input_tree->SetBranchAddress("Jet_btagCSVV2"	,Jet_btagCSVV2);

 Float_t MET_pt; 
 Float_t MET_pt_jerDown; 
 Float_t MET_pt_jerUp;
 Float_t MET_pt_jesDown; 
 Float_t MET_pt_jesUp; 
 Float_t MET_phi; 
 Float_t MET_phi_jerDown; 
 Float_t MET_phi_jerUp;
 Float_t MET_phi_jesDown; 
 Float_t MET_phi_jesUp; 
 the_input_tree->SetBranchAddress("MET_pt",		&MET_pt); 
 the_input_tree->SetBranchAddress("MET_pt_jerUp",	&MET_pt_jerUp); 
 the_input_tree->SetBranchAddress("MET_pt_jerDown",	&MET_pt_jerDown); 
 the_input_tree->SetBranchAddress("MET_pt_jesTotalUp",	&MET_pt_jesUp); 
 the_input_tree->SetBranchAddress("MET_pt_jesTotalDown",&MET_pt_jesDown); 
 the_input_tree->SetBranchAddress("MET_phi",		&MET_phi); 
 the_input_tree->SetBranchAddress("MET_phi_jerUp",	&MET_phi_jerUp); 
 the_input_tree->SetBranchAddress("MET_phi_jerDown",	&MET_phi_jerDown); 
 the_input_tree->SetBranchAddress("MET_phi_jesTotalUp",	&MET_phi_jesUp); 
 the_input_tree->SetBranchAddress("MET_phi_jesTotalDown",&MET_phi_jesDown); 
 
 Int_t PV_npvs;
 the_input_tree->SetBranchAddress("PV_npvs",&PV_npvs); 

 Bool_t hltIsoMu24;
 the_input_tree->SetBranchAddress("HLT_IsoMu24",&hltIsoMu24); 
 Bool_t hltMu17Mu8;
 the_input_tree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",&hltMu17Mu8);
 Bool_t hltEle27;
 the_input_tree->SetBranchAddress("HLT_Ele27_WPTight_Gsf",&hltEle27); 
 Bool_t hltEle35;
 the_input_tree->SetBranchAddress("HLT_Ele35_WPTight_Gsf",&hltEle35);
 Bool_t hltEle23Ele12;
 the_input_tree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",&hltEle23Ele12);

 float genWeight;
 float puWeight;
 if(infilecatv[ifile] != 0){
   the_input_tree->SetBranchAddress("genWeight",&genWeight);
   the_input_tree->SetBranchAddress("puWeight",&puWeight);
   the_input_tree->SetBranchAddress("Muon_genPartFlav",Muon_genPartFlav);
   the_input_tree->SetBranchAddress("Electron_genPartFlav",Electron_genPartFlav);  
 }

//for loop over events for given sample
 cout<<"number of Events: "<<numberOfEvents<<endl;
 
 for (int i=0; i < numberOfEvents; ++i) 
     {
       Bool_t passFilter[11] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};

       the_input_tree->GetEntry(i);
     
       if(i%100000==0) cout<<"Event: "<<i<<endl;
      
       if(!(hltIsoMu24 || hltEle27 || hltEle35 || hltEle23Ele12 || hltMu17Mu8)) continue;

       vector<int> idLep;  vector<int> idLepPdg; vector<int> idTight; unsigned int goodIsTight = 0; vector<TLorentzVector> leptons; int goodIsGenLep = 0; vector<int> isGenLep;
      
	for(unsigned int e=0; e<nElectron; e++)
	  {
	    if(Electron_pt[e]>10. && TMath::Abs(Electron_eta[e])<2.5 && ((TMath::Abs(Electron_eta[e])<1.5 && TMath::Abs(Electron_dxy[e]) < 0.05 && TMath::Abs(Electron_dz[e]) < 0.1 ) || (TMath::Abs(Electron_eta[e])>=1.5 && TMath::Abs(Electron_dxy[e]) < 0.1 && TMath::Abs(Electron_dz[e]) < 0.2)))
	      { 
				if(Electron_cutBasedHLTPreSel[e]!=1) {checkTheID->Fill(Electron_cutBased[e]);}
				if(Electron_cutBasedHLTPreSel[e]==1) {checkTheID->Fill(Electron_cutBased[e]+5.);}
			  	if(Electron_cutBased[e]==4 && Electron_cutBasedHLTPreSel[e]==1) 	{idLep.push_back(e); idTight.push_back(1); goodIsTight++ ; idLepPdg.push_back(Electron_pdgId[e]);}
	  			else if(Electron_cutBasedHLTPreSel[e]==1) 				{idLep.push_back(e); idTight.push_back(0);               ; idLepPdg.push_back(Electron_pdgId[e]);}
	  			else continue;
			  if     ( (Electron_genPartFlav[e] == 1 || Electron_genPartFlav[e] == 15 ) && infilecatv[ifile]!=0) {goodIsGenLep++;isGenLep.push_back(1);}
			  else if(!(Electron_genPartFlav[e] == 1 || Electron_genPartFlav[e] == 15 ) && infilecatv[ifile]!=0) {isGenLep.push_back(0);}
	  			TLorentzVector lepTemp; lepTemp.SetPtEtaPhiM(Electron_pt[e],Electron_eta[e],Electron_phi[e],Electron_mass[e]);
	  			leptons.push_back(lepTemp);	
			}
      }

      for(unsigned int m=0; m<nMuon; m++)
      	{	  
		if(Muon_pt[m]>10. && TMath::Abs(Muon_eta[m])<2.4 && TMath::Abs(Muon_dxy[m]) < 0.02 && TMath::Abs(Muon_dz[m]) < 0.1)
	  		{
	  			if     (Muon_tightId[m] && Muon_pfRelIso04_all[m] < 0.15) 	{idLep.push_back(m); idTight.push_back(1); goodIsTight++ ; idLepPdg.push_back(Muon_pdgId[m]);}
	  			else if(Muon_tightId[m] && Muon_pfRelIso04_all[m] < 0.4) 	{idLep.push_back(m); idTight.push_back(0);               ; idLepPdg.push_back(Muon_pdgId[m]);}
	  			else continue;
				if     ( (Muon_genPartFlav[m] == 1 || Muon_genPartFlav[m] == 15 ) && infilecatv[ifile]!=0) {goodIsGenLep++;isGenLep.push_back(1);}
				else if(!(Muon_genPartFlav[m] == 1 || Muon_genPartFlav[m] == 15 ) && infilecatv[ifile]!=0) {isGenLep.push_back(0);}
	  			TLorentzVector lepTemp; lepTemp.SetPtEtaPhiM(Muon_pt[m],Muon_eta[m],Muon_phi[m],Muon_mass[m]);
	  			leptons.push_back(lepTemp);	
			}
      }

      if(leptons.size() != 3) continue;
      
      bool passLepPtCutScaleUp	=false;
      bool passLepPtCutScaleDown=false;
      bool passLepPtCut		=false;
      if(  leptons.at(0).Pt() > 20. 	   && leptons.at(1).Pt() > 20. 	      && leptons.at(2).Pt() > 20.  ) 	    passLepPtCut=true;
      if(  leptons.at(0).Pt() > (20.*0.99) && leptons.at(1).Pt() > (20.*0.99) && leptons.at(2).Pt() > (20.*0.99)  ) passLepPtCutScaleDown=true;
      if(  leptons.at(0).Pt() > (20.*1.01) && leptons.at(1).Pt() > (20.*1.01) && leptons.at(2).Pt() > (20.*1.01)  ) passLepPtCutScaleUp=true;

       if(!(passLepPtCut || passLepPtCutScaleUp || passLepPtCutScaleDown)) continue;
      //if(!(leptons.at(0).Pt() > 25. || leptons.at(1).Pt() > 25. || leptons.at(2).Pt() > 25.) ) continue;
      

      double systTotLep[2] = {1.0, 1.0}; // m/e
      for(unsigned nl=0; nl<idLep.size(); nl++){
        if(TMath::Abs((int)idLepPdg[nl]) == 13) systTotLep[0] = systTotLep[0] * 1.015;
        if(TMath::Abs((int)idLepPdg[nl]) == 11) systTotLep[1] = systTotLep[1] * 1.02;            
      }


      //Find Z boson, categorize events
      double minMassll = 999.0;
      double minMassZ = 999.0;
      double mass3l = 0.0;
      double eta3l = 0.0;
      double deltaRllMin = 999.0;
      int type3l = 0;
      int tagZ[3] = {-1,-1,-1};
      for(unsigned nl0=0; nl0<idLep.size()-1; nl0++){
        for(unsigned nl1=nl0+1; nl1<idLep.size(); nl1++){
      	  if((int)idLepPdg[nl0] * (int)idLepPdg[nl1] > 0)	 continue;
          TLorentzVector dilepAux( leptons.at(nl0) + (leptons.at(nl1)) );
	  
	  if(TMath::Abs((int)idLepPdg[nl0])==TMath::Abs((int)idLepPdg[nl1]) &&
	     TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassZ-91.1876)) {
	     minMassZ = dilepAux.M();tagZ[0]=nl0;tagZ[1]=nl1;
	     
	     if(TMath::Abs((int)idLepPdg[nl0]) == 13) 			  type3l = 0;
	     else                                                         type3l = 1;
	  }
	  if(minMassll > dilepAux.M()) minMassll = dilepAux.M();
        }
      }
      
      for(unsigned nl0=0; nl0<idLep.size(); nl0++){
        if((int)nl0==tagZ[0]||(int)nl0==tagZ[1]) 		 continue;
        tagZ[2] = nl0;
        break;
      }  
      if(tagZ[0] == -1 || tagZ[1] == -1 || tagZ[2] == -1)	 continue;
      if(TMath::Abs((int)idLepPdg[tagZ[2]]) == 13) 		       type3l += 0;
      else							       type3l += 2;
      TLorentzVector trilep(	( ( leptons.at(0)) )
        		       +( ( leptons.at(1)) )
        		       +( ( leptons.at(2)) ));

      mass3l = trilep.M();
      eta3l  = trilep.Eta();

      //GenLeptonMatching
      


      //Jets
      vector<int> idJet;vector<int> idJet50; vector<int> idJetJESUP; vector<int> idJetJESDOWN; vector<int> idJetJERUP; vector<int> idJetJERDOWN; double btagjet[2] = {0., 0.}; vector<TLorentzVector> jets,jetsJesUp,jetsJesDown,jetsJerUp,jetsJerDown;
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      double theHT = 0;
      //tightCSV2 0.9693

      for(unsigned int nj=0; nj<nJet ; nj++){

        if(Jet_pt[nj] < 10) continue;	
        TLorentzVector jetTemp; 	jetTemp.SetPtEtaPhiM(Jet_pt[nj],Jet_eta[nj],Jet_phi[nj],Jet_mass[nj]);
	TLorentzVector jetTempJesUp; 	jetTempJesUp.SetPtEtaPhiM(Jet_pt_jesUp[nj],Jet_eta[nj],Jet_phi[nj],Jet_mass_jesUp[nj]);
	TLorentzVector jetTempJesDown; 	jetTempJesDown.SetPtEtaPhiM(Jet_pt_jesDown[nj],Jet_eta[nj],Jet_phi[nj],Jet_mass_jesDown[nj]);
	TLorentzVector jetTempJerUp; 	jetTempJerUp.SetPtEtaPhiM(Jet_pt_jerUp[nj],Jet_eta[nj],Jet_phi[nj],Jet_mass_jerUp[nj]);
	TLorentzVector jetTempJerDown; 	jetTempJerDown.SetPtEtaPhiM(Jet_pt_jerDown[nj],Jet_eta[nj],Jet_phi[nj],Jet_mass_jerDown[nj]);
	//jets.push_back(jetTemp);
        
        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(jetTemp.DeltaR(leptons.at(nl)) < 0.4) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

	if(jetTemp.Pt() > 30 && 
	   Jet_btagCSVV2[nj] > bDiscrMax) bDiscrMax = Jet_btagCSVV2[nj];

        if( jetTemp.Pt() >= 50) {
	  theHT = theHT + (jetTemp.Pt());
	  idJet.push_back(nj);
	  jets.push_back(jetTemp);
	}
	if( jetTempJesUp.Pt() >= 50) {
	  jetsJesUp.push_back(jetTempJesUp);
	}
	if( jetTempJesDown.Pt() >= 50) {
	  jetsJesDown.push_back(jetTempJesDown);
	}
	if( jetTempJerUp.Pt() >= 50) {
	  jetsJerUp.push_back(jetTempJerUp);
	}
	if( jetTempJerDown.Pt() >= 50) {
	  jetsJerDown.push_back(jetTempJerDown);
	}
      }

      passFilter[ 5] = minMassll > 4.;
      passFilter[ 6] = MET_pt > 30;
      passFilter[ 7] = minMassZ > minMass && minMassZ < maxMass;
      passFilter[ 8] = leptons.at(tagZ[2]).Pt() > 20;
      passFilter[ 9] = mass3l > 100;
      passFilter[10] = true;
      if(true) passFilter[10] = bDiscrMax < 0.9693;

      bool passAllCuts[1] 		= {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] 	&& passFilter[9] && passFilter[10] && passLepPtCut};
      bool passAllCutsLepScaleUp[1] 	= {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8]   	&& passFilter[9] && passFilter[10] && passLepPtCutScaleUp};
      bool passAllCutsLepScaleDown[1] 	= {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8]   && passFilter[9] && passFilter[10]   && passLepPtCutScaleDown};
      bool passAllCutsJESUP 	= {passFilter[5] &&  MET_pt_jesUp   > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};
      bool passAllCutsJESDOWN   = {passFilter[5] &&  MET_pt_jesDown > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};
      bool passAllCutsJERUP 	= {passFilter[5] &&  MET_pt_jerUp   > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};
      bool passAllCutsJERDOWN   = {passFilter[5] &&  MET_pt_jerDown > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};

      bool fakeEnriched[2]= {passFilter[5] && passFilter[6] && !passFilter[7] && passFilter[8] && passFilter[9] && !passFilter[10] && passLepPtCut,
			     passFilter[5] &&                  !passFilter[7] && passFilter[8] && passFilter[9] && !passFilter[10] && passLepPtCut};
	

     TLorentzVector dijetjesUp,dijetjesDown;
      double dijetDeltaEtaJesUp = 0; double dijetDeltaEtaJesDown = 0;
      if(jetsJesUp.size() >= 2 ) {
	dijetjesUp=jetsJesUp[0]+jetsJesUp[1];				
	dijetDeltaEtaJesUp = TMath::Abs(jetsJesUp[0].Eta()-jetsJesUp[1].Eta());
      }
      if(jetsJesDown.size() >= 2) {
	dijetjesDown=jetsJesDown[0]+jetsJesDown[1];	
        dijetDeltaEtaJesDown = TMath::Abs(jetsJesDown[0].Eta()-jetsJesDown[1].Eta());
      }
      TLorentzVector dijetjerUp,dijetjerDown;
      double dijetDeltaEtaJerUp = 0; double dijetDeltaEtaJerDown = 0;
      if(jetsJerUp.size() >= 2 && infilecatv[ifile]!=0) {
	dijetjerUp=jetsJerUp[0]+jetsJerUp[1];				
	dijetDeltaEtaJerUp = TMath::Abs(jetsJerUp[0].Eta()-jetsJerUp[1].Eta());
      }
      if(jetsJerDown.size() >= 2 && infilecatv[ifile]!=0) {
	dijetjerDown=jetsJerDown[0]+jetsJerDown[1];				
	dijetDeltaEtaJerDown = TMath::Abs(jetsJerDown[0].Eta()-jetsJerDown[1].Eta());
	}

     TLorentzVector metP4; 		metP4.SetPtEtaPhiM(MET_pt,0.,MET_phi,0.);
     TLorentzVector metP4JesUp; 	metP4.SetPtEtaPhiM(MET_pt_jesUp,0.,MET_phi_jesUp,0.);
     TLorentzVector metP4JesDown; 	metP4.SetPtEtaPhiM(MET_pt_jesDown,0.,MET_phi_jesDown,0.);
     TLorentzVector metP4JerUp; 	metP4.SetPtEtaPhiM(MET_pt_jerUp,0.,MET_phi_jerUp,0.);
     TLorentzVector metP4JerDown; 	metP4.SetPtEtaPhiM(MET_pt_jerDown,0.,MET_phi_jerDown,0.);
     TLorentzVector wp4 = leptons.at(tagZ[2]) + metP4;
     TLorentzVector zp4 = leptons.at(tagZ[0]) + leptons.at(tagZ[1]) ;
     double transverseMass2 = TMath::Sqrt( TMath::Power(zp4.Et() + wp4.Et(),2) - TMath::Power(zp4.Px() + wp4.Px(),2) - TMath::Power(zp4.Py() + wp4.Py(),2));

     //JES UNCERTAINTY
     TLorentzVector  wp4JESUP 	= leptons.at(tagZ[2]) + metP4JesUp;
     TLorentzVector  wp4JESDOWN = leptons.at(tagZ[2]) + metP4JesDown;
     TLorentzVector  wp4JERUP 	= leptons.at(tagZ[2]) + metP4JerUp;
     TLorentzVector  wp4JERDOWN = leptons.at(tagZ[2]) + metP4JerDown;
     TLorentzVector wp4LEPScaleUP 	= (leptons.at(tagZ[2])*1.01 	+ metP4);
     TLorentzVector zp4LEPScaleUP 	= (leptons.at(tagZ[0])*1.01 	+ leptons.at(tagZ[1])*1.01);  
     TLorentzVector wp4LEPScaleDOWN 	= (leptons.at(tagZ[2])*0.99 	+ metP4);
     TLorentzVector zp4LEPScaleDOWN 	= (leptons.at(tagZ[0])*0.99 	+ leptons.at(tagZ[1])*0.99);



     TLorentzVector dijet;
     double dijetDeltaEta	= 0;
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
	dijet 		= jets[0]+jets[1];
	dijetDeltaEta 	= TMath::Abs(jets[0].Eta()- jets[1].Eta());
	zep3l 		= TMath::Abs(eta3l - 0.5*(jets[0].Eta() + jets[1].Eta()));

	passVBFLoose 		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>100.;     
      	passVBFLooseJESUP 	= passAllCutsJESUP 		&& jetsJesUp.size()   >= 2 	&& dijetjesUp.M()	>100.;
      	passVBFLooseJESDOWN 	= passAllCutsJESDOWN	 	&& jetsJesDown.size() >= 2 	&& dijetjesDown.M()	>100.;
      	passVBFLooseJERUP 	= passAllCutsJERUP 		&& jetsJerUp.size()   >= 2      && dijetjerUp.M()	>100.;
      	passVBFLooseJERDOWN 	= passAllCutsJERDOWN		&& jetsJerDown.size() >= 2 	&& dijetjerDown.M()	>100.;
      	passVBFLooseLEPScaleUP   = passAllCutsLepScaleUp[0] 	&& idJet.size() >= 2 		&& dijet.M() 		>100.;
      	passVBFLooseLEPScaleDOWN = passAllCutsLepScaleDown[0] 	&& idJet.size() >= 2 		&& dijet.M() 		>100.;

      	passVBF 		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>500.	&& dijetDeltaEta > 2.5;     
        passVBFzep 		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>500.	&& dijetDeltaEta > 2.5	&& zep3l<2.5;     
	//passVBFzep2Jet50 	= passAllCuts[0] 		&& idJet50.size() >= 2 		&&  	(( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M()	>500.	       	&& TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta() - ((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()) > 2.5				&& zep3l<2.5; 
    
      	passVBFJESUP 		= passAllCutsJESUP 		&& jetsJesUp.size()   >= 2 	&& dijetjesUp.M()	>500.	&& dijetDeltaEtaJesUp > 2.5;
      	passVBFJESDOWN 		= passAllCutsJESDOWN	 	&& jetsJesDown.size() >= 2 	&& dijetjesDown.M()	>500.	&& dijetDeltaEtaJesDown > 2.5;
      	passVBFJERUP 		= passAllCutsJERUP 		&& jetsJerUp.size()   >= 2      && dijetjerUp.M()	>500.	&& dijetDeltaEtaJerUp > 2.5;
      	passVBFJERDOWN 		= passAllCutsJERDOWN		&& jetsJerDown.size() >= 2 	&& dijetjerDown.M()	>500.	&& dijetDeltaEtaJerDown > 2.5;
      	passVBFLEPScaleUP  	= passAllCutsLepScaleUp[0] 	&& idJet.size() >= 2 		&& dijet.M() 		>500.	&& dijetDeltaEta > 2.5;
      	passVBFLEPScaleDOWN 	= passAllCutsLepScaleDown[0] 	&& idJet.size() >= 2 		&& dijet.M() 		>500.	&& dijetDeltaEta > 2.5;
	
	passVBFLooseMJJLess500  = passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		<500.;
	passVBFLooseDEtaLess25  = passAllCuts[0] 		&& idJet.size() >= 2 						&& dijetDeltaEta < 2.5;
	pass2JetControl		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>100.   && (dijetDeltaEta < 2.5 || dijet.M()<500. || zep3l > 2.5 ) ;

	
	
	//cout<<zep3l<<endl;

      }

      int theCategory = infilecatv[ifile];
      unsigned int typeFakeLepton[2] = {0,0};
      int nFakeCount = 0;
      double fakeSF = 1.0;
      if(/*usePureMC == false && infilecatv[ifile]<6*/1){
        if     ((infilecatv[ifile] == 0 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add Z+jets from data
	  for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    if(tagZ[0] == (int)nl) nFakeCount = nFakeCount + 1;
	    if(tagZ[1] == (int)nl) nFakeCount = nFakeCount + 2;
	    if(tagZ[2] == (int)nl) nFakeCount = nFakeCount + 4;
	    if(tagZ[2] != (int)nl)
	    fakeSF = fakeSF * fakeRateFactor(leptons.at(nl).Pt(),TMath::Abs(leptons.at(nl).Eta()),TMath::Abs(idLepPdg[nl]),1,"default");
	    else
	    fakeSF = fakeSF * fakeRateFactor(leptons.at(nl).Pt(),TMath::Abs(leptons.at(nl).Eta()),TMath::Abs(idLepPdg[nl]),1,"default");
	    theCategory = 1;
	    if(TMath::Abs((int)idLepPdg[nl]) == 13) 			typeFakeLepton[0]++;
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
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] == 2 && goodIsTight != idTight.size()){ // remove Z+gamma, fakeable objects
          fakeSF = 0.0;
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

      double effSF=1.;
      if(infilecatv[ifile] != 0){
      for(unsigned int l=0;l<leptons.size();l++){
	if(TMath::Abs(idLepPdg[l]) == 11) { 
	  effSF = effSF * histEleRECOSF->GetBinContent(histEleRECOSF->FindBin(leptons[l].Eta(),TMath::Min(leptons[l].Pt(),499.)));
	  effSF = effSF * histEleIDSF->GetBinContent(histEleIDSF->FindBin(leptons[l].Eta(),TMath::Min(leptons[l].Pt(),499.)));		  
	}
        else if(TMath::Abs(idLepPdg[l]) == 13) { 
	  effSF = effSF * histMuIDSF->GetBinContent(histMuIDSF->FindBin(TMath::Abs(leptons[l].Eta()),TMath::Min(leptons[l].Pt(),119.)));
	}
	//if(effSF==0.0) cout<<"eta: "<<leptons[l].Eta()<<"pt: "<<leptons[l].Pt()<<endl;
      }
      //cout<<effSF<<endl;
      }

      

      double totalWeight = 1.;
      if(theCategory==0 && passAllCuts[0]) {dataevents++;}
      if(infilecatv[ifile]==0) {xSecLumiWeight=1.; genWeight=1. ;puWeight=1.; /*fakeSF=1.*/; effSF=1.;}
      totalWeight = xSecLumiWeight*genWeight*puWeight*fakeSF*effSF;
      
      //cout<<i<<"  "<<theCategory<<"   "<<fakeSF<<"   "<<totalWeight<<endl;
      //cout<<i<<"  "<<idLep.size()<<"   "<<goodIsTight<<"   "<<goodIsGenLep<<"   "<<fakeSF<<"   "<<genWeight<<"   "<<puWeight<<"   "<<xSecLumiWeight<<endl;
      //if(idLep.size() != goodIsTight) continue;     
      //cout<<passAllCuts[0]<<endl;

      
      for(int thePlot=0; thePlot<allPlots; thePlot++){
		double theVar = 0.0;
		bool makePlot = false;
 		if     (thePlot ==  0 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(leptons.at(tagZ[0]).Pt(),1199.999);}
		else if(thePlot ==  1 && passAllCuts[0])                      	 {makePlot = true;theVar = TMath::Min(leptons.at(tagZ[1]).Pt(),199.999);}
		else if(thePlot ==  2 && passAllCuts[0])                         {makePlot = true;theVar =  type3l;}
		else if(thePlot ==  3 && fakeEnriched[0])                        {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
		else if(thePlot ==  4 && fakeEnriched[1])                        {makePlot = true;theVar = TMath::Min((double)MET_pt,199.999);}
		else if(thePlot ==  5 && passAllCuts[0])                          {makePlot = true;theVar = TMath::Min((double)MET_pt,199.999);}
		else if(thePlot ==  6 && passAllCuts[0])                          {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
		else if(thePlot ==  7 && passAllCuts[0] && idJet.size()>0.)       {makePlot = true;theVar = TMath::Min((double)jets[0].Pt(),500.);}
		else if(thePlot ==  8 && passAllCuts[0] && idJet.size()>0.)       {makePlot = true;theVar =(double)jets[0].Eta();}
		else if(thePlot ==  9 && passAllCuts[0] )       		{makePlot = true;theVar =(double)leptons.at(tagZ[0]).Eta();}
		else if(thePlot ==  10 && passAllCuts[0] )       		{makePlot = true;theVar =(double)leptons.at(tagZ[1]).Eta();}
		else if(thePlot ==  11 && passAllCuts[0] )       		{makePlot = true;theVar =(double)PV_npvs;}
		
		if( (makePlot && theCategory<7) || (makePlot && theCategory == 7 && infilenamev[ifile].Contains("M900")) ) {histo[     4][thePlot][theCategory]->Fill(theVar,totalWeight);}
		if( (makePlot && theCategory<7) || (makePlot && theCategory == 7 && infilenamev[ifile].Contains("M900")) ) histo[type3l][thePlot][theCategory]->Fill(theVar,totalWeight);
      	  }

	if(1 && (passVBFLoose || passVBFLooseJESDOWN || passVBFLooseJESUP || passVBFLooseJERUP || passVBFLooseJERDOWN || passVBFLooseLEPScaleUP || passVBFLooseLEPScaleDOWN) 
	 && ( theCategory!=0 	/*|| (theCategory==0 	&& !passVBF)*/)) {
	//	if(1 && passVBFLoose && theCategory!=0) {


	double MVAVar 			= -999.;
	double MVAVarLEPScaleUP 	= -999.;
	double MVAVarLEPScaleDOWN 	= -999.;
	double MVAVarJESDOWN 		= -999.;
	double MVAVarJESUP   		= -999.;
	double MVAVarJERUP 		= -999.;
	double MVAVarJERDOWN 		= -999.;
	double bound      	        =  12. ;
	if(isHiggs) bound 		= 1999.9;
      
	if(passVBFLoose){	 
	  MVAVar 		= (double)TMath::Min(superVar(dijet.M(),	dijetDeltaEta,	zp4,	wp4 , zep3l,	isHiggs ),	bound);
      	}          	
	if(passVBFLooseJESDOWN)	 	MVAVarJESDOWN 	= (double)TMath::Min(superVar(dijetjesDown.M()	,dijetDeltaEtaJesDown,	zp4,	wp4JESDOWN	,zep3l, isHiggs)	,bound);
	if(passVBFLooseJESUP)		MVAVarJESUP    	= (double)TMath::Min(superVar(dijetjesUp.M()	,dijetDeltaEtaJesUp,	zp4,	wp4JESUP	,zep3l,	isHiggs)	,bound);
	if(passVBFLooseJERDOWN) 	MVAVarJERDOWN 	= (double)TMath::Min(superVar(dijetjerDown.M()	,dijetDeltaEtaJerDown,	zp4,	wp4JERDOWN	,zep3l,	isHiggs)	,bound);
	if(passVBFLooseJERUP)		MVAVarJERUP    	= (double)TMath::Min(superVar(dijetjerUp.M()	,dijetDeltaEtaJerUp,	zp4,	wp4JESUP	,zep3l,	isHiggs)	,bound);
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

	    /* histo_Zg_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
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
	     for(int npdf=0; npdf<100; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);*/
             histo_Zg_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Zg_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Zg_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Zg_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Zg_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_Zg_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             //histo_Zg_CMS_puUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             //histo_Zg_CMS_puDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
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
	     
	     /*histo_WZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_WZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     /* if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<100; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
             else
	     for(int npdf=0; npdf<100; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);*/
             histo_WZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_WZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_WZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_WZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             //histo_WZ_CMS_puUp  	->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             //histo_WZ_CMS_puDown	->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
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
	     /*histo_ZZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_ZZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     /*if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);*/
             histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             //histo_ZZ_CMS_puUp  	->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             //histo_ZZ_CMS_puDown	->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
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
	     /*histo_VVV_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_VVV_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     /*if(initPDFTag != -1)
	       for(int npdf=0; npdf<100; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<100; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);*/
             histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             //histo_VVV_CMS_puUp  		->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             //histo_VVV_CMS_puDown		->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
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
	     /*histo_EWWZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_EWWZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_EWWZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_EWWZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_EWWZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_EWWZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     /*if(initPDFTag != -1)
	       for(int npdf=0; npdf<100; npdf++) histo_EWWZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
             else if(infilenamev[ifile].Contains("powheg") == true)
	       for(int npdf=0; npdf<100; npdf++) histo_EWWZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
             else
	     for(int npdf=0; npdf<100; npdf++) histo_EWWZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);*/
             histo_EWWZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_EWWZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_EWWZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_EWWZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_EWWZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_EWWZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
	     // histo_EWWZ_CMS_puUp  	->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             //histo_EWWZ_CMS_puDown	->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
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
		   /*histo_Higgs_[i]->Fill(MVAVar,totalWeight);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
		   histo_Higgs_CMS_QCDScaleBoundingPerMass[i][5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
		   /*if(initPDFTag != -1)
		     for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBoundingPerMass[i][npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
		   else if(infilenamev[ifile].Contains("powheg") == true)
		     for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBoundingPerMass[i][npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
		   else
		   for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBoundingPerMass[i][npdf]->Fill(MVAVar,totalWeight);*/
		   histo_Higgs_CMS_MVALepEffMperMassBoundingAvg[i] ->Fill(MVAVar,totalWeight*1.00);
		   histo_Higgs_CMS_MVALepEffEperMassBoundingAvg[i] ->Fill(MVAVar,totalWeight*1.00);
		   histo_Higgs_CMS_MVALepEffMperMassBoundingUp[i]  ->Fill(MVAVar,totalWeight*systTotLep[0]);
		   histo_Higgs_CMS_MVALepEffEperMassBoundingUp[i]  ->Fill(MVAVar,totalWeight*systTotLep[1]);
		   histo_Higgs_CMS_MVALepEffMperMassBoundingDown[i]->Fill(MVAVar,totalWeight/systTotLep[0]);
		   histo_Higgs_CMS_MVALepEffEperMassBoundingDown[i]->Fill(MVAVar,totalWeight/systTotLep[1]);
		   //histo_Higgs_CMS_puUp[i]  			->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
		   //histo_Higgs_CMS_puDown[i]			->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	
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
	  printf("CATEGORY PROBLEM!\n"); return;}
	}//end of datacard writing
    }//end of tree per file
 }//end of all files



for(int thePlot=0; thePlot<allPlots; thePlot++){
    for(int j=0; j<allStates; j++){
      char output[200];
      sprintf(output,"histowz_nice_%d_%d.root",j,thePlot);	  
      TFile* outFilePlotsNote = new TFile(output,"recreate");
      outFilePlotsNote->cd();
      for(int np=0; np<histBins; np++) histo[j][thePlot][np]->Write();
      checkTheID->Write();
      outFilePlotsNote->Close();
    }
  }

  double sumEvents = 0;
  for(int np=1; np<histBins; np++) sumEvents += histo[4][2][np]->GetSumOfWeights();
  printf("yields: %f |",histo[4][2][0]->GetSumOfWeights());
  for(int np=1; np<histBins; np++) printf(" %.3f",histo[4][2][np]->GetSumOfWeights());
  printf(" = %.3f\n",sumEvents);
  double sumEventsType[5] = {0,0,0,0,0};
  double sumEventsTypeE[5] = {0,0,0,0,0};

  sumEventsType[0] = 0.;sumEventsType[1] = 0.;sumEventsType[2] = 0.;sumEventsType[3] = 0.;sumEventsType[4] = 0.;
  sumEventsTypeE[0] = 0.;sumEventsTypeE[1] = 0.;sumEventsTypeE[2] = 0.;sumEventsTypeE[3] = 0.;sumEventsTypeE[4] = 0.;

cout<<endl<<"-------   Inclusive YIELDS    ---------"<<endl;
printf("                  all                 mmm                 eem                 mme                 eee\n");
printf("-----------------------------------------------------------------------------------------------------------\n");
for(int np=0; np<histBins; np++) {
    double sumEvents = 0.; double sumEventsE = 0.;
    for(int i=1; i<=4; i++) {sumEvents = sumEvents + histo[4][2][np]->GetBinContent(i); sumEventsE = sumEventsE + histo[4][2][np]->GetBinError(i)*histo[4][2][np]->GetBinError(i);}
    if(np!=0){
    sumEventsType[0] = sumEventsType[0] + sumEvents;                          sumEventsTypeE[0] = sumEventsTypeE[0] + sumEventsE;
    sumEventsType[1] = sumEventsType[1] + histo[4][2][np]->GetBinContent(1); sumEventsTypeE[1] = sumEventsTypeE[1] + histo[4][2][np]->GetBinError(1) * histo[4][2][np]->GetBinError(1);
    sumEventsType[2] = sumEventsType[2] + histo[4][2][np]->GetBinContent(2); sumEventsTypeE[2] = sumEventsTypeE[2] + histo[4][2][np]->GetBinError(2) * histo[4][2][np]->GetBinError(2);
    sumEventsType[3] = sumEventsType[3] + histo[4][2][np]->GetBinContent(3); sumEventsTypeE[3] = sumEventsTypeE[3] + histo[4][2][np]->GetBinError(3) * histo[4][2][np]->GetBinError(3);
    sumEventsType[4] = sumEventsType[4] + histo[4][2][np]->GetBinContent(4); sumEventsTypeE[4] = sumEventsTypeE[4] + histo[4][2][np]->GetBinError(4) * histo[4][2][np]->GetBinError(4);
    }
    sumEventsE = sqrt(sumEventsE);
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    processName[np].Data(),sumEvents,sumEventsE,
    histo[4][2][np]->GetBinContent(1),histo[4][2][np]->GetBinError(1),histo[4][2][np]->GetBinContent(2),histo[4][2][np]->GetBinError(2),
    histo[4][2][np]->GetBinContent(3),histo[4][2][np]->GetBinError(3),histo[4][2][np]->GetBinContent(4),histo[4][2][np]->GetBinError(4));
    if(np==0)
    printf("-----------------------------------------------------------------------------------------------------------\n");
  }
    printf("-----------------------------------------------------------------------------------------------------------\n");
    printf("(%6s): %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f | %7.2f +/- %5.2f\n",
    "   all",sumEventsType[0],sqrt(sumEventsTypeE[0]),
    sumEventsType[1],sqrt(sumEventsTypeE[1]),sumEventsType[2],sqrt(sumEventsTypeE[2]),
    sumEventsType[3],sqrt(sumEventsTypeE[3]),sumEventsType[4],sqrt(sumEventsTypeE[4]));

    cout<<"dataevents: "<<dataevents<<endl;

  char outputLimits[200]; 
  sprintf(outputLimits,outputDirectory+"sample_%d_wz3l%2s_input_%4s.root",sampleID,finalStateName,ECMsb.Data());
 
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data   ->Write();
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
      newcardShape << Form("shapes * * wz3l3lHig_input_13TeV2017.root histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n");
      newcardShape << Form("shapes data_obs * wz3l3lHig_input_13TeV2017.root histo_Data\n");
      newcardShape << Form("shapes Higgs * wz3l3lHig_input_13TeV2017.root histo_Higgs_M$MASS histo_Higgs_M$MASS_$SYSTEMATIC\n");
      newcardShape << Form("Observation %d\n", -1);
      newcardShape << Form("bin wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
      newcardShape << Form("process Zg VVV WZ ZZ Fake EWWZ Higgs\n");
      newcardShape << Form("process 1 2 5 3 4 6 0 \n");
      newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f %8.5f  \n",-1.,-1.,-1.,-1.,-1.,-1.,-1.) ;
      newcardShape << Form("%s                               lnN  	%7.5f   %7.5f   %7.5f   %7.5f    -    %7.5f    %7.5f   \n","lumi_13TeV"         	,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE);
      newcardShape << Form("%s		                     lnN    	  -      -      %7.5f     -      -     -         -     \n","norm_WZ_13TeV2017"		,1.05);	
      newcardShape << Form("%s				     lnN  	%7.5f   %7.5f     -     %7.5f    -    %7.5f    %7.5f   \n","CMS_eff_b_mistag_13TeV2017" ,0.98,0.98,0.98,0.98,0.98);
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
      newcardShape << Form("shapes * * wz3l3lEW_input_13TeV2017.root histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n");
      newcardShape << Form("shapes data_obs * wz3l3lEW_input_13TeV2017.root histo_Data\n");
      newcardShape << Form("Observation %d\n", -1);
      newcardShape << Form("bin wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
      newcardShape << Form("process Zg VVV WZ ZZ Fake EWWZ\n");
      newcardShape << Form("process 1 2 5 3 4 0 \n");
      newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",-1.,-1.,-1.,-1.,-1.,-1.) ;
      newcardShape << Form("%s                               lnN  	%7.5f   %7.5f   %7.5f   %7.5f    -    %7.5f       \n","lumi_13TeV"         	,lumiE,lumiE,lumiE,lumiE,lumiE);
      newcardShape << Form("%s		                   lnN    	  -      -      %7.5f     -      -     -          \n","norm_WZ_13TeV2017"		,1.05);	
      newcardShape << Form("%s				   lnN  	%7.5f   %7.5f     -     %7.5f    -    %7.5f       \n","CMS_eff_b_mistag_13TeV2017",0.98,0.98,0.98,0.98);
      newcardShape << Form("%s		                   lnN   	  -       -       -       -    %7.5f    -         \n","CMS_wz3l_FakeSys" 		,1.30);	
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n",effMName		        ,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n",effEName			,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n","CMS_pu" 			,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n","CMS_scale_j"	        ,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n","CMS_res_j"		,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f       \n","CMS_scale_l"		,1.,1.,1.,1.,1.);
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
    

}
