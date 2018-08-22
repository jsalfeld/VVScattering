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

#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"


#include "MitAnalysisRunII/macros/80x/factors.h"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"


double superVar(double mjjIn, double jetDetaJJ, TLorentzVector Zp4, TLorentzVector Wp4, double Zep3l, bool isHiggsX=false ){

double bin = -9.;

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


double superVarSep(double mjjIn, double jetDetaJJ, TLorentzVector Zp4, TLorentzVector Wp4, double Zep3l, int channel,bool isHiggsX=false ){

double bin = -9.;

 if(!isHiggsX){
   if((jetDetaJJ <= 2.5 || mjjIn<=500 || Zep3l > 2.5))      bin=0. + channel*13.;

   if(jetDetaJJ <= 3.5 	&& jetDetaJJ > 2.5 && Zep3l < 2.5){
     if(mjjIn>500. 	&& mjjIn<=1000.)       	bin=1.+ channel*13.;
     if(mjjIn>1000. 	&& mjjIn<=1350.) 	bin=2.+ channel*13.;
     if(mjjIn>1350. 	&& mjjIn<=1750.) 	bin=3.+ channel*13.;
     if(mjjIn>1750.)             		bin=4.+ channel*13.;
   }

   if(jetDetaJJ <= 5. 	&& jetDetaJJ > 3.5 && Zep3l < 2.5){
     if(mjjIn>500. 	&& mjjIn<=1000.) 	bin=5.+ channel*13.;
     if(mjjIn>1000. 	&& mjjIn<=1350.) 	bin=6.+ channel*13.;
     if(mjjIn>1350. 	&& mjjIn<=1750.) 	bin=7.+ channel*13.;
     if(mjjIn>1750.)             		bin=8.+ channel*13.;
   }

   if( jetDetaJJ > 5. && Zep3l < 2.5){
     if(mjjIn>500. 	&& mjjIn<=1000.) 	bin=9.+ channel*13.;
     if(mjjIn>1000. 	&& mjjIn<=1350.) 	bin=10.+ channel*13.;
     if(mjjIn>1350. 	&& mjjIn<=1750.) 	bin=11.+ channel*13.;
     if(mjjIn>1750.)             		bin=12.+ channel*13.;
  }
 }

if(isHiggsX){

   if( jetDetaJJ <= 2.5     || mjjIn<=500 )      bin=-1;
   if( jetDetaJJ  > 2.5     && mjjIn >500 ){
    	bin = TMath::Sqrt( TMath::Power(Zp4.Et() + Wp4.Et(),2) - TMath::Power(Zp4.Px() + Wp4.Px(),2) - TMath::Power(Zp4.Py() + Wp4.Py(),2));

   }
}

if( bin == -9 ) cout<<"superVar DANGER: "<<"mJJ: "<<mjjIn<<",deltaJJ: "<<jetDetaJJ<<",zep3l:  "<<Zep3l<<endl;

return bin;
}




enum selType                     { SIGSEL, nSelTypes};
TString selTypeName[nSelTypes]= { "SIGSEL"};

enum systType                     {METUP=0, METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"METUP","METDOWN"};

double mcPrescale = 1.0;
bool usePureMC = false;

int whichSkim = 2;

vector<long int> eventNumbers0;
vector<long int> eventNumbers1;
vector<long int> eventNumbers2;
vector<long int> eventNumbers3;
vector<long int> lsNumbers0;
vector<long int> lsNumbers1;
vector<long int> lsNumbers2;
vector<long int> lsNumbers3;

vector<long int> runNumbers0;
vector<long int> runNumbers1;
vector<long int> runNumbers2;
vector<long int> runNumbers3;

vector<vector<long int>> eventNumbers(4);
vector<vector<long int>> lsNumbers(4);
vector<vector<long int>> runNumbers(4);

TRandom2 r;

TString triggerSuffix = "";

void wzScatteringAnalysisBatchSyncImprovedChannel(
 double minMass =  76,
 double maxMass = 106,
 bool applyBtagging = true,
 TString typeLepSel = "default",
 TString type3rdLepSel = "default",
 TString outputDirectory = "./",
 bool isHiggs = false,
 int sampleID = 7
 ){

  double jetIDpTcut = 50;

  TString SignalSuffix[11] = {"M200","M300","M400","M500","M600","M700","M800","M900","M1000","M1500","M2000"};
  
  Int_t period = 1;
  TString filesPathDA  = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/met_";
  TString filesPathMC  = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/met_";
  TString filesPathMC2 = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/mc/met_";

  Double_t lumi = 35.9;

  cout<<"APPETIZER"<<endl;

  Float_t new_v;
  Float_t jetOnePt;
  Float_t jetTwoPt;
  Float_t jetTwoEta;
  Float_t jetOneEta;
  Float_t jetTwoPhi;
  Float_t jetOnePhi;
  Float_t jetDeltaEta;
  Float_t jetDeltaPhi;
  Float_t jetMJJ;
  Float_t jetZepp;
  Float_t trainWeight;
  
  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString>	infilenamev;  
  vector<Int_t> 	infilecatv;  

  TString puPath = "";
  if      (period==1){

    puPath = "MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root";


    infilenamev.push_back(Form("%sdata_Run2016B.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016C.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016D.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016E.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016F.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016G.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016H.root",filesPathDA.Data())); infilecatv.push_back(0);
  

  if(usePureMC == true){
infilenamev.push_back(Form("%sTT_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(1);
  }



  infilenamev.push_back("/eos/cms/store/caf/user/ceballos/Nero/output_80x/met_WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8.root");   	infilecatv.push_back(6);
  
  for(int i=0; i<11;i++){
    if(i<1) {infilenamev.push_back("/eos/cms/store/user/jsalfeld/Nero/v0.2/met_ChargedHiggsToWZTo3LNu_M200_13TeV-madgraph-pythia8.root")	 ;infilecatv.push_back(7);}
    if(i>=1){   
	infilenamev.push_back("/eos/cms/store/caf/user/ceballos/Nero/output_jsalfeld_80x/met_ChargedHiggsToWZTo3LNu_"+SignalSuffix[i]+"_13TeV-madgraph-pythia8.root")	;infilecatv.push_back(7);
 	cout<<SignalSuffix[i]<<endl;}
     }	

  //infilenamev.push_back("/eos/cms/store/caf/user/ceballos/Nero/output_jsalfeld_80x/met_ChargedHiggsToWZTo3LNu_M400_13TeV-madgraph-pythia8.root");   	infilecatv.push_back(7);
 

  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));                   	infilecatv.push_back(2);
   
 

  //infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root",filesPathMC.Data()));			      	infilecatv.push_back(3);
   
      infilenamev.push_back(Form("%sWZTo3LNu_0Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));			      	infilecatv.push_back(3);
   
   	infilenamev.push_back(Form("%sWZTo3LNu_1Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));			      	infilecatv.push_back(3);
   
   	infilenamev.push_back(Form("%sWZTo3LNu_2Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));			      	infilecatv.push_back(3);
   
   	infilenamev.push_back(Form("%sWZTo3LNu_3Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));			      	infilecatv.push_back(3);
 

  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8.root",filesPathMC.Data()));					      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZJJTo4L_EWK_13TeV-madgraph-pythia8.root",filesPathMC.Data()));					      infilecatv.push_back(4);
  

  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));	 	      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(4);
  
  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data())); 			      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data())); 			      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                           infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                           infilecatv.push_back(5);  
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));      infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));           infilecatv.push_back(5);
   infilenamev.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                                  infilecatv.push_back(5);  
  //infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root",filesPathMC.Data()));			      	infilecatv.push_back(3);
  
  }
  else {assert(0);}
  //std::cout<<"hallo"<<std::endl;
  if(infilenamev.size() != infilecatv.size()) {assert(0); cout<<"cannot be"<<endl;return;}

  /*infilenamev.clear();infilecatv.clear();

  for(int i=0; i<11;i++){
    if(i<1) infilenamev.push_back("/eos/cms/store/user/jsalfeld/Nero/v0.2/met_ChargedHiggsToWZTo3LNu_M200_13TeV-madgraph-pythia8.root")					 ;infilecatv.push_back(7);
    if(i>=1) { 
	infilenamev.push_back("/eos/cms/store/caf/user/ceballos/Nero/output_jsalfeld_80x/met_ChargedHiggsToWZTo3LNu_"+SignalSuffix[i]+"_13TeV-madgraph-pythia8.root");   infilecatv.push_back(7);
 	cout<<SignalSuffix[i]<<endl;}
	}	*/
  //infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8.root",filesPathMC.Data()));					      infilecatv.push_back(4);
//infilenamev.push_back("/eos/cms/store/user/jsalfeld/WllMoriond17/met_WLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8.root");  infilecatv.push_back(6);


  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
  //cout<<"hi"<<endl;

  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));


  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
  delete fPUFile;
  //cout<<"hello1"<<endl;  

  TFile *fTrackElectronReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root"));
  TH2D *fhDeltrksf= (TH2D*)(fTrackElectronReco_SF->Get("scalefactors_Reco_Electron")); assert(fhDeltrksf); fhDeltrksf->SetDirectory(0);
   delete fTrackElectronReco_SF;

   TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root"));
   TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("scalefactors_Medium_Electron"));
   TH2D *fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_Tight_Electron"));
   assert(fhDElMediumSF);
   assert(fhDElTightSF);
   fhDElMediumSF->SetDirectory(0);
   fhDElTightSF->SetDirectory(0);
   delete fElSF;



   //TFile *fTrackMuonReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/trackMuReco_SF.root"));
   TFile *fTrackMuonReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/Tracking_EfficienciesAndSF_BCDEFGH.root"));
   //TH1D *fhDmutrksfptg10 = (TH1D*)(fTrackMuonReco_SF->Get("mutrksfptg10")); assert(fhDmutrksfptg10); fhDmutrksfptg10->SetDirectory(0);
    TH1D *fhDmutrksfptg10 = (TH1D*)(fTrackMuonReco_SF->Get("ratio_eff_eta3_dr030e030_corr")); assert(fhDmutrksfptg10); fhDmutrksfptg10->SetDirectory(0);
   //TH1D *fhDmutrksfptl10 = (TH1D*)(fTrackMuonReco_SF->Get("mutrksfptl10")); assert(fhDmutrksfptl10); fhDmutrksfptl10->SetDirectory(0);
   delete fTrackMuonReco_SF;

 
   TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root"));
   TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("scalefactors_TightId_Muon")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
   delete fMuSF;



   TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root"));
   TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("scalefactors_Iso_MuonTightId")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
   delete fMuIsoSF;



 TString theVeryTightSFName = "MitAnalysisRunII/data/80x/veryTightSF_37ifb.root";

  if(strcmp(typeLepSel.Data(),"veryverytight")==0){
    theVeryTightSFName = "MitAnalysisRunII/data/80x/veryveryTightSF_37ifb.root";
    printf("Using veryverytight SF\n");
  }

  TFile *fElVeryTightSF = TFile::Open(Form("%s",theVeryTightSFName.Data()));
  TH1D *fhDVeryTightSF = (TH1D*)(fElVeryTightSF->Get("veryTightSF"));
  assert(fhDVeryTightSF);
  fhDVeryTightSF->SetDirectory(0);
  delete fElVeryTightSF;


cout<<"Hi"<<endl;


  double totalFakeDataCount[4][9];
  for(int i=0; i<4; i++) for(int j=0; j<9; j++) totalFakeDataCount[i][j] = 0;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 65;
  const int histBins = 8;
  const int allStates = 5;
  TH1D* histo[allStates][allPlots][histBins];
  TString processName[histBins] = {"..Data", ".Fakes", "Zgamma", "....WZ", "....ZZ", "...VVV","..EWWZ" ,".Higgs"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 50; 	xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >=  1 && thePlot <=  2) {nBinPlot = 17; 	xminPlot = 30.0; xmaxPlot = 200.0;}
    else if(thePlot >=  3 && thePlot <=  3) {nBinPlot = 100; 	xminPlot =50.0; xmaxPlot = 150.0;}
    else if(thePlot >=  4 && thePlot <=  4) {nBinPlot = 200; 	xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  5 && thePlot <=  5) {nBinPlot = 200; 	xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >=  6 && thePlot <=  6) {nBinPlot = 100; 	xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =  80; 	xminPlot = 0.0; xmaxPlot =   4.0;}
    else if(thePlot >=  8 && thePlot <= 11) {nBinPlot = 100; 	xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 12 && thePlot <= 12) {nBinPlot =   7; 	xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 13 && thePlot <= 13) {nBinPlot =  40; 	xminPlot = 0.0; xmaxPlot =  40.0;}
    else if(thePlot >= 14 && thePlot <= 14) {nBinPlot =  40; 	xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot >= 15 && thePlot <= 15) {nBinPlot =   4; 	xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >= 16 && thePlot <= 17) {nBinPlot =  90; 	xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 18 && thePlot <= 18) {nBinPlot =   7; 	xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 100; 	xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 20 && thePlot <= 20) {nBinPlot = 100; 	xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot >= 21 && thePlot <= 21) {nBinPlot = 100; 	xminPlot = -1.0; xmaxPlot = 10.0;}
    else if(thePlot >= 22 && thePlot <= 22) {nBinPlot = 100; 	xminPlot = 0.0; xmaxPlot = 1000.0;}
    else if(thePlot >= 23 && thePlot <= 24) {nBinPlot = 100; 	xminPlot = 0.0; xmaxPlot = 1000.0;}
    else if(thePlot >= 25 && thePlot <= 25) {nBinPlot =  50; 	xminPlot =-0.5; xmaxPlot =  49.5;}
    else if(thePlot >= 26 && thePlot <= 31) {nBinPlot = 200; 	xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 32 && thePlot <= 32) {nBinPlot = 100; 	xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot >= 33 && thePlot <= 36) {nBinPlot =   7; 	xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 37 && thePlot <= 37) {nBinPlot =  50; 	xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >= 38 && thePlot <= 38) {nBinPlot = 100; 	xminPlot = 0.0; xmaxPlot =2000.0;}
    else if(thePlot >= 39 && thePlot <= 39) {nBinPlot =  80; 	xminPlot = 0.0; xmaxPlot =   8.0;}
    else if(thePlot >= 40 && thePlot <= 40) {nBinPlot =   4; 	xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >= 41 && thePlot <= 42) {nBinPlot =  50;   xminPlot = 0.0; xmaxPlot =   500;}
    else if(thePlot >= 43 && thePlot <= 44) {nBinPlot =  10;    xminPlot =-5.0; xmaxPlot =   5.0;}
    else if(thePlot >= 45 && thePlot <= 45) {nBinPlot =  40;     xminPlot = 100.0; xmaxPlot =  1000;}
    else if(thePlot >= 46 && thePlot <= 46) {nBinPlot =  10;    xminPlot = 0.0; xmaxPlot =   5.0;}
    else if(thePlot >= 47 && thePlot <= 49) {nBinPlot =   20;    xminPlot = 0.0; xmaxPlot =   200;}
    else if(thePlot >= 50 && thePlot <= 52) {nBinPlot =   50;    xminPlot =-2.5; xmaxPlot =   2.5;}
    else if(thePlot >= 53 && thePlot <= 53) {nBinPlot =   4;     xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >= 54 && thePlot <= 54) {nBinPlot =   200;   xminPlot =0.0; xmaxPlot =   1000;}
    else if(thePlot >= 55 && thePlot <= 55) {nBinPlot =   20;    xminPlot =0.0; xmaxPlot =   200;}
    else if(thePlot >= 56 && thePlot <= 56) {nBinPlot =   30;    xminPlot =0.0; xmaxPlot =   1500;}
    else if(thePlot >= 57 && thePlot <= 57) {nBinPlot =   50;    xminPlot =0.0; xmaxPlot =   10.0;}
    else if(thePlot >= 58 && thePlot <= 58) {nBinPlot =   50;    xminPlot =0.0; xmaxPlot =   1000.0;}
    else if(thePlot >= 59 && thePlot <= 59) {nBinPlot =   15;    xminPlot =0.0; xmaxPlot =   15.0;}
    else if(thePlot >= 60 && thePlot <= 60) {nBinPlot =   4;   xminPlot =-0.5; xmaxPlot =   3.5;}//Jet1 pT 2jet
    else if(thePlot >= 61 && thePlot <= 61) {nBinPlot =   4;    xminPlot =-0.5; xmaxPlot =   3.5;}//jet2 pT 2Jet
    else if(thePlot >= 62 && thePlot <= 62) {nBinPlot =   50;    xminPlot =0.0; xmaxPlot =   500.0;}//l1 pT 2jet
    else if(thePlot >= 63 && thePlot <= 63) {nBinPlot =   50;    xminPlot =0.0; xmaxPlot =   500.0;}//l2 pT 2jet
    else if(thePlot >= 64 && thePlot <= 64) {nBinPlot =   50;    xminPlot =0.0; xmaxPlot =   500.0;}//l3 pT 2jet
    else if(thePlot >= 65 && thePlot <= 65) {nBinPlot =   100;   xminPlot =0.0; xmaxPlot =   1000.0;}//MJJ
    else if(thePlot >= 66 && thePlot <= 66) {nBinPlot =   10;    xminPlot =0.0; xmaxPlot =   10.0;}//deta
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) {
      for(int j=0; j<allStates; j++) histo[j][thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    }
    histos->Reset();histos->Clear();
  }

    TH2D* Higgs_2d = new TH2D("Higgs_2d", "Higgs_2d", 100, 0., 1000, 100, 0,1000);
    Higgs_2d->Sumw2();
    TH2D* WZ_2d = new TH2D("WZ_2d", "WZ_2d", 100, 0., 1000, 100,0,1000);
    Higgs_2d->Sumw2();
    cout<<"Hi"<<endl;

  TString ECMsb  = "13TeV2016";
  //const int nBinMVA = 4; Float_t xbins[nBinMVA+1] = {0, 100, 250, 450, 1500};
  //const int nBinMVA = 5; Float_t xbins[nBinMVA+1] = {250., 500., 1000., 1350., 1750. , 2000. };
  const int nBinMVA = 52; Float_t xbins[nBinMVA+1] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.,51.,52.}; 
  Float_t xbinsX[nBinMVA+1] = {-1.,0.,50.,100.,150.,200.,250.,300.,400.,500.,700.,1000.,1500.,2000.};
  //Float_t xbins[nBinMVA+1] = {-1.,0.,50.,100.,150.,200.,250.,300.,400.,500.,700.,1000.,1250.,1500.,1750.,2000.};
  //Float_t xbinsX[nBinMVA+1] = {-1.,0,100.,200.,300.,400.,500.,700.,800.,1000.,1250.,1500.,1750.,2000.};
  if( isHiggs)  {for(int k=0; k<nBinMVA+1; k++) xbins[k] = xbinsX[k];}
        //for(int k=0; k<nBinMVA+1; k++) cout<<xbins[k]<<endl;
  // const int nBinMVA = 20; Float_t xbins[nBinMVA+1] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900,950,1000};
  //const int nBinMVA = 1; Float_t xbins[nBinMVA+1] = {0, 1};
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
cout<<"Hi"<<endl;

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


  TH1D* histo_Fake_CMS_wz3l_FakeSys0Up            = 	new TH1D( "histo_Fake_CMS_wz3l_FakeSys0Up"  	, "histo_Fake_CMS_wz3l_FakeSys0Up"  		, nBinMVA, xbins);	histo_Fake_CMS_wz3l_FakeSys0Up  ->Sumw2();
  TH1D* histo_Fake_CMS_wz3l_FakeSys0Down          = 	new TH1D( "histo_Fake_CMS_wz3l_FakeSys0Down"  	, "histo_Fake_CMS_wz3l_FakeSys0Down"  		, nBinMVA, xbins); 	histo_Fake_CMS_wz3l_FakeSys0Down  ->Sumw2();
  TH1D* histo_Fake_CMS_wz3l_FakeSys1Up            = 	new TH1D( "histo_Fake_CMS_wz3l_FakeSys1Up"  	, "histo_Fake_CMS_wz3l_FakeSys1Up"  		, nBinMVA, xbins);	histo_Fake_CMS_wz3l_FakeSys1Up  ->Sumw2();
  TH1D* histo_Fake_CMS_wz3l_FakeSys1Down          = 	new TH1D( "histo_Fake_CMS_wz3l_FakeSys1Down"  	, "histo_Fake_CMS_wz3l_FakeSys1Down"  		, nBinMVA, xbins); 	histo_Fake_CMS_wz3l_FakeSys1Down  ->Sumw2();
  TH1D* histo_Fake_CMS_wz3l_FakeSys2Up            = 	new TH1D( "histo_Fake_CMS_wz3l_FakeSys2Up"  	, "histo_Fake_CMS_wz3l_FakeSys2Up"  		, nBinMVA, xbins);	histo_Fake_CMS_wz3l_FakeSys2Up  ->Sumw2();
  TH1D* histo_Fake_CMS_wz3l_FakeSys2Down          = 	new TH1D( "histo_Fake_CMS_wz3l_FakeSys2Down"  	, "histo_Fake_CMS_wz3l_FakeSys2Down"  		, nBinMVA, xbins); 	histo_Fake_CMS_wz3l_FakeSys2Down  ->Sumw2();
  TH1D* histo_Fake_CMS_wz3l_FakeSys3Up            = 	new TH1D( "histo_Fake_CMS_wz3l_FakeSys3Up"  	, "histo_Fake_CMS_wz3l_FakeSys3Up"  		, nBinMVA, xbins);	histo_Fake_CMS_wz3l_FakeSys3Up  ->Sumw2();
  TH1D* histo_Fake_CMS_wz3l_FakeSys3Down          = 	new TH1D( "histo_Fake_CMS_wz3l_FakeSys3Down"  	, "histo_Fake_CMS_wz3l_FakeSys3Down"  		, nBinMVA, xbins); 	histo_Fake_CMS_wz3l_FakeSys3Down  ->Sumw2();
  
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

TFile* fakeFile = TFile::Open("MitAnalysisRunII/data/80x/fakeShapeInclusive.root");  
 TH1F* fakeShape = (TH1F*)(fakeFile->Get("histo_Fake"));assert(fakeShape); fakeShape->SetDirectory(0);
  delete fakeFile;
  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    cout<<"number of samples: "	<<infilenamev.size()<<endl;
    cout<<"sample ID: "		<<ifile<<"  "<<sampleID<<endl;
    cout<<"sample name: "	<<infilenamev[ifile]<<endl;
    if(ifile != sampleID && sampleID > -1) continue;
  
    
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());

    TFile *the_input_file=TFile::Open(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");
    TTree *the_SelBit_tree= (TTree*)the_input_file->FindObjectAny("SelBit_tree");
    TTree *the_PDF_tree   = (TTree*)the_input_file->FindObjectAny("pdfReweight");

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareJets eventJets;
    eventJets.setBranchAddresses(the_input_tree);

    BareLeptons eventLeptons;
    eventLeptons.setBranchAddresses(the_input_tree);

    BareTaus eventTaus;
    eventTaus.setBranchAddresses(the_input_tree);

    BareMet eventMet;
    eventMet.SetExtend();
    eventMet.setBranchAddresses(the_input_tree);

    BareTrigger eventTrigger;
    eventTrigger.setBranchAddresses(the_input_tree);

    BareVertex eventVertex;
    eventVertex.setBranchAddresses(the_input_tree);

    BareMonteCarlo eventMonteCarlo;
    eventMonteCarlo.SetExtend();
    eventMonteCarlo.setBranchAddresses(the_input_tree);
    //cout<<"here"<<endl;

    TNamed *triggerNames = (TNamed*)the_input_file->FindObjectAny("triggerNames");
    char **tokens;
    size_t numtokens;
    tokens = strsplit(triggerNames->GetTitle(), ",", &numtokens);
    if(infilecatv[ifile] == 0){
      for (int i = 0; i < (int)numtokens; i++) {
        printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
      }
    }
    else {
    }

    unsigned int selBit_= 0;
    
    the_SelBit_tree->SetBranchAddress("selBit", &selBit_);
    
    double theMCPrescale = mcPrescale;
    int numberOfEvents = int(the_input_tree->GetEntries()/theMCPrescale);
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    
    for (int i=0; i < numberOfEvents; ++i) {
    
      the_SelBit_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      
      cutFlow->Fill(0.);    
      if((selBit_ & 0x1<<whichSkim) == 0 ) continue;
      cutFlow->Fill(1.);    
      
      the_input_tree->GetEntry(i);

      //cout<<(eventEvent.selBits)<<endl;
      if(!(eventEvent.filterbadPFMuon)){continue;}
      if(!((eventEvent.selBits & BareEvent::FullRecommendation) == BareEvent::FullRecommendation)){continue;}

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
	      (strcmp(tokens[nt],Form("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0)    ||
	      (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0) ||
	      (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0)    ||
	      (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0)  ||
	      (strcmp(tokens[nt],Form("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0)     ||
	      (strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v%s",triggerSuffix.Data()))  == 0)  ||
	      (strcmp(tokens[nt],Form("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v%s",triggerSuffix.Data()))  == 0)
		)   passFilter[1]=kTRUE;	      
	   }

      if(passFilter[0] == kFALSE) continue;
      cutFlow->Fill(2.);    
	
      if(passFilter[1] == kFALSE) continue;
      cutFlow->Fill(3.);    
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {

	if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt())< 10.) continue;

        if(selectIdIsoCut(typeLepSel.Data(),
		TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   	TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),
		(double)(*eventLeptons.iso)[nlep],
		(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))                                                                    
			{idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}

        else if( (((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) ){idTight.push_back(0); idLep.push_back(nlep); if(numberOfLeptons == 4) goodIsTight++;}
        else if( (((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP)){idSoft.push_back(nlep);}
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


      bool passLepPtCutScaleUp	=false;
      bool passLepPtCutScaleDown=false;
      bool passLepPtCut		=false;

      if(goodIsTight == idTight.size()) passFilter[3] = kTRUE;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() >= 25 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() >= 20 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() >= 15) passLepPtCut=true;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt()*0.99 >= 25 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt()*0.99 >= 20 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt()*0.99 >= 15) passLepPtCutScaleDown=true;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt()*1.01 >= 25 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt()*1.01 >= 20 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt()*1.01 >= 15) passLepPtCutScaleUp=true;

      if(!(passLepPtCut || passLepPtCutScaleUp || passLepPtCutScaleDown)) continue;
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

      bool passNMinusOne[6] = {                 passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCut,
                               passFilter[5] &&                  passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCut,
                               passFilter[5] && passFilter[6] &&                  passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCut,
                               passFilter[5] && passFilter[6] && passFilter[7] &&                  passFilter[9] && passFilter[10] && passLepPtCut,
                               passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8]                  && passFilter[10] && passLepPtCut,
			       passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] 		   && passLepPtCut};

      bool passAllCuts[1] 		= {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCut};
      bool passAllCutsLepScaleUp[1] 	= {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCutScaleUp};
      bool passAllCutsLepScaleDown[1] 	= {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCutScaleDown};
     

      bool controlSel[4] = {passFilter[5] && passFilter[6] && !passFilter[7] && passFilter[9] && bDiscrMax > 0.605 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20 && passLepPtCut,
                            passFilter[5] &&                  !passFilter[7] && passFilter[9] && bDiscrMax > 0.605 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20 && passLepPtCut,
			    passFilter[5] && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() > 20 		 && passLepPtCut,
			    passFilter[5] && !passFilter[6] && passFilter[7] && passFilter[9] &&           ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20		 && passLepPtCut};

      bool passTTZSel[2] = {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && idJet.size() >= 4 && btagjet[0] > 0.560 && btagjet[1] > 0.560 && passLepPtCut,
                            passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && idJet.size() == 3 && btagjet[0] > 0.800 && btagjet[1] > 0.560 && passLepPtCut};


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
     
     bool passAllCutsJESUP 	= {passFilter[5] &&  metJESUP   > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCut};
     bool passAllCutsJESDOWN    = {passFilter[5] &&  metJESDOWN > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCut};
     bool passAllCutsJERUP 	= {passFilter[5] &&  metJERUP   > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCut};
     bool passAllCutsJERDOWN    = {passFilter[5] &&  metJERDOWN > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCut};

     bool passAllCutsZEnriched  = {passFilter[5] &&     met_t < 30   && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10] && passLepPtCut};

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
     double zep3lJesUp = 999;
     double zep3lJesDown = 999;
     double zep3lJerUp = 999;
     double zep3lJerDown = 999;
     if(idJet.size()>=2){ 
	dijet 		= 	(( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) )   + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) ));
	dijetDeltaEta 	= 	TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta() - ((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta());
	zep3l 		= 	TMath::Abs(eta3l - 0.5*(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta() + ((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()));

	passVBFLoose 		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>100.;     
       	passVBFLooseLEPScaleUP   = passAllCutsLepScaleUp[0] 		&& idJet.size() >= 2 	&& dijet.M() 		>100.;
      	passVBFLooseLEPScaleDOWN = passAllCutsLepScaleDown[0] 		&& idJet.size() >= 2 	&& dijet.M() 		>100.;

      	passVBF 		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>500.	&& dijetDeltaEta > 2.5;     
        passVBFzep 		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>500.	&& dijetDeltaEta > 2.5	&& zep3l<2.5;     
	passVBFzep2Jet50 	= passAllCuts[0] 		&& idJet50.size() >= 2 		&&  	(( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M()	>500.	       	&& TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()- ((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()) > 2.5				&& zep3l<2.5; 
    
      	passVBFLEPScaleUP  	= passAllCutsLepScaleUp[0] 	&& idJet.size() >= 2 		&& dijet.M() 		>500.	&& dijetDeltaEta > 2.5;
      	passVBFLEPScaleDOWN 	= passAllCutsLepScaleDown[0] 	&& idJet.size() >= 2 		&& dijet.M() 		>500.	&& dijetDeltaEta > 2.5;
	
	passVBFLooseMJJLess500  = passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		<500.;
	passVBFLooseDEtaLess25  = passAllCuts[0] 		&& idJet.size() >= 2 						&&  dijetDeltaEta < 2.5;
	pass2JetControl		= passAllCuts[0] 		&& idJet.size() >= 2 		&& dijet.M() 		>100.   && (dijetDeltaEta < 2.5 || dijet.M()<500. || zep3l > 2.5 ) ;	
      }

      if(idJetJESUP.size()>=2){ 
	zep3lJesUp 		= TMath::Abs(eta3l - 0.5*(((TLorentzVector*)(*eventJets.p4)[idJetJESUP[0]])->Eta() + ((TLorentzVector*)(*eventJets.p4)[idJetJESUP[1]])->Eta()));
       	passVBFLooseJESUP 	= passAllCutsJESUP 		&& idJetJESUP.size()   >= 2 	&& dijetjesUp.M()	>100.;      	
      	passVBFJESUP 		= passAllCutsJESUP 		&& idJetJESUP.size()   >= 2 	&& dijetjesUp.M()	>500.	&& deltaEtaJJjesUp > 2.5;
      	      }
     if(idJetJESDOWN.size()>=2){ 
	zep3lJesDown 		= TMath::Abs(eta3l - 0.5*(((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[0]])->Eta() + ((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[1]])->Eta()));
       	passVBFLooseJESDOWN 	= passAllCutsJESDOWN 		&& idJetJESDOWN.size()   >= 2 	&& dijetjesDown.M()	>100.;      	
      	passVBFJESDOWN 		= passAllCutsJESDOWN 		&& idJetJESDOWN.size()   >= 2 	&& dijetjesDown.M()	>500.	&& deltaEtaJJjesDown > 2.5;
      	      }
      if(idJetJERUP.size()>=2){ 
	zep3lJerUp 		= TMath::Abs(eta3l - 0.5*(((TLorentzVector*)(*eventJets.p4)[idJetJERUP[0]])->Eta() + ((TLorentzVector*)(*eventJets.p4)[idJetJERUP[1]])->Eta()));
       	passVBFLooseJERUP 	= passAllCutsJERUP 		&& idJetJERUP.size()   >= 2 	&& dijetjerUp.M()	>100.;      	
      	passVBFJERUP 		= passAllCutsJERUP 		&& idJetJERUP.size()   >= 2 	&& dijetjerUp.M()	>500.	&& deltaEtaJJjerUp > 2.5;
      	      }
     if(idJetJERDOWN.size()>=2){
	zep3lJerDown 		= TMath::Abs(eta3l - 0.5*(((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[0]])->Eta() + ((TLorentzVector*)(*eventJets.p4)[idJetJERDOWN[1]])->Eta())); 
       	passVBFLooseJERDOWN 	= passAllCutsJERDOWN 		&& idJetJERDOWN.size()   >= 2 	&& dijetjerDown.M()	>100.;      	
      	passVBFJERDOWN 		= passAllCutsJERDOWN 		&& idJetJERDOWN.size()   >= 2 	&& dijetjerDown.M()	>500.	&& deltaEtaJJjerDown > 2.5;
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

      //cout<<effSF<<"   "<<type3l<<endl;
      double totalWeight 			= trigEff*mcWeight*theLumi*effSF*puWeight*fakeSF*theMCPrescale;   
      if(infilecatv[ifile] == 7) totalWeight	= trigEff*mcWeight*theLumi*effSF*puWeight*fakeSF*theMCPrescale/1000.;
      if(infilenamev[ifile].Contains("tZq")) totalWeight = totalWeight*1.25;
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
        else if(thePlot == 53 && passVBF )               				{makePlot = true;theVar = (double)  type3l;}
	else if(thePlot == 54 && passVBF )                                  		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(tagZ[1])) )->Pt();}
	else if(thePlot == 55 && pass1Jet)                                 		{makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[0])) )->Pt();}
	else if(thePlot == 56 && passVBF )                          			{makePlot = true;theVar = (double)transverseMass2;}
	else if(thePlot == 57 && passVBFLooseMJJLess500 )                            	{makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta());}
	else if(thePlot == 58 && passVBFLooseDEtaLess25 )                            	{makePlot = true;theVar = (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M();}
        else if(thePlot == 59 && passVBFLooseMJJLess500 )                            	{makePlot = true;theVar = (double)TMath::Min(superVar(dijet.M(),	dijetDeltaEta,	zp4,	wp4, zep3l,	isHiggs),13.);}
	else if(thePlot == 60 && passVBFLoose && !passVBFzep)                            	{makePlot = true;theVar = (double)type3l;}
	else if(thePlot == 61 && passVBFzep2Jet50 )                            	{makePlot = true;theVar = (double)type3l;}
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
	 && (1/* theCategory!=0 	/*|| (theCategory==0 	&& !passVBF)*/)) {
	//	if(1 && passVBFLoose && theCategory!=0) {


	double MVAVar 			= -999.;
	double MVAVarLEPScaleUP 	= -999.;
	double MVAVarLEPScaleDOWN 	= -999.;
	double MVAVarJESDOWN 		= -999.;
	double MVAVarJESUP   		= -999.;
	double MVAVarJERUP 		= -999.;
	double MVAVarJERDOWN 		= -999.;
	double bound      	        = 53.;
	if(isHiggs) bound 		= 1999.9;
      
	if(passVBFLoose){ 		MVAVar 			= (double)TMath::Min(superVarSep(dijet.M()		, dijetDeltaEta,	zp4,			wp4 		, zep3l, type3l,isHiggs ),		 bound); }          	
	if(passVBFLooseJESDOWN)	 	MVAVarJESDOWN 		= (double)TMath::Min(superVarSep(dijetjesDown.M()	,deltaEtaJJjesDown,	zp4,			wp4JESDOWN	, zep3lJesDown,	type3l, isHiggs)	,bound);
	if(passVBFLooseJESUP)		MVAVarJESUP    		= (double)TMath::Min(superVarSep(dijetjesUp.M()		,deltaEtaJJjesUp,	zp4,			wp4JESUP	, zep3lJesUp,	type3l,	isHiggs)	,bound);
	if(passVBFLooseJERDOWN) 	MVAVarJERDOWN 		= (double)TMath::Min(superVarSep(dijetjerDown.M()	,deltaEtaJJjerDown,	zp4,			wp4JERDOWN	, zep3lJerDown,	type3l,	isHiggs)	,bound);
	if(passVBFLooseJERUP)		MVAVarJERUP    		= (double)TMath::Min(superVarSep(dijetjerUp.M()		,deltaEtaJJjerUp,	zp4,			wp4JESUP	, zep3lJerUp,	type3l,	isHiggs)	,bound);
	if(passVBFLooseLEPScaleDOWN) 	MVAVarLEPScaleDOWN 	= (double)TMath::Min(superVarSep(dijet.M()		,dijetDeltaEta,		zp4LEPScaleDOWN,	wp4LEPScaleDOWN	, zep3l,	type3l,	isHiggs)	,bound);
	if(passVBFLooseLEPScaleUP)	MVAVarLEPScaleUP    	= (double)TMath::Min(superVarSep(dijet.M()		,dijetDeltaEta,		zp4LEPScaleUP,		wp4LEPScaleUP	, zep3l,	type3l,	isHiggs)	,bound);
	

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
	  if(passVBFLooseJESDOWN) 	histo_WZ_CMS_scale_jDown	->Fill(MVAVarJESDOWN,totalWeight);
	  if(passVBFLooseJESUP) 	histo_WZ_CMS_scale_jUp		->Fill(MVAVarJESUP,totalWeight);
	  if(passVBFLooseJERDOWN) 	histo_WZ_CMS_res_jDown		->Fill(MVAVarJERDOWN,totalWeight);
	  if(passVBFLooseJERUP) 	histo_WZ_CMS_res_jUp		->Fill(MVAVarJERUP,totalWeight);
	  if(passVBFLooseLEPScaleUP) 	histo_WZ_CMS_scale_lUp		->Fill(MVAVarLEPScaleUP,totalWeight);
	  if(passVBFLooseLEPScaleDOWN) 	histo_WZ_CMS_scale_lDown	->Fill(MVAVarLEPScaleDOWN,totalWeight);
          
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


  
  
  //cout<<fakeShape->GetBinError(1)/fakeShape->GetBinContent(1)<<endl;
//cout<<fakeShape->GetBinError(0)/fakeShape->GetBinContent(0)<<endl;
  if(sampleID==5){ 
  double fakeArray[4] = {32.65/84.72, 10.75/84.72, 28.9/84.72, 12.42/84.72};
  double fakeArrayUnc[4] = {4.7/32.65, 3.47/10.75, 4.34/28.90, 3.43/12.42};
    for(int i=0; i<4; i++){
      for(int j=1;j<14;j++){
	//if(fakeShape)


	double RelBinError=(double)fakeShape->GetBinError(j)/fakeShape->GetBinContent(j);

	histo_Fake->SetBinContent(j+(i*13.),fakeArray[i]*fakeShape->GetBinContent(j));

	histo_Fake->SetBinError(j+(i*13.),fakeArray[i]*fakeShape->GetBinContent(j)*RelBinError);
	if(i==0){
	  histo_Fake_CMS_wz3l_FakeSys0Up->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j))+(fakeArray[i]*fakeShape->GetBinContent(j)*fakeArrayUnc[i]));
	  histo_Fake_CMS_wz3l_FakeSys0Down->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j))-(fakeArray[i]*fakeShape->GetBinContent(j)*fakeArrayUnc[i]));
	}
	else{
	  histo_Fake_CMS_wz3l_FakeSys0Up->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j)));
	  histo_Fake_CMS_wz3l_FakeSys0Down->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j)));
	}
       if(i==1){
	  histo_Fake_CMS_wz3l_FakeSys1Up->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j))+(fakeArray[i]*fakeShape->GetBinContent(j)*fakeArrayUnc[i]));
	  histo_Fake_CMS_wz3l_FakeSys1Down->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j))-(fakeArray[i]*fakeShape->GetBinContent(j)*fakeArrayUnc[i]));
	}
	else{
	  histo_Fake_CMS_wz3l_FakeSys1Up->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j)));
	  histo_Fake_CMS_wz3l_FakeSys1Down->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j)));
	}
	if(i==2){
	  histo_Fake_CMS_wz3l_FakeSys2Up->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j))+(fakeArray[i]*fakeShape->GetBinContent(j)*fakeArrayUnc[i]));
	  histo_Fake_CMS_wz3l_FakeSys2Down->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j))-(fakeArray[i]*fakeShape->GetBinContent(j)*fakeArrayUnc[i]));
	}
	else{
	  histo_Fake_CMS_wz3l_FakeSys2Up->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j)));
	  histo_Fake_CMS_wz3l_FakeSys2Down->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j)));
	}
	if(i==3){
	  histo_Fake_CMS_wz3l_FakeSys3Up->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j))+(fakeArray[i]*fakeShape->GetBinContent(j)*fakeArrayUnc[i]));
	  histo_Fake_CMS_wz3l_FakeSys3Down->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j))-(fakeArray[i]*fakeShape->GetBinContent(j)*fakeArrayUnc[i]));
	}
	else{
	  histo_Fake_CMS_wz3l_FakeSys3Up->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j)));
	  histo_Fake_CMS_wz3l_FakeSys3Down->SetBinContent(j+(i*13.),(fakeArray[i]*fakeShape->GetBinContent(j)));
	}
	
      	}
	}
    //delete fakeFile;
    
  histo_Fake  				->Write();
  histo_Fake_CMS_wz3l_FakeSys3Up->Write();
  histo_Fake_CMS_wz3l_FakeSys3Down->Write();
  histo_Fake_CMS_wz3l_FakeSys2Up->Write();
  histo_Fake_CMS_wz3l_FakeSys2Down->Write();
  histo_Fake_CMS_wz3l_FakeSys1Up->Write();
  histo_Fake_CMS_wz3l_FakeSys1Down->Write();
  histo_Fake_CMS_wz3l_FakeSys0Up->Write();
  histo_Fake_CMS_wz3l_FakeSys0Down->Write();
  histo_FakeE  				->Write();
  }

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
      newcardShape << Form("%s				     lnN  	%7.5f   %7.5f   %7.5    %7.5f    -    %7.5f    %7.5f   \n","CMS_eff_b_mistag_13TeV2016" ,0.98,0.98,0.98,0.98,0.98,0.98);
      newcardShape << Form("%s		                     lnN   	  -       -       -       -    %7.5f    -        -     \n","CMS_wz3l_FakeSys" 		,1.30);	
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n",effMName			,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n",effEName			,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n","CMS_pu" 			,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n","CMS_scale_j"		,1.,1.,1.,1.,1.,1.);
      newcardShape << Form("%s                               shape  	%7.5f   %7.5f   %7.5f   %7.5f   -     %7.5f    %7.5f   \n","CMS_res_j"			,1.,1.,1.,1.,1.,1.);
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
      newcardShape << Form("%s				   lnN  	%7.5f   %7.5f   %7.5f     %7.5f    -    %7.5f     \n","CMS_eff_b_mistag_13TeV2016",0.98,0.98,0.98,0.98,0.98);
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
      newcardShape << Form("%s		                   shape   	  -       -       -       -    %7.5f    -         \n","CMS_wz3l_FakeSys0" 		,1.);	
      newcardShape << Form("%s		                   shape   	  -       -       -       -    %7.5f    -         \n","CMS_wz3l_FakeSys1" 		,1.);	
      newcardShape << Form("%s		                   shape   	  -       -       -       -    %7.5f    -         \n","CMS_wz3l_FakeSys2" 		,1.);	
      newcardShape << Form("%s		                   shape   	  -       -       -       -    %7.5f    -         \n","CMS_wz3l_FakeSys3" 		,1.);	
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


