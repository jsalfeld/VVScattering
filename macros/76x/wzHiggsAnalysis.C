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

#include "MitAnalysisRunII/macros/76x/factors.h"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"

enum selType                     { SIGSEL, nSelTypes};
TString selTypeName[nSelTypes]= { "SIGSEL"};

enum systType                     {METUP=0, METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"METUP","METDOWN"};

double mcPrescale = 1.0;
bool usePureMC = false;

TRandom2 r;

void wzHiggsAnalysisNEWwithSysRenameSystematics(
 double minMass =  76,
 double maxMass = 106,
 bool applyBtagging = true,
 TString typeLepSel = "medium",
 TString type3rdLepSel = "default"
 ){


  TString SignalSuffix[11] = {"M200","M300","M400","M500","M600","M700","M800","M900","M1000","M1500","M2000"};

  Int_t period = 1;
   TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_76x/met_";
    TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_76x/met_";
  TString filesPathMCSig  = "/scratch/jsalfeld/VBFWZ_76X_13TeV_Signal/met_";


TString fileFakeShape ="MitAnalysisRunII/data/80x/fakeShape2016.root";


//TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_76x/met_";
//  TString filesPathDA  = "root://eoscms//eos/cms/store/user/jsalfeld/HWZ_Study/data/met_";
   //TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_76x/met_";
  //   TString filesPathMC  = "root://eoscms//eos/cms/store/user/jsalfeld/HWZ_Study/met_";
  Double_t lumi = 2.318;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  if      (period==1){
  puPath = "MitAnalysisRunII/data/76x/puWeights_76x.root";
        infilenamev.push_back(Form("%sdata_AOD_Run2015C_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);
     infilenamev.push_back(Form("%sdata_AOD_Run2015D_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);


  //  infilenamev.push_back("/scratch/ceballos/ntuples_weightsDA_80x/met_data.root");  infilecatv.push_back(0);

  //infilenamev.push_back("/scratch/ceballos/ntuples_weightsDA_80x/data.root");                                                                                          
  //infilecatv.push_back(0);


  if(usePureMC == true){
infilenamev.push_back(Form("%sTT_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1+AODSIM.root",filesPathMC.Data()));
    /*
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                                            infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));						   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(1);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                            infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 		   infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  				   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(1);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));*/



			   infilecatv.push_back(1);
  }
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                   infilecatv.push_back(2);

  //infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(3);

  //infilenamev.push_back(Form("%sWZJJ_EWK_QCD_13TeV-madgraph-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); infilecatv.push_back(3);

  //infilenamev.push_back(Form("%sWZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(3);

  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sWZJJ_EWK_QCD_13TeV-madgraph-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	       infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZJJ_EWK_13TeV-madgraph-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	       infilecatv.push_back(3);
  //  infilenamev.push_back(Form("%sWZJJ_QCD_13TeV-madgraph-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	       infilecatv.push_back(3);



  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3+AODSIM.root",filesPathMC.Data())); 			   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));       
  infilecatv.push_back(5);
 infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(5);
    //infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); infilecatv.push_back(5);



  //infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(5);
  //  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(5);
 // infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(5);
  infilenamev.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                     infilecatv.push_back(5);          
  //infilenamev.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                     infilecatv.push_back(5);
 

    for(int i=0; i<11;i++){
    infilenamev.push_back(Form("%sVBFHiggs_Signal"+SignalSuffix[i]+"_New.root",filesPathMCSig.Data()));           infilecatv.push_back(6);
	}
  
  //infilenamev.push_back(Form("%sM600_Signal.root",filesPathMC.Data()));           infilecatv.push_back(6);

  }
  else {assert(0);}
  //std::cout<<"hallo"<<std::endl;
  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}

  //  infilenamev.clear();infilecatv.clear();

  // infilenamev.push_back(Form("%sdata_AOD_Run2015D_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);
  //    infilenamev.push_back(Form("%sdata_AOD_Run2015C_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);
  // infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(5);
  // infilenamev.push_back("/scratch/ceballos/ntuples_weightsDA_80x/data.root");
  // infilecatv.push_back(0);
   
  //infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(3);
	//  infilenamev.push_back(Form("%sWZJJ_EWK_QCD_13TeV-madgraph-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sVBFHiggs_SignalM700_New.root",filesPathMCSig.Data()));           infilecatv.push_back(6);

  //for(int i=0; i<9;i++){
  //infilenamev.push_back(Form("%sVBFHiggs_Signal"+SignalSuffix[i]+"_New.root",filesPathMCSig.Data()));           infilecatv.push_back(6);
  //}

 //infilenamev.push_back(Form("%sWZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(3);

 // infilenamev.push_back(Form("%sM600_Signal2.root",filesPathMC.Data()));           infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sSignal200.root",filesPathMC.Data()));           infilecatv.push_back(6);


  //infilenamev.push_back(Form("%sWZ.root",filesPathMC.Data()));           infilecatv.push_back(4);

  //infilenamev.push_back(Form("%sWZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(3);



   //  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(3);



  //infilenamev.push_back("root://eoscms//eos/cms/store/user/jsalfeld/WZpowheg/NeroNtuples_8.root");infilecatv.push_back(4);

  //infilenamev.push_back(Form("%sWZJJ_EWK_QCD_13TeV-madgraph-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	       infilecatv.push_back(5);
  // infilenamev.push_back(Form("%sWZJJ_EWK_13TeV-madgraph-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	       infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sWZJJ_QCD_13TeV-madgraph-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	       infilecatv.push_back(5);

  //infilenamev.push_back(Form("%sM600_Signal2.root",filesPathMC.Data()));           infilecatv.push_back(6);

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
  cout<<"hi"<<endl;
  TFile *fFakeShape = TFile::Open(fileFakeShape.Data());
  TH1D *fhDfakeShape = (TH1D*) (fFakeShape->Get("histo_Fake")); assert(fhDfakeShape); fhDfakeShape->SetDirectory(0);
  delete fFakeShape;

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
  delete fPUFile;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Medium_ele"));
  TH2D *fhDElTightSF  = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Tight_ele"));
  TH2D *fhDElMediumMVASF = (TH2D*)(fElSF->Get("unfactorized_scalefactors_MediumMVA_ele"));
  TH2D *fhDElTightMVASF  = (TH2D*)(fElSF->Get("unfactorized_scalefactors_TightMVA_ele"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  assert(fhDElMediumMVASF);
  assert(fhDElTightMVASF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF ->SetDirectory(0);
  fhDElMediumMVASF->SetDirectory(0);
  fhDElTightMVASF ->SetDirectory(0);
  delete fElSF;

  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Medium_mu"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Iso_mu"));
  assert(fhDMuMediumSF);
  assert(fhDMuIsoSF);
  fhDMuMediumSF->SetDirectory(0);
  fhDMuIsoSF->SetDirectory(0);
  delete fMuSF;

  double totalFakeDataCount[4][9];
  for(int i=0; i<4; i++) for(int j=0; j<9; j++) totalFakeDataCount[i][j] = 0;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 56;
  const int histBins = 7;
  const int allStates = 6;
  TH1D* histo[allStates][allPlots][histBins];
  TString processName[histBins] = {"..Data", ".Fakes", "Zgamma", "....WZ", "...ZZ", "...VVV","...Signal"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  2) {nBinPlot = 170; xminPlot = 30.0; xmaxPlot = 200.0;}
    else if(thePlot >=  3 && thePlot <=  3) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 150.0;}
    else if(thePlot >=  4 && thePlot <=  4) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  5 && thePlot <=  5) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400.0;}
    else if(thePlot >=  6 && thePlot <=  6) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot =   4.0;}
    else if(thePlot >=  8 && thePlot <= 11) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 12 && thePlot <= 12) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 13 && thePlot <= 13) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot =  40.0;}
    else if(thePlot >= 14 && thePlot <= 14) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot >= 15 && thePlot <= 15) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >= 16 && thePlot <= 17) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 18 && thePlot <= 18) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 20 && thePlot <= 20) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot >= 21 && thePlot <= 22) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 1000.0;}
    else if(thePlot >= 23 && thePlot <= 24) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =1000.0;}
    else if(thePlot >= 25 && thePlot <= 25) {nBinPlot =  50; xminPlot =-0.5; xmaxPlot =  49.5;}
    else if(thePlot >= 26 && thePlot <= 31) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 32 && thePlot <= 32) {nBinPlot = 100; xminPlot =50.0; xmaxPlot = 250.0;}
    else if(thePlot >= 33 && thePlot <= 36) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >= 37 && thePlot <= 37) {nBinPlot = 120; xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >= 38 && thePlot <= 38) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =2000.0;}
    else if(thePlot >= 39 && thePlot <= 39) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot =   8.0;}
    else if(thePlot >= 40 && thePlot <= 40) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >= 41 && thePlot <= 42) {nBinPlot =   50;   xminPlot = 0.0; xmaxPlot =   500;}
    else if(thePlot >= 43 && thePlot <= 44) {nBinPlot =   10;  xminPlot =-5.0; xmaxPlot =   5.0;}
    else if(thePlot >= 45 && thePlot <= 45) {nBinPlot =   40;   xminPlot = 100.0; xmaxPlot =   500;}
    else if(thePlot >= 46 && thePlot <= 46) {nBinPlot =   25;   xminPlot = 0.0; xmaxPlot =   2.5;}
    else if(thePlot >= 47 && thePlot <= 49) {nBinPlot =   20;  xminPlot = 0.0; xmaxPlot =   200;}
    else if(thePlot >= 50 && thePlot <= 52) {nBinPlot =   50;   xminPlot =-2.5; xmaxPlot =   2.5;}
    else if(thePlot >= 53 && thePlot <= 53) {nBinPlot =   4;    xminPlot =-0.5; xmaxPlot =   3.5;}
    else if(thePlot >= 54 && thePlot <= 54) {nBinPlot =   200;    xminPlot =0.0; xmaxPlot =   1000;}
    else if(thePlot >= 55 && thePlot <= 55) {nBinPlot =   20;    xminPlot =0.0; xmaxPlot =   200;}
    else if(thePlot >= 56 && thePlot <= 56) {nBinPlot =   50;    xminPlot =0.0; xmaxPlot =   500;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) {
      for(int j=0; j<allStates; j++) histo[j][thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    }
    histos->Reset();histos->Clear();
  }

    TH1D* MuonPt = new TH1D("MuonPt", "MuonPt", 50, 0., 200);
    MuonPt->Sumw2();
TH1D* ElePt = new TH1D("ElePt", "ElePt", 50, 0., 200);
    ElePt->Sumw2();

    TH2D* Higgs_2d = new TH2D("Higgs_2d", "Higgs_2d", 100, 0., 1000, 100, 0,1000);
    Higgs_2d->Sumw2();
    TH2D* WZ_2d = new TH2D("WZ_2d", "WZ_2d", 100, 0., 1000, 100,0,1000);
    Higgs_2d->Sumw2();

  TString ECMsb  = "13TeV2015";
  //const int nBinMVA = 4; Float_t xbins[nBinMVA+1] = {0, 100, 250, 450, 1500};
  const int nBinMVA = 9; Float_t xbins[nBinMVA+1] = {0, 100, 200, 400, 600, 800, 1000, 1200, 1500, 2200};
   // const int nBinMVA = 20; Float_t xbins[nBinMVA+1] = {0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900,950,1000};

  //const int nBinMVA = 1; Float_t xbins[nBinMVA+1] = {0, 1};
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();

  TH1D *histo_Data   = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_Zg     = (TH1D*) histoMVA->Clone("histo_Zg"); 
  TH1D *histo_VVV    = (TH1D*) histoMVA->Clone("histo_VVV");	 
  TH1D *histo_WZ     = (TH1D*) histoMVA->Clone("histo_WZ");	 
  TH1D *histo_ZZ     = (TH1D*) histoMVA->Clone("histo_ZZ");	 
  TH1D *histo_Fake  = (TH1D*) histoMVA->Clone("histo_Fake");	 
  TH1D *histo_FakeE  = (TH1D*) histoMVA->Clone("histo_FakeE");	 
  TH1D *histo_Higgs  = (TH1D*) histoMVA->Clone("histo_Higgs");	
  TH1D *histo_Higgs_[11];
  for(int i = 0; i<11 ; i++)
	{
    		histo_Higgs_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]);	 
	}

  char finalStateName[2],effMName[10],effEName[10],momMName[10],momEName[10];
  sprintf(effMName,"CMS_eff_m");sprintf(momMName,"CMS_scale_m");
  sprintf(effEName,"CMS_eff_e");sprintf(momEName,"CMS_scale_e");
  sprintf(finalStateName,"3l");


  TH1D* histo_Zg_CMS_scale_jUp            = 	new TH1D( "histo_Zg_CMS_scale_jUp"  , "histo_Zg_CMS_scale_jUp"  , nBinMVA, xbins);		histo_Zg_CMS_scale_jUp  ->Sumw2();
  TH1D* histo_Zg_CMS_scale_jDown          = 	new TH1D( "histo_Zg_CMS_scale_jDown"  , "histo_Zg_CMS_scale_jDown"  , nBinMVA, xbins); 	histo_Zg_CMS_scale_jDown  ->Sumw2();
  TH1D* histo_VVV_CMS_scale_jUp           = 	new TH1D( "histo_VVV_CMS_scale_jUp"  , "histo_VVV_CMS_scale_jUp"  , nBinMVA, xbins); 	histo_VVV_CMS_scale_jUp  ->Sumw2();
  TH1D* histo_VVV_CMS_scale_jDown         = 	new TH1D( "histo_VVV_CMS_scale_jDown"  , "histo_VVV_CMS_scale_jDown"  , nBinMVA, xbins); 	histo_VVV_CMS_scale_jDown  ->Sumw2();
  TH1D* histo_WZ_CMS_scale_jUp           = 		new TH1D( "histo_WZ_CMS_scale_jUp"  , "histo_WZ_CMS_scale_jUp"  , nBinMVA, xbins); 		histo_WZ_CMS_scale_jUp  ->Sumw2();
  TH1D* histo_WZ_CMS_scale_jDown           = 	new TH1D( "histo_WZ_CMS_scale_jDown"  , "histo_WZ_CMS_scale_jDown"  , nBinMVA, xbins); 	histo_WZ_CMS_scale_jDown  ->Sumw2();
  TH1D* histo_ZZ_CMS_scale_jUp           = 		new TH1D( "histo_ZZ_CMS_scale_jUp"  , "histo_ZZ_CMS_scale_jUp"  , nBinMVA, xbins); 		histo_ZZ_CMS_scale_jUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_scale_jDown           = 	new TH1D( "histo_ZZ_CMS_scale_jDown"  , "histo_ZZ_CMS_scale_jDown"  , nBinMVA, xbins); 	histo_ZZ_CMS_scale_jDown  ->Sumw2();
  TH1D *histo_Higgs_CMS_scale_jUp_[11];
  TH1D *histo_Higgs_CMS_scale_jDown_[11];
  for(int i = 0; i<11 ; i++)
	{
    		histo_Higgs_CMS_scale_jUp_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_scale_jUp");	 
		histo_Higgs_CMS_scale_jDown_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_scale_jDown");	 
	}

  TH1D* histo_Zg_CMS_res_metUp            = 	new TH1D( "histo_Zg_CMS_res_metUp"  , "histo_Zg_CMS_res_metUp"  , nBinMVA, xbins);		histo_Zg_CMS_res_metUp  ->Sumw2();
  TH1D* histo_Zg_CMS_res_metDown          = 	new TH1D( "histo_Zg_CMS_res_metDown"  , "histo_Zg_CMS_res_metDown"  , nBinMVA, xbins); 	histo_Zg_CMS_res_metDown  ->Sumw2();
  TH1D* histo_VVV_CMS_res_metUp           = 	new TH1D( "histo_VVV_CMS_res_metUp"  , "histo_VVV_CMS_res_metUp"  , nBinMVA, xbins); 	histo_VVV_CMS_res_metUp  ->Sumw2();
  TH1D* histo_VVV_CMS_res_metDown         = 	new TH1D( "histo_VVV_CMS_res_metDown"  , "histo_VVV_CMS_res_metDown"  , nBinMVA, xbins); 	histo_VVV_CMS_res_metDown  ->Sumw2();
  TH1D* histo_WZ_CMS_res_metUp           = 		new TH1D( "histo_WZ_CMS_res_metUp"  , "histo_WZ_CMS_res_metUp"  , nBinMVA, xbins); 		histo_WZ_CMS_res_metUp  ->Sumw2();
  TH1D* histo_WZ_CMS_res_metDown           = 	new TH1D( "histo_WZ_CMS_res_metDown"  , "histo_WZ_CMS_res_metDown"  , nBinMVA, xbins); 	histo_WZ_CMS_res_metDown  ->Sumw2();
  TH1D* histo_ZZ_CMS_res_metUp           = 		new TH1D( "histo_ZZ_CMS_res_metUp"  , "histo_ZZ_CMS_res_metUp"  , nBinMVA, xbins); 		histo_ZZ_CMS_res_metUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_res_metDown           = 	new TH1D( "histo_ZZ_CMS_res_metDown"  , "histo_ZZ_CMS_res_metDown"  , nBinMVA, xbins); 	histo_ZZ_CMS_res_metDown  ->Sumw2();
  TH1D *histo_Higgs_CMS_res_metUp_[11];
  TH1D *histo_Higgs_CMS_res_metDown_[11];
  for(int i = 0; i<11 ; i++)
	{
    		histo_Higgs_CMS_res_metUp_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_res_metUp");	 
		histo_Higgs_CMS_res_metDown_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_res_metDown");	 
	}

  TH1D* histo_Zg_CMS_scale_lUp            = 	new TH1D( "histo_Zg_CMS_scale_lUp"  , "histo_Zg_CMS_scale_lUp"  , nBinMVA, xbins);		histo_Zg_CMS_scale_lUp  ->Sumw2();
  TH1D* histo_Zg_CMS_scale_lDown          = 	new TH1D( "histo_Zg_CMS_scale_lDown"  , "histo_Zg_CMS_scale_lDown"  , nBinMVA, xbins); 	histo_Zg_CMS_scale_lDown  ->Sumw2();
  TH1D* histo_VVV_CMS_scale_lUp           = 	new TH1D( "histo_VVV_CMS_scale_lUp"  , "histo_VVV_CMS_scale_lUp"  , nBinMVA, xbins); 	histo_VVV_CMS_scale_lUp  ->Sumw2();
  TH1D* histo_VVV_CMS_scale_lDown         = 	new TH1D( "histo_VVV_CMS_scale_lDown"  , "histo_VVV_CMS_scale_lDown"  , nBinMVA, xbins); 	histo_VVV_CMS_scale_lDown  ->Sumw2();
  TH1D* histo_WZ_CMS_scale_lUp           = 		new TH1D( "histo_WZ_CMS_scale_lUp"  , "histo_WZ_CMS_scale_lUp"  , nBinMVA, xbins); 		histo_WZ_CMS_scale_lUp  ->Sumw2();
  TH1D* histo_WZ_CMS_scale_lDown           = 	new TH1D( "histo_WZ_CMS_scale_lDown"  , "histo_WZ_CMS_scale_lDown"  , nBinMVA, xbins); 	histo_WZ_CMS_scale_lDown  ->Sumw2();
  TH1D* histo_ZZ_CMS_scale_lUp           = 		new TH1D( "histo_ZZ_CMS_scale_lUp"  , "histo_ZZ_CMS_scale_lUp"  , nBinMVA, xbins); 		histo_ZZ_CMS_scale_lUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_scale_lDown           = 	new TH1D( "histo_ZZ_CMS_scale_lDown"  , "histo_ZZ_CMS_scale_lDown"  , nBinMVA, xbins); 	histo_ZZ_CMS_scale_lDown  ->Sumw2();
  TH1D *histo_Higgs_CMS_scale_lUp_[11];
  TH1D *histo_Higgs_CMS_scale_lDown_[11];
  for(int i = 0; i<11 ; i++)
	{
    		histo_Higgs_CMS_scale_lUp_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_scale_lUp");	 
		histo_Higgs_CMS_scale_lDown_[i] = (TH1D*) histoMVA->Clone("histo_Higgs_"+SignalSuffix[i]+"_CMS_scale_lDown");	 
	}

  TH1D* histo_WZ_CMS_ewkscaleUp           =             new TH1D( "histo_WZ_CMS_ewkscaleUp"  , "histo_WZ_CMS_ewkscaleUp"  , nBinMVA, xbins);            histo_WZ_CMS_ewkscaleUp  ->Sumw2();
  TH1D* histo_WZ_CMS_ewkscaleDown           =           new TH1D( "histo_WZ_CMS_ewkscaleDown"  , "histo_WZ_CMS_ewkscaleDown"  , nBinMVA, xbins);                histo_WZ_CMS_ewkscaleDown  ->Sumw2();

  TH1D* histo_WZ_CMS_QCDunc_scaleUp           =                 new TH1D( "histo_WZ_CMS_QCDunc_scaleUp"  , "histo_WZ_CMS_QCDunc_scaleUp"  , nBinMVA, xbins);            histo_WZ_CMS_QCDunc_scaleUp  ->Sumw2();
  TH1D* histo_WZ_CMS_QCDunc_scaleDown           =               new TH1D( "histo_WZ_CMS_QCDunc_scaleDown"  , "histo_WZ_CMS_QCDunc_scaleDown"  , nBinMVA, xbins);                histo_WZ_CMS_QCDunc_scaleDown  ->Sumw2();




  TH1D* histo_Zg_CMS_MVAZHStatBoundingUp           = new TH1D( Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVAZHStatBoundingDown         = new TH1D( Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp         = new TH1D( Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown       = new TH1D( Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingUp           = new TH1D( Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingDown         = new TH1D( Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingUp           = new TH1D( Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingDown         = new TH1D( Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingDown->Sumw2();
  TH1D* histo_Fake_CMS_MVAFakeMStatBoundingUp     = new TH1D( Form("histo_Fake_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Fake_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Fake_CMS_MVAFakeMStatBoundingUp  ->Sumw2();
  TH1D* histo_Fake_CMS_MVAFakeMStatBoundingDown   = new TH1D( Form("histo_Fake_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Fake_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Fake_CMS_MVAFakeMStatBoundingDown->Sumw2();
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingUp     = new TH1D( Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->Sumw2();
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingDown   = new TH1D( Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_FakeE_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingDown->Sumw2();

  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingUp     = new TH1D( Form("histo_Higgs_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_wz3l%s_MVAEMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingDown   = new TH1D( Form("histo_Higgs_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_wz3l%s_MVAEMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingDown->Sumw2();

  TH1D* histo_Diff = new TH1D("dummy", "dummy",1000,-1,1); histo_Diff->Sumw2();

  double histo_WZ_CMS_QCDScaleInitial[7] = {0,0,0,0,0,0,0};
  TH1D* histo_Zg_CMS_QCDScaleBounding[6];
  TH1D* histo_VVV_CMS_QCDScaleBounding[6];
  TH1D* histo_WZ_CMS_QCDScaleBounding[6];
  TH1D* histo_ZZ_CMS_QCDScaleBounding[6];
  TH1D* histo_Higgs_CMS_QCDScaleBounding[6];

  for(int nb=0; nb<6; nb++){
    histo_Zg_CMS_QCDScaleBounding[nb]      = new TH1D(Form("histo_Zg_QCDScale_f%d",nb),      Form("histo_Zg_QCDScale_f%d",nb),nBinMVA, xbins);      histo_Zg_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_VVV_CMS_QCDScaleBounding[nb]     = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),     Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WZ_CMS_QCDScaleBounding[nb]	   = new TH1D(Form("histo_WZ_QCDScale_f%d",nb),      Form("histo_WZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_WZ_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ZZ_CMS_QCDScaleBounding[nb]	   = new TH1D(Form("histo_ZZ_QCDScale_f%d",nb),      Form("histo_ZZ_QCDScale_f%d",nb),nBinMVA, xbins);      histo_ZZ_CMS_QCDScaleBounding[nb]->Sumw2();
histo_Higgs_CMS_QCDScaleBounding[nb]	   = new TH1D(Form("histo_Higgs_QCDScale_f%d",nb),      Form("histo_Higgs_QCDScale_f%d",nb),nBinMVA, xbins);      histo_Higgs_CMS_QCDScaleBounding[nb]->Sumw2();
  }
  TH1D* histo_Zg_CMS_PDFBounding[102];
  TH1D* histo_VVV_CMS_PDFBounding[102];
  TH1D* histo_WZ_CMS_PDFBounding[102];
  TH1D* histo_ZZ_CMS_PDFBounding[102];
  TH1D* histo_Higgs_CMS_PDFBounding[102];
  for(int nb=0; nb<102; nb++){
    histo_Zg_CMS_PDFBounding[nb]      = new TH1D(Form("histo_Zg_PDF_f%d",nb),      Form("histo_Zg_PDF_f%d",nb),     nBinMVA, xbins); histo_Zg_CMS_PDFBounding[nb]->Sumw2();
    histo_VVV_CMS_PDFBounding[nb]     = new TH1D(Form("histo_VVV_PDF_f%d",nb),     Form("histo_VVV_PDF_f%d",nb),    nBinMVA, xbins); histo_VVV_CMS_PDFBounding[nb]->Sumw2();
    histo_WZ_CMS_PDFBounding[nb]      = new TH1D(Form("histo_WZ_PDF_f%d",nb),      Form("histo_WZ_PDF_f%d",nb),     nBinMVA, xbins); histo_WZ_CMS_PDFBounding[nb]->Sumw2();
    histo_ZZ_CMS_PDFBounding[nb]      = new TH1D(Form("histo_ZZ_PDF_f%d",nb),      Form("histo_ZZ_PDF_f%d",nb),     nBinMVA, xbins); histo_ZZ_CMS_PDFBounding[nb]->Sumw2();
    histo_Higgs_CMS_PDFBounding[nb]      = new TH1D(Form("histo_Higgs_PDF_f%d",nb),      Form("histo_Higgs_PDF_f%d",nb),     nBinMVA, xbins); histo_Higgs_CMS_PDFBounding[nb]->Sumw2();
  }

  TH1D* histo_Zg_CMS_MVAZHStatBoundingBinUp[nBinMVA];
  TH1D* histo_Zg_CMS_MVAZHStatBoundingBinDown[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinUp[nBinMVA];
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinDown[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinUp[nBinMVA];
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinDown[nBinMVA];
  TH1D* histo_Fake_CMS_MVAFakeMStatBoundingBinUp[nBinMVA];
  TH1D* histo_Fake_CMS_MVAFakeMStatBoundingBinDown[nBinMVA];
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nBinMVA];
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++){
    histo_Zg_CMS_MVAZHStatBoundingBinUp[nb]         = new TH1D(Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dUp"   ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingBinUp[nb]   ->Sumw2();
    histo_Zg_CMS_MVAZHStatBoundingBinDown[nb]       = new TH1D(Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb), Form("histo_Zg_CMS_wz3l%s_MVAZHStatBounding_%s_Bin%dDown" ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Zg_CMS_MVAZHStatBoundingBinDown[nb] ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	    = new TH1D(Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	    = new TH1D(Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_wz3l%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	    = new TH1D(Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WZ_CMS_wz3l%s_MVAWZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	    = new TH1D(Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dUp"        ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	      ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	    = new TH1D(Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb), Form("histo_ZZ_CMS_wz3l%s_MVAZZStatBounding_%s_Bin%dDown"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]      ->Sumw2();
    histo_Fake_CMS_MVAFakeMStatBoundingBinUp[nb]   = new TH1D(Form("histo_Fake_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_Fake_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Fake_CMS_MVAFakeMStatBoundingBinUp[nb]       ->Sumw2();
    histo_Fake_CMS_MVAFakeMStatBoundingBinDown[nb] = new TH1D(Form("histo_Fake_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_Fake_CMS_wz3l%s_MVAFakeMStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Fake_CMS_MVAFakeMStatBoundingBinDown[nb]	 ->Sumw2();
    histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nb]   = new TH1D(Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nb]       ->Sumw2();
    histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nb] = new TH1D(Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_FakeE_CMS_wz3l%s_MVAFakeEStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nb]	 ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]   = new TH1D(Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb), Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%s_Bin%dUp"  ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]       ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb] = new TH1D(Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb), Form("histo_Higgs_CMS_wz3l%s_MVAHiggsStatBounding_%s_Bin%dDown",finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb]	 ->Sumw2();
  }

  TH1D* histo_Zg_CMS_MVALepEffMBoundingUp        = new TH1D( Form("histo_Zg_%sUp",effMName)  , Form("histo_Zg_%sUp",effMName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVALepEffMBoundingDown      = new TH1D( Form("histo_Zg_%sDown",effMName), Form("histo_Zg_%sDown",effMName), nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effMName)  , Form("histo_VVV_%sUp",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effMName), Form("histo_VVV_%sDown",effMName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_WZ_%sUp",effMName)  , Form("histo_WZ_%sUp",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_WZ_%sDown",effMName), Form("histo_WZ_%sDown",effMName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_ZZ_%sUp",effMName)  , Form("histo_ZZ_%sUp",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_ZZ_%sDown",effMName), Form("histo_ZZ_%sDown",effMName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingDown->Sumw2();
TH1D* histo_Higgs_CMS_MVALepEffMBoundingUp    	 = new TH1D( Form("histo_Higgs_%sUp",effMName)  , Form("histo_Higgs_%sUp",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingDown  	 = new TH1D( Form("histo_Higgs_%sDown",effMName), Form("histo_Higgs_%sDown",effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_MVALepEffMBoundingAvg       = new TH1D( Form("histo_Zg_%sAvg",effMName)  , Form("histo_Zg_%sAvg",effMName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effMName)  , Form("histo_VVV_%sAvg",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_WZ_%sAvg",effMName)  , Form("histo_WZ_%sAvg",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_ZZ_%sAvg",effMName)  , Form("histo_ZZ_%sAvg",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
TH1D* histo_Higgs_CMS_MVALepEffMBoundingAvg    	 = new TH1D( Form("histo_Higgs_%sAvg",effMName)  , Form("histo_Higgs_%sAvg",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingAvg  ->Sumw2();

  TH1D* histo_Zg_CMS_MVALepEffEBoundingUp        = new TH1D( Form("histo_Zg_%sUp",effEName)  , Form("histo_Zg_%sUp",effEName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVALepEffEBoundingDown      = new TH1D( Form("histo_Zg_%sDown",effEName), Form("histo_Zg_%sDown",effEName), nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effEName)  , Form("histo_VVV_%sUp",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effEName), Form("histo_VVV_%sDown",effEName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_WZ_%sUp",effEName)  , Form("histo_WZ_%sUp",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_WZ_%sDown",effEName), Form("histo_WZ_%sDown",effEName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_ZZ_%sUp",effEName)  , Form("histo_ZZ_%sUp",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_ZZ_%sDown",effEName), Form("histo_ZZ_%sDown",effEName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingUp    	 = new TH1D( Form("histo_Higgs_%sUp",effEName)  , Form("histo_Higgs_%sUp",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingDown  	 = new TH1D( Form("histo_Higgs_%sDown",effEName), Form("histo_Higgs_%sDown",effEName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_MVALepEffEBoundingAvg       = new TH1D( Form("histo_Zg_%sAvg",effEName)  , Form("histo_Zg_%sAvg",effEName)  , nBinMVA, xbins); histo_Zg_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effEName)  , Form("histo_VVV_%sAvg",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_WZ_%sAvg",effEName)  , Form("histo_WZ_%sAvg",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_ZZ_%sAvg",effEName)  , Form("histo_ZZ_%sAvg",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingAvg    	 = new TH1D( Form("histo_Higgs_%sAvg",effEName)  , Form("histo_Higgs_%sAvg",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingAvg  ->Sumw2();

  TH1D* histo_Zg_CMS_MVAMETBoundingUp           = new TH1D( Form("histo_Zg_CMS_scale_metUp")  , Form("histo_Zg_CMS_scale_metUp")  , nBinMVA, xbins); histo_Zg_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_MVAMETBoundingDown         = new TH1D( Form("histo_Zg_CMS_scale_metDown"), Form("histo_Zg_CMS_scale_metDown"), nBinMVA, xbins); histo_Zg_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingUp   	= new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingDown 	= new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown"), nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_WZ_CMS_scale_metUp")  , Form("histo_WZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_WZ_CMS_scale_metDown"), Form("histo_WZ_CMS_scale_metDown"), nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_ZZ_CMS_scale_metUp")  , Form("histo_ZZ_CMS_scale_metUp")  , nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_ZZ_CMS_scale_metDown"), Form("histo_ZZ_CMS_scale_metDown"), nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingUp    	= new TH1D( Form("histo_Higgs_CMS_scale_metUp")  , Form("histo_Higgs_CMS_scale_metUp")  , nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingDown  	= new TH1D( Form("histo_Higgs_CMS_scale_metDown"), Form("histo_Higgs_CMS_scale_metDown"), nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingDown->Sumw2();

  TH1D* histo_Zg_CMS_PUBoundingUp               = new TH1D( Form("histo_Zg_CMS_puUp")  , Form("histo_Zg_CMS_puUp")  , nBinMVA, xbins); histo_Zg_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_Zg_CMS_PUBoundingDown             = new TH1D( Form("histo_Zg_CMS_puDown"), Form("histo_Zg_CMS_puDown"), nBinMVA, xbins); histo_Zg_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingUp           	= new TH1D( Form("histo_VVV_CMS_puUp")  , Form("histo_VVV_CMS_puUp")  , nBinMVA, xbins); histo_VVV_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingDown         	= new TH1D( Form("histo_VVV_CMS_puDown"), Form("histo_VVV_CMS_puDown"), nBinMVA, xbins); histo_VVV_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_PUBoundingUp            	= new TH1D( Form("histo_WZ_CMS_puUp")  , Form("histo_WZ_CMS_puUp")  , nBinMVA, xbins); histo_WZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_PUBoundingDown  	        = new TH1D( Form("histo_WZ_CMS_puDown"), Form("histo_WZ_CMS_puDown"), nBinMVA, xbins); histo_WZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingUp    	        = new TH1D( Form("histo_ZZ_CMS_puUp")  , Form("histo_ZZ_CMS_puUp")  , nBinMVA, xbins); histo_ZZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingDown  	        = new TH1D( Form("histo_ZZ_CMS_puDown"), Form("histo_ZZ_CMS_puDown"), nBinMVA, xbins); histo_ZZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingUp    	        = new TH1D( Form("histo_Higgs_CMS_puUp")  , Form("histo_Higgs_CMS_puUp")  , nBinMVA, xbins); histo_Higgs_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingDown  	        = new TH1D( Form("histo_Higgs_CMS_puDown"), Form("histo_Higgs_CMS_puDown"), nBinMVA, xbins); histo_Higgs_CMS_PUBoundingDown->Sumw2();
  unsigned int numberOfLeptons = 3;

  double totalEventsProcess[50];
  std::vector<double> sumEventsProcess(infilenamev.size(), 0.0);

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    std::cout<<"hallo1"<<std::endl;
    /*TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");
    TTree *the_PDF_tree   = (TTree*)the_input_file.FindObjectAny("pdfReweight");*/

    TFile *the_input_file=TFile::Open(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");
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
    eventMonteCarlo.setBranchAddresses(the_input_tree);

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
      printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());
    }

    char weightDef[256];
    int initPDFTag = -1;
    if(the_PDF_tree) {
      the_PDF_tree->SetBranchAddress("weightDef", &weightDef);
      for (int i=0; i<the_PDF_tree->GetEntries(); ++i) {
        the_PDF_tree->GetEntry(i);
        char **tokensPDF;
        size_t numtokensPDF;
        tokensPDF = strsplit(weightDef, " = ", &numtokensPDF);
        for (int k = 0; k < (int)numtokensPDF; k++) if(strcmp(tokensPDF[k],"292201") == 0||
	                                               strcmp(tokensPDF[k],"292001") == 0||
						       strcmp(tokensPDF[k],"260001") == 0) {initPDFTag = i; break;}
	if(initPDFTag != -1) break;
      }
    }
    if(infilecatv[ifile] != 0 && initPDFTag == -1 && infilenamev[ifile].Contains("powheg") == false) {
      printf("PDFTAG PROBLEM\n");
      if(the_PDF_tree) {
        printf("PDFTree Entries: %d\n",(int)the_PDF_tree->GetEntries());
        for (int i=0; i<the_PDF_tree->GetEntries(); ++i) {
          the_PDF_tree->GetEntry(i);
	  printf("PDF(%d): %s\n",i,weightDef);
        }
      }
      else {
        printf("PDFTree not available\n");
      }
      //return;
    }

    double theMCPrescale = mcPrescale;
    int numberOfEvents = int(the_input_tree->GetEntries()/theMCPrescale);
    //int nProcess = 0;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    if((infilecatv[ifile] > 5)) numberOfEvents = 50000; 
    for (int i=0; i < numberOfEvents; ++i) {
      the_input_tree->GetEntry(i);
      
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
      //if(infilecatv[ifile]==3 && i%6!=0) continue;
      //nProcess++;
      if(i%1000000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      Bool_t passFilter[11] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
      for (int nt = 0; nt <(int)numtokens; nt++) {
        if((*eventTrigger.triggerFired)[nt] == 0) continue;
        
	if((strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
           (strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
           (strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*")     == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")    == 0) ||
           (strcmp(tokens[nt],"HLT_IsoMu27_v*")     == 0) ||
           (strcmp(tokens[nt],"HLT_IsoMu20_v*")     == 0) ||
           (strcmp(tokens[nt],"HLT_IsoTkMu20_v*")    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele23_WPLoose_Gsf_v*")    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele22_eta2p1_WP75_Gsf_v*")    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele27_WPLoose_Gsf_v*")    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele27_WP85_Gsf_v*")      == 0)
           ) passFilter[1] = kTRUE;

      }
      
      //std::cout<<"hallo2"<<std::endl;
      	if(passFilter[0] == kFALSE) continue;
	//////////////////////////////////
	/////////////// CCCCCCHHHHHHHEEEEEECKKKKKK  /////////////  
	/////////////////////////////////
	if(infilecatv[ifile] != 6){
	if(passFilter[1] == kFALSE) continue;}
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
	/*printf("%s  %d  %f  %f  %f  %i  %f  \n",
	       typeLepSel.Data(),
	TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),
	TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),

	(double)(*eventLeptons.iso)[nlep],
	(int)(*eventLeptons.selBits)[nlep],
	(double)(*eventLeptons.mva)[nlep]);*/

        if(selectIdIsoCut(typeLepSel.Data(),
	TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if( (((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) ){idTight.push_back(0); idLep.push_back(nlep); if(numberOfLeptons == 4) goodIsTight++;}
        else if( (((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP )){idSoft.push_back(nlep);}
      }
      if(idLep.size()!=idTight.size()) {assert(1); return;}
      if(idLep.size()==numberOfLeptons) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;
      //std::cout<<"hallo"<<std::endl;
      if(goodIsTight == idTight.size()) passFilter[3] = kTRUE;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() <= 10) continue;


      for(int i=0; i<3; i++){

	if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[i]])==13){
	  MuonPt->Fill(((TLorentzVector*)(*eventLeptons.p4)[idLep[i]])->Pt());}    
	else{ElePt->Fill(((TLorentzVector*)(*eventLeptons.p4)[idLep[i]])->Pt());}
      }

      double dPhiLepMETMin = 999.;
      int signQ = 0;
      double systTotLep[2] = {1.0, 1.0}; // m/e
      for(unsigned nl=0; nl<idLep.size(); nl++){
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) systTotLep[0] = systTotLep[0] * 1.005;
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 11) systTotLep[1] = systTotLep[1] * 1.015;
        signQ = signQ + (int)(*eventLeptons.pdgId)[idLep[nl]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]);
        if(dPhiLepMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0]))))
           dPhiLepMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));      
      }
      // std::cout<<(double)eventMet.trackMet->Pt()<<"   "<<((TLorentzVector*)(*eventMet.p4)[0])->Pt()<<std::endl;
      double minMET  = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      double minPMET = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      if(dPhiLepMETMin < TMath::Pi()/2) minPMET = minPMET * sin(dPhiLepMETMin);
      // std::cout<<"hi3"<<std::endl;
      passFilter[4] = TMath::Abs(signQ) == 1;
      if(passFilter[4] == kFALSE) continue;

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
 
      vector<int> idJet; vector<int> idJetJESUP; vector<int> idJetJESDOWN; double btagjet[2] = {0., 0.};
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      double theHT = 0;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;        

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(dPhiJetMET   == -1) dPhiJetMET   = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 15 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() >= 30) {
	
        if     ((float)(*eventJets.bDiscr)[nj] > btagjet[0]) {btagjet[1] = btagjet[0]; btagjet[0] = (float)(*eventJets.bDiscr)[nj];}
	else if((float)(*eventJets.bDiscr)[nj] > btagjet[1]) {btagjet[1] = (float)(*eventJets.bDiscr)[nj];}

        theHT = theHT + ((TLorentzVector*)(*eventJets.p4)[nj])->Pt();
	idJet.push_back(nj);}

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()* 0.95 >= 30) idJetJESDOWN.push_back(nj);
	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()* 1.05 >= 30) idJetJESUP.push_back(nj);


      }

      for(unsigned nl0=0; nl0<idLep.size(); nl0++){
        if((int)nl0==tagZ[0]||(int)nl0==tagZ[1]) continue;
        tagZ[2] = nl0;
        break;
      }
      
      if(tagZ[0] == -1 || tagZ[1] == -1 || tagZ[2] == -1) continue;
      
      bool tight3rdLepId = true;
      if(idTight[tagZ[2]] == 1){
        tight3rdLepId = selectIdIsoCut(type3rdLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Eta()),(double)(*eventLeptons.iso)[idLep[tagZ[2]]],(int)(*eventLeptons.selBits)[idLep[tagZ[2]]],(double)(*eventLeptons.mva)[idLep[tagZ[2]]]);
      }


      //cout<<"hi4"<<endl;


      if(tight3rdLepId == false) continue;

      if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[tagZ[2]]]) == 13) type3l += 0;
      else							       type3l += 2;
      TLorentzVector trilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + 
        		      ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) + 
        		      ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[2])) ) ));
      mass3l = trilep.M();
      //std::cout<<type3l<<std::endl;

      passFilter[ 5] = minMassll > 4; //deltaRllMin > 0.1;
      passFilter[ 6] = ((TLorentzVector*)(*eventMet.p4)[0])->Pt() > 30;
      passFilter[ 7] = minMassZ > minMass && minMassZ < maxMass;
      passFilter[ 8] = ((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt() > 20;
      passFilter[ 9] = mass3l > 100;
      passFilter[10] = true;
      if(applyBtagging) passFilter[10] = bDiscrMax < 0.935;

      if(passFilter[7]==kTRUE && (tagZ[0] == tagZ[1] || tagZ[0] == tagZ[2] || tagZ[1] == tagZ[2])) {printf("ZPROBLEM!\n");assert(0);return;}

      bool passNMinusOne[6] = {                 passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10],
                               passFilter[5] &&                  passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10],
                               passFilter[5] && passFilter[6] &&                  passFilter[8] && passFilter[9] && passFilter[10],
                               passFilter[5] && passFilter[6] && passFilter[7] &&                  passFilter[9] && passFilter[10],
                               passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8]                  && passFilter[10],
			       passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9]};

      bool passAllCuts[1] = {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};
     

      bool controlSel[3] = {passFilter[5] && passFilter[6] && !passFilter[7] && passFilter[9] && bDiscrMax > 0.605 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20,
                            passFilter[5] &&                  !passFilter[7] && passFilter[9] && bDiscrMax > 0.605 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20,
			    passFilter[5] && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() > 20};

      bool passTTZSel[2] = {passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && idJet.size() >= 4 && btagjet[0] > 0.560 && btagjet[1] > 0.560,
                            passFilter[5] && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9] && idJet.size() == 3 && btagjet[0] > 0.800 && btagjet[1] > 0.560};

      

      bool passVBFLoose = passAllCuts[0] && idJet.size() >= 2 &&
                     (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M() > 250;

      bool pass1Jet = passAllCuts[0] && idJet.size() >= 1;



      bool pass2JetControl = passAllCuts[0] && idJet.size() >= 2 &&
                     (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M() < 500 && TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()) < 2.5 &&(( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M() > 100;

      double deltaPhiLeptonMet = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtLN = TMath::Sqrt(2.0*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiLeptonMet)));

      double deltaPhiTriLeptonMet = TMath::Abs(trilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
      double mtEvent = TMath::Sqrt(2.0*trilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiTriLeptonMet)));

     bool passSystCuts[nSystTypes] = {
          passFilter[5] && (double)(*eventMet.ptJESUP)[0]   > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10],
	  passFilter[5] && (double)(*eventMet.ptJESDOWN)[0] > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]
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
     

      double transverseMass = TMath::Sqrt( TMath::Power(pT_0+pT_1+pT_2+met_t,2) - TMath::Power(px_0+px_1+px_2+met_x,2) - TMath::Power(py_0+py_1+py_2+met_y,2));

     TLorentzVector wp4 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + *((TLorentzVector*)(*eventMet.p4)[0]);
     TLorentzVector zp4 = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]]) + *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]) ;
     double transverseMass2 = TMath::Sqrt( TMath::Power(zp4.Et() + wp4.Et(),2) - TMath::Power(zp4.Px() + wp4.Px(),2) - TMath::Power(zp4.Py() + wp4.Py(),2));

     //JES UNCERTAINTY
     TLorentzVector  wp4JESUP = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + *((TLorentzVector*)(*eventMet.p4)[0])*((double)(*eventMet.ptJESUP)[0]/met_t);
     TLorentzVector  wp4JESDOWN = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + *((TLorentzVector*)(*eventMet.p4)[0])*((double)(*eventMet.ptJESDOWN)[0]/met_t);
     //double transverseMass2JESUP = TMath::Sqrt( TMath::Power(zp4.Et() + wp4JESUP.Et(),2) - TMath::Power(zp4.Px() + wp4JESUP.Px(),2) - TMath::Power(zp4.Py() + wp4JESUP.Py(),2));
     //double transverseMass2JESDOWN = TMath::Sqrt( TMath::Power(zp4.Et() + wp4JESDOWN.Et(),2) - TMath::Power(zp4.Px() + wp4JESDOWN.Px(),2) - TMath::Power(zp4.Py() + wp4JESDOWN.Py(),2));
     if(infilecatv[ifile] == 6 && idJetJESUP.size() >= 2 && idJetJESDOWN.size() >= 2){     
       		TLorentzVector metJESup = (*(TLorentzVector*)(*eventMet.p4)[0]);
       		TLorentzVector metJESdown = (*(TLorentzVector*)(*eventMet.p4)[0]);
       		for(int i = 0; i < eventJets.p4->GetEntriesFast() ; i++){  
	 		metJESup = metJESup + ((double)(*eventJets.unc)[i])*(*((TLorentzVector*)(*eventJets.p4)[i]));
	 		metJESdown = metJESdown - ((double)(*eventJets.unc)[i])*(*((TLorentzVector*)(*eventJets.p4)[i]));
       		}   

       wp4JESUP = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + (*((TLorentzVector*)(*eventMet.p4)[0]))*((double)metJESup.Pt()/met_t);
       wp4JESDOWN = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + (*((TLorentzVector*)(*eventMet.p4)[0]))*((double)metJESdown.Pt()/met_t);		       
       //if(transverseMass2JESDOWN != transverseMass2JESDOWN ) cout<<"PROBLEM with JES!"<<endl;     
     }   
    double transverseMass2JESUP = TMath::Sqrt( TMath::Power(zp4.Et() + wp4JESUP.Et(),2) - TMath::Power(zp4.Px() + wp4JESUP.Px(),2) - TMath::Power(zp4.Py() + wp4JESUP.Py(),2));
    double transverseMass2JESDOWN = TMath::Sqrt( TMath::Power(zp4.Et() + wp4JESDOWN.Et(),2) - TMath::Power(zp4.Px() + wp4JESDOWN.Px(),2) - TMath::Power(zp4.Py() + wp4JESDOWN.Py(),2));


     //MET RESOLUTION 10%
     double metRESUP = r.Gaus(met_t,0.1*met_t);
     TLorentzVector  wp4METRESUP = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + *((TLorentzVector*)(*eventMet.p4)[0])*(metRESUP/met_t);
     TLorentzVector  wp4METRESDOWN = *((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]) + *((TLorentzVector*)(*eventMet.p4)[0]);
     double transverseMass2METRESUP = TMath::Sqrt( TMath::Power(zp4.Et() + wp4METRESUP.Et(),2) - TMath::Power(zp4.Px() + wp4METRESUP.Px(),2) - TMath::Power(zp4.Py() + wp4METRESUP.Py(),2));
     double transverseMass2METRESDOWN = TMath::Sqrt( TMath::Power(zp4.Et() + wp4METRESDOWN.Et(),2) - TMath::Power(zp4.Px() + wp4METRESDOWN.Px(),2) - TMath::Power(zp4.Py() + wp4METRESDOWN.Py(),2));

     //LEPTON SCALE UNC 1%
     TLorentzVector wp4LEPScaleUP = (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]))*1.01 + *((TLorentzVector*)(*eventMet.p4)[0]);
     TLorentzVector zp4LEPScaleUP = (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]]))*1.01 + (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]))*1.01 ;  
     TLorentzVector wp4LEPScaleDOWN = (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]]))*0.99 + *((TLorentzVector*)(*eventMet.p4)[0]);
     TLorentzVector zp4LEPScaleDOWN = (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]]))*0.99 + (*((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]]))*0.99 ;
     double transverseMass2LEPScaleUP = TMath::Sqrt( TMath::Power(zp4LEPScaleUP.Et() + wp4LEPScaleUP.Et(),2) - TMath::Power(zp4LEPScaleUP.Px() + wp4LEPScaleUP.Px(),2) - TMath::Power(zp4LEPScaleUP.Py() + wp4LEPScaleUP.Py(),2));
     double transverseMass2LEPScaleDOWN = TMath::Sqrt( TMath::Power(zp4LEPScaleDOWN.Et() + wp4LEPScaleDOWN.Et(),2) - TMath::Power(zp4LEPScaleDOWN.Px() + wp4LEPScaleDOWN.Px(),2) - TMath::Power(zp4LEPScaleDOWN.Py() + wp4LEPScaleDOWN.Py(),2));


      bool passAllCutsMETRESUP = {passFilter[5] &&  metRESUP > 30 && passFilter[7] && passFilter[8] && passFilter[9] && passFilter[10]};    
      bool passVBF = passAllCuts[0] && idJet.size() >= 2 &&
                     (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M() > 500 &&
		     TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()) > 2.5;
      //cout<<"hi4a2"<<endl;
      bool passVBFJESUP = passAllCuts[0] && idJetJESUP.size() >= 2 &&
                     ( ( *(TLorentzVector*)(eventJets.p4->At(idJetJESUP[0])) )*1.05 + ( *(TLorentzVector*)(eventJets.p4->At(idJetJESUP[1])) )*1.05 ).M() > 500 &&
		     TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJetJESUP[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJetJESUP[1]])->Eta()) > 2.5;
      bool passVBFJESDOWN = passAllCuts[0] && idJetJESDOWN.size() >= 2 &&
                     ( ( *(TLorentzVector*)(eventJets.p4->At(idJetJESDOWN[0])) )*0.95 + ( *(TLorentzVector*)(eventJets.p4->At(idJetJESDOWN[1])) )*0.95 ).M() > 500 &&
		     TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJetJESDOWN[1]])->Eta()) > 2.5;
      bool passVBFMETRESUP = passAllCutsMETRESUP && idJet.size() >= 2 &&
                     (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M() > 500 &&
		     TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()) > 2.5;
      bool passVBFMETRESDOWN = passAllCuts[0] && idJet.size() >= 2 &&
                     (( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M() > 500 &&
		     TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()) > 2.5;

      bool passVBFLEPScaleUP = passVBF;
      bool passVBFLEPScaleDOWN = passVBF;

      // luminosity
      double theLumi  = 1.0; if(infilecatv[ifile] != 0) theLumi  = lumi;
      // pile-up
      double puWeight     = 1.0; if(infilecatv[ifile] != 0) puWeight     = nPUScaleFactor(fhDPU    , (double)eventMonteCarlo.puTrueInt);
      double puWeightUp   = 1.0; if(infilecatv[ifile] != 0) puWeightUp   = nPUScaleFactor(fhDPUUp  , (double)eventMonteCarlo.puTrueInt);
      double puWeightDown = 1.0; if(infilecatv[ifile] != 0) puWeightDown = nPUScaleFactor(fhDPUDown, (double)eventMonteCarlo.puTrueInt);
      // lepton efficiency
      double effSF = 1.0;
      if(infilecatv[ifile] != 0 && infilecatv[ifile] != 6){
        for(unsigned int nl=0; nl<idLep.size(); nl++){
	  if(tagZ[2] != (int)nl)
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
		typeLepSel.Data(),fhDMuMediumSF,fhDMuIsoSF,fhDElMediumSF,fhDElTightSF,fhDElMediumMVASF,fhDElTightMVASF);
          else
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
		type3rdLepSel.Data(),fhDMuMediumSF,fhDMuIsoSF,fhDElMediumSF,fhDElTightSF,fhDElMediumMVASF,fhDElTightMVASF);
        }
      }

      // fake rate
      nFakeCount = 0;
      unsigned int typeFakeLepton[2] = {0,0};
      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false && infilecatv[ifile]<6){
        if     ((infilecatv[ifile] == 0 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add Z+jets from data
	  for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    if(tagZ[0] == (int)nl) nFakeCount = nFakeCount + 1;
	    if(tagZ[1] == (int)nl) nFakeCount = nFakeCount + 2;
	    if(tagZ[2] == (int)nl) nFakeCount = nFakeCount + 4;
	    if(tagZ[2] != (int)nl)
	    fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    else
	    fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,type3rdLepSel.Data());
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
     
      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      //double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale;
      if(infilecatv[ifile] == 6) totalWeight=theLumi*puWeight*effSF/float(numberOfEvents);//float(the_input_tree->GetEntries());//*(0.101*3*0.1086); Need to Normalize afterwards
      //if(infilecatv[ifile] == 3 && infilenamev[ifile].Contains("WZJJ_QCD")) totalWeight=totalWeight * 0.565*(6)*1.06; // reweight according to data 15% unc. QCD+EWK
      if(infilecatv[ifile] == 3 && infilenamev[ifile].Contains("WZJJ_EWK")) totalWeight=totalWeight; // reweight according to data 15% unc. QCD+EWK	
      //if(infilecatv[ifile] == 3) totalWeight=totalWeight * 0.73; // reweight according to data 15% unc. QCD
      if(totalWeight == 0) continue;
      if(theCategory == 0) totalWeight = 1.0;
      // end event weighting
      if((infilecatv[ifile] != 0 || theCategory == 0) && passAllCuts[0]) sumEventsProcess[ifile] += totalWeight;

      if(passAllCuts[0] && infilecatv[ifile] == 0) totalFakeDataCount[type3l][nFakeCount]++;


      
 
      for(int thePlot=0; thePlot<allPlots; thePlot++){
	double theVar = 0.0;
	bool makePlot = false;
 	if     (thePlot ==  0 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(mtLN,199.999);}
	else if(thePlot ==  1 && passNMinusOne[0])                       {makePlot = true;theVar = TMath::Min(minMassll,199.999);}
	else if(thePlot ==  2 && pass2JetControl)                       {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot ==  3 && passNMinusOne[2])                       {makePlot = true;theVar = TMath::Max(TMath::Min(minMassZ,149.999),50.001);}
	else if(thePlot ==  4 && passNMinusOne[3])                       {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt(),199.999);}
	else if(thePlot ==  5 && passNMinusOne[4])                       {makePlot = true;theVar = TMath::Min(mass3l,399.999);}
	else if(thePlot ==  6 && passNMinusOne[5])                       {makePlot = true;theVar = TMath::Min(bDiscrMax,0.999);}
	else if(thePlot ==  7 && passNMinusOne[0])                       {makePlot = true;theVar = TMath::Min(deltaRllMin,3.999);}
	else if(thePlot ==  8 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(minPMET,199.999);}
	else if(thePlot ==  9 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[0]]])->Pt(),199.999);}
	else if(thePlot == 10 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[1]]])->Pt(),199.999);}
	else if(thePlot == 11 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(minMET,199.999);}
	else if(thePlot == 12 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 13 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min((double)eventEvent.rho,39.999);}
	else if(thePlot == 14 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	else if(thePlot == 15 && passAllCuts[0])                         {makePlot = true;theVar = (double)type3l;}
	else if(thePlot == 16 && passAllCuts[0])                         {makePlot = true;theVar = dPhiJetMET*180/TMath::Pi();}
	else if(thePlot == 17 && passAllCuts[0])                         {makePlot = true;theVar = dPhiLepMETMin*180/TMath::Pi();}
	else if(thePlot == 18 && controlSel[0])                          {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 19 && controlSel[1])                          {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
	else if(thePlot == 20 && controlSel[2])                          {makePlot = true;theVar = TMath::Max(TMath::Min(mass3l,249.999),50.001);}
	else if(thePlot == 21 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(mtEvent,999.999);}
	else if(thePlot == 22 && passVBF)                                {makePlot = true;theVar = TMath::Min(transverseMass2,999.999);}
	else if(thePlot == 23 && passTTZSel[0])                          {makePlot = true;theVar = TMath::Min(theHT,999.999);}
	else if(thePlot == 24 && passTTZSel[1])                          {makePlot = true;theVar = TMath::Min(theHT,999.999);}
	else if(thePlot == 25 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min((double)(numberQuarks[0]+10*numberQuarks[1]),49.499);}
	else if(thePlot == 26 && passNMinusOne[3])	                 {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[tagZ[2]]])->Pt(),199.999);}
	else if(thePlot == 27 && passAllCuts[0])                         {makePlot = true;theVar = TMath::Min(mtLN,199.999);}
	else if(thePlot == 28 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 11 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 29 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 13 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) <  1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 30 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 11 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 31 && controlSel[2] && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[2]]) == 13 && TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Eta()) >= 1.479) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt(),199.999);}
	else if(thePlot == 32 && controlSel[2])                          {makePlot = true;theVar = TMath::Max(TMath::Min(mass3l,249.999),50.001);}
	else if(thePlot == 33 && passNMinusOne[5] && numberQuarks[1] == 0                  ) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 34 && passNMinusOne[5] && numberQuarks[1] == 0 && passFilter[10]) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 35 && passNMinusOne[5] && numberQuarks[1]  > 0                  ) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 36 && passNMinusOne[5] && numberQuarks[1]  > 0 && passFilter[10]) {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
	else if(thePlot == 37 && pass2JetControl)                                {makePlot = true;theVar = TMath::Min(transverseMass2,1199.999);}
	else if(thePlot == 38 && passVBFLoose)                                {makePlot = true;theVar = TMath::Min((( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M(),1999.999);}
	else if(thePlot == 39 && passVBFLoose)                                {makePlot = true;theVar = TMath::Min(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()),7.999);}
        //else if(thePlot == 40 && pass2JetControl)                         {makePlot = true;theVar = (double)type3l;}
        else if(thePlot == 40 && passVBF)                         {makePlot = true;theVar = (double)type3l;}
	else if(thePlot == 41 && pass2JetControl)                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[0])) )->Pt();}
        else if(thePlot == 42 && pass2JetControl)                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[1])) )->Pt();}
        else if(thePlot == 43 && pass2JetControl)                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[0])) )->Eta();}
        else if(thePlot == 44 && pass2JetControl)                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[1])) )->Eta();}
        else if(thePlot == 45 && pass2JetControl)                         {makePlot = true;theVar = (double) TMath::Min((( ( *(TLorentzVector*)(eventJets.p4->At(idJet[0])) ) + ( *(TLorentzVector*)(eventJets.p4->At(idJet[1])) ) )).M(),499.999);}
        else if(thePlot == 46 && pass2JetControl)                         {makePlot = true;theVar = (double) TMath::Min(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()),2.4999);}
	else if(thePlot == 47 && pass2JetControl )                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[0])) )->Pt();}
	else if(thePlot == 48 && pass2JetControl )                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[1])) )->Pt();}
	else if(thePlot == 49 && pass2JetControl )                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[2])) )->Pt();}
        else if(thePlot == 50 && pass2JetControl )                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[0])) )->Eta();}
	else if(thePlot == 51 && pass2JetControl )                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[1])) )->Eta();}
	else if(thePlot == 52 && pass2JetControl )                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[2])) )->Eta();}
        else if(thePlot == 53 && pass2JetControl)                         {makePlot = true;theVar = (double)type3l;}
	else if(thePlot == 54 && passVBF)                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventLeptons.p4->At(idLep[0])) )->Pt();}
	else if(thePlot == 55 && pass1Jet)                         {makePlot = true;theVar = (double)( (TLorentzVector*)(eventJets.p4->At(idJet[0])) )->Pt();}
	else if(thePlot == 56 && passAllCuts[0])                         {makePlot = true;theVar = (double)transverseMass2;}

	if( (makePlot && theCategory<6) || (makePlot && theCategory == 6 && infilenamev[ifile].Contains("M700")) ) histo[     4][thePlot][theCategory]->Fill(theVar,totalWeight);
	if( (makePlot && theCategory<6) || (makePlot && theCategory == 6 && infilenamev[ifile].Contains("M700")) ) histo[type3l][thePlot][theCategory]->Fill(theVar,totalWeight);
      }
      //if(infilenamev[ifile].Contains("M700")){
      if(1 && passVBF ) {

	double MVAVar = TMath::Min(transverseMass2, 2199.99);
	double MVAVarJESUP = TMath::Min(transverseMass2JESUP, 2199.99);
	double MVAVarJESDOWN = TMath::Min(transverseMass2JESDOWN, 2199.99);
	double MVAVarMETRESUP = TMath::Min(transverseMass2METRESUP, 2199.99);
	double MVAVarMETRESDOWN = TMath::Min(transverseMass2METRESDOWN, 2199.99);
	double MVAVarLEPScaleUP = TMath::Min(transverseMass2LEPScaleUP, 2199.99);
	double MVAVarLEPScaleDOWN = TMath::Min(transverseMass2LEPScaleDOWN, 2199.99);
        

        if     (theCategory == 0){
	  if(passAllCuts[SIGSEL]) histo_Data->Fill(MVAVar,totalWeight);
	 
	}
        else if(theCategory == 1){
	  if(passAllCuts[SIGSEL]) {
	    if(typeFakeLepton[0]+typeFakeLepton[1] != goodIsTight == idTight.size()) {printf("PROBLEMFake %d %d %d %d\n",typeFakeLepton[0],typeFakeLepton[1],goodIsTight,(int)idTight.size()); return;}

	    if(MVAVar>50.) cout<<"fake weight"<<totalWeight<<endl;
	   
	   
	     histo_Fake->Fill(MVAVar,totalWeight);
	  }
        }
        else if(theCategory == 2 ){
	  if(passAllCuts[SIGSEL]) {
	  
	     histo_Zg->Fill(MVAVar,totalWeight);

	     histo_Zg_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_Zg_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_Zg_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_Zg_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_Zg_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_Zg_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_Zg_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_Zg_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Zg_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Zg_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Zg_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Zg_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_Zg_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_Zg_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_Zg_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
          if(passSystCuts[METUP])  histo_Zg_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_Zg_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
	  if(passVBFJESUP) histo_Zg_CMS_scale_jUp->Fill(MVAVarJESUP,totalWeight);
	  if(passVBFJESDOWN) histo_Zg_CMS_scale_jDown->Fill(MVAVarJESDOWN,totalWeight);
	  if(passVBFMETRESUP) histo_Zg_CMS_res_metUp->Fill(MVAVarMETRESUP,totalWeight);
	  if(passVBFMETRESDOWN) histo_Zg_CMS_res_metDown->Fill(MVAVarMETRESDOWN,totalWeight);
	  if(passVBFLEPScaleUP) histo_Zg_CMS_scale_lUp->Fill(MVAVarLEPScaleUP,totalWeight);
	  if(passVBFLEPScaleDOWN) histo_Zg_CMS_scale_lDown->Fill(MVAVarLEPScaleDOWN,totalWeight);
	}
        else if(theCategory == 3){
	  if(/*passAllCuts[SIGSEL]*/ 1) {


	     histo_WZ->Fill(MVAVar,totalWeight);

	     if(infilenamev[ifile].Contains("powheg")){
               histo_WZ_CMS_QCDunc_scaleDown->Fill(MVAVar,totalWeight*0.88);
               histo_WZ_CMS_QCDunc_scaleUp->Fill(MVAVar,totalWeight*1.12);}
             else{
               histo_WZ_CMS_QCDunc_scaleDown->Fill(MVAVar,totalWeight);
               histo_WZ_CMS_QCDunc_scaleUp->Fill(MVAVar,totalWeight);}

             if(infilenamev[ifile].Contains("EWK")){
               histo_WZ_CMS_ewkscaleDown->Fill(MVAVar,totalWeight*0.7);
               histo_WZ_CMS_ewkscaleUp->Fill(MVAVar,totalWeight*1.3);}
             else{
               histo_WZ_CMS_ewkscaleDown->Fill(MVAVar,totalWeight);
               histo_WZ_CMS_ewkscaleUp->Fill(MVAVar,totalWeight);}


	     WZ_2d->Fill( (double)((TLorentzVector*)(eventLeptons.p4->At(idLep[0])) )->Pt(),transverseMass2,totalWeight);
	     histo_WZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_WZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_WZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_WZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_WZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_WZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
             else
	     for(int npdf=0; npdf<102; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_WZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_WZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_WZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_WZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_WZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_WZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
	  }
	  if(passVBFJESDOWN) histo_WZ_CMS_scale_jDown->Fill(MVAVarJESDOWN,totalWeight);
	  if(passVBFJESUP) histo_WZ_CMS_scale_jUp->Fill(MVAVarJESUP,totalWeight);
	  if(passVBFMETRESDOWN) histo_WZ_CMS_res_metDown->Fill(MVAVarMETRESDOWN,totalWeight);
	  if(passVBFMETRESUP) histo_WZ_CMS_res_metUp->Fill(MVAVarMETRESUP,totalWeight);
	  if(passVBFLEPScaleUP) histo_WZ_CMS_scale_lUp->Fill(MVAVarLEPScaleUP,totalWeight);
	  if(passVBFLEPScaleDOWN) histo_WZ_CMS_scale_lDown->Fill(MVAVarLEPScaleDOWN,totalWeight);
          if(passSystCuts[METUP])  histo_WZ_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_WZ_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
	}
        else if(theCategory == 4){
	  if(passAllCuts[SIGSEL]) {
	     histo_ZZ->Fill(MVAVar,totalWeight);
	     histo_ZZ_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_ZZ_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_ZZ_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_ZZ_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_ZZ_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_ZZ_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_ZZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ZZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
          }
	  if(passVBFJESDOWN) histo_ZZ_CMS_scale_jDown->Fill(MVAVarJESDOWN,totalWeight);
          if(passVBFJESUP) histo_ZZ_CMS_scale_jUp->Fill(MVAVarJESUP,totalWeight);
	  if(passVBFMETRESDOWN) histo_ZZ_CMS_res_metDown->Fill(MVAVarMETRESDOWN,totalWeight);
          if(passVBFMETRESUP) histo_ZZ_CMS_res_metUp->Fill(MVAVarMETRESUP,totalWeight);
          if(passVBFLEPScaleUP) histo_ZZ_CMS_scale_lUp->Fill(MVAVarLEPScaleUP,totalWeight);
	  if(passVBFLEPScaleDOWN) histo_ZZ_CMS_scale_lDown->Fill(MVAVarLEPScaleDOWN,totalWeight);
          if(passSystCuts[METUP])  histo_ZZ_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_ZZ_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 5){
	  if(passAllCuts[SIGSEL]) {
	     histo_VVV->Fill(MVAVar,totalWeight);
	     histo_VVV_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_VVV_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_VVV_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_VVV_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_VVV_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_VVV_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_VVV_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_VVV_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
          }
	  if(passVBFJESDOWN) histo_VVV_CMS_scale_jDown->Fill(MVAVarJESDOWN,totalWeight);
	  if(passVBFJESUP) histo_VVV_CMS_scale_jUp->Fill(MVAVarJESUP,totalWeight);
	  if(passVBFMETRESDOWN) histo_VVV_CMS_res_metDown->Fill(MVAVarMETRESDOWN,totalWeight);
	  if(passVBFMETRESUP) histo_VVV_CMS_res_metUp->Fill(MVAVarMETRESUP,totalWeight);
	  if(passVBFLEPScaleUP) histo_VVV_CMS_scale_lUp->Fill(MVAVarLEPScaleUP,totalWeight);
	  if(passVBFLEPScaleDOWN) histo_VVV_CMS_scale_lDown->Fill(MVAVarLEPScaleDOWN,totalWeight);
          if(passSystCuts[METUP])  histo_VVV_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_VVV_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
        }
       else if( theCategory == 6 ){
	  if(passAllCuts[SIGSEL]) {
	    if(infilenamev[ifile].Contains("SignalM700")) Higgs_2d->Fill( (double)((TLorentzVector*)(eventLeptons.p4->At(idLep[0])) )->Pt(), transverseMass2,totalWeight);
	    for(int i=0; i<11; i++) { 
	      if(infilenamev[ifile].Contains("Signal"+SignalSuffix[i])){
		histo_Higgs_[i]->Fill(MVAVar,totalWeight);
	      	}
	      }

	    //cout<<"trans after  :"<<transverseMass2<<endl;

	     histo_Higgs_CMS_QCDScaleBounding[0]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2));
	     histo_Higgs_CMS_QCDScaleBounding[1]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5));
	     histo_Higgs_CMS_QCDScaleBounding[2]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1));
	     histo_Higgs_CMS_QCDScaleBounding[3]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2));
	     histo_Higgs_CMS_QCDScaleBounding[4]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1));
	     histo_Higgs_CMS_QCDScaleBounding[5]  ->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5));
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<102; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]));
             else if(infilenamev[ifile].Contains("powheg") == true)
	     for(int npdf=0; npdf<102; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+0]));
	     else
	     for(int npdf=0; npdf<102; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_Higgs_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Higgs_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Higgs_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Higgs_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Higgs_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_Higgs_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_Higgs_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_Higgs_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
          }
          if(passSystCuts[METUP])  histo_Higgs_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_Higgs_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
	  if(passVBFJESDOWN) {
	    //cout<<"transDOWN after:"<<transverseMass2JESDOWN<<endl;
	    for(int i=0; i<11; i++) { 
	      if(infilenamev[ifile].Contains("Signal"+SignalSuffix[i]+"_")){
		histo_Higgs_CMS_scale_jDown_[i]->Fill(MVAVarJESDOWN,totalWeight);
	      }
	    }
	  }
	  if(passVBFJESUP) {
	    //cout<<"transUP after:"<<transverseMass2JESUP<<endl;
            for(int i=0; i<11; i++) {
              if(infilenamev[ifile].Contains("Signal"+SignalSuffix[i]+"_")){
                histo_Higgs_CMS_scale_jUp_[i]->Fill(MVAVarJESUP,totalWeight);
              }
            }
          }
         if(passVBFMETRESDOWN) {
	    for(int i=0; i<11; i++) { 
	      if(infilenamev[ifile].Contains("Signal"+SignalSuffix[i]+"_")){
		histo_Higgs_CMS_res_metDown_[i]->Fill(MVAVarMETRESDOWN,totalWeight);
	      }
	    }
	  }
	  if(passVBFMETRESUP) {
            for(int i=0; i<11; i++) {
              if(infilenamev[ifile].Contains("Signal"+SignalSuffix[i]+"_")){
                histo_Higgs_CMS_res_metUp_[i]->Fill(MVAVarMETRESUP,totalWeight);
              }
            }
          }

        if(passVBFLEPScaleDOWN) {
	    for(int i=0; i<11; i++) { 
	      if(infilenamev[ifile].Contains("Signal"+SignalSuffix[i]+"_")){
		histo_Higgs_CMS_scale_lDown_[i]->Fill(MVAVarLEPScaleDOWN,totalWeight);
	      }
	    }
	  }
	  if(passVBFLEPScaleUP) {
            for(int i=0; i<11; i++) {
              if(infilenamev[ifile].Contains("Signal"+SignalSuffix[i]+"_")){
                histo_Higgs_CMS_scale_lUp_[i]->Fill(MVAVarLEPScaleUP,totalWeight);
              }
            }
          }

	}
	else {
	  printf("CATEGORY PROBLEM!\n"); return;
	}
      } // making data cards
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
    for(int j=0; j<allStates; j++){
      char output[200];
      sprintf(output,"histowz_nice_%d_%d.root",j,thePlot);	  
      TFile* outFilePlotsNote = new TFile(output,"recreate");
      outFilePlotsNote->cd();
      for(int np=0; np<histBins; np++) histo[j][thePlot][np]->Write();
      outFilePlotsNote->Close();
    }
  }
  
  printf("QCD Init: WZ(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_WZ_CMS_QCDScaleInitial[6],histo_WZ_CMS_QCDScaleInitial[0],histo_WZ_CMS_QCDScaleInitial[1],histo_WZ_CMS_QCDScaleInitial[2],histo_WZ_CMS_QCDScaleInitial[3],histo_WZ_CMS_QCDScaleInitial[4],histo_WZ_CMS_QCDScaleInitial[5]);
  for(int nj=0; nj<6; nj++) histo_WZ_CMS_QCDScaleInitial[nj] = histo_WZ_CMS_QCDScaleInitial[nj] / histo_WZ_CMS_QCDScaleInitial[6];
  printf("QCD RateInit: WZ(%f/%f/%f/%f/%f/%f)\n",
    histo_WZ_CMS_QCDScaleInitial[0],histo_WZ_CMS_QCDScaleInitial[1],histo_WZ_CMS_QCDScaleInitial[2],histo_WZ_CMS_QCDScaleInitial[3],histo_WZ_CMS_QCDScaleInitial[4],histo_WZ_CMS_QCDScaleInitial[5]);

  // correcting by initial normalization
  for(int nj=0; nj<6; nj++) histo_WZ_CMS_QCDScaleBounding[nj]->Scale(1.0/histo_WZ_CMS_QCDScaleInitial[nj]);

  printf("QCD Corr: WZ(%f:%f/%f/%f/%f/%f/%f) ZZ(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) ZH(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_WZ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ZZ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_Zg->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_Zg_CMS_QCDScaleBounding[5]->GetSumOfWeights());

  for(int i=1; i<=histo_Zg->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_Zg_CMS_MVAZHStatBoundingUp        ->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorUp  *histo_Zg->GetBinError(i),0.000001));
    histo_Zg_CMS_MVAZHStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorDown*histo_Zg->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp      ->SetBinContent(i,TMath::Max(histo_VVV    ->GetBinContent(i)+factorUp  *histo_VVV    ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_VVV    ->GetBinContent(i)+factorDown*histo_VVV    ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_WZ     ->GetBinContent(i)+factorUp  *histo_WZ     ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_WZ     ->GetBinContent(i)+factorDown*histo_WZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_ZZ     ->GetBinContent(i)+factorUp  *histo_ZZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingDown      ->SetBinContent(i,TMath::Max(histo_ZZ     ->GetBinContent(i)+factorDown*histo_ZZ     ->GetBinError(i),0.000001));
    histo_Fake_CMS_MVAFakeMStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_Fake  ->GetBinContent(i)+factorUp  *histo_Fake  ->GetBinError(i),0.000001));
    histo_Fake_CMS_MVAFakeMStatBoundingDown->SetBinContent(i,TMath::Max(histo_Fake  ->GetBinContent(i)+factorDown*histo_Fake  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorUp  *histo_FakeE  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingDown->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorDown*histo_FakeE  ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_Higgs  ->GetBinContent(i)+factorUp  *histo_Higgs  ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingDown->SetBinContent(i,TMath::Max(histo_Higgs  ->GetBinContent(i)+factorDown*histo_Higgs  ->GetBinError(i),0.000001));

    histo_Zg_CMS_MVAZHStatBoundingBinUp[i-1]        ->Add(histo_Zg     ); histo_Zg_CMS_MVAZHStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorUp  *histo_Zg->GetBinError(i),0.000001));
    histo_Zg_CMS_MVAZHStatBoundingBinDown[i-1]      ->Add(histo_Zg     ); histo_Zg_CMS_MVAZHStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_Zg->GetBinContent(i)+factorDown*histo_Zg->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	    ->Add(histo_VVV    ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]     ->SetBinContent(i,TMath::Max(histo_VVV	->GetBinContent(i)+factorUp  *histo_VVV    ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]    ->Add(histo_VVV    ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]   ->SetBinContent(i,TMath::Max(histo_VVV	->GetBinContent(i)+factorDown*histo_VVV    ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]	    ->Add(histo_WZ     ); histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]       ->SetBinContent(i,TMath::Max(histo_WZ	->GetBinContent(i)+factorUp  *histo_WZ     ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]	    ->Add(histo_WZ     ); histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]     ->SetBinContent(i,TMath::Max(histo_WZ	->GetBinContent(i)+factorDown*histo_WZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]	    ->Add(histo_ZZ     ); histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]       ->SetBinContent(i,TMath::Max(histo_ZZ	->GetBinContent(i)+factorUp  *histo_ZZ     ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]	    ->Add(histo_ZZ     ); histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]     ->SetBinContent(i,TMath::Max(histo_ZZ	->GetBinContent(i)+factorDown*histo_ZZ     ->GetBinError(i),0.000001));
    histo_Fake_CMS_MVAFakeMStatBoundingBinUp[i-1]  ->Add(histo_Fake  ); histo_Fake_CMS_MVAFakeMStatBoundingBinUp[i-1]    ->SetBinContent(i,TMath::Max(histo_Fake  ->GetBinContent(i)+factorUp  *histo_Fake  ->GetBinError(i),0.000001));
    histo_Fake_CMS_MVAFakeMStatBoundingBinDown[i-1]->Add(histo_Fake  ); histo_Fake_CMS_MVAFakeMStatBoundingBinDown[i-1]  ->SetBinContent(i,TMath::Max(histo_Fake  ->GetBinContent(i)+factorDown*histo_Fake  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[i-1]  ->Add(histo_FakeE  ); histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[i-1]    ->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorUp  *histo_FakeE  ->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[i-1]->Add(histo_FakeE  ); histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[i-1]  ->SetBinContent(i,TMath::Max(histo_FakeE  ->GetBinContent(i)+factorDown*histo_FakeE  ->GetBinError(i),0.000001));

   histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[i-1]  ->Add(histo_Higgs  ); histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[i-1]    ->SetBinContent(i,TMath::Max(histo_Higgs  ->GetBinContent(i)+factorUp  *histo_Higgs  ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[i-1]->Add(histo_Higgs  ); histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[i-1]  ->SetBinContent(i,TMath::Max(histo_Higgs  ->GetBinContent(i)+factorDown*histo_Higgs  ->GetBinError(i),0.000001));
  }


  char outputLimits[200];
  sprintf(outputLimits,"wz3l%2s_input_%4s.root",finalStateName,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data   ->Write();
  for (int i=1;i<10;i++){
    if(histo_Zg->GetBinContent(i)<=0.){histo_Zg->SetBinContent(i,0.000001);histo_Zg->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_scale_jUp->GetBinContent(i)<=0.){histo_Zg_CMS_scale_jUp->SetBinContent(i,0.000001);histo_Zg_CMS_scale_jUp->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_scale_jDown->GetBinContent(i)<=0.){histo_Zg_CMS_scale_jDown->SetBinContent(i,0.000001);histo_Zg_CMS_scale_jDown->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_res_metUp->GetBinContent(i)<=0.){histo_Zg_CMS_res_metUp->SetBinContent(i,0.000001);histo_Zg_CMS_res_metUp->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_res_metDown->GetBinContent(i)<=0.){histo_Zg_CMS_res_metDown->SetBinContent(i,0.000001);histo_Zg_CMS_res_metDown->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_scale_lUp->GetBinContent(i)<=0.){histo_Zg_CMS_scale_lUp->SetBinContent(i,0.000001);histo_Zg_CMS_scale_lUp->SetBinError(i,0.2*1.6);}
    if(histo_Zg_CMS_scale_lDown->GetBinContent(i)<=0.){histo_Zg_CMS_scale_lDown->SetBinContent(i,0.000001); histo_Zg_CMS_scale_lDown->SetBinError(i,0.2*1.6);}
    if(histo_ZZ->GetBinContent(i)<=0.){histo_ZZ->SetBinContent(i,0.000001);histo_ZZ->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_scale_jUp->GetBinContent(i)<=0.){histo_ZZ_CMS_scale_jUp->SetBinContent(i,0.000001);histo_ZZ_CMS_scale_jUp->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_scale_jDown->GetBinContent(i)<=0.){histo_ZZ_CMS_scale_jDown->SetBinContent(i,0.000001);histo_ZZ_CMS_scale_jDown->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_res_metUp->GetBinContent(i)<=0.){histo_ZZ_CMS_res_metUp->SetBinContent(i,0.000001);histo_ZZ_CMS_res_metUp->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_res_metDown->GetBinContent(i)<=0.){histo_ZZ_CMS_res_metDown->SetBinContent(i,0.000001);histo_ZZ_CMS_res_metDown->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_scale_lUp->GetBinContent(i)<=0.){histo_ZZ_CMS_scale_lUp->SetBinContent(i,0.000001);histo_ZZ_CMS_scale_lUp->SetBinError(i,0.2*1.6);}
    if(histo_ZZ_CMS_scale_lDown->GetBinContent(i)<=0.){histo_ZZ_CMS_scale_lDown->SetBinContent(i,0.000001);histo_ZZ_CMS_scale_lDown->SetBinError(i,0.2*1.6);}
    if(histo_WZ->GetBinContent(i)<=0.){histo_WZ->SetBinContent(i,0.000001);histo_WZ->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_scale_jUp->GetBinContent(i)<=0.){histo_WZ_CMS_scale_jUp->SetBinContent(i,0.000001);histo_WZ_CMS_scale_jUp->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_scale_jDown->GetBinContent(i)<=0.){histo_WZ_CMS_scale_jDown->SetBinContent(i,0.000001);histo_WZ_CMS_scale_jDown->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_res_metUp->GetBinContent(i)<=0.){histo_WZ_CMS_res_metUp->SetBinContent(i,0.000001);histo_WZ_CMS_res_metUp->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_res_metDown->GetBinContent(i)<=0.){histo_WZ_CMS_res_metDown->SetBinContent(i,0.000001);histo_WZ_CMS_res_metDown->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_scale_lUp->GetBinContent(i)<=0.){histo_WZ_CMS_scale_lUp->SetBinContent(i,0.000001);histo_WZ_CMS_scale_lUp->SetBinError(i,0.2*1.6);}
    if(histo_WZ_CMS_scale_lDown->GetBinContent(i)<=0.){histo_WZ_CMS_scale_lDown->SetBinContent(i,0.000001);histo_WZ_CMS_scale_lDown->SetBinError(i,0.2*1.6);}
    if(histo_VVV->GetBinContent(i)<=0.){histo_VVV->SetBinContent(i,0.000001);histo_VVV->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_scale_jUp->GetBinContent(i)<=0.){histo_VVV_CMS_scale_jUp->SetBinContent(i,0.000001);histo_VVV_CMS_scale_jUp->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_scale_jDown->GetBinContent(i)<=0.){histo_VVV_CMS_scale_jDown->SetBinContent(i,0.000001);histo_VVV_CMS_scale_jDown->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_res_metUp->GetBinContent(i)<=0.){histo_VVV_CMS_res_metUp->SetBinContent(i,0.000001);histo_VVV_CMS_res_metUp->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_res_metDown->GetBinContent(i)<=0.){histo_VVV_CMS_res_metDown->SetBinContent(i,0.000001);histo_VVV_CMS_res_metDown->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_scale_lUp->GetBinContent(i)<=0.){histo_VVV_CMS_scale_lUp->SetBinContent(i,0.000001);histo_VVV_CMS_scale_lUp->SetBinError(i,0.2*1.6);}
    if(histo_VVV_CMS_scale_lDown->GetBinContent(i)<=0.){histo_VVV_CMS_scale_lDown->SetBinContent(i,0.000001);histo_VVV_CMS_scale_lDown->SetBinError(i,0.2*1.6);}
	}

  histo_Zg     ->Write();
  histo_Zg_CMS_scale_jUp     ->Write();
  histo_Zg_CMS_scale_jDown     ->Write();
  histo_Zg_CMS_res_metUp     ->Write();
  histo_Zg_CMS_res_metDown     ->Write();
  histo_Zg_CMS_scale_lUp     ->Write();
  histo_Zg_CMS_scale_lDown     ->Write();

  histo_VVV    ->Write();
  histo_VVV_CMS_scale_jUp     ->Write();
  histo_VVV_CMS_scale_jDown     ->Write();
  histo_VVV_CMS_res_metUp     ->Write();
  histo_VVV_CMS_res_metDown     ->Write();
  histo_VVV_CMS_scale_lUp     ->Write();
  histo_VVV_CMS_scale_lDown     ->Write();


  histo_WZ     ->Write();
  histo_WZ_CMS_scale_jUp     ->Write();
  histo_WZ_CMS_scale_jDown     ->Write();
  histo_WZ_CMS_res_metUp     ->Write();
  histo_WZ_CMS_res_metDown     ->Write();
  histo_WZ_CMS_scale_lUp     ->Write();
  histo_WZ_CMS_scale_lDown     ->Write();


  histo_WZ_CMS_QCDunc_scaleDown->Write();
  histo_WZ_CMS_QCDunc_scaleUp->Write();

  histo_WZ_CMS_ewkscaleDown->Write();
  histo_WZ_CMS_ewkscaleUp->Write();

  histo_ZZ     ->Write();
  histo_ZZ_CMS_scale_jUp     ->Write();
  histo_ZZ_CMS_scale_jDown     ->Write();
  histo_ZZ_CMS_res_metUp     ->Write();
  histo_ZZ_CMS_res_metDown     ->Write();
  histo_ZZ_CMS_scale_lUp     ->Write();
  histo_ZZ_CMS_scale_lDown     ->Write();

  cout<<"before: "<<fhDfakeShape->Integral()<<endl;
  //Get Fake bin error from side band
  cout<<"before error 6: "<<fhDfakeShape->GetBinError(6)<<endl;
   cout<<"before content 6: "<<fhDfakeShape->GetBinContent(6)<<endl;
  fhDfakeShape->Scale(histo_Fake->Integral()/fhDfakeShape->Integral());
  cout<<"after error 6: "<<fhDfakeShape->GetBinError(6)<<endl;
  cout<<"after content 6: "<<fhDfakeShape->GetBinContent(6)<<endl;
    for (int i=1;i<8;i++){
    if(histo_Fake->GetBinContent(i)<=0.){
      histo_Fake->SetBinContent(i,0.000001);
      histo_Fake->SetBinError(i,0.2*1.6);
      //if(i>3){
      //histo_Fake->SetBinContent(i,fhDfakeShape->GetBinContent(i));
      //histo_Fake->SetBinError(i,fhDfakeShape->GetBinError(i));
      //}

    }
    }
    cout<<"after: "<<fhDfakeShape->Integral()<<endl;
  histo_Fake  ->Write();
  histo_FakeE  ->Write();

  for(int i=0; i<11; i++) {
	histo_Higgs_[i]  ->Write();
	histo_Higgs_CMS_scale_jUp_[i]  ->Write();
	histo_Higgs_CMS_scale_jDown_[i]  ->Write();
	histo_Higgs_CMS_res_metUp_[i]  ->Write();
	histo_Higgs_CMS_res_metDown_[i]  ->Write();
	histo_Higgs_CMS_scale_lUp_[i]  ->Write();
	histo_Higgs_CMS_scale_lDown_[i]  ->Write();
	}


 
  // MuonPt->Write();
  //ElePt->Write();
  //WZ_2d->Write();
  //Higgs_2d->Write();

  cout << histo_Data   ->GetSumOfWeights() << " ";
  cout << histo_Zg     ->GetSumOfWeights() << " ";
  cout << histo_VVV    ->GetSumOfWeights() << " ";
  cout << histo_WZ     ->GetSumOfWeights() << " ";
  cout << histo_ZZ     ->GetSumOfWeights() << " ";
  cout << histo_Fake  ->GetSumOfWeights() << " ";
  cout << histo_FakeE  ->GetSumOfWeights() << " ";
  cout << histo_Higgs_[1]  ->GetSumOfWeights() << " ";
  cout << endl;

  /*  printf("uncertainties Stat\n");
  histo_Zg_CMS_MVAZHStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg	->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAZHStatBoundingUp   ->GetBinContent(i)/histo_Zg   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVAZHStatBoundingDown      ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg	->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAZHStatBoundingDown ->GetBinContent(i)/histo_Zg	->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingDown    ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingUp	     ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAWZStatBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ	->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingDown      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingUp	     ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAZZStatBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ	->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingDown      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Fake_CMS_MVAFakeMStatBoundingUp  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Fake->GetBinContent(i)>0)printf("%5.1f ",histo_Fake_CMS_MVAFakeMStatBoundingUp	     ->GetBinContent(i)/histo_Fake   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Fake_CMS_MVAFakeMStatBoundingDown->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Fake->GetBinContent(i)>0)printf("%5.1f ",histo_Fake_CMS_MVAFakeMStatBoundingDown      ->GetBinContent(i)/histo_Fake   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_FakeE->GetBinContent(i)>0)printf("%5.1f ",histo_FakeE_CMS_MVAFakeEStatBoundingUp	     ->GetBinContent(i)/histo_FakeE   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_FakeE_CMS_MVAFakeEStatBoundingDown->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_FakeE->GetBinContent(i)>0)printf("%5.1f ",histo_FakeE_CMS_MVAFakeEStatBoundingDown      ->GetBinContent(i)/histo_FakeE   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffM\n");
  histo_Zg_CMS_MVALepEffMBoundingUp       ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVALepEffMBoundingDown     ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffE\n");
  histo_Zg_CMS_MVALepEffEBoundingUp       ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVALepEffEBoundingDown     ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingUp       ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingDown     ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg     ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties MET\n");
  histo_Zg_CMS_MVAMETBoundingUp           ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_Zg      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_MVAMETBoundingDown         ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_MVAMETBoundingDown   ->GetBinContent(i)/histo_Zg   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties PU\n");
  histo_Zg_CMS_PUBoundingUp               ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_PUBoundingUp      ->GetBinContent(i)/histo_Zg      ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Zg_CMS_PUBoundingDown             ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_Zg        ->GetBinContent(i)>0)printf("%5.1f ",histo_Zg_CMS_PUBoundingDown   ->GetBinContent(i)/histo_Zg   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingUp    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingDown  ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WZ_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_WZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingUp   ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ZZ_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_Zg->GetNbinsX(); i++) {if(histo_ZZ    ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingDown ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  */  

outFileLimits->Close();

//double lumiE = 1.027;
  double lumiE = 1.023;
  double lumiEstab = 1.015;

  double systLepResE[4] = {1.01,1.01,1.01,1.01};
  double systLepResM[4] = {1.01,1.01,1.01,1.01};

  for(int nb=1; nb<=nBinMVA; nb++){
     // QCD study
    double systQCDScale[4] = {TMath::Abs(histo_Zg_CMS_QCDScaleBounding[0] ->GetBinContent(nb)-histo_Zg ->GetBinContent(nb)),
                              TMath::Abs(histo_VVV_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_VVV->GetBinContent(nb)),
                              TMath::Abs(histo_WZ_CMS_QCDScaleBounding[0] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb)),
                              TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[0] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb))};
    for(int nqcd=1; nqcd<6; nqcd++) {
      if(TMath::Abs(histo_Zg_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb) -histo_Zg ->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_Zg_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_Zg ->GetBinContent(nb));
      if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb));
      if(TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb) -histo_WZ ->GetBinContent(nb)) > systQCDScale[2]) systQCDScale[2] = TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb));
      if(TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb) -histo_ZZ ->GetBinContent(nb)) > systQCDScale[3]) systQCDScale[3] = TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb));
    }                 
    if(histo_Zg ->GetBinContent(nb) != 0) systQCDScale[0] = 1 + systQCDScale[0]/histo_Zg ->GetBinContent(nb); else systQCDScale[0] = 1;
    if(histo_VVV->GetBinContent(nb) != 0) systQCDScale[1] = 1 + systQCDScale[1]/histo_VVV->GetBinContent(nb); else systQCDScale[1] = 1;
    if(histo_WZ ->GetBinContent(nb) != 0) systQCDScale[2] = 1 + systQCDScale[2]/histo_WZ ->GetBinContent(nb); else systQCDScale[2] = 1;
    if(histo_ZZ ->GetBinContent(nb) != 0) systQCDScale[3] = 1 + systQCDScale[3]/histo_ZZ ->GetBinContent(nb); else systQCDScale[3] = 1;
    printf("QCDScale(%d): %f %f %f %f\n",nb,systQCDScale[0],systQCDScale[1],systQCDScale[2],systQCDScale[3]);
    
    // PDF study
    double systPDF[5];
    histo_Diff->Reset();
    if(histo_Zg ->GetBinContent(nb) > 0) for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_Zg_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_Zg ->GetBinContent(nb))/histo_Zg ->GetBinContent(nb));
    systPDF[0] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    if(histo_VVV->GetBinContent(nb) > 0) for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_VVV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))/histo_VVV->GetBinContent(nb));
    systPDF[1] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    if(histo_WZ ->GetBinContent(nb) > 0) for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_WZ_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb))/histo_WZ ->GetBinContent(nb));
    systPDF[2] = 1.0+histo_Diff->GetRMS();
    histo_Diff->Reset();
    if(histo_ZZ ->GetBinContent(nb) > 0) for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_ZZ_CMS_PDFBounding[npdf] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb))/histo_ZZ ->GetBinContent(nb));
    systPDF[3] = 1.0+histo_Diff->GetRMS();
    if(histo_Higgs ->GetBinContent(nb) > 0) for(int npdf=1; npdf<102; npdf++) histo_Diff->Fill((histo_Higgs_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_Higgs ->GetBinContent(nb))/histo_Higgs ->GetBinContent(nb));
    systPDF[4] = 1.0+histo_Diff->GetRMS(); 
    printf("PDF(%d): %f %f %f %f %f\n",nb,systPDF[0],systPDF[1],systPDF[2],systPDF[3],systPDF[4]);

    double systLepEffM[4] = {1.0,1.0,1.0,1.0};
    if     (histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_Zg_CMS_MVALepEffMBoundingUp   ->GetBinContent(nb) > 0) systLepEffM[0] = histo_Zg_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_Zg_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[0] = histo_Zg_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_Zg_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffMBoundingDown     ->GetBinContent(nb) > 0) systLepEffM[1] = histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_WZ_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_WZ_CMS_MVALepEffMBoundingDown	   ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_ZZ_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)  > 0 && histo_ZZ_CMS_MVALepEffMBoundingDown	   ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingDown->GetBinContent(nb);

    double systLepEffE[4] = {1.0,1.0,1.0,1.0};
    if     (histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_Zg_CMS_MVALepEffEBoundingUp   ->GetBinContent(nb) > 0) systLepEffE[0] = histo_Zg_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_Zg_CMS_MVALepEffEBoundingDown ->GetBinContent(nb) > 0) systLepEffE[0] = histo_Zg_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_Zg_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[1] = histo_VVV_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffEBoundingDown     ->GetBinContent(nb) > 0) systLepEffE[1] = histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_WZ_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[2] = histo_WZ_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_WZ_CMS_MVALepEffEBoundingDown	   ->GetBinContent(nb) > 0) systLepEffE[2] = histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_WZ_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_ZZ_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[3] = histo_ZZ_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)  > 0 && histo_ZZ_CMS_MVALepEffEBoundingDown	   ->GetBinContent(nb) > 0) systLepEffE[3] = histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingDown->GetBinContent(nb);

    double systMet[4] = {1.0,1.0,1.0,1.0};
    if     (histo_Zg->GetBinContent(nb) > 0 && histo_Zg_CMS_MVAMETBoundingUp   ->GetBinContent(nb) > 0) systMet[0] = histo_Zg_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    else if(histo_Zg->GetBinContent(nb) > 0 && histo_Zg_CMS_MVAMETBoundingDown ->GetBinContent(nb) > 0) systMet[0] = histo_Zg->GetBinContent(nb)/histo_Zg_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[1] = histo_VVV_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAMETBoundingDown     ->GetBinContent(nb) > 0) systMet[1] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_WZ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[2] = histo_WZ_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVAMETBoundingDown      ->GetBinContent(nb) > 0) systMet[2] = histo_WZ->GetBinContent(nb)/histo_WZ_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[3] = histo_ZZ_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVAMETBoundingDown      ->GetBinContent(nb) > 0) systMet[3] = histo_ZZ->GetBinContent(nb)/histo_ZZ_CMS_MVAMETBoundingDown->GetBinContent(nb);

    double systPU[4] = {1.0,1.0,1.0,1.0};
    if     (histo_Zg->GetBinContent(nb) > 0 && histo_Zg_CMS_PUBoundingUp   ->GetBinContent(nb) > 0) systPU[0] = histo_Zg_CMS_PUBoundingUp->GetBinContent(nb)/histo_Zg->GetBinContent(nb);
    else if(histo_Zg->GetBinContent(nb) > 0 && histo_Zg_CMS_PUBoundingDown ->GetBinContent(nb) > 0) systPU[0] = histo_Zg->GetBinContent(nb)/histo_Zg_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_PUBoundingUp  ->GetBinContent(nb) > 0) systPU[1] = histo_VVV_CMS_PUBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_PUBoundingDown->GetBinContent(nb) > 0) systPU[1] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_WZ->GetBinContent(nb) > 0 && histo_WZ_CMS_PUBoundingUp   ->GetBinContent(nb) > 0) systPU[2] = histo_WZ_CMS_PUBoundingUp->GetBinContent(nb)/histo_WZ->GetBinContent(nb);
    else if(histo_WZ->GetBinContent(nb) > 0 && histo_WZ_CMS_PUBoundingDown ->GetBinContent(nb) > 0) systPU[2] = histo_WZ->GetBinContent(nb)/histo_WZ_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_ZZ->GetBinContent(nb) > 0 && histo_ZZ_CMS_PUBoundingUp   ->GetBinContent(nb) > 0) systPU[3] = histo_ZZ_CMS_PUBoundingUp->GetBinContent(nb)/histo_ZZ->GetBinContent(nb);
    else if(histo_ZZ->GetBinContent(nb) > 0 && histo_ZZ_CMS_PUBoundingDown ->GetBinContent(nb) > 0) systPU[3] = histo_ZZ->GetBinContent(nb)/histo_ZZ_CMS_PUBoundingDown->GetBinContent(nb);
    for(int npu=0; npu<4; npu++) if(systPU[npu] > 1.02) systPU[npu] = 1.02;
    for(int npu=0; npu<4; npu++) if(systPU[npu] < 0.98) systPU[npu] = 0.98;

    char outputLimitsShape[200];
    sprintf(outputLimitsShape,"histo_limits_wz3l%2s_%4s_bin%d.txt",finalStateName,ECMsb.Data(),nb-1);
    ofstream newcardShape;
    newcardShape.open(outputLimitsShape);
    newcardShape << Form("imax 1 number of channels\n");
    newcardShape << Form("jmax * number of background\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");

    newcardShape << Form("shapes * * wz3l3l_input_13TeV2015BBB5.root histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n");
    newcardShape << Form("shapes data_obs * wz3l3l_input_13TeV2015BBB5.root histo_Data\n");
    newcardShape << Form("shapes Higgs * wz3l3l_input_13TeV2015BBB5.root histo_Higgs_M$MASS histo_Higgs_M$MASS_$SYSTEMATIC\n");
    newcardShape << Form("Observation %d\n", -1/*(int)histo_Data->GetBinContent(nb)*/);
    newcardShape << Form("bin wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d wz%2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
    newcardShape << Form("process Zg VVV WZ ZZ Fake Higgs\n");
    newcardShape << Form("process 1 2 5 3 4 0\n");
    newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",-1.,-1.,-1.,-1.,-1.,-1.) ;//,histo_Zg->GetBinContent(nb),histo_VVV->GetBinContent(nb),histo_WZ->GetBinContent(nb),histo_ZZ->GetBinContent(nb),histo_Fake->GetBinContent(nb)/*,histo_FakeE->GetBinContent(nb)*/,histo_Higgs->GetBinContent(nb));
    newcardShape << Form("%s                               lnN  %7.5f   %7.5f   -   %7.5f   -    %7.5f  \n","lumi_13TeV",lumiE,lumiE,lumiE,lumiE);	
    newcardShape << Form("%s                               lnN  %7.5f   %7.5f   -   %7.5f   -    %7.5f  \n","lumi_stab_13TeV2015",lumiEstab,lumiEstab,lumiEstab,lumiEstab);    
//    newcardShape << Form("norm_WZ_13TeV2015                                lnN    -      -    %7.5f    -    -    -  \n",1.23);	
    newcardShape << Form("CMS_eff_b_mistag_13TeV2015                       lnN  %7.5f   %7.5f   -   %7.5f   -  %7.5f\n",0.98,0.98,0.98,0.98);
    //    newcardShape << Form("CMS_eff_b_mistag_13TeV2016                       lnN  %7.5f   %7.5f   -   %7.5f   -  %7.5f\n",0.98,0.98,0.98,0.98);
    
    newcardShape << Form("CMS_QCDunc_scale                            shape    -    -    %7.5f   -    -    - \n", 1.);
    newcardShape << Form("CMS_ewkscale                            shape    -    -    %7.5f   -    -    - \n", 1.);
    //newcardShape << Form("CMS_eff_b_bjet                         lnN	-     %7.5f   -     -	  -    -    -  \n",1.07);		
    newcardShape << Form("%s                                     lnN  %7.5f   %7.5f   -   %7.5f   -  %7.5f\n",effMName,1.018,1.018,1.018,1.018);
    newcardShape << Form("%s                                     lnN  %7.5f   %7.5f   -   %7.5f   -  %7.5f\n",effEName,1.017,1.017,1.017,1.017);
    //newcardShape << Form("%s                                     lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",momMName,systLepResM[0],systLepResM[1],systLepResM[2],systLepResM[3]);
    //newcardShape << Form("%s                                     lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",momEName,systLepResE[0],systLepResE[1],systLepResE[2],systLepResE[3]);
    //newcardShape << Form("CMS_pu                                 lnN  %7.5f   %7.5f %7.5f %7.5f   -    -  \n",systPU[0],systPU[1],systPU[2],systPU[3]);
    newcardShape << Form("%s                                     shape  %7.5f   %7.5f %7.5f %7.5f   -   %7.5f \n","CMS_scale_j",1.,1.,1.,1.,1.);
    newcardShape << Form("%s                                     shape  %7.5f   %7.5f %7.5f %7.5f   -   %7.5f \n","CMS_res_met",1.,1.,1.,1.,1.);
    newcardShape << Form("%s                                     shape  %7.5f   %7.5f %7.5f %7.5f   -   %7.5f \n","CMS_scale_l",1.,1.,1.,1.,1.);
    newcardShape << Form("%s                                     lnN       -       -    -      -    -   %7.5f \n","ewk_qqbar",1.08);    
    newcardShape << Form("%s                                     shape     -       -    -      -   %7.5f   -  \n","Fake_2015_bbb_histo_Fake_bin_1",1.);    
    newcardShape << Form("%s                                     shape     -       -    -      -   %7.5f   -  \n","Fake_2015_bbb_histo_Fake_bin_2",1.);    
    newcardShape << Form("%s                                     shape     -       -    -      -   %7.5f   -  \n","Fake_2015_bbb_histo_Fake_bin_3",1.);    
    newcardShape << Form("%s                                     shape     -       -    -      -   %7.5f   -  \n","Fake_2015_bbb_histo_Fake_bin_4",1.);    
    newcardShape << Form("%s                                     shape     -       -    -      -   %7.5f   -  \n","Fake_2015_bbb_histo_Fake_bin_5",1.);    
    newcardShape << Form("%s                                     shape     -       -    -      -   %7.5f   -  \n","Fake_2015_bbb_histo_Fake_bin_6",1.);
    newcardShape << Form("%s                                     shape     -       -    -      -   %7.5f   -  \n","Fake_2015_bbb_histo_Fake_bin_7",1.);
    newcardShape << Form("%s                                     shape     -       -    -      -   %7.5f   -  \n","Fake_2015_bbb_histo_Fake_bin_8",1.);
    newcardShape << Form("%s                                     shape     -       -    -      -   %7.5f   -  \n","Fake_2015_bbb_histo_Fake_bin_9",1.);
    newcardShape << Form("%s                                     shape     -       -    -     %7.5f   -    - \n","ZZ_2015_bbb_histo_ZZ_bin_1",1.);
    newcardShape << Form("%s                                     shape     -       -    -     %7.5f   -    - \n","ZZ_2015_bbb_histo_ZZ_bin_2",1.);
    newcardShape << Form("%s                                     shape     -       -    -     %7.5f   -    - \n","ZZ_2015_bbb_histo_ZZ_bin_3",1.);
    newcardShape << Form("%s                                     shape     -       -    -     %7.5f   -    - \n","ZZ_2015_bbb_histo_ZZ_bin_4",1.);
    newcardShape << Form("%s                                     shape     -       -    -     %7.5f   -    - \n","ZZ_2015_bbb_histo_ZZ_bin_5",1.);
    newcardShape << Form("%s                                     shape     -       -    -     %7.5f   -    - \n","ZZ_2015_bbb_histo_ZZ_bin_6",1.);
    newcardShape << Form("%s                                     shape     -       -    -     %7.5f   -    - \n","ZZ_2015_bbb_histo_ZZ_bin_7",1.);
    newcardShape << Form("%s                                     shape     -       -    -     %7.5f   -    - \n","ZZ_2015_bbb_histo_ZZ_bin_8",1.);
    newcardShape << Form("%s                                     shape     -       -    -     %7.5f   -    - \n","ZZ_2015_bbb_histo_ZZ_bin_9",1.);
    newcardShape << Form("%s                                     shape     -    %7.5f   -       -     -    - \n","VVV_2015_bbb_histo_VVV_bin_1",1.);
    newcardShape << Form("%s                                     shape     -    %7.5f   -       -     -    - \n","VVV_2015_bbb_histo_VVV_bin_2",1.);
    newcardShape << Form("%s                                     shape     -    %7.5f   -       -     -    - \n","VVV_2015_bbb_histo_VVV_bin_3",1.);
    newcardShape << Form("%s                                     shape     -    %7.5f   -       -     -    - \n","VVV_2015_bbb_histo_VVV_bin_4",1.);
    newcardShape << Form("%s                                     shape     -    %7.5f   -       -     -    - \n","VVV_2015_bbb_histo_VVV_bin_5",1.);
    newcardShape << Form("%s                                     shape     -    %7.5f   -       -     -    - \n","VVV_2015_bbb_histo_VVV_bin_6",1.);
    newcardShape << Form("%s                                     shape     -    %7.5f   -       -     -    - \n","VVV_2015_bbb_histo_VVV_bin_7",1.);
    newcardShape << Form("%s                                     shape     -    %7.5f   -       -     -    - \n","VVV_2015_bbb_histo_VVV_bin_8",1.);
    newcardShape << Form("%s                                     shape     -    %7.5f   -       -     -    - \n","VVV_2015_bbb_histo_VVV_bin_9",1.);
    newcardShape << Form("%s                                     shape     -       -   %7.5f    -     -    -  \n","WZ_2015_bbb_histo_WZ_bin_1",1.);
    newcardShape << Form("%s                                     shape     -       -   %7.5f    -     -    -  \n","WZ_2015_bbb_histo_WZ_bin_4",1.);
    newcardShape << Form("%s                                     shape     -       -   %7.5f    -     -    -  \n","WZ_2015_bbb_histo_WZ_bin_5",1.);
    newcardShape << Form("%s                                     shape     -       -   %7.5f    -     -    -  \n","WZ_2015_bbb_histo_WZ_bin_6",1.);
    newcardShape << Form("%s                                     shape     -       -   %7.5f    -     -    -  \n","WZ_2015_bbb_histo_WZ_bin_7",1.);
    newcardShape << Form("%s                                     shape     -       -   %7.5f    -     -    -  \n","WZ_2015_bbb_histo_WZ_bin_8",1.);
    newcardShape << Form("%s                                     shape     -       -   %7.5f    -     -    -  \n","WZ_2015_bbb_histo_WZ_bin_9",1.);
    newcardShape << Form("%s                                     shape   %7.5f       -    -    -     -    -  \n","Zg_2015_bbb_histo_Zg_bin_2",1.);
    newcardShape << Form("%s                                     shape   %7.5f       -    -    -     -    -  \n","Zg_2015_bbb_histo_Zg_bin_3",1.);
    newcardShape << Form("%s                                     shape   %7.5f       -    -    -     -    -  \n","Zg_2015_bbb_histo_Zg_bin_4",1.);
    newcardShape << Form("%s                                     shape   %7.5f       -    -    -     -    -  \n","Zg_2015_bbb_histo_Zg_bin_5",1.);
    newcardShape << Form("%s                                     shape   %7.5f       -    -    -     -    -  \n","Zg_2015_bbb_histo_Zg_bin_6",1.);
    newcardShape << Form("%s                                     shape   %7.5f       -    -    -     -    -  \n","Zg_2015_bbb_histo_Zg_bin_7",1.);
    newcardShape << Form("%s                                     shape   %7.5f       -    -    -     -    -  \n","Zg_2015_bbb_histo_Zg_bin_8",1.);
    newcardShape << Form("%s                                     shape   %7.5f       -    -    -     -    -  \n","Zg_2015_bbb_histo_Zg_bin_9",1.);
    newcardShape << Form("pdf_qqbar                              lnN  %7.5f      %7.5f  -    %7.5f   -   %7.5f \n",systPDF[0],systPDF[1],systPDF[3],1.0129);
    newcardShape << Form("nuisance edit rename * * CMS_eff_m CMS_eff_m_13TeV2015\n");
    newcardShape << Form("nuisance edit rename * * CMS_eff_e CMS_eff_e_13TeV2015\n");
    newcardShape << Form("nuisance edit rename * * CMS_res_met CMS_res_met_13TeV2015\n");
    newcardShape << Form("nuisance edit rename * * CMS_scale_l CMS_scale_l_13TeV2015\n");   
	
     //JES shape
    //MET shape
    //LEPES Shape

    newcardShape.close();
  }
}
