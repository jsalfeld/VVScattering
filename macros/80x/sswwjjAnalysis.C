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

#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"

#include "MitAnalysisRunII/macros/80x/factors.h"
#include "MitAnalysisRunII/macros/80x/BTagCalibrationStandalone.cc"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"

const double effTimesXsFiducial[2] = {0.177093*0.0269642*1000, 0.157640*0.0269642*1000}; // eff x cross section x 1000 (pb -> fb)

//double WSSF[5]  = {1.253783,1.239172,0.944641,0.997775,1.099893};
//double WSSFE[5] = {0.373235,0.163163,0.053062,0.029079,0.037923};
//double WSSF[5]  = {1.159353,1.186796,0.932848,0.996142,1.097076};
//double WSSFE[5] = {0.399227,0.175151,0.056634,0.032720,0.042750};
double WSSF[5]  = {1.131470,1.150168,0.928350,0.982761,1.100479};
double WSSFE[5] = {0.410094,0.177351,0.057828,0.033097,0.043396};

double the_sf_ZLL = 0.80;
const double bTagCuts[2] = {0.8484,0.9535}; // 0.5426/0.8484/0.9535 (check BTagCalibration2Reader!)

void func_ws_sf(double eta, double pt, double theSF[2]);
double func_ws_eff(double eta1, double eta2, TH1D *fhEff);

enum selType                     {SIGSEL=0, TOPSEL,   WZSEL,   DILSEL,   SS2JSEL,   OS2JSEL,   SSZLLSEL, nSelTypes};
TString selTypeName[nSelTypes]= {"SIGSEL", "TOPSEL", "WZSEL", "DILSEL", "SS2JSEL", "OS2JSEL", "SSZLLSEL"};

enum systType                     {JESUP=0, JESDOWN,  JERUP,  JERDOWN,  METUP,  METDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","JERUP","JERDOWN","METUP","METDOWN"};

bool verbose = true;
bool isMINIAOD = true;
int whichSkim1 = 6;
int whichSkim2 = 7;
double mcPrescale = 1.0;
bool usePureMC = false;
int period = 1;
//const TString typeLepSel = "verytight";
const bool usePUPPI = false;
const bool useWSFromData = true;
const bool useWZFromData = false;
const double mjjCut = 500.;

void sswwjjAnalysis(
 int theControlRegion = 0, 
 TString typeLepSel = "verytight", 
 int finalVar = 0, // 0 == mjj 4 x 6, 1 == mll 4 x 6, 2 == mll  5 x 2, 3 == mll 5 x 1, 4 == mjj vs. mll 4 x 5
 bool isBlinded = false,
 bool isMIT = false
 ){

  // File instances on EOS
  TString filesPathDA  = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/met_";
  TString filesPathMC  = "root://eoscms.cern.ch//eos/cms/store/caf/user/ceballos/Nero/output_80x/met_";
  TString filesPathMC2 = "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/ceballos/Nero/output_80x/mc/met_";
  // File instances on T3 hadoop
  if(isMIT){
    filesPathDA   = "/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/data/met_";
    filesPathMC	  = "/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/mc/met_";
    filesPathMC2  = "/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/mc/met_";
  }
  Double_t lumi = 35.9;

  if(typeLepSel == "medium_mva") {
    WSSF[0] = 0.762245; WSSFE[0] = 0.351553;
    WSSF[1] = 1.271616; WSSFE[1] = 0.180645;
    WSSF[2] = 0.915736; WSSFE[2] = 0.059091;
    WSSF[3] = 0.923135; WSSFE[3] = 0.030046;
    WSSF[4] = 1.038516; WSSFE[4] = 0.040860;
  }

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev, signalName_;  
  vector<Int_t> infilecatv, signalIndex_;  

  TString puPath = "";
  TString triggerSuffix = "*";
  if(isMINIAOD) triggerSuffix = "";
  if      (period==1){
  puPath = "MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root";

  //data samples
  if(isMINIAOD) {
    infilenamev.push_back(Form("%sdata_Run2016B.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016C.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016D.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016E.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016F.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016G.root",filesPathDA.Data())); infilecatv.push_back(0);
    infilenamev.push_back(Form("%sdata_Run2016H.root",filesPathDA.Data())); infilecatv.push_back(0);
  } else {
  }

  //MC samples
  //signal: EWK + QCD
  //infilenamev.push_back(Form("%sWpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.root",filesPathMC.Data()));               infilecatv.push_back(1);
  infilenamev.push_back(Form("%sWpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8.root",filesPathMC.Data()));                   infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sWpWpJJ_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));			       infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sWmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));			       infilecatv.push_back(1);

  //QCD to be subtracted from signal
  //infilenamev.push_back(Form("%sWpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.root",filesPathMC.Data()));                   infilecatv.push_back(-1);

  //QCD to be added to background
  infilenamev.push_back(Form("%sWpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.root",filesPathMC.Data()));                   infilecatv.push_back(2);

  //WZ
  //infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root",filesPathMC.Data()));                             infilecatv.push_back(3); 
  //infilenamev.push_back(Form("%sWZTo3LNu_mllmin01_13TeV-powheg-pythia8.root",filesPathMC.Data()));                                 infilecatv.push_back(3); 
  infilenamev.push_back(Form("%sWZTo3LNu_0Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_1Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_2Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_3Jets_MLL-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_0Jets_MLL-4To50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_1Jets_MLL-4To50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_2Jets_MLL-4To50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWZTo3LNu_3Jets_MLL-4To50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	   infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWLLJJ_WToLNu_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8.root",filesPathMC.Data()));           infilecatv.push_back(3);
  infilenamev.push_back(Form("%sWLLJJ_WToLNu_MLL-4To60_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8.root",filesPathMC.Data())); infilecatv.push_back(3);

  //ZZ
  //infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8.root",filesPathMC.Data()));				     infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));	             infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));	             infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));	             infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		     infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));  	             infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data())); 	             infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZJJTo4L_EWK_13TeV-madgraph-pythia8.root",filesPathMC.Data()));                              infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZJJTo4L_QCD_13TeV-madgraph-pythia8.root",filesPathMC.Data()));                              infilecatv.push_back(4);

  //VVV
  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(5);
  infilenamev.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));				      infilecatv.push_back(5);

  //Wrong sign
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg.root",filesPathMC.Data()));                                            infilecatv.push_back(6);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV.root",filesPathMC.Data()));                                        infilecatv.push_back(6);

  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data()));                       infilecatv.push_back(6);

  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8.root",filesPathMC.Data())); 				      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg.root",filesPathMC2.Data()));			                      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));    infilecatv.push_back(6);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));infilecatv.push_back(6);

  //infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));        infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));            infilecatv.push_back(6);
  infilenamev.push_back(Form("%sDY1JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));       infilecatv.push_back(6);
  infilenamev.push_back(Form("%sDY2JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));       infilecatv.push_back(6);
  infilenamev.push_back(Form("%sDY3JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));       infilecatv.push_back(6);
  infilenamev.push_back(Form("%sDY4JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));       infilecatv.push_back(6);
  infilenamev.push_back(Form("%sDY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sDY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sDY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sDY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(6);

  //Wgamma
  //infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));		        infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));		        infilecatv.push_back(7);
  infilenamev.push_back(Form("%sWGJJToLNu_EWK_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8.root",filesPathMC.Data()));             infilecatv.push_back(7);

  //DPS
  infilenamev.push_back(Form("%sWWTo2L2Nu_DoubleScattering_13TeV-pythia8.root",filesPathMC.Data()));                          infilecatv.push_back(8);

  //Non-prompt leptons
  //infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));                infilecatv.push_back(9);
  //infilenamev.push_back(Form("%sTTToSemiLeptonic_13TeV-powheg.root",filesPathMC2.Data()));                                    infilecatv.push_back(9);

  for(int ifile=0; ifile<(int)infilenamev.size(); ifile++) {
    signalIndex_.push_back(-1); // Populate vector of signal indices with -1 for the non-MC-signal files
  }

  {
  int i=0;
  signalName_.push_back(Form("mh%d", 200)); infilenamev.push_back(Form("%sDoublyChargedHiggsGMmodel_HWW_M200_13TeV-madgraph.root",filesPathMC.Data())); infilecatv.push_back(11); signalIndex_.push_back(i); i++;
  signalName_.push_back(Form("mh%d", 600)); infilenamev.push_back(Form("%sDoublyChargedHiggsGMmodel_HWW_M600_13TeV-madgraph.root",filesPathMC.Data())); infilecatv.push_back(11); signalIndex_.push_back(i); i++;
  signalName_.push_back(Form("mh%d", 300)); infilenamev.push_back(Form("%sDoublyChargedHiggsGMmodel_HWW_M300_13TeV-madgraph.root",filesPathMC.Data())); infilecatv.push_back(11); signalIndex_.push_back(i); i++;
  signalName_.push_back(Form("mh%d", 400)); infilenamev.push_back(Form("%sDoublyChargedHiggsGMmodel_HWW_M400_13TeV-madgraph.root",filesPathMC.Data())); infilecatv.push_back(11); signalIndex_.push_back(i); i++;
  signalName_.push_back(Form("mh%d", 500)); infilenamev.push_back(Form("%sDoublyChargedHiggsGMmodel_HWW_M500_13TeV-madgraph.root",filesPathMC.Data())); infilecatv.push_back(11); signalIndex_.push_back(i); i++;
  signalName_.push_back(Form("mh%d", 700)); infilenamev.push_back(Form("%sDoublyChargedHiggsGMmodel_HWW_M700_13TeV-madgraph.root",filesPathMC.Data())); infilecatv.push_back(11); signalIndex_.push_back(i); i++;
  signalName_.push_back(Form("mh%d", 800)); infilenamev.push_back(Form("%sDoublyChargedHiggsGMmodel_HWW_M800_13TeV-madgraph.root",filesPathMC.Data())); infilecatv.push_back(11); signalIndex_.push_back(i); i++;
  signalName_.push_back(Form("mh%d", 900)); infilenamev.push_back(Form("%sDoublyChargedHiggsGMmodel_HWW_M900_13TeV-madgraph.root",filesPathMC.Data())); infilecatv.push_back(11); signalIndex_.push_back(i); i++;
  signalName_.push_back(Form("mh%d",1000)); infilenamev.push_back(Form("%sDoublyChargedHiggsGMmodel_HWW_M1000_13TeV-madgraph.root",filesPathMC.Data()));infilecatv.push_back(11); signalIndex_.push_back(i); i++;
  }

  } // end period == 1
  else {assert(0); return;}

  //infilenamev.clear();infilecatv.clear();
  //infilenamev.push_back(Form("nero.root")); infilecatv.push_back(6);

  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}
  
  int nSigModels=signalName_.size();

  double selectedFiducial[2] = {0, 0};
  double selectedNonFiducial[2] = {0, 0};

  double denBTagging[5][5][3],jetEpsBtagLOOSE[5][5][3],jetEpsBtagTIGHT[5][5][3];
  double numBTaggingLOOSE[5][5][3],numBTaggingTIGHT[5][5][3];
  for(int i0=0; i0<5; i0++) {
    for(int i1=0; i1<5; i1++) {
      for(int i2=0; i2<3; i2++) {
        denBTagging[i0][i1][i2] = 0.0;
        numBTaggingLOOSE[i0][i1][i2] = 0.0;numBTaggingTIGHT[i0][i1][i2] = 0.0;
	if     (i2==BTagEntry::FLAV_B)    jetEpsBtagLOOSE[i0][i1][i2] = jetEpsBtagBMEDIUM[i0][i1];
	else if(i2==BTagEntry::FLAV_C)    jetEpsBtagLOOSE[i0][i1][i2] = jetEpsBtagCMEDIUM[i0][i1];
	else if(i2==BTagEntry::FLAV_UDSG) jetEpsBtagLOOSE[i0][i1][i2] = jetEpsBtagLMEDIUM[i0][i1];
	if     (i2==BTagEntry::FLAV_B)    jetEpsBtagTIGHT[i0][i1][i2] = jetEpsBtagBTIGHT[i0][i1];
	else if(i2==BTagEntry::FLAV_C)    jetEpsBtagTIGHT[i0][i1][i2] = jetEpsBtagCTIGHT[i0][i1];
	else if(i2==BTagEntry::FLAV_UDSG) jetEpsBtagTIGHT[i0][i1][i2] = jetEpsBtagLTIGHT[i0][i1];
      }
    }
  }
  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);

  BTagCalibration2 *btagCalib = new BTagCalibration2("csvv2","MitAnalysisRunII/data/80x/CSVv2_Moriond17_B_H.csv");

  BTagCalibration2Reader btagReaderBCLOOSE(btagCalib,BTagEntry::OP_MEDIUM,"comb","central");
  BTagCalibration2Reader btagReaderLLOOSE(btagCalib,BTagEntry::OP_MEDIUM,"incl","central");
  BTagCalibration2Reader btagReaderBCLOOSEUP(btagCalib,BTagEntry::OP_MEDIUM,"comb","up");
  BTagCalibration2Reader btagReaderLLOOSEUP(btagCalib,BTagEntry::OP_MEDIUM,"incl","up");
  BTagCalibration2Reader btagReaderBCLOOSEDOWN(btagCalib,BTagEntry::OP_MEDIUM,"comb","down");
  BTagCalibration2Reader btagReaderLLOOSEDOWN(btagCalib,BTagEntry::OP_MEDIUM,"incl","down");

  BTagCalibration2Reader btagReaderBCTIGHT(btagCalib,BTagEntry::OP_TIGHT,"comb","central");
  BTagCalibration2Reader btagReaderLTIGHT(btagCalib,BTagEntry::OP_TIGHT,"incl","central");
  BTagCalibration2Reader btagReaderBCTIGHTUP(btagCalib,BTagEntry::OP_TIGHT,"comb","up");
  BTagCalibration2Reader btagReaderLTIGHTUP(btagCalib,BTagEntry::OP_TIGHT,"incl","up");
  BTagCalibration2Reader btagReaderBCTIGHTDOWN(btagCalib,BTagEntry::OP_TIGHT,"comb","down");
  BTagCalibration2Reader btagReaderLTIGHTDOWN(btagCalib,BTagEntry::OP_TIGHT,"incl","down");

  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
  delete fPUFile;

  TFile *fTrackElectronReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root"));
  TH2D *fhDeltrksf= (TH2D*)(fTrackElectronReco_SF->Get("scalefactors_Reco_Electron")); assert(fhDeltrksf); fhDeltrksf->SetDirectory(0);
  delete fTrackElectronReco_SF;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("scalefactors_Medium_Electron"));
  TH2D *fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_Tight_Electron"));
  if(typeLepSel == "medium_mva") fhDElMediumSF = (TH2D*)(fElSF->Get("scalefactors_MediumMVA_Electron"));
  if(typeLepSel == "default_mva") fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_TightMVA_Electron"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF->SetDirectory(0);
  delete fElSF;

  TFile *fElVeryTightSF = TFile::Open(Form("MitAnalysisRunII/data/80x/veryTightSF_37ifb.root"));
  TH1D *fhDVeryTightSF = (TH1D*)(fElVeryTightSF->Get("veryTightSF"));
  assert(fhDVeryTightSF);
  fhDVeryTightSF->SetDirectory(0);
  delete fElVeryTightSF;

  //TFile *fElveryTightWrongSignEff = TFile::Open(Form("MitAnalysisRunII/data/80x/veryTightWrongSignEff_37ifb.root"));
  //TH1D *fhDveryTightWrongSignEff = (TH1D*)(fElveryTightWrongSignEff->Get("veryTightWrongSignEff"));
  //assert(fhDveryTightWrongSignEff);
  //fhDveryTightWrongSignEff->SetDirectory(0);
  //delete fElveryTightWrongSignEff;

  TFile *fTrackMuonReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/Tracking_EfficienciesAndSF_BCDEFGH.root"));
  TH1D *fhDmutrksfptg10 = (TH1D*)(fTrackMuonReco_SF->Get("ratio_eff_eta3_dr030e030_corr")); assert(fhDmutrksfptg10); fhDmutrksfptg10->SetDirectory(0);
  delete fTrackMuonReco_SF;

  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("scalefactors_TightId_Muon")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  //TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonID_Z_RunBCD_prompt80X_7p65.root"));
  //TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuMediumSF); fhDMuMediumSF->SetDirectory(0);
  delete fMuSF;

  TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("scalefactors_Iso_MuonTightId")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
  //TFile *fMuIsoSF = TFile::Open(Form("MitAnalysisRunII/data/80x/MuonIso_Z_RunBCD_prompt80X_7p65.root"));
  //TH2D *fhDMuIsoSF = (TH2D*)(fMuIsoSF->Get("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio")); assert(fhDMuIsoSF); fhDMuIsoSF->SetDirectory(0);
  delete fMuIsoSF;

  TString ECMsb  = "13TeV2016";

  const int nBinMVA = 24; Float_t xbins[nBinMVA+1] = {500, 800, 1100, 1500, 2000,
                                                          2800, 3100, 3500, 4000,
							  4800, 5100, 5500, 6000,
							  6800, 7100, 7500, 8000,
							  8800, 9100, 9500,10000,
							 10800,11100,11500,12000};

  const int nBinWZMVA = 5; Float_t xbinsWZ[nBinWZMVA+1] = {500, 800, 1100, 1500, 2000, 2000.00001};
  int nBinWZMVAModule = 4;

  if     (finalVar == 1){
    xbins[ 0] = 0; xbins[ 1] = 100; xbins[ 2] = 200; xbins[ 3] = 300; xbins[ 4] = 600;
     		   xbins[ 5] =1100; xbins[ 6] =1200; xbins[ 7] =1300; xbins[ 8] =1600;
     		   xbins[ 9] =2100; xbins[10] =2200; xbins[11] =2300; xbins[12] =2600;
     		   xbins[13] =3100; xbins[14] =3200; xbins[15] =3300; xbins[16] =3600;
     		   xbins[17] =4100; xbins[18] =4200; xbins[19] =4300; xbins[20] =4600;
     		   xbins[21] =5100; xbins[22] =5200; xbins[23] =5300; xbins[24] =5600;

    xbinsWZ[ 0] = 0; xbinsWZ[ 1] = 100; xbinsWZ[ 2] = 200; xbinsWZ[ 3] = 300; xbinsWZ[ 4] = 600; xbinsWZ[ 5] = 600.00001;
    nBinWZMVAModule = 4;
  }
  else if(finalVar == 2){
    xbins[ 0] = 0; xbins[ 1] = 100; xbins[ 2] = 200; xbins[ 3] = 300; xbins[ 4] = 400; xbins[ 5] = 600;
     		   xbins[ 6] =1100; xbins[ 7] =1200; xbins[ 8] =1300; xbins[ 9] =1400; xbins[10] =1600;
                   // not used!
     		   xbins[11] =2100; xbins[12] =2200; xbins[13] =2300; xbins[14] =2400; xbins[15] =2600;
     		   xbins[16] =3100; xbins[17] =3200; xbins[18] =3300; xbins[19] =3400; xbins[20] =3600;
     		   xbins[21] =4100; xbins[22] =4200; xbins[23] =4300; xbins[24] =4400;

    xbinsWZ[ 0] = 0; xbinsWZ[ 1] = 100; xbinsWZ[ 2] = 200; xbinsWZ[ 3] = 300; xbinsWZ[ 4] = 400; xbinsWZ[ 5] = 600;
    nBinWZMVAModule = 5;
  }
  else if(finalVar == 3){
    xbins[ 0] = 0; xbins[ 1] = 100; xbins[ 2] = 200; xbins[ 3] = 300; xbins[ 4] = 400; xbins[ 5] = 600;
                   // not used!
     		   xbins[ 6] =1100; xbins[ 7] =1200; xbins[ 8] =1300; xbins[ 9] =1400; xbins[10] =1600;
     		   xbins[11] =2100; xbins[12] =2200; xbins[13] =2300; xbins[14] =2400; xbins[15] =2600;
     		   xbins[16] =3100; xbins[17] =3200; xbins[18] =3300; xbins[19] =3400; xbins[20] =3600;
     		   xbins[21] =4100; xbins[22] =4200; xbins[23] =4300; xbins[24] =4400;

    xbinsWZ[ 0] = 0; xbinsWZ[ 1] = 100; xbinsWZ[ 2] = 200; xbinsWZ[ 3] = 300; xbinsWZ[ 4] = 400; xbinsWZ[ 5] = 600;
    nBinWZMVAModule = 5;
  }
  else if(finalVar == 4){ // just to show nothing needs to be done
  }

  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();
  TH1D* histoWZMVA = new TH1D("histoWZMVA", "histoWZMVA", nBinWZMVA, xbinsWZ);
  histoWZMVA->Sumw2();
  TH1D* histoOneBin = new TH1D("histoOneBin", "histoOneBin", 1, -0.5, 0.5);
  histoOneBin->Sumw2();

  TH1D *histo_Data  = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_EWK   = (TH1D*) histoMVA->Clone("histo_EWK"); 
  TH1D *histo_QCD   = (TH1D*) histoMVA->Clone("histo_QCD");    
  TH1D *histo_WZ    = (TH1D*) histoMVA->Clone("histo_WZ");
  TH1D *histo_ZZ    = (TH1D*) histoMVA->Clone("histo_ZZ");
  TH1D *histo_VVV   = (TH1D*) histoMVA->Clone("histo_VVV");
  TH1D *histo_WS    = (TH1D*) histoMVA->Clone("histo_WS");
  TH1D *histo_WG    = (TH1D*) histoMVA->Clone("histo_WG");
  TH1D *histo_DPS   = (TH1D*) histoMVA->Clone("histo_DPS");
  TH1D *histo_FakeM = (TH1D*) histoMVA->Clone("histo_FakeM");  
  TH1D *histo_FakeE = (TH1D*) histoMVA->Clone("histo_FakeE");  
  TH1D *histo_Higgs[nSigModels];
  for(int nModel=0; nModel<nSigModels; nModel++) {
    histo_Higgs[nModel] = (TH1D*) histoMVA->Clone(Form("histo_Higgs_%s",   signalName_[nModel].Data())); 
  }

  TH1D *histoOneBin_Data  = (TH1D*) histoOneBin->Clone("histoOneBin_Data");
  TH1D *histoOneBin_EWK   = (TH1D*) histoOneBin->Clone("histoOneBin_EWK"); 
  TH1D *histoOneBin_QCD   = (TH1D*) histoOneBin->Clone("histoOneBin_QCD");    
  TH1D *histoOneBin_WZ    = (TH1D*) histoOneBin->Clone("histoOneBin_WZ");
  TH1D *histoOneBin_ZZ    = (TH1D*) histoOneBin->Clone("histoOneBin_ZZ");
  TH1D *histoOneBin_VVV   = (TH1D*) histoOneBin->Clone("histoOneBin_VVV");
  TH1D *histoOneBin_WS    = (TH1D*) histoOneBin->Clone("histoOneBin_WS");
  TH1D *histoOneBin_WG    = (TH1D*) histoOneBin->Clone("histoOneBin_WG");
  TH1D *histoOneBin_DPS   = (TH1D*) histoOneBin->Clone("histoOneBin_DPS");
  TH1D *histoOneBin_FakeM = (TH1D*) histoOneBin->Clone("histoOneBin_FakeM");  
  TH1D *histoOneBin_FakeE = (TH1D*) histoOneBin->Clone("histoOneBin_FakeE");  
  TH1D *histoOneBin_Higgs[nSigModels];
  for(int nModel=0; nModel<nSigModels; nModel++) {
    histoOneBin_Higgs[nModel] = (TH1D*) histoMVA->Clone(Form("histoOneBin_Higgs")); 
  }

  double totalFakeDataCount[6][5];
  for(int i=0; i<6; i++) for(int j=0; j<5; j++) totalFakeDataCount[i][j] = 0;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 41;
  const int histBins = 13;
  TH1D* histo[7][allPlots][histBins];
  TString processName[histBins] = {".Data", "EWKWW", "QCDWW", "...WZ", "...ZZ", "..VVV", "...WS", "...WG", "..DPS", "FakeM", "FakeE", "..Hig1", "..Hig2"};

  const int nBinMJJMVA = 4; Float_t xbinsMJJ[nBinMJJMVA+1] = {500, 800, 1100, 1500, 2000};
  const int nBinMLLMVA = 4; Float_t xbinsMLL[nBinMLLMVA+1] = {0, 100, 180, 300, 600};

  for(int nState=0; nState<7; nState++){
    for(int thePlot=0; thePlot<allPlots; thePlot++){
      if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 2000;}
      else if(thePlot >=  1 && thePlot <=  1) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot = 8;}
      else if(thePlot >=  2 && thePlot <=  2) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot = 6.5;}
      else if(thePlot >=  3 && thePlot <=  4) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400;}
      else if(thePlot >=  5 && thePlot <=  5) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot = 39.5;}
      else if(thePlot >=  6 && thePlot <=  6) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot = 3.5;}
      else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot = 40;}
      else if(thePlot >=  8 && thePlot <=  8) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot = 50;}
      else if(thePlot >=  9 && thePlot <=  9) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 1.0;}
      else if(thePlot >= 10 && thePlot <= 10) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot = 3.5;}
      else if(thePlot >= 11 && thePlot <= 11) {}
      else if(thePlot >= 12 && thePlot <= 12) {}
      else if(thePlot >= 13 && thePlot <= 14) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200;}
      else if(thePlot >= 15 && thePlot <= 15) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot = 8;}
      else if(thePlot >= 16 && thePlot <= 17) {nBinPlot = 100; xminPlot =-5.0; xmaxPlot = 5.0;}
      else if(thePlot >= 18 && thePlot <= 18) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 400;}
      else if(thePlot >= 19 && thePlot <= 19) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 1.0;}
      else if(thePlot >= 20 && thePlot <= 20) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot = 20;}
      else if(thePlot >= 21 && thePlot <= 22) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 1.0;}
      else if(thePlot >= 23 && thePlot <= 24) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200;}
      else if(thePlot >= 25 && thePlot <= 25) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot = 6.5;}
      else if(thePlot >= 26 && thePlot <= 26) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 500;}
      else if(thePlot >= 27 && thePlot <= 27) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 2000;}
      else if(thePlot >= 28 && thePlot <= 29) {nBinPlot =  80; xminPlot = 0.0; xmaxPlot = 8;}
      else if(thePlot >= 30 && thePlot <= 30) {nBinPlot =   4; xminPlot =-0.5; xmaxPlot = 3.5;}
      TH1D* histos;
      if     (thePlot == 11)         histos = new TH1D("histos", "histos", nBinMJJMVA, xbinsMJJ);
      else if(thePlot == 12)         histos = new TH1D("histos", "histos", nBinMLLMVA, xbinsMLL);
      else if(thePlot < allPlots-10) histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
      else if(thePlot < allPlots-7)  histos = new TH1D("histos", "histos", nBinMJJMVA, xbinsMJJ);
      else if(thePlot < allPlots-4)  histos = new TH1D("histos", "histos", nBinMLLMVA, xbinsMLL);
      else                           histos = new TH1D("histos", "histos", nBinMVA, xbins);
      histos->Sumw2();
      for(int i=0; i<histBins; i++) histo[nState][thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
      histos->Reset();histos->Clear();
    }
  }

  double bgdDecay[6][nSelTypes][histBins],weiDecay[6][nSelTypes][histBins];
  for(int nState=0; nState<6; nState++) { for(unsigned int i=0; i<nSelTypes; i++) { for(int j=0; j<histBins; j++) {       
    bgdDecay[nState][i][j] = 0.0; weiDecay[nState][i][j] = 0.0; 
  }}}

  char finalStateName[6],effMName[10],effEName[10],momMName[10],momEName[10],metName[10],jesName[10],jerName[10],puName[10],btagName[20],mistagName[20];
  sprintf(effMName,"CMS_eff2016_m");sprintf(momMName,"CMS_scale2016_m");
  sprintf(effEName,"CMS_eff2016_e");sprintf(momEName,"CMS_scale2016_e");
  sprintf(metName,"CMS_scale_met");sprintf(jesName,"CMS_scale_j");sprintf(jerName,"CMS_jer");
  sprintf(puName,"CMS_pu");sprintf(btagName,"CMS_eff_b_b2016");sprintf(mistagName,"CMS_eff_b_mistag2016");
  sprintf(finalStateName,"wwss");
  if     (theControlRegion == 1) sprintf(finalStateName,"top");
  else if(theControlRegion == 2) sprintf(finalStateName,"wz");

  TH1D* histo_EWK_CMS_MVAEWKStatBoundingUp       = new TH1D( Form("histo_EWK_CMS_wwss%s_MVAEWKStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_EWK_CMS_wwss%s_MVAEWKStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_EWK_CMS_MVAEWKStatBoundingUp  ->Sumw2();
  TH1D* histo_EWK_CMS_MVAEWKStatBoundingDown     = new TH1D( Form("histo_EWK_CMS_wwss%s_MVAEWKStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_EWK_CMS_wwss%s_MVAEWKStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_EWK_CMS_MVAEWKStatBoundingDown->Sumw2();
  TH1D* histo_QCD_CMS_MVAQCDStatBoundingUp       = new TH1D( Form("histo_QCD_CMS_wwss%s_MVAQCDStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_QCD_CMS_wwss%s_MVAQCDStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_QCD_CMS_MVAQCDStatBoundingUp  ->Sumw2();
  TH1D* histo_QCD_CMS_MVAQCDStatBoundingDown     = new TH1D( Form("histo_QCD_CMS_wwss%s_MVAQCDStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_QCD_CMS_wwss%s_MVAQCDStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_QCD_CMS_MVAQCDStatBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingUp         = new TH1D( Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAWZStatBoundingDown       = new TH1D( Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WZ_CMS_MVAWZStatBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingUp         = new TH1D( Form("histo_ZZ_CMS_wwss%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_wwss%s_MVAZZStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingDown       = new TH1D( Form("histo_ZZ_CMS_wwss%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ZZ_CMS_wwss%s_MVAZZStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ZZ_CMS_MVAZZStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp       = new TH1D( Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown     = new TH1D( Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVAWSStatBoundingUp         = new TH1D( Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WS_CMS_MVAWSStatBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAWSStatBoundingDown       = new TH1D( Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WS_CMS_MVAWSStatBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAWGStatBoundingUp         = new TH1D( Form("histo_WG_CMS_wwss%s_MVAWGStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WG_CMS_wwss%s_MVAWGStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WG_CMS_MVAWGStatBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAWGStatBoundingDown       = new TH1D( Form("histo_WG_CMS_wwss%s_MVAWGStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WG_CMS_wwss%s_MVAWGStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WG_CMS_MVAWGStatBoundingDown->Sumw2();
  TH1D* histo_DPS_CMS_MVADPSStatBoundingUp       = new TH1D( Form("histo_DPS_CMS_wwss%s_MVADPSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_DPS_CMS_wwss%s_MVADPSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_DPS_CMS_MVADPSStatBoundingUp  ->Sumw2();
  TH1D* histo_DPS_CMS_MVADPSStatBoundingDown     = new TH1D( Form("histo_DPS_CMS_wwss%s_MVADPSStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_DPS_CMS_wwss%s_MVADPSStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_DPS_CMS_MVADPSStatBoundingDown->Sumw2();
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingUp   = new TH1D( Form("histo_FakeM_CMS_wwss%s_MVAFakeMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_FakeM_CMS_wwss%s_MVAFakeMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeM_CMS_MVAFakeMStatBoundingUp  ->Sumw2();
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingDown = new TH1D( Form("histo_FakeM_CMS_wwss%s_MVAFakeMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_FakeM_CMS_wwss%s_MVAFakeMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeM_CMS_MVAFakeMStatBoundingDown->Sumw2();
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingUp   = new TH1D( Form("histo_FakeE_CMS_wwss%s_MVAFakeEStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_FakeE_CMS_wwss%s_MVAFakeEStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->Sumw2();
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingDown = new TH1D( Form("histo_FakeE_CMS_wwss%s_MVAFakeEStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_FakeE_CMS_wwss%s_MVAFakeEStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_FakeE_CMS_MVAFakeEStatBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingUp[nSigModels];  
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingDown[nSigModels];
  for(int nModel=0; nModel<nSigModels; nModel++) {
    histo_Higgs_CMS_MVAHiggsStatBoundingUp[nModel]      = new TH1D( Form("histo_Higgs_CMS_wwss%s_%s_MVAHiggsStatBounding_%sUp"  ,finalStateName, signalName_[nModel].Data(), ECMsb.Data()), Form("histo_Higgs_CMS_wwss%s_%s_MVAHiggsStatBounding_%sUp"  ,finalStateName, signalName_[nModel].Data(), ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingUp[nModel]  ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingDown[nModel]    = new TH1D( Form("histo_Higgs_CMS_wwss%s_%s_MVAHiggsStatBounding_%sDown",finalStateName, signalName_[nModel].Data(), ECMsb.Data()), Form("histo_Higgs_CMS_wwss%s_%s_MVAHiggsStatBounding_%sDown",finalStateName, signalName_[nModel].Data(), ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingDown[nModel]->Sumw2();
  }

  TH1D* histo_EWK_CMS_MVAEWKStatBoundingBinUp[nBinMVA];      
  TH1D* histo_EWK_CMS_MVAEWKStatBoundingBinDown[nBinMVA];    
  TH1D* histo_QCD_CMS_MVAQCDStatBoundingBinUp[nBinMVA];      
  TH1D* histo_QCD_CMS_MVAQCDStatBoundingBinDown[nBinMVA];    
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinUp[nBinMVA];        
  TH1D* histo_WZ_CMS_MVAWZStatBoundingBinDown[nBinMVA];      
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinUp[nBinMVA];        
  TH1D* histo_ZZ_CMS_MVAZZStatBoundingBinDown[nBinMVA];      
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp[nBinMVA];      
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown[nBinMVA];    
  TH1D* histo_WS_CMS_MVAWSStatBoundingBinUp[nBinMVA];        
  TH1D* histo_WS_CMS_MVAWSStatBoundingBinDown[nBinMVA];      
  TH1D* histo_WG_CMS_MVAWGStatBoundingBinUp[nBinMVA];        
  TH1D* histo_WG_CMS_MVAWGStatBoundingBinDown[nBinMVA];      
  TH1D* histo_DPS_CMS_MVADPSStatBoundingBinUp[nBinMVA];      
  TH1D* histo_DPS_CMS_MVADPSStatBoundingBinDown[nBinMVA];    
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[nBinMVA];  
  TH1D* histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[nBinMVA];
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nBinMVA];  
  TH1D* histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nSigModels][nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nSigModels][nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++) {
    histo_EWK_CMS_MVAEWKStatBoundingBinUp[nb]        = new TH1D(Form("histo_EWK_CMS_wwss%s_MVAEWKStatBounding_%s_Bin%dUp" 	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_EWK_CMS_wwss%s_MVAEWKStatBounding_%s_Bin%dUp"  	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_EWK_CMS_MVAEWKStatBoundingBinDown[nb]      = new TH1D(Form("histo_EWK_CMS_wwss%s_MVAEWKStatBounding_%s_Bin%dDown" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_EWK_CMS_wwss%s_MVAEWKStatBounding_%s_Bin%dDown"  	,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_QCD_CMS_MVAQCDStatBoundingBinUp[nb]        = new TH1D(Form("histo_QCD_CMS_wwss%s_MVAQCDStatBounding_%s_Bin%dUp" 	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_QCD_CMS_wwss%s_MVAQCDStatBounding_%s_Bin%dUp"  	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_QCD_CMS_MVAQCDStatBoundingBinDown[nb]      = new TH1D(Form("histo_QCD_CMS_wwss%s_MVAQCDStatBounding_%s_Bin%dDown" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_QCD_CMS_wwss%s_MVAQCDStatBounding_%s_Bin%dDown"  	,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]          = new TH1D(Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%s_Bin%dUp" 	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%s_Bin%dUp"  	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]        = new TH1D(Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%s_Bin%dDown" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_WZ_CMS_wwss%s_MVAWZStatBounding_%s_Bin%dDown"  	,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]          = new TH1D(Form("histo_ZZ_CMS_wwss%s_MVAZZStatBounding_%s_Bin%dUp" 	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_ZZ_CMS_wwss%s_MVAZZStatBounding_%s_Bin%dUp"  	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]        = new TH1D(Form("histo_ZZ_CMS_wwss%s_MVAZZStatBounding_%s_Bin%dDown" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_ZZ_CMS_wwss%s_MVAZZStatBounding_%s_Bin%dDown"  	,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]        = new TH1D(Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%s_Bin%dUp" 	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%s_Bin%dUp"  	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]      = new TH1D(Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%s_Bin%dDown" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_VVV_CMS_wwss%s_MVAVVVStatBounding_%s_Bin%dDown"  	,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_WS_CMS_MVAWSStatBoundingBinUp[nb]          = new TH1D(Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%s_Bin%dUp" 	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%s_Bin%dUp"  	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_WS_CMS_MVAWSStatBoundingBinDown[nb]        = new TH1D(Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%s_Bin%dDown" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_WS_CMS_wwss%s_MVAWSStatBounding_%s_Bin%dDown"  	,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_WG_CMS_MVAWGStatBoundingBinUp[nb]          = new TH1D(Form("histo_WG_CMS_wwss%s_MVAWGStatBounding_%s_Bin%dUp" 	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_WG_CMS_wwss%s_MVAWGStatBounding_%s_Bin%dUp"  	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_WG_CMS_MVAWGStatBoundingBinDown[nb]        = new TH1D(Form("histo_WG_CMS_wwss%s_MVAWGStatBounding_%s_Bin%dDown" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_WG_CMS_wwss%s_MVAWGStatBounding_%s_Bin%dDown"  	,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_DPS_CMS_MVADPSStatBoundingBinUp[nb]        = new TH1D(Form("histo_DPS_CMS_wwss%s_MVADPSStatBounding_%s_Bin%dUp" 	   ,finalStateName,  ECMsb.Data(),nb), Form("histo_DPS_CMS_wwss%s_MVADPSStatBounding_%s_Bin%dUp"  	  ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_DPS_CMS_MVADPSStatBoundingBinDown[nb]      = new TH1D(Form("histo_DPS_CMS_wwss%s_MVADPSStatBounding_%s_Bin%dDown" 	 ,finalStateName,  ECMsb.Data(),nb), Form("histo_DPS_CMS_wwss%s_MVADPSStatBounding_%s_Bin%dDown"  	,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[nb]    = new TH1D(Form("histo_FakeM_CMS_wwss%s_MVAFakeMStatBounding_%s_Bin%dUp"	       ,finalStateName,  ECMsb.Data(),nb), Form("histo_FakeM_CMS_wwss%s_MVAFakeMStatBounding_%s_Bin%dUp"      ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[nb]  = new TH1D(Form("histo_FakeM_CMS_wwss%s_MVAFakeMStatBounding_%s_Bin%dDown"      ,finalStateName,  ECMsb.Data(),nb), Form("histo_FakeM_CMS_wwss%s_MVAFakeMStatBounding_%s_Bin%dDown"    ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 
    histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nb]    = new TH1D(Form("histo_FakeE_CMS_wwss%s_MVAFakeEStatBounding_%s_Bin%dUp"	       ,finalStateName,  ECMsb.Data(),nb), Form("histo_FakeE_CMS_wwss%s_MVAFakeEStatBounding_%s_Bin%dUp"      ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins);
    histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nb]  = new TH1D(Form("histo_FakeE_CMS_wwss%s_MVAFakeEStatBounding_%s_Bin%dDown"      ,finalStateName,  ECMsb.Data(),nb), Form("histo_FakeE_CMS_wwss%s_MVAFakeEStatBounding_%s_Bin%dDown"    ,finalStateName,  ECMsb.Data(),nb),nBinMVA, xbins); 

    histo_EWK_CMS_MVAEWKStatBoundingBinUp[nb]	   ->Sumw2();
    histo_EWK_CMS_MVAEWKStatBoundingBinDown[nb]    ->Sumw2();
    histo_QCD_CMS_MVAQCDStatBoundingBinUp[nb]	   ->Sumw2();
    histo_QCD_CMS_MVAQCDStatBoundingBinDown[nb]    ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinUp[nb]	   ->Sumw2();
    histo_WZ_CMS_MVAWZStatBoundingBinDown[nb]	   ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[nb]	   ->Sumw2();
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[nb]	   ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	   ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WS_CMS_MVAWSStatBoundingBinUp[nb]	   ->Sumw2();
    histo_WS_CMS_MVAWSStatBoundingBinDown[nb]	   ->Sumw2();
    histo_WG_CMS_MVAWGStatBoundingBinUp[nb]	   ->Sumw2();
    histo_WG_CMS_MVAWGStatBoundingBinDown[nb]	   ->Sumw2();
    histo_DPS_CMS_MVADPSStatBoundingBinUp[nb]	   ->Sumw2();
    histo_DPS_CMS_MVADPSStatBoundingBinDown[nb]    ->Sumw2();
    histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[nb]  ->Sumw2();
    histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[nb]->Sumw2();
    histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[nb]  ->Sumw2();
    histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[nb]->Sumw2();

    for(int nModel=0; nModel<nSigModels; nModel++) { 
      histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nModel][nb]        = new TH1D(Form("histo_Higgs_CMS_wwss%s_%s_MVAHiggsStatBounding_%s_Bin%dUp"       ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb), Form("histo_Higgs_CMS_wwss%s_%s_MVAHiggsStatBounding_%s_Bin%dUp"       ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb),nBinMVA, xbins); 
      histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nModel][nb]      = new TH1D(Form("histo_Higgs_CMS_wwss%s_%s_MVAHiggsStatBounding_%s_Bin%dDown"     ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb), Form("histo_Higgs_CMS_wwss%s_%s_MVAHiggsStatBounding_%s_Bin%dDown"     ,finalStateName, signalName_[nModel].Data(), ECMsb.Data(),nb),nBinMVA, xbins); 
      histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nModel][nb]   ->Sumw2();
      histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nModel][nb] ->Sumw2();
    }
  }

  TH1D* histo_Diff = new TH1D("dummy", "dummy",1000,-1,1); histo_Diff->Sumw2();

  TH1D* histo_EWK_CMS_QCDScaleBounding[6];
  TH1D* histo_QCD_CMS_QCDScaleBounding[6];
  TH1D* histo_WZ_CMS_QCDScaleBounding [6];
  TH1D* histo_ZZ_CMS_QCDScaleBounding [6];
  TH1D* histo_VVV_CMS_QCDScaleBounding[6];
  TH1D* histo_WS_CMS_QCDScaleBounding [6];
  TH1D* histo_WG_CMS_QCDScaleBounding [6];
  TH1D* histo_DPS_CMS_QCDScaleBounding[6];
  TH1D* histo_Higgs_CMS_QCDScaleBounding[nSigModels][6];
  for(int nb=0; nb<6; nb++){
    histo_EWK_CMS_QCDScaleBounding[nb]        = new TH1D(Form("histo_EWK_QCDScale_f%d",nb),	 Form("histo_EWK_QCDScale_f%d",nb),nBinMVA, xbins);	 histo_EWK_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_QCD_CMS_QCDScaleBounding[nb]        = new TH1D(Form("histo_QCD_QCDScale_f%d",nb),	 Form("histo_QCD_QCDScale_f%d",nb),nBinMVA, xbins);	 histo_QCD_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WZ_CMS_QCDScaleBounding [nb]        = new TH1D(Form("histo_WZ_QCDScale_f%d",nb),	 Form("histo_WZ_QCDScale_f%d",nb),nBinMVA, xbins);	 histo_WZ_CMS_QCDScaleBounding [nb]->Sumw2();
    histo_ZZ_CMS_QCDScaleBounding [nb]        = new TH1D(Form("histo_ZZ_QCDScale_f%d",nb),	 Form("histo_ZZ_QCDScale_f%d",nb),nBinMVA, xbins);	 histo_ZZ_CMS_QCDScaleBounding [nb]->Sumw2();
    histo_VVV_CMS_QCDScaleBounding[nb]        = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),	 Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);	 histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WS_CMS_QCDScaleBounding [nb]        = new TH1D(Form("histo_WS_QCDScale_f%d",nb),	 Form("histo_WS_QCDScale_f%d",nb),nBinMVA, xbins);	 histo_WS_CMS_QCDScaleBounding [nb]->Sumw2();
    histo_WG_CMS_QCDScaleBounding [nb]        = new TH1D(Form("histo_WG_QCDScale_f%d",nb),	 Form("histo_WG_QCDScale_f%d",nb),nBinMVA, xbins);	 histo_WG_CMS_QCDScaleBounding [nb]->Sumw2();
    histo_DPS_CMS_QCDScaleBounding[nb]        = new TH1D(Form("histo_DPS_QCDScale_f%d",nb),	 Form("histo_DPS_QCDScale_f%d",nb),nBinMVA, xbins);	 histo_DPS_CMS_QCDScaleBounding[nb]->Sumw2();
    for(int nModel=0; nModel<nSigModels; nModel++) {
      histo_Higgs_CMS_QCDScaleBounding[nModel][nb]   = new TH1D(Form("histo_Higgs_%s_QCDScale_f%d", signalName_[nModel].Data(), nb), Form("histo_Higgs_%s_QCDScale_f%d", signalName_[nModel].Data(), nb),nBinMVA, xbins); histo_Higgs_CMS_QCDScaleBounding[nModel][nb]->Sumw2();
    }
  }

  TH1D* histo_EWK_CMS_PDFBounding[100];
  TH1D* histo_QCD_CMS_PDFBounding[100];
  TH1D* histo_WZ_CMS_PDFBounding [100];
  TH1D* histo_ZZ_CMS_PDFBounding [100];
  TH1D* histo_VVV_CMS_PDFBounding[100];
  TH1D* histo_WS_CMS_PDFBounding [100];
  TH1D* histo_WG_CMS_PDFBounding [100];
  TH1D* histo_DPS_CMS_PDFBounding[100];
  TH1D* histo_Higgs_CMS_PDFBounding[nSigModels][100];
  for(int nb=0; nb<100; nb++){
    histo_EWK_CMS_PDFBounding[nb]        = new TH1D(Form("histo_EWK_PDF_f%d",nb),	 Form("histo_EWK_PDF_f%d",nb),nBinMVA, xbins);	 histo_EWK_CMS_PDFBounding[nb]->Sumw2();
    histo_QCD_CMS_PDFBounding[nb]        = new TH1D(Form("histo_QCD_PDF_f%d",nb),	 Form("histo_QCD_PDF_f%d",nb),nBinMVA, xbins);	 histo_QCD_CMS_PDFBounding[nb]->Sumw2();
    histo_WZ_CMS_PDFBounding [nb]        = new TH1D(Form("histo_WZ_PDF_f%d",nb),	 Form("histo_WZ_PDF_f%d",nb),nBinMVA, xbins);	 histo_WZ_CMS_PDFBounding [nb]->Sumw2();
    histo_ZZ_CMS_PDFBounding [nb]        = new TH1D(Form("histo_ZZ_PDF_f%d",nb),	 Form("histo_ZZ_PDF_f%d",nb),nBinMVA, xbins);	 histo_ZZ_CMS_PDFBounding [nb]->Sumw2();
    histo_VVV_CMS_PDFBounding[nb]        = new TH1D(Form("histo_VVV_PDF_f%d",nb),	 Form("histo_VVV_PDF_f%d",nb),nBinMVA, xbins);	 histo_VVV_CMS_PDFBounding[nb]->Sumw2();
    histo_WS_CMS_PDFBounding [nb]        = new TH1D(Form("histo_WS_PDF_f%d",nb),	 Form("histo_WS_PDF_f%d",nb),nBinMVA, xbins);	 histo_WS_CMS_PDFBounding [nb]->Sumw2();
    histo_WG_CMS_PDFBounding [nb]        = new TH1D(Form("histo_WG_PDF_f%d",nb),	 Form("histo_WG_PDF_f%d",nb),nBinMVA, xbins);	 histo_WG_CMS_PDFBounding [nb]->Sumw2();
    histo_DPS_CMS_PDFBounding[nb]        = new TH1D(Form("histo_DPS_PDF_f%d",nb),	 Form("histo_DPS_PDF_f%d",nb),nBinMVA, xbins);	 histo_DPS_CMS_PDFBounding[nb]->Sumw2();
    for(int nModel=0; nModel<nSigModels; nModel++) {
      histo_Higgs_CMS_PDFBounding[nModel][nb]   = new TH1D(Form("histo_Higgs_%s_PDF_f%d", signalName_[nModel].Data(), nb), Form("histo_Higgs_%s_PDF_f%d", signalName_[nModel].Data(), nb),nBinMVA, xbins); histo_Higgs_CMS_PDFBounding[nModel][nb]->Sumw2();
    }
  }

  TH1D* histo_EWK_CMS_MVALepEffMBoundingUp    = new TH1D( Form("histo_EWK_%sUp",effMName)  , Form("histo_EWK_%sUp",effMName)  , nBinMVA, xbins); histo_EWK_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_EWK_CMS_MVALepEffMBoundingDown  = new TH1D( Form("histo_EWK_%sDown",effMName), Form("histo_EWK_%sDown",effMName), nBinMVA, xbins); histo_EWK_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_QCD_CMS_MVALepEffMBoundingUp    = new TH1D( Form("histo_QCD_%sUp",effMName)  , Form("histo_QCD_%sUp",effMName)  , nBinMVA, xbins); histo_QCD_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_QCD_CMS_MVALepEffMBoundingDown  = new TH1D( Form("histo_QCD_%sDown",effMName), Form("histo_QCD_%sDown",effMName), nBinMVA, xbins); histo_QCD_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingUp     = new TH1D( Form("histo_WZ_%sUp",effMName)  , Form("histo_WZ_%sUp",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingDown   = new TH1D( Form("histo_WZ_%sDown",effMName), Form("histo_WZ_%sDown",effMName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingUp     = new TH1D( Form("histo_ZZ_%sUp",effMName)  , Form("histo_ZZ_%sUp",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingDown   = new TH1D( Form("histo_ZZ_%sDown",effMName), Form("histo_ZZ_%sDown",effMName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingUp    = new TH1D( Form("histo_VVV_%sUp",effMName)  , Form("histo_VVV_%sUp",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingDown  = new TH1D( Form("histo_VVV_%sDown",effMName), Form("histo_VVV_%sDown",effMName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVALepEffMBoundingUp     = new TH1D( Form("histo_WS_%sUp",effMName)  , Form("histo_WS_%sUp",effMName)  , nBinMVA, xbins); histo_WS_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVALepEffMBoundingDown   = new TH1D( Form("histo_WS_%sDown",effMName), Form("histo_WS_%sDown",effMName), nBinMVA, xbins); histo_WS_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffMBoundingUp     = new TH1D( Form("histo_WG_%sUp",effMName)  , Form("histo_WG_%sUp",effMName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffMBoundingDown   = new TH1D( Form("histo_WG_%sDown",effMName), Form("histo_WG_%sDown",effMName), nBinMVA, xbins); histo_WG_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_DPS_CMS_MVALepEffMBoundingUp    = new TH1D( Form("histo_DPS_%sUp",effMName)  , Form("histo_DPS_%sUp",effMName)  , nBinMVA, xbins); histo_DPS_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_DPS_CMS_MVALepEffMBoundingDown  = new TH1D( Form("histo_DPS_%sDown",effMName), Form("histo_DPS_%sDown",effMName), nBinMVA, xbins); histo_DPS_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingUp[nSigModels]; 
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingDown[nSigModels];

  TH1D* histo_EWK_CMS_MVALepEffMBoundingAvg   = new TH1D( Form("histo_EWK_%sAvg",effMName)  , Form("histo_EWK_%sAvg",effMName)  , nBinMVA, xbins); histo_EWK_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_QCD_CMS_MVALepEffMBoundingAvg   = new TH1D( Form("histo_QCD_%sAvg",effMName)  , Form("histo_QCD_%sAvg",effMName)  , nBinMVA, xbins); histo_QCD_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffMBoundingAvg    = new TH1D( Form("histo_WZ_%sAvg",effMName)  , Form("histo_WZ_%sAvg",effMName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffMBoundingAvg    = new TH1D( Form("histo_ZZ_%sAvg",effMName)  , Form("histo_ZZ_%sAvg",effMName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg   = new TH1D( Form("histo_VVV_%sAvg",effMName)  , Form("histo_VVV_%sAvg",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WS_CMS_MVALepEffMBoundingAvg    = new TH1D( Form("histo_WS_%sAvg",effMName)  , Form("histo_WS_%sAvg",effMName)  , nBinMVA, xbins); histo_WS_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffMBoundingAvg    = new TH1D( Form("histo_WG_%sAvg",effMName)  , Form("histo_WG_%sAvg",effMName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_DPS_CMS_MVALepEffMBoundingAvg   = new TH1D( Form("histo_DPS_%sAvg",effMName)  , Form("histo_DPS_%sAvg",effMName)  , nBinMVA, xbins); histo_DPS_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingAvg[nSigModels];

  TH1D* histo_EWK_CMS_MVALepEffEBoundingUp    = new TH1D( Form("histo_EWK_%sUp",effEName)  , Form("histo_EWK_%sUp",effEName)  , nBinMVA, xbins); histo_EWK_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_EWK_CMS_MVALepEffEBoundingDown  = new TH1D( Form("histo_EWK_%sDown",effEName), Form("histo_EWK_%sDown",effEName), nBinMVA, xbins); histo_EWK_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_QCD_CMS_MVALepEffEBoundingUp    = new TH1D( Form("histo_QCD_%sUp",effEName)  , Form("histo_QCD_%sUp",effEName)  , nBinMVA, xbins); histo_QCD_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_QCD_CMS_MVALepEffEBoundingDown  = new TH1D( Form("histo_QCD_%sDown",effEName), Form("histo_QCD_%sDown",effEName), nBinMVA, xbins); histo_QCD_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingUp     = new TH1D( Form("histo_WZ_%sUp",effEName)  , Form("histo_WZ_%sUp",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingDown   = new TH1D( Form("histo_WZ_%sDown",effEName), Form("histo_WZ_%sDown",effEName), nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingUp     = new TH1D( Form("histo_ZZ_%sUp",effEName)  , Form("histo_ZZ_%sUp",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingDown   = new TH1D( Form("histo_ZZ_%sDown",effEName), Form("histo_ZZ_%sDown",effEName), nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingUp    = new TH1D( Form("histo_VVV_%sUp",effEName)  , Form("histo_VVV_%sUp",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingDown  = new TH1D( Form("histo_VVV_%sDown",effEName), Form("histo_VVV_%sDown",effEName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVALepEffEBoundingUp     = new TH1D( Form("histo_WS_%sUp",effEName)  , Form("histo_WS_%sUp",effEName)  , nBinMVA, xbins); histo_WS_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVALepEffEBoundingDown   = new TH1D( Form("histo_WS_%sDown",effEName), Form("histo_WS_%sDown",effEName), nBinMVA, xbins); histo_WS_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffEBoundingUp     = new TH1D( Form("histo_WG_%sUp",effEName)  , Form("histo_WG_%sUp",effEName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffEBoundingDown   = new TH1D( Form("histo_WG_%sDown",effEName), Form("histo_WG_%sDown",effEName), nBinMVA, xbins); histo_WG_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_DPS_CMS_MVALepEffEBoundingUp    = new TH1D( Form("histo_DPS_%sUp",effEName)  , Form("histo_DPS_%sUp",effEName)  , nBinMVA, xbins); histo_DPS_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_DPS_CMS_MVALepEffEBoundingDown  = new TH1D( Form("histo_DPS_%sDown",effEName), Form("histo_DPS_%sDown",effEName), nBinMVA, xbins); histo_DPS_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingUp[nSigModels];  
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingDown[nSigModels];

  TH1D* histo_EWK_CMS_MVALepEffEBoundingAvg   = new TH1D( Form("histo_EWK_%sAvg",effEName)  , Form("histo_EWK_%sAvg",effEName)  , nBinMVA, xbins); histo_EWK_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_QCD_CMS_MVALepEffEBoundingAvg   = new TH1D( Form("histo_QCD_%sAvg",effEName)  , Form("histo_QCD_%sAvg",effEName)  , nBinMVA, xbins); histo_QCD_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WZ_CMS_MVALepEffEBoundingAvg    = new TH1D( Form("histo_WZ_%sAvg",effEName)  , Form("histo_WZ_%sAvg",effEName)  , nBinMVA, xbins); histo_WZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVALepEffEBoundingAvg    = new TH1D( Form("histo_ZZ_%sAvg",effEName)  , Form("histo_ZZ_%sAvg",effEName)  , nBinMVA, xbins); histo_ZZ_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg   = new TH1D( Form("histo_VVV_%sAvg",effEName)  , Form("histo_VVV_%sAvg",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WS_CMS_MVALepEffEBoundingAvg    = new TH1D( Form("histo_WS_%sAvg",effEName)  , Form("histo_WS_%sAvg",effEName)  , nBinMVA, xbins); histo_WS_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffEBoundingAvg    = new TH1D( Form("histo_WG_%sAvg",effEName)  , Form("histo_WG_%sAvg",effEName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_DPS_CMS_MVALepEffEBoundingAvg   = new TH1D( Form("histo_DPS_%sAvg",effEName)  , Form("histo_DPS_%sAvg",effEName)  , nBinMVA, xbins); histo_DPS_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingAvg[nSigModels];

  TH1D* histo_EWK_CMS_MVAMETBoundingUp        = new TH1D( Form("histo_EWK_%sUp",metName)  , Form("histo_EWK_%sUp",metName)  , nBinMVA, xbins); histo_EWK_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_EWK_CMS_MVAMETBoundingDown      = new TH1D( Form("histo_EWK_%sDown",metName), Form("histo_EWK_%sDown",metName), nBinMVA, xbins); histo_EWK_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_QCD_CMS_MVAMETBoundingUp        = new TH1D( Form("histo_QCD_%sUp",metName)  , Form("histo_QCD_%sUp",metName)  , nBinMVA, xbins); histo_QCD_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_QCD_CMS_MVAMETBoundingDown      = new TH1D( Form("histo_QCD_%sDown",metName), Form("histo_QCD_%sDown",metName), nBinMVA, xbins); histo_QCD_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingUp         = new TH1D( Form("histo_WZ_%sUp",metName)  , Form("histo_WZ_%sUp",metName)  , nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAMETBoundingDown       = new TH1D( Form("histo_WZ_%sDown",metName), Form("histo_WZ_%sDown",metName), nBinMVA, xbins); histo_WZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingUp         = new TH1D( Form("histo_ZZ_%sUp",metName)  , Form("histo_ZZ_%sUp",metName)  , nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAMETBoundingDown       = new TH1D( Form("histo_ZZ_%sDown",metName), Form("histo_ZZ_%sDown",metName), nBinMVA, xbins); histo_ZZ_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingUp        = new TH1D( Form("histo_VVV_%sUp",metName)  , Form("histo_VVV_%sUp",metName)  , nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingDown      = new TH1D( Form("histo_VVV_%sDown",metName), Form("histo_VVV_%sDown",metName), nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVAMETBoundingUp         = new TH1D( Form("histo_WS_%sUp",metName)  , Form("histo_WS_%sUp",metName)  , nBinMVA, xbins); histo_WS_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAMETBoundingDown       = new TH1D( Form("histo_WS_%sDown",metName), Form("histo_WS_%sDown",metName), nBinMVA, xbins); histo_WS_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAMETBoundingUp         = new TH1D( Form("histo_WG_%sUp",metName)  , Form("histo_WG_%sUp",metName)  , nBinMVA, xbins); histo_WG_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAMETBoundingDown       = new TH1D( Form("histo_WG_%sDown",metName), Form("histo_WG_%sDown",metName), nBinMVA, xbins); histo_WG_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_DPS_CMS_MVAMETBoundingUp        = new TH1D( Form("histo_DPS_%sUp",metName)  , Form("histo_DPS_%sUp",metName)  , nBinMVA, xbins); histo_DPS_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_DPS_CMS_MVAMETBoundingDown      = new TH1D( Form("histo_DPS_%sDown",metName), Form("histo_DPS_%sDown",metName), nBinMVA, xbins); histo_DPS_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingUp[nSigModels];  
  TH1D* histo_Higgs_CMS_MVAMETBoundingDown[nSigModels];

  TH1D* histo_EWK_CMS_MVAJESBoundingUp        = new TH1D( Form("histo_EWK_%sUp",jesName)  , Form("histo_EWK_%sUp",jesName)  , nBinMVA, xbins); histo_EWK_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_EWK_CMS_MVAJESBoundingDown      = new TH1D( Form("histo_EWK_%sDown",jesName), Form("histo_EWK_%sDown",jesName), nBinMVA, xbins); histo_EWK_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_QCD_CMS_MVAJESBoundingUp        = new TH1D( Form("histo_QCD_%sUp",jesName)  , Form("histo_QCD_%sUp",jesName)  , nBinMVA, xbins); histo_QCD_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_QCD_CMS_MVAJESBoundingDown      = new TH1D( Form("histo_QCD_%sDown",jesName), Form("histo_QCD_%sDown",jesName), nBinMVA, xbins); histo_QCD_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingUp         = new TH1D( Form("histo_WZ_%sUp",jesName)  , Form("histo_WZ_%sUp",jesName)  , nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAJESBoundingDown       = new TH1D( Form("histo_WZ_%sDown",jesName), Form("histo_WZ_%sDown",jesName), nBinMVA, xbins); histo_WZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingUp         = new TH1D( Form("histo_ZZ_%sUp",jesName)  , Form("histo_ZZ_%sUp",jesName)  , nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJESBoundingDown       = new TH1D( Form("histo_ZZ_%sDown",jesName), Form("histo_ZZ_%sDown",jesName), nBinMVA, xbins); histo_ZZ_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingUp        = new TH1D( Form("histo_VVV_%sUp",jesName)  , Form("histo_VVV_%sUp",jesName)  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown      = new TH1D( Form("histo_VVV_%sDown",jesName), Form("histo_VVV_%sDown",jesName), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVAJESBoundingUp         = new TH1D( Form("histo_WS_%sUp",jesName)  , Form("histo_WS_%sUp",jesName)  , nBinMVA, xbins); histo_WS_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAJESBoundingDown       = new TH1D( Form("histo_WS_%sDown",jesName), Form("histo_WS_%sDown",jesName), nBinMVA, xbins); histo_WS_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAJESBoundingUp         = new TH1D( Form("histo_WG_%sUp",jesName)  , Form("histo_WG_%sUp",jesName)  , nBinMVA, xbins); histo_WG_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAJESBoundingDown       = new TH1D( Form("histo_WG_%sDown",jesName), Form("histo_WG_%sDown",jesName), nBinMVA, xbins); histo_WG_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_DPS_CMS_MVAJESBoundingUp        = new TH1D( Form("histo_DPS_%sUp",jesName)  , Form("histo_DPS_%sUp",jesName)  , nBinMVA, xbins); histo_DPS_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_DPS_CMS_MVAJESBoundingDown      = new TH1D( Form("histo_DPS_%sDown",jesName), Form("histo_DPS_%sDown",jesName), nBinMVA, xbins); histo_DPS_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJESBoundingUp[nSigModels];  
  TH1D* histo_Higgs_CMS_MVAJESBoundingDown[nSigModels];

  TH1D* histo_EWK_CMS_MVAJERBoundingUp        = new TH1D( Form("histo_EWK_%sUp",jerName)  , Form("histo_EWK_%sUp",jerName)  , nBinMVA, xbins); histo_EWK_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_EWK_CMS_MVAJERBoundingDown      = new TH1D( Form("histo_EWK_%sDown",jerName), Form("histo_EWK_%sDown",jerName), nBinMVA, xbins); histo_EWK_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_QCD_CMS_MVAJERBoundingUp        = new TH1D( Form("histo_QCD_%sUp",jerName)  , Form("histo_QCD_%sUp",jerName)  , nBinMVA, xbins); histo_QCD_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_QCD_CMS_MVAJERBoundingDown      = new TH1D( Form("histo_QCD_%sDown",jerName), Form("histo_QCD_%sDown",jerName), nBinMVA, xbins); histo_QCD_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVAJERBoundingUp         = new TH1D( Form("histo_WZ_%sUp",jerName)  , Form("histo_WZ_%sUp",jerName)  , nBinMVA, xbins); histo_WZ_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVAJERBoundingDown       = new TH1D( Form("histo_WZ_%sDown",jerName), Form("histo_WZ_%sDown",jerName), nBinMVA, xbins); histo_WZ_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJERBoundingUp         = new TH1D( Form("histo_ZZ_%sUp",jerName)  , Form("histo_ZZ_%sUp",jerName)  , nBinMVA, xbins); histo_ZZ_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVAJERBoundingDown       = new TH1D( Form("histo_ZZ_%sDown",jerName), Form("histo_ZZ_%sDown",jerName), nBinMVA, xbins); histo_ZZ_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJERBoundingUp        = new TH1D( Form("histo_VVV_%sUp",jerName)  , Form("histo_VVV_%sUp",jerName)  , nBinMVA, xbins); histo_VVV_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJERBoundingDown      = new TH1D( Form("histo_VVV_%sDown",jerName), Form("histo_VVV_%sDown",jerName), nBinMVA, xbins); histo_VVV_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVAJERBoundingUp         = new TH1D( Form("histo_WS_%sUp",jerName)  , Form("histo_WS_%sUp",jerName)  , nBinMVA, xbins); histo_WS_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVAJERBoundingDown       = new TH1D( Form("histo_WS_%sDown",jerName), Form("histo_WS_%sDown",jerName), nBinMVA, xbins); histo_WS_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAJERBoundingUp         = new TH1D( Form("histo_WG_%sUp",jerName)  , Form("histo_WG_%sUp",jerName)  , nBinMVA, xbins); histo_WG_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAJERBoundingDown       = new TH1D( Form("histo_WG_%sDown",jerName), Form("histo_WG_%sDown",jerName), nBinMVA, xbins); histo_WG_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_DPS_CMS_MVAJERBoundingUp        = new TH1D( Form("histo_DPS_%sUp",jerName)  , Form("histo_DPS_%sUp",jerName)  , nBinMVA, xbins); histo_DPS_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_DPS_CMS_MVAJERBoundingDown      = new TH1D( Form("histo_DPS_%sDown",jerName), Form("histo_DPS_%sDown",jerName), nBinMVA, xbins); histo_DPS_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJERBoundingUp[nSigModels];  
  TH1D* histo_Higgs_CMS_MVAJERBoundingDown[nSigModels];

  TH1D* histo_EWK_CMS_MVABTAGBoundingUp        = new TH1D( Form("histo_EWK_%sUp",btagName)  , Form("histo_EWK_%sUp",btagName)  , nBinMVA, xbins); histo_EWK_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_EWK_CMS_MVABTAGBoundingDown      = new TH1D( Form("histo_EWK_%sDown",btagName), Form("histo_EWK_%sDown",btagName), nBinMVA, xbins); histo_EWK_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_QCD_CMS_MVABTAGBoundingUp        = new TH1D( Form("histo_QCD_%sUp",btagName)  , Form("histo_QCD_%sUp",btagName)  , nBinMVA, xbins); histo_QCD_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_QCD_CMS_MVABTAGBoundingDown      = new TH1D( Form("histo_QCD_%sDown",btagName), Form("histo_QCD_%sDown",btagName), nBinMVA, xbins); histo_QCD_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_MVABTAGBoundingUp         = new TH1D( Form("histo_WZ_%sUp",btagName)  , Form("histo_WZ_%sUp",btagName)  , nBinMVA, xbins); histo_WZ_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_MVABTAGBoundingDown       = new TH1D( Form("histo_WZ_%sDown",btagName), Form("histo_WZ_%sDown",btagName), nBinMVA, xbins); histo_WZ_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_MVABTAGBoundingUp         = new TH1D( Form("histo_ZZ_%sUp",btagName)  , Form("histo_ZZ_%sUp",btagName)  , nBinMVA, xbins); histo_ZZ_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_MVABTAGBoundingDown       = new TH1D( Form("histo_ZZ_%sDown",btagName), Form("histo_ZZ_%sDown",btagName), nBinMVA, xbins); histo_ZZ_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVABTAGBoundingUp        = new TH1D( Form("histo_VVV_%sUp",btagName)  , Form("histo_VVV_%sUp",btagName)  , nBinMVA, xbins); histo_VVV_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVABTAGBoundingDown      = new TH1D( Form("histo_VVV_%sDown",btagName), Form("histo_VVV_%sDown",btagName), nBinMVA, xbins); histo_VVV_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_MVABTAGBoundingUp         = new TH1D( Form("histo_WS_%sUp",btagName)  , Form("histo_WS_%sUp",btagName)  , nBinMVA, xbins); histo_WS_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_MVABTAGBoundingDown       = new TH1D( Form("histo_WS_%sDown",btagName), Form("histo_WS_%sDown",btagName), nBinMVA, xbins); histo_WS_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVABTAGBoundingUp         = new TH1D( Form("histo_WG_%sUp",btagName)  , Form("histo_WG_%sUp",btagName)  , nBinMVA, xbins); histo_WG_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVABTAGBoundingDown       = new TH1D( Form("histo_WG_%sDown",btagName), Form("histo_WG_%sDown",btagName), nBinMVA, xbins); histo_WG_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_DPS_CMS_MVABTAGBoundingUp        = new TH1D( Form("histo_DPS_%sUp",btagName)  , Form("histo_DPS_%sUp",btagName)  , nBinMVA, xbins); histo_DPS_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_DPS_CMS_MVABTAGBoundingDown      = new TH1D( Form("histo_DPS_%sDown",btagName), Form("histo_DPS_%sDown",btagName), nBinMVA, xbins); histo_DPS_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVABTAGBoundingUp[nSigModels];  
  TH1D* histo_Higgs_CMS_MVABTAGBoundingDown[nSigModels];

  TH1D* histo_EWK_CMS_PUBoundingUp    	      = new TH1D( Form("histo_EWK_%sUp",puName)  , Form("histo_EWK_%sUp",puName)  , nBinMVA, xbins); histo_EWK_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_EWK_CMS_PUBoundingDown  	      = new TH1D( Form("histo_EWK_%sDown",puName), Form("histo_EWK_%sDown",puName), nBinMVA, xbins); histo_EWK_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_QCD_CMS_PUBoundingUp    	      = new TH1D( Form("histo_QCD_%sUp",puName)  , Form("histo_QCD_%sUp",puName)  , nBinMVA, xbins); histo_QCD_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_QCD_CMS_PUBoundingDown  	      = new TH1D( Form("histo_QCD_%sDown",puName), Form("histo_QCD_%sDown",puName), nBinMVA, xbins); histo_QCD_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WZ_CMS_PUBoundingUp    	      = new TH1D( Form("histo_WZ_%sUp",puName)  , Form("histo_WZ_%sUp",puName)  , nBinMVA, xbins); histo_WZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WZ_CMS_PUBoundingDown  	      = new TH1D( Form("histo_WZ_%sDown",puName), Form("histo_WZ_%sDown",puName), nBinMVA, xbins); histo_WZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingUp    	      = new TH1D( Form("histo_ZZ_%sUp",puName)  , Form("histo_ZZ_%sUp",puName)  , nBinMVA, xbins); histo_ZZ_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ZZ_CMS_PUBoundingDown  	      = new TH1D( Form("histo_ZZ_%sDown",puName), Form("histo_ZZ_%sDown",puName), nBinMVA, xbins); histo_ZZ_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingUp    	      = new TH1D( Form("histo_VVV_%sUp",puName)  , Form("histo_VVV_%sUp",puName)  , nBinMVA, xbins); histo_VVV_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingDown  	      = new TH1D( Form("histo_VVV_%sDown",puName), Form("histo_VVV_%sDown",puName), nBinMVA, xbins); histo_VVV_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WS_CMS_PUBoundingUp    	      = new TH1D( Form("histo_WS_%sUp",puName)  , Form("histo_WS_%sUp",puName)  , nBinMVA, xbins); histo_WS_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WS_CMS_PUBoundingDown  	      = new TH1D( Form("histo_WS_%sDown",puName), Form("histo_WS_%sDown",puName), nBinMVA, xbins); histo_WS_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_PUBoundingUp    	      = new TH1D( Form("histo_WG_%sUp",puName)  , Form("histo_WG_%sUp",puName)  , nBinMVA, xbins); histo_WG_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_PUBoundingDown  	      = new TH1D( Form("histo_WG_%sDown",puName), Form("histo_WG_%sDown",puName), nBinMVA, xbins); histo_WG_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_DPS_CMS_PUBoundingUp    	      = new TH1D( Form("histo_DPS_%sUp",puName)  , Form("histo_DPS_%sUp",puName)  , nBinMVA, xbins); histo_DPS_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_DPS_CMS_PUBoundingDown  	      = new TH1D( Form("histo_DPS_%sDown",puName), Form("histo_DPS_%sDown",puName), nBinMVA, xbins); histo_DPS_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingUp[nSigModels];  
  TH1D* histo_Higgs_CMS_PUBoundingDown[nSigModels];

  TH1D* histo_WS_CMS_WSSFUp    	              = new TH1D( Form("histo_WS_CMS_WSSFUp")  , Form("histo_WS_CMS_WSSFUp")  , nBinMVA, xbins); histo_WS_CMS_WSSFUp  ->Sumw2();
  TH1D* histo_WS_CMS_WSSFDown  	              = new TH1D( Form("histo_WS_CMS_WSSFDown"), Form("histo_WS_CMS_WSSFDown"), nBinMVA, xbins); histo_WS_CMS_WSSFDown->Sumw2();

  TH1D* histo_Fake_CMS_SystM    	      = new TH1D( Form("histo_Fake_CMS_SystM"), Form("histo_Fake_CMS_SystM"), nBinMVA, xbins); histo_Fake_CMS_SystM->Sumw2();
  TH1D* histo_Fake_CMS_SystE  	              = new TH1D( Form("histo_Fake_CMS_SystE"), Form("histo_Fake_CMS_SystE"), nBinMVA, xbins); histo_Fake_CMS_SystE->Sumw2();

  for(int nModel=0; nModel<nSigModels; nModel++) { 
    histo_Higgs_CMS_MVALepEffMBoundingUp[nModel]          = new TH1D( Form("histo_Higgs_%s_%sUp",   signalName_[nModel].Data(), effMName), Form("histo_Higgs_%s_%sUp",  signalName_[nModel].Data(), effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingUp[nModel]  ->Sumw2();
    histo_Higgs_CMS_MVALepEffMBoundingDown[nModel]        = new TH1D( Form("histo_Higgs_%s_%sDown", signalName_[nModel].Data(), effMName), Form("histo_Higgs_%s_%sDown",signalName_[nModel].Data(), effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingDown[nModel]->Sumw2();
    histo_Higgs_CMS_MVALepEffMBoundingAvg [nModel]        = new TH1D( Form("histo_Higgs_%s_%sAvg",  	   signalName_[nModel].Data(), effMName), Form("histo_Higgs_%s_%sAvg" ,	   signalName_[nModel].Data(), effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingAvg[nModel]  ->Sumw2();
    histo_Higgs_CMS_MVALepEffEBoundingUp [nModel]         = new TH1D( Form("histo_Higgs_%s_%sUp",		   signalName_[nModel].Data(), effEName), Form("histo_Higgs_%s_%sUp"  ,	   signalName_[nModel].Data(), effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingUp[nModel]  ->Sumw2();
    histo_Higgs_CMS_MVALepEffEBoundingDown [nModel]       = new TH1D( Form("histo_Higgs_%s_%sDown", 	   signalName_[nModel].Data(), effEName), Form("histo_Higgs_%s_%sDown",	   signalName_[nModel].Data(), effEName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingDown[nModel]->Sumw2();
    histo_Higgs_CMS_MVALepEffEBoundingAvg [nModel]        = new TH1D( Form("histo_Higgs_%s_%sAvg",  	   signalName_[nModel].Data(), effEName), Form("histo_Higgs_%s_%sAvg" ,	   signalName_[nModel].Data(), effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingAvg[nModel]  ->Sumw2();
    histo_Higgs_CMS_MVAMETBoundingUp [nModel]             = new TH1D( Form("histo_Higgs_%s_CMS_scale_metUp"  , signalName_[nModel].Data()), 	  Form("histo_Higgs_%s_CMS_scale_metUp"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingUp[nModel]  ->Sumw2();
    histo_Higgs_CMS_MVAMETBoundingDown [nModel]           = new TH1D( Form("histo_Higgs_%s_CMS_scale_metDown", signalName_[nModel].Data()), 	  Form("histo_Higgs_%s_CMS_scale_metDown", signalName_[nModel].Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingDown[nModel]->Sumw2();
    histo_Higgs_CMS_MVAJESBoundingUp [nModel]             = new TH1D( Form("histo_Higgs_%s_CMS_scale_jUp"	 , signalName_[nModel].Data()), 	  Form("histo_Higgs_%s_CMS_scale_jUp"    , signalName_[nModel].Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingUp[nModel]  ->Sumw2();
    histo_Higgs_CMS_MVAJESBoundingDown [nModel]           = new TH1D( Form("histo_Higgs_%s_CMS_scale_jDown"  , signalName_[nModel].Data()), 	  Form("histo_Higgs_%s_CMS_scale_jDown"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingDown[nModel]->Sumw2();
    histo_Higgs_CMS_MVAJERBoundingUp [nModel]             = new TH1D( Form("histo_Higgs_%s_CMS_jerUp"	 , signalName_[nModel].Data()), 	  Form("histo_Higgs_%s_CMS_jerUp"    , signalName_[nModel].Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAJERBoundingUp[nModel]  ->Sumw2();
    histo_Higgs_CMS_MVAJERBoundingDown [nModel]           = new TH1D( Form("histo_Higgs_%s_CMS_jerDown"  , signalName_[nModel].Data()), 	  Form("histo_Higgs_%s_CMS_jerDown"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAJERBoundingDown[nModel]->Sumw2();
    histo_Higgs_CMS_MVABTAGBoundingUp [nModel]            = new TH1D( Form("histo_Higgs_%s_CMS_eff_b_2016Up"	, signalName_[nModel].Data()),  Form("histo_Higgs_%s_CMS_eff_b_2016Up"    , signalName_[nModel].Data()), nBinMVA, xbins); histo_Higgs_CMS_MVABTAGBoundingUp[nModel]  ->Sumw2();
    histo_Higgs_CMS_MVABTAGBoundingDown [nModel]          = new TH1D( Form("histo_Higgs_%s_CMS_eff_b_2016Down"  , signalName_[nModel].Data()), 	 Form("histo_Higgs_%s_CMS_eff_b_2016Down"  , signalName_[nModel].Data()), nBinMVA, xbins); histo_Higgs_CMS_MVABTAGBoundingDown[nModel]->Sumw2();
    histo_Higgs_CMS_PUBoundingUp [nModel]                 = new TH1D( Form("histo_Higgs_%s_CMS_puUp"	 , signalName_[nModel].Data()), 	  Form("histo_Higgs_%s_CMS_puUp"	   , signalName_[nModel].Data()), nBinMVA, xbins); histo_Higgs_CMS_PUBoundingUp[nModel]  ->Sumw2();
    histo_Higgs_CMS_PUBoundingDown [nModel]               = new TH1D( Form("histo_Higgs_%s_CMS_puDown"	 , signalName_[nModel].Data()), 	  Form("histo_Higgs_%s_CMS_puDown"	   , signalName_[nModel].Data()), nBinMVA, xbins); histo_Higgs_CMS_PUBoundingDown[nModel]->Sumw2();
  }

  const int numberCuts = 12;
  TH1D* hDWWLL[7];
  hDWWLL[0] = new TH1D("hDWWLL0", "hDWWLL0", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[1] = new TH1D("hDWWLL1", "hDWWLL1", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[2] = new TH1D("hDWWLL2", "hDWWLL2", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[3] = new TH1D("hDWWLL3", "hDWWLL3", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[4] = new TH1D("hDWWLL4", "hDWWLL4", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[5] = new TH1D("hDWWLL5", "hDWWLL5", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[6] = new TH1D("hDWWLL6", "hDWWLL6", numberCuts, -0.5, numberCuts-0.5);
  TString cutName[numberCuts] = {"#lep, sign, flavor","        lep pT cut","  #jets with pT>30","         btag-veto","          tau veto",
                                 "           mll cut","      loose Z veto","            Z veto","           met cut","           mjj cut",
				 "    deltaetajj cut","    zeppenfeld cut" };
  TString signalName="";
  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());

    int sigModel = (infilecatv[ifile]==11) ? signalIndex_[ifile] : -1;
    if(sigModel>=0) signalName=signalName_[sigModel];

    TFile *the_input_file = TFile::Open(infilenamev[ifile].Data());
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
    TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");
    TTree *the_SelBit_tree= (TTree*)the_input_file->FindObjectAny("SelBit_tree");

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareJets eventJets;
    if(usePUPPI) eventJets.setBranchAddresses(the_input_tree,"puppi");
    else         eventJets.setBranchAddresses(the_input_tree);

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
 
    TNamed *triggerNames = (TNamed*)the_input_file->FindObjectAny("triggerNames");
    char **tokens;
    size_t numtokens;
    tokens = strsplit(triggerNames->GetTitle(), ",", &numtokens);
    if(ifile == 0){
      for (int i = 0; i < (int)numtokens; i++) {
        printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
      }
    }
    else {
    }
    
    bool errorMsgQCDscale = false;
    unsigned int selBit_= 0;
    the_SelBit_tree->SetBranchAddress("selBit", &selBit_);
    hDWWLL[0]->Scale(0.0);hDWWLL[1]->Scale(0.0);hDWWLL[2]->Scale(0.0);hDWWLL[3]->Scale(0.0);hDWWLL[4]->Scale(0.0);hDWWLL[5]->Scale(0.0);hDWWLL[6]->Scale(0.0);
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    if(the_input_tree->GetEntries() != the_SelBit_tree->GetEntries()) {printf("BIG SKIMMING FAILURE\n"); return;}
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_SelBit_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      if((selBit_ & 0x1<<whichSkim1) == 0 && (selBit_ & 0x1<<whichSkim2) == 0) continue;

      the_input_tree->GetEntry(i);
 
      int initPDFTag = 0;
      if((*eventMonteCarlo.pdfRwgt).size() == 0) initPDFTag = -1;

      int sumEvol[6] = {-1, -1, -1, -1, -1, -1};
      Bool_t passFilterSig[numberCuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      Bool_t passFilterCR1[numberCuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      Bool_t passFilterCR2[numberCuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      Bool_t passPresel = kFALSE;
      Bool_t passSignalRegion = kFALSE;
      Bool_t passControlRegionTop      = kFALSE;
      Bool_t passControlRegionWZ       = kFALSE;
      Bool_t passControlRegionDi       = kFALSE;
      Bool_t passControlRegionSS2j     = kFALSE;
      Bool_t passControlRegionOS2j     = kFALSE;
      Bool_t passControlRegionZLL      = kFALSE;
      Bool_t passLooseControlRegionTop = kFALSE;
      Bool_t passLooseControlRegionWZ  = kFALSE;
      Bool_t passMjjLoose = kFALSE;
      Bool_t passOS = kFALSE;

      if(infilecatv[ifile] != 999) {
        for (int nt = 0; nt <(int)numtokens; nt++) {
          if((*eventTrigger.triggerFired)[nt] == 0) continue;
          if((strcmp(tokens[nt],Form("HLT_Ele25_eta2p1_WPTight_Gsf_v%s",triggerSuffix.Data()))  == 0) ||
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
           ) passPresel = kTRUE;
        }
      } else { passPresel = kTRUE;}

      //acess lepton info
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0; vector<int> idLepLoose;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                                   {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake )     {idTight.push_back(0); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP)    {idSoft.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepLoose)== BareLeptons::LepLoose && 
                ((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt() > 10.0)                                  {idLepLoose.push_back(nlep);}
      }

      if(idLep.size()!=idTight.size()) {assert(1); return;}

      if(idLep.size() == 2 || idLep.size() == 3) passPresel = passPresel && kTRUE;
      else                                       passPresel = passPresel && kFALSE;

      if(passPresel == kFALSE) continue; // two or three leptons in the event

      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 25 && 
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20) passPresel = passPresel && kTRUE;
      else                                                           passPresel = passPresel && kFALSE;

      if(passPresel == kFALSE) continue; // ptl1/l2 > 25/20

      // typeSel = 0(m+m+), 1(e+e+), 2(e+m+/m+e+) 3(m-m-), 4(e-e-), 5(e-m-/m-e-) --> ++/-- for finalVar == 2 / one category for  finalVar == 3
      // typePair = 0(mm), 1(ee), 2(em)
      int typeSel = -1; int typePair = -1;
      if(idTight.size() >= 2){
        if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) {typeSel = 0; typePair = 0;}
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) {typeSel = 1; typePair = 1;}
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])!=TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]))         {typeSel = 2; typePair = 2;}
        else {assert(1); return;}
	if((int)(*eventLeptons.pdgId)[idLep[0]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) > 0) typeSel = typeSel + 3;
      }
      else {assert(1); printf("Not possible %d %d %d\n",(int)idLep.size(),(int)idTight.size(),goodIsTight); return;}                                                                                                     ;

      if     (finalVar == 2) {if(typeSel <= 2) typeSel = 0; else typeSel = 1;}
      else if(finalVar == 3) {typeSel = 0;}

     //signQ = 0(opposite sign), +/-2(same sign)
      int signQ = -1;
      signQ = (int)(*eventLeptons.pdgId)[idLep[0]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) + (int)(*eventLeptons.pdgId)[idLep[1]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]);

      //#lep, same sign or opposite sign, flavor combination choice
      if(idTight.size() == 2 && signQ != 0){
	passFilterSig[0] = kTRUE;
	passFilterCR1[0] = kTRUE;
      }
      if(idTight.size() == 2 && signQ == 0 && typePair == 2)  passOS = kTRUE;
      if(idTight.size() == 3 )                                passFilterCR2[0] = kTRUE; 

      //lepton pT cut
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 25 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20){
	passFilterSig[1] = kTRUE;
	passFilterCR1[1] = kTRUE;
	if(idTight.size() == 3 && ((TLorentzVector*)(*eventLeptons.p4)[idLep[2]])->Pt() > 10) passFilterCR2[1] = kTRUE;
      }

      vector<int> idB,idC;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if     (TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 5 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) idB.push_back(ngen0);
        else if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 4 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) idC.push_back(ngen0);
      }

      //access jet info
      vector<int> idJet,idJesUp,idJesDown,idJerUp,idJerDown;
      double bDiscrMax = 0.0;
      double total_bjet_probLOOSE[2] = {1,1};double total_bjet_probLOOSEUP[2] = {1,1};double total_bjet_probLOOSEDOWN[2] = {1,1};
      double total_bjet_probTIGHT[2] = {1,1};double total_bjet_probTIGHTUP[2] = {1,1};double total_bjet_probTIGHTDOWN[2] = {1,1};
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 20) continue;

        //if(((int)(*eventJets.selBits)[nj] & BareJets::JetTight)!= BareJets::JetTight) continue;

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.4) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

	if(infilecatv[ifile] != 0){ // btagging study
          BTagEntry::JetFlavor jetFlavor = BTagEntry::FLAV_UDSG;
	  for(unsigned int ng=0; ng<idB.size(); ng++) {
	    if (((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[idB[ng]])) < 0.4) {jetFlavor = BTagEntry::FLAV_B; break;}
	  }
	  if(jetFlavor == BTagEntry::FLAV_UDSG){
	    for(unsigned int ng=0; ng<idC.size(); ng++) {
	      if (((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[idC[ng]])) < 0.4) {jetFlavor = BTagEntry::FLAV_C; break;}
	    }
          }

          int nJPt = 0;
	  if     (((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) nJPt = 0;
	  else if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 40) nJPt = 1;
	  else if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 60) nJPt = 2;
	  else if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 80) nJPt = 3;
          else                                                       nJPt = 4;
          int nJEta = 0;
	  if     (TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()) < 0.5) nJEta = 0;
	  else if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()) < 1.0) nJEta = 1;
	  else if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()) < 1.5) nJEta = 2;
	  else if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()) < 2.0) nJEta = 3;
          else                                                                     nJEta = 4;
          denBTagging[nJPt][nJEta][jetFlavor]++;
          if((float)(*eventJets.bDiscr)[nj] >= bTagCuts[0]) numBTaggingLOOSE[nJPt][nJEta][jetFlavor]++;
          if((float)(*eventJets.bDiscr)[nj] >= bTagCuts[1]) numBTaggingTIGHT[nJPt][nJEta][jetFlavor]++;

          double bjet_SFLOOSE = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFLOOSE = btagReaderLLOOSE.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
	  else                                  bjet_SFLOOSE = btagReaderBCLOOSE.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
          if(bjet_SFLOOSE == 0) bjet_SFLOOSE = 1;
          double bjet_SFLOOSEUP = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFLOOSEUP = btagReaderLLOOSEUP.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
	  else                                  bjet_SFLOOSEUP = btagReaderBCLOOSEUP.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
          if(bjet_SFLOOSEUP == 0) bjet_SFLOOSEUP = 1;
          double bjet_SFLOOSEDOWN = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFLOOSEDOWN = btagReaderLLOOSEDOWN.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
	  else                                  bjet_SFLOOSEDOWN = btagReaderBCLOOSEDOWN.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
          if(bjet_SFLOOSEDOWN == 0) bjet_SFLOOSEDOWN = 1;

	  if((float)(*eventJets.bDiscr)[nj] >= bTagCuts[0]){
	    total_bjet_probLOOSE[0] = total_bjet_probLOOSE[0] * jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor];
	    total_bjet_probLOOSE[1] = total_bjet_probLOOSE[1] * jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor] * bjet_SFLOOSE;
	    total_bjet_probLOOSEUP[0] = total_bjet_probLOOSEUP[0] * jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor];
	    total_bjet_probLOOSEUP[1] = total_bjet_probLOOSEUP[1] * jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor] * bjet_SFLOOSEUP;
	    total_bjet_probLOOSEDOWN[0] = total_bjet_probLOOSEDOWN[0] * jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor];
	    total_bjet_probLOOSEDOWN[1] = total_bjet_probLOOSEDOWN[1] * jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor] * bjet_SFLOOSEDOWN;
	  } else {
	    total_bjet_probLOOSE[0] = total_bjet_probLOOSE[0] * (1.0 - jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor]);
	    total_bjet_probLOOSE[1] = total_bjet_probLOOSE[1] * (1.0 - jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor] * bjet_SFLOOSE);
	    total_bjet_probLOOSEUP[0] = total_bjet_probLOOSEUP[0] * (1.0 - jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor]);
	    total_bjet_probLOOSEUP[1] = total_bjet_probLOOSEUP[1] * (1.0 - jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor] * bjet_SFLOOSEUP);
	    total_bjet_probLOOSEDOWN[0] = total_bjet_probLOOSEDOWN[0] * (1.0 - jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor]);
	    total_bjet_probLOOSEDOWN[1] = total_bjet_probLOOSEDOWN[1] * (1.0 - jetEpsBtagLOOSE[nJPt][nJEta][jetFlavor] * bjet_SFLOOSEDOWN);
	  }

          double bjet_SFTIGHT = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFTIGHT = btagReaderLTIGHT.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
	  else                                  bjet_SFTIGHT = btagReaderBCTIGHT.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
          if(bjet_SFTIGHT == 0) bjet_SFTIGHT = 1;
          double bjet_SFTIGHTUP = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFTIGHTUP = btagReaderLTIGHTUP.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
	  else                                  bjet_SFTIGHTUP = btagReaderBCTIGHTUP.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
          if(bjet_SFTIGHTUP == 0) bjet_SFTIGHTUP = 1;
          double bjet_SFTIGHTDOWN = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFTIGHTDOWN = btagReaderLTIGHTDOWN.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
	  else                                  bjet_SFTIGHTDOWN = btagReaderBCTIGHTDOWN.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),((TLorentzVector*)(*eventJets.p4)[nj])->Pt());
          if(bjet_SFTIGHTDOWN == 0) bjet_SFTIGHTDOWN = 1;

	  if((float)(*eventJets.bDiscr)[nj] >= bTagCuts[1]){
	    total_bjet_probTIGHT[0] = total_bjet_probTIGHT[0] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor];
	    total_bjet_probTIGHT[1] = total_bjet_probTIGHT[1] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHT;
	    total_bjet_probTIGHTUP[0] = total_bjet_probTIGHTUP[0] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor];
	    total_bjet_probTIGHTUP[1] = total_bjet_probTIGHTUP[1] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHTUP;
	    total_bjet_probTIGHTDOWN[0] = total_bjet_probTIGHTDOWN[0] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor];
	    total_bjet_probTIGHTDOWN[1] = total_bjet_probTIGHTDOWN[1] * jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHTDOWN;
	  } else {
	    total_bjet_probTIGHT[0] = total_bjet_probTIGHT[0] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor]);
	    total_bjet_probTIGHT[1] = total_bjet_probTIGHT[1] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHT);
	    total_bjet_probTIGHTUP[0] = total_bjet_probTIGHTUP[0] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor]);
	    total_bjet_probTIGHTUP[1] = total_bjet_probTIGHTUP[1] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHTUP);
	    total_bjet_probTIGHTDOWN[0] = total_bjet_probTIGHTDOWN[0] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor]);
	    total_bjet_probTIGHTDOWN[1] = total_bjet_probTIGHTDOWN[1] * (1.0 - jetEpsBtagTIGHT[nJPt][nJEta][jetFlavor] * bjet_SFTIGHTDOWN);
	  }

          double sf_jer[2] = {1,1};
          if((float)(*eventJets.ptResUncCentral)[nj] > 0){
	    sf_jer[0] = (float)(*eventJets.ptResUncUp)[nj]   / (float)(*eventJets.ptResUncCentral)[nj];
	    sf_jer[1] = (float)(*eventJets.ptResUncDown)[nj] / (float)(*eventJets.ptResUncCentral)[nj];
          }
          if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*(1+(float)(*eventJets.unc)[nj]) > 30) idJesUp.push_back(nj);
          if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*(1-(float)(*eventJets.unc)[nj]) > 30) idJesDown.push_back(nj);
          if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*sf_jer[0] > 30) idJerUp.push_back(nj);
          if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()*sf_jer[1] > 30) idJerDown.push_back(nj);

        }

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 20) { 
	   if ((float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];
        }

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt()      > 30) {idJet.push_back(nj);}

      }

      //if(TMath::Abs(total_bjet_probLOOSE[1]/total_bjet_probLOOSE[0]-1.0) > 0.05) printf("total_bjet_probLOOSE large correction: %f, bDiscrMax = %f\n",total_bjet_probLOOSE[1]/total_bjet_probLOOSE[0],bDiscrMax);
      //if(TMath::Abs(total_bjet_probTIGHT[1]/total_bjet_probTIGHT[0]-1.0) > 0.05) printf("total_bjet_probTIGHT large correction: %f, bDiscrMax = %f\n",total_bjet_probTIGHT[1]/total_bjet_probTIGHT[0],bDiscrMax);

      //#Jet with pT > 30 GeV
      if(idJet.size() >= 2){
	passFilterSig[2] = kTRUE;
        passFilterCR1[2] = kTRUE;
        passFilterCR2[2] = kTRUE;
      }

      int numberGoodTaus = 0;
      for(int ntau=0; ntau<eventTaus.p4->GetEntriesFast(); ntau++) {
        if(((TLorentzVector*)(*eventTaus.p4)[ntau])->Pt() <= 18.0 ||
           TMath::Abs(((TLorentzVector*)(*eventTaus.p4)[ntau])->Eta()) >= 2.3) continue;
        bool isElMu = false;
        for(unsigned nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventTaus.p4)[ntau])) < 0.3) {
            isElMu = true;
            break;
          }
        }
        if(isElMu == false &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::TauDecayModeFinding	 ) == BareTaus::TauDecayModeFinding &&
           ((int)(*eventTaus.selBits)[ntau] & BareTaus::TauDecayModeFindingNewDMs) == BareTaus::TauDecayModeFindingNewDMs &&
           (double)(*eventTaus.iso)[ntau] < 5.0){
          numberGoodTaus++;
        }
      }

      if(bDiscrMax < bTagCuts[0] && idSoft.size() == 0){
	passFilterSig[3] = kTRUE;
      } else{
        passFilterCR1[3] = kTRUE;
      }
      if(bDiscrMax < bTagCuts[1] && idSoft.size() == 0){
        passFilterCR2[3] = kTRUE;
      }

      //tau veto
      if(numberGoodTaus == 0){
	passFilterSig[4] = kTRUE;
        passFilterCR1[4] = kTRUE;
      }
      passFilterCR2[4] = kTRUE;

      //Mll cut
      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) ));
      if(dilep.M() > 20){
	passFilterSig[5] = kTRUE;
	passFilterCR1[5] = kTRUE;
      }
      passFilterCR2[5] = kTRUE;


      double minMassLooseZ = 999.0;
      for(unsigned nl0=0; nl0<idLep.size(); nl0++){
        for(unsigned nl1=0; nl1<idLepLoose.size(); nl1++){

	  if((int)(*eventLeptons.pdgId)[idLep[nl0]] * (int)(*eventLeptons.pdgId)[idLepLoose[nl1]] > 0) continue;
          TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLepLoose[nl1])) ) ));

	  if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLepLoose[nl1]]) &&
	     TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassLooseZ-91.1876)) {
	     minMassLooseZ = dilepAux.M();
	  }
        }
      }

      // Loose Z veto
      if(TMath::Abs(minMassLooseZ-91.1876) > 15 && idLepLoose.size() == 0){
	passFilterSig[6] = kTRUE;
        passFilterCR1[6] = kTRUE;
      }
      passFilterCR2[6] = kTRUE;
      
      //Z veto
      if((typePair == 1 && TMath::Abs(dilep.M()-91.1876) > 15.0) || typePair != 1){
	passFilterSig[7] = kTRUE;
	passFilterCR1[7] = kTRUE;
      }

      double minMassZ = 999.0;
      if(idTight.size() == 3){
	for(unsigned nl0=0; nl0<idLep.size()-1; nl0++){
          for(unsigned nl1=nl0+1; nl1<idLep.size(); nl1++){

	    if((int)(*eventLeptons.pdgId)[idLep[nl0]] + (int)(*eventLeptons.pdgId)[idLep[nl1]] != 0) continue; // OSSF pairs only
            TLorentzVector dilepAux(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[nl1])) ) ));

	    if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl1]]) &&
	       TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassZ-91.1876)) {
	       minMassZ = dilepAux.M();
	    }
          }
	}
      }
      if(TMath::Abs(minMassZ-91.1876) < 15.0) passFilterCR2[7] = kTRUE;

      //access met info
      TLorentzVector theMET;
      if(usePUPPI){
        theMET.SetPx((double)eventMet.metPuppi->Px());
        theMET.SetPy((double)eventMet.metPuppi->Py());
      } else {
        theMET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px());
        theMET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py());
      }

      //met cut
      double metCut = 40.;
      if(typePair == 1) metCut = 40.;
      if(theMET.Pt() > metCut){
        passFilterSig[8] = kTRUE;
        passFilterCR1[8] = kTRUE;
        passFilterCR2[8] = kTRUE;
      }

      TLorentzVector dijet; double deltaEtaJJ = 0.0;double theLeptonZ[2] = {0.999,0.999};
      if(idJet.size() >= 2){
	//Mjj cut
	dijet = ( *(TLorentzVector*)(*eventJets.p4)[idJet[0]] ) + ( *(TLorentzVector*)(*eventJets.p4)[idJet[1]] );
	if(dijet.M() > mjjCut){
          passFilterSig[9] = kTRUE;
          passFilterCR1[9] = kTRUE;
          passFilterCR2[9] = kTRUE;
	}
	if(dijet.M() > 100){
          passMjjLoose = kTRUE;
	}

	//deltaEta jj cut
	deltaEtaJJ = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta());
	if(deltaEtaJJ > 2.5){
          passFilterSig[10] = kTRUE;
          passFilterCR1[10] = kTRUE;
          passFilterCR2[10] = kTRUE;
	}
	
	theLeptonZ[0] = TMath::Min(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()-(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()+((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta())/2.)/deltaEtaJJ,0.999);
	theLeptonZ[1] = TMath::Min(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()-(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()+((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta())/2.)/deltaEtaJJ,0.999);

        passFilterSig[11] = TMath::Max(theLeptonZ[0],theLeptonZ[1]) < 0.75;
        passFilterCR1[11] = TMath::Max(theLeptonZ[0],theLeptonZ[1]) < 0.75;
        passFilterCR2[11] = kTRUE;
      }

      //                           #lep, sign, flavor  lep pT cut          #jets with pT>30    btag-veto           tau veto            mll cut             loose Z veto        Z veto              met cut             mjj cut             deltaetajj cut       zeppenfeld cut
      passSignalRegion           = passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] && passFilterSig[9] && passFilterSig[10] && passFilterSig[11];
      passControlRegionTop       = passFilterCR1[0] && passFilterCR1[1] && passFilterCR1[2] && passFilterCR1[3] && passFilterCR1[4] && passFilterCR1[5] && passFilterCR1[6] && passFilterCR1[7] && passFilterCR1[8] && passFilterCR1[9] && passFilterCR1[10] && passFilterCR1[11];
      passControlRegionWZ        = passFilterCR2[0] && passFilterCR2[1] && passFilterCR2[2] && passFilterCR2[3] && passFilterCR2[4] && passFilterCR2[5] && passFilterCR2[6] && passFilterCR2[7] && passFilterCR2[8] && passFilterCR2[9] && passFilterCR2[10] && passFilterCR2[11];
      passControlRegionDi        = passOS           && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] && passFilterSig[9] && passFilterSig[10] && passFilterSig[11];
      passControlRegionSS2j      = passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] &&                     passFilterSig[8] && passMjjLoose     &&                      TMath::Max(theLeptonZ[0],theLeptonZ[1]) < 0.999;
      passControlRegionOS2j      = passOS           && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] &&                     passFilterSig[8] && passMjjLoose     &&                      TMath::Max(theLeptonZ[0],theLeptonZ[1]) < 0.999;
      passControlRegionZLL       = passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] &&                     passFilterSig[6] &&!passFilterSig[7] && passFilterSig[8] && passMjjLoose;
      passLooseControlRegionTop  = passFilterCR1[0] && passFilterCR1[1] && passFilterCR1[2] && passFilterCR1[3] && passFilterCR1[4] && passFilterCR1[5] && passFilterCR1[6] && passFilterCR1[7] && passFilterCR1[8]                     && passFilterCR1[10];
      passLooseControlRegionWZ   = passFilterCR2[0] && passFilterCR2[1] && passFilterCR2[2] && passFilterCR2[3] && passFilterCR2[4] && passFilterCR2[5] && passFilterCR2[6] && passFilterCR2[7] && passFilterCR2[8]                     && passFilterCR2[10];

      bool passNMinusOne[8] = {
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] &&                     passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] && passFilterSig[9] && passFilterSig[10] && passFilterSig[11], // btag veto
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] &&                     passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] && passFilterSig[9] && passFilterSig[10] && passFilterSig[11], // tau veto
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] &&                     passFilterSig[6] &&                     passFilterSig[8] && passFilterSig[9] && passFilterSig[10] && passFilterSig[11], // mll cut && Z veto
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] &&                     passFilterSig[7] && passFilterSig[8] && passFilterSig[9] && passFilterSig[10] && passFilterSig[11], // Loose Z veto
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] &&                     passFilterSig[9] && passFilterSig[10] && passFilterSig[11], // met cut
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] &&                     passFilterSig[10] && passFilterSig[11], // mjj cut
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] && passFilterSig[9]                      && passFilterSig[11], // detajj cut
        passFilterSig[0] && passFilterSig[1] && passFilterSig[2] && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] && passFilterSig[9] && passFilterSig[10]                       // zeppenfeld cut
                               };

      bool totalSel = kTRUE;
      for(int isel=0; isel<numberCuts; isel++) {
        totalSel = totalSel && passFilterSig[isel];
	if(totalSel == kTRUE) sumEvol[typeSel]++;
      }

      bool passAllCuts[nSelTypes] = {passSignalRegion,passControlRegionTop,passControlRegionWZ,passControlRegionDi,passControlRegionSS2j,passControlRegionOS2j,passControlRegionZLL};                 

      TLorentzVector dijetjesUp,dijetjesDown;
      double deltaEtaJJjesUp = 0; double deltaEtaJJjesDown = 0;
      if(idJesUp.size() >= 2 ) {
        //dijetjesUp   = ( *(TLorentzVector*)(*eventJets.p4)[idJesUp[0]] )   + ( *(TLorentzVector*)(*eventJets.p4)[idJesUp[1]] );
	dijetjesUp.SetPx(((TLorentzVector*)(*eventJets.p4)[idJesUp[0]])->Px()*(1+(float)(*eventJets.unc)[idJesUp[0]])+((TLorentzVector*)(*eventJets.p4)[idJesUp[1]])->Px()*(1+(float)(*eventJets.unc)[idJesUp[1]]));
	dijetjesUp.SetPy(((TLorentzVector*)(*eventJets.p4)[idJesUp[0]])->Py()*(1+(float)(*eventJets.unc)[idJesUp[0]])+((TLorentzVector*)(*eventJets.p4)[idJesUp[1]])->Py()*(1+(float)(*eventJets.unc)[idJesUp[1]]));
	dijetjesUp.SetPz(((TLorentzVector*)(*eventJets.p4)[idJesUp[0]])->Pz()*(1+(float)(*eventJets.unc)[idJesUp[0]])+((TLorentzVector*)(*eventJets.p4)[idJesUp[1]])->Pz()*(1+(float)(*eventJets.unc)[idJesUp[1]]));
	dijetjesUp.SetE (((TLorentzVector*)(*eventJets.p4)[idJesUp[0]])-> E()*(1+(float)(*eventJets.unc)[idJesUp[0]])+((TLorentzVector*)(*eventJets.p4)[idJesUp[1]])-> E()*(1+(float)(*eventJets.unc)[idJesUp[1]]));
	deltaEtaJJjesUp = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJesUp[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJesUp[1]])->Eta());
      }
      if(idJesDown.size() >= 2) {
        //dijetjesDown = ( *(TLorentzVector*)(*eventJets.p4)[idJesDown[0]] ) + ( *(TLorentzVector*)(*eventJets.p4)[idJesDown[1]] );
	dijetjesDown.SetPx(((TLorentzVector*)(*eventJets.p4)[idJesDown[0]])->Px()*(1-(float)(*eventJets.unc)[idJesDown[0]])+((TLorentzVector*)(*eventJets.p4)[idJesDown[1]])->Px()*(1-(float)(*eventJets.unc)[idJesDown[1]]));
	dijetjesDown.SetPy(((TLorentzVector*)(*eventJets.p4)[idJesDown[0]])->Py()*(1-(float)(*eventJets.unc)[idJesDown[0]])+((TLorentzVector*)(*eventJets.p4)[idJesDown[1]])->Py()*(1-(float)(*eventJets.unc)[idJesDown[1]]));
	dijetjesDown.SetPz(((TLorentzVector*)(*eventJets.p4)[idJesDown[0]])->Pz()*(1-(float)(*eventJets.unc)[idJesDown[0]])+((TLorentzVector*)(*eventJets.p4)[idJesDown[1]])->Pz()*(1-(float)(*eventJets.unc)[idJesDown[1]]));
	dijetjesDown.SetE (((TLorentzVector*)(*eventJets.p4)[idJesDown[0]])-> E()*(1-(float)(*eventJets.unc)[idJesDown[0]])+((TLorentzVector*)(*eventJets.p4)[idJesDown[1]])-> E()*(1-(float)(*eventJets.unc)[idJesDown[1]]));
        deltaEtaJJjesDown = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJesDown[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJesDown[1]])->Eta());
      }

      TLorentzVector dijetjerUp,dijetjerDown;
      double deltaEtaJJjerUp = 0; double deltaEtaJJjerDown = 0;
      if(idJerUp.size() >= 2 ) {
        //dijetjerUp   = ( *(TLorentzVector*)(*eventJets.p4)[idJerUp[0]] )   + ( *(TLorentzVector*)(*eventJets.p4)[idJerUp[1]] );
        double sfjets_jer[2] = {1,1};
        if((float)(*eventJets.ptResUncCentral)[idJerUp[0]] > 0) sfjets_jer[0] = (float)(*eventJets.ptResUncUp)[idJerUp[0]] / (float)(*eventJets.ptResUncCentral)[idJerUp[0]];
        if((float)(*eventJets.ptResUncCentral)[idJerUp[1]] > 0) sfjets_jer[1] = (float)(*eventJets.ptResUncUp)[idJerUp[1]] / (float)(*eventJets.ptResUncCentral)[idJerUp[1]];
	dijetjerUp.SetPx(((TLorentzVector*)(*eventJets.p4)[idJerUp[0]])->Px()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJerUp[1]])->Px()*sfjets_jer[1]);
	dijetjerUp.SetPy(((TLorentzVector*)(*eventJets.p4)[idJerUp[0]])->Py()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJerUp[1]])->Py()*sfjets_jer[1]);
	dijetjerUp.SetPz(((TLorentzVector*)(*eventJets.p4)[idJerUp[0]])->Pz()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJerUp[1]])->Pz()*sfjets_jer[1]);
	dijetjerUp.SetE (((TLorentzVector*)(*eventJets.p4)[idJerUp[0]])-> E()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJerUp[1]])-> E()*sfjets_jer[1]);
	deltaEtaJJjerUp = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJerUp[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJerUp[1]])->Eta());
      }
      if(idJerDown.size() >= 2) {
        //dijetjerDown   = ( *(TLorentzVector*)(*eventJets.p4)[idJerDown[0]] )   + ( *(TLorentzVector*)(*eventJets.p4)[idJerDown[1]] );
        double sfjets_jer[2] = {1,1};
        if((float)(*eventJets.ptResUncCentral)[idJerDown[0]] > 0) sfjets_jer[0] = (float)(*eventJets.ptResUncDown)[idJerDown[0]] / (float)(*eventJets.ptResUncCentral)[idJerDown[0]];
        if((float)(*eventJets.ptResUncCentral)[idJerDown[1]] > 0) sfjets_jer[1] = (float)(*eventJets.ptResUncDown)[idJerDown[1]] / (float)(*eventJets.ptResUncCentral)[idJerDown[1]];
	dijetjerDown.SetPx(((TLorentzVector*)(*eventJets.p4)[idJerDown[0]])->Px()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJerDown[1]])->Px()*sfjets_jer[1]);
	dijetjerDown.SetPy(((TLorentzVector*)(*eventJets.p4)[idJerDown[0]])->Py()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJerDown[1]])->Py()*sfjets_jer[1]);
	dijetjerDown.SetPz(((TLorentzVector*)(*eventJets.p4)[idJerDown[0]])->Pz()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJerDown[1]])->Pz()*sfjets_jer[1]);
	dijetjerDown.SetE (((TLorentzVector*)(*eventJets.p4)[idJerDown[0]])-> E()*sfjets_jer[0]+((TLorentzVector*)(*eventJets.p4)[idJerDown[1]])-> E()*sfjets_jer[1]);
	deltaEtaJJjerDown = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJerDown[0]])->Eta()-((TLorentzVector*)(*eventJets.p4)[idJerDown[1]])->Eta());
      }

      bool passSystCuts[nSystTypes] = {
          passFilterSig[0] && passFilterSig[1] && idJesUp.size() >= 2	&& passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] 						   && dijetjesUp.M() > mjjCut   && deltaEtaJJjesUp > 2.5   && passFilterSig[11],
          passFilterSig[0] && passFilterSig[1] && idJesDown.size() >= 2 && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] 						   && dijetjesDown.M() > mjjCut && deltaEtaJJjesDown > 2.5 && passFilterSig[11],
          passFilterSig[0] && passFilterSig[1] && idJerUp.size() >= 2	&& passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] 						   && dijetjerUp.M() > mjjCut   && deltaEtaJJjerUp > 2.5   && passFilterSig[11],
          passFilterSig[0] && passFilterSig[1] && idJerDown.size() >= 2 && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && passFilterSig[8] 						   && dijetjerDown.M() > mjjCut && deltaEtaJJjerDown > 2.5 && passFilterSig[11],
          passFilterSig[0] && passFilterSig[1] && passFilterSig[2]      && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])->Pt()   > 40 && passFilterSig[9]          && passFilterSig[10]	   && passFilterSig[11],
          passFilterSig[0] && passFilterSig[1] && passFilterSig[2]      && passFilterSig[3] && passFilterSig[4] && passFilterSig[5] && passFilterSig[6] && passFilterSig[7] && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt() > 40 && passFilterSig[9]          && passFilterSig[10]	   && passFilterSig[11]
      };
 
      if     (theControlRegion == 1){
        passSystCuts[JESUP]   = passFilterCR1[0] && passFilterCR1[1] && idJesUp.size() >= 2   && passFilterCR1[3] && passFilterCR1[4] && passFilterCR1[5] && passFilterCR1[6] && passFilterCR1[7] && passFilterCR1[8]							 && dijetjesUp.M() > mjjCut   && deltaEtaJJjesUp > 2.5   && passFilterCR1[11];
        passSystCuts[JESDOWN] = passFilterCR1[0] && passFilterCR1[1] && idJesDown.size() >= 2 && passFilterCR1[3] && passFilterCR1[4] && passFilterCR1[5] && passFilterCR1[6] && passFilterCR1[7] && passFilterCR1[8]							 && dijetjesDown.M() > mjjCut && deltaEtaJJjesDown > 2.5 && passFilterCR1[11];
        passSystCuts[JERUP]   = passFilterCR1[0] && passFilterCR1[1] && idJerUp.size() >= 2   && passFilterCR1[3] && passFilterCR1[4] && passFilterCR1[5] && passFilterCR1[6] && passFilterCR1[7] && passFilterCR1[8]							 && dijetjerUp.M() > mjjCut   && deltaEtaJJjerUp > 2.5   && passFilterCR1[11];
        passSystCuts[JERDOWN] = passFilterCR1[0] && passFilterCR1[1] && idJerDown.size() >= 2 && passFilterCR1[3] && passFilterCR1[4] && passFilterCR1[5] && passFilterCR1[6] && passFilterCR1[7] && passFilterCR1[8]							 && dijetjerDown.M() > mjjCut && deltaEtaJJjerDown > 2.5 && passFilterCR1[11];
        passSystCuts[METUP]   = passFilterCR1[0] && passFilterCR1[1] && passFilterCR1[2]      && passFilterCR1[3] && passFilterCR1[4] && passFilterCR1[5] && passFilterCR1[6] && passFilterCR1[7] && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])->Pt()   > 40 && passFilterCR1[9]	      && passFilterCR1[10]	 && passFilterCR1[11];
        passSystCuts[METDOWN] = passFilterCR1[0] && passFilterCR1[1] && passFilterCR1[2]      && passFilterCR1[3] && passFilterCR1[4] && passFilterCR1[5] && passFilterCR1[6] && passFilterCR1[7] && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt() > 40 && passFilterCR1[9]	      && passFilterCR1[10]	 && passFilterCR1[11];
      }
      else if(theControlRegion == 2){
        passSystCuts[JESUP]   = passFilterCR2[0] && passFilterCR2[1] && idJesUp.size() >= 2   && passFilterCR2[3] && passFilterCR2[4] && passFilterCR2[5] && passFilterCR2[6] && passFilterCR2[7] && passFilterCR2[8]							 && dijetjesUp.M() > mjjCut   && deltaEtaJJjesUp > 2.5   && passFilterCR2[11];
        passSystCuts[JESDOWN] = passFilterCR2[0] && passFilterCR2[1] && idJesDown.size() >= 2 && passFilterCR2[3] && passFilterCR2[4] && passFilterCR2[5] && passFilterCR2[6] && passFilterCR2[7] && passFilterCR2[8]							 && dijetjesDown.M() > mjjCut && deltaEtaJJjesDown > 2.5 && passFilterCR2[11];
        passSystCuts[JERUP]   = passFilterCR2[0] && passFilterCR2[1] && idJerUp.size() >= 2   && passFilterCR2[3] && passFilterCR2[4] && passFilterCR2[5] && passFilterCR2[6] && passFilterCR2[7] && passFilterCR2[8]							 && dijetjerUp.M() > mjjCut   && deltaEtaJJjerUp > 2.5   && passFilterCR2[11];
        passSystCuts[JERDOWN] = passFilterCR2[0] && passFilterCR2[1] && idJerDown.size() >= 2 && passFilterCR2[3] && passFilterCR2[4] && passFilterCR2[5] && passFilterCR2[6] && passFilterCR2[7] && passFilterCR2[8]							 && dijetjerDown.M() > mjjCut && deltaEtaJJjerDown > 2.5 && passFilterCR2[11];
        passSystCuts[METUP]   = passFilterCR2[0] && passFilterCR2[1] && passFilterCR2[2]      && passFilterCR2[3] && passFilterCR2[4] && passFilterCR2[5] && passFilterCR2[6] && passFilterCR2[7] && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])->Pt()   > 40 && passFilterCR2[9]	      && passFilterCR2[10]	 && passFilterCR2[11];
        passSystCuts[METDOWN] = passFilterCR2[0] && passFilterCR2[1] && passFilterCR2[2]      && passFilterCR2[3] && passFilterCR2[4] && passFilterCR2[5] && passFilterCR2[6] && passFilterCR2[7] && ((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt() > 40 && passFilterCR2[9]	      && passFilterCR2[10]       && passFilterCR2[11];
      }

      // begin event weighting
      vector<bool> isGenDupl;
      vector<int>wBoson;
      vector<int>zBoson;
      int numberQuarks[2] = {0,0};
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 4 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[0]++;
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 5 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) numberQuarks[1]++;
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 24) {
	  wBoson.push_back(ngen0);
	}
        if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 23) {
	  zBoson.push_back(ngen0);
	}
        isGenDupl.push_back(0);
	bool isGoodFlags = ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::PromptFinalState) == BareMonteCarlo::PromptFinalState ||
            		   ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::DirectPromptTauDecayProductFinalState) == BareMonteCarlo::DirectPromptTauDecayProductFinalState;
        isGoodFlags = isGoodFlags && (TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 11 || TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 13);
        if(isGoodFlags == false) isGenDupl[ngen0] = 1;
      }

      int genLep = 0;
      int numberGoodGenLep[2] = {0,0};
      for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	if(isGenDupl[ngen] == 1) continue;
        genLep++;
	numberGoodGenLep[0]++;
	if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt() <= 10 ||
	   TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()) >= 2.5) continue;
	numberGoodGenLep[1]++;
      }

      double total_WS_SF[2] = {1.0, 1.0};
      vector<int> isGenLep; unsigned int goodIsGenRSLep = 0; unsigned int goodIsGenWSLep = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        bool isGenRSLepton = false;
        bool isGenWSLepton = false;
        for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	  if(isGenDupl[ngen] == 1) continue;
          if((int)(*eventLeptons.pdgId)[idLep[nl]] == (int)(*eventMonteCarlo.pdgId)[ngen] &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.3) {
	    isGenRSLepton = true;
	  }
          else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.3) {
	    isGenWSLepton = true;
	  }
	}
	if     (isGenRSLepton == true) {isGenLep.push_back(1); goodIsGenRSLep++;}
	else if(isGenWSLepton == true) {isGenLep.push_back(2); goodIsGenWSLep++;func_ws_sf(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),total_WS_SF);}
	else                           {isGenLep.push_back(0);}
      }

      // trigger efficiency
      double trigEff = 1.0;
      //if(infilecatv[ifile] != 0) {
      //  trigEff = trigLookup.GetExpectedTriggerEfficiency(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),
      //  						  ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),
      //  						 TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]));
      //}
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
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
	  	typeLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,true,fhDMuIsoSF,fhDVeryTightSF,true);
        }
      }

      int theCategory = infilecatv[ifile];
      int nFakeCount = 0;

      // wrong sign
      if((theCategory == 3 || theCategory == 4 || theCategory == 5) && goodIsGenWSLep > 0) theCategory  = 6;

      // fake rate
      unsigned int typeFakeLepton[2] = {0,0};
      double fakeSF = 1.0;
      if(usePureMC == false && sigModel == -1){
	if((infilecatv[ifile] == 0 || infilecatv[ifile] == 7 || (goodIsGenRSLep+goodIsGenWSLep) == isGenLep.size()) && goodIsTight != idTight.size()){
            for(unsigned int nl=0; nl<idLep.size(); nl++){
              if(idTight[nl] == 1) continue;
              fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
              theCategory = 9;
	      if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) {typeFakeLepton[0]++; nFakeCount = nFakeCount + 1;}
	      else                                                        {typeFakeLepton[1]++; nFakeCount = nFakeCount + 2;}
            }
            if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) fakeSF =  1.0 * fakeSF; // double fake, MC
            else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) fakeSF = -1.0 * fakeSF; // single fake, MC
            else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) fakeSF = -1.0 * fakeSF; // double fake, data
            else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) fakeSF =  1.0 * fakeSF; // single fake, data
            //if(typeFakeLepton[0] < typeFakeLepton[1] && idLep.size() == 2) {theCategory = 10; typeFakeLepton[1] = 0}
	}
	else if(infilecatv[ifile] != 0 && infilecatv[ifile] != 7 && (goodIsGenRSLep+goodIsGenWSLep) != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
	}
	else if(infilecatv[ifile] != 0 && (goodIsGenRSLep+goodIsGenWSLep) == isGenLep.size()){ // MC with all good leptons
          fakeSF = 1.0;
	}
	else if(infilecatv[ifile] == 0 || infilecatv[ifile] == 7){ // data or W+gamma with all good leptons
          fakeSF = 1.0;
	}
	else {
          printf("PROBLEMFAKES: %d %d %d %d %d %d\n",infilecatv[ifile],goodIsGenRSLep,goodIsGenWSLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
          assert(0);
	}
      }

      if(isBlinded == true && passFilterSig[0] == kTRUE && passFilterSig[9] == kTRUE) continue;

      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale*trigEff;

      // Wrong sign scale factor
      if(theCategory == 6 && useWSFromData  && usePureMC == false) {
        totalWeight = totalWeight * total_WS_SF[0];
      }

      // Wrong sign efficiency in OS events
      //if(passControlRegionDi || passControlRegionOS2j) {
      //  totalWeight = totalWeight * func_ws_eff(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta(),
      //                                          ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta(),
      //                                          fhDveryTightWrongSignEff);
      //}

      // Z->ll scale factor
      if     (infilenamev[ifile].Contains("JetsToLL") && goodIsTight == idTight.size()) {
        totalWeight = totalWeight * the_sf_ZLL;
      }
      else if(infilenamev[ifile].Contains("JetsToLL") && goodIsTight != idTight.size() && (passAllCuts[SIGSEL] || passAllCuts[TOPSEL] || passAllCuts[WZSEL])) {
        totalWeight = 0.0;
      }
      if(infilenamev[ifile].Contains("JetsToLL") && typePair == 0 && signQ != 0) {
        totalWeight = 0.0;
        if(passSignalRegion && goodIsGenWSLep == 0) printf("mm event with no WS candidates\n");
      }

      // Btag scale factor
      if     (idLep.size() == 2) totalWeight = totalWeight * total_bjet_probLOOSE[1]/total_bjet_probLOOSE[0];
      else if(idLep.size() == 3) totalWeight = totalWeight * total_bjet_probTIGHT[1]/total_bjet_probTIGHT[0];
      else                       printf("THIS IS NOT POSSIBLE!: %zu\n",idLep.size());

      double btagCorr[2] = {(total_bjet_probLOOSEUP[1]  /total_bjet_probLOOSEUP[0]  )/(total_bjet_probLOOSE[1]/total_bjet_probLOOSE[0]),
                            (total_bjet_probLOOSEDOWN[1]/total_bjet_probLOOSEDOWN[0])/(total_bjet_probLOOSE[1]/total_bjet_probLOOSE[0])};
      if(idLep.size() == 3){
        btagCorr[0] = (total_bjet_probTIGHTUP[1]  /total_bjet_probTIGHTUP[0]  )/(total_bjet_probTIGHT[1]/total_bjet_probTIGHT[0]);
        btagCorr[1] = (total_bjet_probTIGHTDOWN[1]/total_bjet_probTIGHTDOWN[0])/(total_bjet_probTIGHT[1]/total_bjet_probTIGHT[0]);
      }

      if(theCategory == -1) {theCategory = 1; totalWeight = -1.0 * totalWeight;}

      if(totalWeight == 0) continue;
      // end event weighting

      if(passSignalRegion && infilecatv[ifile] == 0) totalFakeDataCount[typeSel][nFakeCount] = totalFakeDataCount[typeSel][nFakeCount] + 1;

      for(int nl=0; nl <=sumEvol[typeSel]; nl++) if(fakeSF == 1) {hDWWLL[typeSel]->Fill((double)nl,totalWeight);hDWWLL[6]->Fill((double)nl,totalWeight);}

      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i] && sigModel == -1) {
          bgdDecay[typeSel][i][theCategory] += totalWeight;
          weiDecay[typeSel][i][theCategory] += totalWeight*totalWeight;
        }
      }

      // Making histograms for datacards
      double MVAVar = TMath::Min(dijet.M(),1999.999)+2000.*typeSel;
      double MVAVarJESSyst[2] = {TMath::Min(dijetjesUp.M(),1999.999)+2000.*typeSel,TMath::Min(dijetjesDown.M(),1999.999)+2000.*typeSel};
      double MVAVarJERSyst[2] = {TMath::Min(dijetjerUp.M(),1999.999)+2000.*typeSel,TMath::Min(dijetjerDown.M(),1999.999)+2000.*typeSel};

      if     (theControlRegion == 1){
      }
      else if(theControlRegion == 2){
        MVAVar = TMath::Min(dijet.M(),1999.999); 
	MVAVarJESSyst[0] = TMath::Min(dijetjesUp.M(),1999.999); MVAVarJESSyst[1] = TMath::Min(dijetjesDown.M(),1999.999);
	MVAVarJERSyst[0] = TMath::Min(dijetjerUp.M(),1999.999); MVAVarJERSyst[1] = TMath::Min(dijetjerDown.M(),1999.999);
      }

      if     (finalVar == 1 || finalVar == 2 || finalVar == 3){
        MVAVar = TMath::Min(dilep.M(),599.999)+1000.*typeSel;
        if(theControlRegion == 2){
          MVAVar = TMath::Min(dilep.M(),599.999);
        }
        MVAVarJESSyst[0] = MVAVar;
        MVAVarJESSyst[1] = MVAVar;
        MVAVarJERSyst[0] = MVAVar;
        MVAVarJERSyst[1] = MVAVar;
      }
      else if(finalVar == 4){
        int typeSelAux = 0;
        if     (dilep.M() < 100) typeSelAux = 0;
        else if(dilep.M() < 180) typeSelAux = 1;
        else if(dilep.M() < 300) typeSelAux = 2;
        else                     typeSelAux = 3;

        MVAVar = TMath::Min(dijet.M(),1999.999)+2000.*typeSelAux;
        MVAVarJESSyst[0] = TMath::Min(dijetjesUp.M(),1999.999)+2000.*typeSelAux;
        MVAVarJESSyst[1] = TMath::Min(dijetjesDown.M(),1999.999)+2000.*typeSelAux;
        MVAVarJERSyst[0] = TMath::Min(dijetjerUp.M(),1999.999)+2000.*typeSelAux;
        MVAVarJERSyst[1] = TMath::Min(dijetjerDown.M(),1999.999)+2000.*typeSelAux;

        if(theControlRegion == 1 || theControlRegion == 2){
          MVAVar = TMath::Min(dijet.M(),1999.999); 
          MVAVarJESSyst[0] = TMath::Min(dijetjesUp.M(),1999.999);
          MVAVarJESSyst[1] = TMath::Min(dijetjesDown.M(),1999.999);
          MVAVarJERSyst[0] = TMath::Min(dijetjerUp.M(),1999.999);
          MVAVarJERSyst[1] = TMath::Min(dijetjerDown.M(),1999.999);
        }
      }

      for(int thePlot=0; thePlot<allPlots; thePlot++){
        double theVar = 0.0;
        bool makePlot = false;
        if     (thePlot ==  0 && passNMinusOne[5])          {makePlot = true;theVar = TMath::Min(dijet.M(),1999.999);}
        else if(thePlot ==  1 && passNMinusOne[6])          {makePlot = true;theVar = TMath::Min(deltaEtaJJ,7.999);}
        else if(thePlot ==  2 && passSignalRegion)          {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
        else if(thePlot ==  3 && passNMinusOne[4])          {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),399.999);}
        else if(thePlot ==  4 && passNMinusOne[2])          {makePlot = true;theVar = TMath::Min(dilep.M(),399.999);}
        else if(thePlot ==  5 && passSignalRegion)          {makePlot = true;theVar = TMath::Min((double)eventVertex.npv,39.499);}
	else if(thePlot ==  6 && passNMinusOne[1])          {makePlot = true;theVar = TMath::Min((double)numberGoodTaus,3.499);}
	else if(thePlot ==  7 && passNMinusOne[3])          {makePlot = true;theVar = TMath::Min(TMath::Abs(minMassLooseZ-91.1876),39.999);}
	else if(thePlot ==  8 && passSignalRegion)          {makePlot = true;theVar = (double)(numberGoodGenLep[1]+10*numberGoodGenLep[0]);}
	else if(thePlot ==  9 && passNMinusOne[0])          {makePlot = true;theVar = TMath::Min(TMath::Max(bDiscrMax,0.001),0.999);}
	else if(thePlot == 10 && passNMinusOne[0])          {makePlot = true;theVar = TMath::Min((double)idSoft.size(),3.499);}

        else if(thePlot == 11 && passControlRegionDi)       {makePlot = true;theVar = TMath::Min(dijet.M(),1999.999);}
        else if(thePlot == 12 && passControlRegionDi)       {makePlot = true;theVar = TMath::Min(dilep.M(),599.999);}

        else if(thePlot == 13 && passSignalRegion)          {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Pt(),199.999);}
        else if(thePlot == 14 && passSignalRegion)          {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Pt(),199.999);}
        else if(thePlot == 15 && passControlRegionSS2j)     {makePlot = true;theVar = TMath::Min(deltaEtaJJ,7.999);}
        else if(thePlot == 16 && passControlRegionSS2j)     {makePlot = true;theVar = ((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta();}
        else if(thePlot == 17 && passControlRegionSS2j)     {makePlot = true;theVar = ((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta();}
        else if(thePlot == 18 && passControlRegionSS2j)     {makePlot = true;theVar = TMath::Min(dilep.M(),399.999);}
        else if(thePlot == 19 && passControlRegionSS2j)     {makePlot = true;theVar = TMath::Max(theLeptonZ[0],theLeptonZ[1]);}
	else if(thePlot == 20 && passControlRegionWZ)       {makePlot = true;theVar = TMath::Min(TMath::Abs(minMassZ-91.1876),19.999);}
        else if(thePlot == 21 && passNMinusOne[7])          {makePlot = true;theVar = TMath::Min(theLeptonZ[0],theLeptonZ[1]);}
        else if(thePlot == 22 && passNMinusOne[7])          {makePlot = true;theVar = TMath::Max(theLeptonZ[0],theLeptonZ[1]);}
        else if(thePlot == 23 && passSignalRegion)          {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),199.999);}
        else if(thePlot == 24 && passSignalRegion)          {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),199.999);}
        else if(thePlot == 25 && passControlRegionZLL)      {makePlot = true;theVar = TMath::Min((double)idJet.size(),6.499);}
        else if(thePlot == 26 && passSignalRegion)          {makePlot = true;theVar = TMath::Min(dilep.Pt(),499.999);}
        else if(thePlot == 27 && passLooseControlRegionWZ)  {makePlot = true;theVar = TMath::Min(dijet.M(),1999.999);}
        else if(thePlot == 28 && passLooseControlRegionTop) {makePlot = true;theVar = TMath::Min(deltaEtaJJ,7.999);}
        else if(thePlot == 29 && passLooseControlRegionWZ)  {makePlot = true;theVar = TMath::Min(deltaEtaJJ,7.999);}
        else if(thePlot == 30 && passNMinusOne[3])          {makePlot = true;theVar = TMath::Min((double)idLepLoose.size(),3.499);}

        else if(thePlot == 31 && passSignalRegion)          {makePlot = true;theVar = TMath::Min(dijet.M(),1999.999);}
        else if(thePlot == 32 && passControlRegionTop)      {makePlot = true;theVar = TMath::Min(dijet.M(),1999.999);}
        else if(thePlot == 33 && passControlRegionWZ)       {makePlot = true;theVar = TMath::Min(dijet.M(),1999.999);}

        else if(thePlot == 34 && passSignalRegion)          {makePlot = true;theVar = TMath::Min(dilep.M(),599.999);}
        else if(thePlot == 35 && passControlRegionTop)      {makePlot = true;theVar = TMath::Min(dilep.M(),599.999);}
        else if(thePlot == 36 && passControlRegionWZ)       {makePlot = true;theVar = TMath::Min(dilep.M(),599.999);}

        else if(thePlot == 37 && passSignalRegion)          {makePlot = true;theVar = MVAVar;}
        else if(thePlot == 38 && passControlRegionTop)      {makePlot = true;theVar = MVAVar;}
        else if(thePlot == 39 && passControlRegionWZ)       {makePlot = true;theVar = MVAVar;}

        if(makePlot && sigModel <= 0) histo[typeSel][thePlot][theCategory]  ->Fill(theVar,totalWeight);
        if(makePlot && sigModel <= 0) histo[6]      [thePlot][theCategory]  ->Fill(theVar,totalWeight);
        if(makePlot && sigModel == 1) histo[typeSel][thePlot][theCategory+1]->Fill(theVar,totalWeight);
        if(makePlot && sigModel == 1) histo[6]      [thePlot][theCategory+1]->Fill(theVar,totalWeight);
      }

      // Avoid QCD scale and PDF weights that are anomalous high
      double maxQCDscale = (TMath::Abs((double)eventMonteCarlo.r1f2)+TMath::Abs((double)eventMonteCarlo.r1f5)+TMath::Abs((double)eventMonteCarlo.r2f1)+
        		    TMath::Abs((double)eventMonteCarlo.r2f2)+TMath::Abs((double)eventMonteCarlo.r5f1)+TMath::Abs((double)eventMonteCarlo.r5f5))/6.0;
      double PDFAvg = 0.0;
      if(infilecatv[ifile] != 0 && ((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2))){
        if(initPDFTag != -1)
        for(int npdf=0; npdf<100; npdf++) PDFAvg = PDFAvg + TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]);
        PDFAvg = PDFAvg/100.0;
      }

      if     (theCategory == 0){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_Data->Fill(MVAVar,totalWeight);
        }
      }
      else if(theCategory == 1){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_EWK->Fill(MVAVar,totalWeight);
           histo_EWK_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
           histo_EWK_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
           histo_EWK_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
           histo_EWK_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
           histo_EWK_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
           histo_EWK_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
           if(initPDFTag != -1)
           for(int npdf=0; npdf<100; npdf++) histo_EWK_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
           else
           for(int npdf=0; npdf<100; npdf++) histo_EWK_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
           if	  (typePair == 0) histo_EWK_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else if(typePair == 1) histo_EWK_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else                  {histo_EWK_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_EWK_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}
           if	  (typePair == 0) histo_EWK_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else if(typePair == 1) histo_EWK_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else                  {histo_EWK_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);
	                          histo_EWK_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);}
           if	  (typePair == 0) histo_EWK_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else if(typePair == 1) histo_EWK_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else                  {histo_EWK_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);
	                          histo_EWK_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);}
           histo_EWK_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
           histo_EWK_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
           histo_EWK_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
           histo_EWK_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);

          bool passFiducial[2] = {false, false};
	  // Begin gen fiducial selection
	  Int_t countSelectedGenLeptons = 0;
	  vector<int> idGenLep;
	  for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
	    if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 && TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) continue;
	    bool isGoodFlags = ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::PromptFinalState) == BareMonteCarlo::PromptFinalState ||
          		       ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::DirectPromptTauDecayProductFinalState) == BareMonteCarlo::DirectPromptTauDecayProductFinalState;
	    if(!isGoodFlags) continue;
	    idGenLep.push_back(ngen0);
	    bool passSelLepton = 
               //(((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::PromptFinalState)                      == BareMonteCarlo::PromptFinalState ||
               // ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::DirectPromptTauDecayProductFinalState) == BareMonteCarlo::DirectPromptTauDecayProductFinalState) &&
               ((*eventMonteCarlo.flags)[ngen0] & BareMonteCarlo::PromptFinalState)                      == BareMonteCarlo::PromptFinalState &&
               TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Eta()) < 2.5 &&
               ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 20;
	    if(passSelLepton) countSelectedGenLeptons ++;
	  }

	  if(countSelectedGenLeptons >= 2){
	    vector<int> idGenJet30;
	    for(int njetgen=0; njetgen<eventMonteCarlo.jetP4->GetEntriesFast(); njetgen++) {
              if(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[njetgen])->Eta()) >= 5.0) continue;
              bool isGenLepton = false;
              for(unsigned int nglep = 0; nglep<idGenLep.size(); nglep++) {
             	if(((TLorentzVector*)(*eventMonteCarlo.p4)[idGenLep[nglep]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.jetP4)[njetgen])) < 0.3) {
	     	  isGenLepton = true;
             	  break;
	     	}
              }
              if(isGenLepton) continue;
              if(((TLorentzVector*)(*eventMonteCarlo.jetP4)[njetgen])->Pt() > 30) idGenJet30.push_back(njetgen);
	    }

	    if(idGenJet30.size() >= 2){
              TLorentzVector dijet =	   ( *(TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet30[0]] ) +   ( *(TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet30[1]] );
	      double deltaEtaJJ = TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet30[0]])->Eta()-((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet30[1]])->Eta());
              if(deltaEtaJJ > 2.5 && dijet.M() > 300) passFiducial[0] = true;
              if(deltaEtaJJ > 2.5 && dijet.M() > 500) passFiducial[1] = true;
	    }
	  }
	  // End gen fiducial selection
          if(passFiducial[0]) selectedFiducial[0]    = selectedFiducial[0]    + totalWeight;
          else                selectedNonFiducial[0] = selectedNonFiducial[0] + totalWeight;
          if(passFiducial[1]) selectedFiducial[1]    = selectedFiducial[1]    + totalWeight;
          else                selectedNonFiducial[1] = selectedNonFiducial[1] + totalWeight;
        }
        if(passSystCuts[JESUP])  histo_EWK_CMS_MVAJESBoundingUp  ->Fill(MVAVarJESSyst[0],totalWeight);
        if(passSystCuts[JESDOWN])histo_EWK_CMS_MVAJESBoundingDown->Fill(MVAVarJESSyst[1],totalWeight);
        if(passSystCuts[JERUP])  histo_EWK_CMS_MVAJERBoundingUp  ->Fill(MVAVarJERSyst[0],totalWeight);
        if(passSystCuts[JERDOWN])histo_EWK_CMS_MVAJERBoundingDown->Fill(MVAVarJERSyst[1],totalWeight);
        if(passSystCuts[METUP])  histo_EWK_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
        if(passSystCuts[METDOWN])histo_EWK_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
      }
      else if(theCategory == 2){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_QCD->Fill(MVAVar,totalWeight);
           histo_QCD_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
           histo_QCD_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
           histo_QCD_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
           histo_QCD_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
           histo_QCD_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
           histo_QCD_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
           if(initPDFTag != -1)
           for(int npdf=0; npdf<100; npdf++) histo_QCD_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
           else
           for(int npdf=0; npdf<100; npdf++) histo_QCD_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
           if	  (typePair == 0) histo_QCD_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else if(typePair == 1) histo_QCD_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else                  {histo_QCD_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_QCD_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}
           if	  (typePair == 0) histo_QCD_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else if(typePair == 1) histo_QCD_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else                  {histo_QCD_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);
	                          histo_QCD_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);}
           if	  (typePair == 0) histo_QCD_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else if(typePair == 1) histo_QCD_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else                  {histo_QCD_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);
	                          histo_QCD_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);}
           histo_QCD_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
           histo_QCD_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
           histo_QCD_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
           histo_QCD_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
        }
        if(passSystCuts[JESUP])  histo_QCD_CMS_MVAJESBoundingUp  ->Fill(MVAVarJESSyst[0],totalWeight);
        if(passSystCuts[JESDOWN])histo_QCD_CMS_MVAJESBoundingDown->Fill(MVAVarJESSyst[1],totalWeight);
        if(passSystCuts[JERUP])  histo_QCD_CMS_MVAJERBoundingUp  ->Fill(MVAVarJERSyst[0],totalWeight);
        if(passSystCuts[JERDOWN])histo_QCD_CMS_MVAJERBoundingDown->Fill(MVAVarJERSyst[1],totalWeight);
        if(passSystCuts[METUP])  histo_QCD_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
        if(passSystCuts[METDOWN])histo_QCD_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
      }
      else if(theCategory == 3){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_WZ->Fill(MVAVar,totalWeight);
           histo_WZ_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
           histo_WZ_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
           histo_WZ_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
           histo_WZ_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
           histo_WZ_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
           histo_WZ_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
           if(initPDFTag != -1)
           for(int npdf=0; npdf<100; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
           else
           for(int npdf=0; npdf<100; npdf++) histo_WZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
           if	  (typePair == 0) histo_WZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else if(typePair == 1) histo_WZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else                  {histo_WZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_WZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}
           if	  (typePair == 0) histo_WZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else if(typePair == 1) histo_WZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else                  {histo_WZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);
	                          histo_WZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);}
           if	  (typePair == 0) histo_WZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else if(typePair == 1) histo_WZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else                  {histo_WZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);
	                          histo_WZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);}
           histo_WZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
           histo_WZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
           histo_WZ_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
           histo_WZ_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
        }
        if(passSystCuts[JESUP])  histo_WZ_CMS_MVAJESBoundingUp  ->Fill(MVAVarJESSyst[0],totalWeight);
        if(passSystCuts[JESDOWN])histo_WZ_CMS_MVAJESBoundingDown->Fill(MVAVarJESSyst[1],totalWeight);
        if(passSystCuts[JERUP])  histo_WZ_CMS_MVAJERBoundingUp  ->Fill(MVAVarJERSyst[0],totalWeight);
        if(passSystCuts[JERDOWN])histo_WZ_CMS_MVAJERBoundingDown->Fill(MVAVarJERSyst[1],totalWeight);
        if(passSystCuts[METUP])  histo_WZ_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
        if(passSystCuts[METDOWN])histo_WZ_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
      }
      else if(theCategory == 4){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_ZZ->Fill(MVAVar,totalWeight);
           histo_ZZ_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
           histo_ZZ_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
           histo_ZZ_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
           histo_ZZ_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
           histo_ZZ_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
           histo_ZZ_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
           if(initPDFTag != -1)
           for(int npdf=0; npdf<100; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
           else
           for(int npdf=0; npdf<100; npdf++) histo_ZZ_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
           if	  (typePair == 0) histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else if(typePair == 1) histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else                  {histo_ZZ_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_ZZ_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}
           if	  (typePair == 0) histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else if(typePair == 1) histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else                  {histo_ZZ_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);
	                          histo_ZZ_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);}
           if	  (typePair == 0) histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else if(typePair == 1) histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else                  {histo_ZZ_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);
	                          histo_ZZ_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);}
           histo_ZZ_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
           histo_ZZ_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
           histo_ZZ_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
           histo_ZZ_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
        }
        if(passSystCuts[JESUP])  histo_ZZ_CMS_MVAJESBoundingUp  ->Fill(MVAVarJESSyst[0],totalWeight);
        if(passSystCuts[JESDOWN])histo_ZZ_CMS_MVAJESBoundingDown->Fill(MVAVarJESSyst[1],totalWeight);
        if(passSystCuts[JERUP])  histo_ZZ_CMS_MVAJERBoundingUp  ->Fill(MVAVarJERSyst[0],totalWeight);
        if(passSystCuts[JERDOWN])histo_ZZ_CMS_MVAJERBoundingDown->Fill(MVAVarJERSyst[1],totalWeight);
        if(passSystCuts[METUP])  histo_ZZ_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
        if(passSystCuts[METDOWN])histo_ZZ_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
      }
      else if(theCategory == 5){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_VVV->Fill(MVAVar,totalWeight);
           histo_VVV_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
           histo_VVV_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
           histo_VVV_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
           histo_VVV_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
           histo_VVV_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
           histo_VVV_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
           if(initPDFTag != -1)
           for(int npdf=0; npdf<100; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
           else
           for(int npdf=0; npdf<100; npdf++) histo_VVV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
           if	  (typePair == 0) histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else if(typePair == 1) histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else                  {histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}
           if	  (typePair == 0) histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else if(typePair == 1) histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else                  {histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);
	                          histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);}
           if	  (typePair == 0) histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else if(typePair == 1) histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else                  {histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);
	                          histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);}
           histo_VVV_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
           histo_VVV_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
           histo_VVV_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
           histo_VVV_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
        }
        if(passSystCuts[JESUP])  histo_VVV_CMS_MVAJESBoundingUp  ->Fill(MVAVarJESSyst[0],totalWeight);
        if(passSystCuts[JESDOWN])histo_VVV_CMS_MVAJESBoundingDown->Fill(MVAVarJESSyst[1],totalWeight);
        if(passSystCuts[JERUP])  histo_VVV_CMS_MVAJERBoundingUp  ->Fill(MVAVarJERSyst[0],totalWeight);
        if(passSystCuts[JERDOWN])histo_VVV_CMS_MVAJERBoundingDown->Fill(MVAVarJERSyst[1],totalWeight);
        if(passSystCuts[METUP])  histo_VVV_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
        if(passSystCuts[METDOWN])histo_VVV_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
      }
      else if(theCategory == 6){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_WS             ->Fill(MVAVar,totalWeight);
           histo_WS_CMS_WSSFUp  ->Fill(MVAVar,totalWeight*total_WS_SF[1]);
           histo_WS_CMS_WSSFDown->Fill(MVAVar,totalWeight/total_WS_SF[1]);
           histo_WS_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
           histo_WS_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
           histo_WS_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
           histo_WS_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
           histo_WS_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
           histo_WS_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
           if(initPDFTag != -1)
           for(int npdf=0; npdf<100; npdf++) histo_WS_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
           else
           for(int npdf=0; npdf<100; npdf++) histo_WS_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
           if	  (typePair == 0) histo_WS_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else if(typePair == 1) histo_WS_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else                  {histo_WS_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_WS_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}
           if	  (typePair == 0) histo_WS_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else if(typePair == 1) histo_WS_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else                  {histo_WS_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);
	                          histo_WS_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);}
           if	  (typePair == 0) histo_WS_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else if(typePair == 1) histo_WS_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else                  {histo_WS_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);
	                          histo_WS_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);}
           histo_WS_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
           histo_WS_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
           histo_WS_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
           histo_WS_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
        }
        if(passSystCuts[JESUP])  histo_WS_CMS_MVAJESBoundingUp  ->Fill(MVAVarJESSyst[0],totalWeight);
        if(passSystCuts[JESDOWN])histo_WS_CMS_MVAJESBoundingDown->Fill(MVAVarJESSyst[1],totalWeight);
        if(passSystCuts[JERUP])  histo_WS_CMS_MVAJERBoundingUp  ->Fill(MVAVarJERSyst[0],totalWeight);
        if(passSystCuts[JERDOWN])histo_WS_CMS_MVAJERBoundingDown->Fill(MVAVarJERSyst[1],totalWeight);
        if(passSystCuts[METUP])  histo_WS_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
        if(passSystCuts[METDOWN])histo_WS_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
      }
      else if(theCategory == 7){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_WG->Fill(MVAVar,totalWeight);
           histo_WG_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
           histo_WG_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
           histo_WG_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
           histo_WG_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
           histo_WG_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
           histo_WG_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
           if(initPDFTag != -1)
           for(int npdf=0; npdf<100; npdf++) histo_WG_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
           else
           for(int npdf=0; npdf<100; npdf++) histo_WG_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
            if	  (typePair == 0) histo_WG_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else if(typePair == 1) histo_WG_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else                  {histo_WG_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_WG_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}
           if	  (typePair == 0) histo_WG_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else if(typePair == 1) histo_WG_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else                  {histo_WG_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);
	                          histo_WG_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);}
           if	  (typePair == 0) histo_WG_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else if(typePair == 1) histo_WG_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else                  {histo_WG_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);
	                          histo_WG_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);}
           histo_WG_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
           histo_WG_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
           histo_WG_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
           histo_WG_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
        }
        if(passSystCuts[JESUP])  histo_WG_CMS_MVAJESBoundingUp  ->Fill(MVAVarJESSyst[0],totalWeight);
        if(passSystCuts[JESDOWN])histo_WG_CMS_MVAJESBoundingDown->Fill(MVAVarJESSyst[1],totalWeight);
        if(passSystCuts[JERUP])  histo_WG_CMS_MVAJERBoundingUp  ->Fill(MVAVarJERSyst[0],totalWeight);
        if(passSystCuts[JERDOWN])histo_WG_CMS_MVAJERBoundingDown->Fill(MVAVarJERSyst[1],totalWeight);
        if(passSystCuts[METUP])  histo_WG_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
        if(passSystCuts[METDOWN])histo_WG_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
      }
      else if(theCategory == 8){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_DPS->Fill(MVAVar,totalWeight);
           histo_DPS_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
           histo_DPS_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
           histo_DPS_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
           histo_DPS_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
           histo_DPS_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
           histo_DPS_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
           for(int npdf=0; npdf<100; npdf++) histo_DPS_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
           if	  (typePair == 0) histo_DPS_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else if(typePair == 1) histo_DPS_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
           else                  {histo_DPS_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*0.50);
	                          histo_DPS_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*0.50);}
           if	  (typePair == 0) histo_DPS_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else if(typePair == 1) histo_DPS_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*1.02);
           else                  {histo_DPS_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);
	                          histo_DPS_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*0.50*1.02);}
           if	  (typePair == 0) histo_DPS_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else if(typePair == 1) histo_DPS_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.98);
           else                  {histo_DPS_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);
	                          histo_DPS_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight*0.50*0.98);}
           histo_DPS_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
           histo_DPS_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
           histo_DPS_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
           histo_DPS_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
        }
        if(passSystCuts[JESUP])  histo_DPS_CMS_MVAJESBoundingUp  ->Fill(MVAVarJESSyst[0],totalWeight);
        if(passSystCuts[JESDOWN])histo_DPS_CMS_MVAJESBoundingDown->Fill(MVAVarJESSyst[1],totalWeight);
        if(passSystCuts[JERUP])  histo_DPS_CMS_MVAJERBoundingUp  ->Fill(MVAVarJERSyst[0],totalWeight);
        if(passSystCuts[JERDOWN])histo_DPS_CMS_MVAJERBoundingDown->Fill(MVAVarJERSyst[1],totalWeight);
        if(passSystCuts[METUP])  histo_DPS_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
        if(passSystCuts[METDOWN])histo_DPS_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
      }
      else if(theCategory == 9){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
          if(typeFakeLepton[0] < typeFakeLepton[1] && idLep.size() == 2) histo_Fake_CMS_SystE->Fill(MVAVar,totalWeight);
	  else  							 histo_Fake_CMS_SystM->Fill(MVAVar,totalWeight);
          histo_FakeM->Fill(MVAVar,totalWeight);
          //if(dijet.M() > 1500 && dilep.M() > 180) {
          //  printf("HIGHMJJMLL: %f %f %d %f / %f %f %f %f %f %f %f\n",dijet.M(),dilep.M(),typeSel,totalWeight,mcWeight,theLumi,puWeight,effSF,fakeSF,theMCPrescale,trigEff);
          //}
        }
      }
      else if(theCategory == 10){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_FakeE->Fill(MVAVar,totalWeight);
        }
      }
      else if(theCategory == 11){
        if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1) || (passAllCuts[WZSEL] && theControlRegion == 2)) {
           histo_Higgs[sigModel]->Fill(MVAVar,totalWeight);
           histo_Higgs_CMS_QCDScaleBounding[sigModel][0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
           histo_Higgs_CMS_QCDScaleBounding[sigModel][1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
           histo_Higgs_CMS_QCDScaleBounding[sigModel][2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
           histo_Higgs_CMS_QCDScaleBounding[sigModel][3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
           histo_Higgs_CMS_QCDScaleBounding[sigModel][4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
           histo_Higgs_CMS_QCDScaleBounding[sigModel][5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
           if(initPDFTag != -1)
           for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBounding[sigModel][npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
           else
           for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBounding[sigModel][npdf]->Fill(MVAVar,totalWeight);
           if	  (typePair == 0) histo_Higgs_CMS_MVALepEffMBoundingAvg [sigModel]->Fill(MVAVar,totalWeight*1.00);
           else if(typePair == 1) histo_Higgs_CMS_MVALepEffEBoundingAvg [sigModel]->Fill(MVAVar,totalWeight*1.00);
           else                  {histo_Higgs_CMS_MVALepEffMBoundingAvg [sigModel]->Fill(MVAVar,totalWeight*0.50);
	                          histo_Higgs_CMS_MVALepEffEBoundingAvg [sigModel]->Fill(MVAVar,totalWeight*0.50);}
           if	  (typePair == 0) histo_Higgs_CMS_MVALepEffMBoundingUp  [sigModel]->Fill(MVAVar,totalWeight*1.02);
           else if(typePair == 1) histo_Higgs_CMS_MVALepEffEBoundingUp  [sigModel]->Fill(MVAVar,totalWeight*1.02);
           else                  {histo_Higgs_CMS_MVALepEffMBoundingUp  [sigModel]->Fill(MVAVar,totalWeight*0.50*1.02);
	                          histo_Higgs_CMS_MVALepEffEBoundingUp  [sigModel]->Fill(MVAVar,totalWeight*0.50*1.02);}
           if	  (typePair == 0) histo_Higgs_CMS_MVALepEffMBoundingDown[sigModel]->Fill(MVAVar,totalWeight*0.98);
           else if(typePair == 1) histo_Higgs_CMS_MVALepEffEBoundingDown[sigModel]->Fill(MVAVar,totalWeight*0.98);
           else                  {histo_Higgs_CMS_MVALepEffMBoundingDown[sigModel]->Fill(MVAVar,totalWeight*0.50*0.98);
	                          histo_Higgs_CMS_MVALepEffEBoundingDown[sigModel]->Fill(MVAVar,totalWeight*0.50*0.98);}
           histo_Higgs_CMS_PUBoundingUp  [sigModel]->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
           histo_Higgs_CMS_PUBoundingDown[sigModel]->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
           histo_Higgs_CMS_MVABTAGBoundingUp  [sigModel]->Fill(MVAVar,totalWeight*btagCorr[0]);
           histo_Higgs_CMS_MVABTAGBoundingDown[sigModel]->Fill(MVAVar,totalWeight*btagCorr[1]);
        }
        if(passSystCuts[JESUP])  histo_Higgs_CMS_MVAJESBoundingUp  [sigModel]->Fill(MVAVarJESSyst[0],totalWeight);
        if(passSystCuts[JESDOWN])histo_Higgs_CMS_MVAJESBoundingDown[sigModel]->Fill(MVAVarJESSyst[1],totalWeight);
        if(passSystCuts[JERUP])  histo_Higgs_CMS_MVAJERBoundingUp  [sigModel]->Fill(MVAVarJERSyst[0],totalWeight);
        if(passSystCuts[JERDOWN])histo_Higgs_CMS_MVAJERBoundingDown[sigModel]->Fill(MVAVarJERSyst[1],totalWeight);
        if(passSystCuts[METUP])  histo_Higgs_CMS_MVAMETBoundingUp  [sigModel]->Fill(MVAVar,totalWeight);
        if(passSystCuts[METDOWN])histo_Higgs_CMS_MVAMETBoundingDown[sigModel]->Fill(MVAVar,totalWeight);
      }
    }
    printf("                          mm+          ee+          em+          mm-          ee-          em-          all\n");
    for(int nc=0; nc<numberCuts; nc++){
      printf("(%15s): %10.2f   %10.2f   %10.2f   %10.2f   %10.2f   %10.2f   %10.2f\n",cutName[nc].Data(),hDWWLL[0]->GetBinContent(nc+1),hDWWLL[1]->GetBinContent(nc+1),hDWWLL[2]->GetBinContent(nc+1),
                                                                                                hDWWLL[3]->GetBinContent(nc+1),hDWWLL[4]->GetBinContent(nc+1),hDWWLL[5]->GetBinContent(nc+1),
												hDWWLL[6]->GetBinContent(nc+1));
    }
    the_input_file->Close();
  } // end of chain

  printf("----------------------totalFakeDataCount--------------------------------\n");
  for(int ni=0; ni<6; ni++) {
    printf("(%d): ",ni);
    for(int nj=0; nj<5; nj++) printf("%6.1f ",totalFakeDataCount[ni][nj]);
    printf("\n");
  }

  // WZ scale factor from data
  double sfE_WZ[nBinWZMVA]; for(int i=0; i<nBinWZMVA; i++) sfE_WZ[i] = 1.0;;
  printf("WZ SFs applied?: %d\n",useWZFromData);
  printf("WZini:");
  for(int np=1; np<=histo_WZ->GetNbinsX(); np++) {
    printf(" %.2f +/- %.2f",histo_WZ->GetBinContent(np),histo_WZ->GetBinError(np));
  }
  printf("\n");
  for(int nb=1; nb<=histo[6][allPlots-4][0]->GetNbinsX(); nb++) {
    double sumEventsType[3] = {0,0,0}; double sumEventsTypeE[3] = {0,0,0};
    for(int np=0; np<histBins; np++) {
      if     (np==0){
        sumEventsType[0] = sumEventsType[0] + histo[6][allPlots-4][np]->GetBinContent(nb); sumEventsTypeE[0] = sumEventsTypeE[0] + histo[6][allPlots-4][np]->GetBinError(nb)*histo[6][allPlots-4][np]->GetBinError((nb));
      }
      else if(np==3){
        sumEventsType[2] = sumEventsType[2] + histo[6][allPlots-4][np]->GetBinContent(nb); sumEventsTypeE[2] = sumEventsTypeE[2] + histo[6][allPlots-4][np]->GetBinError(nb)*histo[6][allPlots-4][np]->GetBinError((nb));
      }
      else {
        sumEventsType[1] = sumEventsType[1] + histo[6][allPlots-4][np]->GetBinContent(nb); sumEventsTypeE[1] = sumEventsTypeE[1] + histo[6][allPlots-4][np]->GetBinError(nb)*histo[6][allPlots-4][np]->GetBinError((nb));
      }
    }
    double sf_WZ  = (sumEventsType[0]-sumEventsType[1])/sumEventsType[2]; // 0 = data, 1=bg, 2=wz
    double sf_WZE[3] = {sumEventsTypeE[0]/sumEventsType[2]/sumEventsType[2],
  			sumEventsTypeE[1]/sumEventsType[2]/sumEventsType[2],
  			sumEventsTypeE[2]/sumEventsType[2]/sumEventsType[2]*sf_WZ*sf_WZ};
    sfE_WZ[nb-1] = sqrt(sf_WZE[0]+sf_WZE[1]+sf_WZE[2])/sf_WZ;
    printf("(WZ study(%d)): data: %8.2f +/- %6.2f | bkg: %8.2f +/- %6.2f | wz: %8.2f +/- %6.2f -> sf = %4.2f +/- %4.2f (%4.2f/%4.2f/%4.2f)\n",nb-1,
  	   sumEventsType[0],sqrt(sumEventsTypeE[0]),sumEventsType[1],sqrt(sumEventsTypeE[1]),sumEventsType[2],sqrt(sumEventsTypeE[2]),
           sf_WZ,sfE_WZ[nb-1]*sf_WZ,sqrt(sf_WZE[0]),sqrt(sf_WZE[1]),sqrt(sf_WZE[2]));
    if(useWZFromData){
      for(int nt=0; nt<6; nt++) {
        histo_WZ  		       ->SetBinContent(nb+nBinWZMVA*nt,histo_WZ			      ->GetBinContent(nb+nBinWZMVA*nt)*sf_WZ);
        histo_WZ  		       ->SetBinError  (nb+nBinWZMVA*nt,histo_WZ			      ->GetBinError  (nb+nBinWZMVA*nt)*sf_WZ);
        histo_WZ_CMS_MVAJESBoundingUp  ->SetBinContent(nb+nBinWZMVA*nt,histo_WZ_CMS_MVAJESBoundingUp  ->GetBinContent(nb+nBinWZMVA*nt)*sf_WZ);
        histo_WZ_CMS_MVAJESBoundingDown->SetBinError  (nb+nBinWZMVA*nt,histo_WZ_CMS_MVAJESBoundingDown->GetBinError  (nb+nBinWZMVA*nt)*sf_WZ);
        histo_WZ_CMS_MVAJERBoundingUp  ->SetBinContent(nb+nBinWZMVA*nt,histo_WZ_CMS_MVAJERBoundingUp  ->GetBinContent(nb+nBinWZMVA*nt)*sf_WZ);
        histo_WZ_CMS_MVAJERBoundingDown->SetBinError  (nb+nBinWZMVA*nt,histo_WZ_CMS_MVAJERBoundingDown->GetBinError  (nb+nBinWZMVA*nt)*sf_WZ);
        histo_WZ_CMS_MVAMETBoundingUp  ->SetBinContent(nb+nBinWZMVA*nt,histo_WZ_CMS_MVAMETBoundingUp  ->GetBinContent(nb+nBinWZMVA*nt)*sf_WZ);
        histo_WZ_CMS_MVAMETBoundingDown->SetBinError  (nb+nBinWZMVA*nt,histo_WZ_CMS_MVAMETBoundingDown->GetBinError  (nb+nBinWZMVA*nt)*sf_WZ);
      }
    }
  }
  printf("WZend:");
  for(int np=1; np<=histo_WZ->GetNbinsX(); np++) {
    printf(" %.2f +/- %.2f",histo_WZ->GetBinContent(np),histo_WZ->GetBinError(np));
  }
  printf("\n");

  // Z->ll scale factor from data
  { // 0 = data, 1=bg, 2=zll
    double sf_ZLL[2] = {1.0,1.0};
    double sumEventsType[3] = {0,0,0}; double sumEventsTypeE[3] = {0,0,0};
    for(int np=0; np<histBins; np++) {
      for(int nFinalStates=0; nFinalStates<6; nFinalStates++) {
        if     (np==0){
	  sumEventsType[0] = sumEventsType[0] + bgdDecay[nFinalStates][SSZLLSEL][np]; sumEventsTypeE[0] = sumEventsTypeE[0] + weiDecay[nFinalStates][SSZLLSEL][np];
        }
        else if(np==6){
	  sumEventsType[2] = sumEventsType[2] + bgdDecay[nFinalStates][SSZLLSEL][np]; sumEventsTypeE[2] = sumEventsTypeE[2] + weiDecay[nFinalStates][SSZLLSEL][np];
        }
        else {
	  sumEventsType[1] = sumEventsType[1] + bgdDecay[nFinalStates][SSZLLSEL][np]; sumEventsTypeE[1] = sumEventsTypeE[1] + weiDecay[nFinalStates][SSZLLSEL][np];
        }
      }
    }
    sf_ZLL[0] = (sumEventsType[0]-sumEventsType[1])/sumEventsType[2];
    double sf_ZLLE[3] = {sumEventsTypeE[0]/sumEventsType[2]/sumEventsType[2],
                         sumEventsTypeE[1]/sumEventsType[2]/sumEventsType[2],
			 sumEventsTypeE[2]/sumEventsType[2]/sumEventsType[2]*sf_ZLL[0]*sf_ZLL[0]};
    sf_ZLL[1] = sqrt(sf_ZLLE[0]+sf_ZLLE[1]+sf_ZLLE[2]);
    printf("(ZLL study): data: %8.2f +/- %6.2f | bkg: %8.2f +/- %6.2f | zll: %8.2f +/- %6.2f -> sf = %4.2f +/- %4.2f (%4.2f/%4.2f/%4.2f)\n",
    	   sumEventsType[0],sqrt(sumEventsTypeE[0]),sumEventsType[1],sqrt(sumEventsTypeE[1]),sumEventsType[2],sqrt(sumEventsTypeE[2]),
	   sf_ZLL[0],sf_ZLL[1],sqrt(sf_ZLLE[0]),sqrt(sf_ZLLE[1]),sqrt(sf_ZLLE[2]));
  }

  if(histo_EWK->GetSumOfWeights() > 0){
    printf("Signal selected Fiducial0 ==> good: %f / bad: %f -> %f\n",selectedFiducial[0],selectedNonFiducial[0],selectedNonFiducial[0]/histo_EWK->GetSumOfWeights());
    printf("Signal selected Fiducial1 ==> good: %f / bad: %f -> %f\n",selectedFiducial[1],selectedNonFiducial[1],selectedNonFiducial[1]/histo_EWK->GetSumOfWeights());
    printf("Cross section Fiducial0/1 (fb): %f %f\n",effTimesXsFiducial[0],effTimesXsFiducial[1]);
    printf("Efficiency Fiducial0/1 ==> %f %f\n",selectedFiducial[0]/effTimesXsFiducial[0]/lumi,selectedFiducial[1]/effTimesXsFiducial[1]/lumi);
  }

  histo[6][allPlots-1][0] ->Add(histo_Data);
  histo[6][allPlots-1][1] ->Add(histo_EWK);
  histo[6][allPlots-1][2] ->Add(histo_QCD);
  histo[6][allPlots-1][3] ->Add(histo_WZ);
  histo[6][allPlots-1][4] ->Add(histo_ZZ);
  histo[6][allPlots-1][5] ->Add(histo_VVV);
  histo[6][allPlots-1][6] ->Add(histo_WS);
  histo[6][allPlots-1][7] ->Add(histo_WG);
  histo[6][allPlots-1][8] ->Add(histo_DPS);
  histo[6][allPlots-1][9] ->Add(histo_FakeM);
  histo[6][allPlots-1][10]->Add(histo_FakeE);
  if(nSigModels >= 1)
  histo[6][allPlots-1][11]->Add(histo_Higgs[0]);
  if(nSigModels >= 2)
  histo[6][allPlots-1][12]->Add(histo_Higgs[1]);

  for(int nModel=0; nModel<7; nModel++){
    for(int thePlot=0; thePlot<allPlots; thePlot++){
      char output[200];
      sprintf(output,"histowwss%d_%d_%d.root",finalVar,nModel,thePlot);	  
      TFile* outFilePlotsNote = new TFile(output,"recreate");
      outFilePlotsNote->cd();
      for(int np=0; np<histBins; np++) histo[nModel][thePlot][np]->Write();
      outFilePlotsNote->Close();
    }
  }

  for(int ns=0; ns<nSelTypes; ns++) {
    printf("Selection: %s\n",selTypeName[ns].Data());
    double sumEventsType[7] = {0,0,0,0,0,0,0}; double sumEventsTypeE[7] = {0,0,0,0,0,0,0};
    for(int np=0; np<histBins; np++) {
      double tot = 0; double totE = 0;
      for(int nModel=0; nModel<6; nModel++) {
        tot = tot + bgdDecay[nModel][ns][np]; totE = totE + weiDecay[nModel][ns][np];
        if(np!=0){
          sumEventsType[nModel] = sumEventsType[nModel] + bgdDecay[nModel][ns][np]; sumEventsTypeE[nModel] = sumEventsTypeE[nModel] + weiDecay[nModel][ns][np];
	  sumEventsType[     6] = sumEventsType[     6] + bgdDecay[nModel][ns][np]; sumEventsTypeE[     6] = sumEventsTypeE[     6] + weiDecay[nModel][ns][np];
        }
      }
      if(!processName[np].Contains("Hig"))
      printf("(%5s): %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f => %8.2f +/- %6.2f\n",
      processName[np].Data(),bgdDecay[0][ns][np],sqrt(weiDecay[0][ns][np]),bgdDecay[1][ns][np],sqrt(weiDecay[1][ns][np]),
        		     bgdDecay[2][ns][np],sqrt(weiDecay[2][ns][np]),bgdDecay[3][ns][np],sqrt(weiDecay[3][ns][np]),
        		     bgdDecay[4][ns][np],sqrt(weiDecay[4][ns][np]),bgdDecay[5][ns][np],sqrt(weiDecay[5][ns][np]),
			     tot,sqrt(totE));
    }
    printf("(..bkg): %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %8.2f +/- %6.2f => %8.2f +/- %6.2f\n",
    	   sumEventsType[0],sqrt(sumEventsTypeE[0]),sumEventsType[1],sqrt(sumEventsTypeE[1]),
           sumEventsType[2],sqrt(sumEventsTypeE[2]),sumEventsType[3],sqrt(sumEventsTypeE[3]),
           sumEventsType[4],sqrt(sumEventsTypeE[4]),sumEventsType[5],sqrt(sumEventsTypeE[5]),
	   sumEventsType[6],sqrt(sumEventsTypeE[6]));
  }

  printf("QCD Corr: EWK(%f:%f/%f/%f/%f/%f/%f) QCD(%f:%f/%f/%f/%f/%f/%f) WZ(%f:%f/%f/%f/%f/%f/%f) ZZ(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) WS(%f:%f/%f/%f/%f/%f/%f) WG(%f:%f/%f/%f/%f/%f/%f) DPS(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_EWK->GetSumOfWeights(),histo_EWK_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_EWK_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_EWK_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_EWK_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_EWK_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_EWK_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_QCD->GetSumOfWeights(),histo_QCD_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_QCD_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_QCD_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_QCD_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_QCD_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_QCD_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_WZ ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[0] ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[1] ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[2] ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[3] ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[4] ->GetSumOfWeights(),histo_WZ_CMS_QCDScaleBounding[5] ->GetSumOfWeights(),
    histo_ZZ ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[0] ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[1] ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[2] ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[3] ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[4] ->GetSumOfWeights(),histo_ZZ_CMS_QCDScaleBounding[5] ->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_WS ->GetSumOfWeights(),histo_WS_CMS_QCDScaleBounding[0] ->GetSumOfWeights(),histo_WS_CMS_QCDScaleBounding[1] ->GetSumOfWeights(),histo_WS_CMS_QCDScaleBounding[2] ->GetSumOfWeights(),histo_WS_CMS_QCDScaleBounding[3] ->GetSumOfWeights(),histo_WS_CMS_QCDScaleBounding[4] ->GetSumOfWeights(),histo_WS_CMS_QCDScaleBounding[5] ->GetSumOfWeights(),
    histo_WG ->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[0] ->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[1] ->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[2] ->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[3] ->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[4] ->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[5] ->GetSumOfWeights(),
    histo_DPS->GetSumOfWeights(),histo_DPS_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_DPS_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_DPS_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_DPS_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_DPS_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_DPS_CMS_QCDScaleBounding[5]->GetSumOfWeights());

  for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {
    double factorUp = +1.0; double factorDown = -1.0;
    histo_EWK_CMS_MVAEWKStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_EWK  ->GetBinContent(i)+factorUp  *histo_EWK  ->GetBinError(i),0.000001));
    histo_EWK_CMS_MVAEWKStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_EWK  ->GetBinContent(i)+factorDown*histo_EWK  ->GetBinError(i),0.000001));
    histo_QCD_CMS_MVAQCDStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_QCD  ->GetBinContent(i)+factorUp  *histo_QCD  ->GetBinError(i),0.000001));
    histo_QCD_CMS_MVAQCDStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_QCD  ->GetBinContent(i)+factorDown*histo_QCD  ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_WZ   ->GetBinContent(i)+factorUp  *histo_WZ   ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingDown	    ->SetBinContent(i,TMath::Max(histo_WZ   ->GetBinContent(i)+factorDown*histo_WZ   ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_ZZ   ->GetBinContent(i)+factorUp  *histo_ZZ   ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingDown	    ->SetBinContent(i,TMath::Max(histo_ZZ   ->GetBinContent(i)+factorDown*histo_ZZ   ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_VVV  ->GetBinContent(i)+factorUp  *histo_VVV  ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_VVV  ->GetBinContent(i)+factorDown*histo_VVV  ->GetBinError(i),0.000001));
    histo_WS_CMS_MVAWSStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_WS   ->GetBinContent(i)+factorUp  *histo_WS   ->GetBinError(i),0.000001));
    histo_WS_CMS_MVAWSStatBoundingDown	    ->SetBinContent(i,TMath::Max(histo_WS   ->GetBinContent(i)+factorDown*histo_WS   ->GetBinError(i),0.000001));
    histo_WG_CMS_MVAWGStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_WG   ->GetBinContent(i)+factorUp  *histo_WG   ->GetBinError(i),0.000001));
    histo_WG_CMS_MVAWGStatBoundingDown	    ->SetBinContent(i,TMath::Max(histo_WG   ->GetBinContent(i)+factorDown*histo_WG   ->GetBinError(i),0.000001));
    histo_DPS_CMS_MVADPSStatBoundingUp	    ->SetBinContent(i,TMath::Max(histo_DPS  ->GetBinContent(i)+factorUp  *histo_DPS  ->GetBinError(i),0.000001));
    histo_DPS_CMS_MVADPSStatBoundingDown    ->SetBinContent(i,TMath::Max(histo_DPS  ->GetBinContent(i)+factorDown*histo_DPS  ->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_FakeM->GetBinContent(i)+factorUp  *histo_FakeM->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingDown->SetBinContent(i,TMath::Max(histo_FakeM->GetBinContent(i)+factorDown*histo_FakeM->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->SetBinContent(i,TMath::Max(histo_FakeE->GetBinContent(i)+factorUp  *histo_FakeE->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingDown->SetBinContent(i,TMath::Max(histo_FakeE->GetBinContent(i)+factorDown*histo_FakeE->GetBinError(i),0.000001));

    histo_EWK_CMS_MVAEWKStatBoundingBinUp[i-1]	    ->Add(histo_EWK  ); histo_EWK_CMS_MVAEWKStatBoundingBinUp[i-1]	->SetBinContent(i,TMath::Max(histo_EWK  ->GetBinContent(i)+factorUp  *histo_EWK  ->GetBinError(i),0.000001));
    histo_EWK_CMS_MVAEWKStatBoundingBinDown[i-1]    ->Add(histo_EWK  ); histo_EWK_CMS_MVAEWKStatBoundingBinDown[i-1]	->SetBinContent(i,TMath::Max(histo_EWK  ->GetBinContent(i)+factorDown*histo_EWK  ->GetBinError(i),0.000001));
    histo_QCD_CMS_MVAQCDStatBoundingBinUp[i-1]	    ->Add(histo_QCD  ); histo_QCD_CMS_MVAQCDStatBoundingBinUp[i-1]	->SetBinContent(i,TMath::Max(histo_QCD  ->GetBinContent(i)+factorUp  *histo_QCD  ->GetBinError(i),0.000001));
    histo_QCD_CMS_MVAQCDStatBoundingBinDown[i-1]    ->Add(histo_QCD  ); histo_QCD_CMS_MVAQCDStatBoundingBinDown[i-1]	->SetBinContent(i,TMath::Max(histo_QCD  ->GetBinContent(i)+factorDown*histo_QCD  ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]	    ->Add(histo_WZ   ); histo_WZ_CMS_MVAWZStatBoundingBinUp[i-1]	->SetBinContent(i,TMath::Max(histo_WZ	->GetBinContent(i)+factorUp  *histo_WZ   ->GetBinError(i),0.000001));
    histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]	    ->Add(histo_WZ   ); histo_WZ_CMS_MVAWZStatBoundingBinDown[i-1]	->SetBinContent(i,TMath::Max(histo_WZ	->GetBinContent(i)+factorDown*histo_WZ   ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]	    ->Add(histo_ZZ   ); histo_ZZ_CMS_MVAZZStatBoundingBinUp[i-1]	->SetBinContent(i,TMath::Max(histo_ZZ	->GetBinContent(i)+factorUp  *histo_ZZ   ->GetBinError(i),0.000001));
    histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]	    ->Add(histo_ZZ   ); histo_ZZ_CMS_MVAZZStatBoundingBinDown[i-1]	->SetBinContent(i,TMath::Max(histo_ZZ	->GetBinContent(i)+factorDown*histo_ZZ   ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	    ->Add(histo_VVV  ); histo_VVV_CMS_MVAVVVStatBoundingBinUp[i-1]	->SetBinContent(i,TMath::Max(histo_VVV  ->GetBinContent(i)+factorUp  *histo_VVV  ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]    ->Add(histo_VVV  ); histo_VVV_CMS_MVAVVVStatBoundingBinDown[i-1]	->SetBinContent(i,TMath::Max(histo_VVV  ->GetBinContent(i)+factorDown*histo_VVV  ->GetBinError(i),0.000001));
    histo_WS_CMS_MVAWSStatBoundingBinUp[i-1]	    ->Add(histo_WS   ); histo_WS_CMS_MVAWSStatBoundingBinUp[i-1]	->SetBinContent(i,TMath::Max(histo_WS	->GetBinContent(i)+factorUp  *histo_WS   ->GetBinError(i),0.000001));
    histo_WS_CMS_MVAWSStatBoundingBinDown[i-1]	    ->Add(histo_WS   ); histo_WS_CMS_MVAWSStatBoundingBinDown[i-1]	->SetBinContent(i,TMath::Max(histo_WS	->GetBinContent(i)+factorDown*histo_WS   ->GetBinError(i),0.000001));
    histo_WG_CMS_MVAWGStatBoundingBinUp[i-1]	    ->Add(histo_WG   ); histo_WG_CMS_MVAWGStatBoundingBinUp[i-1]	->SetBinContent(i,TMath::Max(histo_WG	->GetBinContent(i)+factorUp  *histo_WG   ->GetBinError(i),0.000001));
    histo_WG_CMS_MVAWGStatBoundingBinDown[i-1]	    ->Add(histo_WG   ); histo_WG_CMS_MVAWGStatBoundingBinDown[i-1]	->SetBinContent(i,TMath::Max(histo_WG	->GetBinContent(i)+factorDown*histo_WG   ->GetBinError(i),0.000001));
    histo_DPS_CMS_MVADPSStatBoundingBinUp[i-1]	    ->Add(histo_DPS  ); histo_DPS_CMS_MVADPSStatBoundingBinUp[i-1]	->SetBinContent(i,TMath::Max(histo_DPS  ->GetBinContent(i)+factorUp  *histo_DPS  ->GetBinError(i),0.000001));
    histo_DPS_CMS_MVADPSStatBoundingBinDown[i-1]    ->Add(histo_DPS  ); histo_DPS_CMS_MVADPSStatBoundingBinDown[i-1]	->SetBinContent(i,TMath::Max(histo_DPS  ->GetBinContent(i)+factorDown*histo_DPS  ->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[i-1]  ->Add(histo_FakeM); histo_FakeM_CMS_MVAFakeMStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_FakeM->GetBinContent(i)+factorUp  *histo_FakeM->GetBinError(i),0.000001));
    histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[i-1]->Add(histo_FakeM); histo_FakeM_CMS_MVAFakeMStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_FakeM->GetBinContent(i)+factorDown*histo_FakeM->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[i-1]  ->Add(histo_FakeE); histo_FakeE_CMS_MVAFakeEStatBoundingBinUp[i-1]  ->SetBinContent(i,TMath::Max(histo_FakeE->GetBinContent(i)+factorUp  *histo_FakeE->GetBinError(i),0.000001));
    histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[i-1]->Add(histo_FakeE); histo_FakeE_CMS_MVAFakeEStatBoundingBinDown[i-1]->SetBinContent(i,TMath::Max(histo_FakeE->GetBinContent(i)+factorDown*histo_FakeE->GetBinError(i),0.000001));

    for(int nModel=0; nModel<nSigModels; nModel++) { 
      histo_Higgs_CMS_MVAHiggsStatBoundingUp[nModel]          ->SetBinContent(i,TMath::Max(histo_Higgs[nModel]->GetBinContent(i)+factorUp  *histo_Higgs[nModel]->GetBinError(i),0.000001));
      histo_Higgs_CMS_MVAHiggsStatBoundingDown[nModel]        ->SetBinContent(i,TMath::Max(histo_Higgs[nModel]->GetBinContent(i)+factorDown*histo_Higgs[nModel]->GetBinError(i),0.000001));
      histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nModel][i-1]  ->Add(histo_Higgs[nModel]);
      histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nModel][i-1]->Add(histo_Higgs[nModel]);
      histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nModel][i-1]  ->SetBinContent(i,TMath::Max(histo_Higgs[nModel]->GetBinContent(i)+factorUp  *histo_Higgs[nModel]->GetBinError(i),0.000001));
      histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nModel][i-1]->SetBinContent(i,TMath::Max(histo_Higgs[nModel]->GetBinContent(i)+factorDown*histo_Higgs[nModel]->GetBinError(i),0.000001));
    }
  }
  for(int nModel=-1; nModel<nSigModels; nModel++) {
    TString theSignalName = "sm";
    if(nModel >= 0)  theSignalName = signalName_[nModel].Data();
    char outputLimits[200];
    sprintf(outputLimits,"wwss%d_%s_%s_input_%s.root",finalVar,finalStateName,theSignalName.Data(),ECMsb.Data());
    TFile* outFileLimits = new TFile(outputLimits,"recreate");
    outFileLimits->cd();
    
    if(verbose) {
      printf("*************** Model: %s ******************\n",theSignalName.Data());
      cout << histo_Data ->GetSumOfWeights() << " ";
      cout << histo_EWK  ->GetSumOfWeights() << " ";
      cout << histo_QCD  ->GetSumOfWeights() << " ";
      cout << histo_WZ   ->GetSumOfWeights() << " ";
      cout << histo_ZZ   ->GetSumOfWeights() << " ";
      cout << histo_VVV  ->GetSumOfWeights() << " ";
      cout << histo_WS   ->GetSumOfWeights() << " ";
      cout << histo_WG   ->GetSumOfWeights() << " ";
      cout << histo_DPS  ->GetSumOfWeights() << " ";
      cout << histo_FakeM->GetSumOfWeights() << " ";
      cout << histo_FakeE->GetSumOfWeights() << " ";
      cout << histo_Higgs[TMath::Max(nModel,0)]  ->GetSumOfWeights() << " ";
      cout << endl;
      printf("uncertainties Stat\n");
      for(int i=1; i<=histo_EWK  ->GetNbinsX(); i++) {if(histo_EWK  ->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVAEWKStatBoundingUp      ->GetBinContent(i)/histo_EWK  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_EWK  ->GetNbinsX(); i++) {if(histo_EWK  ->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVAEWKStatBoundingDown    ->GetBinContent(i)/histo_EWK  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD  ->GetNbinsX(); i++) {if(histo_QCD  ->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVAQCDStatBoundingUp      ->GetBinContent(i)/histo_QCD  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD  ->GetNbinsX(); i++) {if(histo_QCD  ->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVAQCDStatBoundingDown    ->GetBinContent(i)/histo_QCD  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ   ->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingUp	 ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ   ->GetNbinsX(); i++) {if(histo_WZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAWZStatBoundingDown      ->GetBinContent(i)/histo_WZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ   ->GetNbinsX(); i++) {if(histo_ZZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingUp	 ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ   ->GetNbinsX(); i++) {if(histo_ZZ   ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAZZStatBoundingDown      ->GetBinContent(i)/histo_ZZ   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV  ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp      ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV  ->GetNbinsX(); i++) {if(histo_VVV  ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown    ->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS   ->GetNbinsX(); i++) {if(histo_WS   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAWSStatBoundingUp	 ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS   ->GetNbinsX(); i++) {if(histo_WS   ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAWSStatBoundingDown      ->GetBinContent(i)/histo_WS   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG   ->GetNbinsX(); i++) {if(histo_WG   ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAWGStatBoundingUp	 ->GetBinContent(i)/histo_WG   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG   ->GetNbinsX(); i++) {if(histo_WG   ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAWGStatBoundingDown      ->GetBinContent(i)/histo_WG   ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS  ->GetNbinsX(); i++) {if(histo_DPS  ->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVADPSStatBoundingUp      ->GetBinContent(i)/histo_DPS  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS  ->GetNbinsX(); i++) {if(histo_DPS  ->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVADPSStatBoundingDown    ->GetBinContent(i)/histo_DPS  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_FakeM->GetNbinsX(); i++) {if(histo_FakeM->GetBinContent(i)>0)printf("%5.1f ",histo_FakeM_CMS_MVAFakeMStatBoundingUp  ->GetBinContent(i)/histo_FakeM->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_FakeM->GetNbinsX(); i++) {if(histo_FakeM->GetBinContent(i)>0)printf("%5.1f ",histo_FakeM_CMS_MVAFakeMStatBoundingDown->GetBinContent(i)/histo_FakeM->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_FakeE->GetNbinsX(); i++) {if(histo_FakeE->GetBinContent(i)>0)printf("%5.1f ",histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->GetBinContent(i)/histo_FakeE->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_FakeE->GetNbinsX(); i++) {if(histo_FakeE->GetBinContent(i)>0)printf("%5.1f ",histo_FakeE_CMS_MVAFakeEStatBoundingDown->GetBinContent(i)/histo_FakeE->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties LepEffM\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_EWK_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_EWK_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_QCD_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_QCD_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingUp   ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffMBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffMBoundingDown ->GetBinContent(i)/histo_WZ_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingUp   ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffMBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffMBoundingDown ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS_CMS_MVALepEffMBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVALepEffMBoundingUp   ->GetBinContent(i)/histo_WS_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS_CMS_MVALepEffMBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVALepEffMBoundingDown ->GetBinContent(i)/histo_WS_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffMBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffMBoundingUp   ->GetBinContent(i)/histo_WG_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffMBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffMBoundingDown ->GetBinContent(i)/histo_WG_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVALepEffMBoundingUp  ->GetBinContent(i)/histo_DPS_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS_CMS_MVALepEffMBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVALepEffMBoundingDown->GetBinContent(i)/histo_DPS_CMS_MVALepEffMBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties LepEffE\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_EWK_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_EWK_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_QCD_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_QCD_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingUp   ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ_CMS_MVALepEffEBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVALepEffEBoundingDown ->GetBinContent(i)/histo_WZ_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingUp   ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ_CMS_MVALepEffEBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVALepEffEBoundingDown ->GetBinContent(i)/histo_ZZ_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS_CMS_MVALepEffEBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVALepEffEBoundingUp   ->GetBinContent(i)/histo_WS_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS_CMS_MVALepEffEBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVALepEffEBoundingDown ->GetBinContent(i)/histo_WS_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffEBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffEBoundingUp   ->GetBinContent(i)/histo_WG_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffEBoundingAvg ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffEBoundingDown ->GetBinContent(i)/histo_WG_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVALepEffEBoundingUp  ->GetBinContent(i)/histo_DPS_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS_CMS_MVALepEffEBoundingAvg->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVALepEffEBoundingDown->GetBinContent(i)/histo_DPS_CMS_MVALepEffEBoundingAvg->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties MET\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVAMETBoundingUp  ->GetBinContent(i)/histo_EWK->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVAMETBoundingDown->GetBinContent(i)/histo_EWK->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVAMETBoundingUp  ->GetBinContent(i)/histo_QCD->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVAMETBoundingDown->GetBinContent(i)/histo_QCD->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_WZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_WZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_ZZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_ZZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingUp  ->GetBinContent(i)/histo_VVV->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingDown->GetBinContent(i)/histo_VVV->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAMETBoundingUp   ->GetBinContent(i)/histo_WG ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAMETBoundingDown ->GetBinContent(i)/histo_WG ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVAMETBoundingUp  ->GetBinContent(i)/histo_DPS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVAMETBoundingDown->GetBinContent(i)/histo_DPS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties JES\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVAJESBoundingUp  ->GetBinContent(i)/histo_EWK->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVAJESBoundingDown->GetBinContent(i)/histo_EWK->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVAJESBoundingUp  ->GetBinContent(i)/histo_QCD->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVAJESBoundingDown->GetBinContent(i)/histo_QCD->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingUp   ->GetBinContent(i)/histo_WZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJESBoundingDown ->GetBinContent(i)/histo_WZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingUp   ->GetBinContent(i)/histo_ZZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJESBoundingDown ->GetBinContent(i)/histo_ZZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp  ->GetBinContent(i)/histo_VVV->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(i)/histo_VVV->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAJESBoundingUp   ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAJESBoundingDown ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAJESBoundingUp   ->GetBinContent(i)/histo_WG ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAJESBoundingDown ->GetBinContent(i)/histo_WG ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVAJESBoundingUp  ->GetBinContent(i)/histo_DPS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVAJESBoundingDown->GetBinContent(i)/histo_DPS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties JER\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVAJERBoundingUp  ->GetBinContent(i)/histo_EWK->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVAJERBoundingDown->GetBinContent(i)/histo_EWK->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVAJERBoundingUp  ->GetBinContent(i)/histo_QCD->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVAJERBoundingDown->GetBinContent(i)/histo_QCD->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJERBoundingUp   ->GetBinContent(i)/histo_WZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVAJERBoundingDown ->GetBinContent(i)/histo_WZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJERBoundingUp   ->GetBinContent(i)/histo_ZZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVAJERBoundingDown ->GetBinContent(i)/histo_ZZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJERBoundingUp  ->GetBinContent(i)/histo_VVV->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJERBoundingDown->GetBinContent(i)/histo_VVV->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAJERBoundingUp   ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVAJERBoundingDown ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAJERBoundingUp   ->GetBinContent(i)/histo_WG ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAJERBoundingDown ->GetBinContent(i)/histo_WG ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVAJERBoundingUp  ->GetBinContent(i)/histo_DPS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVAJERBoundingDown->GetBinContent(i)/histo_DPS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties BTAG\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVABTAGBoundingUp  ->GetBinContent(i)/histo_EWK->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_MVABTAGBoundingDown->GetBinContent(i)/histo_EWK->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVABTAGBoundingUp  ->GetBinContent(i)/histo_QCD->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_MVABTAGBoundingDown->GetBinContent(i)/histo_QCD->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVABTAGBoundingUp   ->GetBinContent(i)/histo_WZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_MVABTAGBoundingDown ->GetBinContent(i)/histo_WZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVABTAGBoundingUp   ->GetBinContent(i)/histo_ZZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_MVABTAGBoundingDown ->GetBinContent(i)/histo_ZZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVABTAGBoundingUp  ->GetBinContent(i)/histo_VVV->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVABTAGBoundingDown->GetBinContent(i)/histo_VVV->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVABTAGBoundingUp   ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_MVABTAGBoundingDown ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVABTAGBoundingUp   ->GetBinContent(i)/histo_WG ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVABTAGBoundingDown ->GetBinContent(i)/histo_WG ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVABTAGBoundingUp  ->GetBinContent(i)/histo_DPS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_MVABTAGBoundingDown->GetBinContent(i)/histo_DPS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties PU\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_PUBoundingUp   ->GetBinContent(i)/histo_EWK->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_EWK->GetNbinsX(); i++) {if(histo_EWK->GetBinContent(i)>0)printf("%5.1f ",histo_EWK_CMS_PUBoundingDown ->GetBinContent(i)/histo_EWK->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_PUBoundingUp   ->GetBinContent(i)/histo_QCD->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_QCD->GetNbinsX(); i++) {if(histo_QCD->GetBinContent(i)>0)printf("%5.1f ",histo_QCD_CMS_PUBoundingDown ->GetBinContent(i)/histo_QCD->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingUp    ->GetBinContent(i)/histo_WZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WZ ->GetNbinsX(); i++) {if(histo_WZ ->GetBinContent(i)>0)printf("%5.1f ",histo_WZ_CMS_PUBoundingDown  ->GetBinContent(i)/histo_WZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingUp    ->GetBinContent(i)/histo_ZZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_ZZ ->GetNbinsX(); i++) {if(histo_ZZ ->GetBinContent(i)>0)printf("%5.1f ",histo_ZZ_CMS_PUBoundingDown  ->GetBinContent(i)/histo_ZZ ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingUp   ->GetBinContent(i)/histo_VVV->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_VVV->GetNbinsX(); i++) {if(histo_VVV->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingDown ->GetBinContent(i)/histo_VVV->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_PUBoundingUp    ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_PUBoundingDown  ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_PUBoundingUp    ->GetBinContent(i)/histo_WG ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WG ->GetNbinsX(); i++) {if(histo_WG ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_PUBoundingDown  ->GetBinContent(i)/histo_WG ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_PUBoundingUp   ->GetBinContent(i)/histo_DPS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_DPS->GetNbinsX(); i++) {if(histo_DPS->GetBinContent(i)>0)printf("%5.1f ",histo_DPS_CMS_PUBoundingDown ->GetBinContent(i)/histo_DPS->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      printf("uncertainties WSSF\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_WSSFUp          ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
      for(int i=1; i<=histo_WS ->GetNbinsX(); i++) {if(histo_WS ->GetBinContent(i)>0)printf("%5.1f ",histo_WS_CMS_WSSFDown        ->GetBinContent(i)/histo_WS ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
    } 

    histo_Data ->Write();
    histo_EWK  ->Write();
    histo_QCD  ->Write();
    histo_WZ   ->Write();
    histo_ZZ   ->Write();
    histo_VVV  ->Write();
    histo_WS   ->Write();
    histo_WG   ->Write();
    histo_DPS  ->Write();
    histo_FakeM->Write();
    histo_FakeE->Write();
    histo_Higgs[TMath::Max(nModel,0)]->Write();

    histo_EWK_CMS_MVAEWKStatBoundingUp	    ->Write();
    histo_EWK_CMS_MVAEWKStatBoundingDown    ->Write();
    histo_QCD_CMS_MVAQCDStatBoundingUp	    ->Write();
    histo_QCD_CMS_MVAQCDStatBoundingDown    ->Write();
    histo_WZ_CMS_MVAWZStatBoundingUp	    ->Write();
    histo_WZ_CMS_MVAWZStatBoundingDown	    ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingUp	    ->Write();
    histo_ZZ_CMS_MVAZZStatBoundingDown	    ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingUp	    ->Write();
    histo_VVV_CMS_MVAVVVStatBoundingDown    ->Write();
    histo_WS_CMS_MVAWSStatBoundingUp	    ->Write();
    histo_WS_CMS_MVAWSStatBoundingDown	    ->Write();
    histo_WG_CMS_MVAWGStatBoundingUp	    ->Write();
    histo_WG_CMS_MVAWGStatBoundingDown	    ->Write();
    histo_DPS_CMS_MVADPSStatBoundingUp	    ->Write();
    histo_DPS_CMS_MVADPSStatBoundingDown    ->Write();
    histo_FakeM_CMS_MVAFakeMStatBoundingUp  ->Write();
    histo_FakeM_CMS_MVAFakeMStatBoundingDown->Write();
    histo_FakeE_CMS_MVAFakeEStatBoundingUp  ->Write();
    histo_FakeE_CMS_MVAFakeEStatBoundingDown->Write();
    histo_Higgs_CMS_MVAHiggsStatBoundingUp[TMath::Max(nModel,0)]   ->Write();
    histo_Higgs_CMS_MVAHiggsStatBoundingDown[TMath::Max(nModel,0)] ->Write();
    
    histo_EWK_CMS_MVALepEffMBoundingUp  ->Write();
    histo_EWK_CMS_MVALepEffMBoundingDown->Write();
    histo_QCD_CMS_MVALepEffMBoundingUp  ->Write();
    histo_QCD_CMS_MVALepEffMBoundingDown->Write();
    histo_WZ_CMS_MVALepEffMBoundingUp	->Write();
    histo_WZ_CMS_MVALepEffMBoundingDown ->Write();
    histo_ZZ_CMS_MVALepEffMBoundingUp	->Write();
    histo_ZZ_CMS_MVALepEffMBoundingDown ->Write();
    histo_VVV_CMS_MVALepEffMBoundingUp  ->Write();
    histo_VVV_CMS_MVALepEffMBoundingDown->Write();
    histo_WS_CMS_MVALepEffMBoundingUp	->Write();
    histo_WS_CMS_MVALepEffMBoundingDown ->Write();
    histo_WG_CMS_MVALepEffMBoundingUp	->Write();
    histo_WG_CMS_MVALepEffMBoundingDown ->Write();
    histo_DPS_CMS_MVALepEffMBoundingUp  ->Write();
    histo_DPS_CMS_MVALepEffMBoundingDown->Write();
    histo_Higgs_CMS_MVALepEffMBoundingUp[TMath::Max(nModel,0)]  ->Write();
    histo_Higgs_CMS_MVALepEffMBoundingDown[TMath::Max(nModel,0)]->Write();

    histo_EWK_CMS_MVALepEffEBoundingUp  ->Write();
    histo_EWK_CMS_MVALepEffEBoundingDown->Write();
    histo_QCD_CMS_MVALepEffEBoundingUp  ->Write();
    histo_QCD_CMS_MVALepEffEBoundingDown->Write();
    histo_WZ_CMS_MVALepEffEBoundingUp	->Write();
    histo_WZ_CMS_MVALepEffEBoundingDown ->Write();
    histo_ZZ_CMS_MVALepEffEBoundingUp	->Write();
    histo_ZZ_CMS_MVALepEffEBoundingDown ->Write();
    histo_VVV_CMS_MVALepEffEBoundingUp  ->Write();
    histo_VVV_CMS_MVALepEffEBoundingDown->Write();
    histo_WS_CMS_MVALepEffEBoundingUp	->Write();
    histo_WS_CMS_MVALepEffEBoundingDown ->Write();
    histo_WG_CMS_MVALepEffEBoundingUp	->Write();
    histo_WG_CMS_MVALepEffEBoundingDown ->Write();
    histo_DPS_CMS_MVALepEffEBoundingUp  ->Write();
    histo_DPS_CMS_MVALepEffEBoundingDown->Write();
    histo_Higgs_CMS_MVALepEffEBoundingUp[TMath::Max(nModel,0)]  ->Write();
    histo_Higgs_CMS_MVALepEffEBoundingDown[TMath::Max(nModel,0)]->Write();
    
    histo_EWK_CMS_MVAMETBoundingUp  ->Write();
    histo_EWK_CMS_MVAMETBoundingDown->Write();
    histo_QCD_CMS_MVAMETBoundingUp  ->Write();
    histo_QCD_CMS_MVAMETBoundingDown->Write();
    histo_WZ_CMS_MVAMETBoundingUp   ->Write();
    histo_WZ_CMS_MVAMETBoundingDown ->Write();
    histo_ZZ_CMS_MVAMETBoundingUp   ->Write();
    histo_ZZ_CMS_MVAMETBoundingDown ->Write();
    histo_VVV_CMS_MVAMETBoundingUp  ->Write();
    histo_VVV_CMS_MVAMETBoundingDown->Write();
    histo_WS_CMS_MVAMETBoundingUp   ->Write();
    histo_WS_CMS_MVAMETBoundingDown ->Write();
    histo_WG_CMS_MVAMETBoundingUp   ->Write();
    histo_WG_CMS_MVAMETBoundingDown ->Write();
    histo_DPS_CMS_MVAMETBoundingUp  ->Write();
    histo_DPS_CMS_MVAMETBoundingDown->Write();
    histo_Higgs_CMS_MVAMETBoundingUp[TMath::Max(nModel,0)]      ->Write();
    histo_Higgs_CMS_MVAMETBoundingDown[TMath::Max(nModel,0)]    ->Write();

    histo_EWK_CMS_MVAJESBoundingUp  ->Write();
    histo_EWK_CMS_MVAJESBoundingDown->Write();
    histo_QCD_CMS_MVAJESBoundingUp  ->Write();
    histo_QCD_CMS_MVAJESBoundingDown->Write();
    histo_WZ_CMS_MVAJESBoundingUp   ->Write();
    histo_WZ_CMS_MVAJESBoundingDown ->Write();
    histo_ZZ_CMS_MVAJESBoundingUp   ->Write();
    histo_ZZ_CMS_MVAJESBoundingDown ->Write();
    histo_VVV_CMS_MVAJESBoundingUp  ->Write();
    histo_VVV_CMS_MVAJESBoundingDown->Write();
    histo_WS_CMS_MVAJESBoundingUp   ->Write();
    histo_WS_CMS_MVAJESBoundingDown ->Write();
    histo_WG_CMS_MVAJESBoundingUp   ->Write();
    histo_WG_CMS_MVAJESBoundingDown ->Write();
    histo_DPS_CMS_MVAJESBoundingUp  ->Write();
    histo_DPS_CMS_MVAJESBoundingDown->Write();
    histo_Higgs_CMS_MVAJESBoundingUp[TMath::Max(nModel,0)]  ->Write();
    histo_Higgs_CMS_MVAJESBoundingDown[TMath::Max(nModel,0)]->Write(); 


    histo_EWK_CMS_MVAJERBoundingUp  ->Write();
    histo_EWK_CMS_MVAJERBoundingDown->Write();
    histo_QCD_CMS_MVAJERBoundingUp  ->Write();
    histo_QCD_CMS_MVAJERBoundingDown->Write();
    histo_WZ_CMS_MVAJERBoundingUp   ->Write();
    histo_WZ_CMS_MVAJERBoundingDown ->Write();
    histo_ZZ_CMS_MVAJERBoundingUp   ->Write();
    histo_ZZ_CMS_MVAJERBoundingDown ->Write();
    histo_VVV_CMS_MVAJERBoundingUp  ->Write();
    histo_VVV_CMS_MVAJERBoundingDown->Write();
    histo_WS_CMS_MVAJERBoundingUp   ->Write();
    histo_WS_CMS_MVAJERBoundingDown ->Write();
    histo_WG_CMS_MVAJERBoundingUp   ->Write();
    histo_WG_CMS_MVAJERBoundingDown ->Write();
    histo_DPS_CMS_MVAJERBoundingUp  ->Write();
    histo_DPS_CMS_MVAJERBoundingDown->Write();
    histo_Higgs_CMS_MVAJERBoundingUp[TMath::Max(nModel,0)]  ->Write();
    histo_Higgs_CMS_MVAJERBoundingDown[TMath::Max(nModel,0)]->Write(); 

    histo_EWK_CMS_MVABTAGBoundingUp  ->Write();
    histo_EWK_CMS_MVABTAGBoundingDown->Write();
    histo_QCD_CMS_MVABTAGBoundingUp  ->Write();
    histo_QCD_CMS_MVABTAGBoundingDown->Write();
    histo_WZ_CMS_MVABTAGBoundingUp   ->Write();
    histo_WZ_CMS_MVABTAGBoundingDown ->Write();
    histo_ZZ_CMS_MVABTAGBoundingUp   ->Write();
    histo_ZZ_CMS_MVABTAGBoundingDown ->Write();
    histo_VVV_CMS_MVABTAGBoundingUp  ->Write();
    histo_VVV_CMS_MVABTAGBoundingDown->Write();
    histo_WS_CMS_MVABTAGBoundingUp   ->Write();
    histo_WS_CMS_MVABTAGBoundingDown ->Write();
    histo_WG_CMS_MVABTAGBoundingUp   ->Write();
    histo_WG_CMS_MVABTAGBoundingDown ->Write();
    histo_DPS_CMS_MVABTAGBoundingUp  ->Write();
    histo_DPS_CMS_MVABTAGBoundingDown->Write();
    histo_Higgs_CMS_MVABTAGBoundingUp[TMath::Max(nModel,0)]  ->Write();
    histo_Higgs_CMS_MVABTAGBoundingDown[TMath::Max(nModel,0)]->Write(); 

    histo_EWK_CMS_PUBoundingUp  ->Write();
    histo_EWK_CMS_PUBoundingDown->Write();
    histo_QCD_CMS_PUBoundingUp  ->Write();
    histo_QCD_CMS_PUBoundingDown->Write();
    histo_WZ_CMS_PUBoundingUp   ->Write();
    histo_WZ_CMS_PUBoundingDown ->Write();
    histo_ZZ_CMS_PUBoundingUp   ->Write();
    histo_ZZ_CMS_PUBoundingDown ->Write();
    histo_VVV_CMS_PUBoundingUp  ->Write();
    histo_VVV_CMS_PUBoundingDown->Write();
    histo_WS_CMS_PUBoundingUp   ->Write();
    histo_WS_CMS_PUBoundingDown ->Write();
    histo_WG_CMS_PUBoundingUp   ->Write();
    histo_WG_CMS_PUBoundingDown ->Write();
    histo_DPS_CMS_PUBoundingUp  ->Write();
    histo_DPS_CMS_PUBoundingDown->Write();
    histo_Higgs_CMS_PUBoundingUp[TMath::Max(nModel,0)]  ->Write();
    histo_Higgs_CMS_PUBoundingDown[TMath::Max(nModel,0)]->Write();

    histo_WS_CMS_WSSFUp         ->Write();
    histo_WS_CMS_WSSFDown       ->Write();

    outFileLimits->Close();

    double lumiE = 1.026;
    double systLepResE[9] = {1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01};
    double systLepResM[9] = {1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01};
 
    for(int nb=1; nb<=nBinMVA; nb++){
      // QCD study
      double systQCDScale[9] = {TMath::Abs(histo_EWK_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_EWK->GetBinContent(nb)),
                                TMath::Abs(histo_QCD_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_QCD->GetBinContent(nb)),
                                TMath::Abs(histo_WZ_CMS_QCDScaleBounding[0] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb)),
                                TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[0] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb)),
                                TMath::Abs(histo_VVV_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_VVV->GetBinContent(nb)),
                                TMath::Abs(histo_WS_CMS_QCDScaleBounding[0] ->GetBinContent(nb)-histo_WS ->GetBinContent(nb)),
                                TMath::Abs(histo_WG_CMS_QCDScaleBounding[0] ->GetBinContent(nb)-histo_WG ->GetBinContent(nb)),
                                TMath::Abs(histo_DPS_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_DPS->GetBinContent(nb)),
                                TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[TMath::Max(nModel,0)][0]->GetBinContent(nb)-histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb))};
      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(histo_EWK_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_EWK->GetBinContent(nb)) > systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_EWK_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_EWK->GetBinContent(nb));
        if(TMath::Abs(histo_QCD_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_QCD->GetBinContent(nb)) > systQCDScale[1]) systQCDScale[1] = TMath::Abs(histo_QCD_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_QCD->GetBinContent(nb));
        if(TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb)) > systQCDScale[2]) systQCDScale[2] = TMath::Abs(histo_WZ_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_WZ ->GetBinContent(nb));
        if(TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb)) > systQCDScale[3]) systQCDScale[3] = TMath::Abs(histo_ZZ_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_ZZ ->GetBinContent(nb));
        if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb)) > systQCDScale[4]) systQCDScale[4] = TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb));
        if(TMath::Abs(histo_WS_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_WS ->GetBinContent(nb)) > systQCDScale[5]) systQCDScale[5] = TMath::Abs(histo_WS_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_WS ->GetBinContent(nb));
        if(TMath::Abs(histo_WG_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_WG ->GetBinContent(nb)) > systQCDScale[6]) systQCDScale[6] = TMath::Abs(histo_WG_CMS_QCDScaleBounding[nqcd] ->GetBinContent(nb)-histo_WG ->GetBinContent(nb));
        if(TMath::Abs(histo_DPS_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_DPS->GetBinContent(nb)) > systQCDScale[7]) systQCDScale[7] = TMath::Abs(histo_DPS_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_DPS->GetBinContent(nb));
        if(TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[TMath::Max(nModel,0)][nqcd]->GetBinContent(nb)  -histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb)) > systQCDScale[8]) systQCDScale[8] = TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[TMath::Max(nModel,0)][nqcd]->GetBinContent(nb)  -histo_Higgs[TMath::Max(nModel,0)]  ->GetBinContent(nb));
      }
      if(histo_EWK->GetBinContent(nb) > 0) systQCDScale[0] = 1 + systQCDScale[0]/histo_EWK->GetBinContent(nb); else systQCDScale[0] = 1;
      if(histo_QCD->GetBinContent(nb) > 0) systQCDScale[1] = 1 + systQCDScale[1]/histo_QCD->GetBinContent(nb); else systQCDScale[1] = 1;
      if(histo_WZ ->GetBinContent(nb) > 0) systQCDScale[2] = 1 + systQCDScale[2]/histo_WZ ->GetBinContent(nb); else systQCDScale[2] = 1;
      if(histo_ZZ ->GetBinContent(nb) > 0) systQCDScale[3] = 1 + systQCDScale[3]/histo_ZZ ->GetBinContent(nb); else systQCDScale[3] = 1;
      if(histo_VVV->GetBinContent(nb) > 0) systQCDScale[4] = 1 + systQCDScale[4]/histo_VVV->GetBinContent(nb); else systQCDScale[4] = 1;
      if(histo_WS ->GetBinContent(nb) > 0) systQCDScale[5] = 1 + systQCDScale[5]/histo_WS ->GetBinContent(nb); else systQCDScale[5] = 1;
      if(histo_WG ->GetBinContent(nb) > 0) systQCDScale[6] = 1 + systQCDScale[6]/histo_WG ->GetBinContent(nb); else systQCDScale[6] = 1;
      if(histo_DPS->GetBinContent(nb) > 0) systQCDScale[7] = 1 + systQCDScale[7]/histo_DPS->GetBinContent(nb); else systQCDScale[7] = 1;
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systQCDScale[8] = 1 + systQCDScale[8]/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb); else systQCDScale[0] = 1;

      for(int ntype=0; ntype<9; ntype++) if(systQCDScale[ntype] < 0) systQCDScale[ntype] = 1.0;
      if(verbose) printf("QCDScale(%d): %f %f %f %f %f %f %f %f %f\n",nb,systQCDScale[0],systQCDScale[1],systQCDScale[2],systQCDScale[3],systQCDScale[4],systQCDScale[5],systQCDScale[6],systQCDScale[7],systQCDScale[8]);
      for(int ntype=0; ntype<9; ntype++) if(systQCDScale[ntype] > 1.5) systQCDScale[ntype] = 1.5;
  
      // PDF study
      double systPDF[9];
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_EWK_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_EWK->GetBinContent(nb))/histo_EWK->GetBinContent(nb));
      systPDF[0] = 1.0+histo_Diff->GetRMS();
 
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_QCD_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_QCD->GetBinContent(nb))/histo_QCD->GetBinContent(nb));
      systPDF[1] = 1.0+histo_Diff->GetRMS();
 
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_WZ_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_WZ->GetBinContent(nb))/histo_WZ->GetBinContent(nb));
      systPDF[2] = 1.0+histo_Diff->GetRMS();
 
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_ZZ_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_ZZ->GetBinContent(nb))/histo_ZZ->GetBinContent(nb));
      systPDF[3] = 1.0+histo_Diff->GetRMS();
 
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_VVV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))/histo_VVV->GetBinContent(nb));
      systPDF[4] = 1.0+histo_Diff->GetRMS();
 
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_WS_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_WS->GetBinContent(nb))/histo_WS->GetBinContent(nb));
      systPDF[5] = 1.0+histo_Diff->GetRMS();
 
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_WG_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_WG->GetBinContent(nb))/histo_WG->GetBinContent(nb));
      systPDF[6] = 1.0+histo_Diff->GetRMS();
 
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_DPS_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_DPS->GetBinContent(nb))/histo_DPS->GetBinContent(nb));
      systPDF[7] = 1.0+histo_Diff->GetRMS();
   
      histo_Diff->Reset();
      for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_Higgs_CMS_PDFBounding[TMath::Max(nModel,0)][npdf]->GetBinContent(nb)-histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb))/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb));
      systPDF[8] = 1.0+histo_Diff->GetRMS();

      if(verbose) printf("PDF(%d): %f %f %f %f %f %f %f %f %f\n",nb,systPDF[0],systPDF[1],systPDF[2],systPDF[3],systPDF[4],systPDF[5],systPDF[6],systPDF[7],systPDF[8]);
      for(int i=0; i<9; i++) if(systPDF[i] <= 1.0) systPDF[i] = 1.01;

      double systLepEffM[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      if     (histo_EWK_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_EWK_CMS_MVALepEffMBoundingUp  ->GetBinContent(nb) > 0) systLepEffM[0] = histo_EWK_CMS_MVALepEffMBoundingUp ->GetBinContent(nb)/histo_EWK_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb);
      else if(histo_EWK_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_EWK_CMS_MVALepEffMBoundingDown->GetBinContent(nb) > 0) systLepEffM[0] = histo_EWK_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_EWK_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_QCD_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_QCD_CMS_MVALepEffMBoundingUp  ->GetBinContent(nb) > 0) systLepEffM[1] = histo_QCD_CMS_MVALepEffMBoundingUp ->GetBinContent(nb)/histo_QCD_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb);
      else if(histo_QCD_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_QCD_CMS_MVALepEffMBoundingDown->GetBinContent(nb) > 0) systLepEffM[1] = histo_QCD_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_QCD_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_WZ_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVALepEffMBoundingUp	->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingUp  ->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingAvg  ->GetBinContent(nb);
      else if(histo_WZ_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[2] = histo_WZ_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb)/histo_WZ_CMS_MVALepEffMBoundingDown ->GetBinContent(nb);
      if     (histo_ZZ_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVALepEffMBoundingUp	->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingUp  ->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingAvg  ->GetBinContent(nb);
      else if(histo_ZZ_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[3] = histo_ZZ_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffMBoundingDown ->GetBinContent(nb);
      if     (histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffMBoundingUp  ->GetBinContent(nb) > 0) systLepEffM[4] = histo_VVV_CMS_MVALepEffMBoundingUp ->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb);
      else if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb) > 0) systLepEffM[4] = histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_WS_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb) > 0 && histo_WS_CMS_MVALepEffMBoundingUp	->GetBinContent(nb) > 0) systLepEffM[5] = histo_WS_CMS_MVALepEffMBoundingUp  ->GetBinContent(nb)/histo_WS_CMS_MVALepEffMBoundingAvg  ->GetBinContent(nb);
      else if(histo_WS_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb) > 0 && histo_WS_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[5] = histo_WS_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb)/histo_WS_CMS_MVALepEffMBoundingDown ->GetBinContent(nb);
      if     (histo_WG_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffMBoundingUp	->GetBinContent(nb) > 0) systLepEffM[6] = histo_WG_CMS_MVALepEffMBoundingUp  ->GetBinContent(nb)/histo_WG_CMS_MVALepEffMBoundingAvg  ->GetBinContent(nb);
      else if(histo_WG_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[6] = histo_WG_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb)/histo_WG_CMS_MVALepEffMBoundingDown ->GetBinContent(nb);
      if     (histo_DPS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_DPS_CMS_MVALepEffMBoundingUp  ->GetBinContent(nb) > 0) systLepEffM[7] = histo_DPS_CMS_MVALepEffMBoundingUp ->GetBinContent(nb)/histo_DPS_CMS_MVALepEffMBoundingAvg ->GetBinContent(nb);
      else if(histo_DPS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_DPS_CMS_MVALepEffMBoundingDown->GetBinContent(nb) > 0) systLepEffM[7] = histo_DPS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_DPS_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
      if     (histo_Higgs_CMS_MVALepEffMBoundingAvg[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffMBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systLepEffM[8] = histo_Higgs_CMS_MVALepEffMBoundingUp [TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffMBoundingAvg [TMath::Max(nModel,0)]->GetBinContent(nb);
      else if(histo_Higgs_CMS_MVALepEffMBoundingAvg[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffMBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systLepEffM[8] = histo_Higgs_CMS_MVALepEffMBoundingAvg[TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffMBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb);
  
      double systLepEffE[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      if     (histo_EWK_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_EWK_CMS_MVALepEffEBoundingUp  ->GetBinContent(nb) > 0) systLepEffE[0] = histo_EWK_CMS_MVALepEffEBoundingUp ->GetBinContent(nb)/histo_EWK_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb);
      else if(histo_EWK_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_EWK_CMS_MVALepEffEBoundingDown->GetBinContent(nb) > 0) systLepEffE[0] = histo_EWK_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_EWK_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
      if     (histo_QCD_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_QCD_CMS_MVALepEffEBoundingUp  ->GetBinContent(nb) > 0) systLepEffE[1] = histo_QCD_CMS_MVALepEffEBoundingUp ->GetBinContent(nb)/histo_QCD_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb);
      else if(histo_QCD_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_QCD_CMS_MVALepEffEBoundingDown->GetBinContent(nb) > 0) systLepEffE[1] = histo_QCD_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_QCD_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
      if     (histo_WZ_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVALepEffEBoundingUp	->GetBinContent(nb) > 0) systLepEffE[2] = histo_WZ_CMS_MVALepEffEBoundingUp  ->GetBinContent(nb)/histo_WZ_CMS_MVALepEffEBoundingAvg  ->GetBinContent(nb);
      else if(histo_WZ_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVALepEffEBoundingDown ->GetBinContent(nb) > 0) systLepEffE[2] = histo_WZ_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb)/histo_WZ_CMS_MVALepEffEBoundingDown ->GetBinContent(nb);
      if     (histo_ZZ_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVALepEffEBoundingUp	->GetBinContent(nb) > 0) systLepEffE[3] = histo_ZZ_CMS_MVALepEffEBoundingUp  ->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingAvg  ->GetBinContent(nb);
      else if(histo_ZZ_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVALepEffEBoundingDown ->GetBinContent(nb) > 0) systLepEffE[3] = histo_ZZ_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb)/histo_ZZ_CMS_MVALepEffEBoundingDown ->GetBinContent(nb);
      if     (histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffEBoundingUp  ->GetBinContent(nb) > 0) systLepEffE[4] = histo_VVV_CMS_MVALepEffEBoundingUp ->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb);
      else if(histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(nb) > 0) systLepEffE[4] = histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
      if     (histo_WS_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb) > 0 && histo_WS_CMS_MVALepEffEBoundingUp	->GetBinContent(nb) > 0) systLepEffE[5] = histo_WS_CMS_MVALepEffEBoundingUp  ->GetBinContent(nb)/histo_WS_CMS_MVALepEffEBoundingAvg  ->GetBinContent(nb);
      else if(histo_WS_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb) > 0 && histo_WS_CMS_MVALepEffEBoundingDown ->GetBinContent(nb) > 0) systLepEffE[5] = histo_WS_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb)/histo_WS_CMS_MVALepEffEBoundingDown ->GetBinContent(nb);
      if     (histo_WG_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffEBoundingUp	->GetBinContent(nb) > 0) systLepEffE[6] = histo_WG_CMS_MVALepEffEBoundingUp  ->GetBinContent(nb)/histo_WG_CMS_MVALepEffEBoundingAvg  ->GetBinContent(nb);
      else if(histo_WG_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffEBoundingDown ->GetBinContent(nb) > 0) systLepEffE[6] = histo_WG_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb)/histo_WG_CMS_MVALepEffEBoundingDown ->GetBinContent(nb);
      if     (histo_DPS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_DPS_CMS_MVALepEffEBoundingUp  ->GetBinContent(nb) > 0) systLepEffE[7] = histo_DPS_CMS_MVALepEffEBoundingUp ->GetBinContent(nb)/histo_DPS_CMS_MVALepEffEBoundingAvg ->GetBinContent(nb);
      else if(histo_DPS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_DPS_CMS_MVALepEffEBoundingDown->GetBinContent(nb) > 0) systLepEffE[7] = histo_DPS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_DPS_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
      if     (histo_Higgs_CMS_MVALepEffEBoundingAvg[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffEBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systLepEffE[8] = histo_Higgs_CMS_MVALepEffEBoundingUp [TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffEBoundingAvg [TMath::Max(nModel,0)]->GetBinContent(nb);
      else if(histo_Higgs_CMS_MVALepEffEBoundingAvg[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffEBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systLepEffE[8] = histo_Higgs_CMS_MVALepEffEBoundingAvg[TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffEBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb);
  
      double systMetUp  [9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      double systMetDown[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      if(histo_EWK->GetBinContent(nb) > 0 && histo_EWK_CMS_MVAMETBoundingUp  ->GetBinContent(nb) > 0) systMetUp  [0] = histo_EWK_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_EWK->GetBinContent(nb);
      if(histo_EWK->GetBinContent(nb) > 0 && histo_EWK_CMS_MVAMETBoundingDown->GetBinContent(nb) > 0) systMetDown[0] = histo_EWK_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_EWK->GetBinContent(nb);
      if(histo_QCD->GetBinContent(nb) > 0 && histo_QCD_CMS_MVAMETBoundingUp  ->GetBinContent(nb) > 0) systMetUp  [1] = histo_QCD_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_QCD->GetBinContent(nb);
      if(histo_QCD->GetBinContent(nb) > 0 && histo_QCD_CMS_MVAMETBoundingDown->GetBinContent(nb) > 0) systMetDown[1] = histo_QCD_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_QCD->GetBinContent(nb);
      if(histo_WZ ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVAMETBoundingUp   ->GetBinContent(nb) > 0) systMetUp  [2] = histo_WZ_CMS_MVAMETBoundingUp   ->GetBinContent(nb)/histo_WZ ->GetBinContent(nb);
      if(histo_WZ ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVAMETBoundingDown ->GetBinContent(nb) > 0) systMetDown[2] = histo_WZ_CMS_MVAMETBoundingDown ->GetBinContent(nb)/histo_WZ ->GetBinContent(nb);
      if(histo_ZZ ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVAMETBoundingUp   ->GetBinContent(nb) > 0) systMetUp  [3] = histo_ZZ_CMS_MVAMETBoundingUp   ->GetBinContent(nb)/histo_ZZ ->GetBinContent(nb);
      if(histo_ZZ ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVAMETBoundingDown ->GetBinContent(nb) > 0) systMetDown[3] = histo_ZZ_CMS_MVAMETBoundingDown ->GetBinContent(nb)/histo_ZZ ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb) > 0 && histo_VVV_CMS_MVAMETBoundingUp  ->GetBinContent(nb) > 0) systMetUp  [4] = histo_VVV_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb) > 0 && histo_VVV_CMS_MVAMETBoundingDown->GetBinContent(nb) > 0) systMetDown[4] = histo_VVV_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_MVAMETBoundingUp   ->GetBinContent(nb) > 0) systMetUp  [5] = histo_WS_CMS_MVAMETBoundingUp   ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_MVAMETBoundingDown ->GetBinContent(nb) > 0) systMetDown[5] = histo_WS_CMS_MVAMETBoundingDown ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WG ->GetBinContent(nb) > 0 && histo_WG_CMS_MVAMETBoundingUp   ->GetBinContent(nb) > 0) systMetUp  [6] = histo_WG_CMS_MVAMETBoundingUp   ->GetBinContent(nb)/histo_WG ->GetBinContent(nb);
      if(histo_WG ->GetBinContent(nb) > 0 && histo_WG_CMS_MVAMETBoundingDown ->GetBinContent(nb) > 0) systMetDown[6] = histo_WG_CMS_MVAMETBoundingDown ->GetBinContent(nb)/histo_WG ->GetBinContent(nb);
      if(histo_DPS->GetBinContent(nb) > 0 && histo_DPS_CMS_MVAMETBoundingUp  ->GetBinContent(nb) > 0) systMetUp  [7] = histo_DPS_CMS_MVAMETBoundingUp  ->GetBinContent(nb)/histo_DPS->GetBinContent(nb);
      if(histo_DPS->GetBinContent(nb) > 0 && histo_DPS_CMS_MVAMETBoundingDown->GetBinContent(nb) > 0) systMetDown[7] = histo_DPS_CMS_MVAMETBoundingDown->GetBinContent(nb)/histo_DPS->GetBinContent(nb);
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVAMETBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systMetUp  [8] = histo_Higgs_CMS_MVAMETBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVAMETBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systMetDown[8] = histo_Higgs_CMS_MVAMETBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);
      for(int i=0; i<9; i++) if(systMetUp  [i] == 1) systMetUp  [i] = 0.998;
      for(int i=0; i<9; i++) if(systMetDown[i] == 1) systMetDown[i] = 1.002;
      for(int nmet=0; nmet<9; nmet++) if(systMetUp[nmet]   > 1.10) systMetUp[nmet]   = 1.10;
      for(int nmet=0; nmet<9; nmet++) if(systMetUp[nmet]   < 0.90) systMetUp[nmet]   = 0.90;
      for(int nmet=0; nmet<9; nmet++) if(systMetDown[nmet] > 1.10) systMetDown[nmet] = 1.10;
      for(int nmet=0; nmet<9; nmet++) if(systMetDown[nmet] < 0.90) systMetDown[nmet] = 0.90;

      double systJesUp  [9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      double systJesDown[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      if(histo_EWK->GetBinContent(nb) > 0 && histo_EWK_CMS_MVAJESBoundingUp  ->GetBinContent(nb) > 0) systJesUp  [0] = histo_EWK_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_EWK->GetBinContent(nb);
      if(histo_EWK->GetBinContent(nb) > 0 && histo_EWK_CMS_MVAJESBoundingDown->GetBinContent(nb) > 0) systJesDown[0] = histo_EWK_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_EWK->GetBinContent(nb);
      if(histo_QCD->GetBinContent(nb) > 0 && histo_QCD_CMS_MVAJESBoundingUp  ->GetBinContent(nb) > 0) systJesUp  [1] = histo_QCD_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_QCD->GetBinContent(nb);
      if(histo_QCD->GetBinContent(nb) > 0 && histo_QCD_CMS_MVAJESBoundingDown->GetBinContent(nb) > 0) systJesDown[1] = histo_QCD_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_QCD->GetBinContent(nb);
      if(histo_WZ ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVAJESBoundingUp   ->GetBinContent(nb) > 0) systJesUp  [2] = histo_WZ_CMS_MVAJESBoundingUp   ->GetBinContent(nb)/histo_WZ ->GetBinContent(nb);
      if(histo_WZ ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVAJESBoundingDown ->GetBinContent(nb) > 0) systJesDown[2] = histo_WZ_CMS_MVAJESBoundingDown ->GetBinContent(nb)/histo_WZ ->GetBinContent(nb);
      if(histo_ZZ ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVAJESBoundingUp   ->GetBinContent(nb) > 0) systJesUp  [3] = histo_ZZ_CMS_MVAJESBoundingUp   ->GetBinContent(nb)/histo_ZZ ->GetBinContent(nb);
      if(histo_ZZ ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVAJESBoundingDown ->GetBinContent(nb) > 0) systJesDown[3] = histo_ZZ_CMS_MVAJESBoundingDown ->GetBinContent(nb)/histo_ZZ ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb) > 0 && histo_VVV_CMS_MVAJESBoundingUp  ->GetBinContent(nb) > 0) systJesUp  [4] = histo_VVV_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb) > 0 && histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(nb) > 0) systJesDown[4] = histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_MVAJESBoundingUp   ->GetBinContent(nb) > 0) systJesUp  [5] = histo_WS_CMS_MVAJESBoundingUp   ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_MVAJESBoundingDown ->GetBinContent(nb) > 0) systJesDown[5] = histo_WS_CMS_MVAJESBoundingDown ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WG ->GetBinContent(nb) > 0 && histo_WG_CMS_MVAJESBoundingUp   ->GetBinContent(nb) > 0) systJesUp  [6] = histo_WG_CMS_MVAJESBoundingUp   ->GetBinContent(nb)/histo_WG ->GetBinContent(nb);
      if(histo_WG ->GetBinContent(nb) > 0 && histo_WG_CMS_MVAJESBoundingDown ->GetBinContent(nb) > 0) systJesDown[6] = histo_WG_CMS_MVAJESBoundingDown ->GetBinContent(nb)/histo_WG ->GetBinContent(nb);
      if(histo_DPS->GetBinContent(nb) > 0 && histo_DPS_CMS_MVAJESBoundingUp  ->GetBinContent(nb) > 0) systJesUp  [7] = histo_DPS_CMS_MVAJESBoundingUp  ->GetBinContent(nb)/histo_DPS->GetBinContent(nb);
      if(histo_DPS->GetBinContent(nb) > 0 && histo_DPS_CMS_MVAJESBoundingDown->GetBinContent(nb) > 0) systJesDown[7] = histo_DPS_CMS_MVAJESBoundingDown->GetBinContent(nb)/histo_DPS->GetBinContent(nb);
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVAJESBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systJesUp  [8] = histo_Higgs_CMS_MVAJESBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVAJESBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systJesDown[8] = histo_Higgs_CMS_MVAJESBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);
      for(int njes=0; njes<9; njes++) if(systJesUp[njes]   > 1.10) systJesUp[njes]   = 1.10;
      for(int njes=0; njes<9; njes++) if(systJesUp[njes]   < 0.90) systJesUp[njes]   = 0.90;
      for(int njes=0; njes<9; njes++) if(systJesDown[njes] > 1.10) systJesDown[njes] = 1.10;
      for(int njes=0; njes<9; njes++) if(systJesDown[njes] < 0.90) systJesDown[njes] = 0.90;

      double systJerUp  [9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      double systJerDown[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      if(histo_EWK->GetBinContent(nb) > 0 && histo_EWK_CMS_MVAJERBoundingUp  ->GetBinContent(nb) > 0) systJerUp  [0] = histo_EWK_CMS_MVAJERBoundingUp  ->GetBinContent(nb)/histo_EWK->GetBinContent(nb);
      if(histo_EWK->GetBinContent(nb) > 0 && histo_EWK_CMS_MVAJERBoundingDown->GetBinContent(nb) > 0) systJerDown[0] = histo_EWK_CMS_MVAJERBoundingDown->GetBinContent(nb)/histo_EWK->GetBinContent(nb);
      if(histo_QCD->GetBinContent(nb) > 0 && histo_QCD_CMS_MVAJERBoundingUp  ->GetBinContent(nb) > 0) systJerUp  [1] = histo_QCD_CMS_MVAJERBoundingUp  ->GetBinContent(nb)/histo_QCD->GetBinContent(nb);
      if(histo_QCD->GetBinContent(nb) > 0 && histo_QCD_CMS_MVAJERBoundingDown->GetBinContent(nb) > 0) systJerDown[1] = histo_QCD_CMS_MVAJERBoundingDown->GetBinContent(nb)/histo_QCD->GetBinContent(nb);
      if(histo_WZ ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVAJERBoundingUp   ->GetBinContent(nb) > 0) systJerUp  [2] = histo_WZ_CMS_MVAJERBoundingUp   ->GetBinContent(nb)/histo_WZ ->GetBinContent(nb);
      if(histo_WZ ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVAJERBoundingDown ->GetBinContent(nb) > 0) systJerDown[2] = histo_WZ_CMS_MVAJERBoundingDown ->GetBinContent(nb)/histo_WZ ->GetBinContent(nb);
      if(histo_ZZ ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVAJERBoundingUp   ->GetBinContent(nb) > 0) systJerUp  [3] = histo_ZZ_CMS_MVAJERBoundingUp   ->GetBinContent(nb)/histo_ZZ ->GetBinContent(nb);
      if(histo_ZZ ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVAJERBoundingDown ->GetBinContent(nb) > 0) systJerDown[3] = histo_ZZ_CMS_MVAJERBoundingDown ->GetBinContent(nb)/histo_ZZ ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb) > 0 && histo_VVV_CMS_MVAJERBoundingUp  ->GetBinContent(nb) > 0) systJerUp  [4] = histo_VVV_CMS_MVAJERBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb) > 0 && histo_VVV_CMS_MVAJERBoundingDown->GetBinContent(nb) > 0) systJerDown[4] = histo_VVV_CMS_MVAJERBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_MVAJERBoundingUp   ->GetBinContent(nb) > 0) systJerUp  [5] = histo_WS_CMS_MVAJERBoundingUp   ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_MVAJERBoundingDown ->GetBinContent(nb) > 0) systJerDown[5] = histo_WS_CMS_MVAJERBoundingDown ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WG ->GetBinContent(nb) > 0 && histo_WG_CMS_MVAJERBoundingUp   ->GetBinContent(nb) > 0) systJerUp  [6] = histo_WG_CMS_MVAJERBoundingUp   ->GetBinContent(nb)/histo_WG ->GetBinContent(nb);
      if(histo_WG ->GetBinContent(nb) > 0 && histo_WG_CMS_MVAJERBoundingDown ->GetBinContent(nb) > 0) systJerDown[6] = histo_WG_CMS_MVAJERBoundingDown ->GetBinContent(nb)/histo_WG ->GetBinContent(nb);
      if(histo_DPS->GetBinContent(nb) > 0 && histo_DPS_CMS_MVAJERBoundingUp  ->GetBinContent(nb) > 0) systJerUp  [7] = histo_DPS_CMS_MVAJERBoundingUp  ->GetBinContent(nb)/histo_DPS->GetBinContent(nb);
      if(histo_DPS->GetBinContent(nb) > 0 && histo_DPS_CMS_MVAJERBoundingDown->GetBinContent(nb) > 0) systJerDown[7] = histo_DPS_CMS_MVAJERBoundingDown->GetBinContent(nb)/histo_DPS->GetBinContent(nb);
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVAJERBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systJerUp  [8] = histo_Higgs_CMS_MVAJERBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVAJERBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systJerDown[8] = histo_Higgs_CMS_MVAJERBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);
      for(int njer=0; njer<9; njer++) if(systJerUp[njer]   > 1.10) systJerUp[njer]   = 1.10;
      for(int njer=0; njer<9; njer++) if(systJerUp[njer]   < 0.90) systJerUp[njer]   = 0.90;
      for(int njer=0; njer<9; njer++) if(systJerDown[njer] > 1.10) systJerDown[njer] = 1.10;
      for(int njer=0; njer<9; njer++) if(systJerDown[njer] < 0.90) systJerDown[njer] = 0.90;

      double systBtagUp  [9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      double systBtagDown[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      if(histo_EWK->GetBinContent(nb) > 0 && histo_EWK_CMS_MVABTAGBoundingUp  ->GetBinContent(nb) > 0) systBtagUp  [0] = histo_EWK_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_EWK->GetBinContent(nb);
      if(histo_EWK->GetBinContent(nb) > 0 && histo_EWK_CMS_MVABTAGBoundingDown->GetBinContent(nb) > 0) systBtagDown[0] = histo_EWK_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_EWK->GetBinContent(nb);
      if(histo_QCD->GetBinContent(nb) > 0 && histo_QCD_CMS_MVABTAGBoundingUp  ->GetBinContent(nb) > 0) systBtagUp  [1] = histo_QCD_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_QCD->GetBinContent(nb);
      if(histo_QCD->GetBinContent(nb) > 0 && histo_QCD_CMS_MVABTAGBoundingDown->GetBinContent(nb) > 0) systBtagDown[1] = histo_QCD_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_QCD->GetBinContent(nb);
      if(histo_WZ ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVABTAGBoundingUp   ->GetBinContent(nb) > 0) systBtagUp  [2] = histo_WZ_CMS_MVABTAGBoundingUp   ->GetBinContent(nb)/histo_WZ ->GetBinContent(nb);
      if(histo_WZ ->GetBinContent(nb) > 0 && histo_WZ_CMS_MVABTAGBoundingDown ->GetBinContent(nb) > 0) systBtagDown[2] = histo_WZ_CMS_MVABTAGBoundingDown ->GetBinContent(nb)/histo_WZ ->GetBinContent(nb);
      if(histo_ZZ ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVABTAGBoundingUp   ->GetBinContent(nb) > 0) systBtagUp  [3] = histo_ZZ_CMS_MVABTAGBoundingUp   ->GetBinContent(nb)/histo_ZZ ->GetBinContent(nb);
      if(histo_ZZ ->GetBinContent(nb) > 0 && histo_ZZ_CMS_MVABTAGBoundingDown ->GetBinContent(nb) > 0) systBtagDown[3] = histo_ZZ_CMS_MVABTAGBoundingDown ->GetBinContent(nb)/histo_ZZ ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb) > 0 && histo_VVV_CMS_MVABTAGBoundingUp  ->GetBinContent(nb) > 0) systBtagUp  [4] = histo_VVV_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb) > 0 && histo_VVV_CMS_MVABTAGBoundingDown->GetBinContent(nb) > 0) systBtagDown[4] = histo_VVV_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_MVABTAGBoundingUp   ->GetBinContent(nb) > 0) systBtagUp  [5] = histo_WS_CMS_MVABTAGBoundingUp   ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_MVABTAGBoundingDown ->GetBinContent(nb) > 0) systBtagDown[5] = histo_WS_CMS_MVABTAGBoundingDown ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WG ->GetBinContent(nb) > 0 && histo_WG_CMS_MVABTAGBoundingUp   ->GetBinContent(nb) > 0) systBtagUp  [6] = histo_WG_CMS_MVABTAGBoundingUp   ->GetBinContent(nb)/histo_WG ->GetBinContent(nb);
      if(histo_WG ->GetBinContent(nb) > 0 && histo_WG_CMS_MVABTAGBoundingDown ->GetBinContent(nb) > 0) systBtagDown[6] = histo_WG_CMS_MVABTAGBoundingDown ->GetBinContent(nb)/histo_WG ->GetBinContent(nb);
      if(histo_DPS->GetBinContent(nb) > 0 && histo_DPS_CMS_MVABTAGBoundingUp  ->GetBinContent(nb) > 0) systBtagUp  [7] = histo_DPS_CMS_MVABTAGBoundingUp  ->GetBinContent(nb)/histo_DPS->GetBinContent(nb);
      if(histo_DPS->GetBinContent(nb) > 0 && histo_DPS_CMS_MVABTAGBoundingDown->GetBinContent(nb) > 0) systBtagDown[7] = histo_DPS_CMS_MVABTAGBoundingDown->GetBinContent(nb)/histo_DPS->GetBinContent(nb);
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVABTAGBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systBtagUp  [8] = histo_Higgs_CMS_MVABTAGBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVABTAGBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systBtagDown[8] = histo_Higgs_CMS_MVABTAGBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);
      for(int nbtag=0; nbtag<9; nbtag++) if(systBtagUp[nbtag]   > 1.10) systBtagUp[nbtag]   = 1.10;
      for(int nbtag=0; nbtag<9; nbtag++) if(systBtagUp[nbtag]   < 0.90) systBtagUp[nbtag]   = 0.90;
      for(int nbtag=0; nbtag<9; nbtag++) if(systBtagDown[nbtag] > 1.10) systBtagDown[nbtag] = 1.10;
      for(int nbtag=0; nbtag<9; nbtag++) if(systBtagDown[nbtag] < 0.90) systBtagDown[nbtag] = 0.90;

      double systPUUp  [9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      double systPUDown[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
      if(histo_EWK->GetBinContent(nb) > 0 && histo_EWK_CMS_PUBoundingUp  ->GetBinContent(nb) > 0) systPUUp  [0] = histo_EWK_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_EWK->GetBinContent(nb);
      if(histo_EWK->GetBinContent(nb) > 0 && histo_EWK_CMS_PUBoundingDown->GetBinContent(nb) > 0) systPUDown[0] = histo_EWK_CMS_PUBoundingDown->GetBinContent(nb)/histo_EWK->GetBinContent(nb);
      if(histo_QCD->GetBinContent(nb) > 0 && histo_QCD_CMS_PUBoundingUp  ->GetBinContent(nb) > 0) systPUUp  [1] = histo_QCD_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_QCD->GetBinContent(nb);
      if(histo_QCD->GetBinContent(nb) > 0 && histo_QCD_CMS_PUBoundingDown->GetBinContent(nb) > 0) systPUDown[1] = histo_QCD_CMS_PUBoundingDown->GetBinContent(nb)/histo_QCD->GetBinContent(nb);
      if(histo_WZ ->GetBinContent(nb) > 0 && histo_WZ_CMS_PUBoundingUp   ->GetBinContent(nb) > 0) systPUUp  [2] = histo_WZ_CMS_PUBoundingUp   ->GetBinContent(nb)/histo_WZ ->GetBinContent(nb);
      if(histo_WZ ->GetBinContent(nb) > 0 && histo_WZ_CMS_PUBoundingDown ->GetBinContent(nb) > 0) systPUDown[2] = histo_WZ_CMS_PUBoundingDown ->GetBinContent(nb)/histo_WZ ->GetBinContent(nb);
      if(histo_ZZ ->GetBinContent(nb) > 0 && histo_ZZ_CMS_PUBoundingUp   ->GetBinContent(nb) > 0) systPUUp  [3] = histo_ZZ_CMS_PUBoundingUp   ->GetBinContent(nb)/histo_ZZ ->GetBinContent(nb);
      if(histo_ZZ ->GetBinContent(nb) > 0 && histo_ZZ_CMS_PUBoundingDown ->GetBinContent(nb) > 0) systPUDown[3] = histo_ZZ_CMS_PUBoundingDown ->GetBinContent(nb)/histo_ZZ ->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb) > 0 && histo_VVV_CMS_PUBoundingUp  ->GetBinContent(nb) > 0) systPUUp  [4] = histo_VVV_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_VVV->GetBinContent(nb) > 0 && histo_VVV_CMS_PUBoundingDown->GetBinContent(nb) > 0) systPUDown[4] = histo_VVV_CMS_PUBoundingDown->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_PUBoundingUp   ->GetBinContent(nb) > 0) systPUUp  [5] = histo_WS_CMS_PUBoundingUp   ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_PUBoundingDown ->GetBinContent(nb) > 0) systPUDown[5] = histo_WS_CMS_PUBoundingDown ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WG ->GetBinContent(nb) > 0 && histo_WG_CMS_PUBoundingUp   ->GetBinContent(nb) > 0) systPUUp  [6] = histo_WG_CMS_PUBoundingUp   ->GetBinContent(nb)/histo_WG ->GetBinContent(nb);
      if(histo_WG ->GetBinContent(nb) > 0 && histo_WG_CMS_PUBoundingDown ->GetBinContent(nb) > 0) systPUDown[6] = histo_WG_CMS_PUBoundingDown ->GetBinContent(nb)/histo_WG ->GetBinContent(nb);
      if(histo_DPS->GetBinContent(nb) > 0 && histo_DPS_CMS_PUBoundingUp  ->GetBinContent(nb) > 0) systPUUp  [7] = histo_DPS_CMS_PUBoundingUp  ->GetBinContent(nb)/histo_DPS->GetBinContent(nb);
      if(histo_DPS->GetBinContent(nb) > 0 && histo_DPS_CMS_PUBoundingDown->GetBinContent(nb) > 0) systPUDown[7] = histo_DPS_CMS_PUBoundingDown->GetBinContent(nb)/histo_DPS->GetBinContent(nb);
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_PUBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systPUUp  [8] = histo_Higgs_CMS_PUBoundingUp  [TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);
      if(histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb) > 0 && histo_Higgs_CMS_PUBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb) > 0) systPUDown[8] = histo_Higgs_CMS_PUBoundingDown[TMath::Max(nModel,0)]->GetBinContent(nb)/histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);
      for(int npu=0; npu<9; npu++) if(systPUUp[npu]   > 1.02) systPUUp[npu]   = 1.02;
      for(int npu=0; npu<9; npu++) if(systPUUp[npu]   < 0.98) systPUUp[npu]   = 0.98;
      for(int npu=0; npu<9; npu++) if(systPUDown[npu] > 1.02) systPUDown[npu] = 1.02;
      for(int npu=0; npu<9; npu++) if(systPUDown[npu] < 0.98) systPUDown[npu] = 0.98;	

      double systWSSFUp  [1] = {1.0};
      double systWSSFDown[1] = {1.0};
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_WSSFUp         ->GetBinContent(nb) > 0) systWSSFUp  [0] = histo_WS_CMS_WSSFUp	      ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);
      if(histo_WS ->GetBinContent(nb) > 0 && histo_WS_CMS_WSSFDown       ->GetBinContent(nb) > 0) systWSSFDown[0] = histo_WS_CMS_WSSFDown     ->GetBinContent(nb)/histo_WS ->GetBinContent(nb);

      double lumiEWZ = lumiE;
      if(useWZFromData){
        lumiEWZ = 1.0;
	systLepEffM[2] = 1.0; systLepEffE[2] = 1.0; 
	systMetUp[2] = 1.0; systJesUp[2] = 1.0; systJerUp[2] = 1.0; systBtagUp[2] = 1.0; systPUUp[2] = 1.0;
	systMetDown[2] = 1.0; systJesDown[2] = 1.0; systJerDown[2] = 1.0; systBtagDown[2] = 1.0; systPUDown[2] = 1.0;
	systLepResM[2] = 1.0; systLepResE[2] = 1.0; 
	systQCDScale[2] = 1.0; systPDF[2] = 1.0;
      }

      char outputLimits[200];
      sprintf(outputLimits,"wwss%d_%s_%s_input_%s_bin%d.root",finalVar,finalStateName,theSignalName.Data(),ECMsb.Data(),nb-1);
      TFile* outFileLimits = new TFile(outputLimits,"recreate");
      outFileLimits->cd();
    
      histoOneBin_Data ->Reset();
      histoOneBin_EWK  ->Reset();
      histoOneBin_QCD  ->Reset();
      histoOneBin_WZ   ->Reset();
      histoOneBin_ZZ   ->Reset();
      histoOneBin_VVV  ->Reset();
      histoOneBin_WS   ->Reset();
      histoOneBin_WG   ->Reset();
      histoOneBin_DPS  ->Reset();
      histoOneBin_FakeM->Reset();
      histoOneBin_FakeE->Reset();
      histoOneBin_Higgs[TMath::Max(nModel,0)]->Reset();
      histoOneBin_Data ->SetBinContent(1,            histo_Data ->GetBinContent(nb));
      histoOneBin_EWK  ->SetBinContent(1, TMath::Max(histo_EWK  ->GetBinContent(nb),0.0));
      histoOneBin_QCD  ->SetBinContent(1, TMath::Max(histo_QCD  ->GetBinContent(nb),0.0));
      histoOneBin_WZ   ->SetBinContent(1, TMath::Max(histo_WZ   ->GetBinContent(nb),0.0));
      histoOneBin_ZZ   ->SetBinContent(1, TMath::Max(histo_ZZ   ->GetBinContent(nb),0.0));
      histoOneBin_VVV  ->SetBinContent(1, TMath::Max(histo_VVV  ->GetBinContent(nb),0.0));
      histoOneBin_WS   ->SetBinContent(1, TMath::Max(histo_WS   ->GetBinContent(nb),0.0));
      histoOneBin_WG   ->SetBinContent(1, TMath::Max(histo_WG   ->GetBinContent(nb),0.0));
      histoOneBin_DPS  ->SetBinContent(1, TMath::Max(histo_DPS  ->GetBinContent(nb),0.0));
      histoOneBin_FakeM->SetBinContent(1, TMath::Max(histo_FakeM->GetBinContent(nb),0.0));
      histoOneBin_FakeE->SetBinContent(1, TMath::Max(histo_FakeE->GetBinContent(nb),0.0));
      histoOneBin_Higgs[TMath::Max(nModel,0)]  ->SetBinContent(1, TMath::Max(histo_Higgs[TMath::Max(nModel,0)]  ->GetBinContent(nb),0.0));
      histoOneBin_Data ->Write();
      histoOneBin_EWK  ->Write();
      histoOneBin_QCD  ->Write();
      histoOneBin_WZ   ->Write();
      histoOneBin_ZZ   ->Write();
      histoOneBin_VVV  ->Write();
      histoOneBin_WS   ->Write();
      histoOneBin_WG   ->Write();
      histoOneBin_DPS  ->Write();
      histoOneBin_FakeM->Write();
      histoOneBin_FakeE->Write();
      histoOneBin_Higgs[TMath::Max(nModel,0)]  ->Write();

      outFileLimits->Close();

      double sigYield = 0;
      if(nModel >= 0) sigYield = histo_Higgs[TMath::Max(nModel,0)]->GetBinContent(nb);

      char outputLimitsShape[200];                                            
      sprintf(outputLimitsShape,"histo_limits_wwss%d_%s_%s_shape_%s_bin%d.txt",finalVar,finalStateName,theSignalName.Data(),ECMsb.Data(),nb-1);
      ofstream newcardShape;
      newcardShape.open(outputLimitsShape);
      newcardShape << Form("imax 1 number of channels\n");
      newcardShape << Form("jmax * number of background\n");
      newcardShape << Form("kmax * number of nuisance parameters\n");
      newcardShape << Form("shapes *   *   %s  histoOneBin_$PROCESS\n",outputLimits);
      newcardShape << Form("shapes data_obs * %s  histoOneBin_Data \n",outputLimits);
      newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
      newcardShape << Form("bin %2s%4s%d %2s%4s%d %2s%4s%d %2s%4s%d %2s%4s%d %2s%4s%d %2s%4s%d %2s%4s%d %2s%4s%d %2s%4s%d %2s%4s%d\n",finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1,finalStateName,ECMsb.Data(),nb-1);
      newcardShape << Form("process EWK QCD WZ ZZ VVV WS WG DPS FakeM FakeE Higgs\n");
      if(nModel == -1) newcardShape << Form("process  0 1 2 3 4 5 6 7 8 9 10\n");
      else             newcardShape << Form("process 10 1 2 3 4 5 6 7 8 9  0\n");
      newcardShape << Form("rate %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n",TMath::Max(histo_EWK->GetBinContent(nb),0.0),TMath::Max(histo_QCD->GetBinContent(nb),0.0),TMath::Max(histo_WZ->GetBinContent(nb),0.0),TMath::Max(histo_ZZ->GetBinContent(nb),0.0),TMath::Max(histo_VVV->GetBinContent(nb),0.0),TMath::Max(histo_WS->GetBinContent(nb),0.0),TMath::Max(histo_WG->GetBinContent(nb),0.0),TMath::Max(histo_DPS->GetBinContent(nb),0.0),TMath::Max(histo_FakeM->GetBinContent(nb),0.0),TMath::Max(histo_FakeE->GetBinContent(nb),0.0),sigYield);
      newcardShape << Form("lumi_%4s       lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f   -	-   %7.5f\n",ECMsb.Data(),lumiE,lumiE,lumiEWZ,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE);		     
      newcardShape << Form("%s             lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f   -	-   %7.5f\n",effMName,systLepEffM[0],systLepEffM[1],systLepEffM[2],systLepEffM[3],systLepEffM[4],systLepEffM[5],systLepEffM[6],systLepEffM[7],systLepEffM[8]);
      newcardShape << Form("%s             lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f   -	-   %7.5f\n",effEName,systLepEffE[0],systLepEffE[1],systLepEffE[2],systLepEffE[3],systLepEffE[4],systLepEffE[5],systLepEffE[6],systLepEffE[7],systLepEffE[8]);
      newcardShape << Form("%s             lnN  %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -	-   %7.5f/%7.5f\n",jesName,systJesUp[0],systJesDown[0],systJesUp[1],systJesDown[1],systJesUp[2],systJesDown[2],systJesUp[3],systJesDown[3],systJesUp[4],systJesDown[4],systJesUp[5],systJesDown[5],systJesUp[6],systJesDown[6],systJesUp[7],systJesDown[7],systJesUp[8],systJesDown[8]);
      newcardShape << Form("%s             lnN  %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -	-   %7.5f/%7.5f\n",jerName,systJerUp[0],systJerDown[0],systJerUp[1],systJerDown[1],systJerUp[2],systJerDown[2],systJerUp[3],systJerDown[3],systJerUp[4],systJerDown[4],systJerUp[5],systJerDown[5],systJerUp[6],systJerDown[6],systJerUp[7],systJerDown[7],systJerUp[8],systJerDown[8]);
      newcardShape << Form("%s             lnN  %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -	-   %7.5f/%7.5f\n",metName,systMetUp[0],systMetDown[0],systMetUp[1],systMetDown[1],systMetUp[2],systMetDown[2],systMetUp[3],systMetDown[3],systMetUp[4],systMetDown[4],systMetUp[5],systMetDown[5],systMetUp[6],systMetDown[6],systMetUp[7],systMetDown[7],systMetUp[8],systMetDown[8]);
      newcardShape << Form("%s             lnN  %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -	-   %7.5f/%7.5f\n",btagName,systBtagUp[0],systBtagDown[0],systBtagUp[1],systBtagDown[1],systBtagUp[2],systBtagDown[2],systBtagUp[3],systBtagDown[3],systBtagUp[4],systBtagDown[4],systBtagUp[5],systBtagDown[5],systBtagUp[6],systBtagDown[6],systBtagUp[7],systBtagDown[7],systBtagUp[8],systBtagDown[8]);
      newcardShape << Form("%s             lnN  %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f %7.5f/%7.5f   -	-   %7.5f/%7.5f\n",puName,systPUUp[0],systPUDown[0],systPUUp[1],systPUDown[1],systPUUp[2],systPUDown[2],systPUUp[3],systPUDown[3],systPUUp[4],systPUDown[4],systPUUp[5],systPUDown[5],systPUUp[6],systPUDown[6],systPUUp[7],systPUDown[7],systPUUp[8],systPUDown[8]);
      newcardShape << Form("%s             lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f   -	-   %7.5f\n",momMName,systLepResM[0],systLepResM[1],systLepResM[2],systLepResM[3],systLepResM[4],systLepResM[5],systLepResM[6],systLepResM[7],systLepResM[8]);
      newcardShape << Form("%s             lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f   -	-   %7.5f\n",momEName,systLepResE[0],systLepResE[1],systLepResE[2],systLepResE[3],systLepResE[4],systLepResE[5],systLepResE[6],systLepResE[7],systLepResE[8]);
      newcardShape << Form("CMS_WZ_l       lnN    -     -   %7.5f   -     -     -    -     -      -     -     -  \n",1.01);
      newcardShape << Form("QCDscale_EWK   lnN  %7.5f   -     -     -     -     -     -     -     -     -     -  \n",systQCDScale[0]);		
      newcardShape << Form("QCDscale_VV    lnN    -   %7.5f %7.5f %7.5f   -     -     -   %7.5f   -     -     -  \n",systQCDScale[1],systQCDScale[2],systQCDScale[3],systQCDScale[7]);		
      newcardShape << Form("QCDscale_VVV   lnN    -     -     -     -   %7.5f   -     -     -     -     -     -  \n",systQCDScale[4]);		
      newcardShape << Form("QCDscale_tt    lnN    -     -     -     -     -   %7.5f   -     -     -     -     -  \n",systQCDScale[5]);		
      newcardShape << Form("QCDscale_WG    lnN    -     -     -     -     -     -   %7.5f   -     -     -     -  \n",systQCDScale[6]);		
      newcardShape << Form("QCDscale_Higgs lnN    -     -     -     -     -     -     -     -     -     -   %7.5f\n",systQCDScale[8]);		
      newcardShape << Form("pdf_qqbar      lnN  %7.5f %7.5f %7.5f %7.5f %7.5f   -   %7.5f %7.5f   -     -   %7.5f\n",systPDF[0],systPDF[1],systPDF[2],systPDF[3],systPDF[4],systPDF[6],systPDF[7],systPDF[8]);
      newcardShape << Form("pdf_gg         lnN    -     -     -     -     -   %7.5f   -     -     -     -     -  \n",systPDF[5]);
      if(theControlRegion != 2){

      //if(histo_FakeM->GetBinContent(nb)>0 && histo_Fake_CMS_SystM->GetBinContent(nb)>0)
      //newcardShape << Form("CMS_FakeM      lnN    -     -     -     -     -     -    -     -    %7.5f   -     -  \n",1.0+0.3*TMath::Min(histo_Fake_CMS_SystM->GetBinContent(nb)/histo_FakeM->GetBinContent(nb) ,1.0));
      //if(histo_FakeM->GetBinContent(nb)>0 && histo_Fake_CMS_SystE->GetBinContent(nb)>0)
      //newcardShape << Form("CMS_FakeE      lnN    -     -     -     -     -     -    -     -    %7.5f   -     -  \n",1.0+0.3*TMath::Min(histo_Fake_CMS_SystE->GetBinContent(nb)/histo_FakeM->GetBinContent(nb) ,1.0));

      newcardShape << Form("CMS_FakeM      lnN    -     -     -     -     -     -    -     -    %7.5f   -     -  \n",1.30);

      if(histo_FakeE->GetBinContent(nb)>0)
      newcardShape << Form("CMS_FakeE      lnN    -     -     -     -     -     -    -     -      -   %7.5f   -  \n",1.30);

      } else {
      if(histo_FakeM->GetBinContent(nb)>0)
      newcardShape << Form("CMS_Fake3l     lnN    -     -     -     -     -     -    -     -    %7.5f   -     -  \n",1.30);
      }

      newcardShape << Form("WS_Norm        lnN    -     -     -     -     -   %7.5f  -     -      -     -     -  \n",systWSSFUp[0]);		
      if(useWZFromData){
      newcardShape << Form("CMS_WZ_Norm%d  lnN    -     -   %7.5f   -     -     -    -     -      -     -     -  \n",(nb-1)%4,1.0+sfE_WZ[(nb-1)%4]);
      }
      else {
      newcardShape << Form("CMS_wwss_WZnorm_bin%d rateParam  * WZ 1 [0.1,10]\n",(nb-1)%nBinWZMVAModule);         
      }

      //newcardShape << Form("CMS_wwss_wjetsMnorm_bin rateParam  * FakeM 1 [0.1,10]\n");         
      //newcardShape << Form("CMS_wwss_wjetsEnorm_bin rateParam  * FakeE 1 [0.1,10]\n");         

      if(histo_EWK  ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_EWKStatBounding2016_%s_Bin%d    lnN %7.5f   -     -     -     -     -     -     -     -     -   - \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_EWK  ->GetBinError(nb)/histo_EWK  ->GetBinContent(nb),0.999));
      if(histo_QCD  ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_QCDStatBounding2016_%s_Bin%d    lnN   -   %7.5f   -     -     -     -     -     -     -     -   - \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_QCD  ->GetBinError(nb)/histo_QCD  ->GetBinContent(nb),0.999));
      if(histo_WZ   ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_WZStatBounding2016_%s_Bin%d     lnN   -     -   %7.5f   -     -     -     -     -     -     -   - \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_WZ   ->GetBinError(nb)/histo_WZ   ->GetBinContent(nb),0.999));
      if(histo_ZZ   ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_ZZStatBounding2016_%s_Bin%d     lnN   -     -     -   %7.5f   -     -     -     -     -     -   - \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_ZZ   ->GetBinError(nb)/histo_ZZ   ->GetBinContent(nb),0.999));
      if(histo_VVV  ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_VVVStatBounding2016_%s_Bin%d    lnN   -     -     -     -   %7.5f   -     -     -     -     -   - \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_VVV  ->GetBinError(nb)/histo_VVV  ->GetBinContent(nb),0.999));
      if(histo_WS   ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_WSStatBounding2016_%s_Bin%d     lnN   -     -     -     -     -   %7.5f   -     -     -     -   - \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_WS   ->GetBinError(nb)/histo_WS   ->GetBinContent(nb),0.999));
      if(histo_WG   ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_WGStatBounding2016_%s_Bin%d     lnN   -     -     -     -     -     -   %7.5f   -     -     -   - \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_WG   ->GetBinError(nb)/histo_WG   ->GetBinContent(nb),0.999));
      if(histo_DPS  ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_DPSStatBounding2016_%s_Bin%d    lnN   -     -     -     -     -     -     -   %7.5f   -     -   - \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_DPS  ->GetBinError(nb)/histo_DPS  ->GetBinContent(nb),0.999));
      if(histo_FakeM->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_FakeMStatBounding2016_%s_Bin%d  lnN   -     -     -     -     -     -     -     -   %7.5f   -   - \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_FakeM->GetBinError(nb)/histo_FakeM->GetBinContent(nb),0.999));
      if(histo_FakeE->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_FakeEStatBounding2016_%s_Bin%d  lnN   -     -     -     -     -     -     -     -    -   %7.5f  - \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_FakeE->GetBinError(nb)/histo_FakeE->GetBinContent(nb),0.999));
      if(histo_Higgs[TMath::Max(nModel,0)]  ->GetBinContent(nb) > 0) newcardShape << Form("CMS_wwss%s_histo_HiggsStatBounding2016_%s_Bin%d    lnN -   -     -     -     -     -     -     -     -     -    %7.5f \n",finalStateName,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_Higgs[TMath::Max(nModel,0)]  ->GetBinError(nb)/histo_Higgs[TMath::Max(nModel,0)]  ->GetBinContent(nb),0.999));
    }
  }
  for(int iF=0; iF<3; iF++){
    printf("double jetEpsBtagLOOSE[%d] = \n",iF);
    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
        printf("%5.4f",numBTaggingLOOSE[iPt][iEta][iF]/denBTagging[iPt][iEta][iF]);
        if(iPt!=4||iEta!=4) printf(",");
        if(iPt==4) printf("\n");
      }
    }
  }
  for(int iF=0; iF<3; iF++){
    printf("double jetEpsBtagTIGHT[%d] = \n",iF);
    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
        printf("%5.4f",numBTaggingTIGHT[iPt][iEta][iF]/denBTagging[iPt][iEta][iF]);
        if(iPt!=4||iEta!=4) printf(",");
        if(iPt==4) printf("\n");
      }
    }
  }
  delete btagCalib;
}

void func_ws_sf(double eta, double pt, double SF[2]){
  int iEta = -1;
  if	 (TMath::Abs(eta) < 0.5) iEta = 0;
  else if(TMath::Abs(eta) < 1.0) iEta = 1;
  else if(TMath::Abs(eta) < 1.5) iEta = 2;
  else if(TMath::Abs(eta) < 2.0) iEta = 3;
  else                           iEta = 4;
  SF[0] = SF[0] * WSSF[iEta];
  SF[1] = SF[1] * (1.0 + WSSFE[iEta]);
}

double func_ws_eff(double eta1, double eta2, TH1D *fhEff){
  return TMath::Max(fhEff->GetBinContent(fhEff->GetXaxis()->FindBin(eta1)),
                    fhEff->GetBinContent(fhEff->GetXaxis()->FindBin(eta2)));
}
