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

#include "WWAnalysis/resummation/WWpTreweight.h"

enum selType                     { SIGSEL,   SSSEL,   TOPSEL,   DYSEL, nSelTypes};
TString selTypeName[nSelTypes]= { "SIGSEL", "SSSEL", "TOPSEL", "DYSEL"};

enum systType                     {JESUP=0, JESDOWN,  METUP,  METDOWN,   JERUP,   JERDOWN, nSystTypes};
TString systTypeName[nSystTypes]= {"JESUP","JESDOWN","METUP","METDOWN", "JERUP", "JERDOWN"};

const double bTagCuts[1] = {0.5426}; // 0.5426/0.8484/0.9535 (check BTagCalibration2Reader!)

bool isMINIAOD = true;
int whichSkim = 3;
double mcPrescale = 1.0;
bool usePureMC = false;
bool applyGStarVeto = true; 
int period = 1;
const bool useDYMVA = false;
const bool usePUPPI = false;
const bool useTopCR = true;
//const TString typeLepSel = "default";

double topNorm[3]  = {0.81,0.95,1.00};
double topNormE[3] = {0.09,0.11,0.01};

void wwAnalysis(
 unsigned int nJetsType = 0,
 int theControlRegion = 0,
 TString typeLepSel = "defaultTight",
 bool isShapeAna = false,
 bool isMIT = true
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

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  TString triggerSuffix = "*";
  if(isMINIAOD) triggerSuffix = "";
  if      (period==1){
  puPath = "MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root";

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

  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg.root",filesPathMC.Data()));                                            infilecatv.push_back(1);

  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV.root",filesPathMC.Data()));					      infilecatv.push_back(2);

  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg.root",filesPathMC2.Data()));					      infilecatv.push_back(3);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));    infilecatv.push_back(3);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",filesPathMC.Data()));infilecatv.push_back(3);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));                   infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));                       infilecatv.push_back(4);
  infilenamev.push_back(Form("%sDYJetsToTauTau_ForcedMuEleDecay_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data())); infilecatv.push_back(4);

  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8.root",filesPathMC.Data()));  				      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8.root",filesPathMC.Data()));					      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));	 	      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8.root",filesPathMC.Data()));                        infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWZTo3LNu_mllmin01_13TeV-powheg-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(5);

  infilenamev.push_back(Form("%sWWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data())); 			      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data())); 			      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                           infilecatv.push_back(6);
  infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                           infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));	      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(6);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root",filesPathMC.Data()));           infilecatv.push_back(6);
  infilenamev.push_back(Form("%stZq_ll_4f_13TeV-amcatnlo-pythia8.root",filesPathMC.Data()));                                  infilecatv.push_back(6);

  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));                  infilecatv.push_back(7);
  infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",filesPathMC.Data()));                   infilecatv.push_back(7);
  
  //if(applyGStarVeto == true) {
  //infilenamev.push_back(Form("%sWGstarToLNuEE_012Jets_13TeV-madgraph.root",filesPathMC.Data()));                              infilecatv.push_back(8);
  //infilenamev.push_back(Form("%sWGstarToLNuMuMu_012Jets_13TeV-madgraph.root",filesPathMC.Data()));                            infilecatv.push_back(8);
  //}

  //infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",filesPathMC.Data()));		   infilecatv.push_back(9);

  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root",filesPathMC.Data()));		      infilecatv.push_back(11);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgenv628_pythia8.root",filesPathMC.Data()));              infilecatv.push_back(11);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));			      infilecatv.push_back(11);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8.root",filesPathMC.Data()));                            infilecatv.push_back(11);
  }
  else {assert(0); return;}

  if(infilenamev.size() != infilecatv.size()) {assert(0); return;}
  
  //infilenamev.clear();infilecatv.clear();
  //infilenamev.push_back(Form("/mnt/hadoop/scratch/ceballos/Nero/v2.2/output_80x/mc/GluGluWWTo2L2Nu_MCFM_13TeV.root")); infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg-herwigpp.root",filesPathMC.Data())); infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8.root",filesPathMC.Data()));infilecatv.push_back(11);
  //infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_herwigpp.root",filesPathMC.Data()));infilecatv.push_back(11);
  //infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg-CUETP8M1Up.root",filesPathMC.Data())); infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg-CUETP8M1Down.root",filesPathMC.Data())); infilecatv.push_back(1);

  double denBTagging[5][5][3],jetEpsBtagLOOSE[5][5][3];
  double numBTaggingLOOSE[5][5][3];
  for(int i0=0; i0<5; i0++) {
    for(int i1=0; i1<5; i1++) {
      for(int i2=0; i2<3; i2++) {
        denBTagging[i0][i1][i2] = 0.0;
        numBTaggingLOOSE[i0][i1][i2] = 0.0;
	if     (i2==BTagEntry::FLAV_B)    jetEpsBtagLOOSE[i0][i1][i2] = jetEpsBtagBLOOSE[i0][i1];
	else if(i2==BTagEntry::FLAV_C)    jetEpsBtagLOOSE[i0][i1][i2] = jetEpsBtagCLOOSE[i0][i1];
	else if(i2==BTagEntry::FLAV_UDSG) jetEpsBtagLOOSE[i0][i1][i2] = jetEpsBtagLLOOSE[i0][i1];
      }
    }
  }

  //Float_t fMVACut[4][4];
  //InitializeJetIdCuts(fMVACut);
  
  BTagCalibration2 *btagCalib = new BTagCalibration2("csvv2","MitAnalysisRunII/data/80x/CSVv2_Moriond17_B_H.csv");
  BTagCalibration2Reader btagReaderBCLOOSE(btagCalib,BTagEntry::OP_LOOSE,"comb","central");
  BTagCalibration2Reader btagReaderLLOOSE(btagCalib,BTagEntry::OP_LOOSE,"incl","central");
  BTagCalibration2Reader btagReaderBCLOOSEUP(btagCalib,BTagEntry::OP_LOOSE,"comb","up");
  BTagCalibration2Reader btagReaderLLOOSEUP(btagCalib,BTagEntry::OP_LOOSE,"incl","up");
  BTagCalibration2Reader btagReaderBCLOOSEDOWN(btagCalib,BTagEntry::OP_LOOSE,"comb","down");
  BTagCalibration2Reader btagReaderLLOOSEDOWN(btagCalib,BTagEntry::OP_LOOSE,"incl","down");

  WWpTreweight theWWpTreweight(
              "WWAnalysis/resummation/central.dat",
              "WWAnalysis/resummation/resum_up.dat",
              "WWAnalysis/resummation/resum_down.dat",
              "WWAnalysis/resummation/scale_up.dat",
              "WWAnalysis/resummation/scale_down.dat",
              "WWAnalysis/resummation/nnlo_central.dat",
              "WWAnalysis/resummation/powheg_2l2nu_nlo.dat",
              "WWAnalysis/resummation/powheg_2l2nu_qup_nlo.dat",
              "WWAnalysis/resummation/powheg_2l2nu_qdown_nlo.dat",
              "WWAnalysis/resummation/powheg_2l2nu_sup_nlo.dat",
              "WWAnalysis/resummation/powheg_2l2nu_sdown_nlo.dat",
              "WWAnalysis/resummation/powheg_2l2nu_nnlo.dat");

  float dymva_= -999.;
  unsigned int nlep_= -1;
  unsigned int njets_= -1;

  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU     = (TH1D*)(fPUFile->Get("puWeights"));     assert(fhDPU);    fhDPU    ->SetDirectory(0);
  TH1D *fhDPUUp   = (TH1D*)(fPUFile->Get("puWeightsUp"));   assert(fhDPUUp);  fhDPUUp  ->SetDirectory(0);
  TH1D *fhDPUDown = (TH1D*)(fPUFile->Get("puWeightsDown")); assert(fhDPUDown);fhDPUDown->SetDirectory(0);
  delete fPUFile;

  TFile *fTrackElectronReco_SF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root"));
  TH2D *fhDeltrksf= (TH2D*)(fTrackElectronReco_SF->Get("scalefactors_Reco_Electron")); assert(fhDeltrksf); fhDeltrksf->SetDirectory(0);
  delete fTrackElectronReco_SF;

  TFile *fElSF_latinos = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_latinos_37ifb.root"));
  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("scalefactors_Medium_Electron"));
  TH2D *fhDElTightSF;
  if(strcmp(typeLepSel.Data(),"defaultTight")==0){
    printf("Using defaultTight SF\n");
    fhDElTightSF = (TH2D*)(fElSF_latinos->Get("scalefactors_Tight_Electron"));
  } else {
    fhDElTightSF = (TH2D*)(fElSF->Get("scalefactors_Tight_Electron"));
  }
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF->SetDirectory(0);
  delete fElSF;
  delete fElSF_latinos;

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

  TFile *fWWPtRatio = TFile::Open(Form("MitAnalysisRunII/data/74x/MyRatioWWpTHistogramAll.root"));
  TH1D *fhDWWPtRatio           = (TH1D*)(fWWPtRatio->Get("wwpt"));
  TH1D *fhDWWPtRatio_scaleup   = (TH1D*)(fWWPtRatio->Get("wwpt_scaleup"));
  TH1D *fhDWWPtRatio_scaledown = (TH1D*)(fWWPtRatio->Get("wwpt_scaledown"));
  TH1D *fhDWWPtRatio_resumup   = (TH1D*)(fWWPtRatio->Get("wwpt_resumup"));
  TH1D *fhDWWPtRatio_resumdown = (TH1D*)(fWWPtRatio->Get("wwpt_resumdown"));
  assert(fhDWWPtRatio	       );
  assert(fhDWWPtRatio_scaleup  );
  assert(fhDWWPtRatio_scaledown);
  assert(fhDWWPtRatio_resumup  );
  assert(fhDWWPtRatio_resumdown);
  fhDWWPtRatio  	->SetDirectory(0);
  fhDWWPtRatio_scaleup  ->SetDirectory(0);
  fhDWWPtRatio_scaledown->SetDirectory(0);
  fhDWWPtRatio_resumup  ->SetDirectory(0);
  fhDWWPtRatio_resumdown->SetDirectory(0);
  delete fWWPtRatio;

  TString ECMsb  = "13TeV2016";
  const int nBinMVA = 8; Float_t xbins[nBinMVA+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  if(isShapeAna && theControlRegion == 0){
    xbins[0] =  12; xbins[1] = 100; xbins[2] = 150; xbins[3] = 200;
    xbins[4] = 250, xbins[5] = 300; xbins[6] = 400; xbins[7] = 600;
    xbins[8] = 800;
  }
  TH1D* histoMVA = new TH1D("histoMVA", "histoMVA", nBinMVA, xbins);
  histoMVA->Sumw2();
  TH1D* histoOneBin = new TH1D("histoOneBin", "histoOneBin", 1, -0.5, 0.5);
  histoOneBin->Sumw2();

  TH1D *histo_FakeData  = (TH1D*) histoMVA->Clone("histo_FakeData");

  TH1D *histo_Data   = (TH1D*) histoMVA->Clone("histo_Data");
  TH1D *histo_qqWW   = (TH1D*) histoMVA->Clone("histo_qqWW"); 
  TH1D *histo_ggWW   = (TH1D*) histoMVA->Clone("histo_ggWW");
  TH1D *histo_Top    = (TH1D*) histoMVA->Clone("histo_Top");	 
  TH1D *histo_DY     = (TH1D*) histoMVA->Clone("histo_DY");
  TH1D *histo_VV     = (TH1D*) histoMVA->Clone("histo_VV");
  TH1D *histo_VVV    = (TH1D*) histoMVA->Clone("histo_VVV");
  TH1D *histo_WG     = (TH1D*) histoMVA->Clone("histo_WG");
  TH1D *histo_WGS    = (TH1D*) histoMVA->Clone("histo_WGS");
  TH1D *histo_WjetsM = (TH1D*) histoMVA->Clone("histo_WjetsM");	 
  TH1D *histo_WjetsE = (TH1D*) histoMVA->Clone("histo_WjetsE");	 
  TH1D *histo_Higgs  = (TH1D*) histoMVA->Clone("histo_Higgs");	 

  TH1D *histoOneBin_Data   = (TH1D*) histoOneBin->Clone("histoOneBin_Data");
  TH1D *histoOneBin_qqWW   = (TH1D*) histoOneBin->Clone("histoOneBin_qqWW"); 
  TH1D *histoOneBin_ggWW   = (TH1D*) histoOneBin->Clone("histoOneBin_ggWW");
  TH1D *histoOneBin_Top    = (TH1D*) histoOneBin->Clone("histoOneBin_Top");   
  TH1D *histoOneBin_DY     = (TH1D*) histoOneBin->Clone("histoOneBin_DY");
  TH1D *histoOneBin_VV     = (TH1D*) histoOneBin->Clone("histoOneBin_VV");
  TH1D *histoOneBin_VVV    = (TH1D*) histoOneBin->Clone("histoOneBin_VVV");
  TH1D *histoOneBin_WG     = (TH1D*) histoOneBin->Clone("histoOneBin_WG");
  TH1D *histoOneBin_WGS    = (TH1D*) histoOneBin->Clone("histoOneBin_WGS");
  TH1D *histoOneBin_WjetsM = (TH1D*) histoOneBin->Clone("histoOneBin_WjetsM");        
  TH1D *histoOneBin_WjetsE = (TH1D*) histoOneBin->Clone("histoOneBin_WjetsE");        
  TH1D *histoOneBin_Higgs  = (TH1D*) histoOneBin->Clone("histoOneBin_Higgs");	      

  double totalFakeDataCount[4][5];
  for(int i=0; i<4; i++) for(int j=0; j<5; j++) totalFakeDataCount[i][j] = 0;
  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 25;
  const int histBins = 12;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {".Data", ".qqWW", ".ggWW", "..Top", "...DY", "...VV", "..VVV", "...WG", "..WGS", "WjetsM", "WjetsE", "Higgs"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    bool isMVAPlot = false;
    if     (thePlot >=  0 && thePlot <=  2) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  3 && thePlot <=  3) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >=  4 && thePlot <=  4) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >=  5 && thePlot <=  5) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;if(isShapeAna) isMVAPlot = true;}

    else if(thePlot >=  6 && thePlot <=  8) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  9 && thePlot <=  9) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 10 && thePlot <= 10) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >= 11 && thePlot <= 11) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;if(isShapeAna) isMVAPlot = true;}

    else if(thePlot >= 12 && thePlot <= 14) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 15 && thePlot <= 15) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 16 && thePlot <= 16) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >= 17 && thePlot <= 17) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;if(isShapeAna) isMVAPlot = true;}

    else if(thePlot >= 18 && thePlot <= 20) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >= 21 && thePlot <= 21) {nBinPlot =  90; xminPlot = 0.0; xmaxPlot = 180.0;}
    else if(thePlot >= 22 && thePlot <= 22) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >= 23 && thePlot <= 23) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 800.0;if(isShapeAna) isMVAPlot = true;}

    else if(thePlot >= 24 && thePlot <= 24) {nBinPlot = 500; xminPlot = 0.0; xmaxPlot = 500.0;}

    TH1D* histos;
    if(isMVAPlot == false) histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    else                   histos = new TH1D("histos", "histos", nBinMVA, xbins);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Reset();histos->Clear();
  }

  char finalStateName[2],effMName[10],effEName[10],momMName[10],momEName[10];
  sprintf(effMName,"CMS_eff_m");sprintf(momMName,"CMS_scale_m");
  sprintf(effEName,"CMS_eff_e");sprintf(momEName,"CMS_scale_e");
  sprintf(finalStateName,"wwlnln");
  if (theControlRegion == 1) sprintf(finalStateName,"wwtop");

  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingUp       = new TH1D( Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingDown     = new TH1D( Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingUp       = new TH1D( Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingDown     = new TH1D( Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVATopStatBoundingUp         = new TH1D( Form("histo_Top_CMS_ww%s_MVATopStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Top_CMS_ww%s_MVATopStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVATopStatBoundingDown       = new TH1D( Form("histo_Top_CMS_ww%s_MVATopStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Top_CMS_ww%s_MVATopStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVADYStatBoundingUp           = new TH1D( Form("histo_DY_CMS_ww%s_MVADYStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_DY_CMS_ww%s_MVADYStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_DY_CMS_MVADYStatBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVADYStatBoundingDown         = new TH1D( Form("histo_DY_CMS_ww%s_MVADYStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_DY_CMS_ww%s_MVADYStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_DY_CMS_MVADYStatBoundingDown->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingUp     = new TH1D( Form("histo_Zjets_CMS_ww%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_ww%s_MVAZjetsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingUp  ->Sumw2();
  TH1D* histo_Zjets_CMS_MVAZjetsStatBoundingDown   = new TH1D( Form("histo_Zjets_CMS_ww%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Zjets_CMS_ww%s_MVAZjetsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Zjets_CMS_MVAZjetsStatBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVAVVStatBoundingUp           = new TH1D( Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVAVVStatBoundingDown         = new TH1D( Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingUp         = new TH1D( Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingDown       = new TH1D( Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAWGStatBoundingUp           = new TH1D( Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WG_CMS_MVAWGStatBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAWGStatBoundingDown         = new TH1D( Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WG_CMS_MVAWGStatBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVAWGSStatBoundingUp         = new TH1D( Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WGS_CMS_MVAWGSStatBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVAWGSStatBoundingDown       = new TH1D( Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WGS_CMS_MVAWGSStatBoundingDown->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingUp   = new TH1D( Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingUp  ->Sumw2();
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingDown = new TH1D( Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingDown->Sumw2();
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingUp   = new TH1D( Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingUp  ->Sumw2();
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingDown = new TH1D( Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingUp     = new TH1D( Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%sUp"  ,finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingDown   = new TH1D( Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%sDown",finalStateName,ECMsb.Data()), Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%sDown",finalStateName,ECMsb.Data()), nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVAWW_nlo       = new TH1D( Form("histo_qqWW_CMS_MVAWW_nlo"), Form("histo_qqWW_CMS_MVAWW_nlo"), nBinMVA, xbins); histo_qqWW_CMS_MVAWW_nlo->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWW_qup_nlo   = new TH1D( Form("histo_qqWW_CMS_MVAWW_qup_nlo"), Form("histo_qqWW_CMS_MVAWW_qup_nlo"), nBinMVA, xbins); histo_qqWW_CMS_MVAWW_qup_nlo->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWW_qdown_nlo = new TH1D( Form("histo_qqWW_CMS_MVAWW_qdown_nlo"), Form("histo_qqWW_CMS_MVAWW_qdown_nlo"), nBinMVA, xbins); histo_qqWW_CMS_MVAWW_qdown_nlo->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWW_sup_nlo   = new TH1D( Form("histo_qqWW_CMS_MVAWW_sup_nlo"), Form("histo_qqWW_CMS_MVAWW_sup_nlo"), nBinMVA, xbins); histo_qqWW_CMS_MVAWW_sup_nlo->Sumw2();
  TH1D* histo_qqWW_CMS_MVAWW_sdown_nlo = new TH1D( Form("histo_qqWW_CMS_MVAWW_sdown_nlo"), Form("histo_qqWW_CMS_MVAWW_sdown_nlo"), nBinMVA, xbins); histo_qqWW_CMS_MVAWW_sdown_nlo->Sumw2();

  TH1D* histo_Diff = new TH1D("dummy", "dummy",1000,-1,1); histo_Diff->Sumw2();

  TH1D* histo_qqWW_CMS_QCDScaleBounding[6];
  TH1D* histo_ggWW_CMS_QCDScaleBounding[6];
  TH1D* histo_Top_CMS_QCDScaleBounding[6];
  TH1D* histo_DY_CMS_QCDScaleBounding[6];
  TH1D* histo_VV_CMS_QCDScaleBounding[6];
  TH1D* histo_VVV_CMS_QCDScaleBounding[6];
  TH1D* histo_WG_CMS_QCDScaleBounding[6];
  TH1D* histo_WGS_CMS_QCDScaleBounding[6];
  TH1D* histo_Higgs_CMS_QCDScaleBounding[6];
  for(int nb=0; nb<6; nb++){
    histo_qqWW_CMS_QCDScaleBounding[nb]  = new TH1D(Form("histo_qqWW_QCDScale_f%d",nb),     Form("histo_qqWW_QCDScale_f%d",nb),nBinMVA, xbins);     histo_qqWW_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_ggWW_CMS_QCDScaleBounding[nb]  = new TH1D(Form("histo_ggWW_QCDScale_f%d",nb),     Form("histo_ggWW_QCDScale_f%d",nb),nBinMVA, xbins);     histo_ggWW_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_Top_CMS_QCDScaleBounding[nb]   = new TH1D(Form("histo_Top_QCDScale_f%d",nb),     Form("histo_Top_QCDScale_f%d",nb),nBinMVA, xbins);     histo_Top_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_DY_CMS_QCDScaleBounding[nb]    = new TH1D(Form("histo_DY_QCDScale_f%d",nb),     Form("histo_DY_QCDScale_f%d",nb),nBinMVA, xbins);     histo_DY_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_VV_CMS_QCDScaleBounding[nb]    = new TH1D(Form("histo_VV_QCDScale_f%d",nb),     Form("histo_VV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_VVV_CMS_QCDScaleBounding[nb]   = new TH1D(Form("histo_VVV_QCDScale_f%d",nb),     Form("histo_VVV_QCDScale_f%d",nb),nBinMVA, xbins);     histo_VVV_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WG_CMS_QCDScaleBounding[nb]    = new TH1D(Form("histo_WG_QCDScale_f%d",nb),     Form("histo_WG_QCDScale_f%d",nb),nBinMVA, xbins);     histo_WG_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_WGS_CMS_QCDScaleBounding[nb]   = new TH1D(Form("histo_WGS_QCDScale_f%d",nb),     Form("histo_WGS_QCDScale_f%d",nb),nBinMVA, xbins);     histo_WGS_CMS_QCDScaleBounding[nb]->Sumw2();
    histo_Higgs_CMS_QCDScaleBounding[nb] = new TH1D(Form("histo_Higgs_QCDScale_f%d",nb),      Form("histo_Higgs_QCDScale_f%d",nb),nBinMVA, xbins);      histo_Higgs_CMS_QCDScaleBounding[nb]->Sumw2();
  }
  TH1D* histo_qqWW_CMS_PDFBounding[100];
  TH1D* histo_ggWW_CMS_PDFBounding[100];
  TH1D* histo_Top_CMS_PDFBounding[100];
  TH1D* histo_DY_CMS_PDFBounding[100];
  TH1D* histo_VV_CMS_PDFBounding[100];
  TH1D* histo_VVV_CMS_PDFBounding[100];
  TH1D* histo_WG_CMS_PDFBounding[100];
  TH1D* histo_WGS_CMS_PDFBounding[100];
  TH1D* histo_Higgs_CMS_PDFBounding[100];
  for(int nb=0; nb<100; nb++){
    histo_qqWW_CMS_PDFBounding[nb]   = new TH1D(Form("histo_qqWW_PDF_f%d",nb),     Form("histo_qqWW_PDF_f%d",nb),    nBinMVA, xbins); histo_qqWW_CMS_PDFBounding[nb]->Sumw2();
    histo_ggWW_CMS_PDFBounding[nb]   = new TH1D(Form("histo_ggWW_PDF_f%d",nb),     Form("histo_ggWW_PDF_f%d",nb),    nBinMVA, xbins); histo_ggWW_CMS_PDFBounding[nb]->Sumw2();
    histo_Top_CMS_PDFBounding[nb]    = new TH1D(Form("histo_Top_PDF_f%d",nb),     Form("histo_Top_PDF_f%d",nb),    nBinMVA, xbins); histo_Top_CMS_PDFBounding[nb]->Sumw2();
    histo_DY_CMS_PDFBounding[nb]     = new TH1D(Form("histo_DY_PDF_f%d",nb),     Form("histo_DY_PDF_f%d",nb),    nBinMVA, xbins); histo_DY_CMS_PDFBounding[nb]->Sumw2();
    histo_VV_CMS_PDFBounding[nb]     = new TH1D(Form("histo_VV_PDF_f%d",nb),     Form("histo_VV_PDF_f%d",nb),    nBinMVA, xbins); histo_VV_CMS_PDFBounding[nb]->Sumw2();
    histo_VVV_CMS_PDFBounding[nb]    = new TH1D(Form("histo_VVV_PDF_f%d",nb),     Form("histo_VVV_PDF_f%d",nb),    nBinMVA, xbins); histo_VVV_CMS_PDFBounding[nb]->Sumw2();
    histo_WG_CMS_PDFBounding[nb]     = new TH1D(Form("histo_WG_PDF_f%d",nb),     Form("histo_WG_PDF_f%d",nb),    nBinMVA, xbins); histo_WG_CMS_PDFBounding[nb]->Sumw2();
    histo_WGS_CMS_PDFBounding[nb]    = new TH1D(Form("histo_WGS_PDF_f%d",nb),     Form("histo_WGS_PDF_f%d",nb),    nBinMVA, xbins); histo_WGS_CMS_PDFBounding[nb]->Sumw2();
    histo_Higgs_CMS_PDFBounding[nb]  = new TH1D(Form("histo_Higgs_PDF_f%d",nb),     Form("histo_Higgs_PDF_f%d",nb),    nBinMVA, xbins); histo_Higgs_CMS_PDFBounding[nb]->Sumw2();
  }

  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[nBinMVA];
  TH1D* histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[nBinMVA];
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingBinUp[nBinMVA];
  TH1D* histo_ggWW_CMS_MVAggWWStatBoundingBinDown[nBinMVA];
  TH1D* histo_Top_CMS_MVATopStatBoundingBinUp[nBinMVA];
  TH1D* histo_Top_CMS_MVATopStatBoundingBinDown[nBinMVA];
  TH1D* histo_DY_CMS_MVADYStatBoundingBinUp[nBinMVA];
  TH1D* histo_DY_CMS_MVADYStatBoundingBinDown[nBinMVA];
  TH1D* histo_VV_CMS_MVAVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VV_CMS_MVAVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinUp[nBinMVA];
  TH1D* histo_VVV_CMS_MVAVVVStatBoundingBinDown[nBinMVA];
  TH1D* histo_WG_CMS_MVAWGStatBoundingBinUp[nBinMVA];
  TH1D* histo_WG_CMS_MVAWGStatBoundingBinDown[nBinMVA];
  TH1D* histo_WGS_CMS_MVAWGSStatBoundingBinUp[nBinMVA];
  TH1D* histo_WGS_CMS_MVAWGSStatBoundingBinDown[nBinMVA];
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[nBinMVA];
  TH1D* histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[nBinMVA];
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[nBinMVA];
  TH1D* histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nBinMVA];
  TH1D* histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nBinMVA];
  for(int nb=0; nb<nBinMVA; nb++){
    histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[nb]	    = new TH1D(Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingBinUp[nb]      ->Sumw2();
    histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[nb]   = new TH1D(Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_qqWW_CMS_ww%s_MVAqqWWStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_qqWW_CMS_MVAqqWWStatBoundingBinDown[nb]    ->Sumw2();
    histo_ggWW_CMS_MVAggWWStatBoundingBinUp[nb]	    = new TH1D(Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingBinUp[nb]      ->Sumw2();
    histo_ggWW_CMS_MVAggWWStatBoundingBinDown[nb]   = new TH1D(Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_ggWW_CMS_ww%s_MVAggWWStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_ggWW_CMS_MVAggWWStatBoundingBinDown[nb]    ->Sumw2();
    histo_Top_CMS_MVATopStatBoundingBinUp[nb]	    = new TH1D(Form("histo_Top_CMS_ww%s_MVATopStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Top_CMS_ww%s_MVATopStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingBinUp[nb]      ->Sumw2();
    histo_Top_CMS_MVATopStatBoundingBinDown[nb]	    = new TH1D(Form("histo_Top_CMS_ww%s_MVATopStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_Top_CMS_ww%s_MVATopStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Top_CMS_MVATopStatBoundingBinDown[nb]    ->Sumw2();
    histo_DY_CMS_MVADYStatBoundingBinUp[nb]	    = new TH1D(Form("histo_DY_CMS_ww%s_MVADYStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_DY_CMS_ww%s_MVADYStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_DY_CMS_MVADYStatBoundingBinUp[nb]      ->Sumw2();
    histo_DY_CMS_MVADYStatBoundingBinDown[nb]	    = new TH1D(Form("histo_DY_CMS_ww%s_MVADYStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_DY_CMS_ww%s_MVADYStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_DY_CMS_MVADYStatBoundingBinDown[nb]    ->Sumw2();
    histo_VV_CMS_MVAVVStatBoundingBinUp[nb]	    = new TH1D(Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VV_CMS_MVAVVStatBoundingBinDown[nb]	    = new TH1D(Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VV_CMS_ww%s_MVAVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VV_CMS_MVAVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]	    = new TH1D(Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinUp[nb]      ->Sumw2();
    histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]	    = new TH1D(Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_VVV_CMS_ww%s_MVAVVVStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_VVV_CMS_MVAVVVStatBoundingBinDown[nb]    ->Sumw2();
    histo_WG_CMS_MVAWGStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WG_CMS_MVAWGStatBoundingBinUp[nb]      ->Sumw2();
    histo_WG_CMS_MVAWGStatBoundingBinDown[nb]	    = new TH1D(Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_WG_CMS_ww%s_MVAWGStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WG_CMS_MVAWGStatBoundingBinDown[nb]    ->Sumw2();
    histo_WGS_CMS_MVAWGSStatBoundingBinUp[nb]	    = new TH1D(Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WGS_CMS_MVAWGSStatBoundingBinUp[nb]      ->Sumw2();
    histo_WGS_CMS_MVAWGSStatBoundingBinDown[nb]	    = new TH1D(Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_WGS_CMS_ww%s_MVAWGSStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WGS_CMS_MVAWGSStatBoundingBinDown[nb]    ->Sumw2();
    histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[nb]   = new TH1D(Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingBinUp[nb]      ->Sumw2();
    histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[nb] = new TH1D(Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_WjetsM_CMS_ww%s_MVAWjetsMStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsM_CMS_MVAWjetsMStatBoundingBinDown[nb]    ->Sumw2();
    histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[nb]   = new TH1D(Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingBinUp[nb]      ->Sumw2();
    histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[nb] = new TH1D(Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_WjetsE_CMS_ww%s_MVAWjetsEStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_WjetsE_CMS_MVAWjetsEStatBoundingBinDown[nb]    ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]   = new TH1D(Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb), Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%s_Bin%dUp"      ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingBinUp[nb]      ->Sumw2();
    histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb] = new TH1D(Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb), Form("histo_Higgs_CMS_ww%s_MVAHiggsStatBounding_%s_Bin%dDown"    ,finalStateName,ECMsb.Data(),nb),nBinMVA, xbins); histo_Higgs_CMS_MVAHiggsStatBoundingBinDown[nb]    ->Sumw2();
  }

  TH1D* histo_qqWW_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_qqWW_%sUp",effMName)  , Form("histo_qqWW_%sUp",effMName)  , nBinMVA, xbins); histo_qqWW_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_qqWW_%sDown",effMName), Form("histo_qqWW_%sDown",effMName), nBinMVA, xbins); histo_qqWW_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_ggWW_%sUp",effMName)  , Form("histo_ggWW_%sUp",effMName)  , nBinMVA, xbins); histo_ggWW_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_ggWW_%sDown",effMName), Form("histo_ggWW_%sDown",effMName), nBinMVA, xbins); histo_ggWW_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_Top_%sUp",effMName)  , Form("histo_Top_%sUp",effMName)  , nBinMVA, xbins); histo_Top_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_Top_%sDown",effMName), Form("histo_Top_%sDown",effMName), nBinMVA, xbins); histo_Top_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_DY_%sUp",effMName)  , Form("histo_DY_%sUp",effMName)  , nBinMVA, xbins); histo_DY_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_DY_%sDown",effMName), Form("histo_DY_%sDown",effMName), nBinMVA, xbins); histo_DY_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_VV_%sUp",effMName)  , Form("histo_VV_%sUp",effMName)  , nBinMVA, xbins); histo_VV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_VV_%sDown",effMName), Form("histo_VV_%sDown",effMName), nBinMVA, xbins); histo_VV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effMName)  , Form("histo_VVV_%sUp",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effMName), Form("histo_VVV_%sDown",effMName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_WG_%sUp",effMName)  , Form("histo_WG_%sUp",effMName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_WG_%sDown",effMName), Form("histo_WG_%sDown",effMName), nBinMVA, xbins); histo_WG_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_WGS_%sUp",effMName)  , Form("histo_WGS_%sUp",effMName)  , nBinMVA, xbins); histo_WGS_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_WGS_%sDown",effMName), Form("histo_WGS_%sDown",effMName), nBinMVA, xbins); histo_WGS_CMS_MVALepEffMBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingUp   	 = new TH1D( Form("histo_Higgs_%sUp",effMName)  , Form("histo_Higgs_%sUp",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingDown 	 = new TH1D( Form("histo_Higgs_%sDown",effMName), Form("histo_Higgs_%sDown",effMName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_qqWW_%sAvg",effMName)  , Form("histo_qqWW_%sAvg",effMName)  , nBinMVA, xbins); histo_qqWW_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_ggWW_%sAvg",effMName)  , Form("histo_ggWW_%sAvg",effMName)  , nBinMVA, xbins); histo_ggWW_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_Top_%sAvg",effMName)  , Form("histo_Top_%sAvg",effMName)  , nBinMVA, xbins); histo_Top_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_DY_%sAvg",effMName)  , Form("histo_DY_%sAvg",effMName)  , nBinMVA, xbins); histo_DY_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_VV_%sAvg",effMName)  , Form("histo_VV_%sAvg",effMName)  , nBinMVA, xbins); histo_VV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effMName)  , Form("histo_VVV_%sAvg",effMName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_WG_%sAvg",effMName)  , Form("histo_WG_%sAvg",effMName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffMBoundingAvg   	 = new TH1D( Form("histo_WGS_%sAvg",effMName)  , Form("histo_WGS_%sAvg",effMName)  , nBinMVA, xbins); histo_WGS_CMS_MVALepEffMBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffMBoundingAvg    = new TH1D( Form("histo_Higgs_%sAvg",effMName)  , Form("histo_Higgs_%sAvg",effMName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffMBoundingAvg  ->Sumw2();

  TH1D* histo_qqWW_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_qqWW_%sUp",effEName)  , Form("histo_qqWW_%sUp",effEName)  , nBinMVA, xbins); histo_qqWW_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_qqWW_%sDown",effEName), Form("histo_qqWW_%sDown",effEName), nBinMVA, xbins); histo_qqWW_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_ggWW_%sUp",effEName)  , Form("histo_ggWW_%sUp",effEName)  , nBinMVA, xbins); histo_ggWW_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_ggWW_%sDown",effEName), Form("histo_ggWW_%sDown",effEName), nBinMVA, xbins); histo_ggWW_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_Top_%sUp",effEName)  , Form("histo_Top_%sUp",effEName)  , nBinMVA, xbins); histo_Top_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_Top_%sDown",effEName), Form("histo_Top_%sDown",effEName), nBinMVA, xbins); histo_Top_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_DY_%sUp",effEName)  , Form("histo_DY_%sUp",effEName)  , nBinMVA, xbins); histo_DY_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_DY_%sDown",effEName), Form("histo_DY_%sDown",effEName), nBinMVA, xbins); histo_DY_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_VV_%sUp",effEName)  , Form("histo_VV_%sUp",effEName)  , nBinMVA, xbins); histo_VV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_VV_%sDown",effEName), Form("histo_VV_%sDown",effEName), nBinMVA, xbins); histo_VV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_VVV_%sUp",effEName)  , Form("histo_VVV_%sUp",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_VVV_%sDown",effEName), Form("histo_VVV_%sDown",effEName), nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_WG_%sUp",effEName)  , Form("histo_WG_%sUp",effEName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_WG_%sDown",effEName), Form("histo_WG_%sDown",effEName), nBinMVA, xbins); histo_WG_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_WGS_%sUp",effEName)  , Form("histo_WGS_%sUp",effEName)  , nBinMVA, xbins); histo_WGS_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_WGS_%sDown",effEName), Form("histo_WGS_%sDown",effEName), nBinMVA, xbins); histo_WGS_CMS_MVALepEffEBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingUp   	 = new TH1D( Form("histo_Higgs_%sUp",effEName)  , Form("histo_Higgs_%sUp",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingDown 	 = new TH1D( Form("histo_Higgs_%sDown",effEName), Form("histo_Higgs_%sDown",effEName), nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_qqWW_%sAvg",effEName)  , Form("histo_qqWW_%sAvg",effEName)  , nBinMVA, xbins); histo_qqWW_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_ggWW_%sAvg",effEName)  , Form("histo_ggWW_%sAvg",effEName)  , nBinMVA, xbins); histo_ggWW_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_Top_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_Top_%sAvg",effEName)  , Form("histo_Top_%sAvg",effEName)  , nBinMVA, xbins); histo_Top_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_DY_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_DY_%sAvg",effEName)  , Form("histo_DY_%sAvg",effEName)  , nBinMVA, xbins); histo_DY_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VV_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_VV_%sAvg",effEName)  , Form("histo_VV_%sAvg",effEName)  , nBinMVA, xbins); histo_VV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_VVV_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_VVV_%sAvg",effEName)  , Form("histo_VVV_%sAvg",effEName)  , nBinMVA, xbins); histo_VVV_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WG_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_WG_%sAvg",effEName)  , Form("histo_WG_%sAvg",effEName)  , nBinMVA, xbins); histo_WG_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_WGS_CMS_MVALepEffEBoundingAvg   	 = new TH1D( Form("histo_WGS_%sAvg",effEName)  , Form("histo_WGS_%sAvg",effEName)  , nBinMVA, xbins); histo_WGS_CMS_MVALepEffEBoundingAvg  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVALepEffEBoundingAvg    = new TH1D( Form("histo_Higgs_%sAvg",effEName)  , Form("histo_Higgs_%sAvg",effEName)  , nBinMVA, xbins); histo_Higgs_CMS_MVALepEffEBoundingAvg  ->Sumw2();

  TH1D* histo_qqWW_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_qqWW_CMS_scale_metUp")  , Form("histo_qqWW_CMS_scale_metUp")  , nBinMVA, xbins); histo_qqWW_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_qqWW_CMS_scale_metDown"), Form("histo_qqWW_CMS_scale_metDown"), nBinMVA, xbins); histo_qqWW_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_ggWW_CMS_scale_metUp")  , Form("histo_ggWW_CMS_scale_metUp")  , nBinMVA, xbins); histo_ggWW_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_ggWW_CMS_scale_metDown"), Form("histo_ggWW_CMS_scale_metDown"), nBinMVA, xbins); histo_ggWW_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_Top_CMS_scale_metUp")  , Form("histo_Top_CMS_scale_metUp")  , nBinMVA, xbins); histo_Top_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_Top_CMS_scale_metDown"), Form("histo_Top_CMS_scale_metDown"), nBinMVA, xbins); histo_Top_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_DY_CMS_scale_metUp")  , Form("histo_DY_CMS_scale_metUp")  , nBinMVA, xbins); histo_DY_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_DY_CMS_scale_metDown"), Form("histo_DY_CMS_scale_metDown"), nBinMVA, xbins); histo_DY_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_VV_CMS_scale_metUp")  , Form("histo_VV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_VV_CMS_scale_metDown"), Form("histo_VV_CMS_scale_metDown"), nBinMVA, xbins); histo_VV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_VVV_CMS_scale_metUp")  , Form("histo_VVV_CMS_scale_metUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_VVV_CMS_scale_metDown"), Form("histo_VVV_CMS_scale_metDown"), nBinMVA, xbins); histo_VVV_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_WG_CMS_scale_metUp")  , Form("histo_WG_CMS_scale_metUp")  , nBinMVA, xbins); histo_WG_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_WG_CMS_scale_metDown"), Form("histo_WG_CMS_scale_metDown"), nBinMVA, xbins); histo_WG_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_WGS_CMS_scale_metUp")  , Form("histo_WGS_CMS_scale_metUp")  , nBinMVA, xbins); histo_WGS_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_WGS_CMS_scale_metDown"), Form("histo_WGS_CMS_scale_metDown"), nBinMVA, xbins); histo_WGS_CMS_MVAMETBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingUp   	 = new TH1D( Form("histo_Higgs_CMS_scale_metUp")  , Form("histo_Higgs_CMS_scale_metUp")  , nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAMETBoundingDown 	 = new TH1D( Form("histo_Higgs_CMS_scale_metDown"), Form("histo_Higgs_CMS_scale_metDown"), nBinMVA, xbins); histo_Higgs_CMS_MVAMETBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_qqWW_CMS_scale_jUp")  , Form("histo_qqWW_CMS_scale_jUp")  , nBinMVA, xbins); histo_qqWW_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_qqWW_CMS_scale_jDown"), Form("histo_qqWW_CMS_scale_jDown"), nBinMVA, xbins); histo_qqWW_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_ggWW_CMS_scale_jUp")  , Form("histo_ggWW_CMS_scale_jUp")  , nBinMVA, xbins); histo_ggWW_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_ggWW_CMS_scale_jDown"), Form("histo_ggWW_CMS_scale_jDown"), nBinMVA, xbins); histo_ggWW_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_Top_CMS_scale_jUp")  , Form("histo_Top_CMS_scale_jUp")  , nBinMVA, xbins); histo_Top_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_Top_CMS_scale_jDown"), Form("histo_Top_CMS_scale_jDown"), nBinMVA, xbins); histo_Top_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_DY_CMS_scale_jUp")  , Form("histo_DY_CMS_scale_jUp")  , nBinMVA, xbins); histo_DY_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_DY_CMS_scale_jDown"), Form("histo_DY_CMS_scale_jDown"), nBinMVA, xbins); histo_DY_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_VV_CMS_scale_jUp")  , Form("histo_VV_CMS_scale_jUp")  , nBinMVA, xbins); histo_VV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_VV_CMS_scale_jDown"), Form("histo_VV_CMS_scale_jDown"), nBinMVA, xbins); histo_VV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_VVV_CMS_scale_jUp")  , Form("histo_VVV_CMS_scale_jUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_VVV_CMS_scale_jDown"), Form("histo_VVV_CMS_scale_jDown"), nBinMVA, xbins); histo_VVV_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_WG_CMS_scale_jUp")  , Form("histo_WG_CMS_scale_jUp")  , nBinMVA, xbins); histo_WG_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_WG_CMS_scale_jDown"), Form("histo_WG_CMS_scale_jDown"), nBinMVA, xbins); histo_WG_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_WGS_CMS_scale_jUp")  , Form("histo_WGS_CMS_scale_jUp")  , nBinMVA, xbins); histo_WGS_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_WGS_CMS_scale_jDown"), Form("histo_WGS_CMS_scale_jDown"), nBinMVA, xbins); histo_WGS_CMS_MVAJESBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJESBoundingUp   	 = new TH1D( Form("histo_Higgs_CMS_scale_jUp")  , Form("histo_Higgs_CMS_scale_jUp")  , nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJESBoundingDown 	 = new TH1D( Form("histo_Higgs_CMS_scale_jDown"), Form("histo_Higgs_CMS_scale_jDown"), nBinMVA, xbins); histo_Higgs_CMS_MVAJESBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVAJERBoundingUp   	 = new TH1D( Form("histo_qqWW_CMS_jerUp")  , Form("histo_qqWW_CMS_jerUp")  , nBinMVA, xbins); histo_qqWW_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVAJERBoundingDown 	 = new TH1D( Form("histo_qqWW_CMS_jerDown"), Form("histo_qqWW_CMS_jerDown"), nBinMVA, xbins); histo_qqWW_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVAJERBoundingUp   	 = new TH1D( Form("histo_ggWW_CMS_jerUp")  , Form("histo_ggWW_CMS_jerUp")  , nBinMVA, xbins); histo_ggWW_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVAJERBoundingDown 	 = new TH1D( Form("histo_ggWW_CMS_jerDown"), Form("histo_ggWW_CMS_jerDown"), nBinMVA, xbins); histo_ggWW_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVAJERBoundingUp   	 = new TH1D( Form("histo_Top_CMS_jerUp")  , Form("histo_Top_CMS_jerUp")  , nBinMVA, xbins); histo_Top_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVAJERBoundingDown 	 = new TH1D( Form("histo_Top_CMS_jerDown"), Form("histo_Top_CMS_jerDown"), nBinMVA, xbins); histo_Top_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVAJERBoundingUp   	 = new TH1D( Form("histo_DY_CMS_jerUp")  , Form("histo_DY_CMS_jerUp")  , nBinMVA, xbins); histo_DY_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVAJERBoundingDown 	 = new TH1D( Form("histo_DY_CMS_jerDown"), Form("histo_DY_CMS_jerDown"), nBinMVA, xbins); histo_DY_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVAJERBoundingUp   	 = new TH1D( Form("histo_VV_CMS_jerUp")  , Form("histo_VV_CMS_jerUp")  , nBinMVA, xbins); histo_VV_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVAJERBoundingDown 	 = new TH1D( Form("histo_VV_CMS_jerDown"), Form("histo_VV_CMS_jerDown"), nBinMVA, xbins); histo_VV_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVAJERBoundingUp   	 = new TH1D( Form("histo_VVV_CMS_jerUp")  , Form("histo_VVV_CMS_jerUp")  , nBinMVA, xbins); histo_VVV_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVAJERBoundingDown 	 = new TH1D( Form("histo_VVV_CMS_jerDown"), Form("histo_VVV_CMS_jerDown"), nBinMVA, xbins); histo_VVV_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVAJERBoundingUp   	 = new TH1D( Form("histo_WG_CMS_jerUp")  , Form("histo_WG_CMS_jerUp")  , nBinMVA, xbins); histo_WG_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVAJERBoundingDown 	 = new TH1D( Form("histo_WG_CMS_jerDown"), Form("histo_WG_CMS_jerDown"), nBinMVA, xbins); histo_WG_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVAJERBoundingUp   	 = new TH1D( Form("histo_WGS_CMS_jerUp")  , Form("histo_WGS_CMS_jerUp")  , nBinMVA, xbins); histo_WGS_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVAJERBoundingDown 	 = new TH1D( Form("histo_WGS_CMS_jerDown"), Form("histo_WGS_CMS_jerDown"), nBinMVA, xbins); histo_WGS_CMS_MVAJERBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJERBoundingUp   	 = new TH1D( Form("histo_Higgs_CMS_jerUp")  , Form("histo_Higgs_CMS_jerUp")  , nBinMVA, xbins); histo_Higgs_CMS_MVAJERBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVAJERBoundingDown 	 = new TH1D( Form("histo_Higgs_CMS_jerDown"), Form("histo_Higgs_CMS_jerDown"), nBinMVA, xbins); histo_Higgs_CMS_MVAJERBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_MVABTAGBoundingUp   	 = new TH1D( Form("histo_qqWW_CMS_eff_b_b2016Up")  , Form("histo_qqWW_CMS_eff_b_b2016Up")  , nBinMVA, xbins); histo_qqWW_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_MVABTAGBoundingDown 	 = new TH1D( Form("histo_qqWW_CMS_eff_b_b2016Down"), Form("histo_qqWW_CMS_eff_b_b2016Down"), nBinMVA, xbins); histo_qqWW_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_MVABTAGBoundingUp   	 = new TH1D( Form("histo_ggWW_CMS_eff_b_b2016Up")  , Form("histo_ggWW_CMS_eff_b_b2016Up")  , nBinMVA, xbins); histo_ggWW_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_MVABTAGBoundingDown 	 = new TH1D( Form("histo_ggWW_CMS_eff_b_b2016Down"), Form("histo_ggWW_CMS_eff_b_b2016Down"), nBinMVA, xbins); histo_ggWW_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_MVABTAGBoundingUp   	 = new TH1D( Form("histo_Top_CMS_eff_b_b2016Up")  , Form("histo_Top_CMS_eff_b_b2016Up")  , nBinMVA, xbins); histo_Top_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_MVABTAGBoundingDown 	 = new TH1D( Form("histo_Top_CMS_eff_b_b2016Down"), Form("histo_Top_CMS_eff_b_b2016Down"), nBinMVA, xbins); histo_Top_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_MVABTAGBoundingUp   	 = new TH1D( Form("histo_DY_CMS_eff_b_b2016Up")  , Form("histo_DY_CMS_eff_b_b2016Up")  , nBinMVA, xbins); histo_DY_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_MVABTAGBoundingDown 	 = new TH1D( Form("histo_DY_CMS_eff_b_b2016Down"), Form("histo_DY_CMS_eff_b_b2016Down"), nBinMVA, xbins); histo_DY_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_MVABTAGBoundingUp   	 = new TH1D( Form("histo_VV_CMS_eff_b_b2016Up")  , Form("histo_VV_CMS_eff_b_b2016Up")  , nBinMVA, xbins); histo_VV_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_MVABTAGBoundingDown 	 = new TH1D( Form("histo_VV_CMS_eff_b_b2016Down"), Form("histo_VV_CMS_eff_b_b2016Down"), nBinMVA, xbins); histo_VV_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_MVABTAGBoundingUp   	 = new TH1D( Form("histo_VVV_CMS_eff_b_b2016Up")  , Form("histo_VVV_CMS_eff_b_b2016Up")  , nBinMVA, xbins); histo_VVV_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_MVABTAGBoundingDown 	 = new TH1D( Form("histo_VVV_CMS_eff_b_b2016Down"), Form("histo_VVV_CMS_eff_b_b2016Down"), nBinMVA, xbins); histo_VVV_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_MVABTAGBoundingUp   	 = new TH1D( Form("histo_WG_CMS_eff_b_b2016Up")  , Form("histo_WG_CMS_eff_b_b2016Up")  , nBinMVA, xbins); histo_WG_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_MVABTAGBoundingDown 	 = new TH1D( Form("histo_WG_CMS_eff_b_b2016Down"), Form("histo_WG_CMS_eff_b_b2016Down"), nBinMVA, xbins); histo_WG_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_MVABTAGBoundingUp   	 = new TH1D( Form("histo_WGS_CMS_eff_b_b2016Up")  , Form("histo_WGS_CMS_eff_b_b2016Up")  , nBinMVA, xbins); histo_WGS_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_MVABTAGBoundingDown 	 = new TH1D( Form("histo_WGS_CMS_eff_b_b2016Down"), Form("histo_WGS_CMS_eff_b_b2016Down"), nBinMVA, xbins); histo_WGS_CMS_MVABTAGBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_MVABTAGBoundingUp   	 = new TH1D( Form("histo_Higgs_CMS_eff_b_b2016Up")  , Form("histo_Higgs_CMS_eff_b_b2016Up")  , nBinMVA, xbins); histo_Higgs_CMS_MVABTAGBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_MVABTAGBoundingDown 	 = new TH1D( Form("histo_Higgs_CMS_eff_b_b2016Down"), Form("histo_Higgs_CMS_eff_b_b2016Down"), nBinMVA, xbins); histo_Higgs_CMS_MVABTAGBoundingDown->Sumw2();

  TH1D* histo_qqWW_CMS_PUBoundingUp   	         = new TH1D( Form("histo_qqWW_CMS_puUp")  , Form("histo_qqWW_CMS_puUp")  , nBinMVA, xbins); histo_qqWW_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_qqWW_CMS_PUBoundingDown 	         = new TH1D( Form("histo_qqWW_CMS_puDown"), Form("histo_qqWW_CMS_puDown"), nBinMVA, xbins); histo_qqWW_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_ggWW_CMS_PUBoundingUp   	         = new TH1D( Form("histo_ggWW_CMS_puUp")  , Form("histo_ggWW_CMS_puUp")  , nBinMVA, xbins); histo_ggWW_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_ggWW_CMS_PUBoundingDown 	         = new TH1D( Form("histo_ggWW_CMS_puDown"), Form("histo_ggWW_CMS_puDown"), nBinMVA, xbins); histo_ggWW_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_Top_CMS_PUBoundingUp   	         = new TH1D( Form("histo_Top_CMS_puUp")  , Form("histo_Top_CMS_puUp")  , nBinMVA, xbins); histo_Top_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_Top_CMS_PUBoundingDown 	         = new TH1D( Form("histo_Top_CMS_puDown"), Form("histo_Top_CMS_puDown"), nBinMVA, xbins); histo_Top_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_DY_CMS_PUBoundingUp   	         = new TH1D( Form("histo_DY_CMS_puUp")  , Form("histo_DY_CMS_puUp")  , nBinMVA, xbins); histo_DY_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_DY_CMS_PUBoundingDown 	         = new TH1D( Form("histo_DY_CMS_puDown"), Form("histo_DY_CMS_puDown"), nBinMVA, xbins); histo_DY_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_VV_CMS_PUBoundingUp   	         = new TH1D( Form("histo_VV_CMS_puUp")  , Form("histo_VV_CMS_puUp")  , nBinMVA, xbins); histo_VV_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_VV_CMS_PUBoundingDown 	         = new TH1D( Form("histo_VV_CMS_puDown"), Form("histo_VV_CMS_puDown"), nBinMVA, xbins); histo_VV_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingUp   	         = new TH1D( Form("histo_VVV_CMS_puUp")  , Form("histo_VVV_CMS_puUp")  , nBinMVA, xbins); histo_VVV_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_VVV_CMS_PUBoundingDown 	         = new TH1D( Form("histo_VVV_CMS_puDown"), Form("histo_VVV_CMS_puDown"), nBinMVA, xbins); histo_VVV_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WG_CMS_PUBoundingUp   	         = new TH1D( Form("histo_WG_CMS_puUp")  , Form("histo_WG_CMS_puUp")  , nBinMVA, xbins); histo_WG_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WG_CMS_PUBoundingDown 	         = new TH1D( Form("histo_WG_CMS_puDown"), Form("histo_WG_CMS_puDown"), nBinMVA, xbins); histo_WG_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_WGS_CMS_PUBoundingUp   	         = new TH1D( Form("histo_WGS_CMS_puUp")  , Form("histo_WGS_CMS_puUp")  , nBinMVA, xbins); histo_WGS_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_WGS_CMS_PUBoundingDown 	         = new TH1D( Form("histo_WGS_CMS_puDown"), Form("histo_WGS_CMS_puDown"), nBinMVA, xbins); histo_WGS_CMS_PUBoundingDown->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingUp   	         = new TH1D( Form("histo_Higgs_CMS_puUp")  , Form("histo_Higgs_CMS_puUp")  , nBinMVA, xbins); histo_Higgs_CMS_PUBoundingUp  ->Sumw2();
  TH1D* histo_Higgs_CMS_PUBoundingDown 	         = new TH1D( Form("histo_Higgs_CMS_puDown"), Form("histo_Higgs_CMS_puDown"), nBinMVA, xbins); histo_Higgs_CMS_PUBoundingDown->Sumw2();

  double bgdDecay[nSelTypes*4][histBins],weiDecay[nSelTypes*2][histBins];
  for(unsigned int i=0; i<nSelTypes*2; i++) {
    for(int j=0; j<histBins; j++) {
      bgdDecay[i][j] = 0.0; weiDecay[i][j] = 0.0; 
    }
  }

  const int numberCuts = 10;
  TH1D* hDWWLL[3];
  hDWWLL[0] = new TH1D("hDWWLL0", "hDWWLL0", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[1] = new TH1D("hDWWLL1", "hDWWLL1", numberCuts, -0.5, numberCuts-0.5);
  hDWWLL[2] = new TH1D("hDWWLL2", "hDWWLL2", numberCuts, -0.5, numberCuts-0.5);
  TString cutName[numberCuts] = {"ptl>20/20", "3rd lepton veto", "mll>12", "MET>20", "Z veto", "minPMET>20/45", "ptll>30/45", "btag-veto", "Nsoftmu=0", "Njets=0"};
  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infilenamev[ifile].Data());

    TFile *the_input_file = TFile::Open(infilenamev[ifile].Data());
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
    //TTree *the_input_all  = (TTree*)the_input_file->FindObjectAny("all");
    TTree *the_SelBit_tree= (TTree*)the_input_file->FindObjectAny("SelBit_tree");
    TTree *the_PDF_tree   = (TTree*)the_input_file->FindObjectAny("pdfReweight");

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
    eventMonteCarlo.setBranchAddresses(the_input_tree);
 
    if(useDYMVA == true){
      the_input_tree->SetBranchAddress("dymva", &dymva_);
      the_input_tree->SetBranchAddress("nlep", &nlep_ );
      the_input_tree->SetBranchAddress("njets", &njets_);
    }

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

    unsigned int selBit_= 0;
    the_SelBit_tree->SetBranchAddress("selBit", &selBit_);
    hDWWLL[0]->Scale(0.0);hDWWLL[1]->Scale(0.0);hDWWLL[2]->Scale(0.0);
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    if(the_input_tree->GetEntries() != the_SelBit_tree->GetEntries()) {printf("BIG SKIMMING FAILURE\n"); return;}
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_SelBit_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      if((selBit_ & 0x1<<whichSkim) == 0) continue;

      the_input_tree->GetEntry(i);
 
      int initPDFTag = 0;
      if((*eventMonteCarlo.pdfRwgt).size() == 0) initPDFTag = -1;

      int sumEvol[3] = {-1, -1, -1};
      Bool_t passFilter[numberCuts] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      Bool_t passPresel = kFALSE;
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

      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}
      }

      if(idLep.size()!=idTight.size()) {assert(1); return;}

      if(idLep.size() >= 2) passPresel = passPresel && kTRUE;
      else                  passPresel = passPresel && kFALSE;

      if(passPresel == kFALSE) continue;

      int typeSel = -1;
      if(idTight.size() >= 2){
        if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) typeSel = 0;
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) typeSel = 1;
        else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])!=TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]))         typeSel = 2;
        else {assert(1); return;}                                                                                                     ;
      }
      else {assert(1); printf("Not possible %d %d %d\n",(int)idLep.size(),(int)idTight.size(),goodIsTight); return;}                                                                                                     ;

      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 25 &&
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20) passFilter[0] = kTRUE;

      if((goodIsTight == idTight.size() || usePureMC == false) && idTight.size() == 2) passFilter[1] = kTRUE;

      TLorentzVector theMET;
      if(usePUPPI){
        theMET.SetPx((double)eventMet.metPuppi->Px());
        theMET.SetPy((double)eventMet.metPuppi->Py());
      } else {
        theMET.SetPx(((TLorentzVector*)(*eventMet.p4)[0])->Px());
        theMET.SetPy(((TLorentzVector*)(*eventMet.p4)[0])->Py());
      }
      double dPhiLepMETMin = 999.;double dPhiLepTrackMETMin = 999.;
      int signQ = 0;
      double systTotLep[2] = {1.0, 1.0}; // m/e
      for(unsigned nl=0; nl<idLep.size(); nl++){
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) systTotLep[0] = systTotLep[0] * 1.02;
        if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 11) systTotLep[1] = systTotLep[1] * 1.02;
        signQ = signQ + (int)(*eventLeptons.pdgId)[idLep[nl]]/TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]);
        if(dPhiLepMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(theMET)))
           dPhiLepMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(theMET));      
        if(dPhiLepTrackMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*eventMet.trackMet)))
           dPhiLepTrackMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*eventMet.trackMet));      
      }
      double PMET[3] = {theMET.Pt(),
                        (double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesUp])  ->Pt(),
			(double)((TLorentzVector*)(*eventMet.metSyst)[BareMet::JesDown])->Pt()};
      double PTrackMET[3] = {(double)eventMet.trackMet->Pt(),(double)eventMet.trackMet->Pt(),(double)eventMet.trackMet->Pt()};

      if(dPhiLepMETMin < TMath::Pi()/2) {PMET[0] = PMET[0] * sin(dPhiLepMETMin);
                                         PMET[1] = PMET[1] * sin(dPhiLepMETMin);
					 PMET[2] = PMET[2] * sin(dPhiLepMETMin);}
      if(dPhiLepTrackMETMin < TMath::Pi()/2) {PTrackMET[0] = PTrackMET[0] * sin(dPhiLepTrackMETMin);
                                              PTrackMET[1] = PTrackMET[1] * sin(dPhiLepTrackMETMin);
					      PTrackMET[2] = PTrackMET[2] * sin(dPhiLepTrackMETMin);}
      double minPMET[3] = {TMath::Min(PMET[0],PTrackMET[0]),TMath::Min(PMET[1],PTrackMET[1]),TMath::Min(PMET[2],PTrackMET[2])};

      bool goodSameSign = passFilter[0] && signQ != 0;
      passFilter[0] = passFilter[0] * (signQ == 0);

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 
      if(dilep.M() > 12) passFilter[2] = kTRUE;

      vector<int> idB,idC;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        if     (TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 5 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) idB.push_back(ngen0);
        else if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) == 4 && ((TLorentzVector*)(*eventMonteCarlo.p4)[ngen0])->Pt() > 15) idC.push_back(ngen0);
      }

      vector<int> idJet,idJesUp,idJesDown,idJerUp,idJerDown;
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double total_bjet_probLOOSE[2] = {1,1};double total_bjet_probLOOSEUP[2] = {1,1};double total_bjet_probLOOSEDOWN[2] = {1,1};
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      double theHT = ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() +
                     ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() +
                     theMET.Pt();
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() <= 20) continue;
	if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()) >= 4.7) continue;
        //bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;
	//if(((int)(*eventJets.selBits)[nj] & BareJets::JetLoose) != BareJets::JetLoose) continue;

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

	if(infilecatv[ifile] != 0){
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
          denBTagging[nJEta][nJPt][jetFlavor]++;
          if((float)(*eventJets.bDiscr)[nj] >= bTagCuts[0]) numBTaggingLOOSE[nJEta][nJPt][jetFlavor]++;

          double bjet_SFLOOSE = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFLOOSE = btagReaderLLOOSE.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
	  else                                  bjet_SFLOOSE = btagReaderBCLOOSE.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
          if(bjet_SFLOOSE == 0) bjet_SFLOOSE = 1;
          double bjet_SFLOOSEUP = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFLOOSEUP = btagReaderLLOOSEUP.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
	  else                                  bjet_SFLOOSEUP = btagReaderBCLOOSEUP.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
          if(bjet_SFLOOSEUP == 0) bjet_SFLOOSEUP = 1;
          double bjet_SFLOOSEDOWN = 1;
	  if(jetFlavor == BTagEntry::FLAV_UDSG) bjet_SFLOOSEDOWN = btagReaderLLOOSEDOWN.eval (jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
	  else                                  bjet_SFLOOSEDOWN = btagReaderBCLOOSEDOWN.eval(jetFlavor,TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()),TMath::Max(((TLorentzVector*)(*eventJets.p4)[nj])->Pt(),20.0));
          if(bjet_SFLOOSEDOWN == 0) bjet_SFLOOSEDOWN = 1;

	  if((float)(*eventJets.bDiscr)[nj] >= bTagCuts[0]){
	    total_bjet_probLOOSE[0] = total_bjet_probLOOSE[0] * jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor];
	    total_bjet_probLOOSE[1] = total_bjet_probLOOSE[1] * jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor] * bjet_SFLOOSE;
	    total_bjet_probLOOSEUP[0] = total_bjet_probLOOSEUP[0] * jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor];
	    total_bjet_probLOOSEUP[1] = total_bjet_probLOOSEUP[1] * jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor] * bjet_SFLOOSEUP;
	    total_bjet_probLOOSEDOWN[0] = total_bjet_probLOOSEDOWN[0] * jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor];
	    total_bjet_probLOOSEDOWN[1] = total_bjet_probLOOSEDOWN[1] * jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor] * bjet_SFLOOSEDOWN;
	  } else {
	    total_bjet_probLOOSE[0] = total_bjet_probLOOSE[0] * (1.0 - jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor]);
	    total_bjet_probLOOSE[1] = total_bjet_probLOOSE[1] * (1.0 - jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor] * bjet_SFLOOSE);
	    total_bjet_probLOOSEUP[0] = total_bjet_probLOOSEUP[0] * (1.0 - jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor]);
	    total_bjet_probLOOSEUP[1] = total_bjet_probLOOSEUP[1] * (1.0 - jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor] * bjet_SFLOOSEUP);
	    total_bjet_probLOOSEDOWN[0] = total_bjet_probLOOSEDOWN[0] * (1.0 - jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor]);
	    total_bjet_probLOOSEDOWN[1] = total_bjet_probLOOSEDOWN[1] * (1.0 - jetEpsBtagLOOSE[nJEta][nJPt][jetFlavor] * bjet_SFLOOSEDOWN);
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

        if(dPhiJetMET   == -1) dPhiJetMET   = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(theMET))*180./TMath::Pi();
        if(dPhiJetDiLep == -1) dPhiJetDiLep = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventJets.p4)[nj])))*180./TMath::Pi();

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 20 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 30) theHT = theHT + ((TLorentzVector*)(*eventJets.p4)[nj])->Pt();

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 30) idJet.push_back(nj);
      }

      if(useDYMVA == true){
        if(nlep_ != idLep.size()) {printf("PROBLEM nlep %d != %d\n",(int)nlep_,(int)idLep.size()); assert(1); return;}
        if(njets_ != idJet.size()) {printf("PROBLEM njet %d != %d\n",(int)njets_,(int)idJet.size()); assert(1); return;}
      }

      int typePair = 0;
      if     (idLep.size() == 2 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])) typePair = 1;

      if(theMET.Pt() > 20) passFilter[3] = kTRUE;
      passFilter[4] = kTRUE;
      if(minPMET[0] > 20) passFilter[5] = kTRUE;
      if(dilep.Pt() > 30) passFilter[6] = kTRUE;
      if(typePair == 1){
        if(TMath::Abs(dilep.M()-91.1876) <= 15.0) passFilter[4] = kFALSE;
        if     (useDYMVA == false && minPMET[0] <= 45) passFilter[5] = kFALSE;
        else if(useDYMVA == true && dymva_ < 0.3) passFilter[5] = kFALSE;
        if(dilep.Pt() <= 45) passFilter[6] = kFALSE;
      }
      if(bDiscrMax < bTagCuts[0]) passFilter[7] = kTRUE;
      //if(idSoft.size() == 0) passFilter[8] = kTRUE;
      passFilter[8] = kTRUE;
      if(idJet.size() == nJetsType) passFilter[9] = kTRUE;

      bool totalSel = kTRUE;
      for(int isel=0; isel<numberCuts; isel++) {
        totalSel = totalSel && passFilter[isel];
	if(totalSel == kTRUE) sumEvol[typeSel]++;
      }
      /*if(totalSel){
      printf("LLL %d %d %llu %d %f %f %f %f %f %f %f %f %f %f %f ",eventEvent.runNum,eventEvent.lumiNum,eventEvent.eventNum,typeSel,
      ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Phi(),
      ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Phi(),
      theMET.Pt(),(double)eventMet.trackMet->Pt(),dPhiLepMETMin,dPhiLepTrackMETMin,minPMET[0]);
      if(idJet.size()>0) printf("%f %f %f\n",((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Pt(),((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta(),((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Phi());
      else               printf("0.0 0.0 0.0\n");
      }*/
      bool passSameSignRegion = goodSameSign  && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5] &&  passFilter[6] &&   passFilter[7] &&  passFilter[8]  && passFilter[9];
      bool passTopRegion      = passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5] &&  passFilter[6] && (!passFilter[7] || !passFilter[8]) && passFilter[9];
      bool passDYRegion       = passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5] && !passFilter[6] &&   passFilter[7] &&  passFilter[8]  && passFilter[9];

      bool passNoJetCutRegion = passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5] &&  passFilter[6] &&   passFilter[7] &&  passFilter[8];

      bool passAllCuts[nSelTypes] = {totalSel, passSameSignRegion, passTopRegion, passDYRegion};

      bool passSystCuts[nSystTypes] = {
        passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5]                                           && passFilter[6] && passFilter[7] && passFilter[8] && idJesUp.size()   == nJetsType,
        passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5]                                           && passFilter[6] && passFilter[7] && passFilter[8] && idJesDown.size() == nJetsType,
        passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && (minPMET[1] > 20 && (minPMET[1] > 45 || typePair == 0)) && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9],
        passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && (minPMET[2] > 20 && (minPMET[2] > 45 || typePair == 0)) && passFilter[6] && passFilter[7] && passFilter[8] && passFilter[9],
        passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5]                                           && passFilter[6] && passFilter[7] && passFilter[8] && idJerUp.size()   == nJetsType,
        passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5]                                           && passFilter[6] && passFilter[7] && passFilter[8] && idJerDown.size() == nJetsType
      };

      if     (theControlRegion == 1){
        passSystCuts[JESUP]   = passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5]                                           && passFilter[6] && (!passFilter[7] || !passFilter[8]) && idJesUp.size()   == nJetsType;
        passSystCuts[JESDOWN] = passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5]                                           && passFilter[6] && (!passFilter[7] || !passFilter[8]) && idJesDown.size() == nJetsType;
        passSystCuts[METUP]   = passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && (minPMET[1] > 20 && (minPMET[1] > 45 || typePair == 0)) && passFilter[6] && (!passFilter[7] || !passFilter[8]) && passFilter[9];
        passSystCuts[METDOWN] = passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && (minPMET[2] > 20 && (minPMET[2] > 45 || typePair == 0)) && passFilter[6] && (!passFilter[7] || !passFilter[8]) && passFilter[9];
        passSystCuts[JERUP]   = passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5]                                           && passFilter[6] && (!passFilter[7] || !passFilter[8]) && idJerUp.size()   == nJetsType;
        passSystCuts[JERDOWN] = passFilter[0] && passFilter[1] && passFilter[2] && passFilter[3] && passFilter[4] && passFilter[5]                                           && passFilter[6] && (!passFilter[7] || !passFilter[8]) && idJerDown.size() == nJetsType;
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
      for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	  if(isGenDupl[ngen] == 1) continue;
          genLep++;
      }

      double thePtwwWeight[5] = {1.0,1.0,1.0,1.0,1.0};
      if(infilecatv[ifile] == 1 && wBoson.size() == 2){
        TLorentzVector wwSystem(( ( *(TLorentzVector*)(eventMonteCarlo.p4->At(wBoson[0])) ) + ( *(TLorentzVector*)(eventMonteCarlo.p4->At(wBoson[1])) ) )); 
        Int_t nptwwbin[5] = {fhDWWPtRatio	   ->GetXaxis()->FindBin(TMath::Min(wwSystem.Pt(),499.999)),
	                     fhDWWPtRatio_scaleup  ->GetXaxis()->FindBin(TMath::Min(wwSystem.Pt(),499.999)),
	                     fhDWWPtRatio_scaledown->GetXaxis()->FindBin(TMath::Min(wwSystem.Pt(),499.999)),
	                     fhDWWPtRatio_resumup  ->GetXaxis()->FindBin(TMath::Min(wwSystem.Pt(),499.999)),
	                     fhDWWPtRatio_resumdown->GetXaxis()->FindBin(TMath::Min(wwSystem.Pt(),499.999))};
        thePtwwWeight[0] = fhDWWPtRatio 	 ->GetBinContent(nptwwbin[0]);
	thePtwwWeight[1] = fhDWWPtRatio_scaleup  ->GetBinContent(nptwwbin[1]);
	thePtwwWeight[2] = fhDWWPtRatio_scaledown->GetBinContent(nptwwbin[2]);
	thePtwwWeight[3] = fhDWWPtRatio_resumup  ->GetBinContent(nptwwbin[3]);
	thePtwwWeight[4] = fhDWWPtRatio_resumdown->GetBinContent(nptwwbin[4]);
      }
      vector<int> isGenLep; unsigned int goodIsGenLep = 0; unsigned int muGenLep = 0;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        bool isGenLepton = false;
        bool isGenMuon = false;
        for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	  if(isGenDupl[ngen] == 1) continue;
          if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.3) {
	    isGenLepton = true;
            if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) == 13) isGenMuon = true;
	    break;
	  }
	}
	if(isGenLepton == true) {isGenLep.push_back(1); goodIsGenLep++; if(isGenMuon == true) muGenLep++;}
	else                    {isGenLep.push_back(0);}
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

      // fake rate
      int type2l = -1;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == 13 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]) == 13) type2l = 0;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == 13 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]) == 11) type2l = 1;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == 11 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]) == 13) type2l = 2;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == 11 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]) == 11) type2l = 3;
      else {printf("IMPOSSIBLE TYPE2L\n");}
      int nFakeCount = 0;
      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false){
        if     (theCategory == 9){ // remove W+jets from MC
          fakeSF = 0.0;
        }
        else if(theCategory == 4 && goodIsTight != idTight.size()){ // remove Z+jets from MC as fakeable objects
          fakeSF = 0.0;
        }
        else if((infilecatv[ifile] == 0 || infilecatv[ifile] == 7 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add W+jets from data or W+gamma from MC
          unsigned int typeFakeLepton[2] = {0,0};
	  for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) nFakeCount = nFakeCount + 1;
	    else                                                        nFakeCount = nFakeCount + 2;

	    fakeSF = fakeSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    theCategory = 9;
	    if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == 13) typeFakeLepton[0]++;
	    else                                                        typeFakeLepton[1]++;
          }
          if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) fakeSF =  1.0 * fakeSF; // double fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) fakeSF = -1.0 * fakeSF; // single fake, MC
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) fakeSF = -1.0 * fakeSF; // double fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) fakeSF =  1.0 * fakeSF; // single fake, data
          if(typeFakeLepton[0] < typeFakeLepton[1]) theCategory = 10;
        }
        else if(infilecatv[ifile] != 0 && infilecatv[ifile] != 7 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
	  fakeSF = 1.0;
        }
        else if(infilecatv[ifile] == 0 || infilecatv[ifile] == 7){ // data or W+gamma with all good leptons
	  fakeSF = 1.0;
        }
	else {
	  printf("PROBLEM: %d %d %d %d %d\n",infilecatv[ifile],goodIsGenLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
	  assert(0);
	}
      }
      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale*trigEff*thePtwwWeight[0];
      //double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF*theMCPrescale*trigEff;

      // top-quark estimation
      if(infilecatv[ifile] == 3 && useTopCR == false) totalWeight = totalWeight * topNorm[TMath::Min((int)nJetsType,2)];

      // z pt correction (not applied anymore)
      //if(infilecatv[ifile] == 4 && zBoson.size() == 1) totalWeight = totalWeight * zpt_correction(((TLorentzVector*)(*eventMonteCarlo.p4)[zBoson[0]])->Pt(), 0);

      // wg* veto applied in wg events
      if(infilecatv[ifile] == 7 && applyGStarVeto == true && genLep >= 3) totalWeight = 0;
      // scale factor on wg* MC events
      if(infilecatv[ifile] == 8) totalWeight = totalWeight * 1.2;

      // Btag scale factor
      if(useTopCR == true) totalWeight = totalWeight * total_bjet_probLOOSE[1]/total_bjet_probLOOSE[0];

      double btagCorr[2] = {(total_bjet_probLOOSEUP[1]  /total_bjet_probLOOSEUP[0]  )/(total_bjet_probLOOSE[1]/total_bjet_probLOOSE[0]),
                            (total_bjet_probLOOSEDOWN[1]/total_bjet_probLOOSEDOWN[0])/(total_bjet_probLOOSE[1]/total_bjet_probLOOSE[0])};

      //if(passAllCuts[SIGSEL]) printf("AAA %f %f %f %f %f %f %f %f\n",totalWeight,mcWeight,theLumi,puWeight,effSF,fakeSF,thePtwwWeight[0],total_bjet_probLOOSE[1]/total_bjet_probLOOSE[0]);

      if(totalWeight == 0) continue;
      // end event weighting

      //if(theCategory == 9 && muGenLep >= 1) theCategory = 10;

      if(passAllCuts[SIGSEL] && infilecatv[ifile] == 0) totalFakeDataCount[type2l][nFakeCount] = totalFakeDataCount[type2l][nFakeCount] + 1;

      for(int nl=0; nl <=sumEvol[typeSel]; nl++) if(fakeSF == 1) hDWWLL[typeSel]->Fill((double)nl,totalWeight);

      for(unsigned int i=0; i<nSelTypes; i++) {
        if(passAllCuts[i]) {
          bgdDecay[i+typePair*nSelTypes][theCategory] += totalWeight;
          weiDecay[i+typePair*nSelTypes][theCategory] += totalWeight*totalWeight;
        }
      }

      for(int thePlot=0; thePlot<allPlots; thePlot++){
        double theVar = 0.0;
        bool makePlot = false;
        if     (thePlot ==  0 && passAllCuts[SIGSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),199.999);}
        else if(thePlot ==  1 && passAllCuts[SIGSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),199.999);}
        else if(thePlot ==  2 && passAllCuts[SIGSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
        else if(thePlot ==  3 && passAllCuts[SIGSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]]))*180./TMath::Pi();}
        else if(thePlot ==  4 && passAllCuts[SIGSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Min(theHT,499.999);}
        else if(thePlot ==  5 && passAllCuts[SIGSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Min(dilep.M(),799.999);}

        else if(thePlot ==  6 && passAllCuts[SSSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),199.999);}
        else if(thePlot ==  7 && passAllCuts[SSSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),199.999);}
        else if(thePlot ==  8 && passAllCuts[SSSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
        else if(thePlot ==  9 && passAllCuts[SSSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]]))*180./TMath::Pi();}
        else if(thePlot == 10 && passAllCuts[SSSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Min(theHT,499.999);}
        else if(thePlot == 11 && passAllCuts[SSSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Min(dilep.M(),799.999);}

        else if(thePlot == 12 && passAllCuts[TOPSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),199.999);}
        else if(thePlot == 13 && passAllCuts[TOPSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),199.999);}
        else if(thePlot == 14 && passAllCuts[TOPSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
        else if(thePlot == 15 && passAllCuts[TOPSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]]))*180./TMath::Pi();}
        else if(thePlot == 16 && passAllCuts[TOPSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Min(theHT,499.999);}
        else if(thePlot == 17 && passAllCuts[TOPSEL] && typeSel == 2) {makePlot = true;theVar = TMath::Min(dilep.M(),799.999);}

        else if(thePlot == 18 && passAllCuts[DYSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),199.999);}
        else if(thePlot == 19 && passAllCuts[DYSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),199.999);}
        else if(thePlot == 20 && passAllCuts[DYSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);}
        else if(thePlot == 21 && passAllCuts[DYSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->DeltaPhi(*(TLorentzVector*)(*eventLeptons.p4)[idLep[1]]))*180./TMath::Pi();}
        else if(thePlot == 22 && passAllCuts[DYSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Min(theHT,499.999);}
        else if(thePlot == 23 && passAllCuts[DYSEL] && typeSel == 2)  {makePlot = true;theVar = TMath::Min(dilep.M(),799.999);}

        else if(thePlot == 24 && passNoJetCutRegion && typeSel == 2)  {makePlot = true;theVar = TMath::Min(theHT,499.999);}

        if(makePlot) histo[thePlot][theCategory]->Fill(theVar,totalWeight);
      }

      if(1) {
	double MVAVar = (double)typePair;
	if(isShapeAna){
          if(typePair == 1) {passAllCuts[SIGSEL] = false; passAllCuts[TOPSEL] = false;}
          MVAVar =  TMath::Min(dilep.M(),xbins[nBinMVA]-0.001);
        }

        // Avoid QCD scale and PDF weights that are anomalous high
        double maxQCDscale = (TMath::Abs((double)eventMonteCarlo.r1f2)+TMath::Abs((double)eventMonteCarlo.r1f5)+TMath::Abs((double)eventMonteCarlo.r2f1)+
                              TMath::Abs((double)eventMonteCarlo.r2f2)+TMath::Abs((double)eventMonteCarlo.r5f1)+TMath::Abs((double)eventMonteCarlo.r5f5))/6.0;
        double PDFAvg = 0.0;
        if(infilecatv[ifile] != 0 && ((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1))){
          if(initPDFTag != -1)
          for(int npdf=0; npdf<100; npdf++) PDFAvg = PDFAvg + TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag]);
          PDFAvg = PDFAvg/100.0;
        }

	if(infilecatv[ifile] == 0){
          if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {

            int typeDiSel = -1;
            if     (idTight[0] == 1 && idTight[1] == 1) {typeDiSel = 0;}
            else if(idTight[0] == 0 && idTight[1] == 1) {typeDiSel = 1;}
            else if(idTight[0] == 1 && idTight[1] == 0) {typeDiSel = 2;}
            else if(idTight[0] == 0 && idTight[1] == 0) {typeDiSel = 3;}
            double totalFakeWeight = 
            fakePromptRateFactor(
            ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]),
            ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]),
            typeLepSel.Data(),typeDiSel);

            histo_FakeData->Fill(MVAVar,totalFakeWeight);
          }
	}

        if     (theCategory == 0){
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) histo_Data->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 1){ // qqWW
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
	     histo_qqWW->Fill(MVAVar,totalWeight);

             histo_qqWW_CMS_MVAWW_nlo      ->Fill(MVAVar,totalWeight);
             histo_qqWW_CMS_MVAWW_qup_nlo  ->Fill(MVAVar,totalWeight*thePtwwWeight[1]);
             histo_qqWW_CMS_MVAWW_qdown_nlo->Fill(MVAVar,totalWeight*thePtwwWeight[2]);
             histo_qqWW_CMS_MVAWW_sup_nlo  ->Fill(MVAVar,totalWeight*thePtwwWeight[3]);
             histo_qqWW_CMS_MVAWW_sdown_nlo->Fill(MVAVar,totalWeight*thePtwwWeight[4]);

	     histo_qqWW_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_qqWW_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_qqWW_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_qqWW_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_qqWW_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_qqWW_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_qqWW_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_qqWW_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_qqWW_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_qqWW_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_qqWW_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_qqWW_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_qqWW_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_qqWW_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_qqWW_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_qqWW_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_qqWW_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_qqWW_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_qqWW_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_qqWW_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_qqWW_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_qqWW_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERUP])  histo_qqWW_CMS_MVAJERBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERDOWN])histo_qqWW_CMS_MVAJERBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 2){ // ggWW
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
	     histo_ggWW->Fill(MVAVar,totalWeight);

	     histo_ggWW_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_ggWW_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_ggWW_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_ggWW_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_ggWW_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_ggWW_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_ggWW_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_ggWW_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_ggWW_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ggWW_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_ggWW_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_ggWW_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_ggWW_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_ggWW_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_ggWW_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_ggWW_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_ggWW_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_ggWW_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_ggWW_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_ggWW_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_ggWW_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_ggWW_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERUP])  histo_ggWW_CMS_MVAJERBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERDOWN])histo_ggWW_CMS_MVAJERBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 3){ // Top
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
	     histo_Top->Fill(MVAVar,totalWeight);

	     histo_Top_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_Top_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_Top_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_Top_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_Top_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_Top_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_Top_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_Top_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_Top_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Top_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Top_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Top_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Top_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_Top_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_Top_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_Top_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_Top_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_Top_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_Top_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_Top_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_Top_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_Top_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERUP])  histo_Top_CMS_MVAJERBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERDOWN])histo_Top_CMS_MVAJERBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 4){ // DY
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
	     histo_DY->Fill(MVAVar,totalWeight);

	     histo_DY_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_DY_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_DY_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_DY_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_DY_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_DY_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_DY_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_DY_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_DY_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_DY_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_DY_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_DY_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_DY_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_DY_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_DY_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_DY_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_DY_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_DY_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_DY_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_DY_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_DY_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_DY_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERUP])  histo_DY_CMS_MVAJERBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERDOWN])histo_DY_CMS_MVAJERBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 5){ // VV
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
	     histo_VV->Fill(MVAVar,totalWeight);

	     histo_VV_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_VV_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_VV_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_VV_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_VV_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_VV_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_VV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_VV_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_VV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_VV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_VV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_VV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_VV_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_VV_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_VV_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_VV_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_VV_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_VV_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_VV_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_VV_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERUP])  histo_VV_CMS_MVAJERBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERDOWN])histo_VV_CMS_MVAJERBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 6){ // VVV
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
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
             histo_VVV_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_VVV_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_VVV_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_VVV_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_VVV_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_VVV_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_VVV_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_VVV_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_VVV_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_VVV_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_VVV_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_VVV_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERUP])  histo_VVV_CMS_MVAJERBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERDOWN])histo_VVV_CMS_MVAJERBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 7){ // WG
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
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
             histo_WG_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WG_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WG_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_WG_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_WG_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_WG_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_WG_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_WG_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_WG_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_WG_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_WG_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_WG_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_WG_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_WG_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERUP])  histo_WG_CMS_MVAJERBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERDOWN])histo_WG_CMS_MVAJERBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 8){ // WGS
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
	     histo_WGS->Fill(MVAVar,totalWeight);

	     histo_WGS_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_WGS_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_WGS_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_WGS_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_WGS_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_WGS_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_WGS_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_WGS_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_WGS_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WGS_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_WGS_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_WGS_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_WGS_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_WGS_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_WGS_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_WGS_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_WGS_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_WGS_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_WGS_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_WGS_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_WGS_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_WGS_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERUP])  histo_WGS_CMS_MVAJERBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERDOWN])histo_WGS_CMS_MVAJERBoundingDown->Fill(MVAVar,totalWeight);
        }
        else if(theCategory == 9){ // WjetsM
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
	     histo_WjetsM->Fill(MVAVar,totalWeight);
	  }
        }
        else if(theCategory == 10){ // WjetsE
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
	     histo_WjetsE->Fill(MVAVar,totalWeight);
	  }
        }
        else if(theCategory == 11){ // Higgs
	  if((passAllCuts[SIGSEL] && theControlRegion == 0) || (passAllCuts[TOPSEL] && theControlRegion == 1)) {
	     histo_Higgs->Fill(MVAVar,totalWeight);

	     histo_Higgs_CMS_QCDScaleBounding[0]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f2)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[1]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r1f5)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[2]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f1)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[3]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r2f2)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[4]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f1)/maxQCDscale);
	     histo_Higgs_CMS_QCDScaleBounding[5]->Fill(MVAVar,totalWeight*TMath::Abs((double)eventMonteCarlo.r5f5)/maxQCDscale);
	     if(initPDFTag != -1)
	     for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight*TMath::Abs((double)(*eventMonteCarlo.pdfRwgt)[npdf+initPDFTag])/PDFAvg);
	     else
	     for(int npdf=0; npdf<100; npdf++) histo_Higgs_CMS_PDFBounding[npdf]->Fill(MVAVar,totalWeight);
             histo_Higgs_CMS_MVALepEffMBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Higgs_CMS_MVALepEffEBoundingAvg ->Fill(MVAVar,totalWeight*1.00);
             histo_Higgs_CMS_MVALepEffMBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[0]);
             histo_Higgs_CMS_MVALepEffEBoundingUp  ->Fill(MVAVar,totalWeight*systTotLep[1]);
             histo_Higgs_CMS_MVALepEffMBoundingDown->Fill(MVAVar,totalWeight/systTotLep[0]);
             histo_Higgs_CMS_MVALepEffEBoundingDown->Fill(MVAVar,totalWeight/systTotLep[1]);
             histo_Higgs_CMS_PUBoundingUp  ->Fill(MVAVar,totalWeight*puWeightUp  /puWeight);
             histo_Higgs_CMS_PUBoundingDown->Fill(MVAVar,totalWeight*puWeightDown/puWeight);
             histo_Higgs_CMS_MVABTAGBoundingUp  ->Fill(MVAVar,totalWeight*btagCorr[0]);
             histo_Higgs_CMS_MVABTAGBoundingDown->Fill(MVAVar,totalWeight*btagCorr[1]);
	  }
          if(passSystCuts[JESUP])  histo_Higgs_CMS_MVAJESBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JESDOWN])histo_Higgs_CMS_MVAJESBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[METUP])  histo_Higgs_CMS_MVAMETBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[METDOWN])histo_Higgs_CMS_MVAMETBoundingDown->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERUP])  histo_Higgs_CMS_MVAJERBoundingUp  ->Fill(MVAVar,totalWeight);
          if(passSystCuts[JERDOWN])histo_Higgs_CMS_MVAJERBoundingDown->Fill(MVAVar,totalWeight);
        }
	else {
	  printf("CATEGORY PROBLEM!\n"); return;
	}
      } // making data cards

    }
    printf("                   mm           ee           em\n");
    for(int nc=0; nc<numberCuts; nc++){
      printf("(%15s): %10.2f   %10.2f   %10.2f\n",cutName[nc].Data(),hDWWLL[0]->GetBinContent(nc+1),hDWWLL[1]->GetBinContent(nc+1),hDWWLL[2]->GetBinContent(nc+1));
    }
    the_input_file->Close();
  } // end of chain

  printf("----------------------totalFakeDataCount--------------------------------\n");
  for(int ni=0; ni<4; ni++) {
    printf("(%d): ",ni);
    for(int nj=0; nj<5; nj++) printf("%6.1f ",totalFakeDataCount[ni][nj]);
    printf("\n");
  }
  printf("----------------------totalFakeData--------------------------------\n");
  printf("total: %.2f\n",histo_FakeData->GetSumOfWeights());
  double sumDataFakes[2] = {0.0, 0.0};
  for(int np=1; np<=histo_FakeData->GetNbinsX(); np++) {
    printf(" %.2f +/-  %.2f",histo_FakeData->GetBinContent(np),histo_FakeData->GetBinError(np));
    sumDataFakes[0] = sumDataFakes[0] + histo_FakeData->GetBinContent(np);
    sumDataFakes[1] = sumDataFakes[1] + histo_FakeData->GetBinError(np)*histo_FakeData->GetBinError(np);
  }
  printf("\n");
  printf("sumDataFakes: %f +/- %f\n",sumDataFakes[0],sqrt(sumDataFakes[1]));

  printf("                    em                     ee/mm                     ll\n");
  printf("----------------------------------------------------------------------------------\n");
  for(int ns=0; ns<nSelTypes; ns++) {
    printf("Selection: %s\n",selTypeName[ns].Data());
    double sumEventsBckType[3] = {0,0,0}; double sumEventsBckTypeE[3] = {0,0,0};
    double sumEventsSigType[3] = {0,0,0}; double sumEventsSigTypeE[3] = {0,0,0};
    for(int np=0; np<histBins; np++) {       
       bgdDecay[ns+nSelTypes*2][np] = bgdDecay[ns+nSelTypes*0][np] + bgdDecay[ns+nSelTypes*1][np];
       weiDecay[ns+nSelTypes*2][np] = weiDecay[ns+nSelTypes*0][np] + weiDecay[ns+nSelTypes*1][np];
       printf("(%6s): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
       processName[np].Data(),bgdDecay[ns+nSelTypes*0][np],sqrt(weiDecay[ns+nSelTypes*0][np]),bgdDecay[ns+nSelTypes*1][np],sqrt(weiDecay[ns+nSelTypes*1][np]),
                              bgdDecay[ns+nSelTypes*2][np],sqrt(weiDecay[ns+nSelTypes*2][np]));
       if(np!=0 && np!=1 && np !=2){
         sumEventsBckType[0] = sumEventsBckType[0] + bgdDecay[ns+nSelTypes*0][np]; sumEventsBckTypeE[0] = sumEventsBckTypeE[0] + weiDecay[ns+nSelTypes*0][np];
         sumEventsBckType[1] = sumEventsBckType[1] + bgdDecay[ns+nSelTypes*1][np]; sumEventsBckTypeE[1] = sumEventsBckTypeE[1] + weiDecay[ns+nSelTypes*1][np];
         sumEventsBckType[2] = sumEventsBckType[2] + bgdDecay[ns+nSelTypes*2][np]; sumEventsBckTypeE[2] = sumEventsBckTypeE[2] + weiDecay[ns+nSelTypes*2][np];
       }
       if(np==1 || np==2){
         sumEventsSigType[0] = sumEventsSigType[0] + bgdDecay[ns+nSelTypes*0][np]; sumEventsSigTypeE[0] = sumEventsSigTypeE[0] + weiDecay[ns+nSelTypes*0][np];
         sumEventsSigType[1] = sumEventsSigType[1] + bgdDecay[ns+nSelTypes*1][np]; sumEventsSigTypeE[1] = sumEventsSigTypeE[1] + weiDecay[ns+nSelTypes*1][np];
         sumEventsSigType[2] = sumEventsSigType[2] + bgdDecay[ns+nSelTypes*2][np]; sumEventsSigTypeE[2] = sumEventsSigTypeE[2] + weiDecay[ns+nSelTypes*2][np];
       }
    }
    printf("(...bkg): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
           sumEventsBckType[0],sqrt(sumEventsBckTypeE[0]),sumEventsBckType[1],sqrt(sumEventsBckTypeE[1]),sumEventsBckType[2],sqrt(sumEventsBckTypeE[2]));
    printf("(....WW): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
           sumEventsSigType[0],sqrt(sumEventsSigTypeE[0]),sumEventsSigType[1],sqrt(sumEventsSigTypeE[1]),sumEventsSigType[2],sqrt(sumEventsSigTypeE[2]));
    printf("(.total): %9.2f +/- %7.2f | %9.2f +/- %7.2f | %9.2f +/- %7.2f\n",
           sumEventsBckType[0]+sumEventsSigType[0],sqrt(sumEventsBckTypeE[0]+sumEventsSigTypeE[0]),sumEventsBckType[1]+sumEventsSigType[1],sqrt(sumEventsBckTypeE[1]+sumEventsSigTypeE[1]),sumEventsBckType[2]+sumEventsSigType[2],sqrt(sumEventsBckTypeE[2]+sumEventsSigTypeE[2]));
    printf("--------------------------------------------------------------------------------\n");
  }
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    char output[200];
    sprintf(output,"histoww_nice%d_%d.root",nJetsType,thePlot);	  
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }
  printf("QCD Corr: qqWW(%f:%f/%f/%f/%f/%f/%f) ggWW(%f:%f/%f/%f/%f/%f/%f) Top(%f:%f/%f/%f/%f/%f/%f) DY(%f:%f/%f/%f/%f/%f/%f) VV(%f:%f/%f/%f/%f/%f/%f) VVV(%f:%f/%f/%f/%f/%f/%f) WG(%f:%f/%f/%f/%f/%f/%f) WGS(%f:%f/%f/%f/%f/%f/%f) Higgs(%f:%f/%f/%f/%f/%f/%f)\n",
    histo_qqWW->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_qqWW_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_ggWW->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_ggWW_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_Top->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_Top_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_DY->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_DY_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VV->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_VVV->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_VVV_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_WG->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WG_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_WGS->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_WGS_CMS_QCDScaleBounding[5]->GetSumOfWeights(),
    histo_Higgs->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[0]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[1]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[2]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[3]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[4]->GetSumOfWeights(),histo_Higgs_CMS_QCDScaleBounding[5]->GetSumOfWeights());

  printf("ptWW  raw: nlo(%f)up/down(%f/%f/%f/%f)nnlo(%f)\n",histo_qqWW_CMS_MVAWW_nlo->GetSumOfWeights(),
                                     histo_qqWW_CMS_MVAWW_qup_nlo->GetSumOfWeights(),histo_qqWW_CMS_MVAWW_qdown_nlo->GetSumOfWeights(),
                                     histo_qqWW_CMS_MVAWW_sup_nlo->GetSumOfWeights(),histo_qqWW_CMS_MVAWW_sdown_nlo->GetSumOfWeights(),
				     histo_qqWW->GetSumOfWeights());

  double ratio_nnlo_nlo = histo_qqWW->GetSumOfWeights()/histo_qqWW_CMS_MVAWW_nlo->GetSumOfWeights();
  histo_qqWW_CMS_MVAWW_nlo      ->Scale(ratio_nnlo_nlo);
  histo_qqWW_CMS_MVAWW_qup_nlo  ->Scale(ratio_nnlo_nlo);
  histo_qqWW_CMS_MVAWW_qdown_nlo->Scale(ratio_nnlo_nlo);
  histo_qqWW_CMS_MVAWW_sup_nlo  ->Scale(ratio_nnlo_nlo);
  histo_qqWW_CMS_MVAWW_sdown_nlo->Scale(ratio_nnlo_nlo);

  printf("ptWW norm: nlo(%f)up/down(%f/%f/%f/%f)nnlo(%f)\n",histo_qqWW_CMS_MVAWW_nlo->GetSumOfWeights(),
                                     histo_qqWW_CMS_MVAWW_qup_nlo->GetSumOfWeights(),histo_qqWW_CMS_MVAWW_qdown_nlo->GetSumOfWeights(),
                                     histo_qqWW_CMS_MVAWW_sup_nlo->GetSumOfWeights(),histo_qqWW_CMS_MVAWW_sdown_nlo->GetSumOfWeights(),
				     histo_qqWW->GetSumOfWeights());

  for(int i=1; i<=histo_qqWW->GetNbinsX(); i++){
    double factorUp = +1.0; double factorDown = -1.0;
    histo_qqWW_CMS_MVAqqWWStatBoundingUp       ->SetBinContent(i,TMath::Max(histo_qqWW    ->GetBinContent(i)+factorUp  *histo_qqWW    ->GetBinError(i),0.000001));
    histo_qqWW_CMS_MVAqqWWStatBoundingDown     ->SetBinContent(i,TMath::Max(histo_qqWW    ->GetBinContent(i)+factorDown*histo_qqWW    ->GetBinError(i),0.000001));
    histo_ggWW_CMS_MVAggWWStatBoundingUp       ->SetBinContent(i,TMath::Max(histo_ggWW    ->GetBinContent(i)+factorUp  *histo_ggWW    ->GetBinError(i),0.000001));
    histo_ggWW_CMS_MVAggWWStatBoundingDown     ->SetBinContent(i,TMath::Max(histo_ggWW    ->GetBinContent(i)+factorDown*histo_ggWW    ->GetBinError(i),0.000001));
    histo_Top_CMS_MVATopStatBoundingUp         ->SetBinContent(i,TMath::Max(histo_Top	 ->GetBinContent(i)+factorUp  *histo_Top    ->GetBinError(i),0.000001));
    histo_Top_CMS_MVATopStatBoundingDown       ->SetBinContent(i,TMath::Max(histo_Top	 ->GetBinContent(i)+factorDown*histo_Top    ->GetBinError(i),0.000001));
    histo_DY_CMS_MVADYStatBoundingUp           ->SetBinContent(i,TMath::Max(histo_DY    ->GetBinContent(i)+factorUp  *histo_DY    ->GetBinError(i),0.000001));
    histo_DY_CMS_MVADYStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_DY    ->GetBinContent(i)+factorDown*histo_DY    ->GetBinError(i),0.000001));
    histo_VV_CMS_MVAVVStatBoundingUp           ->SetBinContent(i,TMath::Max(histo_VV    ->GetBinContent(i)+factorUp  *histo_VV    ->GetBinError(i),0.000001));
    histo_VV_CMS_MVAVVStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_VV    ->GetBinContent(i)+factorDown*histo_VV    ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingUp         ->SetBinContent(i,TMath::Max(histo_VVV    ->GetBinContent(i)+factorUp  *histo_VVV    ->GetBinError(i),0.000001));
    histo_VVV_CMS_MVAVVVStatBoundingDown       ->SetBinContent(i,TMath::Max(histo_VVV    ->GetBinContent(i)+factorDown*histo_VVV    ->GetBinError(i),0.000001));
    histo_WG_CMS_MVAWGStatBoundingUp           ->SetBinContent(i,TMath::Max(histo_WG    ->GetBinContent(i)+factorUp  *histo_WG    ->GetBinError(i),0.000001));
    histo_WG_CMS_MVAWGStatBoundingDown         ->SetBinContent(i,TMath::Max(histo_WG    ->GetBinContent(i)+factorDown*histo_WG    ->GetBinError(i),0.000001));
    histo_WGS_CMS_MVAWGSStatBoundingUp         ->SetBinContent(i,TMath::Max(histo_WGS    ->GetBinContent(i)+factorUp  *histo_WGS    ->GetBinError(i),0.000001));
    histo_WGS_CMS_MVAWGSStatBoundingDown       ->SetBinContent(i,TMath::Max(histo_WGS    ->GetBinContent(i)+factorDown*histo_WGS    ->GetBinError(i),0.000001));
    histo_WjetsM_CMS_MVAWjetsMStatBoundingUp   ->SetBinContent(i,TMath::Max(histo_WjetsM    ->GetBinContent(i)+factorUp  *histo_WjetsM    ->GetBinError(i),0.000001));
    histo_WjetsM_CMS_MVAWjetsMStatBoundingDown ->SetBinContent(i,TMath::Max(histo_WjetsM    ->GetBinContent(i)+factorDown*histo_WjetsM    ->GetBinError(i),0.000001));
    histo_WjetsE_CMS_MVAWjetsEStatBoundingUp   ->SetBinContent(i,TMath::Max(histo_WjetsE    ->GetBinContent(i)+factorUp  *histo_WjetsE    ->GetBinError(i),0.000001));
    histo_WjetsE_CMS_MVAWjetsEStatBoundingDown ->SetBinContent(i,TMath::Max(histo_WjetsE    ->GetBinContent(i)+factorDown*histo_WjetsE    ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingUp     ->SetBinContent(i,TMath::Max(histo_Higgs    ->GetBinContent(i)+factorUp  *histo_Higgs    ->GetBinError(i),0.000001));
    histo_Higgs_CMS_MVAHiggsStatBoundingDown   ->SetBinContent(i,TMath::Max(histo_Higgs    ->GetBinContent(i)+factorDown*histo_Higgs    ->GetBinError(i),0.000001));
  }
  char outputLimits[200];
  sprintf(outputLimits,"ww%4s_%dj_input_%s.root",finalStateName,nJetsType,ECMsb.Data());
  TFile* outFileLimits = new TFile(outputLimits,"recreate");
  outFileLimits->cd();
  histo_Data ->Write();
  histo_qqWW ->Write();
  histo_ggWW ->Write();
  histo_Top  ->Write();
  histo_DY   ->Write();
  histo_VV   ->Write();
  histo_VVV  ->Write();
  histo_WG   ->Write();
  histo_WGS  ->Write();
  histo_WjetsM->Write();
  histo_WjetsE->Write();
  histo_Higgs->Write();
  cout << histo_Data  ->GetSumOfWeights() << " ";
  cout << histo_qqWW  ->GetSumOfWeights() << " ";
  cout << histo_ggWW  ->GetSumOfWeights() << " ";
  cout << histo_Top   ->GetSumOfWeights() << " ";
  cout << histo_DY    ->GetSumOfWeights() << " ";
  cout << histo_VV    ->GetSumOfWeights() << " ";
  cout << histo_VVV   ->GetSumOfWeights() << " ";
  cout << histo_WG    ->GetSumOfWeights() << " ";
  cout << histo_WGS   ->GetSumOfWeights() << " ";
  cout << histo_WjetsM->GetSumOfWeights() << " ";
  cout << histo_WjetsE->GetSumOfWeights() << " ";
  cout << histo_Higgs ->GetSumOfWeights() << " ";
  cout << endl;
  printf("uncertainties Stat\n");
  histo_qqWW_CMS_MVAqqWWStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW   ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAqqWWStatBoundingUp	   ->GetBinContent(i)/histo_qqWW  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAqqWWStatBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW   ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAqqWWStatBoundingDown    ->GetBinContent(i)/histo_qqWW  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAggWWStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW   ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAggWWStatBoundingUp	   ->GetBinContent(i)/histo_ggWW  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAggWWStatBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW   ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAggWWStatBoundingDown    ->GetBinContent(i)/histo_ggWW  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVATopStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVATopStatBoundingUp	->GetBinContent(i)/histo_Top  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVATopStatBoundingDown    ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVATopStatBoundingDown	->GetBinContent(i)/histo_Top  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVADYStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY     ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVADYStatBoundingUp      ->GetBinContent(i)/histo_DY  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVADYStatBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY     ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVADYStatBoundingDown    ->GetBinContent(i)/histo_DY  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAVVStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV     ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAVVStatBoundingUp      ->GetBinContent(i)/histo_VV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAVVStatBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV     ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAVVStatBoundingDown    ->GetBinContent(i)/histo_VV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingUp	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAVVVStatBoundingDown    ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAVVVStatBoundingDown	->GetBinContent(i)/histo_VVV  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAWGStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG     ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAWGStatBoundingUp      ->GetBinContent(i)/histo_WG  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAWGStatBoundingDown      ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG     ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAWGStatBoundingDown    ->GetBinContent(i)/histo_WG  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAWGSStatBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAWGSStatBoundingUp	->GetBinContent(i)/histo_WGS  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAWGSStatBoundingDown    ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAWGSStatBoundingDown	->GetBinContent(i)/histo_WGS  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsM_CMS_MVAWjetsMStatBoundingUp  ->Write();for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsM ->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsM_CMS_MVAWjetsMStatBoundingUp      ->GetBinContent(i)/histo_WjetsM  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsM_CMS_MVAWjetsMStatBoundingDown->Write();for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsM ->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsM_CMS_MVAWjetsMStatBoundingDown    ->GetBinContent(i)/histo_WjetsM  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsE_CMS_MVAWjetsEStatBoundingUp  ->Write();for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsE ->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsE_CMS_MVAWjetsEStatBoundingUp      ->GetBinContent(i)/histo_WjetsE  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WjetsE_CMS_MVAWjetsEStatBoundingDown->Write();for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WjetsE ->GetBinContent(i)>0)printf("%5.1f ",histo_WjetsE_CMS_MVAWjetsEStatBoundingDown    ->GetBinContent(i)/histo_WjetsE  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAHiggsStatBoundingUp  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAHiggsStatBoundingUp      ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAHiggsStatBoundingDown->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs  ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAHiggsStatBoundingDown    ->GetBinContent(i)/histo_Higgs  ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffM\n");
  histo_qqWW_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_qqWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_qqWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_ggWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_ggWW_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_Top_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_Top_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_DY_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_DY_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_WG_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_WG_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_WGS_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVALepEffMBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_WGS_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffMBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffMBoundingUp      ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffMBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffMBoundingDown    ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffMBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties LepEffE\n");
  histo_qqWW_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_qqWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_qqWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_ggWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_ggWW_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_Top_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_Top_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_DY_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_DY_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_VVV_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_WG_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_WG_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_WGS_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVALepEffEBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_WGS_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffEBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffEBoundingUp      ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVALepEffEBoundingDown  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVALepEffEBoundingDown    ->GetBinContent(i)/histo_Higgs_CMS_MVALepEffEBoundingAvg    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties MET\n");
  histo_qqWW_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAMETBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAMETBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAMETBoundingUp      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAMETBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAMETBoundingDown    ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JES\n");
  histo_qqWW_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVAJESBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAJESBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAJESBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAJESBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJESBoundingUp      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAJESBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJESBoundingDown    ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties JER\n");
  histo_qqWW_CMS_MVAJERBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAJERBoundingUp      ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAJERBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAJERBoundingDown    ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAJERBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAJERBoundingUp      ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVAJERBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVAJERBoundingDown    ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAJERBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAJERBoundingUp      ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVAJERBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVAJERBoundingDown    ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVAJERBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVAJERBoundingUp      ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVAJERBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVAJERBoundingDown    ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAJERBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAJERBoundingUp      ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVAJERBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVAJERBoundingDown    ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJERBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJERBoundingUp      ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVAJERBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVAJERBoundingDown    ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAJERBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAJERBoundingUp      ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVAJERBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVAJERBoundingDown    ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAJERBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAJERBoundingUp      ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVAJERBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVAJERBoundingDown    ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAJERBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJERBoundingUp      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVAJERBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVAJERBoundingDown    ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties BTAG\n");
  histo_qqWW_CMS_MVABTAGBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVABTAGBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVABTAGBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_Top    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVABTAGBoundingUp          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVABTAGBoundingUp          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVABTAGBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_VVV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVABTAGBoundingUp          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVABTAGBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_WGS    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVABTAGBoundingUp	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVABTAGBoundingUp      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_MVABTAGBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_MVABTAGBoundingDown    ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties PU\n");
  histo_qqWW_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_PUBoundingUp	   ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_PUBoundingDown    ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_PUBoundingUp	   ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_ggWW_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_ggWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_ggWW_CMS_PUBoundingDown    ->GetBinContent(i)/histo_ggWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_PUBoundingUp	 ->GetBinContent(i)/histo_Top	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Top_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Top    ->GetBinContent(i)>0)printf("%5.1f ",histo_Top_CMS_PUBoundingDown	 ->GetBinContent(i)/histo_Top	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_PUBoundingUp      ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_DY_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_DY    ->GetBinContent(i)>0)printf("%5.1f ",histo_DY_CMS_PUBoundingDown    ->GetBinContent(i)/histo_DY    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_PUBoundingUp      ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VV_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VV_CMS_PUBoundingDown    ->GetBinContent(i)/histo_VV    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingUp	 ->GetBinContent(i)/histo_VVV	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_VVV_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_VVV    ->GetBinContent(i)>0)printf("%5.1f ",histo_VVV_CMS_PUBoundingDown	 ->GetBinContent(i)/histo_VVV	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_PUBoundingUp      ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WG_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WG    ->GetBinContent(i)>0)printf("%5.1f ",histo_WG_CMS_PUBoundingDown    ->GetBinContent(i)/histo_WG    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_PUBoundingUp	 ->GetBinContent(i)/histo_WGS	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_WGS_CMS_PUBoundingDown	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_WGS    ->GetBinContent(i)>0)printf("%5.1f ",histo_WGS_CMS_PUBoundingDown	 ->GetBinContent(i)/histo_WGS	 ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_PUBoundingUp	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_PUBoundingUp      ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_Higgs_CMS_PUBoundingDown	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_Higgs    ->GetBinContent(i)>0)printf("%5.1f ",histo_Higgs_CMS_PUBoundingDown    ->GetBinContent(i)/histo_Higgs    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  printf("uncertainties WW NNLO\n");
  histo_qqWW_CMS_MVAWW_qup_nlo	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWW_qup_nlo         ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAWW_qdown_nlo	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWW_qdown_nlo       ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAWW_sup_nlo	          ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWW_sup_nlo         ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  histo_qqWW_CMS_MVAWW_sdown_nlo	  ->Write(); for(int i=1; i<=histo_qqWW->GetNbinsX(); i++) {if(histo_qqWW    ->GetBinContent(i)>0)printf("%5.1f ",histo_qqWW_CMS_MVAWW_sdown_nlo       ->GetBinContent(i)/histo_qqWW    ->GetBinContent(i)*100);else printf("100.0 ");} printf("\n");
  outFileLimits->Close();
  double theUEPS[3] = {1.032, 1.036, 1.036};
  if(nJetsType != 0){
    theUEPS[0] = 0.995; theUEPS[1] = 0.947; theUEPS[2] = 0.965;
  }
  double lumiE = 1.025;
  double systLepResE[9] = {1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01};
  double systLepResM[9] = {1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01,1.01};
  double syst_btag = 1.02;
  for(int nb=1; nb<=nBinMVA; nb++){
     // QCD study
    double systQCDScale[9] = {TMath::Abs(histo_qqWW_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_qqWW->GetBinContent(nb)),
                              TMath::Abs(histo_ggWW_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_ggWW->GetBinContent(nb)),
                              TMath::Abs(histo_Top_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_Top->GetBinContent(nb)),
                              TMath::Abs(histo_DY_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_DY->GetBinContent(nb)),
                              TMath::Abs(histo_VV_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_VV->GetBinContent(nb)),
                              TMath::Abs(histo_VVV_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_VVV->GetBinContent(nb)),
                              TMath::Abs(histo_WG_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_WG->GetBinContent(nb)),
                              TMath::Abs(histo_WGS_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_WGS->GetBinContent(nb)),
                              TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[0]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb))};
    for(int nqcd=1; nqcd<6; nqcd++) {
      if(TMath::Abs(histo_qqWW_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_qqWW->GetBinContent(nb)) >   systQCDScale[0]) systQCDScale[0] = TMath::Abs(histo_qqWW_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_qqWW->GetBinContent(nb));
      if(TMath::Abs(histo_ggWW_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_ggWW->GetBinContent(nb)) >   systQCDScale[1]) systQCDScale[1] = TMath::Abs(histo_ggWW_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_ggWW->GetBinContent(nb));
      if(TMath::Abs(histo_Top_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Top->GetBinContent(nb)) >     systQCDScale[2]) systQCDScale[2] = TMath::Abs(histo_Top_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Top->GetBinContent(nb));
      if(TMath::Abs(histo_DY_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_DY->GetBinContent(nb)) >       systQCDScale[3]) systQCDScale[3] = TMath::Abs(histo_DY_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_DY->GetBinContent(nb));
      if(TMath::Abs(histo_VV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VV->GetBinContent(nb)) >       systQCDScale[4]) systQCDScale[4] = TMath::Abs(histo_VV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VV->GetBinContent(nb));
      if(TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb)) >     systQCDScale[5]) systQCDScale[5] = TMath::Abs(histo_VVV_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_VVV->GetBinContent(nb));
      if(TMath::Abs(histo_WG_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_WG->GetBinContent(nb)) >       systQCDScale[6]) systQCDScale[6] = TMath::Abs(histo_WG_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_WG->GetBinContent(nb));
      if(TMath::Abs(histo_WGS_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_WGS->GetBinContent(nb)) >     systQCDScale[7]) systQCDScale[7] = TMath::Abs(histo_WGS_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_WGS->GetBinContent(nb));
      if(TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb)) > systQCDScale[8]) systQCDScale[8] = TMath::Abs(histo_Higgs_CMS_QCDScaleBounding[nqcd]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb));
    }
    if(histo_qqWW ->GetBinContent(nb) != 0) systQCDScale[0] = 1 + systQCDScale[0]/histo_qqWW ->GetBinContent(nb); else systQCDScale[0] = 1;
    if(histo_ggWW ->GetBinContent(nb) != 0) systQCDScale[1] = 1 + systQCDScale[1]/histo_ggWW ->GetBinContent(nb); else systQCDScale[1] = 1;
    if(histo_Top  ->GetBinContent(nb) != 0) systQCDScale[2] = 1 + systQCDScale[2]/histo_Top  ->GetBinContent(nb); else systQCDScale[2] = 1;
    if(histo_DY   ->GetBinContent(nb) != 0) systQCDScale[3] = 1 + systQCDScale[3]/histo_DY   ->GetBinContent(nb); else systQCDScale[3] = 1;
    if(histo_VV   ->GetBinContent(nb) != 0) systQCDScale[4] = 1 + systQCDScale[4]/histo_VV   ->GetBinContent(nb); else systQCDScale[4] = 1;
    if(histo_VVV  ->GetBinContent(nb) != 0) systQCDScale[5] = 1 + systQCDScale[5]/histo_VVV  ->GetBinContent(nb); else systQCDScale[5] = 1;
    if(histo_WG   ->GetBinContent(nb) != 0) systQCDScale[6] = 1 + systQCDScale[6]/histo_WG   ->GetBinContent(nb); else systQCDScale[6] = 1;
    if(histo_WGS  ->GetBinContent(nb) != 0) systQCDScale[7] = 1 + systQCDScale[7]/histo_WGS  ->GetBinContent(nb); else systQCDScale[7] = 1;
    if(histo_Higgs->GetBinContent(nb) != 0) systQCDScale[8] = 1 + systQCDScale[8]/histo_Higgs->GetBinContent(nb); else systQCDScale[8] = 1;
    printf("QCDScale(%d): %f %f %f %f %f %f %f %f %f\n",nb,systQCDScale[0],systQCDScale[1],systQCDScale[2],systQCDScale[3],systQCDScale[4],systQCDScale[5],systQCDScale[6],systQCDScale[7],systQCDScale[8]);
    
    // DY QCDScale change
    systQCDScale[3] = 1.10;

    // PDF study
    double systPDF[9];
    histo_Diff->Reset();
    if(histo_ggWW->GetBinContent(nb) > 0)for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_ggWW_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_ggWW->GetBinContent(nb))/histo_ggWW->GetBinContent(nb));
    systPDF[0] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_ggWW->GetBinContent(nb) > 0)for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_ggWW_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_ggWW->GetBinContent(nb))/histo_ggWW->GetBinContent(nb));
    systPDF[1] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_Top->GetBinContent(nb) > 0)for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_Top_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_Top->GetBinContent(nb))/histo_Top->GetBinContent(nb));
    systPDF[2] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_DY->GetBinContent(nb) > 0)for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_DY_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_DY->GetBinContent(nb))/histo_DY->GetBinContent(nb));
    systPDF[3] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_VV->GetBinContent(nb) > 0)for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_VV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VV->GetBinContent(nb))/histo_VV->GetBinContent(nb));
    systPDF[4] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_VVV->GetBinContent(nb) > 0)for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_VVV_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_VVV->GetBinContent(nb))/histo_VVV->GetBinContent(nb));
    systPDF[5] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_WG->GetBinContent(nb) > 0)for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_WG_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_WG->GetBinContent(nb))/histo_WG->GetBinContent(nb));
    systPDF[6] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_WGS->GetBinContent(nb) > 0)for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_WGS_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_WGS->GetBinContent(nb))/histo_WGS->GetBinContent(nb));
    systPDF[7] = 1.0+histo_Diff->GetRMS();

    histo_Diff->Reset();
    if(histo_Higgs->GetBinContent(nb) > 0)for(int npdf=0; npdf<100; npdf++) histo_Diff->Fill((histo_Higgs_CMS_PDFBounding[npdf]->GetBinContent(nb)-histo_Higgs->GetBinContent(nb))/histo_Higgs->GetBinContent(nb));
    systPDF[8] = 1.0+histo_Diff->GetRMS();

    printf("PDF(%d): %f %f %f %f %f %f %f %f %f\n",nb,systPDF[0],systPDF[1],systPDF[2],systPDF[3],systPDF[4],systPDF[5],systPDF[6],systPDF[7],systPDF[8]);

    double systLepEffM[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_qqWW_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[0] = histo_qqWW_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_qqWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_qqWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_qqWW_CMS_MVALepEffMBoundingDown   ->GetBinContent(nb) > 0) systLepEffM[0] = histo_qqWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_qqWW_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_ggWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_ggWW_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[1] = histo_ggWW_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_ggWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_ggWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_ggWW_CMS_MVALepEffMBoundingDown   ->GetBinContent(nb) > 0) systLepEffM[1] = histo_ggWW_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_ggWW_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_Top_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_Top_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[2] = histo_Top_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_Top_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_Top_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_Top_CMS_MVALepEffMBoundingDown     ->GetBinContent(nb) > 0) systLepEffM[2] = histo_Top_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_Top_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_DY_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_DY_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[3] = histo_DY_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_DY_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_DY_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_DY_CMS_MVALepEffMBoundingDown       ->GetBinContent(nb) > 0) systLepEffM[3] = histo_DY_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_DY_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_VV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VV_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[4] = histo_VV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_VV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VV_CMS_MVALepEffMBoundingDown       ->GetBinContent(nb) > 0) systLepEffM[4] = histo_VV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[5] = histo_VVV_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffMBoundingDown     ->GetBinContent(nb) > 0) systLepEffM[5] = histo_VVV_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_WG_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[6] = histo_WG_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_WG_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_WG_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffMBoundingDown       ->GetBinContent(nb) > 0) systLepEffM[6] = histo_WG_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_WG_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_WGS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_WGS_CMS_MVALepEffMBoundingUp	   ->GetBinContent(nb) > 0) systLepEffM[7] = histo_WGS_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_WGS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_WGS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_WGS_CMS_MVALepEffMBoundingDown     ->GetBinContent(nb) > 0) systLepEffM[7] = histo_WGS_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_WGS_CMS_MVALepEffMBoundingDown->GetBinContent(nb);
    if     (histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffMBoundingUp   ->GetBinContent(nb) > 0) systLepEffM[8] = histo_Higgs_CMS_MVALepEffMBoundingUp->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb);
    else if(histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffMBoundingDown ->GetBinContent(nb) > 0) systLepEffM[8] = histo_Higgs_CMS_MVALepEffMBoundingAvg->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffMBoundingDown->GetBinContent(nb);

    double systLepEffE[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_qqWW_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[0] = histo_qqWW_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_qqWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_qqWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_qqWW_CMS_MVALepEffEBoundingDown   ->GetBinContent(nb) > 0) systLepEffE[0] = histo_qqWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_qqWW_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_ggWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_ggWW_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[1] = histo_ggWW_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_ggWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_ggWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_ggWW_CMS_MVALepEffEBoundingDown   ->GetBinContent(nb) > 0) systLepEffE[1] = histo_ggWW_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_ggWW_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_Top_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_Top_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[2] = histo_Top_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_Top_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_Top_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_Top_CMS_MVALepEffEBoundingDown     ->GetBinContent(nb) > 0) systLepEffE[2] = histo_Top_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_Top_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_DY_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_DY_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[3] = histo_DY_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_DY_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_DY_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_DY_CMS_MVALepEffEBoundingDown       ->GetBinContent(nb) > 0) systLepEffE[3] = histo_DY_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_DY_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_VV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VV_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[4] = histo_VV_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_VV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_VV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VV_CMS_MVALepEffEBoundingDown       ->GetBinContent(nb) > 0) systLepEffE[4] = histo_VV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_VV_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[5] = histo_VVV_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_VVV_CMS_MVALepEffEBoundingDown     ->GetBinContent(nb) > 0) systLepEffE[5] = histo_VVV_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_VVV_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_WG_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[6] = histo_WG_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_WG_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_WG_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_WG_CMS_MVALepEffEBoundingDown       ->GetBinContent(nb) > 0) systLepEffE[6] = histo_WG_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_WG_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_WGS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_WGS_CMS_MVALepEffEBoundingUp	   ->GetBinContent(nb) > 0) systLepEffE[7] = histo_WGS_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_WGS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_WGS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_WGS_CMS_MVALepEffEBoundingDown     ->GetBinContent(nb) > 0) systLepEffE[7] = histo_WGS_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_WGS_CMS_MVALepEffEBoundingDown->GetBinContent(nb);
    if     (histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffEBoundingUp   ->GetBinContent(nb) > 0) systLepEffE[8] = histo_Higgs_CMS_MVALepEffEBoundingUp->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb);
    else if(histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb) > 0 && histo_Higgs_CMS_MVALepEffEBoundingDown ->GetBinContent(nb) > 0) systLepEffE[8] = histo_Higgs_CMS_MVALepEffEBoundingAvg->GetBinContent(nb)/histo_Higgs_CMS_MVALepEffEBoundingDown->GetBinContent(nb);

    double systMet[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[0] = histo_qqWW_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAMETBoundingDown   ->GetBinContent(nb) > 0) systMet[0] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[1] = histo_ggWW_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_ggWW->GetBinContent(nb);
    else if(histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVAMETBoundingDown   ->GetBinContent(nb) > 0) systMet[1] = histo_ggWW->GetBinContent(nb)/histo_ggWW_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[2] = histo_Top_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_Top->GetBinContent(nb);
    else if(histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVAMETBoundingDown     ->GetBinContent(nb) > 0) systMet[2] = histo_Top->GetBinContent(nb)/histo_Top_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[3] = histo_DY_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_DY->GetBinContent(nb);
    else if(histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVAMETBoundingDown       ->GetBinContent(nb) > 0) systMet[3] = histo_DY->GetBinContent(nb)/histo_DY_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[4] = histo_VV_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_VV->GetBinContent(nb);
    else if(histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVAMETBoundingDown       ->GetBinContent(nb) > 0) systMet[4] = histo_VV->GetBinContent(nb)/histo_VV_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[5] = histo_VVV_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAMETBoundingDown     ->GetBinContent(nb) > 0) systMet[5] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[6] = histo_WG_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_WG->GetBinContent(nb);
    else if(histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVAMETBoundingDown       ->GetBinContent(nb) > 0) systMet[6] = histo_WG->GetBinContent(nb)/histo_WG_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVAMETBoundingUp	    ->GetBinContent(nb) > 0) systMet[7] = histo_WGS_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_WGS->GetBinContent(nb);
    else if(histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVAMETBoundingDown     ->GetBinContent(nb) > 0) systMet[7] = histo_WGS->GetBinContent(nb)/histo_WGS_CMS_MVAMETBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAMETBoundingUp   ->GetBinContent(nb) > 0) systMet[8] = histo_Higgs_CMS_MVAMETBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAMETBoundingDown ->GetBinContent(nb) > 0) systMet[8] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_MVAMETBoundingDown->GetBinContent(nb);
    for(int i=0; i<9; i++) if(systMet[i] == 1) systMet[i] = 0.998;
    for(int nmet=0; nmet<9; nmet++) if(systMet[nmet] > 1.10) systMet[nmet] = 1.10;
    for(int nmet=0; nmet<9; nmet++) if(systMet[nmet] < 0.90) systMet[nmet] = 0.90;
 
    double systJes[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[0] = histo_qqWW_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAJESBoundingDown   ->GetBinContent(nb) > 0) systJes[0] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[1] = histo_ggWW_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_ggWW->GetBinContent(nb);
    else if(histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVAJESBoundingDown   ->GetBinContent(nb) > 0) systJes[1] = histo_ggWW->GetBinContent(nb)/histo_ggWW_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[2] = histo_Top_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_Top->GetBinContent(nb);
    else if(histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVAJESBoundingDown     ->GetBinContent(nb) > 0) systJes[2] = histo_Top->GetBinContent(nb)/histo_Top_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[3] = histo_DY_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_DY->GetBinContent(nb);
    else if(histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVAJESBoundingDown       ->GetBinContent(nb) > 0) systJes[3] = histo_DY->GetBinContent(nb)/histo_DY_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[4] = histo_VV_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_VV->GetBinContent(nb);
    else if(histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVAJESBoundingDown       ->GetBinContent(nb) > 0) systJes[4] = histo_VV->GetBinContent(nb)/histo_VV_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[5] = histo_VVV_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAJESBoundingDown     ->GetBinContent(nb) > 0) systJes[5] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[6] = histo_WG_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_WG->GetBinContent(nb);
    else if(histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVAJESBoundingDown       ->GetBinContent(nb) > 0) systJes[6] = histo_WG->GetBinContent(nb)/histo_WG_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVAJESBoundingUp	    ->GetBinContent(nb) > 0) systJes[7] = histo_WGS_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_WGS->GetBinContent(nb);
    else if(histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVAJESBoundingDown     ->GetBinContent(nb) > 0) systJes[7] = histo_WGS->GetBinContent(nb)/histo_WGS_CMS_MVAJESBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAJESBoundingUp   ->GetBinContent(nb) > 0) systJes[8] = histo_Higgs_CMS_MVAJESBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAJESBoundingDown ->GetBinContent(nb) > 0) systJes[8] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_MVAJESBoundingDown->GetBinContent(nb);
    for(int njes=0; njes<9; njes++) if(systJes[njes] > 1.10) systJes[njes] = 1.10;
    for(int njes=0; njes<9; njes++) if(systJes[njes] < 0.90) systJes[njes] = 0.90;
 
    double systJer[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAJERBoundingUp	    ->GetBinContent(nb) > 0) systJer[0] = histo_qqWW_CMS_MVAJERBoundingUp->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAJERBoundingDown   ->GetBinContent(nb) > 0) systJer[0] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_MVAJERBoundingDown->GetBinContent(nb);
    if     (histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVAJERBoundingUp	    ->GetBinContent(nb) > 0) systJer[1] = histo_ggWW_CMS_MVAJERBoundingUp->GetBinContent(nb)/histo_ggWW->GetBinContent(nb);
    else if(histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVAJERBoundingDown   ->GetBinContent(nb) > 0) systJer[1] = histo_ggWW->GetBinContent(nb)/histo_ggWW_CMS_MVAJERBoundingDown->GetBinContent(nb);
    if     (histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVAJERBoundingUp	    ->GetBinContent(nb) > 0) systJer[2] = histo_Top_CMS_MVAJERBoundingUp->GetBinContent(nb)/histo_Top->GetBinContent(nb);
    else if(histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVAJERBoundingDown     ->GetBinContent(nb) > 0) systJer[2] = histo_Top->GetBinContent(nb)/histo_Top_CMS_MVAJERBoundingDown->GetBinContent(nb);
    if     (histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVAJERBoundingUp	    ->GetBinContent(nb) > 0) systJer[3] = histo_DY_CMS_MVAJERBoundingUp->GetBinContent(nb)/histo_DY->GetBinContent(nb);
    else if(histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVAJERBoundingDown       ->GetBinContent(nb) > 0) systJer[3] = histo_DY->GetBinContent(nb)/histo_DY_CMS_MVAJERBoundingDown->GetBinContent(nb);
    if     (histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVAJERBoundingUp	    ->GetBinContent(nb) > 0) systJer[4] = histo_VV_CMS_MVAJERBoundingUp->GetBinContent(nb)/histo_VV->GetBinContent(nb);
    else if(histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVAJERBoundingDown       ->GetBinContent(nb) > 0) systJer[4] = histo_VV->GetBinContent(nb)/histo_VV_CMS_MVAJERBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAJERBoundingUp	    ->GetBinContent(nb) > 0) systJer[5] = histo_VVV_CMS_MVAJERBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVAJERBoundingDown     ->GetBinContent(nb) > 0) systJer[5] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVAJERBoundingDown->GetBinContent(nb);
    if     (histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVAJERBoundingUp	    ->GetBinContent(nb) > 0) systJer[6] = histo_WG_CMS_MVAJERBoundingUp->GetBinContent(nb)/histo_WG->GetBinContent(nb);
    else if(histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVAJERBoundingDown       ->GetBinContent(nb) > 0) systJer[6] = histo_WG->GetBinContent(nb)/histo_WG_CMS_MVAJERBoundingDown->GetBinContent(nb);
    if     (histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVAJERBoundingUp	    ->GetBinContent(nb) > 0) systJer[7] = histo_WGS_CMS_MVAJERBoundingUp->GetBinContent(nb)/histo_WGS->GetBinContent(nb);
    else if(histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVAJERBoundingDown     ->GetBinContent(nb) > 0) systJer[7] = histo_WGS->GetBinContent(nb)/histo_WGS_CMS_MVAJERBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAJERBoundingUp   ->GetBinContent(nb) > 0) systJer[8] = histo_Higgs_CMS_MVAJERBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVAJERBoundingDown ->GetBinContent(nb) > 0) systJer[8] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_MVAJERBoundingDown->GetBinContent(nb);
    for(int njer=0; njer<9; njer++) if(systJer[njer] > 1.10) systJer[njer] = 1.10;
    for(int njer=0; njer<9; njer++) if(systJer[njer] < 0.90) systJer[njer] = 0.90;

    double systBtag[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVABTAGBoundingUp    ->GetBinContent(nb) > 0) systBtag[0] = histo_qqWW_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVABTAGBoundingDown  ->GetBinContent(nb) > 0) systBtag[0] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVABTAGBoundingUp    ->GetBinContent(nb) > 0) systBtag[1] = histo_ggWW_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_ggWW->GetBinContent(nb);
    else if(histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_MVABTAGBoundingDown  ->GetBinContent(nb) > 0) systBtag[1] = histo_ggWW->GetBinContent(nb)/histo_ggWW_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVABTAGBoundingUp	    ->GetBinContent(nb) > 0) systBtag[2] = histo_Top_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_Top->GetBinContent(nb);
    else if(histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_MVABTAGBoundingDown    ->GetBinContent(nb) > 0) systBtag[2] = histo_Top->GetBinContent(nb)/histo_Top_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVABTAGBoundingUp	    ->GetBinContent(nb) > 0) systBtag[3] = histo_DY_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_DY->GetBinContent(nb);
    else if(histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_MVABTAGBoundingDown      ->GetBinContent(nb) > 0) systBtag[3] = histo_DY->GetBinContent(nb)/histo_DY_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVABTAGBoundingUp	    ->GetBinContent(nb) > 0) systBtag[4] = histo_VV_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_VV->GetBinContent(nb);
    else if(histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_MVABTAGBoundingDown      ->GetBinContent(nb) > 0) systBtag[4] = histo_VV->GetBinContent(nb)/histo_VV_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVABTAGBoundingUp	    ->GetBinContent(nb) > 0) systBtag[5] = histo_VVV_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_MVABTAGBoundingDown    ->GetBinContent(nb) > 0) systBtag[5] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVABTAGBoundingUp	    ->GetBinContent(nb) > 0) systBtag[6] = histo_WG_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_WG->GetBinContent(nb);
    else if(histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_MVABTAGBoundingDown      ->GetBinContent(nb) > 0) systBtag[6] = histo_WG->GetBinContent(nb)/histo_WG_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVABTAGBoundingUp	    ->GetBinContent(nb) > 0) systBtag[7] = histo_WGS_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_WGS->GetBinContent(nb);
    else if(histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_MVABTAGBoundingDown    ->GetBinContent(nb) > 0) systBtag[7] = histo_WGS->GetBinContent(nb)/histo_WGS_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVABTAGBoundingUp  ->GetBinContent(nb) > 0) systBtag[8] = histo_Higgs_CMS_MVABTAGBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_MVABTAGBoundingDown->GetBinContent(nb) > 0) systBtag[8] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_MVABTAGBoundingDown->GetBinContent(nb);
    for(int nbtag=0; nbtag<9; nbtag++) if(systBtag[nbtag] > 1.10) systBtag[nbtag] = 1.10;
    for(int nbtag=0; nbtag<9; nbtag++) if(systBtag[nbtag] < 0.90) systBtag[nbtag] = 0.90;
 
    double systPU[9] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_PUBoundingUp	    ->GetBinContent(nb) > 0) systPU[0] = histo_qqWW_CMS_PUBoundingUp->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_PUBoundingDown       ->GetBinContent(nb) > 0) systPU[0] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_PUBoundingUp	    ->GetBinContent(nb) > 0) systPU[1] = histo_ggWW_CMS_PUBoundingUp->GetBinContent(nb)/histo_ggWW->GetBinContent(nb);
    else if(histo_ggWW->GetBinContent(nb)> 0 && histo_ggWW_CMS_PUBoundingDown       ->GetBinContent(nb) > 0) systPU[1] = histo_ggWW->GetBinContent(nb)/histo_ggWW_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_PUBoundingUp	    ->GetBinContent(nb) > 0) systPU[2] = histo_Top_CMS_PUBoundingUp->GetBinContent(nb)/histo_Top->GetBinContent(nb);
    else if(histo_Top->GetBinContent(nb)> 0 && histo_Top_CMS_PUBoundingDown         ->GetBinContent(nb) > 0) systPU[2] = histo_Top->GetBinContent(nb)/histo_Top_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_PUBoundingUp	            ->GetBinContent(nb) > 0) systPU[3] = histo_DY_CMS_PUBoundingUp->GetBinContent(nb)/histo_DY->GetBinContent(nb);
    else if(histo_DY->GetBinContent(nb)> 0 && histo_DY_CMS_PUBoundingDown           ->GetBinContent(nb) > 0) systPU[3] = histo_DY->GetBinContent(nb)/histo_DY_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_PUBoundingUp	            ->GetBinContent(nb) > 0) systPU[4] = histo_VV_CMS_PUBoundingUp->GetBinContent(nb)/histo_VV->GetBinContent(nb);
    else if(histo_VV->GetBinContent(nb)> 0 && histo_VV_CMS_PUBoundingDown           ->GetBinContent(nb) > 0) systPU[4] = histo_VV->GetBinContent(nb)/histo_VV_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_PUBoundingUp	    ->GetBinContent(nb) > 0) systPU[5] = histo_VVV_CMS_PUBoundingUp->GetBinContent(nb)/histo_VVV->GetBinContent(nb);
    else if(histo_VVV->GetBinContent(nb)> 0 && histo_VVV_CMS_PUBoundingDown         ->GetBinContent(nb) > 0) systPU[5] = histo_VVV->GetBinContent(nb)/histo_VVV_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_PUBoundingUp	            ->GetBinContent(nb) > 0) systPU[6] = histo_WG_CMS_PUBoundingUp->GetBinContent(nb)/histo_WG->GetBinContent(nb);
    else if(histo_WG->GetBinContent(nb)> 0 && histo_WG_CMS_PUBoundingDown           ->GetBinContent(nb) > 0) systPU[6] = histo_WG->GetBinContent(nb)/histo_WG_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_PUBoundingUp	    ->GetBinContent(nb) > 0) systPU[7] = histo_WGS_CMS_PUBoundingUp->GetBinContent(nb)/histo_WGS->GetBinContent(nb);
    else if(histo_WGS->GetBinContent(nb)> 0 && histo_WGS_CMS_PUBoundingDown         ->GetBinContent(nb) > 0) systPU[7] = histo_WGS->GetBinContent(nb)/histo_WGS_CMS_PUBoundingDown->GetBinContent(nb);
    if     (histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_PUBoundingUp       ->GetBinContent(nb) > 0) systPU[8] = histo_Higgs_CMS_PUBoundingUp->GetBinContent(nb)/histo_Higgs->GetBinContent(nb);
    else if(histo_Higgs->GetBinContent(nb)> 0 && histo_Higgs_CMS_PUBoundingDown     ->GetBinContent(nb) > 0) systPU[8] = histo_Higgs->GetBinContent(nb)/histo_Higgs_CMS_PUBoundingDown->GetBinContent(nb);
    for(int npu=0; npu<9; npu++) if(systPU[npu] > 1.02) systPU[npu] = 1.02;
    for(int npu=0; npu<9; npu++) if(systPU[npu] < 0.98) systPU[npu] = 0.98;

    double systWWNNLO[2] = {1.0,1.0};
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAWW_qup_nlo	    ->GetBinContent(nb) > 0) systWWNNLO[0] = histo_qqWW_CMS_MVAWW_qup_nlo->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAWW_qdown_nlo      ->GetBinContent(nb) > 0) systWWNNLO[0] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_MVAWW_qdown_nlo->GetBinContent(nb);
    if     (histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAWW_sup_nlo	    ->GetBinContent(nb) > 0) systWWNNLO[1] = histo_qqWW_CMS_MVAWW_sup_nlo->GetBinContent(nb)/histo_qqWW->GetBinContent(nb);
    else if(histo_qqWW->GetBinContent(nb)> 0 && histo_qqWW_CMS_MVAWW_sdown_nlo      ->GetBinContent(nb) > 0) systWWNNLO[1] = histo_qqWW->GetBinContent(nb)/histo_qqWW_CMS_MVAWW_sdown_nlo->GetBinContent(nb);

    char outputLimits[200];
    sprintf(outputLimits,"ww_%s_%dj_input_%s_bin%d.root",finalStateName,nJetsType,ECMsb.Data(),nb-1);
    TFile* outFileLimits = new TFile(outputLimits,"recreate");
    outFileLimits->cd();
    
    histoOneBin_Data  ->Reset();
    histoOneBin_qqWW  ->Reset();
    histoOneBin_ggWW  ->Reset();
    histoOneBin_Top   ->Reset();
    histoOneBin_DY    ->Reset();
    histoOneBin_VV    ->Reset();
    histoOneBin_VVV   ->Reset();
    histoOneBin_WG    ->Reset();
    histoOneBin_WGS   ->Reset();
    histoOneBin_WjetsM->Reset();
    histoOneBin_WjetsE->Reset();
    histoOneBin_Higgs ->Reset();
    histoOneBin_Data  ->SetBinContent(1,	    histo_Data  ->GetBinContent(nb));
    histoOneBin_qqWW  ->SetBinContent(1, TMath::Max(histo_qqWW  ->GetBinContent(nb),0.0));
    histoOneBin_ggWW  ->SetBinContent(1, TMath::Max(histo_ggWW  ->GetBinContent(nb),0.0));
    histoOneBin_Top   ->SetBinContent(1, TMath::Max(histo_Top	->GetBinContent(nb),0.0));
    histoOneBin_DY    ->SetBinContent(1, TMath::Max(histo_DY	->GetBinContent(nb),0.0));
    histoOneBin_VV    ->SetBinContent(1, TMath::Max(histo_VV	->GetBinContent(nb),0.0));
    histoOneBin_VVV   ->SetBinContent(1, TMath::Max(histo_VVV	->GetBinContent(nb),0.0));
    histoOneBin_WG    ->SetBinContent(1, TMath::Max(histo_WG	->GetBinContent(nb),0.0));
    histoOneBin_WGS   ->SetBinContent(1, TMath::Max(histo_WGS	->GetBinContent(nb),0.0));
    histoOneBin_WjetsM->SetBinContent(1, TMath::Max(histo_WjetsM->GetBinContent(nb),0.0));
    histoOneBin_WjetsE->SetBinContent(1, TMath::Max(histo_WjetsE->GetBinContent(nb),0.0));
    histoOneBin_Higgs ->SetBinContent(1, TMath::Max(histo_Higgs ->GetBinContent(nb),0.0));
    histoOneBin_Data  ->Write();
    histoOneBin_qqWW  ->Write();
    histoOneBin_ggWW  ->Write();
    histoOneBin_Top   ->Write();
    histoOneBin_DY    ->Write();
    histoOneBin_VV    ->Write();
    histoOneBin_VVV   ->Write();
    histoOneBin_WG    ->Write();
    histoOneBin_WGS   ->Write();
    histoOneBin_WjetsM->Write();
    histoOneBin_WjetsE->Write();
    histoOneBin_Higgs ->Write();

    outFileLimits->Close();

    char outputLimitsShape[200];
    sprintf(outputLimitsShape,"histo_limits_ww%s_%dj_%s_bin%d.txt",finalStateName,nJetsType,ECMsb.Data(),nb-1);
    ofstream newcardShape;
    newcardShape.open(outputLimitsShape);
    newcardShape << Form("imax 1 number of channels\n");
    newcardShape << Form("jmax * number of background\n");
    newcardShape << Form("kmax * number of nuisance parameters\n");
    newcardShape << Form("shapes *   *   %s  histoOneBin_$PROCESS\n",outputLimits);
    newcardShape << Form("shapes data_obs * %s  histoOneBin_Data \n",outputLimits);
    newcardShape << Form("Observation %d\n",(int)histo_Data->GetBinContent(nb));
    newcardShape << Form("bin ww%2s%dj%d ww%2s%dj%d ww%2s%dj%d ww%2s%dj%d ww%2s%dj%d ww%2s%dj%d ww%2s%dj%d ww%2s%dj%d ww%2s%dj%d ww%2s%dj%d ww%2s%dj%d\n",finalStateName,nJetsType,nb-1,finalStateName,nJetsType,nb-1,finalStateName,nJetsType,nb-1,finalStateName,nJetsType,nb-1,finalStateName,nJetsType,nb-1,finalStateName,nJetsType,nb-1,finalStateName,nJetsType,nb-1,finalStateName,nJetsType,nb-1,finalStateName,nJetsType,nb-1,finalStateName,nJetsType,nb-1,finalStateName,nJetsType,nb-1);
    //                            0    1    2   3  4  5   6  7   8     9      10
    newcardShape << Form("process qqWW ggWW Top DY VV VVV WG WGS Higgs WjetsM WjetsE\n");
    newcardShape << Form("process -1 0 1 2 3 4 5 6 7 8 9\n");
    newcardShape << Form("rate %8.5f %8.5f  %8.5f  %8.5f %8.5f %8.5f %8.5f  %8.5f  %8.5f %8.5f %8.5f\n",TMath::Max(histo_qqWW->GetBinContent(nb),0.0),TMath::Max(histo_ggWW->GetBinContent(nb),0.0),TMath::Max(histo_Top->GetBinContent(nb),0.0),TMath::Max(histo_DY->GetBinContent(nb),0.0),TMath::Max(histo_VV->GetBinContent(nb),0.0),TMath::Max(histo_VVV->GetBinContent(nb),0.0),TMath::Max(histo_WG->GetBinContent(nb),0.0),TMath::Max(histo_WGS->GetBinContent(nb),0.0),TMath::Max(histo_Higgs->GetBinContent(nb),0.0),TMath::Max(histo_WjetsM->GetBinContent(nb),0.0),TMath::Max(histo_WjetsE->GetBinContent(nb),0.0));
    newcardShape << Form("lumi_%4s                               lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",ECMsb.Data(),lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE,lumiE);		     
    newcardShape << Form("%s                                     lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",effMName,systLepEffM[0],systLepEffM[1],systLepEffM[2],systLepEffM[3],systLepEffM[4],systLepEffM[5],systLepEffM[6],systLepEffM[7],systLepEffM[8]);
    newcardShape << Form("%s                                     lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",effEName,systLepEffE[0],systLepEffE[1],systLepEffE[2],systLepEffE[3],systLepEffE[4],systLepEffE[5],systLepEffE[6],systLepEffE[7],systLepEffE[8]);
    newcardShape << Form("%s                                     lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",momMName,systLepResM[0],systLepResM[1],systLepResM[2],systLepResM[3],systLepResM[4],systLepResM[5],systLepResM[6],systLepResM[7],systLepResM[8]);
    newcardShape << Form("%s                                     lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",momEName,systLepResE[0],systLepResE[1],systLepResE[2],systLepResE[3],systLepResE[4],systLepResE[5],systLepResE[6],systLepResE[7],systLepResE[8]);
    //newcardShape << Form("CMS_pu                                 lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",systPU[0],systPU[1],systPU[2],systPU[3],systPU[4],systPU[5],systPU[6],systPU[7],systPU[8]);
    newcardShape << Form("CMS_scale_met                          lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",systMet[0],systMet[1],systMet[2],systMet[3],systMet[4],systMet[5],systMet[6],systMet[7],systMet[8]);
    newcardShape << Form("CMS_scale_j                            lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",systJes[0],systJes[1],systJes[2],systJes[3],systJes[4],systJes[5],systJes[6],systJes[7],systJes[8]);
    newcardShape << Form("CMS_jer                                lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",systJer[0],systJer[1],systJer[2],systJer[3],systJer[4],systJer[5],systJer[6],systJer[7],systJer[8]);
    newcardShape << Form("CMS_eff_b_b2016                        lnN  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f  -    -  \n",systBtag[0],systBtag[1],systBtag[2],systBtag[3],systBtag[4],systBtag[5],systBtag[6],systBtag[7],systBtag[8]);
    newcardShape << Form("pdf_qqbar                              lnN  %7.5f  -      -  %7.5f %7.5f %7.5f %7.5f %7.5f	-    -    -  \n",TMath::Max(systPDF[0],1.01),TMath::Max(systPDF[3],1.01),TMath::Max(systPDF[4],1.01),TMath::Max(systPDF[5],1.01),TMath::Max(systPDF[6],1.01),1.01);
    newcardShape << Form("pdf_gg                                 lnN    -   %7.5f %7.5f   -	-     -     -	  -   %7.5f  -    -  \n",1.01,TMath::Max(systPDF[2],1.01),TMath::Max(systPDF[8],1.01));
    newcardShape << Form("QCDscale_VVV		                 lnN    -     -     -     -     -   %7.5f   -	  -	-    -    -  \n",systQCDScale[5]);	    
    newcardShape << Form("QCDscale_ggVV		                 lnN    -   %7.5f   -     -     -     -     -     -     -    -    -  \n",1.15);  
    newcardShape << Form("QCDscale_qqVV		                 lnN  %7.5f   -     -     -   %7.5f   -   %7.5f	%7.5f	-    -    -  \n",systQCDScale[0],systQCDScale[4],systQCDScale[6],systQCDScale[7]);	    
    newcardShape << Form("QCDscale_ggH		                 lnN    -     -     -     -     -     -     -     -   %7.5f  -    -  \n",systQCDScale[8]);	    
    if(useTopCR == true) {
    newcardShape << Form("QCDscale_Top		                 lnN    -     -   %7.5f   -     -     -     -     -     -    -    -  \n",systQCDScale[2]);	    
    newcardShape << Form("CMS_ww_Topnorm_jet%d rateParam  * Top 1 [0.1,10]\n",nJetsType);         
    }
    else {
    newcardShape << Form("norm_Top_%dj		                 lnN    -     -   %7.5f   -     -     -     -     -     -    -    -  \n",nJetsType,1.0 + topNormE[TMath::Min((int)nJetsType,2)]/topNorm[TMath::Min((int)nJetsType,2)]);	    
    }
    newcardShape << Form("norm_DY_%dj		                 lnN    -     -     -   %7.5f   -     -     -     -     -    -    -  \n",nJetsType,systQCDScale[3]);	    
    newcardShape << Form("norm_WGS		                 lnN    -     -     -     -     -     -     -   %7.5f   -    -    -  \n",1.30);	    
    newcardShape << Form("norm_WjetsM		                 lnN    -     -     -     -     -     -     -     -     -   %7.5f -  \n",1.30);	    
    newcardShape << Form("norm_WjetsE		                 lnN    -     -     -     -     -     -     -     -     -    -  %7.5f\n",1.30);	    
    newcardShape << Form("WWNNLO_resum		                 lnN  %7.5f   -     -     -     -     -     -     -	-    -    -  \n",systWWNNLO[0]);	    
    newcardShape << Form("WWNNLO_scale		                 lnN  %7.5f   -     -     -     -     -     -     -	-    -    -  \n",systWWNNLO[1]);	    
    newcardShape << Form("UEPS                                   lnN  %7.5f %7.5f   -     -     -     -     -     -   %7.5f  -    -  \n",theUEPS[0],theUEPS[1],theUEPS[2]);	    
    if(histo_qqWW->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAqqWWStatBounding_%s_Bin%d    lnN  %7.5f   -     -     -     -     -     -     -     -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_qqWW->GetBinError(nb)/histo_qqWW->GetBinContent(nb),0.999));
    if(histo_ggWW->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAggWWStatBounding_%s_Bin%d    lnN    -   %7.5f   -     -     -     -     -     -	 -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_ggWW->GetBinError(nb)/histo_ggWW->GetBinContent(nb),0.999));
    if(histo_Top->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVATopStatBounding_%s_Bin%d     lnN    -     -   %7.5f   -     -     -     -     -	 -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_Top->GetBinError(nb)/histo_Top->GetBinContent(nb),0.999));
    if(histo_DY->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVADYStatBounding_%s_Bin%d      lnN    -     -     -   %7.5f   -     -     -     -	 -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_DY->GetBinError(nb)/histo_DY->GetBinContent(nb),0.999));
    if(histo_VV->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAVVStatBounding_%s_Bin%d      lnN    -     -     -     -   %7.5f   -     -     -	 -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_VV->GetBinError(nb)/histo_VV->GetBinContent(nb),0.999));
    if(histo_VVV->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAVVVStatBounding_%s_Bin%d     lnN    -     -     -     -     -   %7.5f   -     -     -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_VVV->GetBinError(nb)/histo_VVV->GetBinContent(nb),0.999));
    if(histo_WG->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAWGStatBounding_%s_Bin%d      lnN    -     -     -     -     -     -   %7.5f   -     -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_WG->GetBinError(nb)/histo_WG->GetBinContent(nb),0.999));
    if(histo_WGS->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAWGSStatBounding_%s_Bin%d     lnN    -     -     -     -     -     -     -   %7.5f   -     -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_WGS->GetBinError(nb)/histo_WGS->GetBinContent(nb),0.999));
    if(histo_Higgs->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAHiggsStatBounding_%s_Bin%d   lnN    -     -     -     -     -     -     -     -   %7.5f   -   -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_Higgs->GetBinError(nb)/histo_Higgs->GetBinContent(nb),0.999));
    if(histo_WjetsM ->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAWjetsMStatBounding_%s_Bin%d  lnN    -     -     -     -     -     -     -     -	 -   %7.5f -  \n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_WjetsM ->GetBinError(nb)/histo_WjetsM ->GetBinContent(nb),0.999));
    if(histo_WjetsE ->GetBinContent(nb)	  > 0) newcardShape << Form("CMS_ww%s_%dj_MVAWjetsEStatBounding_%s_Bin%d  lnN    -     -     -     -     -     -     -     -	 -     - %7.5f\n",finalStateName,nJetsType,ECMsb.Data(),nb-1,1.0+TMath::Min(histo_WjetsE ->GetBinError(nb)/histo_WjetsE ->GetBinContent(nb),0.999));
    newcardShape.close();
  }
  for(int iF=0; iF<3; iF++){
    printf("double jetEpsBtagLOOSE[%d] = \n",iF);
    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
        printf("%5.4f",numBTaggingLOOSE[iEta][iPt][iF]/denBTagging[iEta][iPt][iF]);
        if(iPt!=4||iEta!=4) printf(",");
        if(iPt==4) printf("\n");
      }
    }
  }
  delete btagCalib;

  // delete garbage
  if(!isShapeAna || theControlRegion != 0){
    system(Form("rm -f ww_%s_%dj_input_%s_bin[2-%d].root",finalStateName,nJetsType,ECMsb.Data(),nBinMVA));
    system(Form("rm -f histo_limits_ww%s_%dj_%s_bin[2-%d].txt",finalStateName,nJetsType,ECMsb.Data(),nBinMVA));
  }
}
