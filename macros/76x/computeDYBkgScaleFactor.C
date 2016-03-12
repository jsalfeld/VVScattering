#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"
#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"

#include "MitAnalysisRunII/macros/76x/factors.h"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"

const TString typeLepSel = "default";
const bool useDYMVA = true;
double mcPrescale = 1.0;

// compute systematic uncertainty
Double_t computeSyst(const TH1D *hout, const TH1D *hin, Int_t binUsed, Double_t rErrorMax, Bool_t useRFromData = false, Bool_t isDebug = false);

void computeDYBkgScaleFactor(Int_t period = 1, Double_t MassZCut = 15){

  TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_76x/met_";
  TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_76x/met_";
  Double_t lumi = 2.318;

  float dymva_= -999.;
  unsigned int nlep_= -1;
  unsigned int njets_= -1;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  if(period==1){
  puPath = "MitAnalysisRunII/data/76x/puWeights_76x.root";
  infilenamev.push_back(Form("%sdata_AOD_Run2015C_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(1);

  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  				   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	 	   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(2);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(2);
  }
  else {assert(0);}

  if(infilenamev.size() != infilecatv.size()) assert(0);

  const Double_t mZ = 91.1876;
  
  const Int_t nbins = 4;  
  Float_t bins[nbins+1];

  const Int_t nmass = 3;
  const Double_t mH[nmass]     = {0, 125, 200};  
  bool useRFromData[3] = {1, 1, 1};

  //*******************************************************
  //Yields and  histograms
  //*******************************************************
  vector<Double_t> nin_kee_data, nin_kmm_data;
  
  vector<vector<TH1D*> > hNin_ree_mc,   hNout_ree_mc,   hNin_rmm_mc,   hNout_rmm_mc;
  vector<vector<TH1D*> > hNin_ree_da,   hNout_ree_da,   hNin_rmm_da,   hNout_rmm_da;

  vector<vector<Double_t> > nin_ee_dy, nout_ee_dy, nin_ee_vz, nout_ee_vz, nin_ee_data;  
  vector<vector<Double_t> > nin_mm_dy, nout_mm_dy, nin_mm_vz, nout_mm_vz, nin_mm_data;
  vector<vector<Double_t> > varin_ee_dy, varout_ee_dy, varin_ee_vz, varout_ee_vz;
  vector<vector<Double_t> > varin_mm_dy, varout_mm_dy, varin_mm_vz, varout_mm_vz;
  vector<vector<Double_t> > nin_em_data;

  
  for(UInt_t jetIndex = 0; jetIndex < 3; ++jetIndex) {
    Double_t tmp_nin_kee_data=0, tmp_nin_kmm_data=0;
    vector<TH1D*> tmp_hNin_ree_mc,   tmp_hNout_ree_mc,   tmp_hNin_rmm_mc,   tmp_hNout_rmm_mc;
    vector<TH1D*> tmp_hNin_ree_da,   tmp_hNout_ree_da,   tmp_hNin_rmm_da,   tmp_hNout_rmm_da;
    
    vector<Double_t> tmp_nin_ee_dy, tmp_nout_ee_dy, tmp_nin_ee_vz, tmp_nout_ee_vz, tmp_nin_ee_data;  
    vector<Double_t> tmp_nin_mm_dy, tmp_nout_mm_dy, tmp_nin_mm_vz, tmp_nout_mm_vz, tmp_nin_mm_data;
    vector<Double_t> tmp_varin_ee_dy, tmp_varout_ee_dy, tmp_varin_ee_vz, tmp_varout_ee_vz;
    vector<Double_t> tmp_varin_mm_dy, tmp_varout_mm_dy, tmp_varin_mm_vz, tmp_varout_mm_vz;
    vector<Double_t> tmp_nin_em_data;

    char hname[50];
    for(Int_t imass=0; imass<nmass; imass++) {
      if     (useDYMVA == kTRUE){
        if     (jetIndex == 0){
	  bins[0] = -0.20; bins[1] =  0.10; bins[2] = +0.20; bins[3] =  0.30; bins[4] = 1.0; 
	}
	else if(jetIndex == 1){
	  bins[0] = -0.20; bins[1] = +0.10; bins[2] = +0.20; bins[3] =  0.30; bins[4] = 1.0; 
	}
	else if(jetIndex == 2){
	  bins[0] = -0.20; bins[1] = +0.10; bins[2] = +0.20; bins[3] =  0.30; bins[4] = 1.0; 
	}
      } else {
	bins[0] = 20; bins[1] = 25; bins[2] = 35; bins[3] =  45; bins[4] = 70; 
      }

      sprintf(hname,"hNin_%iJet_ree_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_ree_mc.push_back(new TH1D(hname,"",nbins,bins));  tmp_hNin_ree_mc[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_ree_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_ree_mc.push_back(new TH1D(hname,"",nbins,bins)); tmp_hNout_ree_mc[imass]->Sumw2();
      sprintf(hname,"hNin_%iJet_rmm_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_rmm_mc.push_back(new TH1D(hname,"",nbins,bins));  tmp_hNin_rmm_mc[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_rmm_mc_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_rmm_mc.push_back(new TH1D(hname,"",nbins,bins)); tmp_hNout_rmm_mc[imass]->Sumw2();
      
      sprintf(hname,"hNin_%iJet_ree_da_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_ree_da.push_back(new TH1D(hname,"",nbins,bins));  tmp_hNin_ree_da[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_ree_da_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_ree_da.push_back(new TH1D(hname,"",nbins,bins)); tmp_hNout_ree_da[imass]->Sumw2();
      sprintf(hname,"hNin_%iJet_rmm_da_m%i",Int_t(jetIndex),(Int_t)mH[imass]);  tmp_hNin_rmm_da.push_back(new TH1D(hname,"",nbins,bins));  tmp_hNin_rmm_da[imass]->Sumw2();
      sprintf(hname,"hNout_%iJet_rmm_da_m%i",Int_t(jetIndex),(Int_t)mH[imass]); tmp_hNout_rmm_da.push_back(new TH1D(hname,"",nbins,bins)); tmp_hNout_rmm_da[imass]->Sumw2();
      
      tmp_nin_ee_dy.push_back(0), tmp_nout_ee_dy.push_back(0), tmp_nin_ee_vz.push_back(0), tmp_nout_ee_vz.push_back(0), tmp_nin_ee_data.push_back(0);
      tmp_nin_mm_dy.push_back(0), tmp_nout_mm_dy.push_back(0), tmp_nin_mm_vz.push_back(0), tmp_nout_mm_vz.push_back(0), tmp_nin_mm_data.push_back(0);
      tmp_varin_ee_dy.push_back(0), tmp_varout_ee_dy.push_back(0), tmp_varin_ee_vz.push_back(0), tmp_varout_ee_vz.push_back(0);    
      tmp_varin_mm_dy.push_back(0), tmp_varout_mm_dy.push_back(0), tmp_varin_mm_vz.push_back(0), tmp_varout_mm_vz.push_back(0);
      tmp_nin_em_data.push_back(0);
     }

    nin_kee_data.push_back(tmp_nin_kee_data);
    nin_kmm_data.push_back(tmp_nin_kmm_data);

    hNin_ree_mc.push_back(tmp_hNin_ree_mc); 
    hNout_ree_mc.push_back(tmp_hNout_ree_mc); 
    hNin_rmm_mc.push_back(tmp_hNin_rmm_mc); 
    hNout_rmm_mc.push_back(tmp_hNout_rmm_mc);

    hNin_ree_da.push_back(tmp_hNin_ree_da); 
    hNout_ree_da.push_back(tmp_hNout_ree_da); 
    hNin_rmm_da.push_back(tmp_hNin_rmm_da); 
    hNout_rmm_da.push_back(tmp_hNout_rmm_da);

    nin_ee_dy.push_back(tmp_nin_ee_dy);
    nout_ee_dy.push_back(tmp_nout_ee_dy);
    nin_ee_vz.push_back(tmp_nin_ee_vz);
    nout_ee_vz.push_back(tmp_nout_ee_vz);
    nin_ee_data.push_back(tmp_nin_ee_data);
    nin_mm_dy.push_back(tmp_nin_mm_dy);
    nout_mm_dy.push_back(tmp_nout_mm_dy);
    nin_mm_vz.push_back(tmp_nin_mm_vz);
    nout_mm_vz.push_back(tmp_nout_mm_vz);
    nin_mm_data.push_back(tmp_nin_mm_data);
    varin_ee_dy.push_back(tmp_varin_ee_dy);
    varout_ee_dy.push_back(tmp_varout_ee_dy);
    varin_ee_vz.push_back(tmp_varin_ee_vz);
    varout_ee_vz.push_back(tmp_varout_ee_vz);
    varin_mm_dy.push_back(tmp_varin_mm_dy);
    varout_mm_dy.push_back(tmp_varout_mm_dy);
    varin_mm_vz.push_back(tmp_varin_mm_vz);
    varout_mm_vz.push_back(tmp_varout_mm_vz);
    nin_em_data.push_back(tmp_nin_em_data);
  }
  //*******************************************************
  //Systematic Error on VZ Normalization
  //*******************************************************
  const Double_t vzNormSystematic = 0.10;

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
 
  LeptonScaleLookup trigLookup(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Medium_ele"));
  TH2D *fhDElTightSF  = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Tight_ele"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF ->SetDirectory(0);
  delete fElSF;

  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Medium_mu"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Iso_mu"));
  assert(fhDMuMediumSF);
  assert(fhDMuIsoSF);
  fhDMuMediumSF->SetDirectory(0);
  fhDMuIsoSF->SetDirectory(0);
  delete fMuSF;

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    //TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");

    BareMonteCarlo eventMonteCarlo;
    eventMonteCarlo.setBranchAddresses(the_input_tree);

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

    if(useDYMVA == true){
      the_input_tree->SetBranchAddress("dymva", &dymva_);
      the_input_tree->SetBranchAddress("nlep", &nlep_ );
      the_input_tree->SetBranchAddress("njets", &njets_);
    }

    TNamed *triggerNames = (TNamed*)the_input_file.FindObjectAny("triggerNames");
    char **tokens;
    size_t numtokens;
    tokens = strsplit(triggerNames->GetTitle(), ",", &numtokens);
    if(infilecatv[ifile] == 0){
      for (int i = 0; i < (int)numtokens; i++) {
        printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
      }
    }

    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_input_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      Bool_t passFilter[6] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}
      }
      if(idLep.size()!=idTight.size()) assert(0);
      if(idLep.size()==2) passFilter[0] = kTRUE;
      if(passFilter[0] == kFALSE) continue;

      if(idLep.size() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20) passFilter[1] = kTRUE;
      for (int nt = 0; nt <(int)numtokens; nt++) {
        if((*eventTrigger.triggerFired)[nt] == 0) continue;
        if((strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
           (strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
           (strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*") 	    == 0) ||
           (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")	    == 0) ||
           (strcmp(tokens[nt],"HLT_IsoMu27_v*") 				    == 0) ||
           (strcmp(tokens[nt],"HLT_IsoMu20_v*") 				    == 0) ||
           (strcmp(tokens[nt],"HLT_IsoTkMu20_v*")				    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")	    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele23_WPLoose_Gsf_v*")			    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele22_eta2p1_WP75_Gsf_v*")			    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele27_WPLoose_Gsf_v*")			    == 0) ||
           (strcmp(tokens[nt],"HLT_Ele27_WP85_Gsf_v*")  			    == 0)
           ) passFilter[2] = kTRUE;
      }

      if(infilecatv[ifile] != 0) passFilter[2] = kTRUE; // do not apply trigger filters to MC
      if(passFilter[1] == kFALSE) continue;
      if(passFilter[2] == kFALSE) continue;

      if(goodIsTight == idTight.size()) passFilter[3] = kTRUE;
      if(passFilter[3] == kFALSE) continue;

      passFilter[4] = ((int)(*eventLeptons.pdgId)[idLep[0]]*(int)(*eventLeptons.pdgId)[idLep[1]] < 0);
      if(passFilter[4] == kFALSE) continue;

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 

      double dPhiLepMETMin = 999.;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        if(dPhiLepMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0]))))
           dPhiLepMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));      
      }
      double minMET  = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      double minPMET = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      if(dPhiLepMETMin < TMath::Pi()/2) minPMET = minPMET * sin(dPhiLepMETMin);

      unsigned int nJets = 0;
      bool isBtag = kFALSE;
      double bDiscrMax = 0.0;
      double dPhiJetMET = -1.0;
      double dPhiJetDiLep = -1.0;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;        

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(dPhiJetMET   == -1) dPhiJetMET   = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])))*180./TMath::Pi();
        if(dPhiJetDiLep == -1) dPhiJetDiLep = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventJets.p4)[nj])))*180./TMath::Pi();

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 15 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrMax) bDiscrMax = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) continue;

        nJets++;
      }

      if(useDYMVA == true){
        if(nlep_ != idLep.size()) {printf("PROBLEM nlep %d != %d\n",(int)nlep_,(int)idLep.size()); assert(1); return;}
        if(njets_ != nJets) {printf("PROBLEM njet %d != %d\n",(int)njets_,nJets); assert(1); return;}
      }

      if(nJets <= 2) passFilter[5] = kTRUE;  	    
      if(passFilter[5] == kFALSE) continue;

      Int_t typeLep = 2;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) typeLep = 0;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) typeLep = 1;

      if(infilecatv[ifile] == 0 && TMath::Abs(dilep.M()-mZ)<MassZCut) {
        if(typeLep == 0) nin_kmm_data[nJets]++;
        if(typeLep == 1) nin_kee_data[nJets]++;
      }

      if(minPMET <= 20) continue;
      double varMet = minPMET;
      if(varMet>=70) varMet=69;
      if(useDYMVA == kTRUE) {
        varMet = dymva_;
      }

      // trigger efficiency
      double trigEff = 1.0;
      if(infilecatv[ifile] != 0) {
        trigEff = trigLookup.GetExpectedTriggerEfficiency(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),
        						  ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta(),((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),
        						 TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]));
      }
      // luminosity
      double theLumi  = 1.0; if(infilecatv[ifile] != 0) theLumi  = lumi;
      // pile-up
      double puWeight = 1.0; if(infilecatv[ifile] != 0) puWeight = nPUScaleFactor(fhDPU, (double)eventMonteCarlo.puTrueInt);
      // lepton efficiency
      double effSF = 1.0;
      if(infilecatv[ifile] != 0){
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
		period,typeLepSel.Data(),fhDMuMediumSF,fhDMuIsoSF,fhDElMediumSF,fhDElTightSF);
        }
      }

      double totalWeight = eventMonteCarlo.mcWeight*theLumi*puWeight*effSF*theMCPrescale*trigEff;

      vector<bool> isGenDupl;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
        isGenDupl.push_back(0);
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) isGenDupl[ngen0] = 1;
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen0]) != 13) continue;
        for(int ngen1=ngen0+1; ngen1<eventMonteCarlo.p4->GetEntriesFast(); ngen1++) {
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
      if(infilecatv[ifile] != 0 && goodIsGenLep != isGenLep.size()) totalWeight = 0.0;

      //loop over analyses
      for(Int_t imass=0; imass<nmass; imass++) {
        // MET related cuts need to be applied after the Rout/in computation
        double theCutMassHigh = 1000.;
	if(dilep.Pt() <= 45.0) continue;
        if(bDiscrMax > 0.605 || idSoft.size() != 0) continue;
	if(dilep.M() <= 12.0) continue;

        if(infilecatv[ifile] == 1){ // Z MC
	  if     (TMath::Abs(dilep.M()-mZ)<MassZCut) {
	    if     (typeLep == 1) {
              hNin_ree_mc[nJets][imass]->Fill(varMet,totalWeight);		 
            }
	    else if(typeLep == 0) {
              hNin_rmm_mc[nJets][imass]->Fill(varMet,totalWeight);		     
            }
          } 
	  else if(TMath::Abs(dilep.M()-mZ) >= 15 && dilep.M() < theCutMassHigh) {
	    if     (typeLep == 1) {
              hNout_ree_mc[nJets][imass]->Fill(varMet,totalWeight);  
            }
	    else if(typeLep == 0) {
              hNout_rmm_mc[nJets][imass]->Fill(varMet,totalWeight);
            }
          }
        }
 
        double theDataWeight = 0.0;
        if     (infilecatv[ifile] == 0 && (typeLep == 0 || typeLep == 1)) theDataWeight =  1.0;
	else if(infilecatv[ifile] == 0					) theDataWeight = -1.0;

	// In reality ee == ee+mm for data
	if(theDataWeight != 0.0){
	  if     (TMath::Abs(dilep.M()-mZ)<MassZCut) {
            hNin_ree_da[nJets][imass]->Fill(varMet,theDataWeight);		
          } 
	  else if(TMath::Abs(dilep.M()-mZ) >= 15 && dilep.M() < theCutMassHigh) {
            hNout_ree_da[nJets][imass]->Fill(varMet,theDataWeight);     
          }
        }

	bool passMET = varMet > 45;
	if(useDYMVA == kTRUE) passMET = dymva_ > 0.3;
	if(passMET == kFALSE) continue;

        if(infilecatv[ifile] == 1){ // Z MC
	  if     (TMath::Abs(dilep.M()-mZ)<MassZCut) {
	    if     (typeLep == 1) { 
	      nin_ee_dy[nJets][imass]+=totalWeight; 
	      varin_ee_dy[nJets][imass]+=totalWeight*totalWeight; 
	    }	
	    else if(typeLep == 0) { 
	      nin_mm_dy[nJets][imass]+=totalWeight; 
	      varin_mm_dy[nJets][imass]+=totalWeight*totalWeight;
	    }

           } 
	   else if(TMath::Abs(dilep.M()-mZ) >= 15 && dilep.M() < theCutMassHigh) {
	    if     (typeLep == 1) {
	      nout_ee_dy[nJets][imass]+=totalWeight; 
	      varout_ee_dy[nJets][imass]+=totalWeight*totalWeight;
	    }  
	    else if(typeLep == 0) {
	      nout_mm_dy[nJets][imass]+=totalWeight;
	      varout_mm_dy[nJets][imass]+=totalWeight*totalWeight;
	    }
          }
        }

        //In Z peak region
	if(TMath::Abs(dilep.M()-mZ)<MassZCut) {

          if(typeLep == 1) {
            if(infilecatv[ifile] == 0) { 
              nin_ee_data[nJets][imass]++; 
            }
            if(infilecatv[ifile] == 2)  {
              nin_ee_vz[nJets][imass]+=totalWeight;
              varin_ee_vz[nJets][imass]+=totalWeight*totalWeight;
            }
          }
	
          if(typeLep == 0) {
            if(infilecatv[ifile] == 0) { 
              nin_mm_data[nJets][imass]++;
            }	
            if(infilecatv[ifile] == 2) {
              nin_mm_vz[nJets][imass]+=totalWeight;
              varin_mm_vz[nJets][imass]+=totalWeight*totalWeight;
            }		   
          }

	  if(infilecatv[ifile] == 0 && typeLep == 2) {
            nin_em_data[nJets][imass]++; 
          }

        } 
        // Out of Z peak region
        else if(TMath::Abs(dilep.M()-mZ) >= 15 && dilep.M() < theCutMassHigh) {

	  if(typeLep == 0) {
            if (infilecatv[ifile] == 2) {
              nout_ee_vz[nJets][imass]+=totalWeight;
              varout_ee_vz[nJets][imass]+=totalWeight*totalWeight;
            }
	  }
	
	  if(typeLep == 1) {
            if (infilecatv[ifile] == 2) {
              nout_mm_vz[nJets][imass]+=totalWeight;
              varout_mm_vz[nJets][imass]+=totalWeight*totalWeight;
            }
	  }
        }

      }
    }
  } // end of chain

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  char DYEstimateTableName[200];
  sprintf(DYEstimateTableName,"DYEstimateTable.txt");
  ofstream fout(DYEstimateTableName);

  vector<vector<Double_t> > DYBkgScaleFactorHiggsSelection;
  vector<vector<Double_t> > DYBkgScaleFactorHiggsSelectionErr;
  vector<Double_t> DYBkgScaleFactorWWPreselection;
  vector<Double_t> DYBkgScaleFactorWWPreselectionErr;

  for(UInt_t jetIndex = 0; jetIndex < 3; ++jetIndex) {
    
    vector<Double_t> tmpDYBkgScaleFactorHiggsSelection;
    vector<Double_t> tmpDYBkgScaleFactorHiggsSelectionErr;
    Double_t tmpDYBkgScaleFactorWWPreselection;
    Double_t tmpDYBkgScaleFactorWWPreselectionErr;

    fout << "************************************************************************\n";
    fout << jetIndex << "-Jet Bin DY Bkg Scale Factor Computation\n";
    fout << "************************************************************************\n";

    Double_t k    = sqrt(nin_kee_data[jetIndex]/nin_kmm_data[jetIndex]);
    Double_t kerr = 0.5*k*sqrt(1.0/nin_kee_data[jetIndex] + 1.0/nin_kmm_data[jetIndex]);
    char buffer[200];
  
    fout << "Electron to muon efficiency ratio is " << k << " +/- " << kerr << endl;

    fout << endl;
    fout << jetIndex << "-jet bin summary:" << endl;
    fout << setw(4) << "sel" << "   ";
    fout << setw(25) << "R_out/in       ";
    fout << setw(15) << "mm/em/ee";
    fout << setw(20) << "N_in (OF,VZ sub)";
    fout << setw(20) << "N_in (data/MC)";
    fout << setw(20) << "   N_out data   ";
    fout << setw(20) << "N_out MC  ";
    fout << setw(20) << "N_out (data/MC)" << endl;


    for(Int_t imass=0; imass<nmass; imass++) {
      fout << setw(4) << mH[imass] << "   ";
    
      //
      // compute Routin from MC
      //

      Int_t MetBinToComputeRoutin = nbins-1;

      Double_t nout_ee   = 0;
      Double_t errout_ee = 0;
      Double_t nout_mm   = 0;
      Double_t errout_mm = 0;
      Double_t nout_ll   = 0;
      Double_t errout_ll = 0;      
      Double_t nin_ee   = 0;
      Double_t errin_ee = 0;
      Double_t nin_mm   = 0;
      Double_t errin_mm = 0;
      Double_t nin_ll   = 0;
      Double_t errin_ll = 0;      
      Double_t Ree        = 0;
      Double_t ReeErrStat = 0;
      Double_t ReeErrSyst = 0;      
      Double_t Rmm        = 0;
      Double_t RmmErrStat = 0;
      Double_t RmmErrSyst = 0;
      Double_t Rll = 0;
      Double_t RllErrStat = 0;
      Double_t RllErrSyst = 0;

      if(useRFromData[jetIndex] == false){
	nout_ee   = hNout_ree_mc[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errout_ee = hNout_ree_mc[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);
	nout_mm   = hNout_rmm_mc[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errout_mm = hNout_rmm_mc[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);
	nout_ll   = nout_ee+nout_mm;
	errout_ll = sqrt(errout_ee*errout_ee + errout_mm*errout_mm);

	nin_ee   = hNin_ree_mc[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errin_ee = hNin_ree_mc[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);
	nin_mm   = hNin_rmm_mc[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errin_mm = hNin_rmm_mc[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);
	nin_ll   = nin_ee+nin_mm;
	errin_ll = sqrt(errin_ee*errin_ee + errin_mm*errin_mm);

	Ree        = nout_ee/nin_ee;
	ReeErrStat = Ree*sqrt(errin_ee*errin_ee/nin_ee/nin_ee + errout_ee*errout_ee/nout_ee/nout_ee);
	ReeErrSyst = computeSyst(hNout_ree_mc[jetIndex][imass],hNin_ree_mc[jetIndex][imass], MetBinToComputeRoutin, 0.4);    

	Rmm        = nout_mm/nin_mm;
	RmmErrStat = Rmm*sqrt(errin_mm*errin_mm/nin_mm/nin_mm + errout_mm*errout_mm/nout_mm/nout_mm);
	RmmErrSyst = computeSyst(hNout_rmm_mc[jetIndex][imass],hNin_rmm_mc[jetIndex][imass], MetBinToComputeRoutin, 0.4);

	TH1D *hout = (TH1D*)hNout_ree_mc[jetIndex][imass]->Clone("hout");
	hout->Add(hNout_rmm_mc[jetIndex][imass]);
	TH1D *hin = (TH1D*)hNin_ree_mc[jetIndex][imass]->Clone("hin");
	hin->Add(hNin_rmm_mc[jetIndex][imass]);
	Rll        = nout_ll/nin_ll;
	RllErrStat = Rll*sqrt(errin_ll*errin_ll/nin_ll/nin_ll + errout_ll*errout_ll/nout_ll/nout_ll);
	RllErrSyst = computeSyst(hout,hin, MetBinToComputeRoutin, 0.4);

	delete hout;
	delete hin;
      } else {
	nout_ll   = hNout_ree_da[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errout_ll = hNout_ree_da[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);

	nin_ll   = hNin_ree_da[jetIndex][imass]->GetBinContent(MetBinToComputeRoutin);
	errin_ll = hNin_ree_da[jetIndex][imass]->GetBinError(MetBinToComputeRoutin);

	TH1D *hout = (TH1D*)hNout_ree_da[jetIndex][imass]->Clone("hout");
	TH1D *hin = (TH1D*)hNin_ree_da[jetIndex][imass]->Clone("hin");

	Rll        = nout_ll/nin_ll;
	RllErrStat = Rll*sqrt(errin_ll*errin_ll/nin_ll/nin_ll + errout_ll*errout_ll/nout_ll/nout_ll);
	Double_t rErrorMax = 100000.4;
	cout << "M: " << mH[imass] << endl;
	RllErrSyst = computeSyst(hout,hin,MetBinToComputeRoutin,rErrorMax,true,true);
	delete hout;
	delete hin;
      }

      if(Rll <= 0) {Rll = 0.10; RllErrStat = 0.10; RllErrSyst = 0.10;}
      // DO NOT ALLOW FOR MORE THAN 100% uncertainty
      if(RllErrSyst > 1.0*Rll) RllErrSyst = 1.0*Rll;
      // DO NOT ALLOW FOR LESS THAN 30% uncertainty
      if(RllErrSyst < 0.3*Rll) RllErrSyst = 0.3*Rll;

      sprintf(buffer,"%.2f +/- %.2f +/- %.2f",Rll,RllErrStat,RllErrSyst);
//     sprintf(buffer,"%.3f +/- %.3f +/- %.3f",Ree,ReeErrStat,ReeErrSyst);
//     sprintf(buffer,"%.3f +/- %.3f +/- %.3f",Rmm,RmmErrStat,RmmErrSyst);
      fout << setw(25) << buffer; 

      //
      // raw in-yields in data
      //
      sprintf(buffer,"%i/%i/%i",Int_t(nin_mm_data[jetIndex][imass]), Int_t(nin_em_data[jetIndex][imass]), Int_t(nin_ee_data[jetIndex][imass]));
      fout << setw(15) << buffer;

      //
      // in-yields in data after OF and VZ subtraction
      //
      Double_t nof = nin_em_data[jetIndex][imass];
      Double_t nin_ee_sub, err_ee_sub,nin_mm_sub,err_mm_sub,nin_ll_sub,err_ll_sub;

      //use Opposite Flavor data for bkg subtraction
      nin_ee_sub = nin_ee_data[jetIndex][imass] - 0.5*k*nof - nin_ee_vz[jetIndex][imass];
      err_ee_sub = sqrt(nin_ee_data[jetIndex][imass] + 0.5*0.5*k*k*nof*nof*(kerr*kerr/k/k + 1.0/nof) + varin_ee_vz[jetIndex][imass] + pow(nin_ee_vz[jetIndex][imass]*vzNormSystematic,2) );
      nin_mm_sub = nin_mm_data[jetIndex][imass] - 0.5/k*nof - nin_mm_vz[jetIndex][imass];
      err_mm_sub = sqrt(nin_mm_data[jetIndex][imass] + 0.5*0.5/k/k*nof*nof*(kerr*kerr/k/k + 1.0/nof) + varin_mm_vz[jetIndex][imass] + pow(nin_mm_vz[jetIndex][imass]*vzNormSystematic,2));
      nin_ll_sub = nin_ee_sub + nin_mm_sub;
      err_ll_sub = sqrt(nin_mm_data[jetIndex][imass] + nin_ee_data[jetIndex][imass] 
        			 + 0.5*0.5*( (k+1.0/k)*(k+1.0/k)*nof*nof*kerr*kerr + (k+1.0/k)*(k+1.0/k)*nof ) 
        			 + varin_mm_vz[jetIndex][imass] + varin_ee_vz[jetIndex][imass] + pow((nin_ee_vz[jetIndex][imass]+nin_mm_vz[jetIndex][imass])*vzNormSystematic,2));
      if(nin_ee_sub <= 0) nin_ee_sub = 1;
      if(nin_mm_sub <= 0) nin_mm_sub = 1;
      if(nin_ll_sub <= 0) nin_ll_sub = 1;

      sprintf(buffer,"  %.2f +/- %.2f : %.2f : %.2f + %.2f  ",nin_ll_sub,err_ll_sub, 0.5*k*nof + 0.5/k*nof, nin_ee_vz[jetIndex][imass], nin_mm_vz[jetIndex][imass]);
      fout << setw(20) << buffer;
    
      //
      // in-yield data/MC scale factor
      //

      Double_t sfin_ee     = nin_ee_sub/nin_ee_dy[jetIndex][imass];
      Double_t sfin_ee_err = sfin_ee*sqrt(err_ee_sub*err_ee_sub/nin_ee_sub/nin_ee_sub + (varin_ee_dy[jetIndex][imass])*(varin_ee_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass]));
      Double_t sfin_mm     = nin_mm_sub/nin_mm_dy[jetIndex][imass];
      Double_t sfin_mm_err = sfin_mm*sqrt(err_mm_sub*err_mm_sub/nin_mm_sub/nin_mm_sub + (varin_mm_dy[jetIndex][imass])*(varin_mm_dy[jetIndex][imass])/(nin_mm_dy[jetIndex][imass])/(nin_mm_dy[jetIndex][imass]));
      Double_t sfin_ll     = nin_ll_sub/(nin_ee_dy[jetIndex][imass]+nin_mm_dy[jetIndex][imass]);
//       Double_t sfin_ll_err = sfin_ll*sqrt(err_ll_sub*err_ll_sub/nin_ll_sub/nin_ll_sub 
//                                           + (varin_ee_dy[jetIndex][imass]+varin_mm_dy[jetIndex][imass])*(varin_ee_dy[jetIndex][imass]+varin_mm_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass]+nin_mm_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass]+nin_mm_dy[jetIndex][imass]));
      Double_t sfin_ll_err = sfin_ll*sqrt(err_ll_sub*err_ll_sub/nin_ll_sub/nin_ll_sub 
                                          + (varin_ee_dy[jetIndex][imass]+varin_mm_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass]+nin_mm_dy[jetIndex][imass])/(nin_ee_dy[jetIndex][imass]+nin_mm_dy[jetIndex][imass]));
   
      sprintf(buffer,"%.2f +/- %.2f",sfin_ll,sfin_ll_err);
//    sprintf(buffer,"%.2f +/- %.2f",sfin_ee,sfin_ee_err);
//    sprintf(buffer,"%.2f +/- %.2f",sfin_mm,sfin_mm_err);
      fout << setw(20) << buffer;
    
      //
      // out-yield prediction
      //
      Double_t nout_ee_pre = Ree*nin_ee_sub;
      Double_t nout_ee_sta = nout_ee_pre*sqrt(ReeErrStat*ReeErrStat/Ree/Ree + err_ee_sub*err_ee_sub/nin_ee_sub/nin_ee_sub);
      Double_t nout_ee_sys = nout_ee_pre*ReeErrSyst/Ree;
      Double_t nout_ee_err = sqrt(nout_ee_sys*nout_ee_sys + nout_ee_sta*nout_ee_sta);
      Double_t nout_mm_pre = Rmm*nin_mm_sub;
      Double_t nout_mm_sta = nout_mm_pre*sqrt(RmmErrStat*RmmErrStat/Rmm/Rmm + err_mm_sub*err_mm_sub/nin_mm_sub/nin_mm_sub);
      Double_t nout_mm_sys = nout_mm_pre*RmmErrSyst/Rmm;
      Double_t nout_mm_err = sqrt(nout_mm_sys*nout_mm_sys + nout_mm_sta*nout_mm_sta);
      Double_t nout_ll_pre = Rll*nin_ll_sub;
      Double_t nout_ll_sta = nout_ll_pre*sqrt((RllErrStat*RllErrStat)/Rll/Rll + err_ll_sub*err_ll_sub/nin_ll_sub/nin_ll_sub);
      Double_t nout_ll_sys = nout_ll_pre*RllErrSyst/Rll;
      Double_t nout_ll_err = sqrt(nout_ll_sys*nout_ll_sys + nout_ll_sta*nout_ll_sta);

      sprintf(buffer,"%.2f +/- %.2f +/- %.2f",nout_ll_pre,nout_ll_sta,nout_ll_sys);
//    sprintf(buffer,"%.2f +/- %.2f +/- %.2f",nout_ee_pre,nout_ee_sta,nout_ee_sys);
//    sprintf(buffer,"%.2f +/- %.2f +/- %.2f",nout_mm_pre,nout_mm_sta,nout_mm_sys);
      fout << "   " << setw(20) << buffer;
    
      //
      // out-yield in MC
      //
      printf("SSS %f %f %f\n",nout_ee_dy[jetIndex][imass]+nout_mm_dy[jetIndex][imass],nout_ee_dy[jetIndex][imass],nout_mm_dy[jetIndex][imass]);
      sprintf(buffer,"%.2f +/- %.2f",nout_ee_dy[jetIndex][imass]+nout_mm_dy[jetIndex][imass],sqrt(varout_ee_dy[jetIndex][imass] + varout_mm_dy[jetIndex][imass]));
//    sprintf(buffer,"%.2f +/- %.2f",nout_ee_dy[jetIndex][imass],sqrt(varout_ee_dy[jetIndex][imass]));
//    sprintf(buffer,"%.2f +/- %.2f",nout_mm_dy[jetIndex][imass],sqrt(varout_mm_dy[jetIndex][imass]));
      fout << setw(20) << buffer;
    
      //
      // out-yield data/MC scale factor
      //
      Double_t sfout_ee     = nout_ee_pre/nout_ee_dy[jetIndex][imass];
      Double_t sfout_ee_err = sfout_ee*sqrt(nout_ee_err*nout_ee_err/nout_ee_pre/nout_ee_pre );
      Double_t sfout_mm     = nout_mm_pre/nout_mm_dy[jetIndex][imass];
      Double_t sfout_mm_err = sfout_mm*sqrt(nout_mm_err*nout_mm_err/nout_mm_pre/nout_mm_pre );
      Double_t sfout_ll     = nout_ll_pre/(nout_ee_dy[jetIndex][imass]+nout_mm_dy[jetIndex][imass]);
      Double_t sfout_ll_err = sfout_ll*sqrt(nout_ll_err*nout_ll_err/nout_ll_pre/nout_ll_pre);
    
      sprintf(buffer,"%.2f +/- %.2f",sfout_ll,sfout_ll_err);
//    sprintf(buffer,"%.2f +/- %.2f",sfout_ee,sfout_ee_err);
//    sprintf(buffer,"%.2f +/- %.2f",sfout_mm,sfout_mm_err);
      fout << setw(20) << buffer;        
      fout << endl;
      if (imass == 0) {
        tmpDYBkgScaleFactorWWPreselection    = sfout_ll;
        tmpDYBkgScaleFactorWWPreselectionErr = sfout_ll_err;
      } else {
        // BIG CHANGE, we quote
        //tmpDYBkgScaleFactorHiggsSelection.push_back(sfout_ll);
        //tmpDYBkgScaleFactorHiggsSelectionErr.push_back(sfout_ll_err);
        tmpDYBkgScaleFactorHiggsSelection.push_back(nout_ll_pre);
        tmpDYBkgScaleFactorHiggsSelectionErr.push_back(sqrt(nout_ll_sta*nout_ll_sta+nout_ll_sys*nout_ll_sys));
      }
    }

    DYBkgScaleFactorHiggsSelection.push_back(tmpDYBkgScaleFactorHiggsSelection);
    DYBkgScaleFactorHiggsSelectionErr.push_back(tmpDYBkgScaleFactorHiggsSelectionErr);
    DYBkgScaleFactorWWPreselection.push_back(tmpDYBkgScaleFactorWWPreselection);
    DYBkgScaleFactorWWPreselectionErr.push_back(tmpDYBkgScaleFactorWWPreselectionErr);
  }

  fout.close();


  //***************************************************************************
  // Generate DY Scale Factor and Systematics Code for card creation
  //***************************************************************************
  char DYBkgScaleFactorsName[200];
  sprintf(DYBkgScaleFactorsName,"DYBkgScaleFactors.h");
  ofstream outf(DYBkgScaleFactorsName);

  outf << "static Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {" << endl;
  
  outf << "  Int_t mHiggs[" << nmass-1 << "] = {";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << mH[i+1];
    if (i < nmass-1-1) outf << ",";    
  }
  outf << "};" << endl;

  outf << "  Double_t DYBkgScaleFactorWWPreselection[3] = { " 
       << DYBkgScaleFactorWWPreselection[0] << ", "
       << DYBkgScaleFactorWWPreselection[1] << ", "
       << DYBkgScaleFactorWWPreselection[2] << " "
       << " };" << endl;

  outf << "  Double_t DYBkgScaleFactorHiggsSelection[3][" << nmass-1 << "] = { " << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << DYBkgScaleFactorHiggsSelection[0][i];    
    if (i < nmass-1-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << DYBkgScaleFactorHiggsSelection[1][i];    
    if (i < nmass-1-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << DYBkgScaleFactorHiggsSelection[2][i];    
    if (i < nmass-1-1) outf << ",";
  }
  outf << "} };" << endl;

  outf << "  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];" << endl;
  
  outf << "  Int_t massIndex = -1;" << endl;
  outf << "  for (UInt_t m=0; m < " << nmass-1 << " ; ++m) {" << endl;
  outf << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outf << "  }" << endl;
  outf << "  if (massIndex >= 0) {" << endl;  
  outf << "    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];" << endl;
  outf << "  } else {" << endl;
  outf << "    return DYBkgScaleFactorWWPreselection[jetBin];" << endl;
  
  outf << "  }" << endl;
  outf << "}" << endl;
  outf << endl;


  outf << "static Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {" << endl;
  
  outf << "  Int_t mHiggs[" << nmass-1 << "] = {";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << mH[i+1];
    if (i < nmass-1-1) outf << ",";    
  }
  outf << "};" << endl;

  outf << "  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { " 
       << (1.0 + DYBkgScaleFactorWWPreselectionErr[0]/DYBkgScaleFactorWWPreselection[0]) << ", "
       << (1.0 + DYBkgScaleFactorWWPreselectionErr[1]/DYBkgScaleFactorWWPreselection[1]) << ", "
       << (1.0 + DYBkgScaleFactorWWPreselectionErr[2]/DYBkgScaleFactorWWPreselection[2]) << " "
       << " };" << endl;

  outf << "  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][" << nmass-1 << "] = { " << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << (1.0 + DYBkgScaleFactorHiggsSelectionErr[0][i] / DYBkgScaleFactorHiggsSelection[0][i]);
    if (i < nmass-1-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << (1.0 + DYBkgScaleFactorHiggsSelectionErr[1][i] / DYBkgScaleFactorHiggsSelection[1][i]);
    if (i < nmass-1-1) outf << ",";
  }
  outf << "}," << endl;
  outf << "    { ";
  for (UInt_t i = 0; i < nmass-1 ; ++i) {
    outf << (1.0 + DYBkgScaleFactorHiggsSelectionErr[2][i] / DYBkgScaleFactorHiggsSelection[2][i]);
    if (i < nmass-1-1) outf << ",";
  }
  outf << "} };" << endl;

  outf << "  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];" << endl;
  
  outf << "  Int_t massIndex = -1;" << endl;
  outf << "  for (UInt_t m=0; m < " << nmass-1 << " ; ++m) {" << endl;
  outf << "    if (mH == mHiggs[m]) massIndex = m;" << endl;
  outf << "  }" << endl;
  outf << "  if (massIndex >= 0) {" << endl;  
  outf << "    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];" << endl;
  outf << "  } else {" << endl;
  outf << "    return DYBkgScaleFactorWWPreselectionKappa[jetBin];" << endl;
  
  outf << "  }" << endl;
  outf << "}" << endl;
  outf << endl;

  outf.close();
}

//--------------------------------------------------------------------------------------------------
Double_t computeSyst(const TH1D *hout, const TH1D *hin, Int_t binUsed, Double_t rErrorMax, Bool_t useRFromData, Bool_t isDebug)
{
  const Int_t nbins = hout->GetXaxis()->GetNbins();
  Double_t r0=(hout->GetBinContent(binUsed))/(hin->GetBinContent(binUsed));
  Double_t dr=0;
  for(Int_t ibin=1; ibin<=nbins; ibin++) {
    if(hin->GetBinContent(ibin) <= 0) continue;
    if(useRFromData == true && ibin == nbins) continue;
    Double_t r = (hout->GetBinContent(ibin))/(hin->GetBinContent(ibin));

    //only consider last bin for systematics if the uncertainty is less than rErrorMax to protect against statistical fluctuations
    Double_t rError = sqrt( pow(hout->GetBinError(ibin)/hout->GetBinContent(ibin),2) + pow(hin->GetBinError(ibin)/hin->GetBinContent(ibin),2));

    if(r < 0.0) r = 0.0;
    if(isDebug == true) printf("bin(%1d) Nin = %8.3f Nout = %8.3f R0 = %5.3f R = %5.3f dr = %5.3f drMax = %5.3f rError = %5.3f\n",ibin,hin->GetBinContent(ibin),hout->GetBinContent(ibin),r0,r,fabs(r-r0),dr,rError);

    if(rError > rErrorMax) continue;

    if(fabs(r-r0) > dr) dr = fabs(r-r0);
  }
  return dr;
}
