#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"
#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"

#include "VVScattering/macros/74x/factors.h"

const bool usePureMC = false;
const int etaBins = 5;
const TString typeLepSel = "default";

void topAnalysis(
 Int_t typeSel = 6,
 Int_t period = 1
 ){

  TString filesPath  = "/scratch5/ceballos/ntuples_weights_74x/met_";
  Double_t lumi = 0.0715;
  if(period == 1) lumi = 2.263;
  enum { kOther, kTTBAR, kTW, kData };

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  // careful, ttbar myst be category 3 and tW must be category 4
  TString puPath = "";
  if      (period==1){
  puPath = "/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src/VVScattering/data/74x/puWeights_13TeV_25ns.root";
  infilenamev.push_back(Form("%sdata_AOD_Run2015C1_25ns.root",filesPath.Data()));												  infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D3_25ns.root",filesPath.Data()));												  infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D4_25ns.root",filesPath.Data()));												  infilecatv.push_back(0);
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));					          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));                                          infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));         infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(2);
  infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(3);
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));      infilecatv.push_back(4);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sGluGluToZZTo2e2mu_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(7);
  infilenamev.push_back(Form("%sGluGluToZZTo2e2tau_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(7);
  infilenamev.push_back(Form("%sGluGluToZZTo2mu2tau_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(7);
  infilenamev.push_back(Form("%sGluGluToZZTo4e_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				  infilecatv.push_back(7);
  infilenamev.push_back(Form("%sGluGluToZZTo4mu_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				  infilecatv.push_back(7);
  infilenamev.push_back(Form("%sGluGluToZZTo4tau_BackgroundOnly_13TeV_MCFM+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sWZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));                          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));        infilecatv.push_back(7);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));         infilecatv.push_back(7);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(5);
  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  		  infilecatv.push_back(6);
  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3+AODSIM.root",filesPath.Data()));	          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  	          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data())); 		          infilecatv.push_back(7);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			          infilecatv.push_back(7);
  }
  else {assert(0);}

  //infilenamev.push_back(Form("nero.root"));     infilecatv.push_back(0);

  if(infilenamev.size() != infilecatv.size()) assert(0);

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
 
  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  double btag_central_All_2j_den[4],      btag_central_All_2j_num[4],      btag_central_All_1j_den[4],      btag_central_All_1j_num[4];
  double btag_central_All_2j_den_error[4],btag_central_All_2j_num_error[4],btag_central_All_1j_den_error[4],btag_central_All_1j_num_error[4];
  double btag_central_2j_den[4][etaBins],      btag_central_2j_num[4][etaBins];
  double btag_central_2j_den_error[4][etaBins],btag_central_2j_num_error[4][etaBins];
  double btag_highestpt_2j_den[4][1],      btag_highestpt_2j_num[4][1],      btag_highestpt_1j_den[4][1],      btag_highestpt_1j_num[4][1];
  double btag_highestpt_2j_den_error[4][1],btag_highestpt_2j_num_error[4][1],btag_highestpt_1j_den_error[4][1],btag_highestpt_1j_num_error[4][1];
  double btag_lowpt_1j_den[4][1],      btag_lowpt_1j_num[4][1],      btag_lowpt_0j_den[4][1],      btag_lowpt_0j_num[4][1];
  double btag_lowpt_1j_den_error[4][1],btag_lowpt_1j_num_error[4][1],btag_lowpt_0j_den_error[4][1],btag_lowpt_0j_num_error[4][1];
  for (int j = 0; j < 4; j++) {
    btag_central_All_2j_den[j] = 0.0; btag_central_All_2j_den_error[j] = 0.0;
    btag_central_All_2j_num[j] = 0.0; btag_central_All_2j_num_error[j] = 0.0;
    for (int i = 0; i < etaBins; i++) {
      btag_central_2j_den[j][i] = 0.0; btag_central_2j_den_error[j][i] = 0.0;
      btag_central_2j_num[j][i] = 0.0; btag_central_2j_num_error[j][i] = 0.0;
    }
    for (int i = 0; i < 1; i++) {
      btag_highestpt_2j_den[j][i] = 0.0; btag_highestpt_2j_den_error[j][i] = 0.0;
      btag_highestpt_2j_num[j][i] = 0.0; btag_highestpt_2j_num_error[j][i] = 0.0;
      btag_highestpt_1j_den[j][i] = 0.0; btag_highestpt_1j_den_error[j][i] = 0.0;
      btag_highestpt_1j_num[j][i] = 0.0; btag_highestpt_1j_num_error[j][i] = 0.0;
      btag_lowpt_1j_den[j][i] = 0.0; btag_lowpt_1j_den_error[j][i] = 0.0;
      btag_lowpt_1j_num[j][i] = 0.0; btag_lowpt_1j_num_error[j][i] = 0.0;
      btag_lowpt_0j_den[j][i] = 0.0; btag_lowpt_0j_den_error[j][i] = 0.0;
      btag_lowpt_0j_num[j][i] = 0.0; btag_lowpt_0j_num_error[j][i] = 0.0;
    }
  }

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

    TNamed *triggerNames = (TNamed*)the_input_file.FindObjectAny("triggerNames");
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

    for (int i=0; i<(int)the_input_tree->GetEntries()/1; ++i) {
      the_input_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      Bool_t passFilter[6] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 20) passFilter[0] = kTRUE;
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
           ) passFilter[1] = kTRUE;
      }

      if(passFilter[0] == kFALSE) continue;
      if(passFilter[1] == kFALSE) continue;

      vector<int> idLep; vector<int> idTight; vector<int> idSoft; unsigned int goodIsTight = 0;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep); goodIsTight++;}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP)== BareLeptons::LepSoftIP){idSoft.push_back(nlep);}
      }
      if(idLep.size()!=idTight.size()) assert(0);
      if(idLep.size()==2) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 20) continue;

      if(idTight[0]==1&&idTight[1]==1) passFilter[3] = kTRUE;
      if(usePureMC ==  true && passFilter[3] == kFALSE) continue;

      if     (typeSel == 0 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) passFilter[4] = kTRUE;
      else if(typeSel == 3 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) passFilter[4] = kTRUE;
      else if(typeSel == 5 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]))         passFilter[4] = kTRUE;
      else if(typeSel == 6 && TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])!=TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]))         passFilter[4] = kTRUE;
      else if(typeSel == 4)                                                                                                               passFilter[4] = kTRUE;

      passFilter[4] = passFilter[4] * ((int)(*eventLeptons.pdgId)[idLep[0]]*(int)(*eventLeptons.pdgId)[idLep[1]] < 0);

      if(passFilter[4] == kFALSE) continue;

      double dPhiLepMETMin = 999.;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        if(dPhiLepMETMin > TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0]))))
           dPhiLepMETMin = TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));      
      }
      double minMET  = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      double minPMET = TMath::Min(((TLorentzVector*)(*eventMet.p4)[0])->Pt(),(double)eventMet.trackMet->Pt());
      if(dPhiLepMETMin < TMath::Pi()/2) minPMET = minPMET * sin(dPhiLepMETMin);

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 
      if(dilep.M() > 12 && minPMET > 20 && dilep.Pt() > 30) passFilter[5] = kTRUE;
      if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])){
	if(TMath::Abs(dilep.M()-91.1876) <= 15.0 || minPMET <= 45 || dilep.Pt() <= 45) passFilter[5] = kFALSE;
      }
      if(passFilter[5] == kFALSE) continue;

      vector<int> idJet;
      double bDiscrLowPt = -999.;
      int n_jet[2]  = {-999, -999};
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;        

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.4) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

	if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30 && 
	   (float)(*eventJets.bDiscr)[nj] > bDiscrLowPt) bDiscrLowPt = (float)(*eventJets.bDiscr)[nj];

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) continue;
	idJet.push_back(nj);
      }

      vector<int> idSoftNotInJet;
      for(unsigned nl=0; nl<idSoft.size(); nl++){
        bool isInsideJet = kFALSE;
        for(unsigned nj=0; nj<idJet.size(); nj++){
          if(((TLorentzVector*)(*eventJets.p4)[idJet[nj]])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idSoft[nl]])) < 0.4) {isInsideJet = kTRUE; break;}
        }
	if(isInsideJet == kFALSE) idSoftNotInJet.push_back(idSoft[nl]);
      }

      // begin event weighting
      vector<bool> isGenDupl;
      for(int ngen0=0; ngen0<eventMonteCarlo.p4->GetEntriesFast(); ngen0++) {
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
      // luminosity
      double theLumi  = 1.0; if(infilecatv[ifile] != 0) theLumi  = lumi;
      // pile-up
      double puWeight = 1.0; if(infilecatv[ifile] != 0) puWeight = nPUScaleFactor(fhDPU, (double)eventVertex.npv);
      // lepton efficiency
      double effSF = 1.0;
      if(infilecatv[ifile] != 0){
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          effSF = effSF * effScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
        }
      }

      // fake rate
      int theCategory = infilecatv[ifile];
      double fakeSF = 1.0;
      if(usePureMC == false){
        if     (theCategory == 5){ // remove W+jets from MC
          fakeSF = 0.0;
        }
        else if(theCategory == 2 && goodIsTight != idTight.size()){ // remove Z+jets from MC as fakeable objects
          fakeSF = 0.0;
        }
        else if((infilecatv[ifile] == 0 || infilecatv[ifile] == 6 || goodIsGenLep == isGenLep.size()) && goodIsTight != idTight.size()){ // add W+jets from data
          for(unsigned int nl=0; nl<idLep.size(); nl++){
	    if(idTight[nl] == 1) continue;
	    effSF = effSF * fakeRateFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),period,typeLepSel.Data());
	    theCategory = 5;
          }
          if     (infilecatv[ifile] != 0 && goodIsTight == idTight.size()-2) effSF =  1.0 * effSF; // double fake, MC
          else if(infilecatv[ifile] != 0 && goodIsTight == idTight.size()-1) effSF = -1.0 * effSF; // single fake, MC
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-2) effSF = -1.0 * effSF; // double fake, data
          else if(infilecatv[ifile] == 0 && goodIsTight == idTight.size()-1) effSF =  1.0 * effSF; // single fake, data
        }
        else if(infilecatv[ifile] != 0 && infilecatv[ifile] != 6 && goodIsGenLep != isGenLep.size()){ // remove MC dilepton fakes from ll events
          fakeSF = 0.0;
        }
        else if(infilecatv[ifile] != 0 && goodIsGenLep == isGenLep.size()){ // MC with all good leptons
	  fakeSF = 1.0;
        }
        else if(infilecatv[ifile] == 0 || infilecatv[ifile] == 6){ // data or W+gamma with all good leptons
	  fakeSF = 1.0;
        }
	else {
	  printf("PROBLEM: %d %d %d %d %d\n",infilecatv[ifile],goodIsGenLep,(int)isGenLep.size(),goodIsTight,(int)idTight.size());
	  assert(0);
	}
      }
      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*fakeSF;
      // end event weighting
      
      int classType = kOther;
      if(theCategory == 0) classType = kData;
      if(theCategory == 3) classType = kTTBAR;
      if(theCategory == 4) classType = kTW;

      Double_t etaMin = 0; Double_t bTagMax[2] = {0,0};
      if(idJet.size() >= 2){
        etaMin = TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta());
	bTagMax[0] = (float)(*eventJets.bDiscr)[idJet[0]];
	bTagMax[1] = (float)(*eventJets.bDiscr)[idJet[1]];
        if(etaMin > TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta())) {
          etaMin =  TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta());
	  bTagMax[0] = (float)(*eventJets.bDiscr)[idJet[1]];
	  bTagMax[1] = (float)(*eventJets.bDiscr)[idJet[0]];
	}
      }

      if((bDiscrLowPt < 0.605 && idSoftNotInJet.size() == 0) && idJet.size() == 2 && bTagMax[1] < 0.605){
        int nEta = TMath::Min(etaMin,2.499)/2.5*etaBins;
        btag_central_2j_den[classType][nEta] 	     += totalWeight;
        btag_central_2j_den_error[classType][nEta]   += totalWeight*totalWeight;
        if(bTagMax[0] >= 0.605 || idSoftNotInJet.size() != 0){
          btag_central_2j_num[classType][nEta]	     += totalWeight;
          btag_central_2j_num_error[classType][nEta] += totalWeight*totalWeight;
        }
        btag_central_All_2j_den[classType]	     += totalWeight;
        btag_central_All_2j_den_error[classType]     += totalWeight*totalWeight;
        if(bTagMax[0] >= 0.605 || idSoftNotInJet.size() != 0){
          btag_central_All_2j_num[classType]         += totalWeight;
          btag_central_All_2j_num_error[classType]   += totalWeight*totalWeight;
        }
      }

      if((bDiscrLowPt < 0.605 && idSoftNotInJet.size() == 0) && idJet.size() == 2 && (float)(*eventJets.bDiscr)[idJet[1]] > 0.605 && TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()) < 2.5 && TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[1]])->Eta()) < 2.5){
        btag_highestpt_2j_den[classType][0]       += totalWeight;
        btag_highestpt_2j_den_error[classType][0] += totalWeight*totalWeight;
        if((float)(*eventJets.bDiscr)[idJet[0]] > 0.605 || idSoftNotInJet.size() != 0){
          btag_highestpt_2j_num[classType][0]       += totalWeight;
          btag_highestpt_2j_num_error[classType][0] += totalWeight*totalWeight;
        }
      }

      if((bDiscrLowPt < 0.605 && idSoftNotInJet.size() == 0) && idJet.size() == 1 && TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[0]])->Eta()) < 2.5){
        btag_highestpt_1j_den[classType][0]       += totalWeight;
        btag_highestpt_1j_den_error[classType][0] += totalWeight*totalWeight;
        if((float)(*eventJets.bDiscr)[idJet[0]] > 0.605 || idSoftNotInJet.size() != 0){
          btag_highestpt_1j_num[classType][0]       += totalWeight;
          btag_highestpt_1j_num_error[classType][0] += totalWeight*totalWeight;
        }
      }

      if(idJet.size() == 1 && (float)(*eventJets.bDiscr)[idJet[0]] > 0.605){
        btag_lowpt_1j_den[classType][0]       += totalWeight;
        btag_lowpt_1j_den_error[classType][0] += totalWeight*totalWeight;
        if(bDiscrLowPt > 0.605 || idSoftNotInJet.size() != 0){
          btag_lowpt_1j_num[classType][0]       += totalWeight;
          btag_lowpt_1j_num_error[classType][0] += totalWeight*totalWeight;
        }
      }

     if(idJet.size() == 0){
        btag_lowpt_0j_den[classType][0]       += totalWeight;
        btag_lowpt_0j_den_error[classType][0] += totalWeight*totalWeight;
        if(bDiscrLowPt > 0.605 || idSoftNotInJet.size() != 0){
          btag_lowpt_0j_num[classType][0]       += totalWeight;
          btag_lowpt_0j_num_error[classType][0] += totalWeight*totalWeight;
        }
      }
    }
  } // end of chain

  const int numberChan = 1;
  TString classLabel[numberChan] = {"all"};
  double btagSF = 1.0;
  //*******************************************************************************
  //2-Jet Bin : BTag Efficiency for central jet
  //*******************************************************************************
  printf("**********eff central jet 2-j**********\n");
  double effttMC_btag_central_2j[etaBins],effttMC_btag_central_2j_error[etaBins],effttMC_btag_central_tt_2j[etaBins],effttMC_btag_central_tt_2j_error[etaBins];
  double effttDA_btag_central_2j[etaBins],effttDA_btag_central_2j_error[etaBins],effttMC_btag_central_tw_2j[etaBins],effttMC_btag_central_tw_2j_error[etaBins];
  double TopBkgScaleFactor_2Jet_central,TopBkgScaleFactorUncertainty_2Jet_central;
  printf("channel               (data/background/top)-num             (data/background/top)-den\n");
  for(int nj=0; nj<etaBins; nj++) {
    //MC efficienciesw
    effttMC_btag_central_tt_2j[nj] = (btag_central_2j_num[1][nj])/(btag_central_2j_den[1][nj]);
    effttMC_btag_central_tw_2j[nj] = (btag_central_2j_num[2][nj])/(btag_central_2j_den[2][nj]);
    effttMC_btag_central_2j[nj]    = (btag_central_2j_num[1][nj] + btag_central_2j_num[2][nj]) / (btag_central_2j_den[1][nj] + btag_central_2j_den[2][nj]);
    
    effttMC_btag_central_tt_2j_error[nj] = sqrt((1.0-effttMC_btag_central_tt_2j[nj])*effttMC_btag_central_tt_2j[nj]/(btag_central_2j_den[1][nj])*
                                                 (btag_central_2j_den_error[1][nj])/(btag_central_2j_den[1][nj]));    
    effttMC_btag_central_tw_2j_error[nj] = sqrt((1.0-effttMC_btag_central_tw_2j[nj])*effttMC_btag_central_tw_2j[nj]/(btag_central_2j_den[2][nj])*
                                                 (btag_central_2j_den_error[2][nj])/(btag_central_2j_den[2][nj]));
    effttMC_btag_central_2j_error[nj]    = sqrt((1.0-effttMC_btag_central_2j[nj])*effttMC_btag_central_2j[nj]/(btag_central_2j_den[1][nj]+btag_central_2j_den[2][nj])*
                                                 (btag_central_2j_den_error[1][nj]+btag_central_2j_den_error[2][nj])/(btag_central_2j_den[1][nj]+btag_central_2j_den[2][nj]));
    //Data efficiencies
    if(btag_central_2j_den[3][nj] > 0){
    effttDA_btag_central_2j[nj] = (btag_central_2j_num[3][nj]-btag_central_2j_num[0][nj]-btag_central_2j_num[2][nj])/
                                  (btag_central_2j_den[3][nj]-btag_central_2j_den[0][nj]-btag_central_2j_den[2][nj]);    
    //effttDA_btag_central_2j[nj] = (btag_central_2j_num[3][nj]-btag_central_2j_num[0][nj])/
    //                             (btag_central_2j_den[3][nj]-btag_central_2j_den[0][nj]);
    effttDA_btag_central_2j[nj] = TMath::Min(TMath::Abs(effttDA_btag_central_2j[nj]),0.999);
    effttDA_btag_central_2j_error[nj] = sqrt((1-effttDA_btag_central_2j[nj])*effttDA_btag_central_2j[nj]/btag_central_2j_den[3][nj]);
    }
    else {
      effttDA_btag_central_2j[nj] = 0.999;
      effttDA_btag_central_2j_error[nj] = 0.999;
    }
    printf("details_central(%d) --> %5.0f/%7.2f/%7.2f  - %5.0f/%7.2f/%7.2f\n",nj,
           btag_central_2j_num[3][nj],btag_central_2j_num[0][nj],btag_central_2j_num[1][nj],
	   btag_central_2j_den[3][nj],btag_central_2j_den[0][nj],btag_central_2j_den[1][nj]);
  }
  printf("channel                      eff_data              eff_tt            ScaleFactor\n");
  for(int nj=0; nj<etaBins; nj++) {
    TopBkgScaleFactor_2Jet_central = effttDA_btag_central_2j[nj]/effttMC_btag_central_tt_2j[nj];
    TopBkgScaleFactorUncertainty_2Jet_central = effttDA_btag_central_2j_error[nj]/effttMC_btag_central_tt_2j[nj];
    printf("scaleFactor_central(%d) --> %6.3f +/- %6.3f | %6.3f +/- %6.3f | %6.3f +/- %6.3f\n",nj,
                           effttDA_btag_central_2j[nj],effttDA_btag_central_2j_error[nj],effttMC_btag_central_tt_2j[nj],effttMC_btag_central_tt_2j_error[nj],
			   TopBkgScaleFactor_2Jet_central,TopBkgScaleFactorUncertainty_2Jet_central);
  }

  // Overall
  double effttMC_btag_central_All_2j,effttMC_btag_central_All_2j_error,effttMC_btag_central_All_tt_2j,effttMC_btag_central_All_tt_2j_error;
  double effttDA_btag_central_All_2j,effttDA_btag_central_All_2j_error,effttMC_btag_central_All_tw_2j,effttMC_btag_central_All_tw_2j_error;
  effttMC_btag_central_All_tt_2j = (btag_central_All_2j_num[1])/(btag_central_All_2j_den[1]);
  effttMC_btag_central_All_tw_2j = (btag_central_All_2j_num[2])/(btag_central_All_2j_den[2]);
  effttMC_btag_central_All_2j    = (btag_central_All_2j_num[1] + btag_central_All_2j_num[2]) / (btag_central_All_2j_den[1] + btag_central_All_2j_den[2]);
  
  effttMC_btag_central_All_tt_2j_error = sqrt((1.0-effttMC_btag_central_All_tt_2j)*effttMC_btag_central_All_tt_2j/(btag_central_All_2j_den[1])*
  					       (btag_central_All_2j_den_error[1])/(btag_central_All_2j_den[1]));    
  effttMC_btag_central_All_tw_2j_error = sqrt((1.0-effttMC_btag_central_All_tw_2j)*effttMC_btag_central_All_tw_2j/(btag_central_All_2j_den[2])*
  					       (btag_central_All_2j_den_error[2])/(btag_central_All_2j_den[2]));
  effttMC_btag_central_All_2j_error    = sqrt((1.0-effttMC_btag_central_All_2j)*effttMC_btag_central_All_2j/(btag_central_All_2j_den[1]+btag_central_All_2j_den[2])*
  					       (btag_central_All_2j_den_error[1]+btag_central_All_2j_den_error[2])/(btag_central_All_2j_den[1]+btag_central_All_2j_den[2]));
  effttDA_btag_central_All_2j = (btag_central_All_2j_num[3]-btag_central_All_2j_num[0]-btag_central_All_2j_num[2])/
  			       (btag_central_All_2j_den[3]-btag_central_All_2j_den[0]-btag_central_All_2j_den[2]);    
  //effttDA_btag_central_All_2j = (btag_central_All_2j_num[3]-btag_central_All_2j_num[0])/
  //				 (btag_central_All_2j_den[3]-btag_central_All_2j_den[0]);    
  effttDA_btag_central_All_2j_error = sqrt((1-effttDA_btag_central_All_2j)*effttDA_btag_central_All_2j/btag_central_All_2j_den[3]);

  printf("details_central(all)         --> %5.0f/%7.2f/%7.2f  - %5.0f/%7.2f/%7.2f\n",
  	 btag_central_All_2j_num[3],btag_central_All_2j_num[0],btag_central_All_2j_num[1],
         btag_central_All_2j_den[3],btag_central_All_2j_den[0],btag_central_All_2j_den[1]);
  double TopBkgScaleFactor_2Jet_central_All = effttDA_btag_central_All_2j/effttMC_btag_central_All_tt_2j;
  double TopBkgScaleFactorUncertainty_2Jet_central_All = effttDA_btag_central_All_2j_error/effttMC_btag_central_All_tt_2j;
  printf("scaleFactor_central_All(all) --> %6.3f +/- %6.3f | %6.3f +/- %6.3f | %6.3f +/- %6.3f\n",
  			 effttDA_btag_central_All_2j,effttDA_btag_central_All_2j_error,effttMC_btag_central_All_tt_2j,effttMC_btag_central_All_tt_2j_error,
        		 TopBkgScaleFactor_2Jet_central_All,TopBkgScaleFactorUncertainty_2Jet_central_All);

  //*******************************************************************************
  //2-Jet Bin : BTag Efficiency for highest pt jet
  //*******************************************************************************
  printf("**********eff highest pt jet 2-j**********\n");
  double effttMC_btag_highestpt_2j[numberChan],effttMC_btag_highestpt_2j_error[numberChan],effttMC_btag_highestpt_tt_2j[numberChan],effttMC_btag_highestpt_tt_2j_error[numberChan];
  double effttDA_btag_highestpt_2j[numberChan],effttDA_btag_highestpt_2j_error[numberChan],effttMC_btag_highestpt_tw_2j[numberChan],effttMC_btag_highestpt_tw_2j_error[numberChan];

  for(int i=0; i<numberChan; i++) {
    //MC efficiencies
    effttMC_btag_highestpt_tt_2j[i] = (btag_highestpt_2j_num[1][i])/(btag_highestpt_2j_den[1][i]);
    effttMC_btag_highestpt_tw_2j[i] = (btag_highestpt_2j_num[2][i])/(btag_highestpt_2j_den[2][i]);
    effttMC_btag_highestpt_2j[i]    = (btag_highestpt_2j_num[1][i] + btag_highestpt_2j_num[2][i]) / (btag_highestpt_2j_den[1][i] + btag_highestpt_2j_den[2][i]);
    
    effttMC_btag_highestpt_tt_2j_error[i] = sqrt((1.0-effttMC_btag_highestpt_tt_2j[i])*effttMC_btag_highestpt_tt_2j[i]/(btag_highestpt_2j_den[1][i])*
                                                 (btag_highestpt_2j_den_error[1][i])/(btag_highestpt_2j_den[1][i]));    
    effttMC_btag_highestpt_tw_2j_error[i] = sqrt((1.0-effttMC_btag_highestpt_tw_2j[i])*effttMC_btag_highestpt_tw_2j[i]/(btag_highestpt_2j_den[2][i])*
                                                 (btag_highestpt_2j_den_error[2][i])/(btag_highestpt_2j_den[2][i]));
    effttMC_btag_highestpt_2j_error[i]    = sqrt((1.0-effttMC_btag_highestpt_2j[i])*effttMC_btag_highestpt_2j[i]/(btag_highestpt_2j_den[1][i]+btag_highestpt_2j_den[2][i])*
                                                 (btag_highestpt_2j_den_error[1][i]+btag_highestpt_2j_den_error[2][i])/(btag_highestpt_2j_den[1][i]+btag_highestpt_2j_den[2][i]));

    //Data efficiencies
    effttDA_btag_highestpt_2j[i] = (btag_highestpt_2j_num[3][i]-btag_highestpt_2j_num[0][i]-btag_highestpt_2j_num[2][i]*btagSF)/
      (btag_highestpt_2j_den[3][i]-btag_highestpt_2j_den[0][i]-btag_highestpt_2j_den[2][i]*btagSF);    
    effttDA_btag_highestpt_2j_error[i] = sqrt((1-effttDA_btag_highestpt_2j[i])*effttDA_btag_highestpt_2j[i]/btag_highestpt_2j_den[3][i]);
  }

  for(int i=0; i<numberChan; i++) printf("scaleFactor2j(%d) --> %6.3f +/- %6.3f\n",i,
             (btag_highestpt_2j_num[3][i]-btag_highestpt_2j_num[0][i])/(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]),
             sqrt(btag_highestpt_2j_num[3][i])/(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]));
  //double TopBkgScaleFactor_2Jet = (btag_highestpt_2j_num[3][4]-btag_highestpt_2j_num[0][4])/(btag_highestpt_2j_num[1][4]+btag_highestpt_2j_num[2][4]);
  //double TopBkgScaleFactorUncertainty_2Jet = sqrt(btag_highestpt_2j_num[3][4])/(btag_highestpt_2j_num[1][4]+btag_highestpt_2j_num[2][4]);

  for(int i=0; i<numberChan; i++) {
    printf("numerator  (%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i].Data(),btag_highestpt_2j_num[3][i],btag_highestpt_2j_num[0][i],(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]),btag_highestpt_2j_num[1][i],btag_highestpt_2j_num[2][i]);
    printf("denominator(%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i].Data(),btag_highestpt_2j_den[3][i],btag_highestpt_2j_den[0][i],(btag_highestpt_2j_den[1][i]+btag_highestpt_2j_den[2][i]),btag_highestpt_2j_den[1][i],btag_highestpt_2j_den[2][i]);
  }

  printf("channel         eff_tttw             eff_tt                eff_tw               eff_data                         ScaleFactor\n");
  for(int i=0; i<numberChan; i++) {
    printf("eff (%s): %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f --> %6.3f +/- %6.3f      : scaleFactor2j(%s) --> %6.3f +/- %6.3f\n",classLabel[i].Data(),
           effttMC_btag_highestpt_2j[i]   ,effttMC_btag_highestpt_2j_error[i],effttMC_btag_highestpt_tt_2j[i],effttMC_btag_highestpt_tt_2j_error[i],
           effttMC_btag_highestpt_tw_2j[i],effttMC_btag_highestpt_tw_2j_error[i],effttDA_btag_highestpt_2j[i],effttDA_btag_highestpt_2j_error[i],
           classLabel[i].Data(),(btag_highestpt_2j_num[3][i]-btag_highestpt_2j_num[0][i])/(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]),
           sqrt(btag_highestpt_2j_num[3][i])/(btag_highestpt_2j_num[1][i]+btag_highestpt_2j_num[2][i]));
  }

  printf("****************************************************************************************************************************************\n");
  //*******************************************************************************
  //1-Jet Bin : BTag Efficiency for highest pt jet
  //*******************************************************************************
  printf("**********eff highest pt jet 1-j**********\n");
  double effttMC_btag_highestpt_1j[numberChan],effttMC_btag_highestpt_1j_error[numberChan],effttMC_btag_highestpt_tt_1j[numberChan],effttMC_btag_highestpt_tt_1j_error[numberChan];
  double effttDA_btag_highestpt_1j[numberChan],effttDA_btag_highestpt_1j_error[numberChan],effttMC_btag_highestpt_tw_1j[numberChan],effttMC_btag_highestpt_tw_1j_error[numberChan];

  for(int i=0; i<numberChan; i++) {
    //MC btag efficiencies
    effttMC_btag_highestpt_tt_1j[i] = btag_highestpt_1j_num[1][i] / btag_highestpt_1j_den[1][i];
    effttMC_btag_highestpt_tw_1j[i] = btag_highestpt_1j_num[2][i] / btag_highestpt_1j_den[2][i];
    effttMC_btag_highestpt_1j[i]    = (btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i])/(btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i]);
    
    effttMC_btag_highestpt_1j_error[i]    = sqrt((1.0-effttMC_btag_highestpt_1j[i])*effttMC_btag_highestpt_1j[i]/(btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i])*
                                                 (btag_highestpt_1j_den_error[1][i]+btag_highestpt_1j_den_error[2][i])/(btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i]));
    effttMC_btag_highestpt_tt_1j_error[i] = sqrt((1.0-effttMC_btag_highestpt_tt_1j[i])*effttMC_btag_highestpt_tt_1j[i]/(btag_highestpt_1j_den[1][i])*
                                                 (btag_highestpt_1j_den_error[1][i])/(btag_highestpt_1j_den[1][i]));
    effttMC_btag_highestpt_tw_1j_error[i] = sqrt((1.0-effttMC_btag_highestpt_tw_1j[i])*effttMC_btag_highestpt_tw_1j[i]/(btag_highestpt_1j_den[2][i])*
                                                 (btag_highestpt_1j_den_error[2][i])/(btag_highestpt_1j_den[2][i]));
    
    //Data btag efficiencies
    effttDA_btag_highestpt_1j[i] = (btag_highestpt_1j_num[3][i]-btag_highestpt_1j_num[0][i]-btag_highestpt_1j_num[2][i]*btagSF)/(btag_highestpt_1j_den[3][i]-btag_highestpt_1j_den[0][i]-btag_highestpt_1j_den[2][i]*btagSF);
    effttDA_btag_highestpt_1j_error[i] = sqrt((1-effttDA_btag_highestpt_1j[i])*effttDA_btag_highestpt_1j[i]/btag_highestpt_1j_den[3][i]);    
  }

  for(int i=0; i<numberChan; i++) {
    printf("numerator  (%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i].Data(),btag_highestpt_1j_num[3][i],btag_highestpt_1j_num[0][i],(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i]),btag_highestpt_1j_num[1][i],btag_highestpt_1j_num[2][i]);
    printf("denominator(%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i].Data(),btag_highestpt_1j_den[3][i],btag_highestpt_1j_den[0][i],(btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i]),btag_highestpt_1j_den[1][i],btag_highestpt_1j_den[2][i]);
  }

  printf("channel       eff_tttw           eff_tt              eff_tw               eff_data              \n");
  for(int i=0; i<numberChan; i++) {
    printf("eff (%d): %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f --> %6.3f +/- %6.3f  \n",i,
           effttMC_btag_highestpt_1j[i]   ,effttMC_btag_highestpt_1j_error[i],effttMC_btag_highestpt_tt_1j[i],effttMC_btag_highestpt_tt_1j_error[i],
           effttMC_btag_highestpt_tw_1j[i],effttMC_btag_highestpt_tw_1j_error[i],effttDA_btag_highestpt_1j[i],effttDA_btag_highestpt_1j_error[i]);
  }
  
  // we use the combined efficiency obtained in the 2-j bin, instead of the obtained final state by final state
  double estimationMC_btag_highestpt_1j[numberChan];     
  double estimationMC_btag_highestpt_1j_err[numberChan]; 
  double estimationDA_btag_highestpt_1j[numberChan];
  double estimationDA_btag_highestpt_1j_error[numberChan]; 
  const int numberEff = 0;  
  const bool useScaleFactorEff1j = true;

  for(int i=0; i<numberChan; i++) {
    estimationMC_btag_highestpt_1j[i] = (1-effttMC_btag_highestpt_2j[numberEff])/effttMC_btag_highestpt_2j[numberEff]*(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i]);
    estimationMC_btag_highestpt_1j_err[i] = effttMC_btag_highestpt_tt_1j_error[i]/effttMC_btag_highestpt_2j[numberEff]/effttMC_btag_highestpt_2j[numberEff]*(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i]);
    if(useScaleFactorEff1j == false){
      estimationDA_btag_highestpt_1j[i] = (1-effttDA_btag_highestpt_2j[numberEff])/effttDA_btag_highestpt_2j[numberEff]*(btag_highestpt_1j_num[3][i]-btag_highestpt_1j_num[0][i]);
    } else {
      double scaleF = effttDA_btag_highestpt_2j[i]/effttMC_btag_highestpt_tt_2j[i];
      estimationDA_btag_highestpt_1j[i] = (1-effttMC_btag_highestpt_1j[i]*scaleF)/effttMC_btag_highestpt_1j[i]*scaleF*(btag_highestpt_1j_num[3][i]-btag_highestpt_1j_num[0][i]);
    }

    estimationDA_btag_highestpt_1j_error[i] = sqrt(((btag_highestpt_1j_num[3][i]-btag_highestpt_1j_num[0][i])*effttDA_btag_highestpt_2j_error[numberEff]/effttDA_btag_highestpt_2j[numberEff]/effttDA_btag_highestpt_2j[numberEff])*
                                                   ((btag_highestpt_1j_num[3][i]-btag_highestpt_1j_num[0][i])*effttDA_btag_highestpt_2j_error[numberEff]/effttDA_btag_highestpt_2j[numberEff]/effttDA_btag_highestpt_2j[numberEff])+
                                                   (1-effttDA_btag_highestpt_2j[numberEff])/effttDA_btag_highestpt_2j[numberEff]*
                                                   (1-effttDA_btag_highestpt_2j[numberEff])/effttDA_btag_highestpt_2j[numberEff]*btag_highestpt_1j_num[3][i]);    
  }

  double systMC_1j[numberChan];
  printf("Predicted ttbar+tW background for 1jet analysis (fails btag):  top background scale factor\n");
  printf("               MC(tt + tW)     Predicted from MC     | Prediction from Data   |    Scale Factor \n");
  for(int i=0; i<numberChan; i++) {
    printf("top 1-jet(%s):    %7.2f        %7.2f +/- %6.2f       | %7.2f +/- %6.2f        |    %5.3f +/- %5.3f\n",classLabel[i].Data(),
           (btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i])-(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i]),
           estimationMC_btag_highestpt_1j[i],estimationMC_btag_highestpt_1j_err[i],
           estimationDA_btag_highestpt_1j[i],estimationDA_btag_highestpt_1j_error[i],
           estimationDA_btag_highestpt_1j[i]      /((btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i])-(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i])),
           estimationDA_btag_highestpt_1j_error[i]/((btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i])-(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i])));
    systMC_1j[i] = ((btag_highestpt_1j_den[1][i]+btag_highestpt_1j_den[2][i])-(btag_highestpt_1j_num[1][i]+btag_highestpt_1j_num[2][i]))/estimationMC_btag_highestpt_1j[i];
    if(systMC_1j[i] < 1.0) systMC_1j[i] = 1.0/systMC_1j[i];
    systMC_1j[i] = 1.0-systMC_1j[i];
  }

  //*******************************************************************************
  //0-Jet Bin
  //*******************************************************************************

  printf("**********eff low pt jet 1-j**********\n");

  //Single Top Scale Factor for 0-Jet bin events which failed b-tagging
  double btagSFTW = estimationDA_btag_highestpt_1j[numberEff] /((btag_highestpt_1j_den[1][numberEff]+btag_highestpt_1j_den[2][numberEff])-(btag_highestpt_1j_num[1][numberEff]+btag_highestpt_1j_num[2][numberEff]));
  printf("btagSFTW = %f\n",btagSFTW);
  double TopBkgScaleFactor_1Jet = btagSFTW;
  double TopBkgScaleFactorUncertainty_1Jet = estimationDA_btag_highestpt_1j_error[numberEff]/((btag_highestpt_1j_den[1][numberEff]+btag_highestpt_1j_den[2][numberEff])-(btag_highestpt_1j_num[1][numberEff]+btag_highestpt_1j_num[2][numberEff]));
  TopBkgScaleFactorUncertainty_1Jet = sqrt(TopBkgScaleFactorUncertainty_1Jet*TopBkgScaleFactorUncertainty_1Jet+systMC_1j[numberEff]*systMC_1j[numberEff]);

  double effttMC_btag_lowpt_1j[numberChan],effttMC_btag_lowpt_1j_error[numberChan],effttMC_btag_lowpt_tt_1j[numberChan],effttMC_btag_lowpt_tt_1j_error[numberChan];
  double effttDA_btag_lowpt_1j[numberChan],effttDA_btag_lowpt_1j_error[numberChan],effttMC_btag_lowpt_tw_1j[numberChan],effttMC_btag_lowpt_tw_1j_error[numberChan];

  for(int i=0; i<numberChan; i++) {
    //MC btag efficiencies
    effttMC_btag_lowpt_1j[i]    = (btag_lowpt_1j_num[1][i]+btag_lowpt_1j_num[2][i])/(btag_lowpt_1j_den[1][i]+btag_lowpt_1j_den[2][i]);
    effttMC_btag_lowpt_tt_1j[i] = (btag_lowpt_1j_num[1][i])/(btag_lowpt_1j_den[1][i]                          );
    effttMC_btag_lowpt_tw_1j[i] = (btag_lowpt_1j_num[2][i])/(btag_lowpt_1j_den[2][i]);
    effttMC_btag_lowpt_1j_error[i]    = sqrt((1.0-effttMC_btag_lowpt_1j[i])*effttMC_btag_lowpt_1j[i]/(btag_lowpt_1j_den[1][i]+btag_lowpt_1j_den[2][i])*
                                             (btag_lowpt_1j_den_error[1][i]+btag_lowpt_1j_den_error[2][i])/(btag_lowpt_1j_den[1][i]+btag_lowpt_1j_den[2][i]));
    effttMC_btag_lowpt_tt_1j_error[i] = sqrt((1.0-effttMC_btag_lowpt_tt_1j[i])*effttMC_btag_lowpt_tt_1j[i]/(btag_lowpt_1j_den[1][i])*
                                             (btag_lowpt_1j_den_error[1][i])/(btag_lowpt_1j_den[1][i]));
    effttMC_btag_lowpt_tw_1j_error[i] = sqrt((1.0-effttMC_btag_lowpt_tw_1j[i])*effttMC_btag_lowpt_tw_1j[i]/(btag_lowpt_1j_den[2][i])*
                                             (btag_lowpt_1j_den_error[2][i])/(btag_lowpt_1j_den[2][i]));

    //Data btag efficiencies
    effttDA_btag_lowpt_1j[i]       = (btag_lowpt_1j_num[3][i]-btag_lowpt_1j_num[0][i]-btag_lowpt_1j_num[2][i]*btagSF*btagSFTW)/ (btag_lowpt_1j_den[3][i]-btag_lowpt_1j_den[0][i]-btag_lowpt_1j_den[2][i]*btagSF*btagSFTW);
    effttDA_btag_lowpt_1j_error[i] = sqrt((1-effttDA_btag_lowpt_1j[i])*effttDA_btag_lowpt_1j[i]/btag_lowpt_1j_den[3][i]);
  }

  for(int i=0; i<numberChan; i++) {
    printf("numerator(%s)   --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i].Data(),btag_lowpt_1j_num[3][i],btag_lowpt_1j_num[0][i],(btag_lowpt_1j_num[1][i]+btag_lowpt_1j_num[2][i]),btag_lowpt_1j_num[1][i],btag_lowpt_1j_num[2][i]);
    printf("denominator(%s) --> data: %4.0f, background: %7.2f, tt+tw: %7.2f, tt: %7.2f, tw: %7.2f\n",classLabel[i].Data(),btag_lowpt_1j_den[3][i],btag_lowpt_1j_den[0][i],(btag_lowpt_1j_den[1][i]+btag_lowpt_1j_den[2][i]),btag_lowpt_1j_den[1][i],btag_lowpt_1j_den[2][i]);
  }

  printf("channel       eff_tttw         eff_tt            eff_tw           eff_data              \n");
  for(int i=0; i<numberChan; i++) {
    printf("eff (%d): %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f --> %6.3f +/- %6.3f\n",i,
           effttMC_btag_lowpt_1j[i]   ,effttMC_btag_lowpt_1j_error[i],effttMC_btag_lowpt_tt_1j[i],effttMC_btag_lowpt_tt_1j_error[i],
           effttMC_btag_lowpt_tw_1j[i],effttMC_btag_lowpt_tw_1j_error[i],effttDA_btag_lowpt_1j[i],effttDA_btag_lowpt_1j_error[i]);
  }

  printf("**********************************************************\n");

  //*******************************************************************************
  //Closure Test 0-Jet Bin
  //*******************************************************************************

  printf("**********eff low pt jet 0-j**********\n");
  double ftw_b[numberChan]; 
  for(int i=0; i<numberChan; i++) {
    ftw_b[i] = effttMC_btag_lowpt_tw_1j[i];
    printf("ftw_b(%d) = %5.3f \n",i,ftw_b[i]);
  }

  double N_top_expected_0j[numberChan]; 
  double fttbar[numberChan],sigma_ftop[numberChan];

  for(int i=0; i<numberChan; i++) {
    N_top_expected_0j[i] = (btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i])-(btag_lowpt_0j_num[1][i]+btag_lowpt_0j_num[2][i]);
    fttbar[i] = (btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]*ftw_b[i])/(btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]);
    sigma_ftop[i] = 0.17*btag_lowpt_0j_den[1][i]*(1.0-ftw_b[i])/(btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]);
  }

  double effMC_btag_lowpt_tt_0j_expected[numberChan];  
  double effMC_btag_lowpt_tw_0j_expected[numberChan];  
  double effMC_btag_lowpt_tt_0j[numberChan];           
  double effMC_btag_lowpt_tw_0j[numberChan];           

  for(int i=0; i<numberChan; i++) {
    effMC_btag_lowpt_tt_0j_expected[i] = btag_lowpt_0j_num[1][i]/btag_lowpt_0j_den[1][i];
    effMC_btag_lowpt_tw_0j_expected[i] = btag_lowpt_0j_num[2][i]/btag_lowpt_0j_den[2][i];
    effMC_btag_lowpt_tt_0j[i]          = 1-(1-effttMC_btag_lowpt_tt_1j[i])*(1-effttMC_btag_lowpt_tt_1j[i]);
    effMC_btag_lowpt_tw_0j[i]          = effttMC_btag_lowpt_tt_1j[i];    
  }

  printf("channel        ttbar MC ( 0Jet / Extrapolated from 1Jet )            tW MC ( 0Jet / Extrapolated from 1Jet) \n");
  for(int i=0; i<numberChan; i++) { 
    printf("(%s),                  %5.3f/%5.3f                                              %5.3f/%5.3f\n",
           classLabel[i].Data(),effMC_btag_lowpt_tt_0j_expected[i],effMC_btag_lowpt_tt_0j[i],effMC_btag_lowpt_tw_0j_expected[i],effMC_btag_lowpt_tw_0j[i]);
  }

  // begin get closure test closing!!!!!!!!!!!!!!!
  // double ftw_2b[numberChan];
  // for(int i=0; i<5; i++) {
  //   ftw_2b[i] = (effMC_btag_lowpt_tw_0j_expected[i]-effttMC_btag_lowpt_tt_1j[i])/(effMC_btag_lowpt_tt_0j[i]-effttMC_btag_lowpt_tt_1j[i]);
  //   printf("ftw_2b(%d) = %5.3f - ",i,ftw_2b[i]);
  //   fttbar[i] = (btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]*ftw_2b[i])/(btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]);
  //   printf("fttbar(%d) = %5.3f | ",i,fttbar[i]);
  // }
  // printf("\n");
  // end get closure test closing!!!!!!!!!!!!!!!

  double effMC_btag_lowpt_0j[numberChan]; 
  double effDA_btag_lowpt_0j[numberChan]; 
  double effMC_btag_lowpt_0j_error[numberChan];
  double effDA_btag_lowpt_0j_error[numberChan];

  for(int i=0; i<numberChan; i++) {
    effMC_btag_lowpt_0j[i] = fttbar[i]*effMC_btag_lowpt_tt_0j[i]+(1-fttbar[i])*effMC_btag_lowpt_tw_0j[i];
    effDA_btag_lowpt_0j[i] = fttbar[i]*(1-(1-effttDA_btag_lowpt_1j[i])*(1-effttDA_btag_lowpt_1j[i]))+(1-fttbar[i])*effttDA_btag_lowpt_1j[i];

    effMC_btag_lowpt_0j_error[i] = sqrt(effttMC_btag_lowpt_1j_error[i]*(fttbar[i]-2*effttMC_btag_lowpt_1j[i]*fttbar[i]+1)*
					effttMC_btag_lowpt_1j_error[i]*(fttbar[i]-2*effttMC_btag_lowpt_1j[i]*fttbar[i]+1));
    effDA_btag_lowpt_0j_error[i] = sqrt(sigma_ftop[i]*(effttDA_btag_lowpt_1j[i]-effttDA_btag_lowpt_1j[i]*effttDA_btag_lowpt_1j[i])*
                                        sigma_ftop[i]*(effttDA_btag_lowpt_1j[i]-effttDA_btag_lowpt_1j[i]*effttDA_btag_lowpt_1j[i])+
					effttDA_btag_lowpt_1j_error[i]*(fttbar[i]-2*effttDA_btag_lowpt_1j[i]*fttbar[i]+1)*
					effttDA_btag_lowpt_1j_error[i]*(fttbar[i]-2*effttDA_btag_lowpt_1j[i]*fttbar[i]+1));
  }

  printf("top tagging efficiency\n");
  printf("Channel    fttbar        Eff toptag(MC)    Eff toptag(MC extrapolated)       Eff toptab Data \n");
  for(int i=0; i<numberChan; i++) {
    printf("(%s)       %5.3f +/- %5.3f,        : %6.3f                 %6.3f +/- %6.3f             %6.3f +/- %6.3f\n",
           classLabel[i].Data(),fttbar[i],sigma_ftop[i],
           (btag_lowpt_0j_num[1][i]+btag_lowpt_0j_num[2][i])/(btag_lowpt_0j_den[1][i]+btag_lowpt_0j_den[2][i]),
           effMC_btag_lowpt_0j[i],effMC_btag_lowpt_0j_error[i],
           effDA_btag_lowpt_0j[i],effDA_btag_lowpt_0j_error[i]);
  }

  double sigma_0f_bck = 0.20;
  double estimationMC_btag_lowpt_0j[numberChan]; 
  double estimationDA_btag_lowpt_0j[numberChan]; 
  double estimationMC_btag_lowpt_0j_error[numberChan]; 
  double estimationDA_btag_lowpt_0j_error[numberChan]; 
  // we use the combined efficiency obtained in the 1-j bin, instead of the obtained final state by final state
  for(int i=0; i<numberChan; i++) {
    estimationMC_btag_lowpt_0j[i] = (1-effMC_btag_lowpt_0j[numberEff])/effMC_btag_lowpt_0j[numberEff]*(btag_lowpt_0j_num[1][i]+btag_lowpt_0j_num[2][i]);
    estimationDA_btag_lowpt_0j[i] = (1-effDA_btag_lowpt_0j[numberEff])/effDA_btag_lowpt_0j[numberEff]*(btag_lowpt_0j_num[3][i]-btag_lowpt_0j_num[0][i]);
    estimationMC_btag_lowpt_0j_error[i] = sqrt(((btag_lowpt_0j_num[1][i]+btag_lowpt_0j_num[2][i])*effMC_btag_lowpt_0j_error[numberEff]/effMC_btag_lowpt_0j[numberEff]/effMC_btag_lowpt_0j[numberEff])*
                                               ((btag_lowpt_0j_num[1][i]+btag_lowpt_0j_num[2][i])*effMC_btag_lowpt_0j_error[numberEff]/effMC_btag_lowpt_0j[numberEff]/effMC_btag_lowpt_0j[numberEff]));
    
    estimationDA_btag_lowpt_0j_error[i] = sqrt(((btag_lowpt_0j_num[3][i]-btag_lowpt_0j_num[0][i])*effDA_btag_lowpt_0j_error[numberEff]/effDA_btag_lowpt_0j[numberEff]/effDA_btag_lowpt_0j[numberEff])*
                                               ((btag_lowpt_0j_num[3][i]-btag_lowpt_0j_num[0][i])*effDA_btag_lowpt_0j_error[numberEff]/effDA_btag_lowpt_0j[numberEff]/effDA_btag_lowpt_0j[numberEff])+
                                               (1-effDA_btag_lowpt_0j[i])/effDA_btag_lowpt_0j[i]*
                                               (1-effDA_btag_lowpt_0j[i])/effDA_btag_lowpt_0j[i]*btag_lowpt_0j_num[3][i]+
                                               TMath::Power(sigma_0f_bck*btag_lowpt_0j_num[0][i]*(1-effDA_btag_lowpt_0j[numberEff])/effDA_btag_lowpt_0j[numberEff],2));
  }

  printf("0-Jet Top Background\n");
  printf("Channel    top-tagged region    bkgs top-tagged          MC tt+tW              MC top            Data top   \n");
  printf("              Event Count       region (non-top)    non-top-tagged count     estimation         estimation  \n");
  for(int i=0; i<numberChan; i++) {
    printf("(%s)              %3d               %6.3f                %6.3f             %6.3f +/- %6.3f    %6.3f +/- %6.3f\n",
           classLabel[i].Data(),
           (int)btag_lowpt_0j_num[3][i],btag_lowpt_0j_num[0][i],N_top_expected_0j[i],
           estimationMC_btag_lowpt_0j[i],estimationMC_btag_lowpt_0j_error[i],
           estimationDA_btag_lowpt_0j[i],estimationDA_btag_lowpt_0j_error[i]);
  }

  // Writing results in latex format
  printf("*****************************************\n");
  printf("\\begin{table}\n");
  printf("\\begin{center}\n");
  printf("{\\tiny\n");
  printf("\\begin{tabular}{|c|c|c|c|}\n");
  printf("\\hline\n");
  printf(" Sample & 0-jet & 1-jet & 2-jet \\\\\n");
  printf("\\hline\n");
  printf("estimated top events in simulation  & %5.3f $\\pm$ %5.1f &  %5.1f $\\pm$ %5.1f & %5.1f $\\pm$ %5.1f \\\\\n",N_top_expected_0j[numberEff],estimationMC_btag_lowpt_0j_error[numberEff],
  (btag_highestpt_1j_den[1][numberEff]+btag_highestpt_1j_den[2][numberEff])-(btag_highestpt_1j_num[1][numberEff]+btag_highestpt_1j_num[2][numberEff]),estimationMC_btag_highestpt_1j_err[numberEff],
  (btag_central_All_2j_den[1]+btag_central_All_2j_den[2])-(btag_central_All_2j_num[1]+btag_central_All_2j_num[2]),sqrt(btag_central_All_2j_den_error[1]+btag_central_All_2j_den_error[2]));
  printf("tagging efficiency                  & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f & %4.1f $\\pm$ %4.1f \\\\ \n",effDA_btag_lowpt_0j[numberEff]*100,effDA_btag_lowpt_0j_error[numberEff]*100,
  effttDA_btag_highestpt_1j[numberEff]*100,effttDA_btag_highestpt_1j_error[numberEff]*100,effttDA_btag_central_All_2j*100,effttDA_btag_central_All_2j_error*100);
  printf("data events in control region       & %4d & %4d & - \\\\ \n",(int)btag_lowpt_0j_num[3][numberEff],(int)btag_highestpt_1j_num[3][numberEff]);
  printf("background events in control region & %5.1f $\\pm$ %5.1f &  %5.1f $\\pm$ %5.1f & - \\\\ \n",btag_lowpt_0j_num[0][numberEff],btag_lowpt_0j_num[0][numberEff]*0.15,
  btag_highestpt_1j_num[0][numberEff],btag_highestpt_1j_num[0][numberEff]*0.20);
  printf("top estimation in data              &  %5.1f $\\pm$ %5.1f &  %5.1f $\\pm$ %5.1f & -\\\\\n",estimationDA_btag_lowpt_0j[numberEff],estimationDA_btag_lowpt_0j_error[numberEff],
  estimationDA_btag_highestpt_1j[numberEff],estimationDA_btag_highestpt_1j_error[numberEff]);
  printf("data/simulation scale factor        &  %5.2f $\\pm$ %5.2f &  %5.2f $\\pm$ %5.2f & %5.2f $\\pm$ %5.2f \\\\\n",
  estimationDA_btag_lowpt_0j[numberEff] / N_top_expected_0j[numberEff],  estimationDA_btag_lowpt_0j_error[numberEff] / N_top_expected_0j[numberEff],
  TopBkgScaleFactor_1Jet,TopBkgScaleFactorUncertainty_1Jet,TopBkgScaleFactor_2Jet_central_All,TopBkgScaleFactorUncertainty_2Jet_central_All);
  printf("\\hline\n");
  printf("\\end{tabular}\n");
  printf("}\n");
  printf("\\end{center}\n");
  printf("\\end{table}\n");
  printf("*****************************************\n");

  // Writing scale factors in *.h files
  double TopBkgScaleFactor_0Jet = estimationDA_btag_lowpt_0j[numberEff] / N_top_expected_0j[numberEff];
  double TopBkgScaleFactorUncertainty_0Jet = estimationDA_btag_lowpt_0j_error[numberEff] / N_top_expected_0j[numberEff];
}
