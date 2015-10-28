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

#include "MitAnalysisRunII/macros/factors.h"

void triggerAnalysis(Int_t period = 1){

  TString filesPath  = "/scratch5/ceballos/ntuples_weights/";
  Double_t lumi = 0.0715;
  if(period == 1) lumi = 0.5947;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  if     (period==0){
  puPath = "/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src/MitAnalysisRunII/data/puWeights_13TeV_50ns.root";
  infilenamev.push_back(Form("%sdata_AOD_50ns.root",filesPath.Data()));														  infilecatv.push_back(0);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));  	  infilecatv.push_back(1);
  assert(0);
  }
  else if(period==1){
  puPath = "/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src/MitAnalysisRunII/data/puWeights_13TeV_25ns.root";
  infilenamev.push_back(Form("%sdata_AOD_25ns.root",filesPath.Data()));														  infilecatv.push_back(0);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3+AODSIM.root",filesPath.Data()));  	          infilecatv.push_back(1);
  }
  else {assert(0);}

  if(infilenamev.size() != infilecatv.size()) assert(0);

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
 
  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  const int ptBins = 5;
  const int etaBins = 3;

  double yieldsP[2][ptBins][etaBins][2],yieldsF[2][ptBins][etaBins][2],eff[2][ptBins][etaBins][2];
  for (int ch=0; ch<2; ch++){
    for (int pt0=0; pt0<ptBins; pt0++){
      for (int eta0=0; eta0<etaBins; eta0++){
        for (int nlep=0; nlep<2; nlep++){
          yieldsP[ch][pt0][eta0][nlep] = 0; yieldsF[ch][pt0][eta0][nlep] = 0; eff[ch][pt0][eta0][nlep] = 0;
        }
      }
    }
  }
  double yieldsTotal[2][2];
  for (int ch=0; ch<2; ch++){
    for (int nlep=0; nlep<2; nlep++){
      yieldsTotal[ch][nlep] = 0;
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
    if(infilecatv[ifile] != 10){
      for (int i = 0; i < (int)numtokens; i++) {
        printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
      }
    }

    int MAX = the_input_tree->GetEntries();
    if(infilecatv[ifile] == 1) MAX = MAX/10.;
    for (int i=0; i<MAX; ++i) {
      the_input_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      Bool_t passFilter[6] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;

      vector<int> idTrigger;
      for (int iL = 0; iL != eventLeptons.p4->GetEntriesFast(); ++iL) idTrigger.push_back(-1);
      for (int nt = 0; nt <(int)numtokens; nt++) {
	if((strcmp(tokens[nt],"HLT_IsoMu27_v*")                                     == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_IsoMu20_v*")				            == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_IsoTkMu20_v*")				    == 0 && (*eventTrigger.triggerFired)[nt] == 1)
           ) {
	  passFilter[1] = kTRUE;
          for (int iL = 0; iL != eventLeptons.p4->GetEntriesFast(); ++iL) {
            if (((*eventTrigger.triggerLeps)[iL] & (1 << nt)) != 0) idTrigger[iL] = 0;
          }
        }
        else
	if((strcmp(tokens[nt],"HLT_Ele27_eta2p1_WP75_Gsf_v*")                       == 0 && (*eventTrigger.triggerFired)[nt] == 1) ||
           (strcmp(tokens[nt],"HLT_Ele27_eta2p1_WPLoose_Gsf_v*")                    == 0 && (*eventTrigger.triggerFired)[nt] == 1)
           ) {
	  passFilter[1] = kTRUE;
          for (int iL = 0; iL != eventLeptons.p4->GetEntriesFast(); ++iL) {
            if (((*eventTrigger.triggerLeps)[iL] & (1 << nt)) != 0) idTrigger[iL] = 1;
          }
        }
      }

      if(passFilter[0] == kFALSE) continue;
      if(passFilter[1] == kFALSE) continue;

      vector<int> idLep; vector<int> idTight;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut("medium",TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
      }
      if(idLep.size()==2) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;

      if(idTight[0]==1&&idTight[1]==1) passFilter[3] = kTRUE;
      if(passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 30 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 30) continue;

      passFilter[4] = ((int)(*eventLeptons.pdgId)[idLep[0]]*(int)(*eventLeptons.pdgId)[idLep[1]] < 0);
      if(passFilter[4] == kFALSE) continue;

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 

      if(TMath::Abs(dilep.M()-91.1876)<15.0) passFilter[5] = kTRUE;  	    
      if(passFilter[5] == kFALSE) continue;

      Int_t typeLep = 2;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) typeLep = 0;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) typeLep = 1;

      if(typeLep == 2) continue;

      double theLumi = lumi; if(infilecatv[ifile] == 0) theLumi = 1.0;
      double puWeight = nPUScaleFactor(fhDPU, (double)eventVertex.npv); if(infilecatv[ifile] == 0) puWeight = 1.0;
      double effSF = 1.0;

      if(infilecatv[ifile] != 0){
        effSF = effScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]),period)*
                effScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt(),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]),period);
      }

      double totalWeight = eventMonteCarlo.mcWeight*theLumi*puWeight*effSF;

      for(unsigned int nlep0=0; nlep0<idLep.size(); nlep0++) {
        for(unsigned int nlep1=0; nlep1<idLep.size(); nlep1++) {
	  if(idTrigger[idLep[nlep0]] == typeLep) {
            int iPt = -1;
            if     (((TLorentzVector*)(*eventLeptons.p4)[idLep[nlep1]])->Pt() < 35) iPt = 0;
            else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[nlep1]])->Pt() < 40) iPt = 1;
            else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[nlep1]])->Pt() < 45) iPt = 2;
            else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[nlep1]])->Pt() < 50) iPt = 3;
            else                                                                    iPt = 4;

            int iEta = -1;
            if     (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nlep1]])->Eta()) < 1.0) iEta = 0;
            else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[nlep1]])->Eta()) < 1.5) iEta = 1;
            else                                                                                  iEta = 2;
	  
	    if(idTrigger[idLep[nlep1]] == typeLep) yieldsP[infilecatv[ifile]][iPt][iEta][typeLep] = yieldsP[infilecatv[ifile]][iPt][iEta][typeLep] + totalWeight;
	    else                                   yieldsF[infilecatv[ifile]][iPt][iEta][typeLep] = yieldsF[infilecatv[ifile]][iPt][iEta][typeLep] + totalWeight;
	     
	  }
        }
       }

      yieldsTotal[infilecatv[ifile]][typeLep] = yieldsTotal[infilecatv[ifile]][typeLep] + totalWeight;
    }
  } // end of chain

  double kFactor[2] = {1., 1.};
  for (int ch=0; ch<2; ch++){
    kFactor[ch] = (yieldsTotal[ch][0]/yieldsTotal[ch][1]+yieldsTotal[ch][1]/yieldsTotal[ch][0])*0.5;
  }
  printf("kFactors data/MC: %f %f\n",kFactor[0],kFactor[1]);

  for (int pt0=0; pt0<ptBins; pt0++){
    for (int eta0=0; eta0<etaBins; eta0++){
      for (int nl=0; nl<2; nl++){
        eff[0][pt0][eta0][nl] = yieldsP[0][pt0][eta0][nl] / (yieldsP[0][pt0][eta0][nl]+yieldsF[0][pt0][eta0][nl]);
        eff[1][pt0][eta0][nl] = yieldsP[1][pt0][eta0][nl] / (yieldsP[1][pt0][eta0][nl]+yieldsF[1][pt0][eta0][nl]);
      }
    }
  }
  for (int pt0=0; pt0<ptBins; pt0++){
    for (int eta0=0; eta0<etaBins; eta0++){
      printf("yield %1d %1d: %6.3f %6.3f - %6.3f %6.3f\n",pt0,eta0,eff[0][pt0][eta0][0],eff[0][pt0][eta0][1],
                                                                   eff[1][pt0][eta0][0],eff[1][pt0][eta0][1]);
    }
  }

  printf("yieldTotal:   %8.1f %8.1f - %8.1f %8.1f\n",yieldsTotal[0][0],yieldsTotal[0][1],yieldsTotal[1][0],yieldsTotal[1][1]);

}
