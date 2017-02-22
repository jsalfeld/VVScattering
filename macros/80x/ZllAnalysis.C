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

#include "MitAnalysisRunII/macros/80x/factors.h"

#include "MitAnalysisRunII/macros/LeptonScaleLookup.h"

int whichSkim = 0;
Int_t period = 1;

void ZllAnalysis(TString typeLepSel = "default"){

  TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_80x/";
  TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_80x/";
  Double_t lumi = 35.9;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  if     (period==1){
  puPath = "MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root";
  infilenamev.push_back(Form("%sdata_Run2016B.root",filesPathDA.Data())); infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_Run2016C.root",filesPathDA.Data())); infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_Run2016D.root",filesPathDA.Data())); infilecatv.push_back(0);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(1);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1+AODSIM.root",filesPathMC.Data()));      infilecatv.push_back(1);
  }
  else {assert(0);}

  if(infilenamev.size() != infilecatv.size()) assert(0);

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);

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
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF->SetDirectory(0);
  delete fElSF;

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

  const int ptBins = 5;
  const int etaBins = 2;

  double yields[2][ptBins][etaBins][ptBins][etaBins][3],yieldsE[2][ptBins][etaBins][ptBins][etaBins][3];
  for (int ch=0; ch<2; ch++){
    for (int pt0=0; pt0<ptBins; pt0++){
      for (int eta0=0; eta0<etaBins; eta0++){
        for (int pt1=0; pt1<ptBins; pt1++){
          for (int eta1=0; eta1<etaBins; eta1++){
            for (int nlep=0; nlep<3; nlep++){
              yields[ch][pt0][eta0][pt1][eta1][nlep] = 0; yieldsE[ch][pt0][eta0][pt1][eta1][nlep] = 0;
            }
          }
        }
      }
    }
  }
  double yieldsTotal[2][3],yieldsTotal20[2][3];
  for (int ch=0; ch<2; ch++){
    for (int nlep=0; nlep<3; nlep++){
      yieldsTotal[ch][nlep] = 0;
      yieldsTotal20[ch][nlep] = 0;
    }
  }
  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    //TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");
    TTree *the_SelBit_tree= (TTree*)the_input_file.FindObjectAny("SelBit_tree");

    BareMonteCarlo eventMonteCarlo;
    eventMonteCarlo.setBranchAddresses(the_input_tree);

    BareEvent eventEvent;
    eventEvent.setBranchAddresses(the_input_tree);

    BareLeptons eventLeptons;
    eventLeptons.setBranchAddresses(the_input_tree);

    BareTrigger eventTrigger;
    eventTrigger.setBranchAddresses(the_input_tree);

    TNamed *triggerNames = (TNamed*)the_input_file.FindObjectAny("triggerNames");
    char **tokens;
    size_t numtokens;
    tokens = strsplit(triggerNames->GetTitle(), ",", &numtokens);
    if(infilecatv[ifile] == 0){
      for (int i = 0; i < (int)numtokens; i++) {
        printf("triggerNames(%2d): \"%s\"\n",(int)i,tokens[i]);
      }
    }

    unsigned int selBit_= 0;
    the_SelBit_tree->SetBranchAddress("selBit", &selBit_);
    int MAX = the_input_tree->GetEntries();
    //if(infilecatv[ifile] == 1) MAX = MAX/10.;
    for (int i=0; i<MAX; ++i) {
      the_SelBit_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());
      if((selBit_ & 0x1<<whichSkim) == 0) continue;

      the_input_tree->GetEntry(i);

      Bool_t passFilter[6] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
      if(eventLeptons.p4->GetEntriesFast() >= 2 &&
     	 ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 20 && 
     	 ((TLorentzVector*)(*eventLeptons.p4)[1])->Pt() > 10) passFilter[0] = kTRUE;
      if(infilecatv[ifile] != 999) {
        for (int nt = 0; nt <(int)numtokens; nt++) {
          if((*eventTrigger.triggerFired)[nt] == 0) continue;
          if((strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
             (strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
             (strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*") == 0) ||
             (strcmp(tokens[nt],"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*")  == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*") 	              == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*")	      == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*") 	      == 0) ||
             (strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*")	      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoMu20_v*") 				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoTkMu20_v*") 				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoMu22_v*") 				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoTkMu22_v*") 				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoMu24_v*")				      == 0) ||
             (strcmp(tokens[nt],"HLT_IsoTkMu24_v*")				      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")	      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*")	      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele25_eta2p1_WPTight_Gsf_v*")                    == 0) ||
             (strcmp(tokens[nt],"HLT_Ele27_eta2p1_WPLoose_Gsf_v*")                    == 0) ||
             (strcmp(tokens[nt],"HLT_Ele27_WPTight_Gsf_v*")			      == 0) ||
             (strcmp(tokens[nt],"HLT_Ele35_WPLoose_Gsf_v*")			      == 0)
             ) passFilter[1] = kTRUE;
	}
      } else { passFilter[1] = kTRUE;}

      //if(infilecatv[ifile] != 0) passFilter[1] = kTRUE; // do not apply trigger filters to MC
      if(passFilter[0] == kFALSE) continue;
      if(passFilter[1] == kFALSE) continue;

      vector<int> idLep; vector<int> idTight;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
      }
      if(idLep.size()==2) passFilter[2] = kTRUE;
      if(passFilter[2] == kFALSE) continue;

      if(idTight[0]==1&&idTight[1]==1) passFilter[3] = kTRUE;
      if(passFilter[3] == kFALSE) continue;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() <= 20 ||
         ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() <= 10) continue;

      passFilter[4] = ((int)(*eventLeptons.pdgId)[idLep[0]]*(int)(*eventLeptons.pdgId)[idLep[1]] < 0);
      if(passFilter[4] == kFALSE) continue;

      TLorentzVector dilep(( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) )); 

      if(TMath::Abs(dilep.M()-91.1876)<15.0) passFilter[5] = kTRUE;  	    
      if(passFilter[5] == kFALSE) continue;

      int iPt[2] = {-1, -1};
      if     (((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 15) iPt[0] = 0;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 20) iPt[0] = 1;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 25) iPt[0] = 2;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 30) iPt[0] = 3;
      else                                                                iPt[0] = 4;
      if     (((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() < 15) iPt[1] = 0;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() < 20) iPt[1] = 1;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() < 25) iPt[1] = 2;
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() < 30) iPt[1] = 3;
      else                                                                iPt[1] = 4;

      int iEta[2] = {-1, -1};
      if     (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) < 1.5) iEta[0] = 0;
      else                                                                              iEta[0] = 1;
      if     (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Eta()) < 1.5) iEta[1] = 0;
      else                                                                              iEta[1] = 1;

      Int_t typeLep = 2;
      if     (TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==13&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==13) typeLep = 0;
      else if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]])==11&&TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]])==11) typeLep = 1;

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
      double puWeight = 1.0; if(infilecatv[ifile] != 0) puWeight = nPUScaleFactor(fhDPU, (double)eventMonteCarlo.puTrueInt);
      double effSF = 1.0;

      if(infilecatv[ifile] != 0){
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          effSF = effSF * effhDScaleFactor(((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Pt(),
	        ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->Eta(),TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]),
	  	typeLepSel.Data(),fhDMuMediumSF,fhDElMediumSF,fhDElTightSF,fhDmutrksfptg10,fhDeltrksf,eventVertex.npv,true,fhDMuIsoSF,true);
        }
        effSF=1; // SF efficiency == 1!
      }

      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
      double totalWeight = mcWeight*theLumi*puWeight*effSF*trigEff;
   
      yields [infilecatv[ifile]][iPt[0]][iEta[0]][iPt[1]][iEta[1]][typeLep] =  yields[infilecatv[ifile]][iPt[0]][iEta[0]][iPt[1]][iEta[1]][typeLep] + totalWeight;

      yieldsE[infilecatv[ifile]][iPt[0]][iEta[0]][iPt[1]][iEta[1]][typeLep] = yieldsE[infilecatv[ifile]][iPt[0]][iEta[0]][iPt[1]][iEta[1]][typeLep] + totalWeight * totalWeight;

      yieldsTotal[infilecatv[ifile]][typeLep]   = yieldsTotal  [infilecatv[ifile]][typeLep] + totalWeight;
      if(((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 20)
      yieldsTotal20[infilecatv[ifile]][typeLep] = yieldsTotal20[infilecatv[ifile]][typeLep] + totalWeight;
    }
  } // end of chain

  double kFactor[2] = {1., 1.};
  for (int ch=0; ch<2; ch++){
    kFactor[ch] = (yieldsTotal[ch][0]/yieldsTotal[ch][1]+yieldsTotal[ch][1]/yieldsTotal[ch][0])*0.5;
  }
  printf("kFactors data/MC: %f %f\n",kFactor[0],kFactor[1]);

  for (int pt0=0; pt0<ptBins; pt0++){
    for (int eta0=0; eta0<etaBins; eta0++){
      for (int pt1=0; pt1<ptBins; pt1++){
        for (int eta1=0; eta1<etaBins; eta1++){
	  if(pt0>=pt1)
          printf("yield %1d %1d %1d %1d: %6.1f %6.1f %6.1f - %6.1f %6.1f %6.1f\n",pt0,eta0,pt1,eta1,yields[0][pt0][eta0][pt1][eta1][0],yields[0][pt0][eta0][pt1][eta1][1],yields[0][pt0][eta0][pt1][eta1][2],yields[1][pt0][eta0][pt1][eta1][0],yields[1][pt0][eta0][pt1][eta1][1],yields[1][pt0][eta0][pt1][eta1][2]);
        }         
      }
    }
  }

  printf("yieldTotal:   %8.1f %8.1f %8.1f - %8.1f %8.1f %8.1f\n",yieldsTotal  [0][0],yieldsTotal  [0][1],yieldsTotal  [0][2],yieldsTotal  [1][0],yieldsTotal  [1][1],yieldsTotal  [1][2]);
  printf("yieldTotal20: %8.1f %8.1f %8.1f - %8.1f %8.1f %8.1f\n",yieldsTotal20[0][0],yieldsTotal20[0][1],yieldsTotal20[0][2],yieldsTotal20[1][0],yieldsTotal20[1][1],yieldsTotal20[1][2]);

  for (int pt0=0; pt0<ptBins; pt0++){
    for (int eta0=0; eta0<etaBins; eta0++){
      for (int pt1=0; pt1<ptBins; pt1++){
        for (int eta1=0; eta1<etaBins; eta1++){
	  if(pt0>=pt1)
          printf("yield %1d %1d %1d %1d: %6.4f +/- %6.4f - %6.4f +/- %6.4f\n",pt0,eta0,pt1,eta1,(yields[0][pt0][eta0][pt1][eta1][0]-yields[0][pt0][eta0][pt1][eta1][2]*kFactor[0])/(yields[1][pt0][eta0][pt1][eta1][0]-yields[1][pt0][eta0][pt1][eta1][2]*kFactor[1]),
                                                                                                                              sqrt(yieldsE[0][pt0][eta0][pt1][eta1][0])*kFactor[0]/(yields[1][pt0][eta0][pt1][eta1][0]-yields[1][pt0][eta0][pt1][eta1][2]*kFactor[1]),
                                                                                                (yields[0][pt0][eta0][pt1][eta1][1]-yields[0][pt0][eta0][pt1][eta1][2]*kFactor[0])/(yields[1][pt0][eta0][pt1][eta1][1]-yields[1][pt0][eta0][pt1][eta1][2]*kFactor[1]),
                                                                                                                              sqrt(yieldsE[0][pt0][eta0][pt1][eta1][0])*kFactor[0]/(yields[1][pt0][eta0][pt1][eta1][1]-yields[1][pt0][eta0][pt1][eta1][2]*kFactor[1]));
        }
      }
    }
  }
}
