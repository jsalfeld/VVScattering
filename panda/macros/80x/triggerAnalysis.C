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
#include "TLorentzVector.h"

#include "VVScattering/panda/macros/80x/pandaFlat.C"

enum TriggerBits {
    kMETTrig	   =(1<<0),
    kSinglePhoTrig =(1<<1),
    kMuEGTrig	   =(1<<2),
    kMuMuTrig	   =(1<<3),
    kMuTrig	   =(1<<4),
    kEGEGTrig	   =(1<<5),
    kEGTrig	   =(1<<6)
};

enum SelectionBit {
 kLoose   =(1<<0),
 kFake    =(1<<1),
 kMedium  =(1<<2),
 kTight   =(1<<3)
};

const double mass_el = 0.000510998928;
const double mass_mu = 0.10566;
void triggerAnalysis(
int typeAna = 0,
bool applyFullLepSel = false,
TString typeLepSel = "medium"
){

  //*******************************************************
  //Input Files
  //*******************************************************
  TString filesPathOld = "/data/t3home000/ceballos/panda/v_005_0/";
  TString filesPath    = "/data/t3home000/ceballos/panda/v_006_0/";
  vector<TString> infileName_;

  if     (typeAna == 0){
    infileName_.push_back(Form("%sMET.root",filesPath.Data()));
  }
  else if(typeAna == 1){
    infileName_.push_back(Form("%sqqWW.root" ,filesPathOld .Data())); 
    infileName_.push_back(Form("%sggWW.root" ,filesPathOld .Data())); 
    infileName_.push_back(Form("%sDYJetsToLL_M-50_NLO.root",filesPathOld .Data()));
    infileName_.push_back(Form("%sqqZZ.root" ,filesPathOld .Data()));
    infileName_.push_back(Form("%sggZZ.root" ,filesPathOld .Data()));
    infileName_.push_back(Form("%sWZ.root" ,filesPathOld .Data()));  
    infileName_.push_back(Form("%sTT2L.root" ,filesPathOld .Data())); 
  }
  else {
    return;
  }

  // Initializations
  double nPassSel[3][2];
  for(int i=0; i<3; i++){
    for(int j=0; j<2; j++){
      nPassSel[i][j] = 0.0;
    }
  }

  TString triggerSuffix = "";

  double xminPlot = 0.0;
  double xmaxPlot = 1.0;
  int nBinPlot = 200;
  const int histBins = 3;
  const int allPlots = 6;
  TH1D* histo[histBins][allPlots];

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  5) {nBinPlot = 60; xminPlot =-0.5; xmaxPlot = 59.5;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[i][thePlot] = (TH1D*) histos->Clone(Form("histo_%d",i));
    histos->Clear();
  }

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infileName_.size(); ifile++) {
    printf("sampleNames(%d): %s\n",ifile,infileName_[ifile].Data());
    TFile *the_input_file = TFile::Open(infileName_[ifile].Data());
    TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("events");
    
    pandaFlat thePandaFlat(the_input_tree);
    Long64_t nentries = thePandaFlat.fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = thePandaFlat.LoadTree(jentry);
      if (ientry < 0) break;
      nb = thePandaFlat.fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%1000000 == 0) printf("--- reading event %8lld (%8lld) of %8lld\n",jentry,ientry,nentries);

      if(thePandaFlat.nLooseLep != 2) continue;

      bool passFakeLepId = ((thePandaFlat.looseLep1SelBit & kFake) == kFake) && ((thePandaFlat.looseLep2SelBit & kFake) == kFake);
      if(passFakeLepId == false) continue;

      bool passMediumLepId = ((thePandaFlat.looseLep1SelBit & kMedium) == kMedium) && ((thePandaFlat.looseLep2SelBit & kMedium) == kMedium);
      if(applyFullLepSel == true && passMediumLepId == false) continue;

      int typePair = -1;
      if     (TMath::Abs(thePandaFlat.looseLep1PdgId)==13&&TMath::Abs(thePandaFlat.looseLep2PdgId)==13) {typePair = 0;}
      else if(TMath::Abs(thePandaFlat.looseLep1PdgId)==11&&TMath::Abs(thePandaFlat.looseLep2PdgId)==11) {typePair = 1;}
      else if(TMath::Abs(thePandaFlat.looseLep1PdgId)!=TMath::Abs(thePandaFlat.looseLep2PdgId))         {typePair = 2;}
      else {assert(1); return;}

      double thePDGMass[2] = {mass_mu, mass_mu};
      if     (abs(typePair) == 1) {thePDGMass[0] = mass_el; thePDGMass[1] = mass_el;}
      else if(abs(typePair) == 2 && abs(thePandaFlat.looseLep1PdgId)==11) {thePDGMass[0] = mass_el;}
      else if(abs(typePair) == 2 && abs(thePandaFlat.looseLep2PdgId)==11) {thePDGMass[1] = mass_el;}
      TLorentzVector v1,v2,metP4;
      v1.SetPtEtaPhiM(thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      v2.SetPtEtaPhiM(thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      TLorentzVector dilep = v1+v2;

      if(v1.Pt() <= 20 || v2.Pt() <= 20 || dilep.M() <= 12) continue;

      int iPt[2] = {-1, -1};
      if     (v1.Pt() < 25) iPt[0] = 0;
      else if(v1.Pt() < 30) iPt[0] = 1;
      else if(v1.Pt() < 35) iPt[0] = 2;
      else if(v1.Pt() < 50) iPt[0] = 3;
      else		    iPt[0] = 4;
      if     (v2.Pt() < 25) iPt[1] = 0;
      else if(v2.Pt() < 30) iPt[1] = 1;
      else if(v2.Pt() < 35) iPt[1] = 2;
      else if(v2.Pt() < 50) iPt[1] = 3;
      else		    iPt[1] = 4;

      int iEta[2] = {-1, -1};
      if    (TMath::Abs(v1.Eta()) < 1.5) iEta[0] = 0;
      else                               iEta[0] = 1;
      if    (TMath::Abs(v2.Eta()) < 1.5) iEta[1] = 0;
      else                               iEta[1] = 1;

      double theBin = -1;
      if     (iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin = 0;
      else if(iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin = 1;
      else if(iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin = 2;
      else if(iPt[0] == 0 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin = 3;

      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin = 4;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin = 5;
      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin = 6;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin = 7;
      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin = 8;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin = 9;
      else if(iPt[0] == 1 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin =10;
      else if(iPt[0] == 1 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =11;

      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin =12;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin =13;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin =14;
      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin =15;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin =16;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin =17;
      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin =18;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin =19;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin =20;
      else if(iPt[0] == 2 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin =21;
      else if(iPt[0] == 2 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =22;
      else if(iPt[0] == 2 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =23;

      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin =24;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin =25;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin =26;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 0) theBin =27;
      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin =28;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin =29;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin =30;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 1) theBin =31;
      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin =32;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin =33;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin =34;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 0) theBin =35;
      else if(iPt[0] == 3 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin =36;
      else if(iPt[0] == 3 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =37;
      else if(iPt[0] == 3 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =38;
      else if(iPt[0] == 3 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 1) theBin =39;

      else if(iPt[0] == 4 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 0) theBin =40;
      else if(iPt[0] == 4 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 0) theBin =41;
      else if(iPt[0] == 4 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 0) theBin =42;
      else if(iPt[0] == 4 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 0) theBin =43;
      else if(iPt[0] == 4 && iPt[1] == 4 && iEta[0] == 0 && iEta[1] == 0) theBin =44;
      else if(iPt[0] == 4 && iPt[1] == 0 && iEta[0] == 0 && iEta[1] == 1) theBin =45;
      else if(iPt[0] == 4 && iPt[1] == 1 && iEta[0] == 0 && iEta[1] == 1) theBin =46;
      else if(iPt[0] == 4 && iPt[1] == 2 && iEta[0] == 0 && iEta[1] == 1) theBin =47;
      else if(iPt[0] == 4 && iPt[1] == 3 && iEta[0] == 0 && iEta[1] == 1) theBin =48;
      else if(iPt[0] == 4 && iPt[1] == 4 && iEta[0] == 0 && iEta[1] == 1) theBin =49;
      else if(iPt[0] == 4 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 0) theBin =50;
      else if(iPt[0] == 4 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 0) theBin =51;
      else if(iPt[0] == 4 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 0) theBin =52;
      else if(iPt[0] == 4 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 0) theBin =53;
      else if(iPt[0] == 4 && iPt[1] == 4 && iEta[0] == 1 && iEta[1] == 0) theBin =54;
      else if(iPt[0] == 4 && iPt[1] == 0 && iEta[0] == 1 && iEta[1] == 1) theBin =55;
      else if(iPt[0] == 4 && iPt[1] == 1 && iEta[0] == 1 && iEta[1] == 1) theBin =56;
      else if(iPt[0] == 4 && iPt[1] == 2 && iEta[0] == 1 && iEta[1] == 1) theBin =57;
      else if(iPt[0] == 4 && iPt[1] == 3 && iEta[0] == 1 && iEta[1] == 1) theBin =58;
      else if(iPt[0] == 4 && iPt[1] == 4 && iEta[0] == 1 && iEta[1] == 1) theBin =59;
 
      else {printf("IMPOSSIBLE\n");}

      double totalWeight = 1.0;
      if(typeAna != 0){
        totalWeight = thePandaFlat.sf_pu;
      }

      Bool_t passMETFilter = kTRUE;
      bool passFilter = (thePandaFlat.trigger & kMuEGTrig) == kMuEGTrig || (thePandaFlat.trigger & kMuMuTrig) == kMuMuTrig ||
                        (thePandaFlat.trigger & kMuTrig)   == kMuTrig   || (thePandaFlat.trigger & kEGEGTrig) == kEGEGTrig ||
		        (thePandaFlat.trigger & kEGTrig)   == kEGTrig;


      if(v2.Pt() > 25.0) nPassSel[typePair][0] = nPassSel[typePair][0] + totalWeight;

      histo[typePair][0]->Fill(theBin,totalWeight);

      if(passMETFilter == true) histo[typePair][2]->Fill(theBin,totalWeight);
      else                      histo[typePair][4]->Fill(theBin,totalWeight);

      if(passFilter){
        if(v2.Pt() > 25.0) nPassSel[typePair][1] = nPassSel[typePair][1] + totalWeight;
        histo[typePair][1]->Fill(theBin,totalWeight);
        if(passMETFilter == true) histo[typePair][3]->Fill(theBin,totalWeight);
        else                      histo[typePair][5]->Fill(theBin,totalWeight);
      }
    } // end events loop
    the_input_file->Close();
  } // end chain loop

  for(int i=0; i<3; i++){
    if(nPassSel[i][0] > 0){
      printf("trigger_eff(%d): %f +/- %f\n",i,nPassSel[i][1]/nPassSel[i][0],sqrt(nPassSel[i][1]/nPassSel[i][0]*(1-nPassSel[i][1]/nPassSel[i][0])/nPassSel[i][0]));
    }
  } // end loop over typePair (4)

  for(int theType=0; theType<histBins; theType++){
    histo[theType][1]->Divide(histo[theType][0]);
    histo[theType][3]->Divide(histo[theType][2]);
    histo[theType][5]->Divide(histo[theType][4]);
  }
  char output[200];
  sprintf(output,"histo_trigger_study_eff_%d_%s_%d.root",(int)applyFullLepSel,typeLepSel.Data(),typeAna);	
  TFile* outFilePlotsEff = new TFile(output,"recreate");
  outFilePlotsEff->cd();
  for(int theType=0; theType<histBins; theType++){
    histo[theType][1]->Write();
  }
  outFilePlotsEff->Close();

  sprintf(output,"histo_trigger_study_eff_%d_%s_%d_met.root",(int)applyFullLepSel,typeLepSel.Data(),typeAna);	
  TFile* outFilePlotsEffMET = new TFile(output,"recreate");
  outFilePlotsEffMET->cd();
  for(int theType=0; theType<histBins; theType++){
    histo[theType][3]->Write();
  }
  outFilePlotsEffMET->Close();

  sprintf(output,"histo_trigger_study_eff_%d_%s_%d_nomet.root",(int)applyFullLepSel,typeLepSel.Data(),typeAna);	
  TFile* outFilePlotsEffNOMET = new TFile(output,"recreate");
  outFilePlotsEffNOMET->cd();
  for(int theType=0; theType<histBins; theType++){
    histo[theType][5]->Write();
  }
  outFilePlotsEffNOMET->Close();
}
