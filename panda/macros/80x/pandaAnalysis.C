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

#include "MitAnalysisRunII/panda/macros/80x/pandaFlat.C"

enum TriggerBits {
    kMETTrig	   =(1<<0),
    kSinglePhoTrig =(1<<1),
    kMuEGTrig	   =(1<<2),
    kMuMuTrig	   =(1<<3),
    kMuTrig	   =(1<<4),
    kEGEGTrig	   =(1<<5),
    kEGTrig	   =(1<<6)
};

double PT(double pt1, double eta1, double phi1, double m1, double pt2, double eta2, double phi2, double m2) {
    TLorentzVector v1,v2;
    v1.SetPtEtaPhiM(pt1,eta1,phi1,m1);
    v2.SetPtEtaPhiM(pt2,eta2,phi2,m2);
    return (v1+v2).Pt();
}

const double mass_el = 0.000510998928;
const double mass_mu = 0.10566;
void pandaAnalysis()
{
  double lumi = 35.9;
  //*******************************************************
  //Input Files
  //*******************************************************
  TString filesPath = "/export/data/t3home000/ceballos/panda/v_004_0/";
  vector<TString> infileName_;  
  vector<Int_t> infileCat_;  
  infileName_.push_back(Form("%sdata.root",filesPath.Data()));                 infileCat_.push_back(0);
  infileName_.push_back(Form("%sDYJetsToLL_M-50_LO.root" ,filesPath.Data()));  infileCat_.push_back(1);
  infileName_.push_back(Form("%sDYJetsToLL_M-10to50.root" ,filesPath.Data())); infileCat_.push_back(1);
  infileName_.push_back(Form("%sqqWW.root" ,filesPath.Data()));                infileCat_.push_back(2);
  infileName_.push_back(Form("%sggWW.root" ,filesPath.Data()));                infileCat_.push_back(2);
  infileName_.push_back(Form("%sZZ.root" ,filesPath.Data()));                  infileCat_.push_back(3);
  infileName_.push_back(Form("%sWZ.root" ,filesPath.Data()));                  infileCat_.push_back(4);
  infileName_.push_back(Form("%sWGstar.root" ,filesPath.Data()));              infileCat_.push_back(4);
  infileName_.push_back(Form("%sVVV.root" ,filesPath.Data()));                 infileCat_.push_back(5);
  infileName_.push_back(Form("%sTTV.root" ,filesPath.Data()));                 infileCat_.push_back(5);
  infileName_.push_back(Form("%sTT2l.root" ,filesPath.Data()));                infileCat_.push_back(6);
  infileName_.push_back(Form("%sTW.root" ,filesPath.Data()));                  infileCat_.push_back(6);
  infileName_.push_back(Form("%sVG.root" ,filesPath.Data()));                  infileCat_.push_back(7);
  infileName_.push_back(Form("%sH125.root" ,filesPath.Data()));                infileCat_.push_back(8);

  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 10;
  const int histBins = 9;
  TH1D* histo[allPlots][histBins];
  TString processName[histBins] = {".Data", "DY", "WW", "ZZ", "WZ", "VVV", "Top", "VG", "H125"};
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  2) {nBinPlot =  60; xminPlot = 75.0; xmaxPlot = 105;}
    else if(thePlot >=  3 && thePlot <=  5) {nBinPlot = 200; xminPlot =  0.0; xmaxPlot = 200;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Reset();histos->Clear();
  }
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
      
      bool passTrigger = (thePandaFlat.trigger & kMuEGTrig) == kMuEGTrig || (thePandaFlat.trigger & kMuMuTrig) == kMuMuTrig ||
                         (thePandaFlat.trigger & kMuTrig)   == kMuTrig   || (thePandaFlat.trigger & kEGEGTrig) == kEGEGTrig ||
			 (thePandaFlat.trigger & kEGTrig)   == kEGTrig;
      if(passTrigger == false) continue;
      if(thePandaFlat.metFilter == 0) continue;
      
      if(thePandaFlat.nLooseLep != 2) continue;
      if(thePandaFlat.looseLep1Pt <= 25 || thePandaFlat.looseLep2Pt <= 20) continue;
      bool passLepId = ((thePandaFlat.looseLep1SelBit & 4) == 4) && ((thePandaFlat.looseLep2SelBit & 4) == 4);
      if(passLepId == false) continue;

      int lepType = -1;
      if     (abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 0;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 1;
      else if(abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 2;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId<0) lepType = 2;
      else if(abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 3;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 4;
      else if(abs(thePandaFlat.looseLep1PdgId)==13 && abs(thePandaFlat.looseLep2PdgId)==11 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 5;
      else if(abs(thePandaFlat.looseLep1PdgId)==11 && abs(thePandaFlat.looseLep2PdgId)==13 && thePandaFlat.looseLep1PdgId*thePandaFlat.looseLep2PdgId>0) lepType = 5;
      else printf("Impossible dilepton combination: %d %d\n",thePandaFlat.looseLep1PdgId,thePandaFlat.looseLep2PdgId);

      if(lepType > 2) continue;
      double thePDGMass[2] = {mass_mu, mass_mu};
      if     (lepType == 1) {thePDGMass[0] = mass_el; thePDGMass[1] = mass_el;}
      else if(lepType == 2 && abs(thePandaFlat.looseLep1PdgId)==11) {thePDGMass[0] = mass_el;}
      else if(lepType == 2 && abs(thePandaFlat.looseLep2PdgId)==11) {thePDGMass[1] = mass_el;}
      TLorentzVector v1,v2;
      v1.SetPtEtaPhiM(thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      v2.SetPtEtaPhiM(thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);

      if(TMath::Abs((v1+v2).M()-91.1876) >= 15) continue;

      double ZPt = (v1+v2).Pt();

      double totalWeight = 1.0;
      int theCategory = infileCat_[ifile];
      if(theCategory != 0){
        totalWeight = thePandaFlat.normalizedWeight * lumi * thePandaFlat.sf_pu *
	              thePandaFlat.sf_trk1 * thePandaFlat.sf_medium1 *
		      thePandaFlat.sf_trk2 * thePandaFlat.sf_medium2;
      }
      histo[lepType+0][theCategory]->Fill((v1+v2).M(),totalWeight);
      histo[lepType+3][theCategory]->Fill(TMath::Min(ZPt,199.999),totalWeight);
      
    }
    cout << nbytes << " " << nb << endl;
  }

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    char output[200];
    sprintf(output,"histozll_%d.root",thePlot);	
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    double totBck = 0;
    for(int i=1; i<=8; i++) totBck =totBck + histo[thePlot][i]->GetSumOfWeights();
    printf("(%d) %f (%f+%f+%f+%f+%f+%f+%f+%f)=%f\n",thePlot,histo[thePlot][0]->GetSumOfWeights(),
    histo[thePlot][1]->GetSumOfWeights(),histo[thePlot][2]->GetSumOfWeights(),
    histo[thePlot][3]->GetSumOfWeights(),histo[thePlot][4]->GetSumOfWeights(),
    histo[thePlot][5]->GetSumOfWeights(),histo[thePlot][6]->GetSumOfWeights(),
    histo[thePlot][6]->GetSumOfWeights(),histo[thePlot][8]->GetSumOfWeights(),
    totBck);
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }
}
