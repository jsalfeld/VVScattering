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
#include "TRandom.h"
#include "TLorentzVector.h"

#include "MitAnalysisRunII/panda/macros/80x/pandaFlat.C"
#include "MitAnalysisRunII/data/80x/RoccoR.cc"

double nPUScaleFactor(TH1D *fhDPU, float npu){
  double mynpu = TMath::Min(npu,(float)79.999);
  Int_t npuxbin = fhDPU->GetXaxis()->FindBin(mynpu);
  //return TMath::Min(fhDPU->GetBinContent(npuxbin),2.5);
  return fhDPU->GetBinContent(npuxbin);
}

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
void pandaZAnalysis(int whichDY = 0, bool isMIT = true, bool isTopSel = false)
{
  TString dirPathRM = TString(gSystem->Getenv("CMSSW_BASE")) + "/src/MitAnalysisRunII/data/80x/rcdata.2016.v3";
  RoccoR rmcor(dirPathRM.Data());
  double lumi = 35.8;
  double k_eff = 0.5 * sqrt(20285930./12446486.);
  //*******************************************************
  //Input Files
  //*******************************************************
  TString filesPathOld = "/data/t3home000/ceballos/panda/v_004_0/";
  TString filesPath    = "/data/t3home000/ceballos/panda/v_005_0/";
  if(isMIT == false) filesPath = "/afs/cern.ch/work/c/ceballos/public/samples/panda/v_005_0/";
  vector<TString> infileName_;
  vector<Int_t> infileCat_;
  infileName_.push_back(Form("%sdata.root",filesPath.Data()));                 infileCat_.push_back(0);
  infileName_.push_back(Form("%sdata.root",filesPathOld.Data()));              infileCat_.push_back(1);

/*
  //infileName_.push_back(Form("%sqqWW.root" ,filesPath.Data()));                infileCat_.push_back(1);
  //infileName_.push_back(Form("%sggWW.root" ,filesPath.Data()));                infileCat_.push_back(1);
  if     (whichDY == 0)
 {infileName_.push_back(Form("%sDYJetsToLL_M-50_LO.root" ,filesPath.Data()));  infileCat_.push_back(2);}
  else if(whichDY == 1)
 {infileName_.push_back(Form("%sDYJetsToLL_M-50_NLO.root",filesPath.Data()));  infileCat_.push_back(2);}
  else if(whichDY == 2)
 {infileName_.push_back(Form("%sDYJetsToLL_POWHEG.root",filesPath.Data()));    infileCat_.push_back(2);}
  else if(whichDY == 3)
 {infileName_.push_back(Form("%sDYJetsToLL_Pt0To50.root",filesPath.Data()));   infileCat_.push_back(2);
  infileName_.push_back(Form("%sDYJetsToLL_Pt50To100.root",filesPath.Data())); infileCat_.push_back(2);
  infileName_.push_back(Form("%sDYJetsToLL_Pt100To250.root",filesPath.Data()));infileCat_.push_back(2);
  infileName_.push_back(Form("%sDYJetsToLL_Pt250To400.root",filesPath.Data()));infileCat_.push_back(2);
  infileName_.push_back(Form("%sDYJetsToLL_Pt400To650.root",filesPath.Data()));infileCat_.push_back(2);
  infileName_.push_back(Form("%sDYJetsToLL_Pt650ToInf.root",filesPath.Data()));infileCat_.push_back(2);
  }
  infileName_.push_back(Form("%sDYJetsToLL_M-10to50.root" ,filesPath.Data())); infileCat_.push_back(2);
  //infileName_.push_back(Form("%sTT2L.root" ,filesPath.Data()));                infileCat_.push_back(3);
  //infileName_.push_back(Form("%sTW.root" ,filesPath.Data()));                  infileCat_.push_back(3);

  infileName_.push_back(Form("%sqqZZ.root" ,filesPath.Data()));                infileCat_.push_back(4);
  infileName_.push_back(Form("%sggZZ.root" ,filesPath.Data()));                infileCat_.push_back(4);
  infileName_.push_back(Form("%sWZ.root" ,filesPath.Data()));                  infileCat_.push_back(4);
  infileName_.push_back(Form("%sVVV.root" ,filesPath.Data()));                 infileCat_.push_back(4);
  infileName_.push_back(Form("%sTTV.root" ,filesPath.Data()));                 infileCat_.push_back(4);
  //infileName_.push_back(Form("%sWGstar.root" ,filesPath.Data()));              infileCat_.push_back(4);
  //infileName_.push_back(Form("%sVG.root" ,filesPath.Data()));                  infileCat_.push_back(6);
  //infileName_.push_back(Form("%sH125.root" ,filesPath.Data()));                infileCat_.push_back(7);
*/

  TFile *fPUFile = TFile::Open(Form("MitAnalysisRunII/data/90x/puWeights_90x.root"));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights")); assert(fhDPU); fhDPU->SetDirectory(0);
  delete fPUFile;

  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 60;
  const int histBins = 9;
  TH1D* histo[allPlots][histBins];
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  1) {nBinPlot = 120; xminPlot = 91.1876-15; xmaxPlot = 91.1876+15;}
    else if(thePlot >=  2 && thePlot <=  2) {nBinPlot = 200; xminPlot = 20.0; xmaxPlot = 220;}
    else if(thePlot >=  3 && thePlot <=  5) {nBinPlot =  10; xminPlot = -0.5; xmaxPlot = 9.5;}
    else if(thePlot >=  6 && thePlot <=  8) {nBinPlot =   5; xminPlot = -0.5; xmaxPlot = 4.5;}
    else if(thePlot >=  9 && thePlot <= 11) {nBinPlot =1000; xminPlot =  0.0; xmaxPlot =1000;}
    else if(thePlot >= 12 && thePlot <= 26) {nBinPlot = 500; xminPlot =  0.0; xmaxPlot = 500;}
    else if(thePlot >= 27 && thePlot <= 32) {nBinPlot = 100; xminPlot =  0.0; xmaxPlot = TMath::Pi();}
    else if(thePlot >= 33 && thePlot <= 38) {nBinPlot = 100; xminPlot =  0.0; xmaxPlot = 2.0;}
    else if(thePlot >= 39 && thePlot <= 44) {nBinPlot = 100; xminPlot =  0.0; xmaxPlot = 5.0;}
    else if(thePlot >= 45 && thePlot <= 47) {nBinPlot =  80; xminPlot = -0.5; xmaxPlot = 79.5;}
    else if(thePlot >= 48 && thePlot <= 53) {nBinPlot = 100; xminPlot = -2.5; xmaxPlot = 2.5;}
    else if(thePlot >= 54 && thePlot <= 59) {nBinPlot = 200; xminPlot = 25.0; xmaxPlot = 225;}
    
    if(isTopSel == true && (thePlot >=  0 && thePlot <=  1)) {nBinPlot = 200; xminPlot = 20.0; xmaxPlot = 220;}

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
      //if(thePandaFlat.metFilter == 0) continue;
      
      if(thePandaFlat.nLooseLep != 2) continue;

      bool passLepId = ((thePandaFlat.looseLep1SelBit & kMedium) == kMedium) && ((thePandaFlat.looseLep2SelBit & kMedium) == kMedium);
      if(passLepId == false) continue;
      int theCategory = infileCat_[ifile];

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

      if(lepType >= 3 || lepType == -1) continue;

      double thePDGMass[2] = {mass_mu, mass_mu};
      if     (abs(lepType) == 1) {thePDGMass[0] = mass_el; thePDGMass[1] = mass_el;}
      else if(abs(lepType) == 2 && abs(thePandaFlat.looseLep1PdgId)==11) {thePDGMass[0] = mass_el;}
      else if(abs(lepType) == 2 && abs(thePandaFlat.looseLep2PdgId)==11) {thePDGMass[1] = mass_el;}
      TLorentzVector v1,v2,metP4;
      v1.SetPtEtaPhiM(thePandaFlat.looseLep1Pt,thePandaFlat.looseLep1Eta,thePandaFlat.looseLep1Phi,thePDGMass[0]);
      v2.SetPtEtaPhiM(thePandaFlat.looseLep2Pt,thePandaFlat.looseLep2Eta,thePandaFlat.looseLep2Phi,thePDGMass[1]);
      TLorentzVector dilep = v1+v2;
      metP4.SetPtEtaPhiM(thePandaFlat.pfmet,0,thePandaFlat.pfmetphi,0);

      double dPhiDiLepMET = TMath::Abs(dilep.DeltaPhi(metP4));
      double ptFrac = TMath::Abs(dilep.Pt()-metP4.Pt())/dilep.Pt();
      double caloMinusPFMETRel = TMath::Abs(thePandaFlat.calomet-thePandaFlat.pfmet)/thePandaFlat.pfmet;
      double dphill = TMath::Abs(v1.DeltaPhi(v2));
      double detall = TMath::Abs(v1.Eta()-v2.Eta());
      double drll = sqrt(dphill*dphill+detall*detall);
      double mtW = TMath::Sqrt(2.0*dilep.Pt()*metP4.Pt()*(1.0 - TMath::Cos(dPhiDiLepMET)));

      bool passSel = (TMath::Abs((v1+v2).M()-91.1876) < 15 || (lepType == 2 && (v1+v2).M() > 20)) && v1.Pt() > 25 && v2.Pt() > 25;
      if(isTopSel == true) passSel = (TMath::Abs((v1+v2).M()-91.1876) > 15 || lepType == 2) && (v1+v2).M() > 20 && v1.Pt() > 25 && v2.Pt() > 25 && mtW > 50;
      if(passSel == false) continue;

      double totalWeight = 1.0;
      if(theCategory != 0 && theCategory != 1){
        totalWeight = thePandaFlat.normalizedWeight * lumi * thePandaFlat.sf_pu *
	              thePandaFlat.sf_trk1 * thePandaFlat.sf_medium1 *
		      thePandaFlat.sf_trk2 * thePandaFlat.sf_medium2;
      }
      else if(theCategory == 1){
        totalWeight = 0.36335 * nPUScaleFactor(fhDPU,thePandaFlat.npv);
        //totalWeight = 1.0;
      }

      histo[lepType+ 0][theCategory]->Fill((v1+v2).M(),totalWeight);
      histo[lepType+ 3][theCategory]->Fill(TMath::Min((double)thePandaFlat.nJet, 9.4999),totalWeight);
      histo[lepType+ 6][theCategory]->Fill(TMath::Min((double)thePandaFlat.jetNLBtags,4.4999),totalWeight);
      histo[lepType+ 9][theCategory]->Fill(TMath::Min((v1+v2).Pt(), 999.999),totalWeight);
      histo[lepType+12][theCategory]->Fill(TMath::Min((double)thePandaFlat.pfmet, 499.999),totalWeight);
      histo[lepType+15][theCategory]->Fill(TMath::Min((double)thePandaFlat.puppimet, 499.999),totalWeight);
      histo[lepType+18][theCategory]->Fill(TMath::Min((double)thePandaFlat.calomet, 499.999),totalWeight);
      histo[lepType+21][theCategory]->Fill(TMath::Min((double)thePandaFlat.trkmet, 499.999),totalWeight);
      histo[lepType+24][theCategory]->Fill(TMath::Min(mtW, 499.999),totalWeight);
      histo[lepType+27][theCategory]->Fill(dPhiDiLepMET,totalWeight);
      histo[lepType+30][theCategory]->Fill(dphill,totalWeight);
      histo[lepType+33][theCategory]->Fill(TMath::Min(ptFrac,1.999),totalWeight);
      histo[lepType+36][theCategory]->Fill(TMath::Min(caloMinusPFMETRel,1.999),totalWeight);
      histo[lepType+39][theCategory]->Fill(TMath::Min(detall,4.999),totalWeight);
      histo[lepType+42][theCategory]->Fill(TMath::Min(drll,4.999),totalWeight);
      histo[lepType+45][theCategory]->Fill(TMath::Min((double)thePandaFlat.npv,79.499),totalWeight);
      histo[lepType+48][theCategory]->Fill(thePandaFlat.looseLep1Eta,totalWeight);
      histo[lepType+51][theCategory]->Fill(thePandaFlat.looseLep2Eta,totalWeight);
      histo[lepType+54][theCategory]->Fill(TMath::Min((double)thePandaFlat.looseLep1Pt, 224.999),totalWeight);
      histo[lepType+57][theCategory]->Fill(TMath::Min((double)thePandaFlat.looseLep2Pt, 224.999),totalWeight);

    } // end event loop
  } // end samples loop

  char output[200];
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    TString addSuffix = "";
    if(isTopSel == true) addSuffix = "_topsel";
    sprintf(output,"histoDY%dzll_%d%s.root",whichDY,thePlot,addSuffix.Data());	
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
