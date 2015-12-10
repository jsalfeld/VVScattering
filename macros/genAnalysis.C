#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include "NeroProducer/Core/interface/BareEvent.hpp"
#include "NeroProducer/Core/interface/BareJets.hpp"
#include "NeroProducer/Core/interface/BareLeptons.hpp"
#include "NeroProducer/Core/interface/BareTaus.hpp"
#include "NeroProducer/Core/interface/BareMet.hpp"
#include "NeroProducer/Core/interface/BareTrigger.hpp"
#include "NeroProducer/Core/interface/BareVertex.hpp"
#include "NeroProducer/Core/interface/BareMonteCarlo.hpp"

#include "MitAnalysisRunII/macros/factors.h"

void genAnalysis(
 Int_t period = 1
 ){

  TString filesPath  = "/scratch5/ceballos/ntuples_noweights/";
  Double_t lumi = 0.0715;
  if(period == 1) lumi = 2.2;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  if      (period==0){
  puPath = "/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src/MitAnalysisRunII/data/puWeights_13TeV_50ns.root";
  //infilenamev.push_back(Form("%sdata_AOD_50ns.root",filesPath.Data()));														  infilecatv.push_back(0);
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1+AODSIM.root",filesPath.Data()));	  infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));  	  infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data())); infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sWZ_TuneCUETP8M1_13TeV-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));					  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sZZ_TuneCUETP8M1_13TeV-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2+AODSIM.root",filesPath.Data()));					  infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1+AODSIM.root",filesPath.Data()));		  infilecatv.push_back(5);
  }
  else if(period==1){
  puPath = "/home/ceballos/cms/cmssw/042/CMSSW_7_4_6/src/MitAnalysisRunII/data/puWeights_13TeV_25ns.root";
  //infilenamev.push_back(Form("%sdata_AOD_25ns.root",filesPath.Data()));														  infilecatv.push_back(0);
  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(1);
  //infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	  //infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3+AODSIM.root",filesPath.Data()));  	          //infilecatv.push_back(2);
  //infilenamev.push_back(Form("%sTTTo2L2Nu_13TeV-powheg+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));						  infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sST_tW_top_5f_//inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));      //infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sST_tW_antitop_5f_//inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));  //infilecatv.push_back(3);
  //infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));					  //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madsp//in_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));					  //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWZTo1L1Nu2Q_13TeV_amcatnloFXFX_madsp//in_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madsp//in_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2+AODSIM.root",filesPath.Data()));			  //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));			  //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));				  //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madsp//in-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	  //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madsp//in-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));	  //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		  //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madsp//in-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		  //infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		  //infilecatv.push_back(5);
  //infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));                     //infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_amcatnloFXFX_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		  //infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		          //infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sHWm//inusJ_HToWW_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		                  //infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sHWplusJ_HToWW_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		                  //infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sHZJ_HToWW_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		                  //infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sGluGluZH_HToWW_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		                  //infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		                  //infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		                  //infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sWplusHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		                  //infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sWm//inusHToTauTau_M125_13TeV_powheg_pythia8+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));		                  //infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madsp//in_pythia8_mWCutfix+RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1+AODSIM.root",filesPath.Data()));      //infilecatv.push_back(7);
  }
  else {assert(0);}
  
  //infilenamev.clear();infilecatv.clear();
  //infilenamev.push_back(Form("/scratch5/ceballos/test/test.root"));   infilecatv.push_back(1);

  if(infilenamev.size() != infilecatv.size()) assert(0);

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);
 
  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 160;
  TH1D* histo[allPlots];

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <= 39) {nBinPlot =  20; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >= 40 && thePlot <= 79) {nBinPlot =   5; xminPlot = 0.0; xmaxPlot =   2.5;}
    else if(thePlot >= 80 && thePlot <= 80) {nBinPlot =   5; xminPlot =-0.5; xmaxPlot =   4.5;}
    else if(thePlot >= 81 && thePlot <= 90) {nBinPlot = 100; xminPlot = 0.0; xmaxPlot =   1.0;}
    else if(thePlot >=100 && thePlot <=109) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >=110 && thePlot <=119) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot =   5.0;}
    else if(thePlot >=120 && thePlot <=129) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >=130 && thePlot <=139) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot =   5.0;}
    else if(thePlot >=140 && thePlot <=149) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot = 100.0;}
    else if(thePlot >=150 && thePlot <=159) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot =   5.0;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    histo[thePlot] = (TH1D*) histos->Clone(Form("histo%d",thePlot));
    histos->Clear();
  }

  //*******************************************************
  // Chain Loop
  //*******************************************************
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {

    TFile the_input_file(infilenamev[ifile]);
    TTree *the_input_tree = (TTree*)the_input_file.FindObjectAny("events");
    //TTree *the_input_all  = (TTree*)the_input_file.FindObjectAny("all");

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

    double isoCut;
    for (int i=0; i<int(the_input_tree->GetEntries()/1.); ++i) {
      the_input_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      double totalWeight = 1;//eventMonteCarlo.mcWeight*lumi;

      vector<bool> isGenDupl; int nGoodGenLeptons = 0;
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
	if(isGenDupl[ngen0] == 0) nGoodGenLeptons++;
      }
      histo[80]->Fill((double)nGoodGenLeptons,totalWeight);
      
      for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	if(isGenDupl[ngen] == 1) continue;
        if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt() <= 10) continue;
	if(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()) >= 2.5) continue;

	int nCount = 0;
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==11) nCount = 10;
	histo[nCount+ 0]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt(),99.999),totalWeight);
	histo[nCount+40]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()),totalWeight);

        int whichRecoLepton = -1;
        for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
          if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[nlep])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.1) {
	    whichRecoLepton = nlep;
	    break;
	  }
	}
	if(whichRecoLepton >= 0) {
	  if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==13) histo[81]->Fill(TMath::Min((double)(*eventLeptons.iso)   [whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt(),0.999),totalWeight);
	  if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==11) histo[82]->Fill(TMath::Min((double)(*eventLeptons.iso)   [whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt(),0.999),totalWeight);
	  if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==13) histo[83]->Fill(TMath::Min((double)(*eventLeptons.chIso) [whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt(),0.999),totalWeight);
	  if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==11) histo[84]->Fill(TMath::Min((double)(*eventLeptons.chIso) [whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt(),0.999),totalWeight);
	  if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==13) histo[85]->Fill(TMath::Min((double)(*eventLeptons.nhIso) [whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt(),0.999),totalWeight);
	  if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==11) histo[86]->Fill(TMath::Min((double)(*eventLeptons.nhIso) [whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt(),0.999),totalWeight);
	  if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==13) histo[87]->Fill(TMath::Min((double)(*eventLeptons.phoIso)[whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt(),0.999),totalWeight);
	  if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==11) histo[88]->Fill(TMath::Min((double)(*eventLeptons.phoIso)[whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt(),0.999),totalWeight);
	  if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==13) histo[89]->Fill(TMath::Min((double)(*eventLeptons.puIso) [whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt(),0.999),totalWeight);
	  if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen])==11) histo[90]->Fill(TMath::Min((double)(*eventLeptons.puIso) [whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt(),0.999),totalWeight);

	  if(((int)(*eventLeptons.selBits)[whichRecoLepton] & BareLeptons::LepBaseline) == BareLeptons::LepBaseline) { 
	    histo[nCount+ 1]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt(),99.999),totalWeight);
   	    histo[nCount+41]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[whichRecoLepton]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Eta()) ? 0.1260 : 0.1440);
	  if(((int)(*eventLeptons.selBits)[whichRecoLepton] & BareLeptons::LepVeto) == BareLeptons::LepVeto 
	  && (double)(*eventLeptons.iso)[whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt() < isoCut) {
	    histo[nCount+ 2]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt(),99.999),totalWeight);
   	    histo[nCount+42]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()),totalWeight);
	  }

	  if(((int)(*eventLeptons.selBits)[whichRecoLepton] & BareLeptons::LepFake) == BareLeptons::LepFake) {
	    histo[nCount+ 3]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt(),99.999),totalWeight);
   	    histo[nCount+43]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[whichRecoLepton]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Eta()) ? 0.0893 : 0.1210);
	  if(((int)(*eventLeptons.selBits)[whichRecoLepton] & BareLeptons::LepLoose) == BareLeptons::LepLoose 
	  && (double)(*eventLeptons.iso)[whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt() < isoCut) {
	    histo[nCount+ 4]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt(),99.999),totalWeight);
   	    histo[nCount+44]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[whichRecoLepton]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Eta()) ? 0.0766 : 0.0678);
	  if(((int)(*eventLeptons.selBits)[whichRecoLepton] & BareLeptons::LepMedium) == BareLeptons::LepMedium 
	  && (double)(*eventLeptons.iso)[whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt() < isoCut) {
	    histo[nCount+ 5]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt(),99.999),totalWeight);
   	    histo[nCount+45]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[whichRecoLepton]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Eta()) ? 0.0354 : 0.0646);
	  if(((int)(*eventLeptons.selBits)[whichRecoLepton] & BareLeptons::LepTight) == BareLeptons::LepTight 
	  && (double)(*eventLeptons.iso)[whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt() < isoCut) {
	    histo[nCount+ 6]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt(),99.999),totalWeight);
   	    histo[nCount+46]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()),totalWeight);
	  }

	  if(((int)(*eventLeptons.selBits)[whichRecoLepton] & BareLeptons::LepSoftIP) == BareLeptons::LepSoftIP) {
	    histo[nCount+ 7]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt(),99.999),totalWeight);
   	    histo[nCount+47]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[whichRecoLepton]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Eta()) ? 0.0766 : 0.0678);
	  if(((int)(*eventLeptons.selBits)[whichRecoLepton] & BareLeptons::LepMediumIP) == BareLeptons::LepMediumIP 
	  && (double)(*eventLeptons.iso)[whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt() < isoCut) {
	    histo[nCount+ 8]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt(),99.999),totalWeight);
   	    histo[nCount+48]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[whichRecoLepton]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Eta()) ? 0.0354 : 0.0646);
	  if(((int)(*eventLeptons.selBits)[whichRecoLepton] & BareLeptons::LepTightIP) == BareLeptons::LepTightIP 
	  && (double)(*eventLeptons.iso)[whichRecoLepton]/((TLorentzVector*)(*eventLeptons.p4)[whichRecoLepton])->Pt() < isoCut) {
	    histo[nCount+ 9]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt(),99.999),totalWeight);
   	    histo[nCount+49]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()),totalWeight);
	  }
	}
      } // end loop over gen leptons


      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
	int nCount = 20;
	if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep])==11) nCount = 30;
        int whichGenLepton = -1;
        for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
	  if(isGenDupl[ngen] == 1) continue;
          if(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Pt() <= 10) continue;
	  if(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])->Eta()) >= 2.5) continue;
          if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[nlep])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.1) {
	    whichGenLepton = nlep;
	    break;
	  }
	}

	if(whichGenLepton == -1) {
	  histo[nCount+ 0]->Fill(TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt(),99.999),totalWeight);
   	  histo[nCount+40]->Fill(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),totalWeight);

	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepBaseline) == BareLeptons::LepBaseline) { 
	    histo[nCount+ 1]->Fill(TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt(),99.999),totalWeight);
   	    histo[nCount+41]->Fill(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()) ? 0.1260 : 0.1440);
	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepVeto) == BareLeptons::LepVeto 
	  && (double)(*eventLeptons.iso)[nlep]/((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt() < isoCut) {
	    histo[nCount+ 2]->Fill(TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt(),99.999),totalWeight);
   	    histo[nCount+42]->Fill(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),totalWeight);
	  }

	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake) == BareLeptons::LepFake) {
	    histo[nCount+ 3]->Fill(TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt(),99.999),totalWeight);
   	    histo[nCount+43]->Fill(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()) ? 0.0893 : 0.1210);
	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepLoose) == BareLeptons::LepLoose 
	  && (double)(*eventLeptons.iso)[nlep]/((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt() < isoCut) {
	    histo[nCount+ 4]->Fill(TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt(),99.999),totalWeight);
   	    histo[nCount+44]->Fill(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()) ? 0.0766 : 0.0678);
	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepMedium) == BareLeptons::LepMedium 
	  && (double)(*eventLeptons.iso)[nlep]/((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt() < isoCut) {
	    histo[nCount+ 5]->Fill(TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt(),99.999),totalWeight);
   	    histo[nCount+45]->Fill(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()) ? 0.0354 : 0.0646);
	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepTight) == BareLeptons::LepTight 
	  && (double)(*eventLeptons.iso)[nlep]/((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt() < isoCut) {
	    histo[nCount+ 6]->Fill(TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt(),99.999),totalWeight);
   	    histo[nCount+46]->Fill(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),totalWeight);
	  }

	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepSoftIP) == BareLeptons::LepSoftIP) {
	    histo[nCount+ 7]->Fill(TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt(),99.999),totalWeight);
   	    histo[nCount+47]->Fill(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()) ? 0.0766 : 0.0678);
	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepMediumIP) == BareLeptons::LepMediumIP 
	  && (double)(*eventLeptons.iso)[nlep]/((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt() < isoCut) {
	    histo[nCount+ 8]->Fill(TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt(),99.999),totalWeight);
   	    histo[nCount+48]->Fill(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),totalWeight);
	  }

          isoCut = 0.12;
	  if(TMath::Abs((int)(*eventLeptons.pdgId)[nlep]) == 11) isoCut = (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()) ? 0.0354 : 0.0646);
	  if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepTightIP) == BareLeptons::LepTightIP 
	  && (double)(*eventLeptons.iso)[nlep]/((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt() < isoCut) {
	    histo[nCount+ 9]->Fill(TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt(),99.999),totalWeight);
   	    histo[nCount+49]->Fill(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),totalWeight);
	  }
	}

      } // end loop over reco leptons

      vector<int> idLep;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut("medium",TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep]))
	                                                                                               {idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idLep.push_back(nlep);}
      }

      vector<int> idGenJet;
      for(int ngenj=0; ngenj<eventMonteCarlo.jetP4->GetEntriesFast(); ngenj++) {
        Bool_t isGenLepton = kFALSE;
        for(int ngenl=0; ngenl<eventMonteCarlo.p4->GetEntriesFast(); ngenl++) {
	if(TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngenl]) != 12 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngenl]) != 14 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngenl]) != 16 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngenl]) != 11 &&
	   TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngenl]) != 13) continue;
          if(((TLorentzVector*)(*eventMonteCarlo.jetP4)[ngenj])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngenl])) < 0.3) isGenLepton = kTRUE;
        }
        if(isGenLepton == kTRUE) continue;
        if(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[ngenj])->Eta()) >= 5) continue;
        idGenJet.push_back(ngenj);
      }

      vector<int> idJet;
      vector<int> idMVAJet;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 30) continue;        
	idJet.push_back(nj);

        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        if(passId == false) continue;        
	idMVAJet.push_back(nj);
      }

      for(unsigned int ngenj=0; ngenj<idGenJet.size(); ngenj++) {
	int nCount = 100;
	if(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Eta()) < 2.5)
	histo[nCount+ 0]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Pt(),99.999),totalWeight);
	if(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Pt() > 30)
	histo[nCount+10]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Eta()),totalWeight);
	
	bool isRecoJet = kFALSE;
        for(unsigned int nj=0; nj<idJet.size(); nj++){
          if(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->DeltaR(*((TLorentzVector*)(*eventJets.p4)[idJet[nj]])) < 0.4) isRecoJet = kTRUE;
        }
        if(isRecoJet == kTRUE) {
	  if(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Eta()) < 2.5)
	  histo[nCount+ 1]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Pt(),99.999),totalWeight);
	  if(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Pt() > 30)
	  histo[nCount+11]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Eta()),totalWeight);
	}

        isRecoJet = kFALSE;
        for(unsigned int nj=0; nj<idMVAJet.size(); nj++){
          if(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->DeltaR(*((TLorentzVector*)(*eventJets.p4)[idMVAJet[nj]])) < 0.4) isRecoJet = kTRUE;
        }
        if(isRecoJet == kTRUE) {
	  if(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Eta()) < 2.5)
	  histo[nCount+ 2]->Fill(TMath::Min(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Pt(),99.999),totalWeight);
	  if(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Pt() > 30)
	  histo[nCount+12]->Fill(TMath::Abs(((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])->Eta()),totalWeight);
	}
      }

      for(unsigned int nj=0; nj<idJet.size(); nj++) {
	int nCount = 120;
	if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[nj]])->Eta()) < 2.5)
	histo[nCount+ 0]->Fill(TMath::Min(((TLorentzVector*)(*eventJets.p4)[idJet[nj]])->Pt(),99.999),totalWeight);
	if(((TLorentzVector*)(*eventJets.p4)[idJet[nj]])->Pt() > 30)
	histo[nCount+10]->Fill(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[nj]])->Eta()),totalWeight);
	
	bool isGenJet = kFALSE;
        for(unsigned int ngenj=0; ngenj<idGenJet.size(); ngenj++){
          if(((TLorentzVector*)(*eventJets.p4)[idJet[nj]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])) < 0.4) isGenJet = kTRUE;
        }
        if(isGenJet == kFALSE) {
  	  if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[nj]])->Eta()) < 2.5)
  	  histo[nCount+ 1]->Fill(TMath::Min(((TLorentzVector*)(*eventJets.p4)[idJet[nj]])->Pt(),99.999),totalWeight);
	  if(((TLorentzVector*)(*eventJets.p4)[idJet[nj]])->Pt() > 30)
	  histo[nCount+11]->Fill(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idJet[nj]])->Eta()),totalWeight);
	}
      }

      for(unsigned int nj=0; nj<idMVAJet.size(); nj++) {
	int nCount = 140;
        if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idMVAJet[nj]])->Eta()) < 2.5)
	histo[nCount+ 0]->Fill(TMath::Min(((TLorentzVector*)(*eventJets.p4)[idMVAJet[nj]])->Pt(),99.999),totalWeight);
	if(((TLorentzVector*)(*eventJets.p4)[idMVAJet[nj]])->Pt() > 30)
	histo[nCount+10]->Fill(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idMVAJet[nj]])->Eta()),totalWeight);
	
	bool isGenJet = kFALSE;
        for(unsigned int ngenj=0; ngenj<idGenJet.size(); ngenj++){
          if(((TLorentzVector*)(*eventJets.p4)[idMVAJet[nj]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.jetP4)[idGenJet[ngenj]])) < 0.4) isGenJet = kTRUE;
        }
        if(isGenJet == kFALSE) {
  	  if(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idMVAJet[nj]])->Eta()) < 2.5)
  	  histo[nCount+ 1]->Fill(TMath::Min(((TLorentzVector*)(*eventJets.p4)[idMVAJet[nj]])->Pt(),99.999),totalWeight);
	  if(((TLorentzVector*)(*eventJets.p4)[idMVAJet[nj]])->Pt() > 30)
	  histo[nCount+11]->Fill(TMath::Abs(((TLorentzVector*)(*eventJets.p4)[idMVAJet[nj]])->Eta()),totalWeight);
	}
      }
    } // end of loop
  } // end of chain

  for(int np= 1; np<10; np++) histo[np]->Divide(histo[ 0]);
  for(int np=11; np<20; np++) histo[np]->Divide(histo[10]);
  for(int np=21; np<30; np++) histo[np]->Divide(histo[20]);
  for(int np=31; np<40; np++) histo[np]->Divide(histo[30]);
  for(int np=41; np<50; np++) histo[np]->Divide(histo[40]);
  for(int np=51; np<60; np++) histo[np]->Divide(histo[50]);
  for(int np=61; np<70; np++) histo[np]->Divide(histo[60]);
  for(int np=71; np<80; np++) histo[np]->Divide(histo[70]);
  for(int np=80; np<=90; np++) histo[np]->Scale(1./histo[np]->GetSumOfWeights());
  for(int np=101; np<110; np++) histo[np]->Divide(histo[100]);
  for(int np=111; np<120; np++) histo[np]->Divide(histo[110]);
  for(int np=121; np<130; np++) histo[np]->Divide(histo[120]);
  for(int np=131; np<140; np++) histo[np]->Divide(histo[130]);
  for(int np=141; np<150; np++) histo[np]->Divide(histo[140]);
  for(int np=151; np<160; np++) histo[np]->Divide(histo[150]);
  char output[200];
  sprintf(output,"histo_geneff.root");	  
  TFile* outFilePlotsNote = new TFile(output,"recreate");
  outFilePlotsNote->cd();
  for(int np=0; np<160; np++) histo[np]->Write();
  outFilePlotsNote->Close();
}
