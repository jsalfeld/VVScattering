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

#include "MitAnalysisRunII/macros/76x/factors.h"

double mcPrescale = 1.0;

void QCDAnalysis(
 Int_t nsel = 0,
 Int_t typeSel = 4,
 Int_t applyPrescale = 1,
 TString typeLepSel = "medium"
 ){

  Int_t period = 1;
  TString filesPathDA  = "/scratch/ceballos/ntuples_weightsDA_76x/qcd_";
  TString filesPathMC  = "/scratch5/ceballos/ntuples_weightsMC_76x/qcd_";
  Double_t lumi = 2.318;

  Double_t prescale[5];

  double denFRDA[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double numFRDA[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double denFRBG[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double numFRBG[5][5] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  
  vector<Int_t> infilecatv;  

  TString puPath = "";
  if     (period==0){
  }
  else if(period==1){
  if     (typeSel == 11) {prescale[0]=0.00000;prescale[1]=0.00828;prescale[2]=0.00499;prescale[3]=0.00716;prescale[4]=0.00887;}
  else if(typeSel == 13) {prescale[0]=0.00245;prescale[1]=0.06880;prescale[2]=0.09568;prescale[3]=0.09479;prescale[4]=0.09383;}

  puPath = "MitAnalysisRunII/data/76x/puWeights_76x.root";
  infilenamev.push_back(Form("%sdata_AOD_Run2015C_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);
  infilenamev.push_back(Form("%sdata_AOD_Run2015D_25ns.root",filesPathDA.Data()));											      infilecatv.push_back(0);

  infilenamev.push_back(Form("%sWWTo2L2Nu_13TeV-powheg+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                                            infilecatv.push_back(1);
  infilenamev.push_back(Form("%sGluGluWWTo2L2Nu_MCFM_13TeV+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(1);

  infilenamev.push_back(Form("%sDYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(2);
  infilenamev.push_back(Form("%sDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1+AODSIM_0.root",filesPathMC.Data()));	   infilecatv.push_back(2);

  infilenamev.push_back(Form("%sTT_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1+AODSIM_0.root",filesPathMC.Data()));                       infilecatv.push_back(3); // tt->X (not tt->2l)
  infilenamev.push_back(Form("%sST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));    infilecatv.push_back(3);
  infilenamev.push_back(Form("%sST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));infilecatv.push_back(3);

  infilenamev.push_back(Form("%sZZTo2L2Nu_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  				   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sZZTo4L_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));					   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4e_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		 	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sGluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sWZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v2+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sWWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sWZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                           infilecatv.push_back(4);
  //infilenamev.push_back(Form("%sZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                         infilecatv.push_back(4);
  infilenamev.push_back(Form("%sTTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sTTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));	   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sTTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sTTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(4);
  infilenamev.push_back(Form("%sTTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));           infilecatv.push_back(4);

  infilenamev.push_back(Form("%sWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(5);

  infilenamev.push_back(Form("%sWGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(6);
  //infilenamev.push_back(Form("%sZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                 infilecatv.push_back(6);

  infilenamev.push_back(Form("%sGluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));		   infilecatv.push_back(7);
  infilenamev.push_back(Form("%sVBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                  infilecatv.push_back(7);
  infilenamev.push_back(Form("%sGluGluHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));			   infilecatv.push_back(7);
  infilenamev.push_back(Form("%sVBFHToTauTau_M125_13TeV_powheg_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));                            infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sVHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data())); 		   infilecatv.push_back(7);
  //infilenamev.push_back(Form("%sttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix+RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1+AODSIM.root",filesPathMC.Data()));  infilecatv.push_back(7);
  }

  //infilenamev.push_back(Form("nero.root"));     infilecatv.push_back(0);

  if(infilenamev.size() != infilecatv.size()) assert(0);

  Float_t fMVACut[4][4];
  InitializeJetIdCuts(fMVACut);

  TFile *fPUFile = TFile::Open(Form("%s",puPath.Data()));
  TH1D *fhDPU = (TH1D*)(fPUFile->Get("puWeights"));
  assert(fhDPU);
  fhDPU->SetDirectory(0);
  delete fPUFile;

  TFile *fElSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDElMediumSF = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Medium_ele"));
  TH2D *fhDElTightSF  = (TH2D*)(fElSF->Get("unfactorized_scalefactors_Tight_ele"));
  TH2D *fhDElMediumMVASF = (TH2D*)(fElSF->Get("unfactorized_scalefactors_MediumMVA_ele"));
  TH2D *fhDElTightMVASF  = (TH2D*)(fElSF->Get("unfactorized_scalefactors_TightMVA_ele"));
  assert(fhDElMediumSF);
  assert(fhDElTightSF);
  assert(fhDElMediumMVASF);
  assert(fhDElTightMVASF);
  fhDElMediumSF->SetDirectory(0);
  fhDElTightSF ->SetDirectory(0);
  fhDElMediumMVASF->SetDirectory(0);
  fhDElTightMVASF ->SetDirectory(0);
  delete fElSF;

  TFile *fMuSF = TFile::Open(Form("MitAnalysisRunII/data/76x/scalefactors_hww.root"));
  TH2D *fhDMuMediumSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Medium_mu"));
  TH2D *fhDMuIsoSF = (TH2D*)(fMuSF->Get("unfactorized_scalefactors_Iso_mu"));
  assert(fhDMuMediumSF);
  assert(fhDMuIsoSF);
  fhDMuMediumSF->SetDirectory(0);
  fhDMuIsoSF->SetDirectory(0);
  delete fMuSF;

  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 11;
  const int histBins = 8;
  TH1D* histo[allPlots][histBins];

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  5) {nBinPlot = 200; xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  6 && thePlot <=  6) {nBinPlot =   7; xminPlot =-0.5; xmaxPlot =   6.5;}
    else if(thePlot >=  7 && thePlot <=  7) {nBinPlot =  40; xminPlot = 0.0; xmaxPlot =  40.0;}
    else if(thePlot >=  8 && thePlot <=  8) {nBinPlot =  40; xminPlot =-0.5; xmaxPlot =  39.5;}
    else if(thePlot >=  9 && thePlot <=  9) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot =  50.0;}
    else if(thePlot >= 10 && thePlot <= 10) {nBinPlot =  50; xminPlot = 0.0; xmaxPlot =   2.5;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) histo[thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    histos->Clear();
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

    Int_t nPassTrigger[12] = {0,0,0,0,0,0,0,0,0,0,0,0};

    Int_t nPassCuts[10] = {0,0,0,0,0,0,0,0,0,0};
    double theMCPrescale = mcPrescale;
    if(infilecatv[ifile] == 0) theMCPrescale = 1.0;
    for (int i=0; i<int(the_input_tree->GetEntries()/theMCPrescale); ++i) {
      the_input_tree->GetEntry(i);
      if(i%100000==0) printf("event %d out of %d\n",i,(int)the_input_tree->GetEntries());

      if(typeSel != 11 && typeSel != 13) assert(0);

      Bool_t passFilter = kFALSE;
      Bool_t passTrigger = kFALSE;
      for (int nt = 0; nt < (int)numtokens; nt++) {
	if(typeSel == 11 &&
          (strcmp(tokens[nt],"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
           strcmp(tokens[nt],"HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
           strcmp(tokens[nt],"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0 ||
           strcmp(tokens[nt],"HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v*") == 0) &&
	  (*eventTrigger.triggerFired)[nt] == 1) passTrigger = kTRUE;
	if(typeSel == 13 &&
	  (strcmp(tokens[nt],"HLT_Mu8_TrkIsoVVL_v*")  == 0 ||              
	   strcmp(tokens[nt],"HLT_Mu17_TrkIsoVVL_v*") == 0 ||              
	   strcmp(tokens[nt],"HLT_Mu24_TrkIsoVVL_v*") == 0 ||              
	   strcmp(tokens[nt],"HLT_Mu34_TrkIsoVVL_v*") == 0) &&           
	  (*eventTrigger.triggerFired)[nt] == 1) passTrigger = kTRUE;
      }
      if(passTrigger == kTRUE &&
         eventLeptons.p4->GetEntriesFast() >= 1 &&
         ((TLorentzVector*)(*eventLeptons.p4)[0])->Pt() > 10 && 
	 (double)((TLorentzVector*)(*eventMet.p4)[0])->Pt() < 30.0) passFilter = kTRUE;

      if(passTrigger == kTRUE) nPassCuts[0]++;
      if(passFilter  == kTRUE) nPassCuts[1]++;
      if(passFilter == kFALSE) continue;

      Bool_t passSel = kFALSE;
      TLorentzVector dilep;double deltaPhiDileptonMet = -999.0; double mtW = -999.0;
      vector<int> idLep; vector<int> idTight;
      for(int nlep=0; nlep<eventLeptons.p4->GetEntriesFast(); nlep++) {
        if(selectIdIsoCut(typeLepSel.Data(),TMath::Abs((int)(*eventLeptons.pdgId)[nlep]),TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Pt()),
	   TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[nlep])->Eta()),(double)(*eventLeptons.iso)[nlep],(int)(*eventLeptons.selBits)[nlep],(double)(*eventLeptons.mva)[nlep]))
	                                                                                               {idTight.push_back(1); idLep.push_back(nlep);}
        else if(((int)(*eventLeptons.selBits)[nlep] & BareLeptons::LepFake)  == BareLeptons::LepFake ) {idTight.push_back(0); idLep.push_back(nlep);}
      }

      vector<int> isGenLep;
      for(unsigned nl=0; nl<idLep.size(); nl++){
        bool isGenLepton = false;
        for(int ngen=0; ngen<eventMonteCarlo.p4->GetEntriesFast(); ngen++) {
          if(TMath::Abs((int)(*eventLeptons.pdgId)[idLep[nl]]) == TMath::Abs((int)(*eventMonteCarlo.pdgId)[ngen]) &&
	    ((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])->DeltaR(*((TLorentzVector*)(*eventMonteCarlo.p4)[ngen])) < 0.3) {
	    isGenLepton = true;
	    break;
	  }
	}
	if(isGenLepton == true) isGenLep.push_back(nl);
      }
      
      vector<int> idJet20,idJet30;
      for(int nj=0; nj<eventJets.p4->GetEntriesFast(); nj++){
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() < 10) continue;
        bool passId = passJetId(fMVACut, (float)(*eventJets.puId)[nj], ((TLorentzVector*)(*eventJets.p4)[nj])->Pt(), TMath::Abs(((TLorentzVector*)(*eventJets.p4)[nj])->Eta()));
        //if(passId == false) continue;        

        Bool_t isLepton = kFALSE;
        for(unsigned int nl=0; nl<idLep.size(); nl++){
          if(((TLorentzVector*)(*eventJets.p4)[nj])->DeltaR(*((TLorentzVector*)(*eventLeptons.p4)[idLep[nl]])) < 0.3) isLepton = kTRUE;
	}
	if(isLepton == kTRUE) continue;

        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 20)  idJet20.push_back(nj);
        if(((TLorentzVector*)(*eventJets.p4)[nj])->Pt() > 30)  idJet30.push_back(nj);
      }

      if     (nsel == 0){ // Z->ll
        if(idLep.size() == 2) nPassCuts[2]++;
	if(idLep.size() == 2 &&
     	  ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 10 && 
     	  ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 10) nPassCuts[3]++;
        if(idLep.size() == 2 &&
     	  ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 10 && 
     	  ((TLorentzVector*)(*eventLeptons.p4)[idLep[1]])->Pt() > 10 &&
          (int)(*eventLeptons.pdgId)[idLep[0]]*(int)(*eventLeptons.pdgId)[idLep[1]] < 0 &&
	  TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == TMath::Abs((int)(*eventLeptons.pdgId)[idLep[1]]) &&
	  TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == typeSel &&
	 (infilecatv[ifile] == 0 || isGenLep.size() == 2)) {
          nPassCuts[4]++;
          dilep = ( ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) ) + ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[1])) ) ); 
          if(TMath::Abs(dilep.M()-91.1876)<15.0) passSel = kTRUE;	
	}
      }
      else if(nsel == 1){ // fake
        if(idLep.size() == 1) nPassCuts[2]++;
	if(idLep.size() == 1 &&
     	  ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 10) nPassCuts[3]++;
        if(idLep.size() == 1 &&
     	  ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 10 &&
	  ((typeSel == 11 && idJet30.size() >= 1) ||(typeSel == 13 && idJet20.size() >= 1)) &&
	  TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == typeSel &&
	 (infilecatv[ifile] == 0 || isGenLep.size() == 1)) {
          dilep = ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) );
          nPassCuts[4]++;
          deltaPhiDileptonMet = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
          mtW = TMath::Sqrt(2.0*dilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiDileptonMet)));
	  if(mtW < 20) passSel = kTRUE;
	}
      }
      else if(nsel == 2){ // W->ln
        if(idLep.size() == 1 && idTight[0] == 1) nPassCuts[2]++;
	if(idLep.size() == 1 && idTight[0] == 1 &&
     	  ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 10) nPassCuts[3]++;
        if(idLep.size() == 1 && idTight[0] == 1 &&
     	  ((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() > 10 &&
	  TMath::Abs((int)(*eventLeptons.pdgId)[idLep[0]]) == typeSel &&
	 (infilecatv[ifile] == 0 || isGenLep.size() == 1)) {
          dilep = ( *(TLorentzVector*)(eventLeptons.p4->At(idLep[0])) );
          nPassCuts[4]++;
          deltaPhiDileptonMet = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
          mtW = TMath::Sqrt(2.0*dilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiDileptonMet)));
	  if(mtW > 50) passSel = kTRUE;
	}
      }

      if(passSel == kTRUE) nPassCuts[5]++;
      if(passSel == kFALSE) continue;

      if(mtW < 0){
        deltaPhiDileptonMet = TMath::Abs(dilep.DeltaPhi(*((TLorentzVector*)(*eventMet.p4)[0])));
        mtW = TMath::Sqrt(2.0*dilep.Pt()*((TLorentzVector*)(*eventMet.p4)[0])->Pt()*(1.0 - cos(deltaPhiDileptonMet)));
      }

      double mcWeight = eventMonteCarlo.mcWeight;
      if(infilecatv[ifile] == 0) mcWeight = 1.0;
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
		typeLepSel.Data(),fhDMuMediumSF,fhDMuIsoSF,fhDElMediumSF,fhDElTightSF,fhDElMediumMVASF,fhDElTightMVASF);
        }
      }
      
      int iPt = -1;
      if     (((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 15){
        iPt = 0;
      }
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 20){
        iPt = 1;
      }
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 25){
        iPt = 2;
      }
      else if(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt() < 30){
        iPt = 3;
      }
      else {
        iPt = 4;
      }
      int iEta = -1;
      if     (TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) < 0.5){
         iEta = 0;
      }
      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) < 1.0){
         iEta = 1;
      }
      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) < 1.5){
         iEta = 2;
      }
      else if(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()) < 2.0){
         iEta = 3;
      }
      else {
         iEta = 4;
      }
      if(nsel == 0) iEta = 0;

      double thePrescale = 1.0;
      if(infilecatv[ifile] != 0 && applyPrescale == 1) thePrescale = prescale[iPt];

      double totalWeight = mcWeight*theLumi*puWeight*effSF*thePrescale*theMCPrescale;

      if(infilecatv[ifile] == 0) {
                            denFRDA[iPt][iEta] = denFRDA[iPt][iEta] + totalWeight;
        if(idTight[0] == 1) numFRDA[iPt][iEta] = numFRDA[iPt][iEta] + totalWeight;
      }
      else {
                            denFRBG[iPt][iEta] = denFRBG[iPt][iEta] + totalWeight;
        if(idTight[0] == 1) numFRBG[iPt][iEta] = numFRBG[iPt][iEta] + totalWeight;
      }

      for(int thePlot=0; thePlot<allPlots; thePlot++){
        double theVar = 0.0;
        if     (thePlot ==  0) theVar = TMath::Min(dilep.M(),199.999);
        else if(thePlot ==  1) theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);
        else if(thePlot ==  2) theVar = TMath::Min((double)eventMet.metNoMu->Pt(),199.999);
        else if(thePlot ==  3) theVar = TMath::Min((double)((TLorentzVector*)(*eventMet.p4)[0])->Pt(),199.999);
        else if(thePlot ==  4) theVar = TMath::Min(dilep.Pt(),199.999);
        else if(thePlot ==  5) theVar = TMath::Min(mtW,199.999);
        else if(thePlot ==  6) theVar = TMath::Min((double)0.0,6.499);
        else if(thePlot ==  7) theVar = TMath::Min((double)eventEvent.rho,39.999);
        else if(thePlot ==  8) theVar = TMath::Min((double)eventVertex.npv,39.499);
        else if(thePlot ==  9) theVar = TMath::Min(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Pt(),49.999);
        else if(thePlot == 10) theVar = TMath::Min(TMath::Abs(((TLorentzVector*)(*eventLeptons.p4)[idLep[0]])->Eta()),2.499);
        histo[thePlot][infilecatv[ifile]]->Fill(theVar,totalWeight);
      }

    }

    printf("eff_cuts: ");
    for(int nc=0; nc<6; nc++){
      double nminusone = the_input_tree->GetEntries();
      if(nc>0) nminusone = nPassCuts[nc-1];
      printf("(%d): %8.5f(%8.5f) | ",nc,100*(double)nPassCuts[nc]/the_input_tree->GetEntries(),100*(double)nPassCuts[nc]/nminusone);
    }
    printf("\n");
  } // end of chain

  for(int nc=0; nc<8; nc++){
    printf("(%d): %5.2f | ",nc,histo[0][nc]->GetSumOfWeights());
  }
  double sumTot[2] = {0.,0.};
  printf("totalMC: %5.2f\n",histo[0][1]->GetSumOfWeights()+histo[0][2]->GetSumOfWeights()+histo[0][3]->GetSumOfWeights()+
                            histo[0][4]->GetSumOfWeights()+histo[0][5]->GetSumOfWeights()+histo[0][6]->GetSumOfWeights()+histo[0][7]->GetSumOfWeights());
  for(int iEta=0; iEta<5; iEta++){
    for(int iPt=0; iPt<5; iPt++){
      sumTot[0] = sumTot[0] + denFRDA[iPt][iEta];
      sumTot[1] = sumTot[1] + denFRBG[iPt][iEta];
      printf("(%d,%d): (%6.1f-%6.1f)/(%6.1f-%6.1f)=%4.3f | ",iPt,iEta,numFRDA[iPt][iEta],numFRBG[iPt][iEta] , denFRDA[iPt][iEta],denFRBG[iPt][iEta],
                                                                     (numFRDA[iPt][iEta]-numFRBG[iPt][iEta])/(denFRDA[iPt][iEta]-denFRBG[iPt][iEta]));
      if(iPt==4) printf("\n");
    }
  }
  printf("sumTot(da/bg) = %f / %f = %f\n",sumTot[0],sumTot[1],sumTot[0]/sumTot[1]);
  if(nsel == 0){
    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
        printf("(%d,%d): (%6.1f)/(%6.1f)=%7.5f | ",iPt,iEta,denFRDA[iPt][iEta],denFRBG[iPt][iEta],denFRDA[iPt][iEta]/denFRBG[iPt][iEta]);
        if(iPt==4) printf("\n");
      }
    }
  }
  else {
    for(int iEta=0; iEta<5; iEta++){
      for(int iPt=0; iPt<5; iPt++){
        printf("(%d,%d): (%6.1f)/(%6.1f)=%4.3f | ",iPt,iEta,numFRDA[iPt][iEta] , denFRDA[iPt][iEta],
                                                           (numFRDA[iPt][iEta])/(denFRDA[iPt][iEta]));
        if(iPt==4) printf("\n");
      }
    }
  }

  if(typeSel == 11) printf("double fake_rate_e[%d][%d] = {\n",5,5);
  else              printf("double fake_rate_m[%d][%d] = {\n",5,5);
  for(int iEta=0; iEta<5; iEta++){
    for(int iPt=0; iPt<5; iPt++){
      printf("%4.3f",TMath::Abs((numFRDA[iPt][iEta]-numFRBG[iPt][iEta])/(denFRDA[iPt][iEta]-denFRBG[iPt][iEta])));
      if(iPt!=4||iEta!=4) printf(",");
      if(iPt==4) printf("\n");
    }
  }
  printf("};\n");
  
  for(int thePlot=0; thePlot<allPlots; thePlot++){
    char output[200];
    sprintf(output,"histo_fake_%d_%d_%d.root",thePlot,nsel,typeSel);	  
    TFile* outFilePlotsNote = new TFile(output,"recreate");
    outFilePlotsNote->cd();
    for(int np=0; np<histBins; np++) histo[thePlot][np]->Write();
    outFilePlotsNote->Close();
  }

}
