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
#include "TRandom2.h"
#include "TLorentzVector.h"

void start2017X(){


  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString>	infilenamev;
  vector<Int_t> 	infilecatv;
  vector<Int_t> 	infilexsecv;
  float                 lumi=41290.;//in pb

  //xSecs in pb!
	
  //infilenamev.push_back("/eos/cms/store/user/jsalfeld/mergedDATAnanoALLSingleEleSingleMuFiltered.root"); infilecatv.push_back(0); infilexsecv.push_back(1);
  //infilenamev.push_back("/eos/cms/store/user/jsalfeld/chekci.root"); infilecatv.push_back(0); infilexsecv.push_back(1);
  
  // infilenamev.push_back("/eos/cms/store/user/jsalfeld/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/mergedX.root"); infilecatv.push_back(3); infilexsecv.push_back(4.715*1.109);
 infilenamev.push_back("/eos/cms/store/user/jsalfeld/GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8/mergedX.root"); infilecatv.push_back(4); infilexsecv.push_back(4.715*1.109);




  double xmin = 0.0;
  double xmax = 1.0;
  int nBinPlot      = 200;
  double xminPlot   = 0.0;
  double xmaxPlot   = 200.0;
  const int allPlots = 3;
  const int histBins = 8;
  const int allStates = 5;

  TH1D* histo[allStates][allPlots][histBins];
  TString processName[histBins] = {"..Data", ".Fakes", "Zgamma", "....WZ", "....ZZ", "...VVV","..EWWZ" ,".Higgs"};

  for(int thePlot=0; thePlot<allPlots; thePlot++){
    if     (thePlot >=  0 && thePlot <=  0) {nBinPlot = 50; 	xminPlot = 0.0; xmaxPlot = 500.0;}
    else if(thePlot >=  1 && thePlot <=  1) {nBinPlot = 40; 	xminPlot = 0.0; xmaxPlot = 200.0;}
    else if(thePlot >=  2 && thePlot <=  2) {nBinPlot = 4; 	xminPlot = -0.5; xmaxPlot = 3.5;}
    TH1D* histos = new TH1D("histos", "histos", nBinPlot, xminPlot, xmaxPlot);
    histos->Sumw2();
    for(int i=0; i<histBins; i++) {
      for(int j=0; j<allStates; j++) histo[j][thePlot][i] = (TH1D*) histos->Clone(Form("histo%d",i));
    }
    histos->Reset();histos->Clear();
  }



  
  //for loop over samples
for(unsigned int ifile=0;ifile<infilenamev.size(); ifile++){
  cout<<infilenamev[ifile]<<endl;
  TFile *the_input_file=TFile::Open((TString)infilenamev[ifile]);
  TTree *the_input_tree = (TTree*)the_input_file->FindObjectAny("Events");
  TTree *the_input_runs  = (TTree*)the_input_file->FindObjectAny("Runs");


  //cross section weight
  float xSecLumiWeight= 1.;
  Long64_t genEventCountn;
  Long64_t genEventCount=1.;
  if(infilecatv[ifile] != 0){
    the_input_runs->SetBranchAddress("genEventCount",&genEventCountn);
    unsigned int numberOfSubSamples = int(the_input_runs->GetEntries());
    cout<<"numberOfSubSamples: "<<numberOfSubSamples<<endl;
    for(unsigned int n=0; n<numberOfSubSamples; n++){
      the_input_runs->GetEntry(n);
      genEventCount = genEventCount + genEventCountn;
    }
    xSecLumiWeight=infilexsecv[ifile]*lumi/genEventCount;
  }

 int numberOfEvents = int(the_input_tree->GetEntries());
 cout<<"number of Events: "<<numberOfEvents<<endl;



 UInt_t nElectron;
 Float_t Electron_pt[nElectron];
 Float_t Electron_eta[nElectron];
 Float_t Electron_phi[nElectron]; 
 Float_t Electron_mass[nElectron]; 
 Int_t Electron_pdgId[nElectron]; 
 Int_t Electron_cutBased[nElectron]; 
 Int_t Electron_cutBasedHLTPreSel[nElectron]; 
 the_input_tree->SetBranchAddress("nElectron",&nElectron); 
 the_input_tree->SetBranchAddress("Electron_pt",Electron_pt); 
 the_input_tree->SetBranchAddress("Electron_eta",Electron_eta); 
 the_input_tree->SetBranchAddress("Electron_phi",Electron_phi); 
 the_input_tree->SetBranchAddress("Electron_mass",Electron_mass); 
 the_input_tree->SetBranchAddress("Electron_pdgId",Electron_pdgId); 
 the_input_tree->SetBranchAddress("Electron_cutBased",Electron_cutBased); 
 the_input_tree->SetBranchAddress("Electron_cutBased_HLTPreSel",Electron_cutBasedHLTPreSel); 

 cout<<"number of Events: "<<numberOfEvents<<endl;
 UInt_t nMuon;
cout<<"number of Events: "<<numberOfEvents<<endl;
 Float_t Muon_pt[nMuon];
cout<<"number of Events: "<<numberOfEvents<<endl;
 Float_t Muon_eta[nMuon];
 Float_t Muon_phi[nMuon]; 
 Float_t Muon_mass[nMuon];
cout<<"number of Events: "<<numberOfEvents<<endl;
 Int_t Muon_pdgId[nMuon]; 
cout<<"number of Events: "<<numberOfEvents<<endl; 
bool  Muon_tightId[nMuon]; 
 
 float Muon_pfRelIso03_all[nMuon]; 
 the_input_tree->SetBranchAddress("nMuon",&nMuon); 
 the_input_tree->SetBranchAddress("Muon_pt",Muon_pt); 
 the_input_tree->SetBranchAddress("Muon_eta",Muon_eta); 
 the_input_tree->SetBranchAddress("Muon_phi",Muon_phi); 
 the_input_tree->SetBranchAddress("Muon_mass",Muon_mass); 
 the_input_tree->SetBranchAddress("Muon_pdgId",Muon_pdgId); 
 the_input_tree->SetBranchAddress("Muon_tightId",Muon_tightId); 
 the_input_tree->SetBranchAddress("Muon_pfRelIso03_all",Muon_pfRelIso03_all); 

cout<<"number of Events: "<<numberOfEvents<<endl;

 float genWeight;
 if(infilecatv[ifile] != 0){
   the_input_tree->SetBranchAddress("genWeight",&genWeight);
 }

//for loop over events for given sample
 
 for (int i=0; i < 10000./*numberOfEvents*/; ++i) 
     {

      the_input_tree->GetEntry(i);
      if(i%100000==0) cout<<"Event: "<<i<<endl;
      vector<int> idLep;  vector<int> idLepPdg; vector<int> idTight; unsigned int goodIsTight = 0; vector<TLorentzVector> leptons;
      
	for(unsigned int e=0; e<nElectron; e++)
	  {
		if(Electron_pt[e]>10 && TMath::Abs(Electron_eta[e])<2.5)
			{
			  if	 (Electron_cutBased[e]         ==4) 	{idLep.push_back(e); idTight.push_back(1); goodIsTight++ ; idLepPdg.push_back(Electron_pdgId[e]);}
	  			else if(Electron_cutBasedHLTPreSel[e]==1) 	{idLep.push_back(e); idTight.push_back(0);               ; idLepPdg.push_back(Electron_pdgId[e]);}
	  			else {continue;}
	  			TLorentzVector lepTemp; lepTemp.SetPtEtaPhiM(Electron_pt[e],Electron_eta[e],Electron_phi[e],Electron_mass[e]);
	  			leptons.push_back(lepTemp);	
	}
      }

      for(unsigned int m=0; m<nMuon; m++)
      	{
	  cout<<i<<"   "<<Muon_pt[m]<<"  "<<Muon_pdgId[m]<<"   "<<m<<"   "<<nMuon<<endl;
		if(Muon_pt[m]>10 && TMath::Abs(Muon_eta[m])<2.4)
	  		{
	  			if     (Muon_tightId[m] && Muon_pfRelIso03_all[m] < 0.15) 	{idLep.push_back(m); idTight.push_back(1); goodIsTight++ ; idLepPdg.push_back(Muon_pdgId[m]);}
	  			else if(Muon_tightId[m] && Muon_pfRelIso03_all[m] < 0.4) 	{idLep.push_back(m); idTight.push_back(0);               ; idLepPdg.push_back(Muon_pdgId[m]);}
	  			else {continue;}
	  			TLorentzVector lepTemp; lepTemp.SetPtEtaPhiM(Muon_pt[m],Muon_eta[m],Muon_phi[m],Muon_mass[m]);
	  			leptons.push_back(lepTemp);	
	}
      }

      //Only 3 leptons
      if(leptons.size()!=3) continue;
      //lepton pT Cuts
      if(leptons.at(0).Pt() < 20. || leptons.at(1).Pt() < 20. || leptons.at(2).Pt() < 20. ) continue;
      //cout<<"  "<<idLepPdg[0]<<"  "<<idLepPdg[1]<<"  "<<idLepPdg[2]<<endl;
      //Find Z boson, categorize events
      double minMassll = 999.0;
      double minMassZ = 999.0;
      double mass3l = 0.0;
      double deltaRllMin = 999.0;
      int type3l = 0;
      int tagZ[3] = {-1,-1,-1};
      for(unsigned nl0=0; nl0<idLep.size()-1; nl0++){
        for(unsigned nl1=nl0+1; nl1<idLep.size(); nl1++){
      	  if((int)idLepPdg[nl0] * (int)idLepPdg[nl1] > 0)	 continue;
          TLorentzVector dilepAux( leptons.at(nl0) + (leptons.at(nl1)) );
	  
	  if(TMath::Abs((int)idLepPdg[nl0])==TMath::Abs((int)idLepPdg[nl1]) &&
	     TMath::Abs(dilepAux.M()-91.1876) < TMath::Abs(minMassZ-91.1876)) {
	     minMassZ = dilepAux.M();tagZ[0]=nl0;tagZ[1]=nl1;
	     //cout<<TMath::Abs((int)idLepPdg[nl0])<<endl;
	     if(TMath::Abs((int)idLepPdg[nl0]) == 13) 			  type3l = 0;
	     else                                                         type3l = 1;
	  }
	  if(minMassll > dilepAux.M()) minMassll = dilepAux.M();
        }
      }
      //cout<<"type3l 1: "<<type3l<<"  "<<idLepPdg[0]<<"  "<<idLepPdg[1]<<"  "<<idLepPdg[2]<<endl;
      for(unsigned nl0=0; nl0<idLep.size(); nl0++){
        if((int)nl0==tagZ[0]||(int)nl0==tagZ[1]) 		 continue;
        tagZ[2] = nl0;
        break;
      }  
      if(tagZ[0] == -1 || tagZ[1] == -1 || tagZ[2] == -1)	 continue;
      if(TMath::Abs((int)idLepPdg[tagZ[2]]) == 13) 		       type3l += 0;
      else							       type3l += 2;
      TLorentzVector trilep(	( ( leptons.at(0)) )
        		       +( ( leptons.at(1)) )
        		       +( ( leptons.at(2)) ));


      //cout<<"type3l 2: "<<type3l<<endl;

      int theCategory = infilecatv[ifile];

      double totalWeight = 1.;
      if(theCategory != 0) totalWeight = xSecLumiWeight*genWeight;


      for(int thePlot=0; thePlot<allPlots; thePlot++){

		double theVar = 0.0;
		bool makePlot = false;
 		if     (thePlot ==  0 )                          				{makePlot = true;theVar = TMath::Min(leptons.at(0).Pt(),1199.999);}
		else if(thePlot ==  1 )                      				 	{makePlot = true;theVar = TMath::Min(leptons.at(1).Pt(),199.999);}
		else if(thePlot ==  2 )                        					{makePlot = true;theVar =  type3l;}
		if( (makePlot && theCategory<7) || (makePlot && theCategory == 7 && infilenamev[ifile].Contains("M900")) ) {histo[     4][thePlot][theCategory]->Fill(theVar,totalWeight);}
     		 //if( (makePlot && theCategory<6) || (makePlot && theCategory == 6 && infilenamev[ifile].Contains("M900")) ) histo[     4][thePlot][theCategory]->Fill(theVar,totalWeight);
		if( (makePlot && theCategory<7) || (makePlot && theCategory == 7 && infilenamev[ifile].Contains("M900")) ) histo[type3l][thePlot][theCategory]->Fill(theVar,totalWeight);
      	  }


     }//end of tree per file
  
 }//end of all files

for(int thePlot=0; thePlot<allPlots; thePlot++){
    for(int j=0; j<allStates; j++){
      char output[200];
      sprintf(output,"histowz_nice_%d_%d.root",j,thePlot);	  
      TFile* outFilePlotsNote = new TFile(output,"recreate");
      outFilePlotsNote->cd();
      for(int np=0; np<histBins; np++) histo[j][thePlot][np]->Write();
      outFilePlotsNote->Close();
    }
  }


}
