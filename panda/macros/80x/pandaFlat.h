//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 29 19:14:55 2017 by ROOT version 6.06/01
// from TTree events/events
// found on file: /afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/finerptbin/DYJetsToLL_M-50_NLO.root
//////////////////////////////////////////////////////////

#ifndef pandaFlat_h
#define pandaFlat_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class pandaFlat {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         scale[6];
   Float_t         sf_btag0;
   Float_t         sf_btag1;
   Float_t         sf_btag2;
   Float_t         sf_btagGT0;
   Float_t         sf_btag0BUp;
   Float_t         sf_btag1BUp;
   Float_t         sf_btag2BUp;
   Float_t         sf_btagGT0BUp;
   Float_t         sf_btag0BDown;
   Float_t         sf_btag1BDown;
   Float_t         sf_btag2BDown;
   Float_t         sf_btagGT0BDown;
   Float_t         sf_btag0MUp;
   Float_t         sf_btag1MUp;
   Float_t         sf_btag2MUp;
   Float_t         sf_btagGT0MUp;
   Float_t         sf_btag0MDown;
   Float_t         sf_btag1MDown;
   Float_t         sf_btag2MDown;
   Float_t         sf_btagGT0MDown;
   Int_t           runNumber;
   Int_t           lumiNumber;
   ULong64_t       eventNumber;
   Int_t           npv;
   Int_t           pu;
   Float_t         mcWeight;
   Int_t           trigger;
   Int_t           metFilter;
   Int_t           egmFilter;
   Int_t           nLooseLep;
   Int_t           looseGenLep1PdgId;
   Int_t           looseGenLep2PdgId;
   Int_t           looseGenLep3PdgId;
   Int_t           looseGenLep4PdgId;
   Int_t           looseLep1PdgId;
   Int_t           looseLep2PdgId;
   Int_t           looseLep3PdgId;
   Int_t           looseLep4PdgId;
   Int_t           looseLep1SelBit;
   Int_t           looseLep2SelBit;
   Int_t           looseLep3SelBit;
   Int_t           looseLep4SelBit;
   Float_t         looseLep1Pt;
   Float_t         looseLep2Pt;
   Float_t         looseLep3Pt;
   Float_t         looseLep4Pt;
   Float_t         looseLep1Eta;
   Float_t         looseLep2Eta;
   Float_t         looseLep3Eta;
   Float_t         looseLep4Eta;
   Float_t         looseLep1Phi;
   Float_t         looseLep2Phi;
   Float_t         looseLep3Phi;
   Float_t         looseLep4Phi;
   Int_t           nJet;
   Int_t           jetNLBtags;
   Int_t           jetNMBtags;
   Int_t           jetNTBtags;
   Float_t         jet1Pt;
   Float_t         jet2Pt;
   Float_t         jet3Pt;
   Float_t         jet4Pt;
   Float_t         jet1Eta;
   Float_t         jet2Eta;
   Float_t         jet3Eta;
   Float_t         jet4Eta;
   Float_t         jet1Phi;
   Float_t         jet2Phi;
   Float_t         jet3Phi;
   Float_t         jet4Phi;
   Float_t         jet1BTag;
   Float_t         jet2BTag;
   Float_t         jet3BTag;
   Float_t         jet4BTag;
   Float_t         jet1GenPt;
   Float_t         jet2GenPt;
   Float_t         jet3GenPt;
   Float_t         jet4GenPt;
   Int_t           jet1Flav;
   Int_t           jet2Flav;
   Int_t           jet3Flav;
   Int_t           jet4Flav;
   Int_t           jet1SelBit;
   Int_t           jet2SelBit;
   Int_t           jet3SelBit;
   Int_t           jet4SelBit;
   Float_t         jet1PtUp;
   Float_t         jet2PtUp;
   Float_t         jet3PtUp;
   Float_t         jet4PtUp;
   Float_t         jet1PtDown;
   Float_t         jet2PtDown;
   Float_t         jet3PtDown;
   Float_t         jet4PtDown;
   Float_t         jet1EtaUp;
   Float_t         jet2EtaUp;
   Float_t         jet3EtaUp;
   Float_t         jet4EtaUp;
   Float_t         jet1EtaDown;
   Float_t         jet2EtaDown;
   Float_t         jet3EtaDown;
   Float_t         jet4EtaDown;
   Float_t         pfmet;
   Float_t         pfmetphi;
   Float_t         pfmetRaw;
   Float_t         pfmetUp;
   Float_t         pfmetDown;
   Float_t         pfmetnomu;
   Float_t         puppimet;
   Float_t         puppimetphi;
   Float_t         calomet;
   Float_t         calometphi;
   Float_t         trkmet;
   Float_t         trkmetphi;
   Float_t         dphipfmet;
   Float_t         dphipuppimet;
   Float_t         genLep1Pt;
   Float_t         genLep2Pt;
   Float_t         genLep1Eta;
   Float_t         genLep2Eta;
   Float_t         genLep1Phi;
   Float_t         genLep2Phi;
   Int_t           genLep1PdgId;
   Int_t           genLep2PdgId;
   Int_t           nTau;
   Float_t         pdfUp;
   Float_t         pdfDown;
   Int_t           nLoosePhoton;
   Float_t         loosePho1Pt;
   Float_t         loosePho1Eta;
   Float_t         loosePho1Phi;
   Float_t         sf_pu;
   Float_t         sf_puUp;
   Float_t         sf_puDown;
   Float_t         sf_zz;
   Float_t         sf_zzUnc;
   Float_t         sf_wz;
   Float_t         sf_zh;
   Float_t         sf_zhUp;
   Float_t         sf_zhDown;
   Float_t         sf_tt;
   Float_t         sf_trk1;
   Float_t         sf_trk2;
   Float_t         sf_trk3;
   Float_t         sf_trk4;
   Float_t         sf_loose1;
   Float_t         sf_loose2;
   Float_t         sf_loose3;
   Float_t         sf_loose4;
   Float_t         sf_medium1;
   Float_t         sf_medium2;
   Float_t         sf_medium3;
   Float_t         sf_medium4;
   Float_t         sf_tight1;
   Float_t         sf_tight2;
   Float_t         sf_tight3;
   Float_t         sf_tight4;
   Float_t         sf_unc1;
   Float_t         sf_unc2;
   Float_t         sf_unc3;
   Float_t         sf_unc4;
   Float_t         normalizedWeight;

   // List of branches
   TBranch        *b_scale;   //!
   TBranch        *b_sf_btag0;   //!
   TBranch        *b_sf_btag1;   //!
   TBranch        *b_sf_btag2;   //!
   TBranch        *b_sf_btagGT0;   //!
   TBranch        *b_sf_btag0BUp;   //!
   TBranch        *b_sf_btag1BUp;   //!
   TBranch        *b_sf_btag2BUp;   //!
   TBranch        *b_sf_btagGT0BUp;   //!
   TBranch        *b_sf_btag0BDown;   //!
   TBranch        *b_sf_btag1BDown;   //!
   TBranch        *b_sf_btag2BDown;   //!
   TBranch        *b_sf_btagGT0BDown;   //!
   TBranch        *b_sf_btag0MUp;   //!
   TBranch        *b_sf_btag1MUp;   //!
   TBranch        *b_sf_btag2MUp;   //!
   TBranch        *b_sf_btagGT0MUp;   //!
   TBranch        *b_sf_btag0MDown;   //!
   TBranch        *b_sf_btag1MDown;   //!
   TBranch        *b_sf_btag2MDown;   //!
   TBranch        *b_sf_btagGT0MDown;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_pu;   //!
   TBranch        *b_mcWeight;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_metFilter;   //!
   TBranch        *b_egmFilter;   //!
   TBranch        *b_nLooseLep;   //!
   TBranch        *b_looseGenLep1PdgId;   //!
   TBranch        *b_looseGenLep2PdgId;   //!
   TBranch        *b_looseGenLep3PdgId;   //!
   TBranch        *b_looseGenLep4PdgId;   //!
   TBranch        *b_looseLep1PdgId;   //!
   TBranch        *b_looseLep2PdgId;   //!
   TBranch        *b_looseLep3PdgId;   //!
   TBranch        *b_looseLep4PdgId;   //!
   TBranch        *b_looseLep1SelBit;   //!
   TBranch        *b_looseLep2SelBit;   //!
   TBranch        *b_looseLep3SelBit;   //!
   TBranch        *b_looseLep4SelBit;   //!
   TBranch        *b_looseLep1Pt;   //!
   TBranch        *b_looseLep2Pt;   //!
   TBranch        *b_looseLep3Pt;   //!
   TBranch        *b_looseLep4Pt;   //!
   TBranch        *b_looseLep1Eta;   //!
   TBranch        *b_looseLep2Eta;   //!
   TBranch        *b_looseLep3Eta;   //!
   TBranch        *b_looseLep4Eta;   //!
   TBranch        *b_looseLep1Phi;   //!
   TBranch        *b_looseLep2Phi;   //!
   TBranch        *b_looseLep3Phi;   //!
   TBranch        *b_looseLep4Phi;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetNLBtags;   //!
   TBranch        *b_jetNMBtags;   //!
   TBranch        *b_jetNTBtags;   //!
   TBranch        *b_jet1Pt;   //!
   TBranch        *b_jet2Pt;   //!
   TBranch        *b_jet3Pt;   //!
   TBranch        *b_jet4Pt;   //!
   TBranch        *b_jet1Eta;   //!
   TBranch        *b_jet2Eta;   //!
   TBranch        *b_jet3Eta;   //!
   TBranch        *b_jet4Eta;   //!
   TBranch        *b_jet1Phi;   //!
   TBranch        *b_jet2Phi;   //!
   TBranch        *b_jet3Phi;   //!
   TBranch        *b_jet4Phi;   //!
   TBranch        *b_jet1BTag;   //!
   TBranch        *b_jet2BTag;   //!
   TBranch        *b_jet3BTag;   //!
   TBranch        *b_jet4BTag;   //!
   TBranch        *b_jet1GenPt;   //!
   TBranch        *b_jet2GenPt;   //!
   TBranch        *b_jet3GenPt;   //!
   TBranch        *b_jet4GenPt;   //!
   TBranch        *b_jet1Flav;   //!
   TBranch        *b_jet2Flav;   //!
   TBranch        *b_jet3Flav;   //!
   TBranch        *b_jet4Flav;   //!
   TBranch        *b_jet1SelBit;   //!
   TBranch        *b_jet2SelBit;   //!
   TBranch        *b_jet3SelBit;   //!
   TBranch        *b_jet4SelBit;   //!
   TBranch        *b_jet1PtUp;   //!
   TBranch        *b_jet2PtUp;   //!
   TBranch        *b_jet3PtUp;   //!
   TBranch        *b_jet4PtUp;   //!
   TBranch        *b_jet1PtDown;   //!
   TBranch        *b_jet2PtDown;   //!
   TBranch        *b_jet3PtDown;   //!
   TBranch        *b_jet4PtDown;   //!
   TBranch        *b_jet1EtaUp;   //!
   TBranch        *b_jet2EtaUp;   //!
   TBranch        *b_jet3EtaUp;   //!
   TBranch        *b_jet4EtaUp;   //!
   TBranch        *b_jet1EtaDown;   //!
   TBranch        *b_jet2EtaDown;   //!
   TBranch        *b_jet3EtaDown;   //!
   TBranch        *b_jet4EtaDown;   //!
   TBranch        *b_pfmet;   //!
   TBranch        *b_pfmetphi;   //!
   TBranch        *b_pfmetRaw;   //!
   TBranch        *b_pfmetUp;   //!
   TBranch        *b_pfmetDown;   //!
   TBranch        *b_pfmetnomu;   //!
   TBranch        *b_puppimet;   //!
   TBranch        *b_puppimetphi;   //!
   TBranch        *b_calomet;   //!
   TBranch        *b_calometphi;   //!
   TBranch        *b_trkmet;   //!
   TBranch        *b_trkmetphi;   //!
   TBranch        *b_dphipfmet;   //!
   TBranch        *b_dphipuppimet;   //!
   TBranch        *b_genLep1Pt;   //!
   TBranch        *b_genLep2Pt;   //!
   TBranch        *b_genLep1Eta;   //!
   TBranch        *b_genLep2Eta;   //!
   TBranch        *b_genLep1Phi;   //!
   TBranch        *b_genLep2Phi;   //!
   TBranch        *b_genLep1PdgId;   //!
   TBranch        *b_genLep2PdgId;   //!
   TBranch        *b_nTau;   //!
   TBranch        *b_pdfUp;   //!
   TBranch        *b_pdfDown;   //!
   TBranch        *b_nLoosePhoton;   //!
   TBranch        *b_loosePho1Pt;   //!
   TBranch        *b_loosePho1Eta;   //!
   TBranch        *b_loosePho1Phi;   //!
   TBranch        *b_sf_pu;   //!
   TBranch        *b_sf_puUp;   //!
   TBranch        *b_sf_puDown;   //!
   TBranch        *b_sf_zz;   //!
   TBranch        *b_sf_zzUnc;   //!
   TBranch        *b_sf_wz;   //!
   TBranch        *b_sf_zh;   //!
   TBranch        *b_sf_zhUp;   //!
   TBranch        *b_sf_zhDown;   //!
   TBranch        *b_sf_tt;   //!
   TBranch        *b_sf_trk1;   //!
   TBranch        *b_sf_trk2;   //!
   TBranch        *b_sf_trk3;   //!
   TBranch        *b_sf_trk4;   //!
   TBranch        *b_sf_loose1;   //!
   TBranch        *b_sf_loose2;   //!
   TBranch        *b_sf_loose3;   //!
   TBranch        *b_sf_loose4;   //!
   TBranch        *b_sf_medium1;   //!
   TBranch        *b_sf_medium2;   //!
   TBranch        *b_sf_medium3;   //!
   TBranch        *b_sf_medium4;   //!
   TBranch        *b_sf_tight1;   //!
   TBranch        *b_sf_tight2;   //!
   TBranch        *b_sf_tight3;   //!
   TBranch        *b_sf_tight4;   //!
   TBranch        *b_sf_unc1;   //!
   TBranch        *b_sf_unc2;   //!
   TBranch        *b_sf_unc3;   //!
   TBranch        *b_sf_unc4;   //!
   TBranch        *b_normalizedWeight;   //!

   pandaFlat(TTree *tree=0);
   virtual ~pandaFlat();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pandaFlat_cxx
pandaFlat::pandaFlat(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/finerptbin/DYJetsToLL_M-50_NLO.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/finerptbin/DYJetsToLL_M-50_NLO.root");
      }
      f->GetObject("events",tree);

   }
   Init(tree);
}

pandaFlat::~pandaFlat()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pandaFlat::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pandaFlat::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void pandaFlat::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("scale", scale, &b_scale);
   fChain->SetBranchAddress("sf_btag0", &sf_btag0, &b_sf_btag0);
   fChain->SetBranchAddress("sf_btag1", &sf_btag1, &b_sf_btag1);
   fChain->SetBranchAddress("sf_btag2", &sf_btag2, &b_sf_btag2);
   fChain->SetBranchAddress("sf_btagGT0", &sf_btagGT0, &b_sf_btagGT0);
   fChain->SetBranchAddress("sf_btag0BUp", &sf_btag0BUp, &b_sf_btag0BUp);
   fChain->SetBranchAddress("sf_btag1BUp", &sf_btag1BUp, &b_sf_btag1BUp);
   fChain->SetBranchAddress("sf_btag2BUp", &sf_btag2BUp, &b_sf_btag2BUp);
   fChain->SetBranchAddress("sf_btagGT0BUp", &sf_btagGT0BUp, &b_sf_btagGT0BUp);
   fChain->SetBranchAddress("sf_btag0BDown", &sf_btag0BDown, &b_sf_btag0BDown);
   fChain->SetBranchAddress("sf_btag1BDown", &sf_btag1BDown, &b_sf_btag1BDown);
   fChain->SetBranchAddress("sf_btag2BDown", &sf_btag2BDown, &b_sf_btag2BDown);
   fChain->SetBranchAddress("sf_btagGT0BDown", &sf_btagGT0BDown, &b_sf_btagGT0BDown);
   fChain->SetBranchAddress("sf_btag0MUp", &sf_btag0MUp, &b_sf_btag0MUp);
   fChain->SetBranchAddress("sf_btag1MUp", &sf_btag1MUp, &b_sf_btag1MUp);
   fChain->SetBranchAddress("sf_btag2MUp", &sf_btag2MUp, &b_sf_btag2MUp);
   fChain->SetBranchAddress("sf_btagGT0MUp", &sf_btagGT0MUp, &b_sf_btagGT0MUp);
   fChain->SetBranchAddress("sf_btag0MDown", &sf_btag0MDown, &b_sf_btag0MDown);
   fChain->SetBranchAddress("sf_btag1MDown", &sf_btag1MDown, &b_sf_btag1MDown);
   fChain->SetBranchAddress("sf_btag2MDown", &sf_btag2MDown, &b_sf_btag2MDown);
   fChain->SetBranchAddress("sf_btagGT0MDown", &sf_btagGT0MDown, &b_sf_btagGT0MDown);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("lumiNumber", &lumiNumber, &b_lumiNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("pu", &pu, &b_pu);
   fChain->SetBranchAddress("mcWeight", &mcWeight, &b_mcWeight);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("metFilter", &metFilter, &b_metFilter);
   fChain->SetBranchAddress("egmFilter", &egmFilter, &b_egmFilter);
   fChain->SetBranchAddress("nLooseLep", &nLooseLep, &b_nLooseLep);
   fChain->SetBranchAddress("looseGenLep1PdgId", &looseGenLep1PdgId, &b_looseGenLep1PdgId);
   fChain->SetBranchAddress("looseGenLep2PdgId", &looseGenLep2PdgId, &b_looseGenLep2PdgId);
   fChain->SetBranchAddress("looseGenLep3PdgId", &looseGenLep3PdgId, &b_looseGenLep3PdgId);
   fChain->SetBranchAddress("looseGenLep4PdgId", &looseGenLep4PdgId, &b_looseGenLep4PdgId);
   fChain->SetBranchAddress("looseLep1PdgId", &looseLep1PdgId, &b_looseLep1PdgId);
   fChain->SetBranchAddress("looseLep2PdgId", &looseLep2PdgId, &b_looseLep2PdgId);
   fChain->SetBranchAddress("looseLep3PdgId", &looseLep3PdgId, &b_looseLep3PdgId);
   fChain->SetBranchAddress("looseLep4PdgId", &looseLep4PdgId, &b_looseLep4PdgId);
   fChain->SetBranchAddress("looseLep1SelBit", &looseLep1SelBit, &b_looseLep1SelBit);
   fChain->SetBranchAddress("looseLep2SelBit", &looseLep2SelBit, &b_looseLep2SelBit);
   fChain->SetBranchAddress("looseLep3SelBit", &looseLep3SelBit, &b_looseLep3SelBit);
   fChain->SetBranchAddress("looseLep4SelBit", &looseLep4SelBit, &b_looseLep4SelBit);
   fChain->SetBranchAddress("looseLep1Pt", &looseLep1Pt, &b_looseLep1Pt);
   fChain->SetBranchAddress("looseLep2Pt", &looseLep2Pt, &b_looseLep2Pt);
   fChain->SetBranchAddress("looseLep3Pt", &looseLep3Pt, &b_looseLep3Pt);
   fChain->SetBranchAddress("looseLep4Pt", &looseLep4Pt, &b_looseLep4Pt);
   fChain->SetBranchAddress("looseLep1Eta", &looseLep1Eta, &b_looseLep1Eta);
   fChain->SetBranchAddress("looseLep2Eta", &looseLep2Eta, &b_looseLep2Eta);
   fChain->SetBranchAddress("looseLep3Eta", &looseLep3Eta, &b_looseLep3Eta);
   fChain->SetBranchAddress("looseLep4Eta", &looseLep4Eta, &b_looseLep4Eta);
   fChain->SetBranchAddress("looseLep1Phi", &looseLep1Phi, &b_looseLep1Phi);
   fChain->SetBranchAddress("looseLep2Phi", &looseLep2Phi, &b_looseLep2Phi);
   fChain->SetBranchAddress("looseLep3Phi", &looseLep3Phi, &b_looseLep3Phi);
   fChain->SetBranchAddress("looseLep4Phi", &looseLep4Phi, &b_looseLep4Phi);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetNLBtags", &jetNLBtags, &b_jetNLBtags);
   fChain->SetBranchAddress("jetNMBtags", &jetNMBtags, &b_jetNMBtags);
   fChain->SetBranchAddress("jetNTBtags", &jetNTBtags, &b_jetNTBtags);
   fChain->SetBranchAddress("jet1Pt", &jet1Pt, &b_jet1Pt);
   fChain->SetBranchAddress("jet2Pt", &jet2Pt, &b_jet2Pt);
   fChain->SetBranchAddress("jet3Pt", &jet3Pt, &b_jet3Pt);
   fChain->SetBranchAddress("jet4Pt", &jet4Pt, &b_jet4Pt);
   fChain->SetBranchAddress("jet1Eta", &jet1Eta, &b_jet1Eta);
   fChain->SetBranchAddress("jet2Eta", &jet2Eta, &b_jet2Eta);
   fChain->SetBranchAddress("jet3Eta", &jet3Eta, &b_jet3Eta);
   fChain->SetBranchAddress("jet4Eta", &jet4Eta, &b_jet4Eta);
   fChain->SetBranchAddress("jet1Phi", &jet1Phi, &b_jet1Phi);
   fChain->SetBranchAddress("jet2Phi", &jet2Phi, &b_jet2Phi);
   fChain->SetBranchAddress("jet3Phi", &jet3Phi, &b_jet3Phi);
   fChain->SetBranchAddress("jet4Phi", &jet4Phi, &b_jet4Phi);
   fChain->SetBranchAddress("jet1BTag", &jet1BTag, &b_jet1BTag);
   fChain->SetBranchAddress("jet2BTag", &jet2BTag, &b_jet2BTag);
   fChain->SetBranchAddress("jet3BTag", &jet3BTag, &b_jet3BTag);
   fChain->SetBranchAddress("jet4BTag", &jet4BTag, &b_jet4BTag);
   fChain->SetBranchAddress("jet1GenPt", &jet1GenPt, &b_jet1GenPt);
   fChain->SetBranchAddress("jet2GenPt", &jet2GenPt, &b_jet2GenPt);
   fChain->SetBranchAddress("jet3GenPt", &jet3GenPt, &b_jet3GenPt);
   fChain->SetBranchAddress("jet4GenPt", &jet4GenPt, &b_jet4GenPt);
   fChain->SetBranchAddress("jet1Flav", &jet1Flav, &b_jet1Flav);
   fChain->SetBranchAddress("jet2Flav", &jet2Flav, &b_jet2Flav);
   fChain->SetBranchAddress("jet3Flav", &jet3Flav, &b_jet3Flav);
   fChain->SetBranchAddress("jet4Flav", &jet4Flav, &b_jet4Flav);
   fChain->SetBranchAddress("jet1SelBit", &jet1SelBit, &b_jet1SelBit);
   fChain->SetBranchAddress("jet2SelBit", &jet2SelBit, &b_jet2SelBit);
   fChain->SetBranchAddress("jet3SelBit", &jet3SelBit, &b_jet3SelBit);
   fChain->SetBranchAddress("jet4SelBit", &jet4SelBit, &b_jet4SelBit);
   fChain->SetBranchAddress("jet1PtUp", &jet1PtUp, &b_jet1PtUp);
   fChain->SetBranchAddress("jet2PtUp", &jet2PtUp, &b_jet2PtUp);
   fChain->SetBranchAddress("jet3PtUp", &jet3PtUp, &b_jet3PtUp);
   fChain->SetBranchAddress("jet4PtUp", &jet4PtUp, &b_jet4PtUp);
   fChain->SetBranchAddress("jet1PtDown", &jet1PtDown, &b_jet1PtDown);
   fChain->SetBranchAddress("jet2PtDown", &jet2PtDown, &b_jet2PtDown);
   fChain->SetBranchAddress("jet3PtDown", &jet3PtDown, &b_jet3PtDown);
   fChain->SetBranchAddress("jet4PtDown", &jet4PtDown, &b_jet4PtDown);
   fChain->SetBranchAddress("jet1EtaUp", &jet1EtaUp, &b_jet1EtaUp);
   fChain->SetBranchAddress("jet2EtaUp", &jet2EtaUp, &b_jet2EtaUp);
   fChain->SetBranchAddress("jet3EtaUp", &jet3EtaUp, &b_jet3EtaUp);
   fChain->SetBranchAddress("jet4EtaUp", &jet4EtaUp, &b_jet4EtaUp);
   fChain->SetBranchAddress("jet1EtaDown", &jet1EtaDown, &b_jet1EtaDown);
   fChain->SetBranchAddress("jet2EtaDown", &jet2EtaDown, &b_jet2EtaDown);
   fChain->SetBranchAddress("jet3EtaDown", &jet3EtaDown, &b_jet3EtaDown);
   fChain->SetBranchAddress("jet4EtaDown", &jet4EtaDown, &b_jet4EtaDown);
   fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
   fChain->SetBranchAddress("pfmetphi", &pfmetphi, &b_pfmetphi);
   fChain->SetBranchAddress("pfmetRaw", &pfmetRaw, &b_pfmetRaw);
   fChain->SetBranchAddress("pfmetUp", &pfmetUp, &b_pfmetUp);
   fChain->SetBranchAddress("pfmetDown", &pfmetDown, &b_pfmetDown);
   fChain->SetBranchAddress("pfmetnomu", &pfmetnomu, &b_pfmetnomu);
   fChain->SetBranchAddress("puppimet", &puppimet, &b_puppimet);
   fChain->SetBranchAddress("puppimetphi", &puppimetphi, &b_puppimetphi);
   fChain->SetBranchAddress("calomet", &calomet, &b_calomet);
   fChain->SetBranchAddress("calometphi", &calometphi, &b_calometphi);
   fChain->SetBranchAddress("trkmet", &trkmet, &b_trkmet);
   fChain->SetBranchAddress("trkmetphi", &trkmetphi, &b_trkmetphi);
   fChain->SetBranchAddress("dphipfmet", &dphipfmet, &b_dphipfmet);
   fChain->SetBranchAddress("dphipuppimet", &dphipuppimet, &b_dphipuppimet);
   fChain->SetBranchAddress("genLep1Pt", &genLep1Pt, &b_genLep1Pt);
   fChain->SetBranchAddress("genLep2Pt", &genLep2Pt, &b_genLep2Pt);
   fChain->SetBranchAddress("genLep1Eta", &genLep1Eta, &b_genLep1Eta);
   fChain->SetBranchAddress("genLep2Eta", &genLep2Eta, &b_genLep2Eta);
   fChain->SetBranchAddress("genLep1Phi", &genLep1Phi, &b_genLep1Phi);
   fChain->SetBranchAddress("genLep2Phi", &genLep2Phi, &b_genLep2Phi);
   fChain->SetBranchAddress("genLep1PdgId", &genLep1PdgId, &b_genLep1PdgId);
   fChain->SetBranchAddress("genLep2PdgId", &genLep2PdgId, &b_genLep2PdgId);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("pdfUp", &pdfUp, &b_pdfUp);
   fChain->SetBranchAddress("pdfDown", &pdfDown, &b_pdfDown);
   fChain->SetBranchAddress("nLoosePhoton", &nLoosePhoton, &b_nLoosePhoton);
   fChain->SetBranchAddress("loosePho1Pt", &loosePho1Pt, &b_loosePho1Pt);
   fChain->SetBranchAddress("loosePho1Eta", &loosePho1Eta, &b_loosePho1Eta);
   fChain->SetBranchAddress("loosePho1Phi", &loosePho1Phi, &b_loosePho1Phi);
   fChain->SetBranchAddress("sf_pu", &sf_pu, &b_sf_pu);
   fChain->SetBranchAddress("sf_puUp", &sf_puUp, &b_sf_puUp);
   fChain->SetBranchAddress("sf_puDown", &sf_puDown, &b_sf_puDown);
   fChain->SetBranchAddress("sf_zz", &sf_zz, &b_sf_zz);
   fChain->SetBranchAddress("sf_zzUnc", &sf_zzUnc, &b_sf_zzUnc);
   fChain->SetBranchAddress("sf_wz", &sf_wz, &b_sf_wz);
   fChain->SetBranchAddress("sf_zh", &sf_zh, &b_sf_zh);
   fChain->SetBranchAddress("sf_zhUp", &sf_zhUp, &b_sf_zhUp);
   fChain->SetBranchAddress("sf_zhDown", &sf_zhDown, &b_sf_zhDown);
   fChain->SetBranchAddress("sf_tt", &sf_tt, &b_sf_tt);
   fChain->SetBranchAddress("sf_trk1", &sf_trk1, &b_sf_trk1);
   fChain->SetBranchAddress("sf_trk2", &sf_trk2, &b_sf_trk2);
   fChain->SetBranchAddress("sf_trk3", &sf_trk3, &b_sf_trk3);
   fChain->SetBranchAddress("sf_trk4", &sf_trk4, &b_sf_trk4);
   fChain->SetBranchAddress("sf_loose1", &sf_loose1, &b_sf_loose1);
   fChain->SetBranchAddress("sf_loose2", &sf_loose2, &b_sf_loose2);
   fChain->SetBranchAddress("sf_loose3", &sf_loose3, &b_sf_loose3);
   fChain->SetBranchAddress("sf_loose4", &sf_loose4, &b_sf_loose4);
   fChain->SetBranchAddress("sf_medium1", &sf_medium1, &b_sf_medium1);
   fChain->SetBranchAddress("sf_medium2", &sf_medium2, &b_sf_medium2);
   fChain->SetBranchAddress("sf_medium3", &sf_medium3, &b_sf_medium3);
   fChain->SetBranchAddress("sf_medium4", &sf_medium4, &b_sf_medium4);
   fChain->SetBranchAddress("sf_tight1", &sf_tight1, &b_sf_tight1);
   fChain->SetBranchAddress("sf_tight2", &sf_tight2, &b_sf_tight2);
   fChain->SetBranchAddress("sf_tight3", &sf_tight3, &b_sf_tight3);
   fChain->SetBranchAddress("sf_tight4", &sf_tight4, &b_sf_tight4);
   fChain->SetBranchAddress("sf_unc1", &sf_unc1, &b_sf_unc1);
   fChain->SetBranchAddress("sf_unc2", &sf_unc2, &b_sf_unc2);
   fChain->SetBranchAddress("sf_unc3", &sf_unc3, &b_sf_unc3);
   fChain->SetBranchAddress("sf_unc4", &sf_unc4, &b_sf_unc4);
   fChain->SetBranchAddress("normalizedWeight", &normalizedWeight, &b_normalizedWeight);
   Notify();
}

Bool_t pandaFlat::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pandaFlat::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pandaFlat::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pandaFlat_cxx
