 // $Id $

#include "MitAnalysisRunII/SelMods/interface/LeptonExampleMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "TFile.h"
#include "TTree.h"

using namespace mithep;
ClassImp(mithep::LeptonExampleMod)

//--------------------------------------------------------------------------------------------------
LeptonExampleMod::LeptonExampleMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fVertexName("NoDefaultNameSet"),
  fMuonName("NoDefaultNameSet"),
  fElectronName("NoDefaultNameSet"),
  fNEventsRun(0),
  fNEventsSelected(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void LeptonExampleMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void LeptonExampleMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

}

//--------------------------------------------------------------------------------------------------
void LeptonExampleMod::Process()
{
  fNEventsRun++;

  auto* OriginalVertex = GetObject<mithep::Collection<Vertex> >(fVertexName, true);
  auto* OriginalCleanMuons = GetObject<mithep::Collection<Muon> >(fMuonName, true);
  auto* OriginalCleanElectrons  = GetObject<mithep::Collection<Electron> >(fElectronName, true);

  mithep::NFArrBool const* ids[nLeptons][32]{};
  for (unsigned iL(0); iL != nLeptons; ++iL) {
    for (unsigned iSel(0); iSel != 32; ++iSel) {
      if (idName_[iL][iSel].Length() != 0)
        ids[iL][iSel] = GetObject<mithep::NFArrBool>(idName_[iL][iSel]);
    }
  }

  Int_t eventLeptonCount = 0;

  unsigned iE(0);
  unsigned iM(0);

  while (true) {
    mithep::Electron* ele = 0;
    if (iE != OriginalCleanElectrons->GetEntries())
      ele = OriginalCleanElectrons->At(iE);

    mithep::Muon* mu = 0;
    if (iM != OriginalCleanMuons->GetEntries())
      mu = OriginalCleanMuons->At(iM);

    if (!ele && !mu)
      break;

    if (ele && mu) {
      if (ele->Pt() > mu->Pt())
        mu = 0;
      else
        ele = 0;
    }

    if (ele) {
      // at least one lepton Id should be true
      unsigned iSel(0);
      for (; iSel != 32; ++iSel) {
        if (ids[kEl][iSel] && ids[kEl][iSel]->At(iE))
          break;
      }

      if (iSel != 32 && ele->Pt() > 10) {
        eventLeptonCount++;
      }
      ++iE;
    }
    else {
      // at least one lepton Id should be true
      unsigned iSel(0);
      for (; iSel != 32; ++iSel) {
        if (ids[kMu][iSel] && ids[kMu][iSel]->At(iM))
          break;
      }

      if (iSel != 32 && mu->Pt() > 10) {
        eventLeptonCount++;
      }
      ++iM;
    }
  }

  bool passAllCuts = OriginalVertex->GetEntries() > 0 &&
                     eventLeptonCount >= 2;

  if(passAllCuts) fNEventsSelected++;
  else            SkipEvent();

  return;
}

//--------------------------------------------------------------------------------------------------
void LeptonExampleMod::SlaveTerminate()
{
  
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here
  if(fNEventsRun > 0)
    std::cout << "run/selected events on LeptonExampleMod(" << fNEventsSelected << "/" << fNEventsRun << ") = " << (double)fNEventsSelected/fNEventsRun << std::endl;
  else
    std::cout << "0 run events on LeptonExampleMod" << std::endl;

}
//--------------------------------------------------------------------------------------------------
void LeptonExampleMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.

}
