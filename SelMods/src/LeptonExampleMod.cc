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

  auto* OriginalVertex = GetObject<mithep::Collection<Vertex> >(fVertexName, true);
  auto* OriginalCleanMuons = GetObject<mithep::Collection<Muon> >(fMuonName, true);
  auto* OriginalCleanElectrons  = GetObject<mithep::Collection<Electron> >(fElectronName, true);

  bool passAllCuts = OriginalVertex->GetEntries() > 0 &&
                     ((OriginalCleanMuons    ->GetEntries() > 0 && OriginalCleanMuons    ->At(0)->Pt() > 10) ||
                      (OriginalCleanElectrons->GetEntries() > 0 && OriginalCleanElectrons->At(0)->Pt() > 10));

  if(passAllCuts) fNEventsSelected++;
  else            SkipEvent();

  return;
}

//--------------------------------------------------------------------------------------------------
void LeptonExampleMod::SlaveTerminate()
{
  
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.
  std::cout << "selected events on LeptonExampleMod: " << fNEventsSelected << std::endl;

} 
//--------------------------------------------------------------------------------------------------
void LeptonExampleMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.

}
