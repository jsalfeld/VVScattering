//--------------------------------------------------------------------------------------------------
// $Id $
//
// LeptonExampleMod
//
// A Module for Selecting leptons
// and produces some distributions.
//
//--------------------------------------------------------------------------------------------------

#ifndef ANALYSISRUNII_SELMODS_LEPTONEXAMPLEMOD_H
#define ANALYSISRUNII_SELMODS_LEPTONEXAMPLEMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class LeptonExampleMod : public BaseMod
  {
    public:
    LeptonExampleMod(const char *name="LeptonExampleMod", 
		 const char *title="Example analysis module with all branches");
      ~LeptonExampleMod() {}

      const char *GetVertexName()             const { return fVertexName;   }
      const char *GetMuonName()               const { return fMuonName;     }
      const char *GetElectronName()           const { return fElectronName; }

      void  SetVertexName(const char *name)         { fVertexName   = name; }
      void  SetMuonName(const char *name)           { fMuonName     = name; }
      void  SetElectronName(const char *name)       { fElectronName = name; }

    protected:
      TString fVertexName;	//name of vertex collection
      TString fMuonName;	//name of muon collection
      TString fElectronName;	//name of electron collection
      Int_t   fNEventsSelected; //selected events

      void         Begin();
      void         Process();
      void         SlaveBegin();
      void         SlaveTerminate();
      void         Terminate();      

      ClassDef(LeptonExampleMod,1) // TAM example analysis module
  };
}
#endif
