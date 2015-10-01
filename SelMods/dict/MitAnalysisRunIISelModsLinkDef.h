#ifndef MITPHYSICS_SELMODS_LINKDEF_H
#define MITPHYSICS_SELMODS_LINKDEF_H

#include "MitAnalysisRunII/SelMods/interface/LeptonExampleMod.h"
#endif
 
#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::LeptonExampleMod+;
#endif
