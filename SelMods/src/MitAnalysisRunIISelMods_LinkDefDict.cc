// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIceballosdIcmsdIcmsswdI044dICMSSW_8_0_7_patch3dItmpdIslc6_amd64_gcc530dIMitAnalysisRunIISelMods_LinkDefDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "MitAnalysisRunII/SelMods/interface/LeptonExampleMod.h"

// Header files passed via #pragma extra_include

namespace mithep {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *mithep_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("mithep", 0 /*version*/, "MitAna/DataCont/interface/BaseCollection.h", 14,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &mithep_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *mithep_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace ROOT {
   static void *new_mithepcLcLLeptonExampleMod(void *p = 0);
   static void *newArray_mithepcLcLLeptonExampleMod(Long_t size, void *p);
   static void delete_mithepcLcLLeptonExampleMod(void *p);
   static void deleteArray_mithepcLcLLeptonExampleMod(void *p);
   static void destruct_mithepcLcLLeptonExampleMod(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::mithep::LeptonExampleMod*)
   {
      ::mithep::LeptonExampleMod *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::mithep::LeptonExampleMod >(0);
      static ::ROOT::TGenericClassInfo 
         instance("mithep::LeptonExampleMod", ::mithep::LeptonExampleMod::Class_Version(), "MitAnalysisRunII/SelMods/interface/LeptonExampleMod.h", 25,
                  typeid(::mithep::LeptonExampleMod), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::mithep::LeptonExampleMod::Dictionary, isa_proxy, 4,
                  sizeof(::mithep::LeptonExampleMod) );
      instance.SetNew(&new_mithepcLcLLeptonExampleMod);
      instance.SetNewArray(&newArray_mithepcLcLLeptonExampleMod);
      instance.SetDelete(&delete_mithepcLcLLeptonExampleMod);
      instance.SetDeleteArray(&deleteArray_mithepcLcLLeptonExampleMod);
      instance.SetDestructor(&destruct_mithepcLcLLeptonExampleMod);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::mithep::LeptonExampleMod*)
   {
      return GenerateInitInstanceLocal((::mithep::LeptonExampleMod*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::mithep::LeptonExampleMod*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace mithep {
//______________________________________________________________________________
atomic_TClass_ptr LeptonExampleMod::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LeptonExampleMod::Class_Name()
{
   return "mithep::LeptonExampleMod";
}

//______________________________________________________________________________
const char *LeptonExampleMod::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::mithep::LeptonExampleMod*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LeptonExampleMod::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::mithep::LeptonExampleMod*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LeptonExampleMod::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::mithep::LeptonExampleMod*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LeptonExampleMod::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::mithep::LeptonExampleMod*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace mithep
namespace mithep {
//______________________________________________________________________________
void LeptonExampleMod::Streamer(TBuffer &R__b)
{
   // Stream an object of class mithep::LeptonExampleMod.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(mithep::LeptonExampleMod::Class(),this);
   } else {
      R__b.WriteClassBuffer(mithep::LeptonExampleMod::Class(),this);
   }
}

} // namespace mithep
namespace ROOT {
   // Wrappers around operator new
   static void *new_mithepcLcLLeptonExampleMod(void *p) {
      return  p ? new(p) ::mithep::LeptonExampleMod : new ::mithep::LeptonExampleMod;
   }
   static void *newArray_mithepcLcLLeptonExampleMod(Long_t nElements, void *p) {
      return p ? new(p) ::mithep::LeptonExampleMod[nElements] : new ::mithep::LeptonExampleMod[nElements];
   }
   // Wrapper around operator delete
   static void delete_mithepcLcLLeptonExampleMod(void *p) {
      delete ((::mithep::LeptonExampleMod*)p);
   }
   static void deleteArray_mithepcLcLLeptonExampleMod(void *p) {
      delete [] ((::mithep::LeptonExampleMod*)p);
   }
   static void destruct_mithepcLcLLeptonExampleMod(void *p) {
      typedef ::mithep::LeptonExampleMod current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::mithep::LeptonExampleMod

namespace {
  void TriggerDictionaryInitialization_MitAnalysisRunIISelMods_LinkDefDict_Impl() {
    static const char* headers[] = {
"MitAnalysisRunII/SelMods/interface/LeptonExampleMod.h",
0
    };
    static const char* includePaths[] = {
"/home/ceballos/cms/cmssw/044/CMSSW_8_0_7_patch3/src",
"/cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.06.00-ikhhed3/include",
"/home/ceballos/cms/cmssw/044/CMSSW_8_0_7_patch3/src/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MitAnalysisRunIISelMods_LinkDefDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace mithep{class __attribute__((annotate(R"ATTRDUMP(TAM example analysis module)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$MitAnalysisRunII/SelMods/dict/MitAnalysisRunIISelModsLinkDef.h")))  LeptonExampleMod;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MitAnalysisRunIISelMods_LinkDefDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "MitAnalysisRunII/SelMods/interface/LeptonExampleMod.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"mithep::LeptonExampleMod", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MitAnalysisRunIISelMods_LinkDefDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MitAnalysisRunIISelMods_LinkDefDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MitAnalysisRunIISelMods_LinkDefDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MitAnalysisRunIISelMods_LinkDefDict() {
  TriggerDictionaryInitialization_MitAnalysisRunIISelMods_LinkDefDict_Impl();
}
