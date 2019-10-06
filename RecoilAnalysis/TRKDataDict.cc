//
// File generated by /afs/slac.stanford.edu/g/babar/package/root/5.14-00e/Linux24SL3_i386_gcc323/bin/rootcint at Tue Mar  2 05:31:34 2010

// Do NOT change. Changes will be lost next time file is generated
//

#include "RConfig.h"
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "TRKDataDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TStreamerInfo.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TCollectionProxy.h"
#include "TIsAProxy.h"
// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void TRKData_ShowMembers(void *obj, TMemberInspector &R__insp, char *R__parent);
   static void *new_TRKData(void *p = 0);
   static void *newArray_TRKData(Long_t size, void *p);
   static void delete_TRKData(void *p);
   static void deleteArray_TRKData(void *p);
   static void destruct_TRKData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRKData*)
   {
      ::TRKData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRKData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRKData", ::TRKData::Class_Version(), "TRKData.hh", 8,
                  typeid(::TRKData), DefineBehavior(ptr, ptr),
                  &::TRKData::Dictionary, isa_proxy, 0,
                  sizeof(::TRKData) );
      instance.SetNew(&new_TRKData);
      instance.SetNewArray(&newArray_TRKData);
      instance.SetDelete(&delete_TRKData);
      instance.SetDeleteArray(&deleteArray_TRKData);
      instance.SetDestructor(&destruct_TRKData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRKData*)
   {
      return GenerateInitInstanceLocal((::TRKData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TRKData*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *TRKData::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *TRKData::Class_Name()
{
   return "TRKData";
}

//______________________________________________________________________________
const char *TRKData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRKData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRKData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRKData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void TRKData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRKData*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *TRKData::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRKData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void TRKData::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRKData.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> _pmin;
      R__b >> _pmax;
      R__b >> _mmin;
      R__b >> _mmax;
      R__b >> _tmin;
      R__b >> _tmax;
      R__b >> _fmin;
      R__b >> _fmax;
      R__b >> _e;
      R__b >> _s;
      R__b.CheckByteCount(R__s, R__c, TRKData::IsA());
   } else {
      R__c = R__b.WriteVersion(TRKData::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << _pmin;
      R__b << _pmax;
      R__b << _mmin;
      R__b << _mmax;
      R__b << _tmin;
      R__b << _tmax;
      R__b << _fmin;
      R__b << _fmax;
      R__b << _e;
      R__b << _s;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

//______________________________________________________________________________
void TRKData::ShowMembers(TMemberInspector &R__insp, char *R__parent)
{
      // Inspect the data members of an object of class TRKData.
      TClass *R__cl = ::TRKData::IsA();
      Int_t R__ncp = strlen(R__parent);
      if (R__ncp || R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__parent, "_pmin", &_pmin);
      R__insp.Inspect(R__cl, R__parent, "_pmax", &_pmax);
      R__insp.Inspect(R__cl, R__parent, "_mmin", &_mmin);
      R__insp.Inspect(R__cl, R__parent, "_mmax", &_mmax);
      R__insp.Inspect(R__cl, R__parent, "_tmin", &_tmin);
      R__insp.Inspect(R__cl, R__parent, "_tmax", &_tmax);
      R__insp.Inspect(R__cl, R__parent, "_fmin", &_fmin);
      R__insp.Inspect(R__cl, R__parent, "_fmax", &_fmax);
      R__insp.Inspect(R__cl, R__parent, "_e", &_e);
      R__insp.Inspect(R__cl, R__parent, "_s", &_s);
      TObject::ShowMembers(R__insp, R__parent);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRKData(void *p) {
      return  p ? new(p) ::TRKData : new ::TRKData;
   }
   static void *newArray_TRKData(Long_t nElements, void *p) {
      return p ? new(p) ::TRKData[nElements] : new ::TRKData[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRKData(void *p) {
      delete ((::TRKData*)p);
   }
   static void deleteArray_TRKData(void *p) {
      delete [] ((::TRKData*)p);
   }
   static void destruct_TRKData(void *p) {
      typedef ::TRKData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::TRKData

/********************************************************
* TRKDataDict.cc
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

extern "C" void G__cpp_reset_tagtableTRKDataDict();

extern "C" void G__set_cpp_environmentTRKDataDict() {
  G__add_compiledheader("TROOT.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("TRKData.hh");
  G__cpp_reset_tagtableTRKDataDict();
}
#include <new>
extern "C" int G__cpp_dllrevTRKDataDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* TRKData */
static int G__TRKDataDict_131_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TRKData* p = NULL;
   long gvp = G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == G__PVOID) || (gvp == 0)) {
       p = new TRKData[n];
     } else {
       p = new((void*) gvp) TRKData[n];
     }
   } else {
     if ((gvp == G__PVOID) || (gvp == 0)) {
       p = new TRKData;
     } else {
       p = new((void*) gvp) TRKData;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   result7->type = 'u';
   result7->tagnum = G__get_linked_tagnum(&G__TRKDataDictLN_TRKData);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TRKData* p = NULL;
   long gvp = G__getgvp();
   switch (libp->paran) {
   case 10:
     //m: 10
     if ((gvp == G__PVOID) || (gvp == 0)) {
       p = new TRKData(
(double) G__double(libp->para[0]), (double) G__double(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])
, (double) G__double(libp->para[4]), (double) G__double(libp->para[5])
, (double) G__double(libp->para[6]), (double) G__double(libp->para[7])
, (double) G__double(libp->para[8]), (double) G__double(libp->para[9]));
     } else {
       p = new((void*) gvp) TRKData(
(double) G__double(libp->para[0]), (double) G__double(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])
, (double) G__double(libp->para[4]), (double) G__double(libp->para[5])
, (double) G__double(libp->para[6]), (double) G__double(libp->para[7])
, (double) G__double(libp->para[8]), (double) G__double(libp->para[9]));
     }
     break;
   case 9:
     //m: 9
     if ((gvp == G__PVOID) || (gvp == 0)) {
       p = new TRKData(
(double) G__double(libp->para[0]), (double) G__double(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])
, (double) G__double(libp->para[4]), (double) G__double(libp->para[5])
, (double) G__double(libp->para[6]), (double) G__double(libp->para[7])
, (double) G__double(libp->para[8]));
     } else {
       p = new((void*) gvp) TRKData(
(double) G__double(libp->para[0]), (double) G__double(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])
, (double) G__double(libp->para[4]), (double) G__double(libp->para[5])
, (double) G__double(libp->para[6]), (double) G__double(libp->para[7])
, (double) G__double(libp->para[8]));
     }
     break;
   case 8:
     //m: 8
     if ((gvp == G__PVOID) || (gvp == 0)) {
       p = new TRKData(
(double) G__double(libp->para[0]), (double) G__double(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])
, (double) G__double(libp->para[4]), (double) G__double(libp->para[5])
, (double) G__double(libp->para[6]), (double) G__double(libp->para[7]));
     } else {
       p = new((void*) gvp) TRKData(
(double) G__double(libp->para[0]), (double) G__double(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])
, (double) G__double(libp->para[4]), (double) G__double(libp->para[5])
, (double) G__double(libp->para[6]), (double) G__double(libp->para[7]));
     }
     break;
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   result7->type = 'u';
   result7->tagnum = G__get_linked_tagnum(&G__TRKDataDictLN_TRKData);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   TRKData* p = NULL;
   long gvp = G__getgvp();
   //m: 1
   if ((gvp == G__PVOID) || (gvp == 0)) {
     p = new TRKData(*(TRKData*) libp->para[0].ref);
   } else {
     p = new((void*) gvp) TRKData(*(TRKData*) libp->para[0].ref);
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   result7->type = 'u';
   result7->tagnum = G__get_linked_tagnum(&G__TRKDataDictLN_TRKData);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      {
         TRKData* pobj;
         TRKData xobj = ((TRKData*) G__getstructoffset())->operator=(*(TRKData*) libp->para[0].ref);
         pobj = new TRKData(xobj);
         result7->obj.i = (long) ((void*) pobj);
         result7->ref = result7->obj.i;
         G__store_tempobject(*result7);
      }
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 103, (long) ((const TRKData*) G__getstructoffset())->operator>=(*(TRKData*) libp->para[0].ref));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 103, (long) ((const TRKData*) G__getstructoffset())->operator<(*(TRKData*) libp->para[0].ref));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 103, (long) ((const TRKData*) G__getstructoffset())->operator==(*(TRKData*) libp->para[0].ref));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 103, (long) ((TRKData*) G__getstructoffset())->isCell((double) G__double(libp->para[0]), (double) G__double(libp->para[1])
, (double) G__double(libp->para[2]), (double) G__double(libp->para[3])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->print();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getPctr());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getMctr());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getTctr());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getFctr());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getPmax());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getMmax());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getTmax());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getFmax());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getPmin());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getMmin());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getTmin());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getFmin());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getE());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 100, (double) ((const TRKData*) G__getstructoffset())->getS());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_24(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->setPmax((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_25(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->setMmax((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_26(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->setTmax((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_27(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->setFmax((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_28(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->setPmin((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_29(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->setMmin((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_30(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->setTmin((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_31(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->setFmin((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_32(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->setE((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_33(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->setS((double) G__double(libp->para[0]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_34(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) TRKData::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_35(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TRKData::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_36(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) TRKData::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_37(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      TRKData::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_38(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((const TRKData*) G__getstructoffset())->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_39(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->ShowMembers(*(TMemberInspector*) libp->para[0].ref, (char*) G__int(libp->para[1]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_40(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->Streamer(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_41(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((TRKData*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_42(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TRKData::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_43(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TRKData::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_44(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) TRKData::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__TRKDataDict_131_0_45(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) TRKData::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef TRKData G__TTRKData;
static int G__TRKDataDict_131_0_46(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   long gvp = G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == G__PVOID) {
       delete[] (TRKData*) soff;
     } else {
       G__setgvp(G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((TRKData*) (soff+(sizeof(TRKData)*i)))->~G__TTRKData();
       }
       G__setgvp(gvp);
     }
   } else {
     if (gvp == G__PVOID) {
       delete (TRKData*) soff;
     } else {
       G__setgvp(G__PVOID);
       ((TRKData*) (soff))->~G__TTRKData();
       G__setgvp(gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* TRKData */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncTRKDataDict {
 public:
  G__Sizep2memfuncTRKDataDict() {p=&G__Sizep2memfuncTRKDataDict::sizep2memfunc;}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncTRKDataDict::*p)();
};

size_t G__get_sizep2memfuncTRKDataDict()
{
  G__Sizep2memfuncTRKDataDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceTRKDataDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__TRKDataDictLN_TRKData))) {
     TRKData *G__Lderived;
     G__Lderived=(TRKData*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__TRKDataDictLN_TRKData),G__get_linked_tagnum(&G__TRKDataDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableTRKDataDict() {

   /* Setting up typedef entry */
   G__search_typename2("Bool_t",103,-1,0,-1);
   G__setnewtype(-1,"Boolean (0=false, 1=true) (bool)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<TStreamerInfo*>",117,G__get_linked_tagnum(&G__TRKDataDictLN_vectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__TRKDataDictLN_reverse_iteratorlEvectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TRKDataDictLN_vectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__TRKDataDictLN_reverse_iteratorlEvectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__TRKDataDictLN_vectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* TRKData */
static void G__setup_memvarTRKData(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__TRKDataDictLN_TRKData));
   { TRKData *p; p=(TRKData*)0x1000; if (p) { }
   G__memvar_setup((void*)NULL,100,0,0,-1,-1,-1,4,"_pmin=",0,(char*)NULL);
   G__memvar_setup((void*)NULL,100,0,0,-1,-1,-1,4,"_pmax=",0,(char*)NULL);
   G__memvar_setup((void*)NULL,100,0,0,-1,-1,-1,4,"_mmin=",0,(char*)NULL);
   G__memvar_setup((void*)NULL,100,0,0,-1,-1,-1,4,"_mmax=",0,(char*)NULL);
   G__memvar_setup((void*)NULL,100,0,0,-1,-1,-1,4,"_tmin=",0,(char*)NULL);
   G__memvar_setup((void*)NULL,100,0,0,-1,-1,-1,4,"_tmax=",0,(char*)NULL);
   G__memvar_setup((void*)NULL,100,0,0,-1,-1,-1,4,"_fmin=",0,(char*)NULL);
   G__memvar_setup((void*)NULL,100,0,0,-1,-1,-1,4,"_fmax=",0,(char*)NULL);
   G__memvar_setup((void*)NULL,100,0,0,-1,-1,-1,4,"_e=",0,(char*)NULL);
   G__memvar_setup((void*)NULL,100,0,0,-1,-1,-1,4,"_s=",0,(char*)NULL);
   G__memvar_setup((void*)NULL,85,0,0,G__get_linked_tagnum(&G__TRKDataDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarTRKDataDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncTRKData(void) {
   /* TRKData */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__TRKDataDictLN_TRKData));
   G__memfunc_setup("TRKData",619,G__TRKDataDict_131_0_1, 105, G__get_linked_tagnum(&G__TRKDataDictLN_TRKData), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TRKData",619,G__TRKDataDict_131_0_2, 105, G__get_linked_tagnum(&G__TRKDataDictLN_TRKData), -1, 0, 10, 1, 1, 0, 
"d - - 0 - p d - - 0 - pmax "
"d - - 0 - mult d - - 0 - mmax "
"d - - 0 - th d - - 0 - tmax "
"d - - 0 - phi d - - 0 - fmax "
"d - - 0 0. e d - - 0 0. s", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("TRKData",619,G__TRKDataDict_131_0_3, 105, G__get_linked_tagnum(&G__TRKDataDictLN_TRKData), -1, 0, 1, 1, 1, 0, "u 'TRKData' - 11 - orig", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("operator=",937,G__TRKDataDict_131_0_4, 117, G__get_linked_tagnum(&G__TRKDataDictLN_TRKData), -1, 0, 1, 1, 1, 0, "u 'TRKData' - 11 - cmp", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("operator>=",999,G__TRKDataDict_131_0_5, 103, -1, G__defined_typename("Bool_t"), 0, 1, 1, 1, 8, "u 'TRKData' - 1 - cmp", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("operator<",936,G__TRKDataDict_131_0_6, 103, -1, G__defined_typename("Bool_t"), 0, 1, 1, 1, 8, "u 'TRKData' - 11 - cmp", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("operator==",998,G__TRKDataDict_131_0_7, 103, -1, G__defined_typename("Bool_t"), 0, 1, 1, 1, 8, "u 'TRKData' - 11 - cmp", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("isCell",604,G__TRKDataDict_131_0_8, 103, -1, G__defined_typename("Bool_t"), 0, 4, 1, 1, 0, 
"d - - 0 - p d - - 0 - m "
"d - - 0 - t d - - 0 - f", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("print",557,G__TRKDataDict_131_0_9, 121, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getPctr",729,G__TRKDataDict_131_0_10, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getMctr",726,G__TRKDataDict_131_0_11, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getTctr",733,G__TRKDataDict_131_0_12, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getFctr",719,G__TRKDataDict_131_0_13, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getPmax",726,G__TRKDataDict_131_0_14, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getMmax",723,G__TRKDataDict_131_0_15, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getTmax",730,G__TRKDataDict_131_0_16, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getFmax",716,G__TRKDataDict_131_0_17, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getPmin",724,G__TRKDataDict_131_0_18, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getMmin",721,G__TRKDataDict_131_0_19, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getTmin",728,G__TRKDataDict_131_0_20, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getFmin",714,G__TRKDataDict_131_0_21, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getE",389,G__TRKDataDict_131_0_22, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getS",403,G__TRKDataDict_131_0_23, 100, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setPmax",738,G__TRKDataDict_131_0_24, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - p", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setMmax",735,G__TRKDataDict_131_0_25, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - m", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setTmax",742,G__TRKDataDict_131_0_26, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - th", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setFmax",728,G__TRKDataDict_131_0_27, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - phi", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setPmin",736,G__TRKDataDict_131_0_28, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - p", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setMmin",733,G__TRKDataDict_131_0_29, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - m", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setTmin",740,G__TRKDataDict_131_0_30, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - th", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setFmin",726,G__TRKDataDict_131_0_31, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - phi", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setE",401,G__TRKDataDict_131_0_32, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - e", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("setS",415,G__TRKDataDict_131_0_33, 121, -1, -1, 0, 1, 1, 1, 0, "d - - 0 - s", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__TRKDataDict_131_0_34, 85, G__get_linked_tagnum(&G__TRKDataDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (TClass* (*)())(&TRKData::Class), 0);
   G__memfunc_setup("Class_Name",982,G__TRKDataDict_131_0_35, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) (const char* (*)())(&TRKData::Class_Name), 0);
   G__memfunc_setup("Class_Version",1339,G__TRKDataDict_131_0_36, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (Version_t (*)())(&TRKData::Class_Version), 0);
   G__memfunc_setup("Dictionary",1046,G__TRKDataDict_131_0_37, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (void (*)())(&TRKData::Dictionary), 0);
   G__memfunc_setup("IsA",253,G__TRKDataDict_131_0_38, 85, G__get_linked_tagnum(&G__TRKDataDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,G__TRKDataDict_131_0_39, 121, -1, -1, 0, 2, 1, 1, 0, 
"u 'TMemberInspector' - 1 - insp C - - 0 - parent", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,G__TRKDataDict_131_0_40, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__TRKDataDict_131_0_41, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__TRKDataDict_131_0_42, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) (const char* (*)())(&TRKData::DeclFileName), 0);
   G__memfunc_setup("ImplFileLine",1178,G__TRKDataDict_131_0_43, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (int (*)())(&TRKData::ImplFileLine), 0);
   G__memfunc_setup("ImplFileName",1171,G__TRKDataDict_131_0_44, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) (const char* (*)())(&TRKData::ImplFileName), 0);
   G__memfunc_setup("DeclFileLine",1152,G__TRKDataDict_131_0_45, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) (int (*)())(&TRKData::DeclFileLine), 0);
   // automatic destructor
   G__memfunc_setup("~TRKData", 745, G__TRKDataDict_131_0_46, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncTRKDataDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalTRKDataDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcTRKDataDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__TRKDataDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__TRKDataDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__TRKDataDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__TRKDataDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__TRKDataDictLN_vectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgR = { "vector<TStreamerInfo*,allocator<TStreamerInfo*> >" , 99 , -1 };
G__linked_taginfo G__TRKDataDictLN_reverse_iteratorlEvectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TStreamerInfo*,allocator<TStreamerInfo*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__TRKDataDictLN_TRKData = { "TRKData" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableTRKDataDict() {
  G__TRKDataDictLN_TClass.tagnum = -1 ;
  G__TRKDataDictLN_TBuffer.tagnum = -1 ;
  G__TRKDataDictLN_TMemberInspector.tagnum = -1 ;
  G__TRKDataDictLN_TObject.tagnum = -1 ;
  G__TRKDataDictLN_vectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgR.tagnum = -1 ;
  G__TRKDataDictLN_reverse_iteratorlEvectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__TRKDataDictLN_TRKData.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableTRKDataDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum(&G__TRKDataDictLN_TClass);
   G__get_linked_tagnum(&G__TRKDataDictLN_TBuffer);
   G__get_linked_tagnum(&G__TRKDataDictLN_TMemberInspector);
   G__get_linked_tagnum(&G__TRKDataDictLN_TObject);
   G__get_linked_tagnum(&G__TRKDataDictLN_vectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgR);
   G__get_linked_tagnum(&G__TRKDataDictLN_reverse_iteratorlEvectorlETStreamerInfomUcOallocatorlETStreamerInfomUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum(&G__TRKDataDictLN_TRKData),sizeof(TRKData),-1,65280,"Testing TRKData",G__setup_memvarTRKData,G__setup_memfuncTRKData);
}
extern "C" void G__cpp_setupTRKDataDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupTRKDataDict()");
  G__set_cpp_environmentTRKDataDict();
  G__cpp_setup_tagtableTRKDataDict();

  G__cpp_setup_inheritanceTRKDataDict();

  G__cpp_setup_typetableTRKDataDict();

  G__cpp_setup_memvarTRKDataDict();

  G__cpp_setup_memfuncTRKDataDict();
  G__cpp_setup_globalTRKDataDict();
  G__cpp_setup_funcTRKDataDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncTRKDataDict();
  return;
}
class G__cpp_setup_initTRKDataDict {
  public:
    G__cpp_setup_initTRKDataDict() { G__add_setup_func("TRKDataDict",(G__incsetup)(&G__cpp_setupTRKDataDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initTRKDataDict() { G__remove_setup_func("TRKDataDict"); }
};
G__cpp_setup_initTRKDataDict G__cpp_setup_initializerTRKDataDict;
