# ----------------------------------------------------------------------
# -- USAGE:
# ---------
# gmake lib    # to build the so-lib to make a library with common helper classes
#
# NOTE: 1. gmake, not make
# ----------------------------------------------------------------------

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
 
ifeq ($(BFARCH),Linux24SL3_i386_gcc323)
# -- Linux2
CXX           = g++
CXXFLAGS      = -ggdb  -Wall -fPIC
#CXXFLAGS      = -O2 -Wall -fPIC
LD            = g++
LDFLAGS       = -ggdb
#LDFLAGS       = -O2
SOFLAGS       = -shared
endif
ifeq ($(BFARCH),Linux2)
# -- Linux2
CXX           = g++-2.95.2
CXXFLAGS      = -g  -Wall -fPIC
#CXXFLAGS      = -O2 -Wall -fPIC
LD            = g++-2.95.2
LDFLAGS       = -g
#LDFLAGS       = -O2
SOFLAGS       = -shared
endif
ifeq ($(BFARCH),Linux24)
# -- Linux24
CXX           = g++-2.95.3
CXXFLAGS      = -g  -Wall -fPIC
#CXXFLAGS      = -O2 -Wall -fPIC
LD            = g++-2.95.3
LDFLAGS       = -g
#LDFLAGS       = -O2
SOFLAGS       = -shared
endif
ifeq ($(BFARCH),SunOS58)
# -- Solaris
CXX           = CC
CXXFLAGS      = -g -KPIC
#CXXFLAGS      = -O2 -KPIC
LD            = CC
#LDFLAGS       = -g
LDFLAGS       = -O2
SOFLAGS       = -G
endif

CXXFLAGS      += $(ROOTCFLAGS)
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
NGLIBS        += ../shlib/$(BFARCH)/libRecoilAnalysis.so 
GLIBS          = $(filter-out -lNew, $(NGLIBS))

VCBNTUPLES = mesData.o  mesDataDict.o \
          recoilAnalysis.o recoilAnalysisDict.o  \
          PIDData.o PIDTable.o  PIDTableDict.o  PIDDataDict.o \
          TRKData.o TRKTable.o TRKTableDict.o  TRKDataDict.o \
          p_energy_loss_corrected.o  recoilDSys.o recoilDSysDict.o recoilBuSys.o

NTUPLES = $(VCBNTUPLES) mesFIT.o   mesFITDict.o  \
	  recoilBase.o  baseClass.o exampleAnalysis.o


ICCFILES = bookHist.icc  init.icc  mcTruth.icc splitOff.icc

$(addprefix ../tmp/$(BFARCH)/RecoilAnalysis/, %.o) : %.cc 
	$(CXX) $(CXXFLAGS) -c $< -o $@


# ================================================================================
lib:  $(addprefix ../tmp/$(BFARCH)/RecoilAnalysis/, $(NTUPLES))
	$(CXX) $(LDFLAGS) $(SOFLAGS) $(addprefix ../tmp/$(BFARCH)/RecoilAnalysis/, $(NTUPLES)) -o ../shlib/$(BFARCH)/libRecoilAnalysis.so
#---------------------------------------------------------------
libVcb:  $(addprefix ../tmp/$(BFARCH)/RecoilAnalysis/, $(VCBNTUPLES))
	$(CXX) $(SOFLAGS) $(addprefix ../tmp/$(BFARCH)/RecoilAnalysis/, $(VCBNTUPLES)) -o ../shlib/$(BFARCH)/libVcbRecoilAnalysis.so
#---------------------------------------------------------------

recoilAnalysisDict.cc:  mesData.hh recoilAnalysis.hh ntuplesLinkDef.h
	$(ROOTSYS)/bin/rootcint -f recoilAnalysisDict.cc -c  mesData.hh recoilAnalysis.hh ntuplesLinkDef.h

recoilBaseDict.cc:  mesData.hh recoilBase.hh ntuplesLinkDef.h
	$(ROOTSYS)/bin/rootcint -f recoilBaseDict.cc -c  mesData.hh recoilBase.hh ntuplesLinkDef.h

baseClassDict.cc:  mesData.hh baseClass.hh ntuplesLinkDef.h
	$(ROOTSYS)/bin/rootcint -f baseClassDict.cc -c  mesData.hh baseClass.hh ntuplesLinkDef.h

recoilDSysDict.cc:  recoilDSys.hh  recoilDSysLinkDef.h
	$(ROOTSYS)/bin/rootcint -f recoilDSysDict.cc -c -I../  recoilDSys.hh recoilDSysLinkDef.h

recoilBuSysDict.cc:  recoilBuSys.hh  ntuplesLinkDef.h
	$(ROOTSYS)/bin/rootcint -f recoilBuSysDict.cc -c -I../  recoilBuSys.hh  ntuplesLinkDef.h

mesDataDict.cc:  mesData.hh 
	$(ROOTSYS)/bin/rootcint -f mesDataDict.cc -c mesData.hh 

mesFITDict.cc:  mesFIT.hh 
	$(ROOTSYS)/bin/rootcint -f mesFITDict.cc -c mesFIT.hh 

PIDTableDict.cc: PIDTable.hh PIDTableLinkDef.h
	$(ROOTSYS)/bin/rootcint -f PIDTableDict.cc -c PIDTable.hh PIDTableLinkDef.h

TRKTableDict.cc: TRKTable.hh TRKTableLinkDef.h
	$(ROOTSYS)/bin/rootcint -f TRKTableDict.cc -c TRKTable.hh TRKTableLinkDef.h

PIDDataDict.cc: PIDData.hh
	$(ROOTSYS)/bin/rootcint -f $@ -c $^

TRKDataDict.cc: TRKData.hh
	$(ROOTSYS)/bin/rootcint -f $@ -c $^

baseClass.o: baseClass.cc baseClass.hh
	$(CXX) $(CXXFLAGS) -DFAST -I../ -c baseClass.cc -o ../tmp/$(BFARCH)/RecoilAnalysis/baseClass.o

# ================================================================================

exampleAnalysisApp: ../tmp/$(BFARCH)/RecoilAnalysis/exampleAnalysisApp.o lib
# ---------------------------------------------------------
	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/exampleAnalysisApp ../tmp/$(BFARCH)/RecoilAnalysis/exampleAnalysisApp.o $(GLIBS)

# ================================================================================


clean:
	rm -f  *Dict*
	rm -rf SunWS_cache
	rm -f  ../tmp/$(BFARCH)/RecoilAnalysis/*.o 
	rm -f  ../shlib/$(BFARCH)/libRecoilAnalysis.so 
	rm -f  ../shlib/$(BFARCH)/libVcbRecoilAnalysis.so 
	rm -f  ../bin/$(BFARCH)/exampleAnalysisApp

