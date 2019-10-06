# Makefile for the thecomparison
# # using the ROOT libraries on Linux platform
#
# Author: Antonio Petrella, 26 jan 2006
#
# needs to be upgraded to other arch.

#------------------------------------------------------------------------------
ifeq ($(BFARCH),Linux24SL3_i386_gcc323) # -- Linux24

CXX           = g++
#CXXFLAGS      = -g -fPIC -Wall 
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared

endif

SrcSuf        = cxx
ObjSuf	      = o
ExeSuf        =
LibSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

 ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
 ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
 ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#ROOTCFLAGS   := $(shell root-config --cflags)
#ROOTLIBS     := $(shell root-config --libs)
#ROOTGLIBS    := $(shell root-config --glibs)
ROOTGLIBS     += $(ROOTSYS)/lib/libHtml.so 

COMPOBJ         = thecomparison.$(ObjSuf) thecomparisonDict.$(ObjSuf)
COMPSRC         = thecomparison.$(SrcSuf) thecomparisonDict.$(SrcSuf)	
COMPSO          = libthecomparison.$(LibSuf) 

RECOILDSYSOBJ  = recoilDSysDict.$(ObjSuf)
RECOILDSYSSRC  = recoilDSys.$(SrcSuf) recoilDSysDict.$(SrcSuf)
RECOILDSYSSO   = librecoilDSys.$(LibSuf)

RECOILOBJ  = recoilAnalysisDict.$(ObjSuf)
RECOILSRC  = recoilAnalysis.$(SrcSuf) recoilAnalysisDict.$(SrcSuf)
RECOILSO   = librecoilAnalysis.$(LibSuf)

MESDATAOBJ = mesDataDict.$(ObjSuf)
MESDATASRC = mesData.$(SrcSuf) mesDataDict.$(SrcSuf)
MESDATASO  = libmesData.$(LibSuf)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)
SYSLIBS =     -L../../shlib/$(BFARCH) -lrecoilAnalysis -lmesData -lrecoilDSys

.SUFFIXES: .cc .o .so
#Dependencies of libthecomparison.so
lib:  $(COMPOBJ) $(MESDATASO) $(RECOILSO) $(RECOILDSYSSO) 
#Rule to create libthecomparison.so	
	$(LD) $(SOFLAGS) $(LDFLAGS) $(COMPOBJ) $(OutPutOpt) ../../shlib/$(BFARCH)/libthecomparison.$(LibSuf) $(SYSLIBS)
	rm -f *.$(ObjSuf)

thecomparison.o: thecomparison.C util.cc thecomparison.h util.hh
	$(CXX) $(CXXFLAGS) -c thecomparison.C util.cc
thecomparisonDict.o: thecomparisonDict.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -c $<
thecomparisonDict.$(SrcSuf): thecomparison.h thecomparisonLinkDef.h
	@echo "Generating dictionary $@..."
	$(ROOTSYS)/bin/rootcint -f $@ -c $^

#Dependencies of librecoilDSys.so
$(RECOILDSYSSO): $(RECOILDSYSOBJ) 
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) ../../shlib/$(BFARCH)/$@ 

recoilDSys.o:  ../../RecoilAnalysis/recoilDSys.cc  ../../RecoilAnalysis/recoilDSys.hh  ../../RecoilAnalysis/mesData.hh 
recoilDSysDict.o: recoilDSysDict.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -c $<
recoilDSysDict.$(SrcSuf):  ../../RecoilAnalysis/recoilDSys.hh  ../../RecoilAnalysis/mesData.hh  ../../RecoilAnalysis/recoilDSysLinkDef.h
	@echo "Generating dictionary $@..."
	$(ROOTSYS)/bin/rootcint -f $@ -c  ../../RecoilAnalysis/recoilDSys.cc  ../../RecoilAnalysis/mesData.cc
#Dependencies of librecoilAnalysis.so
$(RECOILSO): $(RECOILOBJ) 
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) ../../shlib/$(BFARCH)/$@

recoilAnalysis.o: ../../RecoilAnalysis/recoilAnalysis.cc ../../RecoilAnalysis/recoilAnalysis.hh ../../RecoilAnalysis/mesData.hh
recoilAnalysisDict.o: recoilAnalysisDict.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -c $<
recoilAnalysisDict.$(SrcSuf): ../../RecoilAnalysis/recoilAnalysis.hh ../../RecoilAnalysis/mesData.hh ../../RecoilAnalysis/recoilAnalysisLinkDef.h
	@echo "Generating dictionary $@..."
	$(ROOTSYS)/bin/rootcint -f $@ -c ../../RecoilAnalysis/recoilAnalysis.cc ../../RecoilAnalysis/mesData.cc

#Dependencies of libmesData.so
$(MESDATASO): $(MESDATAOBJ) 
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) ../../shlib/$(BFARCH)/$@

mesData.o:   ../../RecoilAnalysis/mesData.cc ../../RecoilAnalysis/mesData.hh
mesDataDict.o: mesDataDict.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -c $<
mesDataDict.$(SrcSuf): ../../RecoilAnalysis/mesData.hh ../../RecoilAnalysis/mesDataLinkDef.h
	@echo "Generating dictionary $@..."
	$(ROOTSYS)/bin/rootcint -f $@ -c ../../RecoilAnalysis/mesData.cc


%.cc %.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.$(ObjSuf)
	rm -f *.cxx
	rm -f *Dict.h
cleanall:
	rm -f *.$(ObjSuf) 
	rm -f ../../shlib/$(BFARCH)/$(MESDATASO) 
	rm -f ../../shlib/$(BFARCH)/$(RECOILDSYSSO)
	rm -f ../../shlib/$(BFARCH)/$(RECOILSO)
	rm -f ../../shlib/$(BFARCH)/$(COMPSO)
	rm -f *Dict.h
	rm -f *.cxx