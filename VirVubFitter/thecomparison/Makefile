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

COMPOBJ         = thecomparison.$(ObjSuf) thecomparisonDict.$(ObjSuf) util.$(ObjSuf) RooCCB.$(ObjSuf) RooThorstenSig.$(ObjSuf) recoilAnalysis.$(ObjSuf) recoilDSys.$(ObjSuf) mesData.$(ObjSuf)
COMPSRC         = thecomparison.$(SrcSuf) thecomparisonDict.$(SrcSuf)
COMPSO          = libthecomparison.$(LibSuf) 

CXXFLAGS     += $(ROOTCFLAGS) $(ROOFITCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)
SYSLIBS =     -L../../shlib/$(BFARCH)

ROOFITCFLAGS = -I. -I../../RooFitCore/tmp -I../../RooFitModels/tmp -I../..
ROOFITLIBS   = -L../../RooFitCore/tmp -L../../RooFitModels/tmp  -lRooFitCore -lRooFitModels

.SUFFIXES: .cc .o .so
#Dependencies of libthecomparison.so
#lib:  $(COMPOBJ) $(MESDATASO) $(RECOILSO) $(RECOILDSYSSO) 
lib:  $(COMPOBJ)
#Rule to create libthecomparison.so	
	$(LD) $(SOFLAGS) $(LDFLAGS) $(COMPOBJ) $(OutPutOpt) ../../shlib/$(BFARCH)/libthecomparison.$(LibSuf) $(SYSLIBS) $(ROOFITLIBS) 
#	rm -f *.$(ObjSuf)

thecomparison.o: thecomparison.C thecomparison.h util.hh RooCCB.hh RooThorstenSig.hh
thecomparisonDict.o: thecomparisonDict.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

thecomparisonDict.$(SrcSuf): thecomparison.h RooCCB.hh RooThorstenSig.hh ../../RecoilAnalysis/recoilAnalysis.hh thecomparisonLinkDef.h
	@echo "Generating dictionary $@..."
	$(ROOTSYS)/bin/rootcint -f $@ -c $(ROOFITCFLAGS) $^

recoilDSys.o:  ../../RecoilAnalysis/recoilDSys.cc  ../../RecoilAnalysis/recoilDSys.hh  ../../RecoilAnalysis/mesData.hh 
	$(CXX) $(CXXFLAGS) -o $@ -c $<
recoilAnalysis.o: ../../RecoilAnalysis/recoilAnalysis.cc ../../RecoilAnalysis/recoilAnalysis.hh ../../RecoilAnalysis/mesData.hh
	$(CXX) $(CXXFLAGS) -o $@ -c $<
mesData.o:   ../../RecoilAnalysis/mesData.cc ../../RecoilAnalysis/mesData.hh
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -f *.$(ObjSuf)
	rm -f *.cxx
	rm -f *Dict.h
cleanall:
	rm -f *.$(ObjSuf) 
	rm -f ../../shlib/$(BFARCH)/$(COMPSO)
	rm -f *Dict.h
	rm -f *.cxx
