# ----------------------------------------------------------------------
# -- USAGE:
# ---------
# gmake -f Makefile lib       # to build the so-lib to analyse the ntuples
# gmake -f Makefile anaRecoil # to build the executable for running on the ntuples
#
# NOTE: 1. gmake, not make
# ----------------------------------------------------------------------

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
ROOTGLIBS     += $(ROOTSYS)/lib/libHtml.so 
ifeq ($(BFARCH),Linux24SL3_i386_gcc323)
# -- Linux24
CXX           = g++
#CXXFLAGS      = -g  -fPIC -Wall 
CXXFLAGS      = -ggdb  -fPIC  
##CXXFLAGS      = -O2  -Wall -fPIC
LD            = g++
LDFLAGS       = -ggdb
#LDFLAGS       = -O2
SOFLAGS       = -ggdb -shared
endif
ifeq ($(BFARCH),Linux24)
# -- Linux24
CXX           = g++-2.95.3 
CXXFLAGS      = -g  -Wall -fPIC
CXXFLAGS      = -O2  -Wall -fPIC
LD            = g++-2.95.3 
LDFLAGS       = -g
LDFLAGS       = -O2
SOFLAGS       = -shared
endif
ifeq ($(BFARCH),SunOS58)
# -- Solaris
CXX           = CC
CXXFLAGS      = -g -KPIC 
CXXFLAGS      = -O2 -KPIC 
LD            = CC
#LD            = purify CC
LDFLAGS       = -g
LDFLAGS       = -O2
SOFLAGS       = -G
endif

CXXFLAGS       += $(ROOTCFLAGS)
NLIBS            = $(ROOTGLIBS) 

NGLIBSFIT       = $(ROOTGLIBS) 
NGLIBSFIT      += -lMinuit
NGLIBSFIT      += ../shlib/$(BFARCH)/libVirFit.so
NGLIBSFIT      += ../shlib/$(BFARCH)/libRecoilAnalysis.so
NGLIBSFIT      += ../RooFitCore/tmp/libRooFitCore.so
NGLIBSFIT      += ../RooFitModels/tmp/libRooFitModels.so

NGLIBSEXCLFIT       = $(ROOTGLIBS) 
NGLIBSEXCLFIT      += -lMinuit
NGLIBSEXCLFIT      += ../shlib/$(BFARCH)/libExclFit.so
#NGLIBSEXCLFIT      += ../shlib/$(BFARCH)/libVirFit.so
NGLIBSEXCLFIT      += ../shlib/$(BFARCH)/libRecoilAnalysis.so
NGLIBSEXCLFIT      += ../RooFitCore/tmp/libRooFitCore.so
NGLIBSEXCLFIT      += ../RooFitModels/tmp/libRooFitModels.so

NGLIBSMXFIT       = $(ROOTGLIBS) 
NGLIBSMXFIT      += -lMinuit
NGLIBSMXFIT      += ../shlib/$(BFARCH)/libmXFit.so
NGLIBSMXFIT      += ../shlib/$(BFARCH)/libRecoilAnalysis.so
NGLIBSMXFIT      += ../RooFitCore/tmp/libRooFitCore.so
NGLIBSMXFIT      += ../RooFitModels/tmp/libRooFitModels.so

NCLIBS         = $(NGLIBS)
NCLIBS        += $(CERNLIB)

LIBS           = $(filter-out -lNew, $(NLIBS))
GLIBS          = $(filter-out -lNew, $(NGLIBS))
CLIBS          = $(filter-out -lNew, $(NCLIBS))
GLIBSFIT       = $(filter-out -lNew, $(NGLIBSFIT))
GLIBSEXCLFIT   = $(filter-out -lNew, $(NGLIBSEXCLFIT))
GLIBSMXFIT     = $(filter-out -lNew, $(NGLIBSMXFIT))

# ======================================================================
# files which needs to be copied and passed through rootcint
RDLLIST = $(wildcard *.rdl)
HHLIST  = $(patsubst %.rdl,%.hh,$(RDLLIST))
CCLIST  = $(patsubst %.rdl,%.cc,$(RDLLIST))
OOLIST  = $(patsubst %.cc,%.o,$(CCLIST))

VIR = VirClass.o VirClassDict.o CMClass.o VirHelper.o PidCorrectMesMean.o XSLBToDstrlnu_DstrToDpi_LinearQ2.o XSLBToDstrlnu_DstrToDpi_CLN.o XSLKin.o XSLItgSimpsonIntegrator.o XSLItgPtrFunction.o XSLItgAbsIntegrator.o XSLItgAbsFunction.o XSLBall04_pilnu.o XSLBall04_etalnu.o XSLBall05.o XSLPseudoScalarFF.o XSLPseudoScalarISGW2.o XSLBToDlnu_CLN.o XSLVectorFF.o XSLVectorISGW2.o $(OOLIST)

EXCLFIT = exclfitNtp.o exclfitNtpDict.o 
MX = mXClass.o mXClassDict.o CMClass.o PidCorrectMesMean.o XSLBToDstrlnu_DstrToDpi_LinearQ2.o XSLBToDstrlnu_DstrToDpi_CLN.o XSLKin.o XSLItgSimpsonIntegrator.o XSLItgPtrFunction.o XSLItgAbsIntegrator.o XSLItgAbsFunction.o

ICCFILES = VirUtil.icc VirInit.icc exclfitInit.icc exclfitUtil.icc
MXICCFILES = mXUtil.icc mXInit.icc

# ======================================================================
$(addprefix ../tmp/$(BFARCH)/VirVubFitter/, %.o) : %.cc 
	$(CXX) $(CXXFLAGS) -I../ -I../RooFitCore/tmp -I../RooFitModels/tmp -I../tmp/$(BFARCH) -I../tmp/$(BFARCH)/VirVubFitter -c $< -o $@

# -- Uncomment this if you want recoilNtp to have a dependency on icc files. 
#    Will lead to automatic recompilation ...
VirClass.cc:  VirClass.hh CMClass.hh PidCorrectMesMean.hh $(ICCFILES) $(addprefix ../tmp/$(BFARCH)/VirVubFitter/, $(HHLIST))

#mXClass.cc:  mXClass.hh CMClass.hh PidCorrectMesMean.hh $(MXICCFILES) $(addprefix ../tmp/$(BFARCH)/VirVubFitter/, $(HHLIST))

# implicit rule to make header files for ROOT based classes
$(addprefix ../tmp/$(BFARCH)/VirVubFitter/, %.hh): %.rdl
	@cp $< $@

# ================================================================================
fit:  $(addprefix ../tmp/$(BFARCH)/VirVubFitter/, $(VIR))
#-------------------------------------------------------
	$(CXX) $(SOFLAGS) $(addprefix ../tmp/$(BFARCH)/VirVubFitter/, $(VIR)) -o ../shlib/$(BFARCH)/libVirFit.so

VirClassDict.cc:  VirClass.hh VirLinkDef.h  
	$(ROOTSYS)/bin/rootcint -f VirClassDict.cc -c -I../RooFitCore/tmp -I../RooFitModels/tmp -I../tmp/$(BFARCH)/VirVubFitter -I../ VirClass.hh CMClass.hh PidCorrectMesMean.hh $(HHLIST) VirLinkDef.h 


# ================================================================================
exclfit:  $(addprefix ../tmp/$(BFARCH)/VirVubFitter/, $(EXCLFIT))
#-------------------------------------------------------
	$(CXX) $(SOFLAGS) $(addprefix ../tmp/$(BFARCH)/VirVubFitter/, $(EXCLFIT)) -o ../shlib/$(BFARCH)/libExclFit.so 

exclfitNtpDict.cc:  exclfitNtp.hh  
	$(ROOTSYS)/bin/rootcint -f exclfitNtpDict.cc -c -I../  -I../RooFitCore/tmp -I../RooFitModels/tmp exclfitNtp.hh 

# ================================================================================
# ================================================================================
mxfit:  $(addprefix ../tmp/$(BFARCH)/VirVubFitter/, $(MX))
#-------------------------------------------------------
	$(CXX) $(SOFLAGS) $(addprefix ../tmp/$(BFARCH)/VirVubFitter/, $(MX)) -o ../shlib/$(BFARCH)/libmXFit.so

mXClassDict.cc:  mXClass.hh mXLinkDef.h  
	$(ROOTSYS)/bin/rootcint -f mXClassDict.cc -c -I../RooFitCore/tmp -I../RooFitModels/tmp -I../tmp/$(BFARCH)/VirVubFitter -I../ mXClass.hh CMClass.hh PidCorrectMesMean.hh mXLinkDef.h 


# ================================================================================
VirFit:  ../tmp/$(BFARCH)/VirVubFitter/VirFit.o fit
# --------------------------------
	$(LD) $(LDFLAGS) -o ../bin/$(BFARCH)/VirFit ../tmp/$(BFARCH)/VirVubFitter/VirFit.o $(GLIBSFIT)

# ================================================================================
FitOnly:  ../tmp/$(BFARCH)/VirVubFitter/FitOnly.o fit
# --------------------------------
	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/FitOnly ../tmp/$(BFARCH)/VirVubFitter/FitOnly.o $(GLIBSFIT)


# ================================================================================
ExclFit:  ../tmp/$(BFARCH)/VirVubFitter/ExclFit.o exclfit
# --------------------------------
	$(LD) $(LDFLAGS)   -o ../bin/$(BFARCH)/ExclFit ../tmp/$(BFARCH)/VirVubFitter/ExclFit.o $(GLIBSEXCLFIT) 

# ================================================================================
b2uFit:  
# --------------------------------
	$(CXX) $(CXXFLAGS) -I../ -I../RooFitCore/tmp -I../RooFitModels/tmp -c b2uFit.cc -o ../tmp/$(BFARCH)/VirVubFitter/b2uFit.o
	$(CXX) $(CXXFLAGS) -I../ -I../RooFitCore/tmp -I../RooFitModels/tmp -c b2uClass.cc -o ../tmp/$(BFARCH)/VirVubFitter/b2uClass.o
	$(LD) $(LDFLAGS)   -o ../bin/$(BFARCH)/b2uFit \
	../tmp/$(BFARCH)/VirVubFitter/b2uFit.o ../tmp/$(BFARCH)/VirVubFitter/b2uClass.o \
	$(ROOTGLIBS) -lMinuit \
	../RooFitCore/tmp/libRooFitCore.so ../RooFitModels/tmp/libRooFitModels.so \
	../shlib/$(BFARCH)/libVirFit.so \
	../shlib/$(BFARCH)/libRecoilAnalysis.so 


# ================================================================================
mXFit:  ../tmp/$(BFARCH)/VirVubFitter/mXFit.o mxfit
# --------------------------------
	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/mXFit ../tmp/$(BFARCH)/VirVubFitter/mXFit.o $(GLIBSMXFIT)



clean:
	rm -f  *Dict*
	rm -rf SunWS_cache
	rm -f  ../tmp/$(BFARCH)/VirVubFitter/*.o 
	rm -f  $(addprefix ../tmp/$(BFARCH)/VirVubFitter/, $(OOLIST))
	rm -f  ../shlib/$(BFARCH)/libVirFit.so 
	rm -f  ../shlib/$(BFARCH)/libExclFit.so 
	rm -f  ../shlib/$(BFARCH)/libmXFit.so 
	rm -f  ../bin/$(BFARCH)/VirFit
	rm -f  ../bin/$(BFARCH)/FitOnly
	rm -f  ../bin/$(BFARCH)/ExclFit 
	rm -f  ../bin/$(BFARCH)/b2uFit
	rm -f  ../bin/$(BFARCH)/mXFit

cleanexcl:
	rm -f  {excl,Excl}*Dict*
	rm -f  ../tmp/$(BFARCH)/VirVubFitter/excl*.o
	rm -f  ../shlib/$(BFARCH)/libExclFit.so
	rm -f  ../bin/$(BFARCH)/ExclFit
