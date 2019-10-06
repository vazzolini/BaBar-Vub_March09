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
 
ifeq ($(BFARCH),Linux24SL3_i386_gcc323)
# -- RHEL 3
CXX           = g++
CXXFLAGS      = -g  -Wall -fPIC
#CXXFLAGS      = -O2 -Wall -fPIC
LD            = g++
LDFLAGS       = -g
#LDFLAGS       = -O2
SOFLAGS       = -shared
endif
ifeq ($(BFARCH),Linux2)
# -- Linux2
CXX           = g++-2.95.2 
CXXFLAGS      = -g  -Wall -fPIC 
CXXFLAGS      = -O2  -Wall -fPIC -w
LD            = g++-2.95.2
LDFLAGS       = -g
LDFLAGS       = -O2
SOFLAGS       = -shared
endif
ifeq ($(BFARCH),Linux24)
# -- Linux24
CXX           = g++-2.95.3 
CXXFLAGS      = -g  -Wall -fPIC
CXXFLAGS      = -O2  -Wall -fPIC -w
LD            = g++-2.95.3 
LDFLAGS       = -O2
LDFLAGS       = -g
SOFLAGS       = -shared
endif
ifeq ($(BFARCH),SunOS58)
# -- Solaris
CXX           = CC
CXXFLAGS      = -O2 -KPIC 
CXXFLAGS      = -g -KPIC 
LD            = CC
#LD            = purify CC
LDFLAGS       = -O2
LDFLAGS       = -g
SOFLAGS       = -G
endif
ifeq ($(BFARCH),SunOS58_sparc_WS6U1)
# -- Solaris
CXX           = CC
CXXFLAGS      = -O2 -KPIC 
CXXFLAGS      = -g -KPIC 
LD            = CC
#LD            = purify CC
LDFLAGS       = -O2
LDFLAGS       = -g
SOFLAGS       = -G
endif


CXXFLAGS       += $(ROOTCFLAGS)
NLIBS           = $(ROOTGLIBS) 
NLIBS          += -lMinuit

NGLIBS          = $(ROOTGLIBS) 
NGLIBS         += ../shlib/$(BFARCH)/libVubAnalysis.so
NGLIBS         += ../shlib/$(BFARCH)/libRecoilAnalysis.so

NGLIBSFIT       = $(ROOTGLIBS) 
NGLIBSFIT      += -lMinuit
NGLIBSFIT      += ../shlib/$(BFARCH)/libFit.so
NGLIBSFIT      += ../shlib/$(BFARCH)/libRecoilAnalysis.so

NGLIBSCOMP       = $(ROOTGLIBS) 
NGLIBSCOMP      += -lMinuit
NGLIBSCOMP      += ../shlib/$(BFARCH)/libComp.so
NGLIBSCOMP      += ../shlib/$(BFARCH)/libRecoilAnalysis.so

NGLIBSCMP       = $(ROOTGLIBS) 
NGLIBSCMP      += -lMinuit
NGLIBSCMP      += ../shlib/$(BFARCH)/libCmp.so
NGLIBSCMP      += ../shlib/$(BFARCH)/libRecoilAnalysis.so

# CERNLIB STUFF
CERNLIB       = /cern/99/lib/libmathlib.a 
CERNLIB      += /cern/99/lib/libpacklib.a 
CERNLIB      += /cern/99/lib/libkernlib.a 

ifeq ($(BFARCH),Linux24SL3_i386_gcc323)
# -- RHEL 3
CERNLIB       = /cern/2002/lib/libmathlib.a 
CERNLIB      += /cern/2002/lib/libpacklib.a 
CERNLIB      += /cern/2002/lib/libkernlib.a 
CERNLIB      += -lg2c -lfrtbegin
endif

ifeq ($(BFARCH),Linux2)
CERNLIB      += -lg2c 
endif
ifeq ($(BFARCH),Linux24)
CERNLIB      += -lg2c 
endif
ifeq ($(BFARCH),SunOS58)
CERNLIB      += -B dynamic
CERNLIB      += -L/opt/SUNWspro/lib -lM77 -lF77 -lsunmath
endif
ifeq ($(BFARCH),SunOS58_sparc_WS6U1)
CERNLIB      += -B dynamic
CERNLIB      += -L/opt/SUNWspro/lib -lM77 -lF77 -lsunmath
endif
NCLIBS         = $(NGLIBS)
NCLIBS        += $(CERNLIB)

LIBS           = $(filter-out -lNew, $(NLIBS))
GLIBS          = $(filter-out -lNew, $(NGLIBS))
CLIBS          = $(filter-out -lNew, $(NCLIBS))
GLIBSFIT       = $(filter-out -lNew, $(NGLIBSFIT))
GLIBSCOMP      = $(filter-out -lNew, $(NGLIBSCOMP))
GLIBSCMP       = $(filter-out -lNew, $(NGLIBSCMP))

# ======================================================================
NTUPLES = valHist.o recoilDSys.o recoilDSysDict.o b2uQA.o sHist.o mesFit.o util.o functions.o
NTUPLES = b2uQA.o b2uQADict.o sHist.o sHistDict.o mesFit.o mesFitDict.o util.o functions.o

FIT = fitNtp.o fitNtpDict.o recoilDSys.o recoilDSysDict.o

COMP = theanal.o theanalDict.o 

CMP = thecomp.o thecompDict.o 

ABCFILES = abcfit_smear.o abcfit_aibi_evol.o abcfit_bmatrix.o abcfit_interface_vub.o ntparam.o

ICCFILES = bookHist.icc  fillHist.icc  fitUtil.icc   init.icc     splitOff.icc    util.icc \
           fastHist.icc  fitInit.icc   gamStudy.icc  mcTruth.icc  trackStudy.icc


# ======================================================================
$(addprefix ../tmp/$(BFARCH)/VubAnalysis/, %.o) : %.cc 
	$(CXX) $(CXXFLAGS) -I../ -c $< -o $@


# -- Uncomment this if you want recoilNtp to have a dependency on icc files. 
#    Will lead to automatic recompilation ...
#recoilNtp.cc:  recoilNtp.hh $(ICCFILES)
#	$(CXX) $(CXXFLAGS) -I../ -c $@ -o ../tmp/$(BFARCH)/VubAnalysis/recoilNtp.o

# ================================================================================
b2uQA:  $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(NTUPLES)) 
#-----------------------------------------------------------------------
	$(CXX) $(SOFLAGS) $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(NTUPLES)) -o ../shlib/$(BFARCH)/libVubAnalysis.so
	$(CXX) $(CXXFLAGS)  -c runB2uQA.cc -o ../tmp/$(BFARCH)/VubAnalysis/runB2uQA.o
	$(CXX) $(CXXFLAGS) -DFAST -I../ -c ../RecoilAnalysis/recoilDSys.cc -o ../tmp/$(BFARCH)/RecoilAnalysis/recoilDSys.o
	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/runB2uQA \
	../tmp/$(BFARCH)/VubAnalysis/runB2uQA.o \
	../tmp/$(BFARCH)/RecoilAnalysis/recoilDSys.o $(LIBS) $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(NTUPLES))

b2uQADict.cc: b2uQA.hh 
	$(ROOTSYS)/bin/rootcint  -f b2uQADict.cc -c b2uQA.hh 

mesFitDict.cc: mesFit.hh 
	$(ROOTSYS)/bin/rootcint  -f mesFitDict.cc -c mesFit.hh 

sHistDict.cc: sHist.hh 
	$(ROOTSYS)/bin/rootcint  -f sHistDict.cc -c sHist.hh 



# ================================================================================
fit:  $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(FIT))
#-------------------------------------------------------
	$(CXX) $(SOFLAGS) $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(FIT)) -o ../shlib/$(BFARCH)/libFit.so

fitNtpDict.cc:  fitNtp.hh  
	$(ROOTSYS)/bin/rootcint -f fitNtpDict.cc -c -I../ fitNtp.hh 


# ================================================================================
comp:  $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(COMP))
#-------------------------------------------------------
	$(CXX) $(SOFLAGS) $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(COMP)) -o ../shlib/$(BFARCH)/libComp.so

theanalDict.cc:  theanal.hh
	$(ROOTSYS)/bin/rootcint -f theanalDict.cc -c -I../ theanal.hh

cmp:  $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(CMP))
#-------------------------------------------------------
	$(CXX) $(SOFLAGS) $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(CMP)) -o ../shlib/$(BFARCH)/libCmp.so

thecompDict.cc:  thecomp.hh
	$(ROOTSYS)/bin/rootcint -f thecompDict.cc -c -I../ thecomp.hh

# ================================================================================
#anaRecoil: ../tmp/$(BFARCH)/VubAnalysis/anaRecoil.o lib
# ---------------------------------------------------------
#	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/anaRecoil ../tmp/$(BFARCH)/VubAnalysis/anaRecoil.o $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(ABCFILES)) $(CLIBS) 

# ================================================================================
anaQA: anaQA.cc valHist.cc recoilNtp.cc 
# ----------------------------------------------------------------------------------
#	rm -f ../tmp/$(BFARCH)/VubAnalysis/valHist.o
#	rm -f ../tmp/$(BFARCH)/VubAnalysis/anaRecoil.o
	rm -f ../tmp/$(BFARCH)/VubAnalysis/recoilNtp.o 
#	rm -f ../tmp/$(BFARCH)/VubAnalysis/recoilDSys.o 
	gmake -f main.mk
	$(CXX) $(CXXFLAGS) -DFAST -I../ -c valHist.cc   -o ../tmp/$(BFARCH)/VubAnalysis/valHist.o  
	$(CXX) $(CXXFLAGS) -DFAST -I../ -c anaQA.cc -o ../tmp/$(BFARCH)/VubAnalysis/anaQA.o
	$(CXX) $(CXXFLAGS) -DFAST -I../ -c recoilNtp.cc -o ../tmp/$(BFARCH)/VubAnalysis/recoilNtp.o
	$(CXX) $(CXXFLAGS) -DFAST -I../ -c recoilDSys.cc -o ../tmp/$(BFARCH)/VubAnalysis/recoilDSys.o
	$(LD) $(LDFLAGS) -o ../bin/$(BFARCH)/anaQA \
	../tmp/$(BFARCH)/VubAnalysis/anaQA.o \
	../tmp/$(BFARCH)/VubAnalysis/recoilNtp.o ../tmp/$(BFARCH)/VubAnalysis/valHist.o ../tmp/$(BFARCH)/VubAnalysis/recoilDSys.o \
	$(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(ABCFILES)) ../shlib/$(BFARCH)/libRecoilAnalysis.so $(LIBS) $(CERNLIB)


# ================================================================================
#abcRecoil: ../tmp/$(BFARCH)/VubAnalysis/anaRecoil.o lib 
# ---------------------------------------------------------
#	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/anaRecoil ../tmp/$(BFARCH)/VubAnalysis/anaRecoil.o $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(ABCFILES)) $(CLIBS) 


# ================================================================================
VubFit:  ../tmp/$(BFARCH)/VubAnalysis/VubFit.o fit
# --------------------------------
	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/VubFit ../tmp/$(BFARCH)/VubAnalysis/VubFit.o $(GLIBSFIT) 

# ================================================================================
compplot:  ../tmp/$(BFARCH)/VubAnalysis/compplot.o cmp
# --------------------------------
#	rm -f ../tmp/$(BFARCH)/VubAnalysis/compplot.o 
	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/compplot ../tmp/$(BFARCH)/VubAnalysis/compplot.o $(GLIBSCMP) 
# ================================================================================
anaplot:  ../tmp/$(BFARCH)/VubAnalysis/anaplot.o comp
# --------------------------------
#	rm -f ../tmp/$(BFARCH)/VubAnalysis/anaplot.o 
	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/anaplot ../tmp/$(BFARCH)/VubAnalysis/anaplot.o $(GLIBSCOMP) 

# ================================================================================
VubAnalysisApp: ../tmp/$(BFARCH)/VubAnalysis/VubAnalysisApp.o lib
# ---------------------------------------------------------
	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/VubAnalysisApp ../tmp/$(BFARCH)/VubAnalysis/VubAnalysisApp.o $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(ABCFILES)) $(CLIBS) 

# ================================================================================
b2u: ../tmp/$(BFARCH)/VubAnalysis/b2u.o
# ---------------------------------------------------------------
	gmake -f main.mk
#	$(CXX) $(CXXFLAGS) -I../ -c b2u.cc -o ../tmp/$(BFARCH)/VubAnalysis/b2u.o
	$(CXX) $(CXXFLAGS) -I../ -c b2uNtp.cc -o ../tmp/$(BFARCH)/VubAnalysis/b2uNtp.o
	$(LD) $(LDFLAGS)  -o ../bin/$(BFARCH)/b2u ../tmp/$(BFARCH)/VubAnalysis/b2u.o ../tmp/$(BFARCH)/VubAnalysis/b2uNtp.o ../shlib/$(BFARCH)/libRecoilAnalysis.so $(addprefix ../tmp/$(BFARCH)/VubAnalysis/, $(ABCFILES)) $(LIBS) $(CERNLIB)


# ================================================================================
Killer:  Killer.o  lib
# --------------------------------
#	$(CXX) $(CXXFLAGS) -c $< 
	$(LD) $(LDFLAGS)  -o Killer Killer.o $(ABCFILES)  $(CLIBS) 

# ================================================================================
all:  
#---------------------------------------------------------------
	cd ../RecoilAnalysis; gmake -f Makefile lib; cd ../VubAnalysis
	gmake -f main.mk
	gmake -f Makefile lib
#	gmake -f Makefile anaRecoil
	gmake -f Makefile VubAnalysisApp
#	gmake -f Makefile anaQA


# ================================================================================
clean:
	rm -f  *Dict*
	rm -rf SunWS_cache
	rm -f  ../tmp/$(BFARCH)/VubAnalysis/*.o 
	rm -f  ../shlib/$(BFARCH)/libVubAnalysis.so 
	rm -f  ../shlib/$(BFARCH)/libVubFastAnalysis.so 
	rm -f  ../shlib/$(BFARCH)/libComp.so 
	rm -f  ../bin/$(BFARCH)/test.exe 
	rm -f  ../bin/$(BFARCH)/VubFit 
	rm -f  ../bin/$(BFARCH)/Killer 
#	rm -f  ../bin/$(BFARCH)/anaRecoil
	rm -f  ../bin/$(BFARCH)/anaQA
	rm -f  ../bin/$(BFARCH)/VubAnalysisApp
	rm -f  ../bin/$(BFARCH)/anaplot
	rm -f  ../bin/$(BFARCH)/compplot
	rm -f  ../bin/$(BFARCH)/runB2uQA
	rm -f  ../bin/$(BFARCH)/b2u

# ================================================================================
cleanall:
	gmake -f Makefile clean
	cd ../RecoilAnalysis; gmake -f Makefile clean; cd ../VubAnalysis
