# Makefile for dataMCcomp
# Antonio 15-09-2008
#-------------------------------------

# ROOT Include Path and libs

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)
ROOTGLIBS     += $(ROOTSYS)/lib/libHtml.so 
ROOTGLIBS     += $(ROOTSYS)/lib/libMinuit.so 

# RooFit	

ROOFITCOREINCLUDE 	= -I../../RooFitCore/tmp/
ROOFITMODELSINCLUDE 	= -I../../RooFitModels/tmp/
ROOFITCORESHLIB 	= -L../../RooFitCore/tmp/
ROOFITMODELSSHLIB 	= -L../../RooFitModels/tmp/
ROOFITSHLIB		+= $(ROOFITCORESHLIB)
ROOFITSHLIB		+= $(ROOFITMODELSSHLIB)

# -- Linux so far --
CXX 		= g++
CXXFLAGS	= -ggdb -Wall
LD		= g++
LDFLAGS		= -ggdb
SOFLAGS		= -ggdb -shared

CXXFLAGS += $(ROOTCFLAGS)
CXXFLAGS += $(ROOFITCOREINCLUDE)
CXXFLAGS += $(ROOFITMODELSINCLUDE)

# ----

OBJECTS		= main.o toyclass.o RooThorstenSig.o Dictionary.o 

all:	$(OBJECTS) 
	$(CXX) -o toygen $(CXXFLAGS) $(ROOTGLIBS) $(ROOFITSHLIB) -lRooFitCore -lRooFitModels $(OBJECTS)

shlib:	$(OBJECTS)
	$(CXX) $(SOFLAGS) $(ROOTGLIBS) $(OBJECTS)

main.o:	main.C 
	$(CXX) -c $(CXXFLAGS) main.C

# util.o:	util.hh util.cc
# 	$(CXX) -c $(CXXFLAGS) util.cc

toyclass.o: toyclass.h toyclass.C
	$(CXX) -c $(CXXFLAGS) toyclass.C

RooThorstenSig.o: RooThorstenSig.hh RooThorstenSig.cc
	$(CXX) -c $(CXXFLAGS) RooThorstenSig.cc

# RooCCB.o: RooCCB.hh RooCCB.cc
# 	$(CXX) -c $(CXXFLAGS) RooCCB.cc

Dictionary.cxx:  RooThorstenSig.hh toyclass.h
	$(ROOTSYS)/bin/rootcint Dictionary.cxx -c $(CXXFLAGS) RooThorstenSig.hh toyclass.h

Dictionary.o: Dictionary.cxx	
	$(CXX) -c $(CXXFLAGS) Dictionary.cxx

clean:
	rm -f *.o
	rm -f Dictionary.cxx
	rm -f *Dict*
	rm toygen