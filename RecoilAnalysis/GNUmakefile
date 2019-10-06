# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile,v 1.3 2006/09/01 16:58:12 menges Exp $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := $(wildcard *Dict.cc)
LIBTMPLFILES := 
LIBDDLORDERED := 

#################################################################
# hack for makeing package rootcint compatible without *.hh -> *.rdl
override ROOTCINTTARGETS := mesData.hh mesFIT.hh PIDTable.hh TRKTable.hh recoilDsys recoilAnalysis.hh PIDData.hh TRKData.hh

#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS := 

$(addprefix $(bindir),$(EXTRABINS)): $(bindir)% : %.o

#################################################################
#++ binary products				[build it with 'bin']

BINS := exampleAnalysisApp 
BINCCFILES := exampleAnalysisApp.cc $(EXTRABINS:=.cc)

#++ Binary rules		 [in production, build it with 'bin']

$(bindir)exampleAnalysisApp: exampleAnalysisApp.o

#++ shell script products.. 			[build it with 'bin']
#BINSCRIPTS := testscript

#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
