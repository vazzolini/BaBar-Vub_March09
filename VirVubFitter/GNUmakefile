# Please see SoftRelTools/HOWTO-GNUmakefile for documentation
# $Id: GNUmakefile,v 1.4 2006/05/18 17:15:37 menges Exp $
#################################################################
#++ library products				[build it with 'lib']

LIBREMOVEFILES := b2uClass.cc b2uFit.cc exclfitNtp.cc $(wildcard *Dict.cc)
LIBTMPLFILES := 
LIBDDLORDERED := 

#################################################################
# hack for makeing package rootcint compatible without *.hh -> *.rdl
override ROOTCINTTARGETS := VirClass.hh CMClass.hh $(wildcard *.rdl)

#################################################################
#++ extra binary products	[not in production, build it with extrabin]

EXTRABINS := FitOnly ExclFit

$(addprefix $(bindir),$(EXTRABINS)): $(bindir)% : %.o

#################################################################
#++ binary products				[build it with 'bin']

BINS := VirFit
BINCCFILES := VirFit.cc $(EXTRABINS:=.cc)

#++ Binary rules		 [in production, build it with 'bin']

$(bindir)VirFit: VirFit.o

#++ shell script products.. 			[build it with 'bin']
#BINSCRIPTS := testscript

#################################################################
#++ include standard makefile from SoftRelTools.
include SoftRelTools/standard.mk
