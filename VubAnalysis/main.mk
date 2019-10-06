# ================================================================
# Default make file for non alpha fortan
# To run: gmake -f filename.mk
#
# libraries are linked in the following order:
# ULIB CERNLIBS
# ================================================================
# name of executable
NAME = test.exe
ifeq ($(BFARCH),Linux24SL3_i386_gcc323)
FC	=g77
else
FC	=f77
endif
# executable is stored in EXEDIR (current working directory)
EXEDIR = $(PWD)
EXEDIR = ../bin/$(BFARCH)
# *.o files are stored in OBJDIR (current working directory)
OBJDIR  = $(PWD)
OBJDIR  =  ../tmp/$(BFARCH)/VubAnalysis

#Build file and directory lists



SRCDIR1 = .
FFILES1 =  abcfit_bmatrix.F abcfit_aibi_evol.F ntparam.F
UINC1 = $(SRCDIR1)


SRCDIR2 =  .
FFILES2 = abcfit_smear.F abcfit_interface_vub.F main.F
UINC2 = $(SRCDIR2)

# Additional include directory
UINCADD = .
#Make full lists
FFILES = $(FFILES1) $(FFILES2) 
SRCDIR = $(SRCDIR1) $(SRCDIR2) 
UINC =   -I$(UINC1) -I$(UINC2)
#
# ===============================================================
# set FORTRAN flags

FFLAGS = $(FCOPT) $(UINC)

## Architechture dependency
ifeq ($(BFARCH),Linux2)
FCOPT = -g -fno-automatic -fvxt -fdollar-ok -fno-backslash \
        -ffixed-line-length-132 -fno-second-underscore
FCOPT = -g -O -fno-automatic -fvxt -fdollar-ok -fno-backslash \
        -ffixed-line-length-132 -fno-second-underscore
endif
ifeq ($(BFARCH),Linux24)
FCOPT = -g -fno-automatic -fvxt -fdollar-ok -fno-backslash \
        -ffixed-line-length-132 -fno-second-underscore
FCOPT = -g -O -fno-automatic -fvxt -fdollar-ok -fno-backslash \
        -ffixed-line-length-132 -fno-second-underscore
endif
ifeq ($(BFARCH),Linux24SL3_i386_gcc323)
FCOPT = -g -fno-automatic -fvxt -fdollar-ok -fno-backslash \
        -ffixed-line-length-132 -fno-second-underscore
FCOPT = -g -O -fno-automatic -fvxt -fdollar-ok -fno-backslash \
        -ffixed-line-length-132 -fno-second-underscore
endif
ifeq ($(BFARCH),SunOS58)
FCOPT = -lm -lsocket -lnsl
FCOPT = -lm -lsocket -lnsl -O2
endif
ifeq ($(BFARCH),SunOS58_sparc_WS6U1)
FCOPT = -lm -lsocket -lnsl
FCOPT = -lm -lsocket -lnsl -O2
endif

# Search path for source and object files:
vpath %.F $(SRCDIR)
vpath %.f $(SRCDIR)
vpath %.o $(OBJDIR)

O_FILES  := $(FFILES:%.F=%.o) 
#OBJFILES := $(addprefix $(OBJDIR)/,$(notdir $(O_FILES)))
OBJFILES = $(addprefix $(OBJDIR)/, $(O_FILES))

# Make full list of source files
#SRCFILES := $(addprefix $(SRCDIR)/, $(FFILES))

$(addprefix $(OBJDIR)/, %.o) : %.F
	$(FC) $(FFLAGS) $< -c -o $@

# Make executable
$(EXEDIR)/$(NAME):  $(OBJFILES) 
	@echo " "
	@echo " L I N K I N G ...."
	@echo " "
#	$(FC) -o $@  $(OBJFILES) \
#        /cern/99/lib/libmathlib.a /cern/99/lib/libpacklib.a\
#        /cern/99/lib/libkernlib.a  $(FFLAGS)









