# ================================================================
# Default make file for non alpha fortan
# To run: gmake -f filename.mk
#
# libraries are linked in the following order:
# ULIB CERNLIBS
# ================================================================
# name of executable
NAME = fillNtparam
ifeq ($(BFARCH),SunOS58)
FC	=f77
endif
ifeq ($(BFARCH),Linux24)
FC	=g77
endif
ifeq ($(BFARCH),Linux2)
FC	=g77
endif

# executable is stored in EXEDIR (current working directory)
EXEDIR = ../bin/$(BFARCH)
EXEDIR = .
# *.o files are stored in OBJDIR (current working directory)
OBJDIR  =  ../tmp/$(BFARCH)/VubAnalysis
OBJDIR  =  .


SRCDIR2 =  $(PWD)
FFILES2 = fill_ntparam.F ntparam.F
UINC2 = $(SRCDIR2)


# Additional include directory
UINCADD = .
#Make full lists
FFILES = $(FFILES2) 
SRCDIR = $(SRCDIR2)
UINC =   -I$(UINC2)
#
# ===============================================================
# set FORTRAN flags



ifeq ($(BFARCH),SunOS58)
FCOPT = -lm -lsocket -lnsl
endif

ifeq ($(BFARCH),Linux24)
FCOPT = -g -fno-automatic -fvxt -fdollar-ok -fno-backslash \
        -ffixed-line-length-132 -fno-second-underscore
endif

ifeq ($(BFARCH),Linux2)
FCOPT = -g -fno-automatic -fvxt -fdollar-ok -fno-backslash \
        -ffixed-line-length-132 -fno-second-underscore
endif

FFLAGS = $(FCOPT) $(UINC)


# Search path for source and object files:
vpath %.F $(SRCDIR)
vpath %.f $(SRCDIR)
vpath %.o $(OBJDIR)

OBJFILES =  fill_ntparam.o ntparam.o

# Make full list of source files
SRCFILES := $(addprefix $(SRCDIR)/, $(FFILES))

# Rules...
# Make object files
#
#$(OBJDIR)/%.o : %.F
#	$(FC) $(FFLAGS) $< -o $@

# Make executable
$(EXEDIR)/$(NAME):  $(OBJFILES) 
	@echo " "
	@echo " L I N K I N G ...."
	@echo " "
	$(FC) -o $@  $(OBJFILES) \
        /cern/99/lib/libmathlib.a /cern/99/lib/libpacklib.a\
        /cern/99/lib/libkernlib.a  $(FFLAGS)


clean:
	/bin/rm -f fillNtparam fill_ntparam.o ntparam.o






