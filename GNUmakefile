ifndef ATMPD_ROOT
	ATMPD_ROOT = ../../..

# define the following if you want to generate version of fitqun that
# does not use sk libraries (for hyperk)
NOSKLIBRARIES = 1
endif

# Uncomment if you want to enable speed optimization
YORI_HAYAKU = 1

FMATH_OPT = -O3
#FMATH_OPT += -DNDEBUG -D_FILE_OFFSET_BITS=64 -msse2 -mfpmath=sse -fomit-frame-pointer -fno-operator-names -ffast-math -march=core2
GCC_VER=$(shell gcc -dumpversion)

ifeq ($(FC),gfortran)
  ifeq ($(shell expr $(GCC_VER) \>= 4.2),1)
    FMATH_OPT += -mtune=core2 #-std=c++11
  endif
  ifeq ($(shell expr $(GCC_VER) \>= 4.5),1)
    FMATH_OPT += -fexcess-precision=fast
  endif
endif

ifndef NOSKLIBRARIES
include $(ATMPD_ROOT)/config.gmk
LOCAL_LIBS = -laplib -lseplib -lodlib
else
include $(FITQUN_ROOT)/config.gmk
LOCAL_LIBS = -L$(WCSIMDIR) -lWCSim
endif

LOCAL_LIBS += -L$(FITQUN_ROOT)
CXXINCLUDES += -I$(FITQUN_ROOT)
CINCLUDES += -I$(FITQUN_ROOT)

#ifeq ($(FC),g77)
#  SITE_DEFINES += -DFQ_G77
#endif

ROOTLIBS     := $(shell root-config --glibs) -lTreePlayer -lMinuit

SITE_DEFINES += -DHEMI_CUDA_DISABLE

ifdef YORI_HAYAKU
  SITE_DEFINES += -DYORI_HAYAKU
#  SKDEBUGFLAGS = -g -O
  CXXDEBUGFLAGS =
  CXXFLAGS = $(CXXDEBUGFLAGS) $(FMATH_OPT) $(CXXOPTIONS) $(CXXINCLUDES) $(CXXDEFINES)
endif




#
#  Objects
#

LIBNAME = fiTQun

INCFILES = fitqunout.h fitqunoutC.h spliTChanOut.h spliTChanOutC.h

ifndef NOSKLIBRARIES
OBJS = $(patsubst %.o,$(FITQUN_ROOT)/%.o,TRuntimeParameters.o TParametersOptionManager.o TUnitsTableParser.o fiTQun.o fiTQun_shared.o fastmath.o fQEventManager.o fQChrgPDF.o spliTChan.o SKTVarConsts.o fQROOTOut.o clrfqcmns.o TScatTable.o TScatTableF.o scattabledict.o scattabledictF.o fortconsts.o fortinit.o fortread.o restoreTQBanks.o fillfqzbsbank.o readfqzbsbank.o trginfo.o getfakepi0.o filldstvars.o fillnt_util.o r_dum.o PDK_MuGamma.o apflscndprt.o)
OBJS += fq_NLL_scanner.o
else
OBJS = $(patsubst %.o,$(FITQUN_ROOT)/%.o, TRuntimeParameters.o TParametersOptionManager.o TUnitsTableParser.o fiTQun.o fiTQun_shared.o fastmath.o fQEventManager.o fQChrgPDF.o spliTChan.o fQROOTOut.o clrfqcmns.o TScatTable.o TScatTableF.o scattabledict.o scattabledictF.o WCSimWrap.o sortzv.o PDK_MuGamma.o)
OBJS += fq_NLL_scanner.o
endif

ifndef NOSKLIBRARIES
EXES = run_fq_NLL_scanner
else
EXES = run_fq_NLL_scannerWC
endif

#
#  Rules for building library 
#

ifndef NOSKLIBRARIES
.PHONY:  lib$(LIBNAME).a $(LIBDIR)lib$(LIBNAME)
else
.PHONY:  lib$(LIBNAME).a
endif

.PHONY:  clean setup includes install.includes depend lib install.lib exec install.exec

all: ${EXES}

install.all::  setup includes install.includes depend lib install.lib exec install.exec

lib$(LIBNAME).a : $(OBJS)
	$(RM) $@
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

ifndef NOSKLIBRARIES
$(LIBDIR)lib$(LIBNAME).a : lib$(LIBNAME).a
	$(RM) $@
	$(INSTALL_LIB) $< $@
endif

run_fq_NLL_scanner: run_fq_NLL_scanner.o lib$(LIBNAME).a
	$(CXX) -g -Wall -o $@ $@.o lib$(LIBNAME).a $(LDLIBS)

ifdef NOSKLIBRARIES
run_fq_NLL_scannerWC: run_fq_NLL_scanner
	mv run_fq_NLL_scanner run_fq_NLL_scannerWC


endif

functest: functest.o lib$(LIBNAME).a
	$(CXX) -g -Wall -o $@ $@.o lib$(LIBNAME).a $(LDLIBS)

speedtest: speedtest.o lib$(LIBNAME).a
	$(CXX) -g -Wall -o $@ $@.o lib$(LIBNAME).a $(LDLIBS)

scattabledict.cc: TScatTable.h
	@echo "Generating ScatTable Dictionary..."
#	@rootcint scattabledict.cc -c TScatTable.h
	@rootcint -f scattabledict.cc -c TScatTable.h TScatTableLinkDef.h

scattabledictF.cc: TScatTableF.h
	@echo "Generating ScatTable Dictionary..."
#	@rootcint scattabledictF.cc -c TScatTableF.h
	@rootcint -f scattabledictF.cc -c TScatTableF.h TScatTableFLinkDef.h

sortzv.o: sortzv.F
	gfortran -c sortzv.F 

clean::
	$(RM) -f *.o *.a *~ ${EXES} scattabledict*.cc scattabledict*.h

# Setup will download the fiTQun tunes for SK.
ifndef NOSKLIBRARIES
setup::
#	sh sourceme
else
setup::
endif

lib:: lib$(LIBNAME).a

install.lib:: $(LIBDIR)lib$(LIBNAME).a

exec:: ${EXES}

includes:: $(INCFILES)

install.includes:: $(INCFILES)
	$(INSTALL_INC) -t $(FORTRAN_INCDIR) $(INCFILES) 

ifndef NOSKLIBRARIES
clean::
	list='$(INCFILES)'; for file in $$list; do \
	( $(RM) $(FORTRAN_INCDIR)/"$$file" ); \
	done
endif

# Get tune files on compilation.
#.PHONY: get_tune

#get_tune
