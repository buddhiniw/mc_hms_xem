## This makefile must be executed with gmake (gnu make).

## These will have to be modified when setting up your own simc.  They
## point to the software necessary to run simc.

## ARGONNE DEFAULT SETUP FLAGS:
#simcdir = .
#Csoft = /disk1/users/reinhold/Csoft/05Dec1996
#CERN_ROOT = /disk1/lib/cern/97a

## CEBAF DEFAULT SETUP FLAGS: ALSO NEED TO UNCOMMENT CERN_ROOT DEFNS FOR HP/SUN/...
simcdir = .
Csoft = /group/hallc/Csoft/Analyzer


## THE REST SHOULD BE OK WITHOUT MODIFICATION.

RM        = rm -f 
SHELL     = /bin/sh

SPEC      = ./hms_simc/

SHARE	  = ./shared/

my_objs = $(SPEC)mc_hms.o $(SPEC)mc_hms_hut.o $(SPEC)mc_hms_recon.o \
	  $(SHARE)project.o $(SHARE)transp.o $(SHARE)rotate_haxis.o \
	  $(SHARE)locforunt.o $(SHARE)musc.o $(SHARE)musc_ext.o \
	  $(SHARE)stringlib.o mt19937.o gauss1.o loren.o

sub_objs = $(SPEC)mc_hms.o $(SPEC)mc_hms_hut.o $(SPEC)mc_hms_recon_sub.o \
	  $(SHARE)project.o $(SHARE)transp_sub.o $(SHARE)rotate_haxis.o \
	  $(SHARE)locforunt.o $(SHARE)musc.o $(SHARE)musc_ext.o \
	  $(SHARE)stringlib.o mt19937.o gauss1.o loren.o

MYOS := $(subst -,,$(shell uname))
#mess $(MYOS)

CERNLIBS = -lgeant$(GEANTVER) -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib

ifeq ($(MYOS),HPUX)
  ifneq (,$(findstring 09,$(shell uname -r)))
    HPUXVERSION := 09
  else
    HPUXVERSION := 10
  endif
  LIBROOT = $(Csoft)/../$(MYOS)$(HPUXVERSION)/lib
else
  LIBROOT = $(Csoft)/../$(MYOS)/lib
endif

ifeq ($(MYOS),HPUX)
  CERN_ROOT = /site/cernlib/hp700_ux90/96a
  FFLAGS=+U77 +ppu -C +e +es +FPVZOU -O +Onolimit -R8
  LDFLAGS=-Wl,-a archive
  OTHERLIBS = \
	-Wl,-L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS) \
	-Wl,-L/usr/lib/X11R5 -lX11 -lm
endif

#ifeq ($(MYOS),Linux)
#  FFLAGS=-I$(Csoft)/SRC/INCLUDE
#  OTHERLIBS = -L$(CERN_ROOT)/lib -lpacklib -lc -lm
#  F77=Af77
#endif                                                                                        
# need to define: setenv ABSOFT /apps/absoft/PRO/usr/absoft
#
# meaning of flags for ABSOFT compiler
#
# -O  (-g) generate optimized(debuggable) code
# -V  accept source code in VAX TAB format
# -f  fold all symbolic names to lower case
# -W  accept wide format, 132 columns (but not more tha that
# -s  all program storage is treated as static
# -B108 apply a single underscore to the names of subroutines and functions
# -B100 (B101) generate code optimized for Pentium II or higher that
#              will not (will still) run on a Pentium
# -N1 all variables initalized with DATA statements are stored as static
# -N2 use real*8 or complex*16 intrinsic functions, overriding real*4 and
#     complex specifications
# -N22 append an undersocre to the names of common blocks instead of prepending
#      _C
# -N90 pass arguments in the same manner as g77 and f2c(g2c)
# -N113 promote real and complex data types without explicit type/length
#       declaration to real*8 and complex*16
#
ifeq ($(MYOS),Linux)
# Uncomment these two lines for redhat 7.2
#  ABSOFT=/apps/absoft/PRO/usr/absoft
#  CERN_ROOT=/apps/cernlib/i386_redhat72/2001
# Uncomment these two lines for Enterprise Linux 3.
#  ABSOFT=/apps/absoft/PRO/opt/absoft
  ABSOFT=/apps/absoft/absoft-8.2/opt/absoft
  CERN_ROOT = /apps/cernlib/i386_rhel3/2003
  FABSFLAGS=-O -V -W -f -s -N1 -B108 -B100 -N90 -N22 -N2 -N113
  INCLUDES=-I.,..,./sos,$(SPEC),./hrsr,./hrsl,./shms,$(Csoft)/SRC/INCLUDE
  EXTRAFLAGS=-DABSOFTFORTRAN
  FFLAGS= $(INCLUDES) $(FABSFLAGS) $(EXTRAFLAGS)
  FFLAG1=$(FFLAGS) -c
  OTHERLIBS = -L$(LIBROOT) -lctp \
        -L$(CERN_ROOT)/lib $(CERNLIBS) -lV77 -lU77 -lg2c -lc -lm \
        -lnsl -lcrypt
# Uncomment the last line above (-lnsl -lcrypt) if you are using cernlib 2001.
  FC  := $(ABSOFT)/bin/f77
  F77 :=$(ABSOFT)/bin/f77
endif

ifeq ($(MYOS),ULTRIX)
  FFLAGS=-check_bounds
  LDFLAGS=
  OTHERLIBS = -L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS)
endif

ifeq ($(MYOS),SunOS)
  CERN_ROOT = /site/cernlib/sun4_solaris2/97a
  FFLAGS=-e -O -I$(Csoft)/SRC/INCLUDE -ftrap=underflow
  ifeq ($(OSTYPE),SunOS4)
    OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -lnsl -lX11
  else
    OTHERLIBS = -L$(CERN_ROOT)/lib $(CERNLIBS) -lnsl -lsocket -lX11
  endif
endif

ifeq ($(MYOS),AIX)
  F77=f77
  FFLAGS=-g -qfixed=132 -qextname -O -I$(Csoft)/SRC/INCLUDE
  OTHERLIBS = -L$(CERN_ROOT)/lib -lpacklib $(CERNLIBS) -lX11
endif

ifeq ($(MYOS),OSF1)
  F77=f77
  CERN_ROOT = /disk1/lib/cern/new
  LIBROOT = $(Csoft)/OSF1/lib
  FFLAGS= -r8 -extend_source -Wl,-taso -I -warn argument_checking \
        -warn declarations -warn truncated_source -warn unused
  LDFLAGS= 
  OTHERLIBS = -Wl,-L$(CERN_ROOT)/lib \
        -lpacklib $(CERNLIBS) -Wl,-L/usr/lib/X11R5 -lX11 -lm 
endif

%.o: %.f
	$(F77) $(FFLAGS) -c $< -o $@

all: mc_hms_single

mc_hms_single: $(my_objs) Makefile mc_hms_single.o
	$(F77) -o $@ $(FFLAGS) mc_hms_single.o $(my_objs) $(OTHERLIBS)

clean:
	rm -f *.o $(SPEC)*.o $(SHARE)*.o mc_hms_single
