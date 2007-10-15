# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

#=============================================================
#  LAPS environment:
#=============================================================

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

RM = rm -f
FLAG = $(FFLAGS) -bnoquiet #-bloadmap

LAPSLIBS = $(SRCROOT)/src/lib/wind/libwind.a \
           $(SRCROOT)/src/wind/libwind.a \
           $(SRCROOT)/src/lib/temp/libtemp.a \
           $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a

DBLAPSLIBS = $(SRCROOT)/src/lib/wind/libwind.a \
             $(SRCROOT)/src/wind/libwind.a \
             $(SRCROOT)/src/lib/temp/libdebugtemp.a \
             $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a

#LAPSFC = $(SRCROOT)/src/sfc/libsfc.a

BUFRLIB = $(SRCROOT)/src/var/bufr/libbufr.a

#WINDLIB = $(SRCROOT)/src/wind/libwind.a

#HUMIDLIB = $(SRCROOT)/src/humid/liblq3.a

.SUFFIXES:
.SUFFIXES: .o .c .f90 .a .f
.f90.o:;    $(FC) $(FLAG) $(FREE) -c $*.f90
.f.o:;      $(FC) $(FLAG) -c $*.f
.c.o:;	    $(CC) $(CPPFLAGS) -c $*.c

#=============================================================
# STMAS_4D:
#=============================================================

ANA = ANALYSIS.exe

STMAS= \
        para_stmas.f90 \
        general_tools.f90 \
        input_bg_and_obs.f90 \
        prep_stmas.f90 \
        post_stmas.f90 \
        lbfgsb.o \
        stmas_core.f90 \
        output_analysis.f90 \
        main.f90

all: $(ANA) 
#$(POST)

$(ANA): $(STMAS) $(LAPSLIBS)
	$(FC) $(FFLAGS) $(FREE) $(STMAS) $(LAPSFC) $(LAPSLIBS) \
	$(BUFRLIB) $(OTHERLIBS) -o $@
	rm -f *.o

$(POST): $(GSIPO) $(LAPSLIBS) $(HUMIDLIB)
	$(FC) $(FFLAGS) $(FREE) $(GSIPO) $(LAPSFC) $(LAPSLIBS) \
	$(BUFRLIB) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(STMAS)
	$(RM) $(ANA);
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS) $(BUFRLIB)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: $(ANA)
	$(INSTALL) $(ANA) $(INSTALLROOT)/bin/$(ANA)

clean:
	rm -f $(ANA) $(POST)  *.o *~ *# *.mod

distclean: clean

