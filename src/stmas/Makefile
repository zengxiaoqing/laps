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
FLAG = $(FFLAGS)

LAPSLIBS = $(SRCROOT)/src/lib/wind/libwind.a \
           $(SRCROOT)/src/wind/libwind.a \
	   $(SRCROOT)/src/balance/libqbalpe.a \
           $(SRCROOT)/src/lib/temp/libtemp.a \
	   $(SRCROOT)/src/humid/liblq3.a \
           $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a \
	   $(SRCROOT)/src/lib/radar/moving/libairad.a 

DBLAPSLIBS = $(SRCROOT)/src/lib/wind/libdebugwind.a \
	     $(SRCROOT)/src/balance/libqbalpe.a \
             $(SRCROOT)/src/lib/temp/libdebugtemp.a \
	     $(SRCROOT)/src/humid/liblq3.a \
             $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a \
	     $(SRCROOT)/src/lib/radar/moving/libdebugairad.a

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

#ANA = ANALYSIS.exe
ANA = STMAS3D.exe
ANA_CG=STMAS3D_CG.exe

STMAS= \
        prmtrs_stmas.f90 \
        GPSwetdelays.f90 \
        uniform_interpolation.f90 \
        generaltools.f90 \
        drawcountour.f90 \
        read_backgrd.f90 \
        readobserves.f90 \
        read_laps_radar.f90 \
        input_bg_obs.f90 \
        prep_stmas4d.f90 \
        post_stmas4d.f90 \
        wcompt_gradt.f90 \
        smcostf_grad.f90 \
        gsbcost_grad.f90 \
        hydcost_grad.f90 \
        costfun_grad.f90 \
        stmas4d_core_lbfgs.f90 \
        output_anals.f90 \
        main.f90

STMAS_CG= \
        prmtrs_stmas.f90 \
        generaltools.f90 \
        drawcountour.f90 \
        read_backgrd.f90 \
        readobserves.f90 \
        input_bg_obs.f90 \
        prep_stmas4d.f90 \
        post_stmas4d.f90 \
        wcompt_gradt.f90 \
        smcostf_grad.f90 \
        gsbcost_grad.f90 \
        hydcost_grad.f90 \
        costfun_grad.f90 \
        cg_descent.o \
        stmas4d_core_cg.f90 \
        output_anals.f90 \
        main.f90

all: $(ANA) 
#$(ANA_CG): compilation problem for cg_descent.f on toro.
#$(POST)

$(ANA): $(STMAS) $(LAPSLIBS)
	$(FC) $(FFLAGS) $(FREE) $(STMAS) $(LAPSFC) $(LAPSLIBS) \
	$(BUFRLIB) $(OTHERLIBS) -o $@
	rm -f *.o

$(ANA_CG): $(STMAS_CG) $(LAPSLIBS)
	$(FC) $(FFLAGS) $(FREE) $(STMAS_CG) $(LAPSFC) $(LAPSLIBS) \
	$(BUFRLIB) $(OTHERLIBS) -o $@
	rm -f *.o

$(POST): $(GSIPO) $(LAPSLIBS) $(HUMIDLIB)
	$(FC) $(FFLAGS) $(FREE) $(GSIPO) $(LAPSFC) $(LAPSLIBS) \
	$(BUFRLIB) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(STMAS) #$(STMAS_CG)
	$(RM) $(ANA) $(ANA_CG);
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS) $(BUFRLIB)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: $(ANA) # $(ANA_CG): compilation problem for cg_descent.f on toro.
	$(INSTALL) $(ANA) $(INSTALLROOT)/bin/$(ANA)
#	$(INSTALL) $(ANA_CG) $(INSTALLROOT)/bin/$(ANA_CG)

clean:
	rm -f $(ANA) $(ANA_CG) $(POST)  *.o *~ *# *.mod

distclean: clean

