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

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a 

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a

LAPSFC = $(SRCROOT)/src/sfc/libsfc.a

BUFRLIB = $(SRCROOT)/src/var/bufr/libbufr.a

WINDLIB = $(SRCROOT)/src/wind/libwind.a


.SUFFIXES:
.SUFFIXES: .o .c .f90 .a .f
.f90.o:;    $(FC) $(FLAG) $(FREE) -c $*.f90
.f.o:;      $(FC) $(FLAG) -c $*.f
.c.o:;	    $(CC) $(CPPFLAGS) -c $*.c

#=============================================================
# 3DVAR:
#=============================================================

PREP = GSI_Prep.exe

GSI= GSI_Prep.f90 \
       LAPS_Parm.f90 \
       LAPS_Conf.f90 \
       LAPS_Bkgd.f90 \
       LAPS_Obsv.f90 \
         $(WINDLIB) \
       LAPS_Memo.f90 \
       GSI_BkObs.f90 \
       LAPS2Mass.f90 \
       wrfbkgout.f90

all: $(PREP) # $(POST)

$(PREP): $(GSI) $(LAPSLIBS) $(BUFRLIB) $(WINDLIB)
	$(FC) $(FFLAGS) $(FREE) $(GSI) $(LAPSFC) $(LAPSLIBS) \
	$(BUFRLIB) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(GSI)
	$(RM) $(PREP) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS) $(BUFRLIB)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: $(PREP)
	$(INSTALL) $(PREP) $(INSTALLROOT)/bin/$(PREP)

clean:
	rm -f $(PREP)  *.o *~ *# *.mod

distclean: clean

