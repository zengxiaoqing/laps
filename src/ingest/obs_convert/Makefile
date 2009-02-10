# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

#=============================================================
#  LAPS environment:
#=============================================================

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

RM = rm -f
FLAG = $(FFLAGS)

LAPSLIBS = $(SRCROOT)/src/lib/temp/libtemp.a \
           $(SRCROOT)/src/lib/wind/libwind.a  \
           $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a \
           $(SRCROOT)/src/var/bufr/libbufr.a

DBLAPSLIBS = $(SRCROOT)/src/lib/temp/libdebugtemp.a \
             $(SRCROOT)/src/lib/wind/libdebugwind.a  \
             $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a \
             $(SRCROOT)/src/var/bufr/libdebugbufr.a

#LAPSFC = $(SRCROOT)/src/sfc/libsfc.a

#HUMIDLIB = $(SRCROOT)/src/humid/liblq3.a

.SUFFIXES:
.SUFFIXES: .o .c .f90 .a .f
.f90.o:;    $(FC) $(FLAG) $(FREE) -c $*.f90
.f.o:;      $(FC) $(FLAG) -c $*.f
.c.o:;	    $(CC) $(CPPFLAGS) -c $*.c

#=============================================================
# 3DVAR:
#=============================================================

CONV = obs_convert.exe

FILES= LAPS_Params.f90 \
       Main_Convrt.f90 \
       LAPS_Dividr.f90 \
       LAPS_Config.f90 \
       LAPS_Ingest.f90 \
       Read_RASS_Data.f90 \
       LAPS_Systim.f90 \
       BUFR_Convrt.f90 \
       WRFD_Convrt.f90

all: $(CONV) 

$(CONV): $(FILES) $(LAPSLIBS)
	$(FC) $(FFLAGS) $(FREE) $(FILES) $(LAPSFC) $(LAPSLIBS) \
	$(OTHERLIBS) -o $@
	rm -f *.o

debug: $(FILES)
	$(RM) $(CONV);
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: $(CONV) 
	$(INSTALL) $(CONV) $(INSTALLROOT)/bin/$(CONV)

clean:
	rm -f $(CONV)  *.o *~ *# *.mod

distclean: clean

