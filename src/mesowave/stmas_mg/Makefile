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

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a 

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a

LAPSFC = $(SRCROOT)/src/sfc/libsfc.a


.SUFFIXES:
.SUFFIXES: .o .c .f90 .a .f
.f90.o:;    $(FC) $(FLAG) $(FREE) -c $*.f90
.f.o:;      $(FC) $(FLAG) -c $*.f

#=============================================================
# STMAS:
#=============================================================

EXE = stmas_mg.x

SRC=Definition.f90 \
    LAPSDatSrc.f90 \
    MemoryMngr.f90 \
    PrePostPrc.f90 \
    STMASAnalz.f90 \
    STMAS_MAIN.f90
    

all: exe

exe: $(EXE)

$(EXE): $(SRC) $(LAPSLIBS)
	$(FC) $(FFLAGS) $(FREE) $(SRC) $(LAPSFC) $(LAPSLIBS) \
	$(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	echo "Hello " $(RM)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE)  *.o *~ *# *.mod
#	rm -f LAPSDatSrc.f90 Definition.f90 STMASAnalz.f90

distclean: clean

# Source code generator:
#LAPSDatSrc.f90: LAPSInfo.f90 LAPSBKGD.f90 LAPSOBSV.f90 LAPS_QCs.f90 \
#	LAPSConf.f90 LAPSIntp.f90
#	$(CPP) $(CPPFLAGS) LAPSDatSrc.F90
#	mv LAPSDatSrc.i LAPSDatSrc.f90
#STMASAnalz.f90: STMASAna.f90
#	$(CPP) $(CPPFLAGS) STMASAnalz.F90 
#	mv STMASAnalz.i STMASAnalz.f90
