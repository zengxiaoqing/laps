# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a   \
           $(SRCROOT)/src/lib/util/libutil.a \
           $(SRCROOT)/src/lib/nav/libnav.a

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a   \
             $(SRCROOT)/src/lib/util/libdebugutil.a \
             $(SRCROOT)/src/lib/nav/libdebugnav.a

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .x .exe

EXE = ingest_rrv.exe

SRC =  genv01lut_sub.f \
       gen_llij_lut_v01.f \
       readradar_v01_cdf.f \
       rrvdat2laps_v01.f \
       buildcomment.f \
       config_rrv_common.f

MAIN=ingest_rrv.f $(SRCROOT)/src/lib/block_data.f

LIB=libgenrrv.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(FC) $(MAIN) $(FFLAGS) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

clean:
	rm -f $(EXE) $(LIB) *.o *~ *# *.i


distclean: clean

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

%.f : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f

%.o : %.F
