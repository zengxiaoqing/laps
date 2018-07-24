# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/radar/wsi_ingest/libwsi.a \
           $(SRCROOT)/src/lib/liblaps.a \
           $(SRCROOT)/src/lib/satellite/libsatellite.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a \
           $(SRCROOT)/src/lib/util/libutil.a \
           $(SRCROOT)/src/lib/nav/libnav.a

DBLAPSLIBS = $(SRCROOT)/src/lib/radar/wsi_ingest/libdebugwsi.a \
             $(SRCROOT)/src/lib/libdebuglaps.a   \
             $(SRCROOT)/src/lib/satellite/libdebugsatellite.a   \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/util/libdebugutil.a \
             $(SRCROOT)/src/lib/nav/libdebugnav.a

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .x .exe

EXE = vrc_driver.x

SRC= nowradwsi2laps.f \
     procnowradZ.f \
     aoml_radar.f \
     read_pass_head.f \
     put_vrc.f

MAIN=vrc_driver.f

LIB=libvrc.a

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


