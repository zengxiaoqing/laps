# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../../../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a   \
           $(LAPSROOT)/src/lib/radar/wsi_ingest/libwsi.a 

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a   \
             $(LAPSROOT)/src/lib/radar/wsi_ingest/libdebugwsi.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .x .exe

EXE = ln3_driver.exe

SRC = nowradwsi2laps.f \
     wait_for_radar.f \
     checkradardata.f

MAIN=ln3_driver.f

LIB=libln3.a

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


