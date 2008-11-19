# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a \
           $(SRCROOT)/src/lib/modules/libmodules.a  

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a

.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = ingest_aircraft.exe

SRC=ingest_aircraft.f \
    ingest_pireps.f \
    ingest_pireps_sub.f \
    ingest_pireps_WFO_sub.f \
    skycover_to_frac.f \
    ingest_acars.f \
    ingest_acars_sub.f \
    ingest_acars_sub_af.f \
    ingest_wisdom.f \
    get_wisdom_data.f \
    write_aircraft_sub.f

all: exe

exe: $(EXE)

$(EXE): $(SRC) $(LAPSLIBS) 
	$(FC) $(FFLAGS) $(SRC) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE)  *.o *~ *#


distclean: clean
