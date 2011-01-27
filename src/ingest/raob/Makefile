# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a  
  
DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a


.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = ingest_sounding.exe

SRC=ingest_sounding.f \
    ingest_raob.f \
    ingest_raob_sub.f \
    ingest_raob_sub_af.f \
    ingest_raob_sub_cwb.f \
    ingest_drpsnd.f \
    ingest_drpsnd_sub_cwb.f \
    avapsread_sub.f \
    get_drpsnd_data.f \
    ingest_goessnd.f \
    process_goes_snd.f \
    get_poes_data.f \
    ingest_satsnd.f \
    ingest_satsnd_af.f \
    tower_driver.f \
    get_rtamps_data.f \
    read_rtamps_data.f \
    get_radiometer_data.f

all: exe

exe: $(EXE)

$(EXE): $(SRC) $(LAPSLIBS) 
	$(FC) $(LITTLEFLAGS) $(FIXED) $(SRC) $(LAPSLIBS) $(OTHERLIBS) -o $@
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
