# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/modules/libmodules.a  \
           $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a

DBLAPSLIBS = $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a

.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = deriv.exe

SRC=laps_deriv.f  \
    laps_deriv_sub.f \
    cloud_deriv_subs.f \
    fireweatherindex.f \
    hainesindex.f \
    fire_fields.f \
    put_derived_wind_prods.f \
    rfill_evap.f \
    schultz.f

all: exe

exe: $(EXE)

$(EXE): $(SRC) $(LAPSLIBS)
	$(FC) $(FFLAGS) $(FIXED) $(SRC) $(LAPSLIBS) $(OTHERLIBS) -o $@
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
