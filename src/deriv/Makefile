# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(LAPSROOT)/src/lib/cloud/libcloud.a  \
           $(LAPSROOT)/src/lib/liblaps.a  \
           $(LAPSROOT)/src/lib/mthermo/libmthermo.a  \
           $(LAPSROOT)/src/lib/util/libutil.a  

DBLAPSLIBS = $(LAPSROOT)/src/lib/cloud/libdebugcloud.a  \
             $(LAPSROOT)/src/lib/libdebuglaps.a  \
             $(LAPSROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(LAPSROOT)/src/lib/util/libdebugutil.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = deriv.exe

SRC=laps_deriv.f  laps_deriv_sub.f   cloud_deriv_subs.f  put_derived_wind_prods.f 

all: exe

exe: $(EXE)

$(EXE): $(SRC)
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
