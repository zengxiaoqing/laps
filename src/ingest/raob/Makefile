# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a  

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a  

.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = ingest_sounding.exe

SRC=ingest_sounding.f \
    ingest_raob.f \
    ingest_raob_sub.f \
    ingest_raob_sub_af.f \
    ingest_satsnd.f \
    ingest_satsnd_af.f

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
