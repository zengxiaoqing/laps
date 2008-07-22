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

EXE = ingest_pro.exe

SRC=ingest_pro.f \
    ingest_sub_pro.f \
    ingest_sub_blppro.f \
    ingest_sub_rsapro.f \
    get_sodar_data.f \
    ingest_vad.f \
    ingest_vad_sub.f \
    write_pro.f

all: exe

exe: $(EXE)

$(EXE): $(SRC) $(LAPSLIBS) 
	$(FC) $(FFLAGS) $(SRC) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

clean:
	rm -f $(EXE)  *.o *~ *#


distclean: clean

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)









