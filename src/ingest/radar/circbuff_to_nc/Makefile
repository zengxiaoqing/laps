# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(LAPSROOT)/src/lib/radar/a2io/liba2rt.a  \
           $(LAPSROOT)/src/lib/radar/a2io/liba2tp.a  \
           $(LAPSROOT)/src/lib/radar/nexrad_nssl/libnexrad.a \
           $(LAPSROOT)/src/lib/radar/nexrad_nssl/libnssl.a \
           $(LAPSROOT)/src/lib/radar/remap_ftn/libremap_ftn.a \
           $(LAPSROOT)/src/lib/liblaps.a  

DBLAPSLIBS = $(LAPSROOT)/src/lib/radar/a2io/liba2rt.a  \
             $(LAPSROOT)/src/lib/radar/a2io/liba2tp.a  \
             $(LAPSROOT)/src/lib/radar/nexrad_nssl/libdebugnexrad.a \
             $(LAPSROOT)/src/lib/radar/nexrad_nssl/libdebugnssl.a \
             $(LAPSROOT)/src/lib/radar/remap_ftn/libdebugremap_ftn.a  \
             $(LAPSROOT)/src/lib/libdebuglaps.a  

.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = circbuff_to_nc.exe

SRC=circbuffer_to_netcdf.c 

all: exe

exe: $(EXE)

$(EXE): $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LAPSLIBS) $(OTHERLIBS) -o $@
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
