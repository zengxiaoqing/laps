# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/radar/a2io/liba2rt.a  \
           $(SRCROOT)/src/lib/radar/a2io/liba2tp.a  \
           $(SRCROOT)/src/lib/radar/nexrad_nssl/libnexrad.a \
           $(SRCROOT)/src/lib/radar/nexrad_nssl/libnssl.a \
           $(SRCROOT)/src/lib/radar/remap_ftn/libremap_ftn.a \
           $(SRCROOT)/src/lib/liblaps.a  

DBLAPSLIBS = $(SRCROOT)/src/lib/radar/a2io/liba2rt.a  \
             $(SRCROOT)/src/lib/radar/a2io/liba2tp.a  \
             $(SRCROOT)/src/lib/radar/nexrad_nssl/libdebugnexrad.a \
             $(SRCROOT)/src/lib/radar/nexrad_nssl/libdebugnssl.a \
             $(SRCROOT)/src/lib/radar/remap_ftn/libdebugremap_ftn.a  \
             $(SRCROOT)/src/lib/libdebuglaps.a  

CMAIN_FSUB=-lxlf90
OTHERLIBS2=$(OTHERLIBS) $(CMAIN_FSUB)

.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = circbuff_to_nc.exe

SRC=circbuffer_to_netcdf.c 

all: exe

exe: $(EXE)

$(EXE): $(SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(SRCROOT)/src/include $(SRC) $(LAPSLIBS) $(OTHERLIBS2) -o $@
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
