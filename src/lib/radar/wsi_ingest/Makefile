# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#

SRCROOT=../../../..

include $(SRCROOT)/src/include/makefile.inc

KRCFLAGS=-g $(INC)

.SUFFIXES:
.SUFFIXES: .o .c .f .a

LIB      = libwsi.a
DEBUGLIB = libdebugwsi.a

SRC=rd_wsi_nowrad_cdf.f   rd_wsi_nexrad_cdf.f \
    get_nl_parms.f        genllijlut.f

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))

all: lib

lib: $(LIB)

$(LIB): $(ARCHIVE) 
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: lib

clean:
	rm -f $(LIB) $(DEBUGLIB) *.o *~ *#

distclean: clean
