# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../..


include $(SRCROOT)/src/include/makefile.inc

KRCFLAGS=-g $(INC)









.SUFFIXES:
.SUFFIXES: .o .c .f .a

LIB      = libremap_ftn.a
DEBUGLIB = libdebugremap_ftn.a




SRC=fill_common.f\
    get_azimuths_deg.f\
    get_scandata.f\
    ld_ray.f\
    lut_gen.f\
    rcptoremote.f\
    readdata.f\
    ref_fill_horz.f\
    timconv.f

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
