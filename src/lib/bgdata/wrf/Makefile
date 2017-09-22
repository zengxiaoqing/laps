# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../..


include $(SRCROOT)/src/include/makefile.inc




.SUFFIXES:
.SUFFIXES: .o .f90 .a

LIB      = ../libbgdata.a
DEBUGLIB = ../libdebugbgdata.a



SRC=lga_driver_wrfarw.f90 \
    get_acceptable_wrf.f90 \
    module_wrf_lga.f90 \
    module_constants.f90 \
    met_functions.f90 \
    smooth2.f90

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))

.f90.o:
	$(RM) $@
	$(FC) -c $(FFLAGS) $(FREE) $(INC) $*.f90

all: lib

lib: $(LIB)

$(LIB): $(ARCHIVE) 
	$(RANLIB) $@

lga_driver_wrfarw.o:	module_wrf_lga.o

module_wrf_lga.o:	module_constants.o

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: lib

cleanup:
	rm -f *.o

clean:
	rm -f $(LIB) $(DEBUGLIB) *.o *~ *#


distclean: clean
