# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../..


include $(LAPSROOT)/src/include/makefile.inc




.SUFFIXES:
.SUFFIXES: .o .c .f .a

LIB      = libbgdata.a
DEBUGLIB = libdebugbgdata.a



SRC= get_bkgnd_dims.f   readeta48_pub.f readconusnws.f \
     readruc2_pub.f get_acceptable_files.f read_netcdf_field.f\
     read_conus211.f readdgprep.f readnogaps.f readruc60native.f \
     es_table.f rdmodelindex.f sfcbkgd.f unpkfgge_a.c



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
