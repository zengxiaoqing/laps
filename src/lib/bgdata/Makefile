# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..


include $(SRCROOT)/src/include/makefile.inc




.SUFFIXES:
.SUFFIXES: .o .c .f .a

LIB      = libbgdata.a
DEBUGLIB = libdebugbgdata.a



SRC= bkgd_stn_interp.f   bkgd_util.f            es_table.f             \
     get_acceptable_files.f get_bkgd_mdl_info.f get_attr_sbn.f         \
     getdims_lapsprd.f   lga_driver.f                                  \
     lga_interp_subs.f   lprep_eta_conusc.f     lprep_ruc2_pub.f       \
     readavnpublic.f     readbgdata.f           readeta48_pub.f        \
     readconusnws.f      readruc2_pub.f         read_fuafsf.f          \
     read_netcdf_field.f read_sbn_grids.f       readdgprep.f           \
     readnogaps.f        readruc60native.f      rdmodelindex.f         \
     sfcbkgd.f           sfcbkgd_sfc.f          tcbogus.f              \
     writelgalgb.f


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
