# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../..


include $(LAPSROOT)/src/include/makefile.inc

KRCFLAGS=-g $(INC)

.SUFFIXES:
.SUFFIXES: .o .c .f .a

LIB      = libsatellite.a
DEBUGLIB = libdebugsatellite.a




SRC=bcd_to_int.f \
    bcdtoint.f \
    byteswp2.f \
    byteswp4.f \
    c_time_fnc.f \
    cnvtglfl.f \
    cnvtint.f \
    cnvtvxfl.f \
    cnvvxdb.f \
    cmpsatres.f \
    convert_to_real.f \
    convert_to_double.f \
    expand_domain.f \
    gen_lut_gvar.f \
    gen_lut_lambert.f \
    genlvdlut_sub.f \
    genlvdlut_lvd.f \
    get_attr_satellite.f \
    get_gwc_oa.f \
    getsat_boundary.f \
    lvd_file_specifier.f \
    rdgvarimg_cdf_header.f \
    read_gwc_header_sub.f \
    read_orb_att.f \
    rewrite_satellite_lvd_nl.f \
    satsublatlon.f \
    stats.f \
    update_gvarimg_parms.f \
    write_orb_att.f \
    write_att_c.c

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
