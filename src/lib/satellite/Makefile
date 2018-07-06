# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..


include $(SRCROOT)/src/include/makefile.inc

KRCFLAGS=-g $(INC)

LIB      = libsatellite.a
DEBUGLIB = libdebugsatellite.a

SRC=bcdtoint.F \
    bcd_to_int.f \
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
    fc01_conv_pts_el_to_met_1.f \
    gen_lut_fx.f \
    gen_lut_gvar.f \
    gen_lut_lambert.f \
    gen_lut_ll.f \
    gen_lut_mercator.f \
    genlutps.f \
    genlvdlut_sub.f \
    genlvdlut_lvd.f \
    geos_transform_pix.c \
    get_attr_gnp.f \
    get_attr_satellite.f \
    getsat_boundary.f \
    lvd_file_specifier.f \
    rdcdfhead.f \
    rdgvarimg_cdf_header.f \
    read_gwc_header_sub.f \
    read_orb_att.f \
    readnet.f \
    satsublatlon.f \
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

ifdef USECPP
%.o : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f
	$(FC) $(FFLAGS) -c $*.f
endif
