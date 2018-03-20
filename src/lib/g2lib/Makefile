# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

# CDEGRIBFLAGS are -DLINUX, -DAIX, etc; 
# DEGRIBFLAGS are -DUSE_PNG, -DUSE_JPEG2000 only; 
CFLAGS += $(CDEGRIBFLAGS) $(DEGRIBFLAGS)
ifndef USECPP
 FFLAGS += $(DEGRIBFLAGS)
 DBFLAGS += $(DEGRIBFLAGS)
endif

LIB      = libg2_4.a
DEBUGLIB = libdebugg2_4.a


SRC=	gridtemplates.f pdstemplates.f drstemplates.f gribmod.f realloc.f addfield.F  \
	addgrid.f addlocal.f getfield.f gb_info.f gf_getfld.f gf_free.f gf_unpack1.f  \
	gf_unpack2.f gf_unpack3.f gf_unpack4.f gf_unpack5.f gf_unpack6.f gf_unpack7.F \
	gettemplates.f getlocal.f getdim.f getpoly.f gribcreate.f gribend.f gribinfo.f \
	mkieee.f rdieee.f simpack.f simunpack.f cmplxpack.f compack.f misspack.f pack_gp.f \
	reduce.f comunpack.f specpack.f specunpack.f jpcpack.F jpcunpack.F enc_jpeg2000.c \
	dec_jpeg2000.c pngpack.F pngunpack.F enc_png.c dec_png.c gbytesc.f skgb.f ixgb2.f \
	getidx.f getg2i.f getg2ir.f getgb2s.f getgb2r.f getgb2l.f getgb2.f getgb2p.f getgb2rp.f \
	putgb2.f g2grids.f gdt2gds.f params.f mova2i.c

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))


all: lib

lib: $(LIB)

$(LIB): $(ARCHIVE) 
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS) $(CDEGRIBFLAGS) $(DEGRIBFLAGS)' \
	FFLAGS='$(DBFLAGS)'

install: lib

ifdef USECPP
%.o :  %.F
	$(CPP) $(CPPFLAGS) $(DEGRIBFLAGS) $< 
	mv $*.i $*.f
	$(FC) -c $(FFLAGS) $*.f
	$(RM) $*.f
endif

clean:
	$(RM) $(LIB) $(DEBUGLIB) *.o *.mod *~ *#

distclean: clean

