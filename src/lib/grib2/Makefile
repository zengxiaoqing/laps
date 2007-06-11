SRCROOT=../..

include $(SRCROOT)/include/makefile.inc

LIB = libgrib2.a
DEBUGLIB = libdebuggrib2.a
SRC  = boustro_flt.f \
	boustro_int.f \
	check_flt.f \
	check_int.f \
	earth.f \
	flr.f \
	flt_map.f \
	fmkieee.f \
	int_map.f \
	length.f \
	pack_gp.f \
	pack_opt.f \
	pk_azimuth.f \
	pkbg.f \
	pk_bmap.f \
	pk_c7.f \
	pk_cmplx.f \
	pk_cylinder.f \
	pk_endian.f \
	pk_equator.f \
	pk_grib2.f \
	pk_lambert.f \
	pk_mercator.f \
	pk_missp.f \
	pk_nomiss.f \
	pk_orthographic.f \
	pk_polster.f \
	pk_s7.f \
	pk_sect0.f \
	pk_sect1.f \
	pk_sect2.f \
	pk_sect3.f \
	pk_sect4.f \
	pk_sect5.f \
	pk_sect6.f \
	pk_sect7.f \
	pk_sect8.f \
	pk_smple.f \
	pk_temp40.f \
	pk_temp41.f \
	pk_temp420.f \
	pk_temp42.f \
	pk_temp430.f \
	pk_temp48.f \
	pk_temp49.f \
	pk_trace.f \
	prep_flt.f \
	prep_int.f \
	prep_noval.f \
	prepr.f \
	prep_sect2_int.f \
	prep_sect2_real.f \
	reduce.f \
	timpr.f \
	io_utils.c \
	pk_swap.c


ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))

all:	lib

lib:	$(LIB)

$(LIB):	$(ARCHIVE)
	$(RANLIB) $@

debug:
	$(MAKE) all LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: 	lib

clean:
	$(RM) $(LIB) $(DEBUGLIB) *.o *~ *#

distclean: 	clean

ifdef USECPP
%.o : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $<
	mv $*.i $*.f
	$(FC) $(FFLAGS) -c $*.f
endif

