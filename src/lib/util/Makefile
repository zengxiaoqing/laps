# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..


include $(SRCROOT)/src/include/makefile.inc

KRCFLAGS=-g $(INC)









.SUFFIXES:
.SUFFIXES: .o .c .f .a

LIB      = libutil.a
DEBUGLIB = libdebugutil.a




SRC=bilinear_laps.f   filter_2dx.f      get_static_field_interp.f\
    helicity.f        interp2sfc.f\
    put_wind_2d.f     trilinear_laps.f  cpt_liw.f         laps_slwc_revb.f\
    stats.f           steer_grid.f      windcnvrt.f       degfrom.f         meanwind.f\
    storm_cent_rt.f   xytomet.f         windfields.f



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
