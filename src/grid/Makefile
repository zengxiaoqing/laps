# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a \
           $(SRCROOT)/src/lib/satellite/libsatellite.a  \
           $(SRCROOT)/src/lib/util/libutil.a \
           $(SRCROOT)/src/lib/modules/libmodules.a

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a \
             $(SRCROOT)/src/lib/satellite/libdebugsatellite.a   \
             $(SRCROOT)/src/lib/util/libdebugutil.a \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a

.SUFFIXES: .o .f .a .F .i .exe

EXE = gridgen_model.exe

SRC   = alt_all.f adjust_geog.f gen_egrid_latlon.f \
        get_tile_list.f gdswiz_comb.f gridgen_utils.f \
        proc_geodat.f read_wrfstatic.f put_laps_static.f \
        rd_static_attr_sub.f

MAIN=gridgen_model.f 

LIB=libgrid.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(FC) $(MAIN) $(FFLAGS) $(FMOD) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@ 
	rm -f *.o

debug: $(SRC) ;
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	  FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE) $(LIB) *.o *~ *# *.i *.mod


distclean: clean
