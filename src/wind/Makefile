# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a  \
           $(LAPSROOT)/src/lib/mthermo/libmthermo.a  \
           $(LAPSROOT)/src/lib/lapack/liblapack.a  \
           $(LAPSROOT)/src/lib/blas/libblas.a  \
           $(LAPSROOT)/src/lib/util/libutil.a  \

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a  \
             $(LAPSROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(LAPSROOT)/src/lib/lapack/libdebuglapack.a  \
             $(LAPSROOT)/src/lib/blas/libdebugblas.a  \
             $(LAPSROOT)/src/lib/util/libdebugutil.a  \

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = wind.exe

SRC= main_sub.f \
     get_wind_obs.f \
     weights.f \
     comp_vr_prof.f \
     comp_laps_meso.f \
     comp_laps_sao.f \
     comp_maps_sao.f \
     comp_laps_vr.f \
     comp_laps_prof.f \
     comp_maps_prof.f \
     comp_laps_pirep.f \
     comp_maps_pirep.f \
     comp_laps1_laps2.f \
     comp_laps_maps.f \
     vertwind.f \
     rdpirep.f \
     rdsao.f \
     read_profiles.f \
     interp_prof.f \
     qcradar.f \
     get_fg_wind.f \
     windanal.f

MAIN=main.f

LIB=libwind.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(FC) $(MAIN) $(FFLAGS) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE) $(LIB) *.o *~ *# *.i


distclean: clean

%.f : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f

%.o : %.F





















