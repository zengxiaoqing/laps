# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a \
           $(LAPSROOT)/src/lib/mthermo/libmthermo.a  
           

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a \
             $(LAPSROOT)/src/lib/mthermo/libdebugmthermo.a              

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = obs_driver.x

SRC   = obs_driver_sub_af.f \
        get_metar_obs.f \
        get_metar_obs_af.f \
        read_metar.f \
        read_metar_cwb.f \
        read_metar_af.f \
        read_synop_cwb.f \
	get_local_obs.f \
	read_local.f \
	get_local_cwb.f \
	read_local_cwb.f \
        get_buoy_obs.f \
	read_buoy.f \
	read_buoy_cwb.f \
	read_ship_cwb.f \
        get_gps_obs.f \
        read_gps.f

MAIN=obs_driver.f 

LIB=libsao.a

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















