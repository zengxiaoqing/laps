# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  


DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a 


.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = wind.exe

SRC= main_sub.f \
     get_wind_obs.f \
     weights.f \
     comparisons.f \
     comp_vr_prof.f \
     comp_laps_vr.f \
     comp_laps1_laps2.f \
     comp_laps_maps.f \
     dualwind_noZ.f \
     dualwind_yesZ.f \
     vertwind.f \
     rdpoint.f \
     rdsfc.f \
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





















