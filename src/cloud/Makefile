# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/cloud/libcloud.a \
           $(SRCROOT)/src/lib/liblaps.a \
           $(SRCROOT)/src/lib/modules/libmodules.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a  

DBLAPSLIBS = $(SRCROOT)/src/lib/cloud/libdebugcloud.a \
             $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a  
             

.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = cloud.exe

SRC=array_plot.f \
    barnes_r5.f \
    compare_analysis_to_rad.f \
    get_fg_cloud.f \
    get_istat_39.f \
    get_sat_data.f \
    getvis.f \
    insertradar.f \
    insertsat.f \
    insertvis.f \
    insert_co2ctp.f \
    laps_cloud.f \
    laps_cloud_sub.f \
    slice.f

all: exe

exe: $(EXE)

$(EXE): $(SRC) $(LAPSLIBS)
	$(FC) $(FFLAGS) $(FIXED) $(SRC) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE)  *.o *~ *#


distclean: clean
