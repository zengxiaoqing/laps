# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../..


include $(LAPSROOT)/src/include/makefile.inc










LAPSLIBS = $(LAPSROOT)/src/lib/cloud/libcloud.a  \
           $(LAPSROOT)/src/lib/liblaps.a  \
           $(LAPSROOT)/src/lib/mthermo/libmthermo.a  \
           $(LAPSROOT)/src/lib/goeslib/libgoeslib.a  \
           $(LAPSROOT)/src/lib/util/libutil.a  

DBLAPSLIBS = $(LAPSROOT)/src/lib/cloud/libdebugcloud.a  \
             $(LAPSROOT)/src/lib/libdebuglaps.a  \
             $(LAPSROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(LAPSROOT)/src/lib/goeslib/libdebuggoeslib.a  \
             $(LAPSROOT)/src/lib/util/libdebugutil.a  


OTHERLIBS = /usr/local/apps/netcdf-3.3.1/lib/libnetcdf.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = cloud.exe





SRC=array_plot.f \
    barnes_r5.f \
    insertpireps.f \
    insertradar.f \
    insertsao.f \
    insertsat.f \
    insertvis.f \
    laps_cloud.f \
    laps_cloud_sub.f \
    slice.f \
    spread.f 

all: exe

exe: $(EXE)

$(EXE): $(SRC)
	$(FC) $(FFLAGS) $(SRC) $(LAPSLIBS) $(OTHERLIBS) -o $@
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
