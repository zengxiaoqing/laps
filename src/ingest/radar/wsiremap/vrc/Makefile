# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../../../..


include $(LAPSROOT)/src/include/makefile.inc










LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a   \
           $(LAPSROOT)/src/lib/util/libutil.a 


DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a   \
             $(LAPSROOT)/src/lib/util/libdebugutil.a 

OTHERLIBS = /usr/local/apps/netcdf-3.3.1/lib/libnetcdf.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .x .exe

EXE = vrc_driver.x




SRC= nowrad2laps.f \
     read_table.f \
     read_wsi_cdf_wsi.f \
     procnowradZ.f \
     read_wsi_cdf_wfo.f \
     readvrclut.f

MAIN=vrc_driver.f

LIB=libvrc.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(FC) $(MAIN) $(FFLAGS) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

#$(EXE.F): $(MAIN) $(ARCHIVE)
#	rm -f $*.f
#	$(CPP) $(CPPFLAGS) $*.F 
#	mv $*.i $*.f
#	$(FC) $(FFLAGS) $*.f $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
#	rm -f *.o $*.f

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

clean:
	rm -f $(EXE) $(LIB) *.o *~ *# *.i


distclean: clean
install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

%.f : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f

%.o : %.F


