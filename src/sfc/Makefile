# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../..


include $(LAPSROOT)/src/include/makefile.inc










LAPSLIBS = $(LAPSROOT)/src/lib/temp/libtemp.a  \
           $(LAPSROOT)/src/lib/liblaps.a  \
           $(LAPSROOT)/src/lib/mthermo/libmthermo.a

DBLAPSLIBS = $(LAPSROOT)/src/lib/temp/libdebugtemp.a  \
             $(LAPSROOT)/src/lib/libdebuglaps.a  \
             $(LAPSROOT)/src/lib/mthermo/libdebugmthermo.a  



OTHERLIBS = /usr/local/apps/netcdf-3.3.1/lib/libnetcdf.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = laps_sfc.x





SRC   = mdatlaps.f \
	lapsvanl.f \
	lapsvanl_sub.f \
	laps_routines.f \
	qc_data.f 

MAIN=laps_sfc.f 



LIB=libsfc.a

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















