# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/temp/libtemp.a  \
           $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a \
           $(SRCROOT)/src/lib/util/libutil.a \
           $(SRCROOT)/src/lib/modules/libmodules.a  

DBLAPSLIBS = $(SRCROOT)/src/lib/temp/libdebugtemp.a  \
             $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a \
             $(SRCROOT)/src/lib/util/libdebugutil.a \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = laps_sfc.x

SRC   = laps_sfc_sub.f \
        mdatlaps.f \
	lapsvanl.f \
	lapsvanl_sub.f \
        barnes_multivariate_sfc.f \
	qc_data.f \
        visibility.f 

MAIN=laps_sfc.f 

LIB=libsfc.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(FC) $(MAIN) $(FFLAGS) $(FIXED) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
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















