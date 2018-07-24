# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/bgdata/libbgdata.a \
           $(SRCROOT)/src/lib/liblaps.a \
           $(SRCROOT)/src/lib/satellite/libsatellite.a  \
           $(SRCROOT)/src/lib/util/libutil.a \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a \
           $(SRCROOT)/src/lib/degrib/libdegrib.a \
           $(SRCROOT)/src/lib/modules/libmodules.a \
           $(SRCROOT)/src/lib/w3lib/libw3.a \
           $(SRCROOT)/src/lib/g2lib/libg2_4.a


DBLAPSLIBS = $(SRCROOT)/src/lib/bgdata/libdebugbgdata.a \
             $(SRCROOT)/src/lib/libdebuglaps.a \
             $(SRCROOT)/src/lib/satellite/libdebugsatellite.a   \
             $(SRCROOT)/src/lib/util/libdebugutil.a \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a \
             $(SRCROOT)/src/lib/degrib/libdebugdegrib.a \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/w3lib/libdebugw3.a \
             $(SRCROOT)/src/lib/g2lib/libdebugg2_4.a

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe .x

EXE = qbalpe.exe

SRC= writeballaps.f qbal_utils.f profile.f var3dsub.f

MAIN=qbalpe.f 

LIB=libqbalpe.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(FC) $(MAIN) $(FFLAGS) $(LIB) $(LAPSLIBS) \
	$(OTHERLIBS) $(DEGRIBLIBS) -o $@
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
