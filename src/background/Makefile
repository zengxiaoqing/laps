# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

#GFORTRAN LAPSLIBS = $(SRCROOT)/src/lib/bgdata/libbgdata.a \
#GFORTRAN           $(SRCROOT)/src/lib/util/libutil.a \
#GFORTRAN           $(SRCROOT)/src/lib/liblaps.a  \
#GFORTRAN           $(SRCROOT)/src/lib/mthermo/libmthermo.a \
#GFORTRAN           $(SRCROOT)/src/lib/modules/libmodules.a \
#GFORTRAN           $(SRCROOT)/src/lib/degrib/libdegrib.a \
#GFORTRAN           $(SRCROOT)/src/lib/w3lib/libw3.a \
#GFORTRAN           $(SRCROOT)/src/lib/g2lib/libg2_4.a 
#GFORTRAN modifications begin
LAPSLIBS = $(SRCROOT)/src/lib/bgdata/libbgdata.a \
           $(SRCROOT)/src/lib/util/libutil.a \
           $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a \
           $(SRCROOT)/src/lib/modules/libmodules.a \
           $(SRCROOT)/src/lib/degrib/libdegrib.a \
           $(SRCROOT)/src/lib/modules/libmodules.a \
           $(SRCROOT)/src/lib/w3lib/libw3.a \
           $(SRCROOT)/src/lib/g2lib/libg2_4.a 
#GFORTRAN modifications end

DBLAPSLIBS = $(SRCROOT)/src/lib/bgdata/libdebugbgdata.a \
             $(SRCROOT)/src/lib/util/libdebugutil.a \
             $(SRCROOT)/src/lib/libdebuglaps.a   \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/degrib/libdebugdegrib.a \
             $(SRCROOT)/src/lib/w3lib/libdebugw3.a \
             $(SRCROOT)/src/lib/g2lib/libdebugg2_4.a 

.SUFFIXES: .o .c .f .a .F .i .exe

EXE = lga.exe

SRC=lga.f

LIB=liblga.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(SRC)
	$(FC) $(SRC) $(FFLAGS) $(LIB) $(LAPSLIBS) \
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



