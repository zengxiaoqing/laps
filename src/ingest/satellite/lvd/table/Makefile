# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../../../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(LAPSROOT)/src/lib/satellite/libsatellite.a \
           $(LAPSROOT)/src/lib/liblaps.a  \
           $(LAPSROOT)/src/lib/goesinav/libgoesinav.a  \
           $(LAPSROOT)/src/lib/nav/libnav.a  

DBLAPSLIBS = $(LAPSROOT)/src/lib/satellite/libdebugsatellite.a \
             $(LAPSROOT)/src/lib/libdebuglaps.a  \
             $(LAPSROOT)/src/lib/goesinav/libdebuggoesinav.a  \
             $(LAPSROOT)/src/lib/nav/libdebugnav.a  

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = genlvdlut.exe

SRC = compute_rirj.f \
      gen_lut_ascii.f \
      read_ascii_satdat.f \
      get_ascii_dimensions.f

MAIN=genlvdlut.f

LIB=libgenlvdlut.a

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
