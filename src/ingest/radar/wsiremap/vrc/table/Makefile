# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../../../../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a   \
           $(LAPSROOT)/src/lib/nav/libnav.a \
           $(LAPSROOT)/src/lib/goesinav/libgoesinav.a 

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a   \
             $(LAPSROOT)/src/lib/nav/libdebugnav.a \
             $(LAPSROOT)/src/lib/goesinav/libdebuggoesinav.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .x .exe

EXE = genvrclut.exe

SRC = gen_llij_lut_wfo.f \
      gen_llij_lut_wsi.f \
      gen_llij_lut_ll.f \
      get_radar_bounds.f 

MAIN=genvrclut.f

LIB=libvrclut.a

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

