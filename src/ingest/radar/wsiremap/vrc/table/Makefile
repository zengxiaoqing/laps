# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/radar/wsi_ingest/libwsi.a \
           $(SRCROOT)/src/lib/nav/libnav.a \
           $(SRCROOT)/src/lib/liblaps.a   



DBLAPSLIBS = $(SRCROOT)/src/lib/radar/wsi_ingest/libdebugwsi.a \
             $(SRCROOT)/src/lib/nav/libdebugnav.a  \
             $(SRCROOT)/src/lib/libdebuglaps.a   

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .x .exe

EXE = genvrclut.exe

SRC=genvrclut.f

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(SRC)
	$(FC) $(FFLAGS) $(SRC) $(LAPSLIBS) $(OTHERLIBS) -o $@
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

