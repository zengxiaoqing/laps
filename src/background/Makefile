# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(BGDATALIB) $(UTILLIB) $(LAPSLIB) $(MTHERMOLIB) 

DBLAPSLIBS = $(DBBGDATALIB) $(DBUTILLIB) $(DBLAPSLIB) $(DBMTHERMOLIB) 


.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = lga.exe

SRC=lga_interp.f  lga_util.f           \
    lprep_eta_conusc.f  lprep_ruc2_pub.f writelgalgb.f

MAIN=lga.f

LIB=liblga.a

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



