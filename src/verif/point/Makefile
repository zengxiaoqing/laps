# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a

RM=rm -f

.SUFFIXES: .o .f .a .exe

LIB = libverif.a
DEBUGLIB = libdebugverif.a
ARCHIVE = ar rv
.f.o:
		$(RM) $@
		$(FC) -c $(FFLAGS) $(INC) $*.f
		$(ARCHIVE) $(LIB) $(<:.f=.o)


EXE1= gen_verif_raob.exe
EXE2= gen_verif_prof.exe

SRC = bkgd_wgi.f \
    gen_verif_prof_sub.f \
    get_prof_pairs.f \
    gen_verif_raob_sub.f \
    get_raob_data_a.f \
    get_raob_pairs.f \
    raob_subs.f

MAIN1 = gen_verif_raob.f
MAIN2 = gen_verif_prof.f

OBJS=$(SRC:.f=.o)

all:	$(EXE1) $(EXE2)

$(EXE1):$(OBJS)
	$(FC) -o $(EXE1) $(FFLAGS) $(MAIN1) -I$(SRCROOT)/src/include $(LIB) $(LAPSLIBS) $(OTHERLIBS) $(MTHERMOLIB) $(UTILLIB);
	rm -f *.o

$(EXE2):$(OBJS)
	$(FC) -o $(EXE2) $(FFLAGS) $(MAIN2) -I$(SRCROOT)/src/include $(LIB) $(LAPSLIBS) $(OTHERLIBS) $(MTHERMOLIB) $(UTILLIB);
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE1) $(RM) $(EXE2) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' LIB='$(DEBUGLIB)' MTHERMOLIB='$(DBMTHERMOLIB)' UTILLIB='$(DBUTILLIB)' \
	FFLAGS='$(DBFLAGS)'

install: all
	$(INSTALL) $(EXE1) $(INSTALLROOT)/bin/$(EXE1); \
	$(INSTALL) $(EXE2) $(INSTALLROOT)/bin/$(EXE2);

clean:
	rm -f $(EXE1) $(EXE2)  *.o *.a *~ *\#

distclean: clean
