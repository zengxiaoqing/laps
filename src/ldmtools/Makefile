# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = 

DBLAPSLIBS = 

LDMBASE=/usr/local/ldm

LOSLIBS=
OTHERLIBS.=-L$(LDMBASE)/lib -lldm $(LOSLIBS)
CFLAGS.=-I$(LDMBASE)/include

.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = fsl_pqinsert

SRC=fsl_pqinsert.c

all: exe

exe: $(EXE)

$(EXE): $(SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SRC) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE)  *.o *~ *#


distclean: clean
