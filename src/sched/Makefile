# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a  

.SUFFIXES:
.SUFFIXES: .o .c .f .a .exe

EXE = systime.exe 
SRC = systime.f 

all: exe

exe: $(EXE)

%.exe: %.f
	$(FC) $(FFLAGS) $*.f $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	@for file in $(EXE);\
	  do \
	  $(INSTALL) $$file $(INSTALLROOT)/bin/$$file; \
	done

clean:
	rm -f $(EXE)  *.o *~ *#

distclean: clean
