# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = 

DBLAPSLIBS = 

.SUFFIXES:
.SUFFIXES: .o .c .f .a .exe

EXE = xfer_laps.exe \
      xfer_model.exe \
      fname2unixtime.exe \
      in_domain.exe

SRC = fill_laps.c \
      fill_bigfile_cmn.c \
      fill_model.c 
      
LIB = lapscdf.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

%.exe: %.c $(LIB)
	$(CC) $(CPPFLAGS) $(CFLAGS) $*.c $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

%.exe: %.f $(LIB)
	$(FC) $(FFLAGS) $*.f $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
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
	rm -f $(EXE)  *.o *~ *\# $(LIB)

distclean: clean
