# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a  

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a  

.SUFFIXES:
.SUFFIXES: .o .c .f .a .exe

EXE = xfer_laps.exe \
      xfer_model.exe \
      fname2unixtime.exe

SRC = xfer_laps.c \
      fill_laps.c \
      fill_bigfile_cmn.c \
      xfer_model.c \
      fill_model.c \
      fname2unixtime.c  
LIB = lapscdf.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

%.exe: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $*.c $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

%.exe: %.f
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
	rm -f $(EXE)  *.o *~ *#

distclean: clean
