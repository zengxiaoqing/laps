# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../..


include $(LAPSROOT)/src/include/makefile.inc










LAPSLIBS = 


DBLAPSLIBS = 


OTHERLIBS = /usr/local/apps/netcdf-3.3.1/lib/libnetcdf.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a .exe

EXE = xfer_laps.exe
SRC = xfer_laps.c fill_laps.c






all: exe

exe: $(EXE)

%.exe: $(SRC)
	$(CC) $(CFLAGS) $(SRC)  $(OTHERLIBS) -o $@
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
