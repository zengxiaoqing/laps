# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT = /usr/nfs/lapb/parallel/laps
DOMAIN   = /data/lapb/parallel/laps/data


include $(LAPSROOT)/src/include/makefile.inc










LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a  


DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a  


OTHERLIBS = /usr/local/apps/netcdf-3.3.1/lib/libnetcdf.a 

.SUFFIXES:
.SUFFIXES: .o .c .f .a .exe

EXE = build_sfc_static.exe  checkbi.exe 
SRC = build_sfc_static.f  checkbi.f






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
	  $(INSTALL) $$file $(LAPSROOT)/bin/$$file; \
	done

clean:
	rm -f $(EXE)  *.o *~ *#


distclean: clean
