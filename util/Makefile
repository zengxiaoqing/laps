# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=..


include $(LAPSROOT)/src/include/makefile.inc









LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a


DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a



OTHERLIBS = 

.SUFFIXES:
.SUFFIXES: .o .c .f .a .exe

EXE = dump_parms.exe
SRC = dump_parms.f 






all: exe

exe: $(EXE)

%.exe: $(SRC)
	$(FC) $(FFLAGS) $(SRC)  $(OTHERLIBS) -o $@
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
