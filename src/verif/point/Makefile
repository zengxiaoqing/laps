# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a \
           $(SRCROOT)/src/lib/util/libutil.a \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a \
             $(SRCROOT)/src/lib/util/libdebugutil.a \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a

.SUFFIXES:
.SUFFIXES: .o .c .f .a .exe

EXE = gen_verif_raob.exe \
      gen_verif_prof.exe

SRC = gen_verif_prof.f \
    bkgd_wgi.f \
    gen_verif_prof_sub.f \
    get_prof_pairs.f \
    gen_verif_raob.f \
    gen_verif_raob_sub.f \
    get_raob_data_a.f \
    get_raob_pairs.f \
    raob_subs.f
      
LIB = lapsverif.a

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
