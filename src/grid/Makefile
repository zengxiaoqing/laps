# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  $(UTILLIB) $(MODULESLIB)

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a $(DBUTILLIB) $(DBMODULESLIB) 

.SUFFIXES: .o .c .f .a .F .i .exe

EXE = gridgen_model.exe

SRC   = get_tile_list.f gridgen_utils.f put_laps_static.f \
        proc_geodat.f showgrid.f \
        create_bcd_bkgnd.c setsup.c \
        getarealim.c \
        xy_to_ll.c ll_to_xy.c 

MAIN=gridgen_model.f 

LIB=libgrid.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(FC) $(MAIN) $(FFLAGS) $(FMOD) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@ 
	rm -f *.o

debug: $(SRC) ;
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE) $(LIB) *.o *~ *# *.i *.mod


distclean: clean
