# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

FFLAGS=$(DBFLAGS)

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = lapsplot.exe

SRC= lapsplot_main.f \
     lapsplot_sub.f \
     array_range.f \
     aski4t.f \
     hsect.f \
     xsect.f \
     divergence.f \
     lapsplot.f \
     frame2.f \
     getuv2d.f \
     logit.f \
     plotobs.f \
     plottype.f \
     plotwinds_2d.f \
     plotgrid_2d.f \
     plotwindob.f \
     plot_mesoob.f \
     plot_contour.f \
     plot_sounding.f \
     supmap.f \
     drawctymap.f \
     weight_vertical.f \
     ccpfil.f \
     callsup_gp.f

MAIN=lapsplot_main.f

LIB=liblapsplot.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(NCARGFC) $(MAIN) $(DBFLAGS) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE) $(LIB) *.o *~ *# *.i


distclean: clean

%.f : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f

%.o : %.F

