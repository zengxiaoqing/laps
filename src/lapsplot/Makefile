# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../..

include $(LAPSROOT)/src/include/makefile.inc

NCARGFC=ncargf77

FFLAGS=$(DBFLAGS)

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a  \
           $(LAPSROOT)/src/lib/mthermo/libmthermo.a  \
           $(LAPSROOT)/src/lib/util/libutil.a  \

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a  \
             $(LAPSROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(LAPSROOT)/src/lib/util/libdebugutil.a  \

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = lapsplot.exe

SRC= lapsplot_main.f \
     lapsplot_sub.f \
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
     supmap.f \
     drawctymap.f \
     weight_vertical.f \
     callsup_gp.f

MAIN=lapsplot_main.f

LIB=liblapsplot.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(NCARGFC) $(MAIN) $(DBFFLAGS) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
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

