# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/cloud/libcloud.a  \
           $(SRCROOT)/src/lib/temp/libtemp.a  \
           $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/bgdata/libbgdata.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a  

DBLAPSLIBS = $(SRCROOT)/src/lib/cloud/libdebugcloud.a  \
             $(SRCROOT)/src/lib/temp/libdebugtemp.a  \
             $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/bgdata/libdebugbgdata.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = lapsplot.exe

SRC= lapsplot_main.f \
     lapsplot_sub.f \
     array_range.f \
     aski4t.f \
     get_lapsplot_parms.f \
     hsect.f \
     hsect_img.f \
     xsect.f \
     divergence.f \
     lapsplot.f \
     frame2.f \
     getuv2d.f \
     hsv_to_rgb.f \
     logit.f \
     plotobs.f \
     plot_il_obs.f \
     read_tdsnd.f \
     plot_allsky.f \
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
     get_fim_data.f \
     get_rams_data.f \
     callsup_gp.f

MAIN=lapsplot_main.f

LIB=liblapsplot.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(NCARGFC) $(MAIN) $(FFLAGS) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
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

