# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../..

include $(SRCROOT)/src/include/makefile.inc

RM = rm -f
FLAG = $(FFLAGS)

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/modules/libmodules.a  \
           $(SRCROOT)/src/lib/mthermo/libmthermo.a  \
           $(SRCROOT)/src/lib/util/libutil.a  \
           $(SRCROOT)/src/lib/goeslib/libgoeslib.a 

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a  \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a  \
             $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a  \
             $(SRCROOT)/src/lib/util/libdebugutil.a  \
             $(SRCROOT)/src/lib/goeslib/libdebuggoeslib.a

LAPSFC = $(SRCROOT)/src/sfc/libsfc.a


.SUFFIXES:
.SUFFIXES: .o .c .f .a

EXE = recurs_iter.x
MOD = Definition.o \
      Util_Tools.o \
      Initialize.o \
      Minimizatn.o \
      ConfigLaps.o

INC =

SRC=MAIN.f90 \
    Writeout.f90 \
    routines.o \
    lso_reader_meso.o \
    get_laps_info.o 
#\
#    laps_routines.o

all: exe

exe: $(EXE)

$(EXE): $(MOD) $(SRC) $(LAPSLIBS)
	$(FC) $(FFLAGS) $(FREE) $(SRC) $(MOD) $(LAPSFC) $(LAPSLIBS) \
	$(OTHERLIBS) -o $@
	rm -f *.o

debug: $(SRC)
	echo "Hello " $(RM)
	$(RM) $(EXE) ;
	$(MAKE) LAPSLIBS='$(DBLAPSLIBS)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: exe
	$(INSTALL) $(EXE) $(INSTALLROOT)/bin/$(EXE)

clean:
	rm -f $(EXE)  *.o *~ *#


distclean: clean


#-------------------- SC MODULES:

DEFINITION = \
         Definition.f90 \
         Constants.f90 \
         Variables.f90 \
         OBSCommon.f90
Definition.o: $(DEFINITION)
	$(FC) -c $(FLAG) $(FREE) Definition.f90

UTIL_TOOLS = \
        Util_Tools.f90 \
        Intplt3d.f90
Util_Tools.o: $(UTIL_TOOLS)
	$(FC) -c $(FLAG) $(FREE) Util_Tools.f90

INITIALIZE = \
        Initialize.f90 \
        Namelist.f90 \
        ReadObsn.f90 \
        Grid2Obs.f90
Initialize.o: Definition.o $(INITIALIZE)
	$(FC) -c $(FLAG) $(INC) $(FREE) Initialize.f90

MINIMIZATN = \
        Minimizatn.f90 \
        Functn.f90 \
        Functn_ad.f90 \
        Iterates.f90 \
        Minimize.f90 \
        RF1D.f90 \
        RF1D_ad.f90 \
        RF3D.f90 \
        RF3D_ad.f90
Minimizatn.o: Definition.o $(MINIMIZATN)
	$(FC) -c $(FLAG) $(INC) $(FREE) Minimizatn.f90

CONFIGLAPS = \
	ConfigLaps.f90 \
	LapsInfo.f90 \
	LSO_Data_QC.f90 \
	WriteAnalysis.f90 \
	GridBarnes.f90
ConfigLaps.o: Definition.o $(CONFIGLAPS)
	$(FC) -c $(FLAG) -I$(SRCROOT)/src/lib/mthermo/ $(FREE) ConfigLaps.f90

routines.o: routines.f
	$(FC) -c $(FLAG) routines.f

lso_reader_meso.o: lso_reader_meso.f
	$(FC) -c $(FFLAGS) $(FIXED) lso_reader_meso.f
QC_tools.o: QC_tools.f
	$(FC) -c $(FFLAGS) $(FIXED) QC_tools.f
get_laps_info: get_laps_info.f
	$(FC) -c $(FFLAGS) $(FIXED) get_laps_info.f
