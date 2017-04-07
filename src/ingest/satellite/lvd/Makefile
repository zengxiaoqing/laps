# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../../../..

include $(SRCROOT)/src/include/makefile.inc

LAPSLIBS = $(SRCROOT)/src/lib/satellite/libsatellite.a \
           $(SRCROOT)/src/lib/liblaps.a \
           $(SRCROOT)/src/lib/modules/libmodules.a \
           $(SRCROOT)/src/lib/goesinav/libgoesinav.a \
           $(SRCROOT)/src/lib/nav/libnav.a \
           $(SRCROOT)/src/lib/util/libutil.a

DBLAPSLIBS = $(SRCROOT)/src/lib/satellite/libdebugsatellite.a \
             $(SRCROOT)/src/lib/libdebuglaps.a \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/goesinav/libdebuggoesinav.a \
             $(SRCROOT)/src/lib/nav/libdebugnav.a \
             $(SRCROOT)/src/lib/util/libdebugutil.a

.SUFFIXES:
.SUFFIXES: .o .c .f .f90 .a .F .i .exe
.f90.o:;    $(FC) $(FLAG) $(FREE) -c $*.f90

EXE = lvd_sat_ingest.exe

SRC = lvd_driver_sub.f \
  readcdf.f \
  readsatdat.f \
  read_ascii_satdat.f \
  read_gwc_satdat.f \
  read_gms_hko.f \
  read_gms_image.c \
  read_gms_taiwan.f \
  readlut.f \
  read_gvarimg_cnt2btemp_lut.f \
  rdcdfblock.f \
  COMScount2tbNrad_sub.f90 \
  get_cldtop_p_sub.f \
  getcdf_satdat.f \
  getsat_attributes.f \
  getafgwc_satdat.f \
  compute_gvarimg_res.f \
  c_time_fnc.f \
  ldlapsdata.f \
  Process_SDHS_GVAR_sub.f \
  satir_process.f \
  satvis_process.f \
  satir2laps.f \
  satvis2laps.f \
  set_missing_sat.f \
  setmissingflag.f \
  satdatfill.f \
  satfill1.f \
  imagecomp.f \
  check_data.f \
  compute_nav_llij.f \
  latlon_to_grij.f \
  genbtemplut.f \
  btempconvert.f \
  btempconv_asc.f \
  wait_for_satdat.f \
  wait_for_wfo_satdat.f \
  setvisflag.f \
  decellularize_gwc_image.f

MAIN=lvd_sat_ingest.f

LIB=liblvd.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)

$(EXE): $(ARCHIVE) $(MAIN)
	$(FC) $(MAIN) $(FFLAGS) $(LIB) $(LAPSLIBS) $(OTHERLIBS) -o $@
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
