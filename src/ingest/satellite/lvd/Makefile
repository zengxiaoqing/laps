# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../../../..

include $(LAPSROOT)/src/include/makefile.inc

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a \
           $(LAPSROOT)/src/lib/satellite/libsatellite.a \
           $(LAPSROOT)/src/lib/goesinav/libgoesinav.a \
           $(LAPSROOT)/src/lib/nav/libnav.a \
           $(LAPSROOT)/src/lib/util/libutil.a

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a \
             $(LAPSROOT)/src/lib/satellite/libdebugsatellite.a \
             $(LAPSROOT)/src/lib/goesinav/libdebuggoesinav.a \
             $(LAPSROOT)/src/lib/nav/libdebugnav.a \
             $(LAPSROOT)/src/lib/util/libdebugutil.a

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = lvd_sat_ingest.exe

SRC = lvd_driver_sub.f \
  readcdf.f \
  readsatdat.f \
  read_ascii_satdat.f \
  read_gwc_satdat.f \
  readlut.f \
  read_gvarimg_cnt2btemp_lut.f \
  rdcdfblock.f \
  getcdf_satdat.f \
  getsat_attributes.f \
  getafgwc_satdat.f \
  bcd_to_int.f \
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
  check_nav_lut.f \
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
