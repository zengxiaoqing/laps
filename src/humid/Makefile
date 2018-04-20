# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..

include $(SRCROOT)/src/include/makefile.inc

#FFLAGS+=-pg

DBFFLAGS+=-pg

LAPSLIBS = $(SRCROOT)/src/lib/liblaps.a  \
           $(SRCROOT)/src/lib/fm/libfm.a \
	   $(SRCROOT)/src/lib/powell/libpowell.a \
           $(SRCROOT)/src/lib/modules/libmodules.a \
           $(SRCROOT)/src/lib/opt90/libopt90.a \
	   $(SRCROOT)/src/lib/mthermo/libmthermo.a   

DBLAPSLIBS = $(SRCROOT)/src/lib/libdebuglaps.a   \
             $(SRCROOT)/src/lib/fm/libdebugfm.a \
             $(SRCROOT)/src/lib/powell/libdebugpowell.a \
             $(SRCROOT)/src/lib/modules/libdebugmodules.a \
             $(SRCROOT)/src/lib/opt90/libdebugopt90.a \
	     $(SRCROOT)/src/lib/mthermo/libdebugmthermo.a   

.SUFFIXES:
.SUFFIXES: .o .c .f .a .i .exe

EXE = lq3driver.x

SRC=module_sfc_structure.f    analq.f            lq3_driver1a.f     rsr.f \
    glsp.f             sh2mr.f            tiros.f \
                         glst.f             lsin.f             slv_laplc.f \
    glstd.f            mak_cld_grid.f     smooth_grid2.f \
    get_laps_gp.f      variational.f         plotfield.f        weight_func.f \
    get_latest_file.f  int_tpw.f          prep_grid.f        write_lh4.f \
    get_rad.f          lh3_compress.f     writef_sp.f \
    getmapsdf.f        locate.f           read_lvd_3_4_5.f   writefile.f \
    moment_b.f         two_d_stats.f      process_gvap.f     read_gvap.f \
    analz_gvap.f       hypso.f            report_change.f \
    process_gps.f      read_gps.f         int_layerpw.f      cloud_sat.f \
    get_newest_file.f  weight_field.f     analz_gps.f        sigma_to_laps.f \
    int_ipw.f          test_cloud.f	  read_madis_gps.f  \
    precip_water_profile.f \
    bias_correction.f  snd_step.f         analz_snd.f       impact_assess.f

MAIN=lq3driver.f

LIB=liblq3.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)



$(EXE): $(ARCHIVE) $(MAIN) $(LAPSLIBS)
	$(RANLIB) $(LIB)
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

























