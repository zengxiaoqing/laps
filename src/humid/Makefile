# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../..

include $(LAPSROOT)/src/include/makefile.inc

#FFLAGS+=-pg

DBFFLAGS+=-pg

LAPSLIBS = $(LAPSROOT)/src/lib/liblaps.a  \
           $(LAPSROOT)/src/lib/fm/libfm.a \
	   $(LAPSROOT)/src/lib/powell/libpowell.a \
	   $(LAPSROOT)/src/lib/mthermo/libmthermo.a   

DBLAPSLIBS = $(LAPSROOT)/src/lib/libdebuglaps.a   \
             $(LAPSROOT)/src/lib/fm/libdebugfm.a \
             $(LAPSROOT)/src/lib/powell/libdebugpowell.a \
	     $(LAPSROOT)/src/lib/mthermo/libdebugmthermo.a   

.SUFFIXES:
.SUFFIXES: .o .c .f .a .F .i .exe

EXE = lq3driver.x

SRC=analq.f            ghbry.f            lq3_driver1a.f     rsr.f \
    glsp.f             sh2mr.f            tiros.f \
    func_o.f             glst.f             lsin.f             slv_laplc.f \
    gen_bl_file.f      glstd.f            mak_cld_grid.f     smooth_grid2.f \
    get_laps_gp.f      goes_sbn.f         plotfield.f        weight_func.f \
    get_latest_file.f  int_tpw.f          prep_grid.f        write_lh4.f \
    get_rad.f          lh3_compress.f     raob_step.f        writef_sp.f \
    getmapsdf.f        locate.f           read_lvd_3_4_5.f   writefile.f \
    moment_b.f         two_d_stats.f      process_gvap.f     read_gvap.f \
    analz_gvap.f       report_change.f    ofm.f              oh_zone.f  \
    precip_water_profile.f                optrans.f	     bias_correction.f

MAIN=lq3driver.f

LIB=liblq3.a

ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))

all: exe

exe: $(ARCHIVE) $(EXE)



$(EXE): $(ARCHIVE) $(MAIN)
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
























