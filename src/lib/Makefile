# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

LAPSROOT=../..


include $(LAPSROOT)/src/include/makefile.inc








.SUFFIXES:
.SUFFIXES: .o .c .f .a

LIB      = liblaps.a
DEBUGLIB = libdebuglaps.a




SRC	= cio.c \
	  rwl_v3.c \
          fort2c_str.c \
          static_routines.c \
          waiting_c.c  \
          fsystem_dum.f \
          add.f \
          addcon.f \
          advect.f \
          aski4tim_gg.f \
          array_diagnosis.f \
          atan3.f \
          atan3d.f \
          barnes_multivariate.f \
          barnes_univariate.f \
          bilinear_interp.f \
          block_data.f \
          c_time2fname.c \
          check_nan.f \
          check_nan1.f \
          check_nan2.f \
          check_nan3.f \
          constant.f \
          convc2f.f \
          conversions.f \
          convf2c.f \
          convf2k.f \
          convk2f.f \
          convkt2ms.f \
          convm2miles.f \
          convmiles2m.f \
          convms2kt.f \
          cv_asc_i4time.f \
          cv_i4tim_asc_lp.f \
          diff.f \
          downcase.f \
          esat1.f \
          fflxc.f \
          filename13.f \
          filter_filenames.f \
          fname_conversions.f \
          get_cloud_deriv.f \
          get_dir_length.f \
          get_directory.f \
          get_domain_perimeter.f \
          get_file_names.f \
          get_file_time.f \
          get_heights_hydrostatic.f \
          get_lowref.f \
          get_maps_lapsgrid.f \
          get_maxref.f \
          get_maxtops.f \
          get_pres_3d.f \
          get_radar_max_pd.f \
          get_static_info.f \
          get_systime.f \
          getcld3d.f \
          getradar.f \
          gettmp2d.f \
          gettmp3d.f \
          getuv3d.f \
          getw3d.f \
          gridconv.f \
          ground_temp.f \
          hum.f \
          i4time_asc_gg.f \
          i4time_fname_lp.f \
          interp.f \
          interpolate_3dfield.f \
          laps_io.f \
          laps_thermo.f \
          lapsgrid.f \
          latlon_to_rlaps.f \
          libdummies.f \
          machine_dependent.F \
          make_rh.f \
          make_ssh.f \
          make_td.f \
          make_fnam_lp.f \
          map_proj.f \
          meanlapse.f \
          meanpres.f \
          move.f \
          multar.f \
          multcon.f \
          my_show_timer.f \
          newreadobs.f \
          normalize.f \
          pcpcnc.f \
          prof_cdf.f \
          psatoz.f \
          qc_field_3d.f \
          radar_qc.f \
          radar_to_latlon.f \
          razm_lat_lon_gm.f \
          readlaps.f \
          rd_laps_static.f \
          read-laps-static.f \
          read_surface_obs.f \
          read_surface_old.f \
          readlapsdata.f \
          readlapsheader.f \
          rfill.f \
          rootfind.f \
          smf_wisp91.f \
          smf_wisp92.f \
          smooth_cross_laps.f \
          snooze_gg.f \
          ssh2.f \
          subcon.f \
          sunfuncs.f \
          time_now.c \
          trig_dum.f \
          upcase.f \
          vectors.f \
          vis2albedo.f \
          visg2log.f \
          viss2log.f \
          vlog2vis.f \
          vv.f \
          wait_for_data.f \
          wfo_fname_conversions.f \
          windspeed.f \
          writelaps.f \
          writelapsdata.f \
          wrt_laps_static.f \
          zero.f \
          zr.f \
          zs.f \
          ztopsa.f \
          getfilenames_c.c \
          i4time_now_gg.c \
          regex.c \
          alloca.c \
          xmalloc.c \
          nan.c



ARCHIVE = $(LIB)($(addsuffix .o, $(basename $(SRC))))
DBARCHIVE = $(DEBUGLIB)($(addsuffix .o, $(basename $(SRC))))


all: lib

lib: $(LIB)

$(LIB): $(ARCHIVE) 
	$(RANLIB) $@

debug: $(SRC)
	$(MAKE) lib LIB='$(DEBUGLIB)' ARCHIVE='$(DBARCHIVE)' \
	CFLAGS='$(DBCFLAGS)' FFLAGS='$(DBFLAGS)'

install: lib

clean:
	$(RM) $(LIB) $(DEBUGLIB) *.o *~ *#

distclean: clean
%.f : %.F
	$(RM) $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f






































