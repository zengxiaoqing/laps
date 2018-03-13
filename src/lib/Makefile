# User options are in src/include/makefile.inc
# This Makefile is designed for gnu make version 3.75 or newer
# Please report problems to lapsbugs@fsl.noaa.gov
#
#

SRCROOT=../..


include $(SRCROOT)/src/include/makefile.inc


LIB      = liblaps.a
DEBUGLIB = libdebuglaps.a


#
# machine_dependent needs to be first
#

SRC	= machine_dependent.F \
          barnes_multivariate.F \
          cio.c \
	  rwl_v3.c \
          fort2c_str.c \
          static_routines.c \
          waiting_c.c  \
          fsystem_dum.f \
          add.f \
          addcon.f \
          advect.f \
          apply_qc_info.f \
          array_diagnosis.f \
          aski4tim_gg.f \
          atan3.f \
          atan3d.f \
          bilinear_interp.f \
          block_data.f \
          calc_contable.f \
          climate_sm.f \
          constant.f \
          contingency_table.f \
          convc2f.f \
          conversions.f \
          convf2c.f \
          convf2k.f \
          convk2f.f \
          convkt2ms.f \
          convm2miles.f \
          convmiles2m.f \
          convms2kt.f \
          csplit.f \
          diff.f \
          downcase.f \
          equtoaltazd.f \
          esat1.f \
          fflxc.f \
          filename13.f \
          filter_filenames.f \
          fname_conversions.f \
          gaussj.f \
          gen_bl_file.f \
          gen_btemps.f \
          get_barnes_superobs.f \
          get_cloud_deriv.f \
          get_dir_length.f \
          get_directory.f \
          get_domain_perimeter.f \
          get_file_names.f \
          get_file_time.f \
          get_heights_hydrostatic.f \
          get_ht_3d.f \
          get_interval_precip.f \
          get_local_towerobs.f \
          get_lowref.f \
          get_maps_lapsgrid.f \
          get_maxref.f \
          get_maxtops.f \
          get_pres_3d.f \
          get_radar_deriv.f \
          get_radar_max_pd.f \
          get_static_info.f \
          get_systime.f \
          getcld3d.f \
          getradar.f \
          get_time_window.f \
          gettmp2d.f \
          gettmp3d.f \
          getuv3d.f \
          getw3d.f \
          ghbry.f \
          gridconv.f \
          ground_temp.f \
          hum.f \
          i4time_asc_gg.f \
          i4time_fname_lp.f \
          ingest_madis_map.f \
          ingest_routines.f \
          insertpireps.f \
          insertsao.f \
          interp.f \
          interpolate_3dfield.f \
          land_albedo.f \
          land_topo.f \
          laps_io.f \
          laps_routines.f \
          laps_thermo.f \
          lapsgrid.f \
          latlon_db_rlaps.f \
          latlon_to_rlaps.f \
          libdummies.f \
          make_rh.f \
          make_ssh.f \
          make_td.f \
          make_optran_layers.f \
          map_proj.f \
          meanpres.f \
          move.f \
          multar.f \
          multcon.f \
          my_show_timer.f \
          nan_routines.f \
          newreadobs.f \
          normalize.f \
          ofm.f \
          oh_zone.f \
          pcpcnc.f \
          prof_cdf.f \
          psatoz.f \
          put_sfc_bal.f \
          qc.f \
          qc_field_3d.f \
          radar_qc.f \
          radar_to_latlon.f \
          razm_lat_lon_gm.f \
          readlaps.f \
          rd_laps_static.f \
          read-laps-static.f \
          read_acars_ob.f \
          read_gps_obs.f \
          read_ldad_prof.f \
          read_namelist.f \
          read_pro_data.f \
          read_snd_data.f \
          read_surface_obs.f \
          read_surface_old.f \
          read_surface_sa.f \
          read_wrfsi_nl.f \
          readlapscompressed.f \
          readlapsdata.f \
          readlapsheader.f \
          rfill.f \
          rh2cld.f \
          rootfind.f \
          rotate.f \
          satgeom.f \
          sfc_routines.f \
          skill_scores.f \
          smf_wisp91.f \
          smf_wisp92.f \
          smooth_box.f \
          smooth_cross_laps.f \
          snooze_gg.f \
          solar_normal.f \
          spread.f \
          xlbfgsb.f \
          xlinpack.f \
          xblas.f \
          xtimer_lbfgs.f \
          ssh2.f \
          stats_1d.f \
          stats_2d.f \
          subcon.f \
          sunfuncs.f \
          time_conversions.f \
          time_now.c \
          trig_dum.f \
          upcase.f \
          up_mflux.f \
          vectors.f \
          vis2albedo.f \
          vv.f \
          wait_for_data.f \
          wfo_fname_conversions.f \
          windspeed.f \
          writelaps.f \
          writelapscompressed.f \
          writelapsdata.f \
          writelapslvls.f \
          writelapslfm.f \
          wrt_laps_static.f \
          write_snd.f \
          xyztopolard.f \
          zero.f \
          zr.f \
          zs.f \
          ztopsa.f \
          ztotsa.f \
          getfilenames_c.c \
          i4time_now_gg.c \
          regex.c \
          alloca.c \
          xmalloc.c 



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

ifdef USECPP
%.o : %.F
	rm -f $*.f
	$(CPP) $(CPPFLAGS) $< 
	mv $*.i $*.f
	$(FC) $(FFLAGS) -c $*.f
endif

