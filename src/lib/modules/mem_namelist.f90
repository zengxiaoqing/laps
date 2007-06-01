Module mem_namelist

! nest7grid.parms variables

include 'lapsparms.for'

integer      iflag_lapsparms
real*4       max_radar_files_nl
real*4       PRESSURE_INTERVAL_L
real*4       PRESSURE_0_L
integer      nk_laps
real*4       standard_latitude
real*4       standard_latitude2
real*4       standard_longitude
integer      NX_L
integer      NY_L
integer      I_PERIMETER
real*4       grid_spacing_m
real*4       grid_cen_lat
real*4       grid_cen_lon
integer      laps_cycle_time
integer      min_to_wait_for_metars
integer      i_delta_sat_t_sec
integer      i2_missing_data
real*4       r_missing_data
integer      MAX_RADARS
real*4       ref_base
real*4       ref_base_useable
real*4       r_hybrid_first_gate
real*4       aircraft_time_window
integer      maxstns
integer      N_PIREP
integer      vert_rad_meso
integer      vert_rad_sao
integer      vert_rad_pirep
integer      vert_rad_prof
real*4       silavwt_parm
real*4       toptwvl_parm
integer      maxstations, maxobs


character*50 c50_lowres_directory
character*40 vertical_grid
character*6  c6_maproj
character*8  radarext_3d
character*8  radarext_3d_accum
character*200 path_to_raw_pirep
character*200 path_to_raw_rass
character*200 path_to_raw_profiler
character*200 path_to_raw_blprass
character*200 path_to_raw_blpprofiler
character*200 path_to_raw_satellite_cdf
character*200 path_to_raw_satellite_gvr
character*200 path_to_raw_sat_wfo_vis
character*200 path_to_raw_sat_wfo_i39
character*200 path_to_raw_sat_wfo_iwv
character*200 path_to_raw_sat_wfo_i11
character*200 path_to_raw_sat_wfo_i12
character*200 path_to_wsi_2d_radar
character*200 path_to_wsi_3d_radar
character*200 path_to_qc_acars
character*200 path_to_raw_raob
character*200 path_to_metar_data
character*200 path_to_local_data
character*200 path_to_buoy_data
character*9   fdda_model_source(maxbgmodels)
character*8   c8_project_common
character*8   c8_blpfmt_common
character*3   c_raddat_type
character*80  c80_description
character*200 path_to_topt30s
character*200 path_to_topt10m
character*200 path_to_pctl10m
character*200 path_to_soil2m
character*200 path_to_landuse30s
character*200 path_to_soiltype_top30s
character*200 path_to_soiltype_bot30s 
character*200 path_to_greenfrac 
character*200 path_to_soiltemp1deg
character*200 path_to_albedo
character*200 path_to_maxsnoalb
character*200 path_to_islope
character*200 path_to_sst


!character*100 path_to_background_model(MAX_BG_MODELS)
!integer laps_background_model(MAX_BG_MODELS)

logical*1    l_compress_radar,l_use_tamdar,l_3dvar,l_pad1

! wind_nl variables
logical :: l_use_raob, l_use_cdw, l_use_radial_vel
real    :: thresh_2_radarobs_lvl_unfltrd  &
          ,thresh_4_radarobs_lvl_unfltrd  &
          ,thresh_9_radarobs_lvl_unfltrd  &
          ,weight_bkg_const_wind  &
          ,weight_radar  &
          ,rms_thresh_wind  
integer :: max_pr,max_pr_levels,i_3d

Contains

!-------------------------------------------------------------------
subroutine read_namelists (namelist_name, filename)

implicit none

character(len=*), intent(in) :: namelist_name, filename


namelist /wind_nl/ l_use_raob, l_use_cdw, l_use_radial_vel  &
                  ,thresh_2_radarobs_lvl_unfltrd  &
                  ,thresh_4_radarobs_lvl_unfltrd  &
                  ,thresh_9_radarobs_lvl_unfltrd  &
                  ,weight_bkg_const_wind  &
                  ,weight_radar  &
                  ,rms_thresh_wind  &
                  ,max_pr,max_pr_levels,i_3d

namelist /lapsparms_NL/ iflag_lapsparms &
                  ,max_radar_files_nl,PRESSURE_INTERVAL_L &
                  ,nk_laps,standard_latitude,standard_latitude2    &    
                  ,standard_longitude,NX_L, NY_L, I_PERIMETER &
                  ,l_compress_radar,l_use_tamdar,l_3dvar &
                  ,grid_spacing_m,grid_cen_lat,grid_cen_lon &
                  ,laps_cycle_time, min_to_wait_for_metars &
                  ,i2_missing_data, r_missing_data, MAX_RADARS &
                  ,ref_base,ref_base_useable,r_hybrid_first_gate &
                  ,maxstns,N_PIREP &
                  ,vert_rad_meso,vert_rad_sao &
                  ,vert_rad_pirep,vert_rad_prof      &
                  ,silavwt_parm,toptwvl_parm &
                  ,vertical_grid,c50_lowres_directory,c6_maproj &
                  ,radarext_3d,radarext_3d_accum &
                  ,aircraft_time_window &
                  ,path_to_raw_pirep &
                  ,path_to_raw_rass,path_to_raw_profiler &
                  ,path_to_raw_blprass,path_to_raw_blpprofiler &
                  ,path_to_wsi_2d_radar,path_to_wsi_3d_radar &
                  ,path_to_qc_acars &
                  ,c8_project_common,c8_blpfmt_common &
                  ,c_raddat_type,c80_description &
                  ,path_to_topt30s ,path_to_topt10m &
                  ,path_to_soiltype_top30s, path_to_soiltype_bot30s &
                  ,path_to_landuse30s,path_to_greenfrac &
                  ,path_to_soiltemp1deg,path_to_albedo,path_to_maxsnoalb &
                  ,path_to_islope,path_to_sst,fdda_model_source


! open the namelist file name

open (12, file=filename, status='old')

! read the requested namelist

if (namelist_name == 'lapsparms') then

   read (12, lapsparms_nl)
   
   ! QC the input variables if desired
   !  .
   !  .
   !  .
   !  .

elseif (namelist_name == 'wind') then

   read (12, wind_nl)
   
   ! QC the input variables if desired
   !  .
   !  .
   !  .
   !  .

endif

close (12)

return
end subroutine

!---------------------------------------------------------------------

end Module
