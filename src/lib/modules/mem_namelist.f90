Module mem_namelist

include 'lapsparms.for'

!       Globally used variables that are independent of the namelists
!       character(len=200) :: generic_data_root, cstaticdir, grid_fnam_common

!       Declarations for namelist variables
        integer    iflag_lapsparms
        real       max_radar_files_nl   ! max_radar_files is in lapsparms.for
        real       PRESSURE_INTERVAL_L
        real       PRESSURE_0_L
        integer    nk_laps
        real       standard_latitude
        real       standard_latitude2
        real       standard_longitude
        integer    NX_L
        integer    NY_L
        integer    I_PERIMETER
        real grid_spacing_m
        real grid_cen_lat
        real grid_cen_lon
        real earth_radius 
        integer    laps_cycle_time, precip_cycle_time
        integer    model_cycle_time, model_fcst_intvl, model_fcst_len
        real purge_time
        integer i2_missing_data
        integer iverbose
        real    r_missing_data
        integer  MAX_RADARS
        real aod,aod_bin(3),aod_asy(3,3),fcterm,aod_ha,ht_ha(4),ssa(3)
        real fcterm_a(3),angstrom_exp_a
        real alpha_ha
        real o3_du,h_o3,d_o3
        real aero_scaleht
        real ref_base
        real ref_base_useable
        real r_hybrid_first_gate
        real aircraft_time_window
        real hydrometeor_scale_pcp
        real hydrometeor_scale_cld
        real solalt_thr_vis
        integer  maxstns
        integer  N_PIREP
        integer max_snd_grid
        integer max_snd_levels
        real redp_lvl
        real prtop
        integer vert_rad_meso
        integer vert_rad_sao
        integer vert_rad_pirep
        integer vert_rad_prof     
        real silavwt_parm
        real toptwvl_parm
        integer iwrite_output
        integer i_offset_radar

        character*40  vertical_grid
        character*10  lvl_coord_cdf
        character*50  c50_lowres_directory
        character*6   c6_maproj
        character*8   radarext_3d
        character*8   radarext_3d_accum
        character*200 path_to_raw_pirep
        character*200 path_to_raw_rass
        character*200 path_to_raw_profiler
        character*200 path_to_raw_blprass
        character*200 path_to_raw_blpprofiler
        character*200 path_to_wsi_2d_radar
        character*200 path_to_wsi_3d_radar
        character*200 path_to_qc_acars
        character*200 path_to_wisdom
        character*30  fdda_model_source(maxbgmodels)
        character*8   c8_project
        character*8   c8_blpfmt
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

        logical*1    l_compress_radar,l_dual_pol,l_use_tamdar,l_3dvar,l_pad1
        logical*1    l_accum_fg, l_accum_radar, l_accum_gauge
        logical*1    l_superob_barnes, l_mosaic_sat, l_fsf_gridgen

        common  /lapsparms/ iflag_lapsparms &
        ,max_radar_files_nl,PRESSURE_INTERVAL_L,PRESSURE_0_L &
        ,nk_laps,standard_latitude,standard_latitude2        &
        ,standard_longitude,NX_L, NY_L, I_PERIMETER          &
        ,grid_spacing_m,grid_cen_lat,grid_cen_lon            &
        ,earth_radius                                        &
        ,laps_cycle_time                                     &
        ,i2_missing_data, r_missing_data, MAX_RADARS         &
        ,ref_base,ref_base_useable,r_hybrid_first_gate       &
        ,maxstns,N_PIREP                                     &
        ,max_snd_grid,max_snd_levels                         &
        ,redp_lvl,prtop                                      &
        ,vert_rad_meso,vert_rad_sao                          &
        ,vert_rad_pirep,vert_rad_prof                        &
        ,silavwt_parm,toptwvl_parm                           &
        ,iwrite_output                                       &
        ,aircraft_time_window                                &
        ,vertical_grid,lvl_coord_cdf                         &
        ,c50_lowres_directory,c6_maproj                      &
        ,radarext_3d,radarext_3d_accum                       &
        ,path_to_raw_pirep                                   &
        ,path_to_raw_rass,path_to_raw_profiler               &
        ,path_to_raw_blprass,path_to_raw_blpprofiler         &
        ,path_to_wsi_2d_radar,path_to_wsi_3d_radar           &
        ,path_to_qc_acars,path_to_wisdom                     &
        ,c8_project,c8_blpfmt                                &
        ,c_raddat_type, c80_description                      &
        ,path_to_topt30s                                     &
        ,path_to_topt10m, path_to_pctl10m, path_to_soil2m    &
        ,path_to_landuse30s,path_to_soiltype_top30s          &
        ,path_to_soiltype_bot30s,path_to_greenfrac           &
        ,path_to_soiltemp1deg,path_to_albedo,path_to_sst     &
        ,path_to_maxsnoalb,path_to_islope                    &
        ,fdda_model_source                                   &
        ,l_compress_radar,l_dual_pol,l_use_tamdar,l_3dvar    &
        ,l_pad1,l_accum_fg,l_accum_radar,l_accum_gauge

! wind_nl variables
logical :: l_use_raob, l_use_cdw, l_use_radial_vel
real    :: weight_bkg_const_wind  &
          ,weight_radar  &
          ,rms_thresh_wind &
          ,stdev_thresh_radial &
          ,r0_barnes_max_m &
          ,brns_conv_rate_wind &
          ,qc_thresh_wind_def &
          ,qc_thresh_wind_pin &
          ,qc_thresh_wind_cdw &
          ,qc_thresh_wind_pro 
integer :: thresh_2_radarobs_lvl_unfltrd  &
          ,thresh_4_radarobs_lvl_unfltrd  &
          ,thresh_9_radarobs_lvl_unfltrd  &
          ,thresh_25_radarobs_lvl_unfltrd  
integer :: max_pr,max_pr_levels,max_wind_obs

! pressures_nl variables
integer, parameter   :: max_p=150
real                 :: pressures(max_p)
integer              :: nplevs

! surface_analysis variables
integer  ::  use_lso_qc,skip_internal_qc, itheta
logical  ::  l_require_lso, l_dev_ck
real     ::  del, gam, ak &
            ,bad_t, bad_td, bad_u, bad_v, bad_p  &
            ,bad_mp, bad_th, bad_the &
            ,bad_tgd_land, bad_tgd_water, bad_vis, bad_tb8  &
            ,thresh_t, thresh_td, thresh_mslp &
            ,rms_wind, rms_temp, rms_dewpoint, rms_pres

          !  redp_lvl utilized in background code also

! temp_nl variables
logical  :: l_read_raob_t, l_use_raob_t
real     :: weight_bkg_const_temp, pres_mix_thresh, rms_thresh_temp, radiometer_ht_temp
integer  :: max_obs, mode_adjust_heights

! cloud_nl variables
logical  :: l_use_vis,l_use_vis_add,l_use_vis_partial,l_use_s8a,l_use_39 &
           ,l_use_pireps ,l_use_metars, l_use_radar, l_corr_parallax
integer  :: latency_co2,i4_sat_window,i4_sat_window_offset,i_varadj
real     :: pct_req_lvd_s8a, echotop_thr_a(3), cld_weight_modelfg

! moisture_switch_nl variables
integer :: covar_switch, print_switch, raob_switch, radiometer_switch, raob_lookback, endian, goes_switch &
          ,cloud_switch, cloud_d, sounder_switch, tiros_switch, sat_skip &
          ,gvap_switch, ihop_flag, time_diff, gps_switch, sfc_mix &
          ,mod_4dda_1
real    :: raob_radius, mod_4dda_factor, t_ref
real    :: max_cdelrh_nl
real    :: cf_set_nl
real    :: cloud_weight_nl
real    :: radio_wt_nl
character(len=256) ::  path_to_gvap12, path_to_gvap10, path_to_gps, path2covar

! lapsprep_nl variables (not yet used)
logical ::  hotstart, balance, make_sfc_uv 
character(len=16) :: output_format(10)
real ::  snow_thresh, lwc2vapor_thresh,rai_frac, sno_frac
         ! RAMS processing variables
real ::  sfcinf
character(len=256) :: var_prefix




Contains

!-------------------------------------------------------------------

subroutine read_namelist_laps (namelist_name, filename)

!use control_coms
!use mem_grid
!use lapsparms_com
!use isan_coms

implicit none

character(len=*), intent(in) :: namelist_name, filename


namelist /lapsparms_NL/ iflag_lapsparms &
                  ,max_radar_files_nl,PRESSURE_INTERVAL_L &
                  ,nk_laps,standard_latitude,standard_latitude2    &    
                  ,standard_longitude,NX_L, NY_L, I_PERIMETER &
                  ,l_compress_radar,l_dual_pol,l_use_tamdar,l_3dvar &
                  ,l_accum_fg, l_accum_radar, l_accum_gauge &
                  ,l_superob_barnes, l_mosaic_sat &
                  ,grid_spacing_m,grid_cen_lat,grid_cen_lon &
                  ,earth_radius &
                  ,laps_cycle_time,precip_cycle_time &
                  ,model_cycle_time, model_fcst_intvl, model_fcst_len &
                  ,purge_time &
                  ,i2_missing_data, iverbose, r_missing_data, MAX_RADARS, i_offset_radar &
                  ,aod,aod_bin,aod_asy,fcterm,aod_ha,ht_ha,ssa,aero_scaleht,o3_du,h_o3,d_o3 &
                  ,fcterm_a,angstrom_exp_a &
                  ,ref_base,ref_base_useable,r_hybrid_first_gate &
                  ,maxstns,N_PIREP &
                  ,max_snd_grid,max_snd_levels,redp_lvl,prtop &
                  ,hydrometeor_scale_pcp,hydrometeor_scale_cld,solalt_thr_vis &
                  ,vert_rad_meso,vert_rad_sao &
                  ,vert_rad_pirep,vert_rad_prof      &
                  ,silavwt_parm,toptwvl_parm &
                  ,iwrite_output &
                  ,aircraft_time_window &
                  ,vertical_grid,lvl_coord_cdf &
                  ,c50_lowres_directory,c6_maproj &
                  ,radarext_3d,radarext_3d_accum &
                  ,path_to_raw_pirep &
                  ,path_to_raw_rass,path_to_raw_profiler &
                  ,path_to_raw_blprass,path_to_raw_blpprofiler &
                  ,path_to_wsi_2d_radar,path_to_wsi_3d_radar &
                  ,path_to_qc_acars,path_to_wisdom &
                  ,c8_project,c8_blpfmt &
                  ,c_raddat_type, c80_description &
                  ,l_fsf_gridgen &
                  ,path_to_topt30s ,path_to_topt10m &
                  ,path_to_soiltype_top30s, path_to_soiltype_bot30s &
                  ,path_to_landuse30s,path_to_greenfrac &
                  ,path_to_soiltemp1deg,path_to_albedo,path_to_maxsnoalb &
                  ,path_to_islope,path_to_sst,fdda_model_source

namelist/pressures_nl/ pressures

namelist /wind_nl/ l_use_raob, l_use_cdw, l_use_radial_vel  &
                  ,thresh_2_radarobs_lvl_unfltrd  &
                  ,thresh_4_radarobs_lvl_unfltrd  &
                  ,thresh_9_radarobs_lvl_unfltrd  &
                  ,thresh_25_radarobs_lvl_unfltrd  &
                  ,stdev_thresh_radial &
                  ,weight_bkg_const_wind  &
                  ,weight_radar  &
                  ,rms_thresh_wind  &
                  ,max_pr,max_pr_levels,max_wind_obs &
                  ,r0_barnes_max_m &
                  ,brns_conv_rate_wind &
                  ,qc_thresh_wind_def &
                  ,qc_thresh_wind_pin &
                  ,qc_thresh_wind_cdw &
                  ,qc_thresh_wind_pro 

namelist /surface_analysis/  &
                  use_lso_qc,skip_internal_qc, itheta  &
                  ,l_require_lso, l_dev_ck  &
                  ,del, gam, ak &
                  ,bad_t, bad_td, bad_u, bad_v, bad_p  &
                  ,bad_mp, bad_th, bad_the &
                  ,bad_tgd_land, bad_tgd_water, bad_vis, bad_tb8  &
                  ,thresh_t, thresh_td, thresh_mslp  &
                  ,rms_wind, rms_temp, rms_dewpoint, rms_pres
                  
namelist /temp_nl/ l_read_raob_t, l_use_raob_t, mode_adjust_heights  &
                  ,weight_bkg_const_temp, pres_mix_thresh, rms_thresh_temp &
                  ,max_obs, radiometer_ht_temp

namelist /cloud_nl/ l_use_vis, l_use_vis_add, l_use_vis_partial &
                   ,l_use_s8a, l_use_pireps, l_use_39, l_use_metars, l_use_radar, l_corr_parallax &
                   ,latency_co2 &
                   ,pct_req_lvd_s8a, echotop_thr_a &
                   ,cld_weight_modelfg &
                   ,i4_sat_window,i4_sat_window_offset &
                   ,i_varadj

namelist /moisture_switch_nl/ &
                   covar_switch, print_switch, raob_switch, radiometer_switch, raob_lookback, endian, raob_radius  &
                  ,goes_switch ,cloud_switch, cloud_d, sounder_switch  &
                  ,tiros_switch, sat_skip, gvap_switch  &
                  ,max_cdelrh_nl,cf_set_nl,cloud_weight_nl,radio_wt_nl & 
                  ,ihop_flag, time_diff, gps_switch, sfc_mix, mod_4dda_1  &
                  ,raob_radius, mod_4dda_factor, t_ref  &
                  ,path_to_gvap12, path_to_gvap10, path_to_gps, path2covar

namelist/lapsprep_nl/ var_prefix, sfcinf, hotstart, balance  &
                     ,output_format, snow_thresh, lwc2vapor_thresh  &
                     ,make_sfc_uv,rai_frac, sno_frac
                  

print*
print*,'======> Read_namelist_laps: ',trim(namelist_name)
print*,'======>          File_name: ',trim(filename)

! open the namelist file name

open (12, file=filename, status='old')

! read the requested namelist
if (namelist_name == 'ilaps_control') then
   ! Set some default values
!  proc_grids = 0

   ! Read ILAPS options information
!  read(12,ilaps_control)

elseif (namelist_name == 'RAMS') then
   ! Read RAMS grid point information
!  open(1,status='OLD',file=rams_grid_input_file)
!  read(1,model_grids)
!  close(1)
   
   ! Since there is more stuff in the LAPS nest7parm file, read this also. 
   !  Then overwrite LAPS grid params with RAMS grid params later.
   !open(1,status='OLD',file=laps_grid_input_file)
   !open(1,status='OLD',file=filename)
   print*,'======> Reading namelist: lapsparms_nl',nk_laps,nx_l,ny_l
   read(12,lapsparms_nl)
   print*,'======> Reading namelist: lapsparms_nl',nk_laps,nx_l,ny_l
   !close(1)
   
elseif (namelist_name == 'LAPS') then
   ! Set some default values
!  proc_grids = 0 ; proc_grids(1) = 1 
   
   ! Read LAPS grid point information
   !open(1,status='OLD',file=laps_grid_input_file)
   read(12,lapsparms_nl)
   !close(1)
   
elseif (namelist_name == 'pressures') then
   read(12,pressures_nl)
   
elseif (namelist_name == 'lapsparms') then

!  default values
   earth_radius = 6370000. ! WRF value
   lvl_coord_cdf = 'HPA'
   prtop = 10000.          ! 100mb for SIGMA_P grid
   l_superob_barnes = .false.
   l_mosaic_sat     = .false.
   l_dual_pol       = .false.
   l_fsf_gridgen    = .false.
   iverbose = 0
   i_offset_radar = -1
   precip_cycle_time = -1
   solalt_thr_vis = 15.
   aod = 0.05              ! default column aerosol optical depth (550nm)
   angstrom_exp_a = -99.   ! aerosol angstrom exponent (flag value)
   aero_scaleht = 1500.    ! default aerosol scale height (m)
   fcterm = 0.05           ! range from 0.00 to 0.09 for large aerosols
                           ! corresponding phase function peak from 20-110
   fcterm_a(:) = -99.      ! coarse model term for pair of DHG functions (each wavelength - flag value)
   aod_ha = .004           ! high altitude aerosol optical depth
   ht_ha(1)=13000.; ht_ha(2)=13000.; ht_ha(3)=25000.; ht_ha(4)=31000.
   o3_du = 300.
   h_o3 = 22000.
   d_o3 = 11000.

!  Single scattering albedo for aerosols
   ssa(1) = .90  ; ssa(2) = .90 ; ssa(3) = .90 ! non-dust
!  ssa(1) = .92  ; ssa(2) = .92 ; ssa(3) = .92 ! urban
!  ssa(1) = .99  ; ssa(2) = .99 ; ssa(3) = .99 ! non-absorbing (sea salt)
!  ssa(1) = .935 ; ssa(2) = .92 ; ssa(3) = .86 ! hematite dust
!  ssa(1) = .95  ; ssa(2) = .85 ; ssa(3) = .75 ! smoke

!  fraction of aerosols in each bin & asymmetry factor
!  Factor of 2 back scatter increase from minimum, peak of 50
   aod_bin(1) = 0.000 
   aod_bin(2) = 0.987 
   aod_bin(3) = 0.013 

!  Coarse mode aerosols (rgb)
   aod_asy(1,1) = +0.957 ; aod_asy(1,2) = +0.962 ; aod_asy(1,3) = +0.967

!  Fine mode aerosols (rgb)
   aod_asy(2,1) = +0.49  ; aod_asy(2,2) = +0.50  ; aod_asy(2,3) = +0.51 

!  Backscattering aerosols
   aod_asy(3,:) = +0.55 

!  Factor of 10 back scatter increase from minimum
!  aod_bin(1) = 0.70   ;   aod_asy(1) = +0.65  
!  aod_bin(2) = 0.12   ;   aod_asy(2) = +0.95      
!  aod_bin(3) = 0.18   ;   aod_asy(3) = -0.65      

!  Factor of ~3 back scatter increase from minimum
!  aod_bin(1) = 0.80   ;   aod_asy(1) = +0.65  
!  aod_bin(2) = 0.14   ;   aod_asy(2) = +0.95      
!  aod_bin(3) = 0.06   ;   aod_asy(3) = -0.65      

!  Original
!  aod_bin(1) = 0.85   ;   aod_asy(1) = +0.65  
!  aod_bin(2) = 0.15   ;   aod_asy(2) = +0.95      
!  aod_bin(3) = 0.00   ;   aod_asy(3) = -0.65      

   read (12, lapsparms_nl, err=901)
   write(6,1)nk_laps,nx_l,ny_l
1  format(' ======> Read from namelist: lapsparms_nl',3i8)
!  print*,'======> Read from namelist: lapsparms_nl',nk_laps,nx_l,ny_l
   
   ! QC the input variables if desired
   !  .
   !  .
   !  .
   !  .
   if(precip_cycle_time .eq. -1)precip_cycle_time = laps_cycle_time

   alpha_ha = aod_ha / ((ht_ha(3)-ht_ha(2))+0.5*(ht_ha(4)-ht_ha(3)))

elseif (namelist_name == 'wind') then

   thresh_25_radarobs_lvl_unfltrd = 450 ! default used mainly for transition

   read (12, wind_nl, err=905)
   
   ! QC the input variables if desired
   if(r0_barnes_max_m .le. 0)then
       write(6,*)' Error in value of r0_barnes_max_m ',r0_barnes_max_m
       goto 905
   endif

elseif (namelist_name == 'sfc_anal') then

   read (12, surface_analysis, err=906)
   
   ! QC the input variables if desired
   !  .
   !  .
   !  .
   !  .

elseif (namelist_name == 'temp_anal') then

   radiometer_ht_temp = 2000.

   read (12, temp_nl, err=907)
   
   ! QC the input variables if desired
   !  .
   !  .
   !  .
   !  .


elseif (namelist_name == 'cloud_anal') then

   i_varadj = 1                    ! keep this on as a default          
   l_corr_parallax = .true.
   l_use_s8a = .true.
   l_use_pireps = .true.

   read (12, cloud_nl, err=908)
   
   ! QC the input variables if desired
   !  .
   !  .
   !  .
   !  .

elseif (namelist_name == 'moisture_anal') then

   read (12, moisture_switch_nl)
   
   ! QC the input variables if desired
   !  .
   !  .
   !  .
   !  .

elseif (namelist_name == 'lapsprep') then
   ! Set some default values
   output_format = ' '
   
   ! Read LAPSPREP info for RAMS-vfile processing
   read(12,lapsprep_nl)
   
else
   print*,'Illegal namelist_name in read_namelist_laps:', trim(namelist_name)
   stop 'read_namelist_laps: illegal namelist_name'
endif

close (12)

print*

return

901  print*,'error reading lapsparms_nl'
     write(*,lapsparms_nl)
     stop

905  print*,'error reading wind_nl'
     write(*,wind_nl)
     stop

906  print*,'error reading surface_analysis'
     write(*,surface_analysis)
     stop

907  print*,'error reading temp_anal'
     write(*,temp_nl)
     stop

908  print*,'error reading cloud_anal'
     write(*,cloud_nl)
     stop

909  print*,'error reading moisture_anal'
     write(*,moisture_switch_nl)
     stop

end subroutine

!---------------------------------------------------------------------

end Module
