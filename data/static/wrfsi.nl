&project_id
 simulation_name = 'wrfsi development'
 user_desc = 'FSL WRFSI developer - LAPB Group'
/

&filetimespec
 start_year = 2000
 start_month = 09
 start_day = 13
 start_hour = 12
 start_minute   =    00                              
 start_second   =    00
 end_year = 2000
 end_month = 09
 end_day = 14
 end_hour = 03
 end_minute     =    00
 end_second     =    00
 interval       =  10800                       
/

&hgridspec
 num_domains = 1
 xdim =  125,100,200,
 ydim =  105,100,200,
 parent_id = 1,1,2,
 ratio_to_parent = 1,4,4,
 domain_origin_parent_x = 1,36,36,
 domain_origin_parent_y = 1,36,36,
 stagger_type = 'A-c',
 map_proj_name = 'polar',
 latlon_grid = '/data/fxa/laps/lapsprd/static/latlon.dat',
 moad_known_lat = 41.0,
 moad_known_lon = -105.5 ,
 moad_known_loc = 'center',
 moad_stand_lats = 40.0, 90.0,
 moad_stand_lons = -105.5,
 moad_delta_x = 10000.,
 moad_delta_y = 10000.,
 silavwt_parm = 0.,
 toptwvl_parm = 2.,
/

&sfcfiles
 topo_30s = '/projects/oplapb/geog/model/topo_30s',
 topo_10m = '/projects/oplapb/geog/model/topo_10m',
 pctland_10m = '/projects/oplapb/geog/model/land_10m',
/

&hinterp_control
  num_domains_to_proc = 1,
  domain_id_list = 1,
  ptop_in_Pa = 10000,
  new_levels_in_Pa =  ,
  sst_to_ice_threshold = -9999,
  linear_interpolation = .false., 
  input_root = 'AVN_FILE',
  terrain_file_name = 'static.wrfsi',
  constants_full_name = 'SST_FILE:LATEST','SNOW_FILE:LATEST'
  output_prefix = 'hinterp'
  verbose = .true.,    
  static_in_output = .true.,
/

&vinterp_control
 input_prefix = 'hinterp'
 output_prefix = 'wrf_input'
 output_coord = 'ZETA',
 stagger_w = .true. 
 num_levels = 31,
 vertical_increment = 100., 
 vertical_stretch = 1.1,
 max_vertical_inc = 1000.,
 max_top = 15000.,
 use_specified_levels = .false.
 levels = 00000., 00050., 00100., 00150., 00200., 
          00250., 00300., 00350., 00400., 00450.,
          00500., 00600., 00700., 00800., 00900.,
          01000., 01250., 01500., 01750., 02000.,
          02500., 03000., 03500., 04000., 04500.,
          05000., 06000., 07000., 08000., 09000.,
          10000., 11000., 12000., 13000., 14000.,
          15000., 
 rh_wrt_liquid = .true.
 verbose = .true. 
 output_nonwrf =.false.                          
 print_setup_only = .false.
/

&paths_to_raw_data
  path_to_raw_pirep = '/public/data/pirep/netcdf/',
  path_to_raw_rass = '/public/data/profiler/rass/noaanet/netcdf/',
  path_to_raw_profiler = '/public/data/profiler/wind/noaanet/netcdf/',
  path_to_raw_blprass = '/public/data/profiler/rass/external/netcdf/',
  path_to_raw_blpprofiler = '/public/data/profiler/wind/external/netcdf/',
  path_to_wsi_2d_radar = '/public/data/radar/wsi/nowrad/netcdf/'
  path_to_wsi_3d_radar = '/public/data/radar/wsi/nexrad/netcdf/',
  path_to_qc_acars = '/public/data/acars/qc/netcdf/'
/

&analysis_control
  c_analysis_type = 'laps',
/
 
&laps_analysis_control
  nx_l=125,
  ny_l=105,
  nz_l=21,
  c_vcoordinate='pressure',
  grid_spacing_m=10000.,
  PRESSURE_BOTTOM=110000.,          
  PRESSURE_INTERVAL=5000.,
  laps_cycle_time=3600,           
  l_highres_laps=.false.,                  
  I_PERIMETER=10,                  
  c50_lowres_dir='NULL',
  craddat_type='wsi',
  radarext_3d='vrc',
  radarext_3d_accum='vrc',
  i2_missing_data=-99,
  r_missing_data=+1e37,
  MAX_RADARS=3,        
  ref_base=-10.,       
  ref_base_useable=0., 
  maxstns=900,         
  N_PIREP=3000,         
  vert_rad_meso=1,     
  vert_rad_sao=1,      
  vert_rad_pirep=1,    
  vert_rad_prof=1,     
  silavwt_parm=0.,     
  toptwvl_parm=2.,     
  c8_project='NIMBUS', 
  fdda_model_source=' ',
/

&si_paths
 extdataroot = '/scratch/oplapb'
 analpath  = '/projects/oplapb/ANAL'
 path_to_AVN = '/public/data/grids/avn/global-65160/grib'
 path_to_ETA = '/public/data/grids/eta/40km_eta212_isobaric/grib'
/

&si_controls
 avail_after_AVN = 6       ! hrs between time of fcst start and availability
 cyc_AVN = 6               ! hrs between model runs
 avail_after_ETA = 4
 cyc_ETA = 6
 fddalen = 0
 fcstlen = 15
 nestinit = 1
/
