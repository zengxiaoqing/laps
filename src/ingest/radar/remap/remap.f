cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
      
      program remap

      include 'remap_constants.dat'      

      character path_to_radar*150, laps_radar_ext*4
     1         ,radar_subdir_dum*3, path_to_vrc*15

      character c_radar_start*10,c_radar_end*10,agroup*10

      character dataroot*255

      real, allocatable, dimension(:,:) :: lat,lon,topo,dum_2d
       
      logical l_realtime

      call get_grid_dim_xy(NX_L,NY_L,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting horizontal domain dimensions'
          go to 999
      endif

      call get_laps_dimensions(NZ_L,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting vertical domain dimensions'
          go to 999
      endif

      allocate(lat(NX_L,NY_L),STAT=istat_alloc)       
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate lat'
          stop
      endif

      allocate(lon(NX_L,NY_L),STAT=istat_alloc)       
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate lon'
          stop
      endif

      allocate(topo(NX_L,NY_L),STAT=istat_alloc)       
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate topo'
          stop
      endif

      allocate(dum_2d(NX_L,NY_L),STAT=istat_alloc)       
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate dum_2d'
          stop
      endif

!     This first call returns only 'n_radars_remap' and 'namelist_parms'
      call get_remap_parms(0,n_radars_remap,max_times,path_to_radar     ! I/O
     1       ,laps_radar_ext,radar_subdir_dum,path_to_vrc,ref_min       ! O
     1       ,min_ref_samples,min_vel_samples,dgr,namelist_parms        ! O
     1       ,istatus)                                                  ! O
      if(istatus .ne. 1)then
          write(6,*)'Warning: bad status return from get_remap_parms'       
          go to 999
      endif

      if(n_radars_remap .eq. -1)then
          call get_max_radars(n_radars_remap,istatus)
          if(istatus .ne. 1)goto999
          write(6,*)' setting n_radars_remap = max_radars_cmn = '
     1             ,n_radars_remap
      endif

      if(namelist_parms%n_groups .gt. 1)then
          n_groups = namelist_parms%n_groups
          call getarg(1,agroup)
          read(agroup,*)i_group
          radars_per_group = float(n_radars_remap) / 
     1                       float(n_groups)
          write(6,*)' i_group / radars_per_group = ',i_group
     1                                              ,radars_per_group
          i_radar_start = 1 + nint(float(i_group-1) * radars_per_group)
          i_radar_end   =     nint(float(i_group  ) * radars_per_group)
          write(6,*)' Processing radar group ',i_radar_start,i_radar_end
      else
          write(6,*)' Processing all radars'
          i_radar_start = 1
          i_radar_end = n_radars_remap
      endif

      call get_grid_center(grid_cen_lat_orig,grid_cen_lon_orig,istatus)

      write(6,*)' REMAP > call get_laps_domain_95'
      call get_laps_domain_95(NX_L,NY_L,lat,lon,topo
     1                       ,dum_2d,grid_spacing_cen_m
     1                       ,istatus)
      if(istatus .ne. 1)then
          write(6,*)' ERROR return from get_laps_domain_95'
          goto 999
      endif

      do i_radar = i_radar_start,i_radar_end

!        Update grid namelist parameters in case domain has relocalized
         call get_directory('root',dataroot,len_root)
         call force_get_laps_config(dataroot,istatus)
         call get_grid_center(grid_cen_lat,grid_cen_lon,istatus)

!        Update lat/lon grid if grid location has changed
         if(grid_cen_lat .ne. grid_cen_lat_orig .OR.
     1      grid_cen_lon .ne. grid_cen_lon_orig      )then
            write(6,*)' REMAP > recall get_laps_domain_95'
            call get_laps_domain_95(NX_L,NY_L,lat,lon,topo
     1                             ,dum_2d,grid_spacing_cen_m
     1                             ,istatus)
            if(istatus .ne. 1)then
                write(6,*)' ERROR return from get_laps_domain_95'
                goto 999
            endif
          endif

          write(6,*)
          write(6,*)' Obtaining parameters for radar # ',i_radar
          call get_remap_parms(i_radar,n_radars_remap,max_times       ! I/O
     1                  ,path_to_radar                                ! O
     1                  ,laps_radar_ext,radar_subdir_dum              ! O
     1                  ,path_to_vrc                                  ! O
     1                  ,ref_min,min_ref_samples,min_vel_samples,dgr  ! O
     1                  ,namelist_parms,istatus)                      ! O
          if(istatus .ne. 1)then
              write(6,*)
     1            'Warning: bad status return from get_remap_parms'       
              go to 999
          endif

          ntimes_radar = 0

          do itimes = 1,max_times
              write(6,*)
              write(6,*)' Looping through radar/time # ',i_radar,itimes       
     1                 ,laps_radar_ext
              call remap_sub(i_radar,itimes
     1                      ,ntimes_radar                              ! I/O 
     1                      ,l_realtime                                ! O
     1                      ,laps_radar_ext       
     1                      ,radar_subdir_dum       
     1                      ,path_to_vrc,path_to_radar,ref_min
     1                      ,min_ref_samples,min_vel_samples,dgr
     1                      ,namelist_parms
     1                      ,NX_L,NY_L,NZ_L
     1                      ,lat,lon,topo                              ! I
     1                      ,istatus)       
              if(istatus .ne. 1)then
                  write(6,*)' remap: istatus returned from remap_sub = '
     1                                             ,istatus
                  if(l_realtime)then
!                     Note that missing tilts are given the same 0 value of
!                     istatus, and we assume that all radar times are experiencing
!                     the missing tilts
                      go to 900 ! saves computer time searching directories
                  else
                      write(6,*)' skip rest of times, go to next radar'
                      go to 900 
                  endif
              endif
          enddo ! itimes

 900  enddo ! i_radar

 999  end

      subroutine remap_sub(i_radar,itimes
     1                    ,ntimes_radar                                ! I/O
     1                    ,l_realtime                                  ! O
     1                    ,laps_radar_ext
     1                    ,c3_radar_subdir,path_to_vrc,path_to_radar
     1                    ,ref_min,min_ref_samples,min_vel_samples,dgr       
     1                    ,namelist_parms
     1                    ,NX_L,NY_L,NZ_L
     1                    ,lat,lon,topo                                ! I
     1                    ,istatus)

c
c     Velocity Obs
c
!     real grid_rvel(NX_L,NY_L,NZ_L) 
!     real grid_rvel_sq(NX_L,NY_L,NZ_L)
!     real grid_nyq(NX_L,NY_L,NZ_L)
!     integer ngrids_vel(NX_L,NY_L,NZ_L)
!     integer n_pot_vel(NX_L,NY_L,NZ_L)

      real, allocatable, dimension(:,:,:) :: grid_rvel ! Radial radar velocities
      real, allocatable, dimension(:,:,:) :: grid_rvel_sq
      real, allocatable, dimension(:,:,:) :: grid_nyq
      integer, allocatable, dimension(:,:,:) :: ngrids_vel
      integer, allocatable, dimension(:,:,:) :: n_pot_vel

c
c     Reflectivity Obs
c
!     real grid_ref (NX_L,NY_L,NZ_L)  !  Radar reflectivities
!     integer ngrids_ref (NX_L,NY_L,NZ_L)
!     integer n_pot_ref (NX_L,NY_L,NZ_L)

      real, allocatable, dimension(:,:,:) :: grid_ref
      integer, allocatable, dimension(:,:,:) :: ngrids_ref
      integer, allocatable, dimension(:,:,:) :: n_pot_ref

      include 'remap_dims.inc'
      include 'remap_constants.dat'      

      integer MAX_REF_TILT
      integer MAX_VEL_TILT

      parameter (MAX_REF_TILT = MAX_REF_GATES * MAX_RAY_TILT)
      parameter (MAX_VEL_TILT = MAX_VEL_GATES * MAX_RAY_TILT)

!     Variables used only in remap_process
      integer i_last_scan,i_first_scan
      integer i_tilt_proc_curr
      integer i4time_vol,i_num_finished_products

!     Variables used for data access and in fill_common 
   
      real b_ref(MAX_REF_TILT)
      real b_vel(MAX_VEL_TILT)
      real b_missing_data

      real v_nyquist_ray_a(MAX_RAY_TILT)
      real azim(MAX_RAY_TILT)
      real eleva 

!     Declarations for call to 'Read_data_88d'
      Real  Slant_ranges_m (max_gates),
     :        Elevation_deg,
     :        Az_array(MAX_RAY_TILT) 
      real, allocatable, dimension(:,:) :: Velocity
      real, allocatable, dimension(:,:) :: Reflect

      integer ref_index, vel_index, io_stat
      integer n_rays, i_scan, i_tilt, n_ref_gates, n_vel_gates

!     Radar Location variables

      integer i_lat,i_lon,i_alt
      real radar_lat
      real radar_lon
      real radar_alt
      character*4 rname_ptr
      character*4 laps_radar_ext
      character*3 c3_radar_subdir
      character*(*) path_to_vrc, path_to_radar

      Integer       ioffset                 
      Integer       joffset                 
      Logical       l_offset_radar

      real lat(NX_L,NY_L)      
      real lon(NX_L,NY_L)      
      real topo(NX_L,NY_L)     

!     Misc Local variables

      character string_time(9)
      character full_fname(91)
      integer initial_ray                ! flag for first ray in volume 
      integer alls_well, knt_bad_stat
      integer i_angle, past_angle
      integer past_scan, past_tilt
      integer len_fname
      integer write_and_exit
      integer i_vcp
      integer ref_ptr, vel_ptr

      integer VERBOSE

      logical l_realtime

!     Function call declarations
      integer get_field_num
      integer get_altitude
      integer get_latitude
      integer get_longitude
      integer get_fixed_angle
      integer get_scan
      integer get_tilt
      integer get_num_rays
      integer get_number_of_gates
      integer get_vcp
      integer get_azi
      integer get_data_field
      integer get_status
      real    get_nyquist

!     Beginning of Executable Code 
!     Some initializations  
      ISTAT = INIT_TIMER()
      VERBOSE = 1

      call get_r_missing_data(r_missing_data,istatus)

      b_missing_data = 255. ! flag value (corresponds to an unsigned integer)

      i_tilt_proc_curr = 1

      call radar_init(i_radar,path_to_radar,path_to_vrc,itimes         ! I
     1               ,b_missing_data                                   ! I
     1               ,i_tilt_proc_curr                                 ! I/O
     1               ,l_realtime                                       ! O
     1               ,i_last_scan,istatus)                             ! O
      write(6,*)' remap_sub: istatus from radar_init (1st call) ='
     1         ,istatus    
      if(istatus .ne. 1)then
          return
      endif

      call get_radarname(rname_ptr,istatus)
      if(istatus .ne. 1)then
          write(6,*)' remap_sub: ERROR returned from get_radarname'      
          return
      endif

      i_first_scan = 1

      i_alt=get_altitude() 
      i_lat=get_latitude() 
      i_lon=get_longitude() 

      radar_alt=  float(i_alt)
      radar_lat=  0.00001 * float(i_lat)
      radar_lon=  0.00001 * float(i_lon) 
      write(6,*)' Radar altitude (m): ',radar_alt  
      write(6,*)' Radar latitude (degrees): ',radar_lat  
      write(6,*)' Radar longitude (degrees): ',radar_lon  

      if(radar_alt .eq. 0. .OR. 
     1   radar_lat .eq. 0. .or. radar_lon .eq. 0. .OR.
     1   abs(radar_lat) .gt. 90.)then
          write(6,*)' ERROR: radar coords not initialized'
          istatus = 0
          return
      endif

      ntimes_radar = ntimes_radar + 1

      call get_grid_spacing_cen(grid_spacing_m, i_status)
      if(i_status .ne. 1)then
          write(6,*)' Error in obtaining grid_spacing_m'
          goto 900 ! return
      endif

      call get_l_offset_radar(nx_l,ny_l,grid_spacing_m,         ! I
     1                        nx_r,ny_r,igrid_r,l_offset_radar) ! O

!     Turn off l_offset for testing
!     nx_r = nx_l
!     ny_r = ny_l
!     igrid_r = 0
!     l_offset_radar = .false.

      call get_ij_offset_radars(nx_l,ny_l,1,                      ! I
     1                          igrid_r,l_offset_radar,           ! I   
     1                          lat,lon,radar_lat,radar_lon,      ! I
     1                          ioffset,joffset)                  ! O

      write(6,*)' Returned from get_ij_offset_radars:',ioffset,joffset

!     call lut_gen FORTRAN routine 
      if(ntimes_radar .eq. 1)then ! first time we have data from this radar
          call lut_gen(rname_ptr,radar_lat,radar_lon,radar_alt
     1                ,ioffset,joffset                            ! I
     1                ,NX_L,NY_L,NZ_L)                            ! I
      endif

      I4_elapsed = ishow_timer()

!     get data indices needed for other a2io library calls  

      ref_index=get_field_num('DBZ') 
      vel_index=get_field_num('VEL')
      write(6,*)
      write(6,*)'  Retrieved reflectivity index as ',ref_index  
      write(6,*)'  Retrieved velocity index as ',vel_index  

!     Misc initializations  

      initial_ray  = 1
      n_rays=0 
      write_and_exit = 0
      alls_well = 1

!     Allocate Velocity arrays
      allocate(grid_rvel(NX_R,NY_R,NZ_L),STAT=istat_alloc)      
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate grid_rvel'
          stop
      endif

      allocate(grid_rvel_sq(NX_R,NY_R,NZ_L),STAT=istat_alloc)      
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate grid_rvel_sq'
          stop
      endif

      allocate(grid_nyq(NX_R,NY_R,NZ_L),STAT=istat_alloc)      
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate grid_nyq'
          stop
      endif

      allocate(ngrids_vel(NX_R,NY_R,NZ_L),STAT=istat_alloc)      
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate ngrids_vel'
          stop
      endif

      allocate(n_pot_vel(NX_R,NY_R,NZ_L),STAT=istat_alloc)      
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate n_pot_vel'
          stop
      endif

!     Allocate Reflectivity arrays
      allocate(grid_ref(NX_R,NY_R,NZ_L),STAT=istat_alloc)      
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate grid_ref'
          stop
      endif

      allocate(ngrids_ref(NX_R,NY_R,NZ_L),STAT=istat_alloc)      
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate ngrids_ref'
          stop
      endif

      allocate(n_pot_ref(NX_R,NY_R,NZ_L),STAT=istat_alloc)      
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate n_pot_ref'
          stop
      endif

!     Begin infinite loop to continuously read radar data  

      do while(alls_well .eq. 1) 


!         Begin loop to fill buffer arrays with data from the circular buffer.
!         Call remap routines and reset pointers at the end of a volume scan  

! Test for existence of velocity data.
! Do we also need to test for reflectivity data?   

          if ( get_status(ref_index) .eq. 0 .or.
     1         get_status(vel_index) .eq. 0 ) then
            knt_bad_stat = 0 
            i_angle = get_fixed_angle() 
            if(i_angle .eq. i_missing_data)then
                write(6,*)' Warning: invalid i_angle in remap_sub'
                istatus = 0
                goto 900 ! return
            endif

            i_scan = get_scan() 
            i_tilt = get_tilt() 
            num_rays = get_num_rays() 
            if(num_rays .lt. 0 .or. num_rays .gt. 10000)then
                write(6,*)' Warning: num_rays out of bounds',num_rays
                istatus = 0
                goto 900 ! return
            endif

            if(n_rays .eq. n_rays/10 * 10)then
!             write(6,*)'  Good status received'
              write(6,*)'  i_angle, i_tilt = ', i_angle, i_tilt
            endif

            if ( initial_ray .eq. 1 ) then
              past_scan = i_scan 
              past_tilt = i_tilt 
              past_angle = i_angle 
              eleva = 0.01 * float(i_angle)

              call get_volume_time(i4time_vol)
              call make_fnam_lp (i4time_vol,string_time,istatus) 

              i_vcp=get_vcp() 
              write(6,*)'   VCP number for this volume: ',i_vcp

              if(VERBOSE .eq. 1)then
                write(6,*)' i4time_vol returned ',i4time_vol
              endif

              write(6,*)' Time is ',string_time
              initial_ray = 0

              ref_ptr = 1
              vel_ptr = 1

            endif ! initial_ray = 1

!           Test for end of tilt
            if( i_tilt .eq. past_tilt .and. i_scan .eq. past_scan 
     1                                .and. n_rays .lt. MAX_RAY_TILT 
     1                                .and. n_rays .lt. num_rays) then

!             Not end of tilt
              n_rays = n_rays + 1
              iarg = get_azi(n_rays)
              if(iarg .eq. i_missing_data)then
                  write(6,*)' Warning: invalid azi in remap_sub',n_rays       

!                 We set this to missing and continue to try remainder of tilt
                  azim(n_rays) = r_missing_data

              else
                  azim(n_rays) = 0.01 * iarg

              endif

              v_nyquist_ray_a(n_rays) = get_nyquist() 

              if(n_rays-1 .eq. n_rays/10 * 10)then
                write(6,*)'    n_rays = ',n_rays  
     1                   ,'    ref_ptr / vel_ptr = ' 
     1                   ,     ref_ptr,vel_ptr  
              endif

              if(VERBOSE .eq. 1)then
                n_ref_gates = get_number_of_gates(ref_index) 
                call get_first_gate(ref_index,first_ref_m,gsp_ref) 
                n_vel_gates = get_number_of_gates(vel_index) 
                call get_first_gate(vel_index,first_vel_m,gsp_vel) 
              endif

              io_stat = get_data_field(ref_index, b_ref(ref_ptr)
     1                                ,ref_ptr  , MAX_REF_GATES
     1                                , b_missing_data) 
              if(io_stat .ne. 1)then
                  write(6,*)' WARNING: incorrect ref data'
                  return
              endif

              io_stat = get_data_field(vel_index, b_vel(vel_ptr)
     1                                ,vel_ptr  , MAX_VEL_GATES
     1                                , b_missing_data) 
              if(io_stat .ne. 1)then
                  write(6,*)' WARNING: incorrect vel data'
                  return
              endif

              ref_ptr = ref_ptr + n_ref_gates 
              vel_ptr = vel_ptr + n_vel_gates 

            else ! end of tilt

              if( i_angle .lt. past_angle .or. i_scan .ne. past_scan )
     1             i_last_scan = 1

! call FILL_COMMON to fill up the common data area for current tilt   

              write(6,*)'  Calling fill_common, i_angle, past_angle ',
     1                                          i_angle, past_angle

              write(6,*)'  n_rays, past_tilt, b_missing_data ',
     1                     n_rays, past_tilt, b_missing_data

              n_ref_gates = get_number_of_gates(ref_index) 
              n_vel_gates = get_number_of_gates(vel_index) 

              if(VERBOSE .eq. 1)then
                write(6,*)'  n_ref_gates, n_vel_gates ',
     1                       n_ref_gates, n_vel_gates  

                write(6,*)'  first_ref_m, first_vel_m ',
     1                       first_ref_m, first_vel_m  

                write(6,*)'  gsp_ref,gsp_vel ',
     1                       gsp_ref,gsp_vel  

              endif

              call fill_common(
     1               b_ref,b_vel,n_rays,i_tilt,
     1               n_ref_gates,n_vel_gates,
     1               first_ref_m,gsp_ref,first_vel_m,gsp_vel,
     1               azim,v_nyquist_ray_a,eleva,b_missing_data) 

! call RADAR_INIT to assess the current tilt (first/last) and read next tilt
              i_tilt_proc_next = i_tilt_proc_curr + 1
              call radar_init(i_radar,path_to_radar,path_to_vrc,itimes   ! I
     1                       ,b_missing_data                             ! I
     1                       ,i_tilt_proc_next                           ! I/O
     1                       ,l_realtime,i_last_scan                     ! O
     1                       ,istatus)                                   ! O
              write(6,*)
     1        ' remap_sub: istatus from radar_init (2nd call) =',istatus    
              if(istatus .ne. 1)then
                  return
              endif

! call REMAP_PROCESS with the current tilt
              write(6,*)
     1        '  Calling remap_process: past_tilt/i_tilt_proc_curr',
     1                                  past_tilt,i_tilt_proc_curr

              if(min_ref_samples .eq. -1)then
                  if(grid_spacing_m .le. 5500.)then
                      min_ref_samples = 1
                  else
                      min_ref_samples = 4
                  endif
              endif

              if(min_vel_samples .eq. -1)then
                  if(grid_spacing_m .le. 3500.)then
                      min_vel_samples = 1
                  else
                      min_vel_samples = 4
                  endif
              endif

              write(6,*)' min_samples ref/vel = ',min_ref_samples
     1                                           ,min_vel_samples
              write(6,*)'  i_last, i_first ',i_last_scan,i_first_scan
              write(6,*)'  i4time_vol, i_num,  istatus',
     1                i4time_vol,i_num_finished_products,istatus

c             Get radar data from the storage area (formerly in remap_process)
c
              allocate(Velocity(max_gates,MAX_RAY_TILT)
     1                                   ,STAT=istat_alloc)      
              if(istat_alloc .ne. 0)then
                  write(6,*)' ERROR: Could not allocate Velocity'
                  stop
              endif

              allocate(Reflect(max_gates,MAX_RAY_TILT)
     1                                  ,STAT=istat_alloc)       
              if(istat_alloc .ne. 0)then
                  write(6,*)' ERROR: Could not allocate Reflect'
                  stop
              endif

              call Read_Data_88D(
     :               i_tilt_proc_curr,
     :               vel_thr_rtau,
     :               r_missing_data,       ! Input
     :               namelist_parms,
     :               gate_spacing_m_ret,   ! Output
     :               i_scan_mode,
     :               Num_sweeps,
     :               Elevation_deg,
     :               n_rays_88d,
     :               n_gates_88d,     ! Ref and Vel are on the same # of gates
     :               Slant_ranges_m,
     :               Velocity,
     :               Reflect,
     :               Az_Array,
     :               vel_nyquist,
     :               istatus_tilt)

!             QC check that elevation angle is within lookup table range
              if(Elevation_deg .lt. MIN_ELEV)then
                 write(6,*)
     1               ' Warning: Elevation_deg < MIN_ELEV in remap_sub'
     1                         ,Elevation_deg,MIN_ELEV
                 istatus = 0
                 goto 900 ! return
              endif

              call remap_process(
     1            i_tilt_proc_curr,i_last_scan,i_first_scan,             ! I
     :            grid_rvel,grid_rvel_sq,grid_nyq,ngrids_vel,n_pot_vel,  ! O
     :            grid_ref,ngrids_ref,n_pot_ref,                         ! O
     1            NX_L,NY_L,NZ_L,NX_R,NY_R,                              ! I
     1            l_offset_radar,ioffset,joffset,                        ! I
     1            lat,lon,topo,                                          ! I
     1            i_scan_mode,                                           ! I
     :            Slant_ranges_m,                                        ! I
     :            n_rays_88d,                                            ! I
     :            n_gates_88d,                                           ! I
     1            Velocity,Reflect,                                      ! I
     1            Az_Array,MAX_RAY_TILT,Elevation_deg,                   ! I
     1            vel_nyquist,                                           ! I
     1            ref_min,min_ref_samples,min_vel_samples,dgr,           ! I
     1            laps_radar_ext,c3_radar_subdir,path_to_vrc,            ! I
     1            namelist_parms,                                        ! I
     1            i4time_vol,                                            ! I
     1            i_num_finished_products,istatus_tilt,istatus)          ! O

              deallocate(Velocity)
              deallocate(Reflect)

              if(istatus .ne. 1)then
                  write(6,*)
     1            ' remap_sub: remap_process returned istatus =',istatus     
                  goto 900 ! return
              endif

              if(i_last_scan .eq. 1)then
                  write(6,*)' Volume completed, return from remap_sub'
                  istatus = 1
                  goto 900 ! return
              endif

              i_tilt_proc_curr = i_tilt_proc_next
              i_first_scan = 0

              if( i_angle .lt. past_angle .or. 
     1            i_scan  .ne. past_scan ) then 
                write(6,*)' Reset to beginning of volume'
                i_first_scan = 1
                i_tilt_proc_curr = 0  
                past_angle= i_angle 
              endif

              n_rays = 0 
              initial_ray = 1

            endif ! test for end of tilt 
          endif   ! close velocity status block 
      enddo       ! close infinite while loop    (increment tilt)

 900  continue

!     Deallocate arrays
      deallocate (grid_rvel)
      deallocate (grid_rvel_sq)
      deallocate (grid_nyq)
      deallocate (ngrids_vel)
      deallocate (n_pot_vel)

      deallocate (grid_ref)
      deallocate (ngrids_ref)
      deallocate (n_pot_ref)

      return
      end

 
       subroutine get_remap_parms(i_radar,n_radars_remap,max_times    ! I/O
     1            ,path_to_radar,laps_radar_ext                       ! O
     1            ,c3_radar_subdir,path_to_vrc                        ! O
     1            ,ref_min,min_ref_samples,min_vel_samples,dgr        ! O
     1            ,namelist_parms,istatus)                            ! O 

       include 'remap_constants.dat'      
       include 'radar_mosaic_dim.inc'      

       integer MAX_RADARS_REMAP
       parameter (MAX_RADARS_REMAP=max_radars_mosaic)

       character*150 path_to_radar_a(MAX_RADARS_REMAP),path_to_vrc_nl       
       character*(*) path_to_radar,path_to_vrc

!      character*4 c4_radarname_a(MAX_RADARS_REMAP)
!      character*4 c4_radarname

       character*4   laps_radar_ext_a(MAX_RADARS_REMAP)
       character*(*) laps_radar_ext

       character*3 c3_radar_subdir

       logical l_line_ref_qc,l_hybrid_first_gate,l_unfold,l_ppi_mode

       namelist /remap_nl/ n_radars_remap,max_times,path_to_radar_a ! ,c4_radarname_a
     1                    ,laps_radar_ext_a,path_to_vrc_nl
     1                    ,ref_min,min_ref_samples,min_vel_samples,dgr
     1                    ,abs_vel_min,l_line_ref_qc,l_hybrid_first_gate       
     1                    ,l_unfold,l_ppi_mode,n_groups
       character*150 static_dir,filename

!      Default values
       l_ppi_mode = .false. ! for backwards compatability with the namelist

       call get_directory('nest7grid',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/remap.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,remap_nl,err=901)
       close(1)

!      Fill namelist parms
       namelist_parms%abs_vel_min = abs_vel_min
       namelist_parms%l_line_ref_qc = l_line_ref_qc 
       namelist_parms%l_hybrid_first_gate = l_hybrid_first_gate
       namelist_parms%l_unfold = l_unfold
       namelist_parms%l_ppi_mode = l_ppi_mode
       namelist_parms%n_groups   = n_groups

       if(i_radar .eq. 0)then
           return
       endif

       if(i_radar        .gt. MAX_RADARS_REMAP  
     1                   .OR. 
     1    n_radars_remap .gt. MAX_RADARS_REMAP)then       
           write(6,*)' ERROR: too many radars in get_remap_parms'
     1              ,i_radar,n_radars_remap,MAX_RADARS_REMAP
           istatus = 0
           return
       endif

       path_to_radar  = path_to_radar_a(i_radar)
!      c4_radarname   = c4_radarname_a(i_radar)

       call s_len(laps_radar_ext_a(i_radar),len_ext)
       laps_radar_ext = laps_radar_ext_a(i_radar)(1:len_ext)

       length = min(len(path_to_vrc),len(path_to_vrc_nl))
       path_to_vrc    = path_to_vrc_nl(1:length)

!      Determine name of radar_subdir if any
       i_ext = 0
       do i = 1,i_radar
           if(laps_radar_ext_a(i) .eq. 'vrc')then
               i_ext = i_ext + 1
           endif
       enddo ! i

       write(c3_radar_subdir,801)i_ext
 801   format(i3.3)

!      write(6,*)' c4_radarname    = ',c4_radarname
       call s_len(laps_radar_ext,len_ext)
       write(6,*)' n_radars_remap  = ',n_radars_remap
       write(6,*)' max_times       = ',max_times
       write(6,*)' laps_radar_ext  = ',laps_radar_ext(1:len_ext)
       write(6,*)' c3_radar_subdir = ',c3_radar_subdir
       write(6,*)' path_to_vrc     = ',path_to_vrc
       write(6,*)' ref_min         = ',ref_min
       write(6,*)' min_ref_samples = ',min_ref_samples
       write(6,*)' min_vel_samples = ',min_vel_samples
       write(6,*)' dgr             = ',dgr
       write(6,*)' abs_vel_min     = ',abs_vel_min
       write(6,*)' l_line_ref_qc   = ',l_line_ref_qc
       write(6,*)' l_hybrid_first_gate = ',l_hybrid_first_gate
       write(6,*)' l_unfold        = ',l_unfold
       write(6,*)' l_ppi_mode      = ',l_ppi_mode

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading remap_nl in ',filename
       write(*,remap_nl)
       istatus = 0
       return

       end

