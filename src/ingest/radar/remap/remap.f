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

      character path_to_radar*150, ext_dum*3
     1         ,radar_subdir_dum*3, path_to_vrc*15
       
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

!     This first call returns only 'n_radars_remap'
      call get_remap_parms(0,n_radars_remap,path_to_radar
     1       ,ext_dum,radar_subdir_dum,path_to_vrc,istatus)       
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

      max_times = 3

      do i_radar = 1,n_radars_remap
          call get_remap_parms(i_radar,n_radars_remap,path_to_radar
     1                  ,ext_dum,radar_subdir_dum
     1                  ,path_to_vrc,istatus)       
          if(istatus .ne. 1)then
              write(6,*)
     1            'Warning: bad status return from get_remap_parms'       
              go to 999
          endif

          do itimes = 1,max_times
              write(6,*)
              write(6,*)' Looping through radar/time # ',i_radar,itimes       
              call remap_sub(i_radar,itimes,ext_dum,radar_subdir_dum       
     1                      ,path_to_vrc,NX_L,NY_L,NZ_L,istatus)
          enddo ! itimes

      enddo ! i_radar

 999  end

      subroutine remap_sub(i_radar,itimes,laps_radar_ext
     1                    ,c3_radar_subdir,path_to_vrc,NX_L,NY_L,NZ_L
     1                    ,istatus)

      include 'remap.inc'
      include 'remap_dims.inc'

      integer MAX_REF_TILT
      integer MAX_VEL_TILT

      parameter (MAX_REF_TILT = MAX_REF_GATES * MAX_RAY_TILT)
      parameter (MAX_VEL_TILT = MAX_VEL_GATES * MAX_RAY_TILT)

!     Variables used only in remap_process
      integer i_last_scan,i_first_scan
      integer i_tilt_proc_curr
      integer i4time_vol,i_num_finished_products,i_status

!     Variables used for data access and in fill_common 
   
      real*4 b_ref(MAX_REF_TILT)
      real*4 b_vel(MAX_VEL_TILT)
      real*4 b_missing_data

      real v_nyquist_ray_a(MAX_RAY_TILT)
      real azim(MAX_RAY_TILT)
      real eleva 

      integer ref_index, vel_index, io_stat
      integer n_rays, i_scan, i_tilt, n_ref_gates, n_vel_gates

!     Radar Location variables

      integer i_lat,i_lon,i_alt
      real radar_lat
      real radar_lon
      real radar_alt
      character*4 rname_ptr
      character*3 laps_radar_ext, c3_radar_subdir
      character*(*) path_to_vrc

!     Misc Local variables

      character string_time(9)
      character full_fname(91)
      integer initial_ray                ! flag for first ray in volume 
      integer alls_well, read_next, knt_bad_stat
      integer i_angle, past_angle
      integer past_scan, past_tilt
      integer len_fname
      integer write_and_exit
      integer i_vcp
      integer ref_ptr, vel_ptr

      integer VERBOSE

!     Function call declarations
      integer get_field_num
      integer get_altitude
      integer get_latitude
      integer get_longitude
      integer get_fixed_angle
      integer get_scan
      integer get_tilt
      integer get_num_rays
      integer get_vcp
      integer get_azi
      integer get_data_field
      integer get_status
      real    get_nyquist

!     Beginning of Executable Code 
!     Some initializations  
      ISTAT = INIT_TIMER()
      VERBOSE = 1

      n_vel_gates = MAX_VEL_GATES
      n_ref_gates = MAX_REF_GATES
      b_missing_data = 255.

      i_tilt_proc_curr = 1

      call radar_init(i_radar,i_tilt_proc_curr,i_last_scan,istatus)
      if(istatus .ne. 1)then
          write(6,*)' remap_sub: istatus returned from radar_init ='
     1              ,istatus      
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

      if(radar_alt .eq. 0. .and. 
     1   radar_lat .eq. 0. .and. radar_lon .eq. 0.)then
          write(6,*)' ERROR: radar coords not initialized'
          istatus = 0
          return
      endif

!     call lut_gen FORTRAN routine 
      if(itimes .eq. 1)then
          call lut_gen(rname_ptr,radar_lat,radar_lon,radar_alt
     1                ,NX_L,NY_L,NZ_L)
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
      read_next = 1
      alls_well = 1

!     Begin infinite loop to continuously read radar data  

      do while(alls_well .eq. 1) 


!       Begin loop to fill buffer arrays with data from the circular buffer.
!       Call remap routines and reset pointers at the end of a volume scan  

!       if (read_next .eq. 1) then
        if (.false.) then
          if(VERBOSE .eq. 1)then
            write(6,*)'  Calling read_radial '
          endif

          io_stat=read_radial() 
          if(io_stat .eq. 1) then
            write(6,*)'  Read_radial returned double eof '
            write_and_exit = 1
          endif

          if(VERBOSE .eq. 1)then
            write(6,*)'  Back from read_radial '
          endif

        else
          read_next = 1

! Test for existence of velocity data.
! Do we also need to test for reflectivity data?   

          if ( get_status(ref_index) .eq. 0 .or.
     1         get_status(vel_index) .eq. 0 ) then
            knt_bad_stat = 0 
            i_angle = get_fixed_angle() 
            i_scan = get_scan() 
            i_tilt = get_tilt() 
            num_rays = get_num_rays() 

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
              call make_fnam_lp (i4time_vol,string_time,i_status) 

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
              azim(n_rays) = 0.01 * get_azi(n_rays)
              v_nyquist_ray_a(n_rays) = get_nyquist() 

              if(n_rays-1 .eq. n_rays/10 * 10)then
                write(6,*)'    n_rays = ',n_rays  
     1                   ,'    ref_ptr / vel_ptr = ' 
     1                   ,     ref_ptr,vel_ptr  
              endif

!             if ( (n_rays-1) % 60 .eq. 0)
!               write(6,*)'  eleva = %f  azim = %f  Nyqst = %f\n",
!                   eleva,azim(n_rays-1),v_nyquist_ray_a(n_rays-1)) 

              if(VERBOSE .eq. 1)then
                ng_ref = get_number_of_gates(ref_index) 
                call get_first_gate(ref_index,first_ref_m,gsp_ref) 
                ng_vel = get_number_of_gates(vel_index) 
                call get_first_gate(vel_index,first_vel_m,gsp_vel) 
              endif

              io_stat = get_data_field(ref_index, b_ref(ref_ptr)
     1                                ,ref_ptr  , MAX_REF_GATES
     1                                , b_missing_data) 
              io_stat = get_data_field(vel_index, b_vel(vel_ptr)
     1                                ,vel_ptr  , MAX_VEL_GATES
     1                                , b_missing_data) 

!             parameter (MAX_REF_TILT = MAX_REF_GATES * MAX_RAY_TILT)
!             parameter (MAX_VEL_TILT = MAX_VEL_GATES * MAX_RAY_TILT)

              ref_ptr = ref_ptr + MAX_REF_GATES 
              vel_ptr = vel_ptr + MAX_VEL_GATES 

            else ! end of tilt

              if( i_angle .lt. past_angle .or. i_scan .ne. past_scan )
     1             i_last_scan = 1

! call FILL_COMMON to fill up the common data area for current tilt   

              write(6,*)'  Calling fill_common, i_angle, past_angle ',
     1                                          i_angle, past_angle

              write(6,*)'  n_rays, past_tilt, b_missing_data ',
     1                     n_rays, past_tilt, b_missing_data

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
              call radar_init(i_radar,i_tilt_proc_next,i_last_scan
     1                                                ,istatus)

! call REMAP_PROCESS with the current tilt
              write(6,*)'  Calling remap_process past_tilt '
     1                                         , past_tilt  
              write(6,*)'  Calling remap_process i_tilt_proc_curr '
     1                                         , i_tilt_proc_curr  

              write(6,*)'  i_last, i_first ',i_last_scan,i_first_scan
              write(6,*)'  i4time_vol, i_num,  i_status',
     1                i4time_vol,i_num_finished_products,i_status

              call remap_process(
     1            i_tilt_proc_curr,i_last_scan,i_first_scan,
     :            grid_rvel,grid_rvel_sq,grid_nyq,ngrids_vel,n_pot_vel,
     :            grid_ref,ngrids_ref,n_pot_ref,
     1            NX_L,NY_L,NZ_L,
     1            laps_radar_ext,c3_radar_subdir,path_to_vrc,
     1            i4time_vol,full_fname,len_fname,
     1            i_num_finished_products,i_status) 


              if(i_last_scan .eq. 1)then
                  write(6,*)' Volume completed, return from remap_sub'
                  istatus = 1
                  return
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
              read_next = 0

            endif ! test for end of tilt 
          endif   ! close velocity status block 
        endif     ! close read_next block
      enddo       ! close infinite while loop    (increment tilt)

      end

 
       subroutine get_remap_parms(i_radar,n_radars_remap
     1            ,path_to_radar,laps_radar_ext
     1            ,c3_radar_subdir,path_to_vrc,istatus) 

       include 'radar_mosaic_dim.inc'      

       integer*4 MAX_RADARS_REMAP
       parameter (MAX_RADARS_REMAP=max_radars_mosaic)

       character*150 path_to_radar_a(MAX_RADARS_REMAP),path_to_vrc_nl       
       character*(*) path_to_radar,path_to_vrc

!      character*4 c4_radarname_a(MAX_RADARS_REMAP)
!      character*4 c4_radarname

       character*4   laps_radar_ext_a(MAX_RADARS_REMAP)
       character*(*) laps_radar_ext

       character*3 c3_radar_subdir

       namelist /remap_nl/ n_radars_remap,path_to_radar_a     ! ,c4_radarname_a
     1                    ,laps_radar_ext_a,path_to_vrc_nl      
       character*150 static_dir,filename
 
       call get_directory('nest7grid',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/remap.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,remap_nl,err=901)
       close(1)

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
       write(6,*)' laps_radar_ext  = ',laps_radar_ext(1:len_ext)
       write(6,*)' c3_radar_subdir = ',c3_radar_subdir
       write(6,*)' path_to_vrc     = ',path_to_vrc

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

