
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

      call remap_sub(NX_L,NY_L,NZ_L,istatus)

 999  end

      subroutine remap_sub(NX_L,NY_L,NZ_L,istatus)

      include 'remap.inc'
      include 'remap_dims.inc'

      integer MAX_REF_TILT
      integer MAX_VEL_TILT

      parameter (MAX_REF_TILT = MAX_REF_GATES * MAX_RAY_TILT)
      parameter (MAX_VEL_TILT = MAX_VEL_GATES * MAX_RAY_TILT)

!     Variables used only in remap_process
      integer i_last_scan,i_first_scan
      integer i_tilt_proc
      integer i4time_vol,i_num_finished_products,i_status

!     Variables used for data access and in fill_common 
   
      real*4 b_ref(MAX_REF_TILT)
      real*4 b_vel(MAX_VEL_TILT)
      real*4 b_missing_data, i4_to_byte

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

!     Misc Local variables

      character sw
      integer iyr, iday, imon, ihour, imin, isec ! time variables 
      character source(80)
      character string_time(9)
      character full_fname(91)
      integer initial_ray                ! flag for first ray in volume 
      integer alls_well, read_next, knt_bad_stat
      integer i_angle, past_angle
      integer past_scan, past_tilt
      integer len_fname
      integer compr_on, xmit_on, write_and_exit
      integer i_mode, i
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
      integer get_azi
      integer get_nyquist
      integer get_data_field

!     Beginning of Executable Code 
!     Some initializations  
      ISTAT = INIT_TIMER()
      VERBOSE = 1

      n_vel_gates = MAX_VEL_GATES
      n_ref_gates = MAX_REF_GATES
      b_missing_data = 255.

!     strtm_ptr = string_time
!     fnm_ptr = full_fname
!     ref0_ptr= b_ref
!     vel0_ptr= b_vel
!     azi0_ptr= azim
!     nyq0_ptr= v_nyquist_ray_a

!     radar_lat=&radar_lat 
!     radar_lon=&radar_lon 
!     radar_alt=&radar_alt 


!     Get Radar name environment variable 
      call getenv('RADARNAME',rname_ptr) 

      if (rname_ptr .eq. ' ')then
        write(6,*)' Could not evaluate RADARNAME environment variable.'
        write(6,*)' Set the 4-character radar name using:'
        write(6,*)'  setenv RADARNAME Kxxx'
        write(6,*)'     before running remap.'
        call exit(1) 
      endif

!     call Archive II initialization routine  
      i_tilt_proc = 1

      call radar_init(i_tilt_proc,i_last_scan)

      i_first_scan = 1

      i_alt=389. 
      i_alt=get_altitude() 
      i_lat=get_latitude() 
      i_lon=get_longitude() 

      radar_alt=  float(i_alt)
      radar_lat=  0.00001 * float(i_lat)
      radar_lon= -0.00001 * float(i_lon) 
      write(6,*)' Radar altitude (m): ',radar_alt  
      write(6,*)' Radar latitude (degrees): ',radar_lat  
      write(6,*)' Radar longitude (degrees): ',radar_lon  

!     call lut_gen FORTRAN routine 
      call lut_gen(rname_ptr,radar_lat,radar_lon,radar_alt,NX_L,NY_L)        

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

          if ( get_status(ref_index) .eq. GOOD_STATUS .or.
     1         get_status(vel_index) .eq. GOOD_STATUS ) then
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
              v_nyquist_ray_a(n_rays) = 0.01 * get_nyquist() 

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
                gsp_ref = get_first_gate(ref_index) 
                ng_vel = get_number_of_gates(vel_index) 
                gsp_vel = get_first_gate(vel_index) 
!               write(6,*)'  ref: Number of gates = %i,  first gate = %i\n",
!                     ng_ref,gsp_ref) 
!               write(6,*)'  vel: Number of gates = %i,  first gate = %i\n",
!                     ng_vel,gsp_vel) 
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

! call the FORTRAN routine to fill up the common data area   

              write(6,*)'  Calling fill_common, i_angle, past_angle ',
     1                                          i_angle, past_angle

              write(6,*)'  n_rays, past_tilt, b_missing_data ',
     1                     n_rays, past_tilt, b_missing_data

              if(VERBOSE .eq. 1)then
                write(6,*)'  n_ref_gates, n_vel_gates ',
     1                       n_ref_gates, n_vel_gates  
              endif

              call fill_common(
     1               b_ref,b_vel,n_rays,i_tilt,
     1               n_ref_gates,n_vel_gates,
     1               azim,v_nyquist_ray_a,eleva,b_missing_data) 

! call the FORTRAN remapper module   
!             Read next tilt
              i_tilt_proc_new = i_tilt_proc + 1
              call radar_init(i_tilt_proc_new,i_last_scan)

              write(6,*)'  Calling remap_process past_tilt '
     1                                         , past_tilt  
              write(6,*)'  Calling remap_process i_tilt_proc '
     1                                         , i_tilt_proc  

              write(6,*)'  i_last, i_first ',i_last_scan,i_first_scan
              write(6,*)'  i4time_vol, i_num,  i_status',
     1                i4time_vol,i_num_finished_products,i_status

              call remap_process(
     1            i_tilt_proc,i_last_scan,i_first_scan,
     :            grid_rvel,grid_rvel_sq,grid_nyq,ngrids_vel,n_pot_vel,
     :            grid_ref,ngrids_ref,n_pot_ref,
     1            NX_L,NY_L,NZ_L,
     1            i4time_vol,full_fname,len_fname,
     1            i_num_finished_products,i_status) 


              if(i_last_scan .eq. 1)then
                  write(6,*)' Volume completed, exit program'
                  call exit(0)
              endif

              i_tilt_proc = i_tilt_proc_new
              i_first_scan = 0

              if( i_angle .lt. past_angle .or. 
     1            i_scan  .ne. past_scan ) then 
                write(6,*)' Reset to beginning of volume'
                i_first_scan = 1
                i_tilt_proc = 0  
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
