
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

      integer NUM_REF_GATES
      integer NUM_VEL_GATES
      integer N_REF_TILT
      integer N_VEL_TILT
      integer N_RAY_TILT
      integer NSIZE

      parameter (NUM_REF_GATES = 460)
      parameter (NUM_VEL_GATES = 920)
      parameter (N_REF_TILT = 174800)
      parameter (N_VEL_TILT = 349600)
      parameter (N_RAY_TILT = 380)
      parameter (NSIZE = 920)

!     Variables used only in remap_process
      integer i_last_scan,i_first_scan
      integer i_tilt_proc
      integer i4time_vol,i_num_finished_products,i_status

!     Variables used for data access and in fill_common 
   
      real b_ref(N_REF_TILT)
      real b_vel(N_VEL_TILT)

      real v_nyquist_ray_a(N_RAY_TILT)
      real azim(N_RAY_TILT)
      real eleva 

      integer ref_index, vel_index, io_stat
      integer n_rays, i_scan, i_tilt, n_ref_gates, n_vel_gates
      real b_missing_data

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

      integer VERBOSE

!     Beginning of Executable Code 
!     Some initializations  
      VERBOSE = 1

      n_vel_gates = NUM_VEL_GATES
      n_ref_gates = NUM_REF_GATES
      b_missing_data = MISSING

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

      call radar_init() 

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

!     get data indices needed for other a2io library calls  

      ref_index=get_field_num('DBZ') 
      write(6,*)'  Retrieved reflectivity index as ',ref_index  
      vel_index=get_field_num('VEL')
      write(6,*)'  Retrieved velocity index as ',vel_index  

!     Misc initializations  

      initial_ray  = 1
      n_rays=0 
      i_first_scan = 1
      i_tilt_proc = 0
      i_last_scan = 0
      write_and_exit = 0
      read_next = 1
      alls_well = 1

!     Begin infinite loop to continuously read radar data  

      do while(alls_well .eq. 1) 

!       Begin loop to fill buffer arrays with data from the circular buffer.
!       Call remap routines and reset pointers at the end of a volume scan  

        if (read_next .eq. 1) then
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

            if(VERBOSE .eq. 1)then
              write(6,*)'  Good status received'
              write(6,*)'  i_angle, i_tilt = ', i_angle, i_tilt
            endif

            if ( initial_ray .eq. 1 ) then
              past_scan = i_scan 
              past_tilt = i_tilt 
              past_angle = i_angle 
              eleva = 0.01 * float(i_angle)

              iyr = get_year() 
              imon = get_month() 
              iday = get_day() 
              ihour = get_hour() 
              imin = get_min() 
              isec = get_sec() 
              i4time_vol = int_to_i4time(iyr,imon,iday,ihour,imin,isec) 
              call make_fnam_lp (i4time_vol,string_time,i_status) 

              i_vcp=get_vcp() 
              write(6,*)'   VCP number for this volume: ',i_vcp

              if(VERBOSE .eq. 1)then
                write(6,*)'   iyr, imon, iday ',iyr,imon,iday
                write(6,*)' i4time_vol returned ',i4time_vol
              endif

              write(6,*)' ihour, imin, isec ',ihour,imin,isec
              write(6,*)' Time is ',string_time
              initial_ray = 0
            endif ! initial_ray = 1

            if( i_tilt .eq. past_tilt .and. i_scan .eq. past_scan .and.
     1          n_rays .lt. N_RAY_TILT ) then

              n_rays = n_rays + 1
              azim(n_rays-1) = 0.01 * get_azi()
              v_nyquist_ray_a(n_rays-1) = 0.01 * get_nyquist() 

!             if(VERBOSE .eq. 1)then
!               write(6,*)'    INFO FOR n_rays = %i \n",n_rays  
!               write(6,*)'    ref_ptr = %i    vel_ptr = %i\n",
!                              ref_ptr,vel_ptr) 
!             endif

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

              io_stat = get_data_field(ref_index, ref_ptr, n_ref_gates) 
              io_stat = get_data_field(vel_index, vel_ptr, n_vel_gates) 

!             ref_ptr += NUM_REF_GATES 
!             vel_ptr += NUM_VEL_GATES 

            else

              if( i_angle .lt. past_angle .or. i_scan .ne. past_scan )
     1             i_last_scan = 1

! call the FORTRAN routine to fill up the common data area   

              write(6,*)'  Calling fill_common, i_angle, past_angle ',
     1                                          i_angle,past_angle

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

              i_tilt_proc = i_tilt_proc + 1

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

              i_last_scan = 0
              i_first_scan = 0

              if( i_angle .lt. past_angle .or. 
     1            i_scan  .ne. past_scan ) then
                i_first_scan = 1
                i_tilt_proc = 0  
                past_angle= i_angle 
              endif

              n_rays = 0 
              initial_ray = 1
              read_next = 0

            endif

!         For bad status, increment bad status counter and try again.  

          elseif( knt_bad_stat .lt. MAX_BAD_STAT .and. 
     1            write_and_exit .eq. 0 ) then

            if(VERBOSE .eq. 1)then
              write(6,*)'  Bad status received for data'
            endif

            knt_bad_stat = knt_bad_stat + 1

!           Once 1000 consecutive bad stati have been received, assume end
!           of data and dump what might be in the buffer.  

          else
            write(6,*)knt_bad_stat, 'bad read status reports received'
                   
            knt_bad_stat = 0 

            if (n_rays .gt. 0) then

              write(6,*)' Transferring available radials ', n_rays  

              if( i_angle .lt. past_angle )i_last_scan = 1

!             call the FORTRAN routine to fill up the common data area   
              call fill_common(
     1             b_ref,b_vel,n_rays,past_tilt,
     1             n_ref_gates,n_vel_gates,
     1             azim,v_nyquist_ray_a,eleva,b_missing_data) 

!             call the FORTRAN remapper module   

              i_tilt_proc = i_tilt_proc + 1

              write(6,*)' Calling remap_process past_tilt = ', past_tilt
              write(6,*)' Calling remap_process i_tilt_proc = '
     1                                                  , i_tilt_proc  
              write(6,*)' i_last, i_first',i_last_scan,i_first_scan
              write(6,*)' i4time_vol, i_num, i_status',
     1              i4time_vol,i_num_finished_products,i_status  

              call remap_process(
     1            i_tilt_proc,i_last_scan,i_first_scan,
     :            grid_rvel,grid_rvel_sq,grid_nyq,ngrids_vel,n_pot_vel,
     :            grid_ref,ngrids_ref,n_pot_ref,
     1            NX_L,NY_L,NZ_L,
     1            i4time_vol,full_fname,len_fname,
     1            i_num_finished_products,i_status) 

              if(write_and_exit .eq. 1) call exit(0) 

              i_last_scan = 0
              i_first_scan = 0

              if( i_angle .lt. past_angle ) then
                i_first_scan = 1
                i_tilt_proc = 0 
                past_angle = i_angle 
              endif

              n_rays = 0 
              initial_ray = 1

            endif ! close n_rays .gt. 0 block  
          endif   ! close velocity status block   
        endif     ! close read_next block
      enddo       ! close infinite while loop  

      end
