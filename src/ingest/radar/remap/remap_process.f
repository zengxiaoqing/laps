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
      Subroutine Remap_process(
     :         i_tilt,                                     ! Integer (input)
     :         i_last_scan,                                ! Integer (input)
     :         i_first_scan,                               ! Integer (input)
     :         grid_rvel,grid_rvel_sq,grid_nyq,ngrids_vel,n_pot_vel, ! (output)
     :         grid_ref,ngrids_ref,n_pot_ref,                        ! (output)
     :         NX_L,NY_L,NZ_L,NX_R,NY_R,                   ! Integer   (input)
     :         l_offset_radar,ioffset,joffset,             !           (input)
     1         lat,lon,topo,                               !           (input)
     1         i_scan_mode,                                !           (input)
     :         Slant_ranges_m,                             !           (input)
     :         n_rays,                                     !           (input)
     :         n_gates,                                    !           (input)
     1         Velocity,Reflect,                           !           (input)
     1         Az_Array,MAX_RAY_TILT,Elevation_deg,        !           (input)
     1         vel_nyquist,                                !           (input)
     :         ref_min,min_ref_samples,min_vel_samples,dgr,! Integer (input)
     :         laps_radar_ext,c3_radar_subdir,             ! Char      (input)
     :         path_to_vrc,                                ! Char      (input)
     :         namelist_parms,                             ! Struct    (input)
     :         i_product_i4time,                           ! Integer (input)
     :         i_num_finished_products,                    ! Integer (output)
     :         i_status_tilt,i_status)                     ! Integer (output)
c
c     Subroutine remap_process
c
c     PURPOSE:
c       Main process routine for the REMAPPING algorithm
c
c **************************** History Section ****************************
c
c       Windsor, C. R.  10-JUL-1985     Original version
c       Albers, Steve    7-APR-1986     Update for version 03
c       Albers, Steve      JUN-1987     Map velocities to LAPS grid
c       Albers, Steve      MAR-1988     Streamlined and converted to 87 data
c       Albers, Steve      DEC-1988     FURTHER conversions for RT87 cartesian
c       Albers, Steve      MAY-1992     Turn off range unfolding for velocities
c       Albers, Steve      FEB-1993     MIN#, 40% FRAC QC for Reflectivity added
c       Albers, Steve      MAY-1994     88D Version for SUN RISC BOX
c       Brewster, Keith    AUG-1994     Clean-out of artifacts
c       Brewster, Keith    APR-1995     Added INITIAL_GATE parameter
c                                       Modified volume nyquist determination.
c       Brewster, Keith    SEP-1995     Added point-by-point Nyquist calc.
c       Albers, Steve   19-DEC-1995     Changed gate_spacing_m variable to 
c                                       gate_spacing_m_ret to prevent 
c                                       reassigning a parameter value. The 
c                                       location is the call to read_data_88D 
c                                       and a new declaration.
c                                       Environment variable evaluations added
c                                       for FTPing and purging the output.
c                                       New streamlined purging function.
c       Albers, Steve      FEB-1996     Linear reflectivity averaging (via lut)
c       Albers, Steve      MAY-1996     New igate_lut to reduce processing 
c       Albers, Steve          1998     More flexibility added

*********************** Declaration Section **************************
c
      include 'trigd.inc'
      implicit none
c
c     Input variables
c
      integer i_tilt, i_status_tilt
      integer i_last_scan
      integer i_first_scan
      integer i_product_i4time
      integer MAX_RAY_TILT
      integer NX_L,NY_L,NZ_L,NX_R,NY_R
      integer ioffset,joffset,io,jo,iomin,iomax,jomin,jomax

      integer min_ref_samples,min_vel_samples

      real ref_min,dgr
c
c     LAPS Grid Dimensions
c
      include 'remap_constants.dat'
      include 'remap.cmn'
c
c     Velocity Obs
c
      real grid_rvel(NX_R,NY_R,NZ_L)  !  Radial radar velocities
      real grid_rvel_sq(NX_R,NY_R,NZ_L)
      real grid_nyq(NX_R,NY_R,NZ_L)
      integer ngrids_vel(NX_R,NY_R,NZ_L)
      integer n_pot_vel(NX_R,NY_R,NZ_L)
c
c     Reflectivity Obs
c
      real grid_ref (NX_R,NY_R,NZ_L)  !  Radar reflectivities
      integer ngrids_ref (NX_R,NY_R,NZ_L)
      integer n_pot_ref (NX_R,NY_R,NZ_L)
c
c     Output variables
c
      integer i_num_finished_products
      integer i_status
c
c     Processing parameters
c
      real re43
      parameter (re43 = 8503700.) ! 4/3 radius of the earth in meters
c
c
      integer max_fields
      parameter (max_fields = 10)
c
c     Variables for NetCDF I/O
c
      character*150 dir

      character*9 a9time
      character*31 ext,ext_in
      character*3 var_a(max_fields)
      character*125 comment_a(max_fields)
      character*10  units_a(max_fields)
      character*4 laps_radar_ext
      character*3 c3_radar_subdir
      character*(*) path_to_vrc

c
c     Functions
c
      real height_to_zcoord
      real height_of_level
c
c     Misc local variables
c
      integer igate,i_scan_mode,jray,end_ext,ilut_ref
c
      Real  Slant_ranges_m (max_gates),
     :        Elevation_deg,
     :        Az_array(MAX_RAY_TILT),
     :        Velocity(max_gates,MAX_RAY_TILT),
     :        Reflect(max_gates,MAX_RAY_TILT)

      real, allocatable, dimension(:,:,:,:) :: out_array_4d

      real r_missing_data

      real lat(NX_L,NY_L)      
      real lon(NX_L,NY_L)      
      real topo(NX_L,NY_L)     
      real dum_2d(NX_L,NY_L)   ! Local
      integer k_eff(NX_L,NY_L) ! Local
c
      logical l_unfold, l_compress_output, l_domain_read, l_offset_radar       
      save l_domain_read
      data l_domain_read /.false./
c 
      real avgvel,vel_nyquist,vel_value,ref_value,lat_dum,lon_dum
      real v_nyquist_tilt(max_tilts)
      real v_nyquist_vol
      real gate_spacing_m_ret,grid_spacing_cen_m
      real height_grid,range_dum,range_new,azimuth,elevation_dum
      real height_guess
c
      integer i,j,k,k_low,ielev,igate_lut,iter,klut,idebug
      integer nazi,iran
      integer num_sweeps,n_rays,n_gates,n_obs_vel,n_output_data,nf
      integer igate_max
      integer igate_interval
      integer n_vel_grids_final,n_vel_grids_prelim
      integer n_ref_grids,n_ref_grids_qc_fail,nycor
      integer istatus,istatus_qc,istat_alloc,len_ext
      integer ishow_timer,i4_elapsed
      integer i_purge
      integer init_ref_gate_hyb,init_ref_gate_actual
      integer mingate_valid_ref,maxgate_valid_ref

      real rvel,azimuth_interval
      real rmax,height_max,rlat_radar,rlon_radar,rheight_radar
      real vknt,rknt,variance,hybrid_range

      character*4 c4_radarname
      character*7 c7_laps_xmit
      character*7 c7_laps_purge
      character*7 c7_laps_sleep

      save height_max, k_low, n_obs_vel, n_output_data
c
c     Beginning of executable code
c
      idebug = 0 ! more verbose output (0-2 range)

      write(6,*)
      write(6,805) i_first_scan,i_last_scan,i_tilt
  805 format(' REMAP_PROCESS > ifirst,ilast,tilt',4i5)

      call get_l_compress_radar(l_compress_output,i_status)
      if(i_status .ne. 1)then
          write(6,*)' Error in obtaining l_compress_radar'
          return
      endif

      rlat_radar = rlat_radar_cmn
      rlon_radar = rlon_radar_cmn
      rheight_radar = rheight_radar_cmn
      c4_radarname = c4_radarname_cmn

      if(rlat_radar .eq. 0.)then
          write(6,*)' ERROR: no radar coords in remap_process'
          i_status = 0
          return
      endif

      i_num_finished_products = 0

      call get_r_missing_data(r_missing_data, i_status)
      if(i_status .ne. 1)then
          write(6,*)' Error in obtaining r_missing_data'
          return
      endif
c
c     For first scan, initialize sums and counters to zero.
c
      IF (i_first_scan .eq. 1 .or. i_first_scan .eq. 999) THEN

        I4_elapsed = ishow_timer()

        write(6,806)
  806   format
     1  (' REMAP_PROCESS > 1st sweep - Initializing vel/ref arrays')

        n_obs_vel = 0
        n_output_data = 0

        grid_rvel(:,:,:) = 0.
        grid_rvel_sq(:,:,:) = 0.
        grid_nyq(:,:,:) = 0.
        grid_ref(:,:,:) = 0.
        ngrids_vel(:,:,:) = 0
        ngrids_ref(:,:,:) = 0
        n_pot_vel(:,:,:) = 0
        n_pot_ref(:,:,:) = 0
c
c       Compute maximum height of data needed.
c
        height_max = height_of_level(NZ_L)
c
c       Define Lower Limit of Radar Coverage in LAPS grid
c
        k_low = int(height_to_zcoord(rheight_radar,i_status))
        k_low = max(k_low,1)

!       True will map one tilt in each LAPS level (for testing only)
        if(namelist_parms%l_ppi_mode)then 
            k_low = 1 
        endif

      END IF ! initialize for 1st scan

c     Former location of 'read_data_88d' call
      IF (i_status_tilt .ne. 1) GO TO 998 ! abnormal return
c
      v_nyquist_tilt(i_tilt) = vel_nyquist
c
      write(6,*)' REMAP_PROCESS > vel_nyquist for this tilt = '
     :        ,i_tilt,vel_nyquist
c
c     Compute max range from elevation angle
c
      rmax = -re43 * sind(elevation_deg)
     :  + sqrt(re43*re43*sind(elevation_deg)**2 +
     :          height_max * (2.*re43 + height_max))


      print *, ' rmax,height_max= ',rmax,height_max

      print *, ' gate_spacing_m,gate_interval= ',gate_spacing_m,
     :           gate_interval

      igate_max = min(int(rmax/gate_spacing_m) , n_gates)

      write(6,809) i_scan_mode,n_gates,igate_max,elevation_deg
 809  format
     :(' REMAP_PROCESS > i_scan_mode,n_gates,igate_max,elevation = '
     :                                                 ,i3,2i5,f5.1)

!     Calculate effective range/height at grid point centers (done iteratively)
      
!     First read domain grid info if needed
!     if(.not. l_domain_read)then
      if(.false.)then
          write(6,*)' REMAP_PROCESS > call get_laps_domain_95'
          call get_laps_domain_95(NX_L,NY_L,lat,lon,topo
     1                           ,dum_2d,grid_spacing_cen_m
     1                           ,istatus)
          if(istatus .ne. 1)then
              write(6,*)' ERROR return from get_laps_domain_95'
              return
          endif
          l_domain_read = .true.
      endif

      I4_elapsed = ishow_timer()

      write(6,*)' REMAP_PROCESS > calculate k_eff array'

      height_guess = height_max

      DO j = 1,NY_L
      DO i = 1,NX_L
!         Use guessed height to obtain initial range value
          call latlon_to_radar(lat(i,j),lon(i,j),height_guess           ! I
     1                        ,azimuth,range_new,elevation_dum          ! O
     1                        ,rlat_radar,rlon_radar,rheight_radar)     ! I  

          iter = 0

          if(range_new .le. rmax)then

            range_dum = r_missing_data

!           Converge on height where radar beam hits the grid point
            do while (abs(range_dum-range_new).gt.10. .and. iter.lt.10)          
              range_dum = range_new

!             Obtain height where radar beam hits grid point
              call radar_to_latlon(lat_dum,lon_dum,height_grid          ! O
     1                            ,azimuth,range_dum,elevation_deg      ! I
     1                            ,rlat_radar,rlon_radar,rheight_radar) ! I

!             Use corrected height to obtain more accurate range
              call latlon_to_radar(lat(i,j),lon(i,j),height_grid        ! I
     1                            ,azimuth,range_new,elevation_dum      ! O
     1                            ,rlat_radar,rlon_radar,rheight_radar) ! I  

              iter = iter + 1

              height_guess = height_grid

            enddo 

            k_eff(i,j) = nint(height_to_zcoord(height_grid,i_status))

          else
            k_eff(i,j) = 0

          endif

      ENDDO
      ENDDO

      I4_elapsed = ishow_timer()
c
c     Find elevation index in look up table
c
      ielev = nint(((elevation_deg-MIN_ELEV) * LUT_ELEVS)
     1                  /(MAX_ELEV-MIN_ELEV)                 )
      ielev = max(ielev,0)
      ielev = min(ielev,LUT_ELEVS)
      write(6,*)' REMAP_PROCESS > elev index = ',ielev

      I4_elapsed = ishow_timer()

      write(6,*)' REMAP_PROCESS > Looping through rays and gates'

      if(namelist_parms%l_hybrid_first_gate)then
          if(elevation_deg .lt. 1.0)then     ! < 1.0
              hybrid_range = 50000.
          elseif(elevation_deg .lt. 2.0)then ! Between 1.0 and 2.0
              hybrid_range = 28000.
          elseif(elevation_deg .lt. 3.0)then ! Between 2.0 and 3.0
              hybrid_range = 18000.
          else
              hybrid_range = 0.
          endif

          init_ref_gate_hyb = hybrid_range / gate_spacing_m
          init_ref_gate_actual = max(INITIAL_REF_GATE,init_ref_gate_hyb)       

          write(6,*)
     1        ' l_hybrid_first_gate flag is set, first range/gate = '      
     1        ,hybrid_range,init_ref_gate_actual

      else
          init_ref_gate_actual = INITIAL_REF_GATE

      endif

      azimuth_interval = 360. / float(lut_azimuths)
      write(6,*)' azimuth_interval = ',azimuth_interval

      write(6,*)' first azimuth = ',az_array(1)

      mingate_valid_ref = 99999
      maxgate_valid_ref = 0

      DO 200 jray=1, n_rays

        if(az_array(jray) .ne. r_missing_data)then
            nazi = nint(az_array(jray) / azimuth_interval)
            nazi = mod(nazi,lut_azimuths)
        else
            goto200
        endif

        igate_interval=1

        DO 180 igate=init_ref_gate_actual,igate_max,igate_interval       

          if(lgate_lut(igate))then ! we'll process this gate, it may have data

            igate_lut = igate/GATE_INTERVAL

            iran = gate_elev_to_projran_lut(igate_lut,ielev)

            i = azran_to_igrid_lut(nazi,iran)
            j = azran_to_jgrid_lut(nazi,iran)

            IF (i .eq. 0 .OR. j .eq. 0) GO TO 180

            io = i - ioffset
            jo = j - joffset 

            if(k_eff(i,j) .ne. 0)then
                k = k_eff(i,j)
                klut = gate_elev_to_z_lut(igate_lut,ielev)
                if(jray .eq. 1 .and. idebug .ge. 1)then
                    write(6,5)igate,igate_lut,iran,i,j,k,klut
5                   format(' igate,igate_lut,iran,i,j,k,klut= ',7i6)
                endif
            else
                k = gate_elev_to_z_lut(igate_lut,ielev)
            endif

            IF (k .eq. 0) GO TO 180

            IF (k .lt. k_low .and. elevation_deg .ge. 0.)then
                write(6,*)' Error: inconsistent k values - ',k,k_low
                stop
            ENDIF

!           True will map one tilt in each LAPS level (for testing only)
            if(namelist_parms%l_ppi_mode)then 
                k = i_tilt 
            endif

            IF( lgate_vel_lut(igate) ) THEN

!           IF( igate .ge. INITIAL_VEL_GATE) THEN
c
c      Velocity Data
c
              n_pot_vel(io,jo,k) = n_pot_vel(io,jo,k) + 1
c
c      Map velocity if data present and abs value of velocity is
c      more than 2 ms-1.
c
              vel_value = Velocity(igate,jray)

              IF (abs(vel_value) .lt. VEL_MIS_CHECK .and.
     :            abs(vel_value) .gt. namelist_parms%abs_vel_min ) THEN

                IF(ngrids_vel(io,jo,k).eq.0 .or. 
     1             vel_nyquist .eq. r_missing_data) THEN
                  rvel =  vel_value

                ELSEIF(vel_nyquist .ne. r_missing_data) THEN ! .and. ngrids > 0
                  avgvel=grid_rvel(io,jo,k)/float(ngrids_vel(io,jo,k))
                  nycor=nint(0.5*(avgvel-vel_value)/
     :                     vel_nyquist)
                  rvel=vel_value+((2*nycor)*vel_nyquist)

                END IF

                n_obs_vel = n_obs_vel + 1
!               if(n_obs_vel .le. 100)then
!                   write(6,*)'n_obs_vel = ',n_obs_vel,i,j,k
!               endif
                grid_rvel(io,jo,k) = grid_rvel(io,jo,k) + rvel
                grid_rvel_sq(io,jo,k) =
     :          grid_rvel_sq(io,jo,k) + rvel*rvel

                ngrids_vel(io,jo,k) = ngrids_vel(io,jo,k) + 1

                if(vel_nyquist .ne. r_missing_data)then
                    grid_nyq(io,jo,k)=grid_nyq(io,jo,k)+vel_nyquist
                endif

              END IF

            END IF
c
c     Map reflectivity
c
            IF( lgate_ref_lut(igate) ) THEN

              n_pot_ref(io,jo,k) = n_pot_ref(io,jo,k) + 1

              ref_value = Reflect(igate,jray)

              IF (abs(ref_value) .lt. REF_MIS_CHECK) THEN ! datum is present

c               grid_ref(io,jo,k) =
c    :          grid_ref(io,jo,k) + ref_value

                if(jray .eq. 1 .and. idebug .ge. 2)then
                    write(6,170)ref_value
 170                format(45x,'Ref = ',f8.2)
                endif

                ilut_ref = nint(ref_value * 10.) ! tenths of a dbz
                grid_ref(io,jo,k) =
     :          grid_ref(io,jo,k) + dbz_to_z_lut(ilut_ref)
                ngrids_ref(io,jo,k) = ngrids_ref(io,jo,k) + 1

                mingate_valid_ref = min(mingate_valid_ref,igate)
                maxgate_valid_ref = max(maxgate_valid_ref,igate)

              END IF

            ENDIF ! l_gate_ref(igate) = .true. and we process the reflectivity

          ENDIF ! l_gate(igate) = .true. and we need to process this gate

  180   CONTINUE ! igate
  200 CONTINUE ! jray

      write(6,815,err=816)elevation_deg,n_obs_vel
     1                   ,mingate_valid_ref,maxgate_valid_ref
  815 format(' REMAP_PROCESS > End Ray/Gate Loop: Elev= ',F10.2
     :      ,'  n_obs_vel = ',I9,' Min/Max Ref Gates ',2I7)

  816 I4_elapsed = ishow_timer()

      IF (i_last_scan .eq. 1) THEN
        write(6,820)
  820   format(
     :  ' REMAP_PROCESS > Last Sweep - Dividing velocity & ref arrays')       


        n_vel_grids_prelim = 0
        n_vel_grids_final = 0
        n_ref_grids = 0
        n_ref_grids_qc_fail = 0

c
c     Diagnostic print-out
c
        write(6,825)
  825   format(' REMAP_PROCESS > Prepare reflectivity Output')

        DO 480 k = 1, k_low-1
        DO 480 jo = 1, NY_R
        DO 480 io = 1, NX_R
          grid_ref(io,jo,k)=r_missing_data
          grid_rvel(io,jo,k)=r_missing_data
          grid_nyq(io,jo,k)=r_missing_data
  480   CONTINUE

        DO 500 k = k_low,NZ_L

          write(6,826) k
  826     format(' REMAP_PROCESS > Dividing: k = ',i2)

          DO 400 jo = 1,NY_R
          DO 400 io = 1,NX_R
c
c     NOTE MIN_VEL_SAMPLES MUST BE GREATER THAN 1
c
            IF(ngrids_vel(io,jo,k) .ge. MIN_VEL_SAMPLES) THEN ! Good gates
              vknt=float(ngrids_vel(io,jo,k))

              IF (vknt .ge. float(n_pot_vel(io,jo,k))
     1                                          * COVERAGE_MIN_VEL) THEN       

                n_vel_grids_prelim = n_vel_grids_prelim + 1
                variance=( grid_rvel_sq(io,jo,k) - 
     :                    (grid_rvel(io,jo,k)*grid_rvel(io,jo,k)/vknt) )
     :                     /(vknt-1.)

                IF (variance .lt. RV_VAR_LIM) THEN ! increment good counter

                  n_vel_grids_final = n_vel_grids_final + 1
                  grid_rvel(io,jo,k) = grid_rvel(io,jo,k)/vknt
                  if(vel_nyquist .ne. r_missing_data)then
                      grid_nyq(io,jo,k) = grid_nyq(io,jo,k)/vknt
                  else
                      grid_nyq(io,jo,k) = r_missing_data
                  endif

                ELSE ! Failed VEL QC test

                  grid_rvel(io,jo,k) = r_missing_data
                  grid_nyq(io,jo,k) = r_missing_data
    
                END IF ! VEL QC test

              ELSE ! Insufficient coverage

                grid_rvel(io,jo,k) = r_missing_data
                grid_nyq(io,jo,k) = r_missing_data

              END IF ! Velocity Coverage check

            ELSE ! Insufficient velocity count

              grid_rvel(io,jo,k) = r_missing_data
              grid_nyq(io,jo,k) = r_missing_data

            END IF ! First check of velocity count
c
c     Reflectivity data
c
!        QC flags within the gridded reflectivities are as follows...
!
!            r_missing_data     Insufficient number of "potential gates" in
!                               the grid volume
!
!            -101.              Insufficient fractional coverage of "actual"
!                               gates in grid volume
!
!            -102.              Reflectivity less than threshold value


            IF(ngrids_ref(io,jo,k) .ge. MIN_REF_SAMPLES) THEN ! Good gates
              rknt=float(ngrids_ref(io,jo,k))
              IF (rknt .ge. float(n_pot_ref(io,jo,k)) 
     1                                          * COVERAGE_MIN_REF) THEN       

!               Calculate mean value of Z
                grid_ref(io,jo,k) = grid_ref(io,jo,k)/rknt

!               Convert from Z to dbZ
                grid_ref(io,jo,k) = alog10(grid_ref(io,jo,k)) * 10.

                IF (grid_ref(io,jo,k) .ge. REF_MIN) THEN

                  n_ref_grids = n_ref_grids + 1
                  IF(n_ref_grids .lt. 200 .and. idebug .ge. 1)
     :               write(6,835) io,jo,k,grid_ref(io,jo,k)
  835                format(' Grid loc: ',3(i4,','),'  Refl: ',f6.1)

                ELSE        ! Failed REF QC test
 
                  n_ref_grids_qc_fail = n_ref_grids_qc_fail + 1
                  grid_ref(io,jo,k) = -102.

                END IF      ! Passed REF QC test

              ELSE       ! Insufficent coverage

                grid_ref(io,jo,k) = -101.

              END IF     ! coverage check of count

            ELSE       ! Insufficent data count

              grid_ref(io,jo,k) = r_missing_data

            END IF     ! first check of count

  400     CONTINUE ! io,jo
  500   CONTINUE ! k

        I4_elapsed = ishow_timer()
c
c     Call QC routine (Now Disabled)
c
        istatus_qc = 1
c       call radar_qc(NX_L,NY_L,NZ_L,grid_rvel,istatus_qc)
        IF (istatus_qc .ne. 1) THEN
          i_num_finished_products = 0
          write(6,840)
  840     format(' REMAP_PROCESS > Bad data detected, no data written')       
          GO TO 998 ! abnormal return
        END IF

        write(6,842) n_ref_grids_qc_fail,n_ref_grids
  842   format(' REMAP_PROCESS > N_REF_QC_FAIL/N_REF = ',2I12)

        IF (n_ref_grids .lt. REF_GRIDS_CHECK) THEN
          i_num_finished_products = 0
          write(6,845) n_ref_grids,REF_GRIDS_CHECK
  845     format(' REMAP_PROCESS > ',i4,' ref grids < ',i4
     :                                 ,'no data file written...')
          GO TO 999 ! normal return
        END IF

        write(6,851)n_ref_obs_old(1),n_ref_grids,i4time_old(1)
     1                                          ,i_product_i4time

  851   format(' REMAP_PROCESS > Ref Obs: Old/New',2i6
     :        ,' I4time: Old/New',2i11)

        i4time_old(1) = i_product_i4time
        n_ref_obs_old(1) = n_ref_grids
c
c     Determine filename extension
        ext = laps_radar_ext
        write(6,*)' REMAP_PROCESS > laps_ext = ',laps_radar_ext
c
c     Prepare to write out data
c
        I4_elapsed = ishow_timer()
c
        var_a(1) = 'REF'
        var_a(2) = 'VEL'
        var_a(3) = 'NYQ'
        units_a(1) = 'dBZ'
        units_a(2) = 'M/S'
        units_a(3) = 'M/S'
        comment_a(1) = 'Doppler Reflectivity'
        comment_a(2) = 'Doppler Velocity'
        comment_a(3) = 'Nyquist Velocity'
        nf = 3

c       DO 555 k=7,9
c       print *, 'sample data on level ',k
c       DO 555 jo=1,NY_R
c       DO 555 io=60,60
c         print *,io,jo,grid_ref(io,jo,k),grid_rvel(io,jo,k)
c 555   CONTINUE
c

        v_nyquist_vol = -999.
        write(6,875) i_tilt
  875   format(' Determine v_nyquist for the ',i4,' tilt volume')
c
        DO 600 i = 1,i_tilt
          write(6,880) i,v_nyquist_tilt(i)
  880     format(' i_tilt:',I6,'  v_nyquist_tilt:',e12.4)
          IF (v_nyquist_tilt(i) .gt. 0.) THEN
            IF (v_nyquist_vol .gt. 0.) THEN
              IF (v_nyquist_tilt(i) .ne. v_nyquist_vol) THEN
                v_nyquist_vol = r_missing_data
                write(6,886)
  886           format(' Nyquist has changed for the tilt',
     1                 ', set v_nyquist_vol to missing.')
                GO TO 601
              END IF
            ELSE
              v_nyquist_vol = v_nyquist_tilt(i)
            END IF
          END IF
  600   CONTINUE
  601   CONTINUE
c
c     Write out header type info into the comment array
c
        l_unfold=.false.
        write(comment_a(1),888)rlat_radar,rlon_radar,rheight_radar
     1       ,n_ref_grids,c4_radarname
        write(comment_a(2),889)rlat_radar,rlon_radar,rheight_radar
     1       ,n_vel_grids_final,c4_radarname,v_nyquist_vol,l_unfold
        write(comment_a(3),888)rlat_radar,rlon_radar,rheight_radar
     1       ,n_vel_grids_final,c4_radarname

  888   format(2f9.3,f8.0,i7,a4,3x)
  889   format(2f9.3,f8.0,i7,a4,3x,e12.4,l2)

        write(6,890)comment_a(1)(1:80)
        write(6,890)comment_a(2)(1:80)
        write(6,890)comment_a(3)(1:80)
  890   format(a80)

        I4_elapsed = ishow_timer()

        if(laps_radar_ext(1:3) .ne. 'vrc')then ! vxx output

            I4_elapsed = ishow_timer()

            call ref_fill_horz(grid_ref,NX_L,NY_L,NZ_L
     1                        ,lat,lon,dgr
     1                        ,NX_R,NY_R,ioffset,joffset
     1                        ,rlat_radar,rlon_radar,rheight_radar
     1                        ,istatus)       
            if(istatus .ne. 1)then
                write(6,*)' Error calling ref_fill_horz'          
                return
            endif

            allocate( out_array_4d(NX_L,NY_L,NZ_L,3), STAT=istat_alloc )       
            if(istat_alloc .ne. 0)then
                write(6,*)' ERROR: Could not allocate out_array_4d'
                stop
            endif

            I4_elapsed = ishow_timer()

            write(6,*)' Filling output arrays'

            if(l_offset_radar)then

                out_array_4d = r_missing_data ! Initialize

                iomin = max((1-ioffset),1)
                iomax = min((NX_L-ioffset),NX_R)

                jomin = max((1-joffset),1)
                jomax = min((NY_L-joffset),NY_R)

                write(6,*)' io range: ',iomin,iomax,
     1                    ' jo range: ',jomin,jomax

                do k = 1,NZ_L
                    do jo = jomin,jomax ! 1,NY_R
                        j = jo + joffset
                        do io = iomin,iomax ! 1,NX_R
                            i = io + ioffset
                            out_array_4d(i,j,k,1) = grid_ref(io,jo,k)
                            out_array_4d(i,j,k,2) = grid_rvel(io,jo,k)       
                            out_array_4d(i,j,k,3) = grid_nyq(io,jo,k)
                        enddo ! io
                    enddo ! jo
                enddo ! k

            else
                out_array_4d(:,:,:,1) = grid_ref(:,:,:)
                out_array_4d(:,:,:,2) = grid_rvel(:,:,:)
                out_array_4d(:,:,:,3) = grid_nyq(:,:,:)

            endif

            I4_elapsed = ishow_timer()

            call make_fnam_lp(i_product_i4time,a9time,istatus)
            if(istatus .ne. 1)return

            call s_len(ext,len_ext)

            if(l_compress_output)then
                write(6,865) c4_radarname,ext(1:len_ext),a9time
 865            format(
     1             ' REMAP_PROCESS > Calling put_compressed_multi_3d'          
     1             ,1x,a4,2x,a,2x,a9)          

                call put_compressed_multi_3d(i_product_i4time,ext,var_a       
     1                                ,units_a,comment_a,out_array_4d
     1                                ,NX_L,NY_L,NZ_L,nf,istatus)

            else
                write(6,866) c4_radarname,ext(1:len_ext),a9time
 866            format(' REMAP_PROCESS > Calling put_laps_multi_3d'
     1                 ,1x,a4,2x,a,2x,a9)          

                call put_laps_multi_3d(i_product_i4time,ext,var_a
     1                                ,units_a,comment_a,out_array_4d
     1                                ,NX_L,NY_L,NZ_L,nf,istatus)

            endif

            deallocate(out_array_4d)

        else ! Single level of data (as per WFO)
            call put_remap_vrc(i_product_i4time,comment_a(1)
     1             ,rlat_radar,rlon_radar,rheight_radar
     1             ,dgr
     1             ,grid_ref,NX_L,NY_L,NZ_L
     1             ,c3_radar_subdir,path_to_vrc,r_missing_data,istatus)       

        endif

        I4_elapsed = ishow_timer()

!       go to 900

900     continue

      END IF ! i_last_scan

      go to 999 ! normal return

!     Return section

998   i_status = 0
      write(6,*) ' WARNING: Return from remap_process with 0 status'
      RETURN

999   i_status = 1
      RETURN

      END


        subroutine purge(ext,nfiles,ntime_min,i4time_now)

!       Keeps number of files according to nfiles or time span according to
!       ntime_min, whichever is greater

        integer MAX_FILES
        parameter (MAX_FILES = 1000)

        character*9 asc_tim_9
        character*31 ext
        character*255 c_filespec
        character c_fnames(MAX_FILES)*80

        call s_len(ext,len_ext)

        c_filespec = '../lapsprd/'//ext(1:len_ext)//'/*.'
     1                            //ext(1:len_ext)//'*'         

        write(6,*)c_filespec

        call    Get_file_names(  c_filespec,
     1			 i_nbr_files_ret,
     1			 c_fnames,MAX_FILES,
     1			 i_status )

        if(i_nbr_files_ret .gt. 0)then
            call get_directory_length(c_fnames(1),lenf)
            write(6,*)i_nbr_files_ret,' file(s) in directory'
        else ! Error Condition
            write(6,*)' No files in directory'
            istatus = 0
            return
        endif

        ntime_sec = ntime_min * 60

10      do i=1,i_nbr_files_ret-nfiles ! Loop through excess versions
            asc_tim_9 = c_fnames(i)(lenf+1:lenf+9)
            call i4time_fname_lp(asc_tim_9,I4time_file,istatus)
            if(i4time_now - i4time_file .gt. ntime_sec)then ! File is too old

!               Delete the file
!               call rm_file(c_fnames(i)(1:lenf+13),istatus)
                call rm_file(c_fnames(i),istatus)


            endif
        enddo

        return
        end


        subroutine rm_file(c_filename,istatus)

        character*(*) c_filename

        integer istatus

        lun = 151

        write(6,*)' rm_file ',c_filename

        open(lun,file=c_filename,status='unknown')

        close(lun,status='delete')

        istatus = 1
 
        return
        end



        subroutine put_remap_vrc(i4time,comment_2d 
     1                         ,rlat_radar,rlon_radar,rheight_radar
     1                         ,dgr
     1                         ,field_3d,imax,jmax,kmax,c3_radar_subdir        
     1                         ,path_to_vrc,r_missing_data,istatus)

!       Stuff from 'put_laps_2d' except how do we handle radar subdir?

        character*150 DIRECTORY
        character*150 DIRECTORY1
        character*31 EXT

        character*125 comment_2d
        character*125 comments_2d(2)
        character*10 units(2)
        character*3 vars(2)

        integer LVL_2D(2)

        real field_3d(imax,jmax,kmax)
        real fields_2d(imax,jmax,2)
        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)
        real dist(imax,jmax)

        character*9 a9time
        character*8 radar_subdir
        character*3 c3_radar_subdir
        character*(*) path_to_vrc
        character*4 LVL_COORD_2D(2)

        call make_fnam_lp(i4time,a9time,istatus)
        if(istatus .ne. 1)return

        write(6,*)' Subroutine put_remap_vrc for ',a9time

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)return

!       Get column max reflectivity (now passing in r_missing_data)
        call get_max_reflect(field_3d,imax,jmax,kmax,r_missing_data
     1                      ,fields_2d(1,1,1) )

        call get_laps_domain(imax,jmax,'nest7grid',lat,lon,topo,istatus)       
        if(istatus .ne. 1)then
            write(6,*)' Error calling get_laps_domain'
            return
        endif

!       Calculate closest radar array
        write(6,*)' Calculating closest radar array (dist to vrc radar)'       
        do i = 1,imax
        do j = 1,jmax
            call latlon_to_radar(lat(i,j),lon(i,j),topo(i,j)
     1                          ,azimuth,dist(i,j),elev
     1                          ,rlat_radar,rlon_radar,rheight_radar)       
        enddo ! j
        enddo ! i

!       call vrc_clutter_thresh(      fields_2d(1,1,1)                   ! I/O
!    1                               ,dist                               ! I
!    1                               ,imax,jmax,ref_base,r_missing_data) ! I

        call ref_fill_horz(fields_2d(1,1,1),imax,jmax,1,lat,lon,dgr
     1                    ,imax,jmax,0,0
     1                    ,rlat_radar,rlon_radar,rheight_radar,istatus)       
        if(istatus .ne. 1)then
            write(6,*)' Error calling ref_fill_horz'          
            return
        endif

!       Utilize closest radar array
        write(6,*)' Utilizing closest radar array (dist to vrc radar)'       
        do i = 1,imax
        do j = 1,jmax
            if(fields_2d(i,j,1) .ne. r_missing_data)then
                fields_2d(i,j,2) = dist(i,j)
            else
                fields_2d(i,j,2) = r_missing_data
            endif
        enddo ! j
        enddo ! i

        ext = 'vrc'

        vars(1) = 'REF'
        vars(2) = 'DIS'

        units(1) = 'DBZ'
        units(2) = 'M'

        lvl_2d(1) = 0
        lvl_2d(2) = 0
        
        lvl_coord_2d(1) = 'MSL'
        lvl_coord_2d(2) = 'MSL'

        comments_2d(1) = comment_2d
        comments_2d(2) = comment_2d

        write(6,*)'path_to_vrc = ',path_to_vrc

        if(path_to_vrc .eq. 'rdr')then
            radar_subdir = c3_radar_subdir
            write(6,*)' radar_subdir = ',radar_subdir

            call get_directory('rdr',directory1,len_dir1)

            directory = directory1(1:len_dir1)//radar_subdir(1:3)
     1                                        //'/vrc/'  
            call s_len(directory,len_dir)

        else ! 'lapsprd'
            call get_directory('vrc',directory,len_dir)

        endif            

        write(6,11)directory(1:len_dir),ext(1:5),vars
11      format(' Writing 2d ',a,1x,a5,2(1x,a3))

        CALL WRITE_LAPS_DATA(I4TIME,DIRECTORY,EXT,imax,jmax,
     1                       2,2,vars,LVL_2D,LVL_COORD_2D,units,
     1                       comments_2d,fields_2d,ISTATUS)
        if(istatus .eq. 1)then
            write(6,*)' VRC successfully written'
        else
            write(6,*)' VRC not successfully written', istatus
        endif


        return
        end

        subroutine vrc_clutter_thresh(ref                             ! I/O
     1                               ,dist                            ! I
     1                               ,ni,nj,ref_base,r_missing_data)  ! I

!       Apply a range dependent reflectivity threshold to filter ground clutter

        real ref(ni,nj), dist(ni,nj)

        do i = 1,ni
        do j = 1,nj
            if(dist(i,j) .lt. 80000.)then
                thresh = 10.
            elseif(dist(i,j) .gt. 100000.)then
                thresh = 0.
            else
                thresh = 10. * (100000. - dist(i,j)) / 20000.
            endif             

            if(      ref(i,j) .lt. thresh 
     1         .and. ref(i,j) .gt. ref_base
     1         .and. ref(i,j) .ne. r_missing_data )then       
                ref(i,j) = ref_base
            endif
        
        enddo ! j
        enddo ! i

        return
        end
