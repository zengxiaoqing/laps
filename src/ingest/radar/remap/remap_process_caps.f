cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
      Subroutine Remap_process(
     :                  i_tilt,                   ! Integer*4 (input)
     :                  i_last_scan,              ! Integer*4 (input)
     :                  i_first_scan,             ! Integer*4 (input)
     :                  i_product_i4time,         ! Integer*4 (input)
     :                  full_fname,               ! Character*91
     :                  fn_length,
     :                  i_num_finished_products,  ! Integer*4 (output)
     :                  i_status )                ! Integer*4 (output)
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
c       Brewster, Keith    AUG_1994     Clean-out of artifacts
c
c ************************** Declaration Section **************************
c
      implicit none
c
c     LAPS Grid Dimensions
c
      include 'lapsparms.for'
      include 'remap_constants.dat'
      include 'remap.cmn'
c
c     Input variables
c
      integer*4 i_tilt
      integer*4 i_last_scan
      integer*4 i_first_scan
      integer*4 i_product_i4time
c
c     Output variables
c
      character*91 full_fname
      integer*4 fn_length
      integer*4 i_num_finished_products
      integer*4 i_status
c
c     Processing parameters
c
      real re43
      parameter (re43 = 8503700.) ! 4/3 radius of the earth in meters
c
c
      integer*4 max_fields
      parameter (max_fields = 10)
c
c     Variables for NetCDF I/O
c
      character*50 dir
      
      character*9 gtime
      character*31 ext,ext_in
      character*3 var_a(max_fields)
      character*125 comment_a(max_fields)
      character*10  units_a(max_fields)

      real*4 out_array_4d(NX_L,NY_L,NZ_L,2)
c
c     Functions
c
      real height_to_zcoord
      real height_of_level
c
c     Misc local variables
c
      integer igate,i_scan_mode,jray,end_ext
c
      Real*4  Slant_ranges_m (max_gates),
     :        Elevation_deg,
     :        Az_array(max_rays),
     :        Velocity(max_gates,max_rays),
     :        Reflect(max_gates,max_rays)
c 
      real*4 vel_nyquist
      real*4 v_nyquist_tilt(max_tilts)
      real*4 v_nyquist_vol
c
      integer i,j,k,k_low,ielev,igate_lut
      integer nazi,iran
      integer num_sweeps,n_rays,n_gates,n_obs_vel,n_output_data,nf
      integer igate_start,igate_max
      integer igate_interval
      integer n_vel_grids_final,n_vel_grids_prelim
      integer n_ref_grids,n_ref_grids_qc
      integer istatus,istatus_qc
      real vel_thr_rtau,rvel
      real rmax,height_max,rlat_radar,rlon_radar,rheight_radar
      real vknt,rknt,std_dev,r_meansq
c
c     Beginning of executable code
c
      rlat_radar = plat_radar
      rlon_radar = plon_radar
      rheight_radar = pheight_radar

      i_num_finished_products = 0

      write(6,805) i_first_scan,i_last_scan,i_tilt
  805 format(' REMAP > PROCESS V940802  ifirst,ilast,tilt'
     :                                           ,4i5)
c
c     For first scan, initialize sums and counters to zero.
c
      IF (i_first_scan .eq. 1 .or. i_first_scan .eq. 999) THEN

        write(6,806)
  806 format(' REMAP > 1st sweep - Initializing vel/ref arrays')

        n_obs_vel = 0
        n_output_data = 0

        DO 100 k = 1,NZ_L
        DO 100 j = 1,NY_L
        DO 100 i = 1,NX_L
          grid_rvel(i,j,k) = 0.
          grid_ref(i,j,k) = 0.
          ngrids_vel(i,j,k) = 0
          ngrids_ref(i,j,k) = 0
          n_pot_vel(i,j,k) = 0
          n_pot_ref(i,j,k) = 0
  100   CONTINUE
c
c     Compute maximum height of data needed.
c
        height_max = height_of_level(NZ_L)
c
c     Define Lower Limit of Radar Coverage in LAPS grid
c
        k_low = int(height_to_zcoord(rheight_radar,i_status))

      END IF
c
c     Get radar data from the storage area.
c
      call Read_Data_88D(
     :               i_tilt,
     :               vel_thr_rtau,
     :               r_missing_data,   ! Input
     :               Gate_spacing_m,   ! Output
     :               Num_sweeps,
     :               Elevation_deg,
     :               n_rays,
     :               n_gates,          ! Ref and Vel are on the same # of gates
     :               Slant_ranges_m,
     :               Velocity,
     :               Reflect,
     :               Az_Array,
     :               vel_nyquist,
     :               i_status)
c
      IF (i_status .ne. 1) GO TO 998
c
      v_nyquist_tilt(i_tilt) = vel_nyquist
c
      write(6,*)' REMAP > PROCESS: vel_nyquist for this tilt = '
     :        ,i_tilt,vel_nyquist
c
c     Find elevation index in look up table
c
      ielev = nint((elevation_deg * LUT_ELEVS)/MAX_ELEV)
      ielev = max(ielev,0)
      ielev = min(ielev,MAX_ELEV)
c
c     Compute max range from elevation angle
c
      rmax = -re43 * sind(elevation_deg)
     :  + sqrt(re43*re43*sind(elevation_deg)**2 +
     :          height_max * (2.*re43 + height_max))

      igate_max = min(int(rmax/gate_spacing_m) , n_gates)

      write(6,809) i_scan_mode,igate_max,elevation_deg
 809  format(' REMAP > i_scan_mode,igate_max,elevation = ',i3,i5,f5.1)

      IF (igate_max .ge. 256) THEN
        igate_interval = 1
      ELSE
        igate_interval = 1
      END IF

      DO 200 jray=1, n_rays

        nazi = nint(az_array(jray))
        IF (nazi .eq. 360) nazi = 0

        IF (mod(nazi,5) .eq. 0) THEN
          igate_start = 2
        ELSE
          igate_start = 30
        END IF

        DO 180 igate=igate_start,igate_max,igate_interval

          igate_lut = igate/GATE_INTERVAL

          iran = gate_elev_to_projran_lut(igate_lut,ielev)

          i = azran_to_igrid_lut(nazi,iran)
          j = azran_to_jgrid_lut(nazi,iran)
          k = gate_elev_to_z_lut(igate_lut,ielev)

          IF (i .eq. 0 .OR. j.eq.0 .OR. k.eq.0 ) GO TO 180
c
c      Velocity Data
c
          n_pot_vel(i,j,k) = n_pot_vel(i,j,k) + 1
c
c      Map velocity if data present and abs value of velocity is
c      more than 2 ms-1.
c
          IF (abs(Velocity(igate,jray)) .lt. VEL_MIS_CHECK .and.
     :        abs(Velocity(igate,jray)) .gt. ABS_VEL_MIN ) THEN

            rvel =  Velocity(igate,jray)
            n_obs_vel = n_obs_vel + 1
            grid_rvel(i,j,k) =
     :          grid_rvel(i,j,k) + rvel
            grid_rvel_sq(i,j,k) =
     :          grid_rvel_sq(i,j,k) + rvel*rvel
            ngrids_vel(i,j,k) =
     :          ngrids_vel(i,j,k) + 1
          END IF
c
c     Map reflectivity
c
          n_pot_ref(i,j,k) = n_pot_ref(i,j,k) + 1

          IF (abs(Reflect(igate,jray)) .lt. REF_MIS_CHECK) THEN

            grid_ref(i,j,k) =
     :          grid_ref(i,j,k) + Reflect(igate,jray)
            ngrids_ref(i,j,k) = ngrids_ref(i,j,k) + 1

          END IF

  180   CONTINUE ! igate
  200 CONTINUE ! jray

      write(6,815) elevation_deg,n_obs_vel
  815 format(' REMAP > Accumulated Sweep ',F10.2,'   n_obs_vel = ',I12)

      IF (i_last_scan .eq. 1) THEN
        write(6,820)
  820   format(' REMAP > Last Sweep - Dividing velocity & ref arrays')

        n_vel_grids_prelim = 0
        n_vel_grids_final = 0
        n_ref_grids = 0
        n_ref_grids_qc = 0

c
c     Diagnostic print-out
c
        write(6,825)
  825   format(' REMAP > PROCESS: Reflectivity Output')

        DO 500 k = k_low,NZ_L

          write(6,826) k
  826     format(' REMAP > k = ',i2)

          DO 400 j = 1,NY_L
          DO 400 i = 1,NX_L

            IF(ngrids_vel(i,j,k) .ge. MIN_VEL_SAMPLES) THEN ! Good gates
              vknt=float(ngrids_vel(i,j,k))

              IF (vknt .ge. float(n_pot_vel(i,j,k))*COVERAGE_MIN) THEN

                n_vel_grids_prelim = n_vel_grids_prelim + 1
                grid_rvel(i,j,k) = grid_rvel(i,j,k)/vknt
                r_meansq = grid_rvel_sq(i,j,k)/vknt
                std_dev = r_meansq - grid_rvel(i,j,k)*grid_rvel(i,j,k)

c               IF (std_dev .lt. 0.) THEN
c                 write(6,830) std_dev,r_meansq,grid_rvel(i,j,k) ** 2 ,
c    1                               r_meansq - grid_rvel(i,j,k) ** 2
c 830             format(' less than zero: ',4F10.2)
c               END IF

                IF (std_dev .lt. RA_RMS_SQ) THEN ! increment good counter

                  n_vel_grids_final = n_vel_grids_final + 1

                ELSE ! Failed VEL QC test

                  grid_rvel(i,j,k) = r_missing_data
    
                END IF ! VEL QC test

              ELSE ! Insufficient coverage

                grid_rvel(i,j,k) = r_missing_data

              END IF ! Velocity Coverage check

            ELSE ! Insufficient velocity count

              grid_rvel(i,j,k) = r_missing_data

            END IF ! First check of velocity count
c
c     Reflectivity data
c
            IF(ngrids_ref(i,j,k) .ge. MIN_REF_SAMPLES) THEN ! Good gates
              rknt=float(ngrids_ref(i,j,k))
              IF (rknt .ge. float(n_pot_ref(i,j,k)) * COVERAGE_MIN) THEN

                grid_ref(i,j,k) = grid_ref(i,j,k)/rknt

                IF (grid_ref(i,j,k) .ge. REF_MIN) THEN

                  n_ref_grids = n_ref_grids + 1
                  IF(n_ref_grids .lt. 20)
     :               write(6,835) i,j,k,grid_ref(i,j,k)
  835                format(' Grid loc: ',3(i4,','),'  Refl: ',f6.1)

                ELSE ! Failed REF QC test
 
                  n_ref_grids_qc = n_ref_grids_qc + 1
                  grid_ref(i,j,k) = r_missing_data

                END IF ! Passed REF QC test

              ELSE ! Insufficent coverage

                grid_ref(i,j,k) = r_missing_data

              END IF   ! coverage check of count

            ELSE ! Insufficent data count

              grid_ref(i,j,k) = r_missing_data

            END IF   ! first check of count

  400     CONTINUE ! i,j
  500   CONTINUE ! k
c
c     Call QC routine (Now Disabled)
c
        istatus_qc = 1
c       call radar_qc(NX_L,NY_L,NZ_L,grid_rvel,istatus_qc)
        IF (istatus_qc .ne. 1) THEN
          i_num_finished_products = 0
          write(6,840)
  840     format(' REMAP > Bad data detected, no data written')
          GO TO 900
        END IF

        IF (n_ref_grids .lt. 10) THEN
          i_num_finished_products = 0
          write(6,845) n_ref_grids
  845     format(' REMAP > Insuf Ref Coverage, no data written: ',I12)
          GO TO 900
        END IF

        write(6,850) n_ref_grids_qc,n_ref_grids
  850   format(' REMAP > N_REF_QC/N_REF = ',2I12)

        write(6,851)n_ref_obs_old(1),n_ref_grids,i4time_old(1)
     1                                          ,i_product_i4time

  851   format(' REMAP > Ref Obs: Old/New',2i6,' I4time: Old/New',2i11)

        i4time_old(1) = i_product_i4time
        n_ref_obs_old(1) = n_ref_grids
c
c     Write out data
c
        write(6,865)
  865   format(' REMAP > Calling write_laps_data')
        ext = 'vrd'
        var_a(1) = 'REF'
        var_a(2) = 'VEL'
        units_a(1) = 'dBZ'
        units_a(2) = 'M/S'
        comment_a(1) = 'Doppler Reflectivity'
        comment_a(2) = 'Doppler Velocity'
        nf = 2

        DO 550 k = 1,NZ_L
        DO 550 j = 1,NY_L
        DO 550 i = 1,NX_L
          out_array_4d(i,j,k,1) = grid_ref(i,j,k)
          out_array_4d(i,j,k,2) = grid_rvel(i,j,k)
  550   CONTINUE

        v_nyquist_vol = v_nyquist_tilt(1)
        write(6,875)
  875   format(' Determine v_nyquist for the volume')
c
        DO 600 i = 1,i_tilt
          write(6,880) i,v_nyquist_tilt(i)
  880     format(' i_tilt: ',I12,'  v_nyquist_tilt: ',F10.2)
          IF (v_nyquist_tilt(i) .eq. r_missing_data) THEN
            v_nyquist_vol = r_missing_data
            write(6,885)
  885       format(' Nyquist is missing for the tilt',
     1             ', set v_nyquist_vol to missing.')
            GO TO 601
          ELSE IF (v_nyquist_tilt(i) .ne. v_nyquist_vol) THEN
            v_nyquist_vol = r_missing_data
            write(6,886)
  886       format(' Nyquist has changed for the tilt',
     1             ', set v_nyquist_vol to missing.')
            GO TO 601
          END IF
  600   CONTINUE
  601   CONTINUE
c
c     Write out header type info into the comment array
c
        write(comment_a(1),888)rlat_radar,rlon_radar,rheight_radar
     1                        ,n_ref_grids,    'KTLX   '
        write(comment_a(2),889)rlat_radar,rlon_radar,rheight_radar
     1                  ,n_vel_grids_final,'KTLX   ',v_nyquist_vol

  888   format(2f9.3,f8.0,i7,a7)
  889   format(2f9.3,f8.0,i7,a7,e10.2)

        write(6,890)comment_a(1)(1:80)
        write(6,890)comment_a(2)(1:80)
  890   format(a80)

        call put_laps_multi_3d(i_product_i4time,ext,var_a,units_a,
     1              comment_a,out_array_4d,NX_L,NY_L,NZ_L,nf,istatus)

        call make_fnam_lp(i_product_i4time,gtime,istatus)
        call downcase(ext,ext_in)
        call s_len(ext_in,end_ext)
        dir='../lapsprd/'//ext_in(1:end_ext)//'/'
        call cvt_fname_data(dir,gtime,ext_in,full_fname,istatus)
        call s_len(full_fname,fn_length)


        write(6,895) n_vel_grids_prelim,n_vel_grids_final,
     1               n_ref_grids,v_nyquist_vol,full_fname
  895   format(
     1    ' REMAP > File written: n_grids prelim/final/ref/v_nyq = '
     1          ,3i7,e10.2,/,
     1    ' REMAP > File name: ',a91)

        i_num_finished_products = 1

      END IF ! i_last_scan

900   i_status = 1
      RETURN

998   i_status = 0
      RETURN

      END
