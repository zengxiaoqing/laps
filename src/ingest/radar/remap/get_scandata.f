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
        subroutine  get_scandata(
     :          i_tilt,             ! Input  (Integer*4)
     :          r_missing_data,     ! Input  (Real*4)
     :          max_rays,           ! Input  (Integer*4)
     :          n_rays,             ! Output (Integer*4)
     :          n_gates,            ! Output (Integer*4)
     :          gate_spacing_m,     ! Output (Real*4)
     :          elevation_deg,      ! Output (Real*4)
     :          i_scan_mode,        ! Output (Integer*4)
     :          v_nyquist_tilt,     ! Output (Real*4)
     :          v_nyquist_ray,      ! Output (Real*4 array)
     :          istatus )           ! Output (Integer*4)

c
c     PURPOSE:
c
c       Get info for radar data.
c
      implicit none
c
c     Input variables
c
      integer*4 i_tilt
      real*4 r_missing_data
      integer*4 max_rays
c
c     Output variables
c
      integer*4 n_rays
      integer*4 n_gates
      real*4 gate_spacing_m
      real*4 elevation_deg
      integer*4 i_scan_mode
      real*4 v_nyquist_tilt
      real*4 v_nyquist_ray(max_rays)
      integer*4 istatus
c
c     Include file
c
      include 'remap_buffer.cmn'
c
c     Misc internal variables
c
      integer i
c
      n_rays = n_rays_cmn
      n_gates = 1840
      gate_spacing_m = 250.
      elevation_deg = elev_cmn
      i_scan_mode = 1
c
c     Transfer nyquist info from coomon to output array
c
      DO 100 i = 1,n_rays
        IF (v_nyquist_ray_a_cmn(i) .eq. 0.) THEN
          v_nyquist_ray(i) = r_missing_data
        ELSE
          v_nyquist_ray(i) = v_nyquist_ray_a_cmn(i)
        END IF
  100 CONTINUE
c
c     Check all nyquists in this tilt.
c     If they are all the same, report that as the tilt nyquist,
c     else report r_missing
c
      v_nyquist_tilt = v_nyquist_ray(1)
      DO 200 i = 1,n_rays
        IF (v_nyquist_ray(i) .eq. r_missing_data) THEN
          v_nyquist_tilt = r_missing_data
          GO TO 999
        ELSE IF (v_nyquist_ray(i) .ne. v_nyquist_tilt) THEN
          v_nyquist_tilt = r_missing_data
          GO TO 999
        END IF
  200 CONTINUE

  999 CONTINUE

      istatus = 1
      RETURN
      END
