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
        subroutine  get_scandata(
     :          i_tilt,             ! Input  (Integer)
     :          r_missing_data,     ! Input  (Real)
     :          max_rays_in,        ! Input  (Integer)
     :          n_rays,             ! Output (Integer)
     :          n_gates,            ! Output (Integer)
     :          gate_spacing_m_ret, ! Output (Real)
     :          elevation_deg,      ! Output (Real)
     :          i_scan_mode,        ! Output (Integer)
     :          v_nyquist_tilt,     ! Output (Real)
     :          v_nyquist_ray,      ! Output (Real array)
     :          istatus )           ! Output (Integer)

c
c     PURPOSE:
c
c       Get info for radar data.
c
      implicit none
c
c     Input variables
c
      integer i_tilt
      real r_missing_data
      integer max_rays_in
c
c     Output variables
c
      integer n_rays
      integer n_gates
      real gate_spacing_m_ret
      real elevation_deg
      integer i_scan_mode
      real v_nyquist_tilt
      real v_nyquist_ray(max_rays_in)
      integer istatus
c
c     Include file
c
      include 'remap_dims.inc'
      include 'remap_constants.dat'
      include 'remap_buffer.cmn'
c
c     Misc internal variables
c
      integer i
c
      n_rays = n_rays_cmn
      n_gates = MAX_GATES
      gate_spacing_m_ret = GATE_SPACING_M
      elevation_deg = elev_cmn
      i_scan_mode = 1
c
c     Transfer nyquist info from common to output array
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
