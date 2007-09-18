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
      subroutine fill_common(ref,vel,n_rays,i_tilt,n_ref_gates,
     :     n_vel_gates,
     :     first_ref_m,gsp_ref_m,first_vel_m,gsp_vel_m,
     :     azim,v_nyquist,eleva,rmissing_data)
c
c     Subroutine fill_common
c
c     PURPOSE:
c     Transfer a tilt of radar data passed-in to data array in common block.
c
c     AUTHOR:
c     Steve Albers, FSL
c
c     MODIFICATIONS:
c     Clean-up    Keith Brewster, CAPS, July, 1994
c
c
c     Sizing parameters
c
      include 'remap_dims.inc'
      include 'remap_buffer.cmn'
!     include 'remap_constants.dat' ! for debugging only
!     include 'remap.cmn' ! for debugging only

      integer n_ref,n_vel
      parameter (n_ref=(max_ray_tilt*max_ref_gates),
     :           n_vel=(max_ray_tilt*max_vel_gates))
c
c     Input declarations
c
      real ref(n_ref)
      real vel(n_vel)
      integer n_rays
      integer i_tilt
      integer n_ref_gates
      integer n_vel_gates
      real azim(max_ray_tilt)
      real v_nyquist(max_ray_tilt)
      real eleva
      real rmissing_data
c
c     Misc internal variables
c
      integer i,j,iref_start,ivel_start
c
      write(6,*)' fill_common: i_tilt = ',i_tilt
      print *, ' ref(1..4) ',
     :             ref(1),ref(2),ref(3),ref(4)
      print *, ' vel(1..4) ',
     :             vel(1),vel(2),vel(3),vel(4)
      print *, ' n_rays = ',n_rays
      print *, ' n_ref_gates = ',n_ref_gates
      print *, ' n_vel_gates = ',n_vel_gates
      print *, ' azim(1..4) = ',azim(1),azim(2),azim(3),azim(4)
      print *, ' Nyquis(1..4) = ',v_nyquist(1),
     :   v_nyquist(2),v_nyquist(3),v_nyquist(4)
      print *, ' eleva = ',eleva
      print *, ' rmissing_data = ',rmissing_data
c
c     Fill in housekeeping stuff in common
c
      n_rays_cmn = n_rays
      i_tilt_cmn = i_tilt
      n_ref_gates_cmn = n_ref_gates
      n_vel_gates_cmn = n_vel_gates

      first_ref_m_cmn = first_ref_m
      first_vel_m_cmn = first_vel_m
      gsp_ref_m_cmn = gsp_ref_m
      gsp_vel_m_cmn = gsp_vel_m
c
      DO 100 i = 1,n_rays
        azim_cmn(i) = azim(i)
        v_nyquist_ray_a_cmn(i) = v_nyquist(i)
  100 CONTINUE

      elev_cmn = eleva
      b_missing_cmn = nint(2.*rmissing_data)
c
c     Fill the common arrays with the scan data
c
      DO 200 j = 1,n_rays
c
c     Reflectivity
c
        iref_start = (j-1) * n_ref_gates
        DO 140 i = 1,n_ref_gates
          b_ref_cmn(i,j) = nint(2.*ref(iref_start + i))
  140   CONTINUE
c
c     Velocity
c
        ivel_start = (j-1) * n_vel_gates
        DO 150 i = 1,n_vel_gates
          b_vel_cmn(i,j) = nint(2.*vel(ivel_start + i))
  150   CONTINUE
  200 CONTINUE

      RETURN
      END
