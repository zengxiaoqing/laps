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
      subroutine fill_common(ref,vel,n_rays,i_tilt,n_ref_gates,
     :     n_vel_gates,azim,v_nyquist,eleva,rmissing_data)
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
      implicit none
c
c     Sizing parameters
c
      integer max_ray,max_ref_gates,max_vel_gates,n_ref,n_vel
      parameter (max_ray=380,max_ref_gates=460,max_vel_gates=920)
      parameter (n_ref=(max_ray*max_ref_gates),
     :           n_vel=(max_ray*max_vel_gates))
c
c     Input declarations
c
      real*4 ref(n_ref)
      real*4 vel(n_vel)
      integer*4 n_rays
      integer*4 i_tilt
      integer*4 n_ref_gates
      integer*4 n_vel_gates
      real*4 azim(max_ray)
      real*4 v_nyquist(max_ray)
      real*4 eleva
      real*4 rmissing_data
c
c     Common block declarations
c
      include 'remap_buffer.cmn'
c
c     Misc internal variables
c
      integer i,j,iref_start,ivel_start
c
      print *, ' REMAP (FILL_COMMON)'
      print *, ' ref(1..4) ',
     :             ref(1),ref(2),ref(3),ref(4)
      print *, ' vel(1..4) ',
     :             vel(1),vel(2),vel(3),vel(4)
      print *, ' n_rays = ',n_rays
      print *, ' i_tilt = ',i_tilt
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
