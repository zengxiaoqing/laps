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
      subroutine read_data_88D(
     :               i_tilt,           ! Input
     :               vel_thr_rtau,
     :               r_missing_data,
     :               gate_spacing_m_ret,   ! Output
     :               num_sweeps,
     :               elevation_deg,
     :               n_rays,
     :               n_gates,
     :               slant_ranges_m,
     :               velocity,
     :               reflect,
     :               az_Array,
     :               v_nyquist_tilt,
     :               istatus)
c
c     Subroutine read_data_xm
c
c       Steve Albers    15-MAR-1988
c       SA                 DEC-1988   Modified more for RT purposes
c       SA                     1994   WSR-88D Remapper
c       Keith Brewster     JUN-1994   Clean-up
c
      implicit None
c
c     Include files
c
      include         'remap_constants.dat'
      include         'remap_dims.inc'
      include         'remap_buffer.cmn'
c
c     Input variables
c
      integer*4 i_tilt
      real*4 vel_thr_rtau
      real*4 r_missing_data
c
c     Output variables
c
      real*4 gate_spacing_m_ret
      integer*4 num_sweeps
      real*4 elevation_deg
      integer*4 n_rays
      integer*4 n_gates
      real*4 slant_ranges_m(max_gates)
      real*4 velocity(max_gates,max_rays)
      real*4 reflect(max_gates,max_rays)
      real*4 az_array(max_rays)
      real*4 v_nyquist_tilt
      integer*4 istatus
c
c     Misc local variables
c
      integer i,i_scan_mode,i_ray
      real v_nyquist_ray(max_rays)
c
c     Get housekeeping info for radar scan.
c
      call get_scandata(
     :          i_tilt,
     :          r_missing_data,
     :          max_rays,
     :          n_rays,
     :          n_gates,
     :          gate_spacing_m_ret,
     :          elevation_deg,
     :          i_scan_mode,
     :          v_nyquist_tilt,
     :          v_nyquist_ray,
     :          istatus)
      IF ( istatus .ne. 1 ) GO TO 998
c
      IF (gate_spacing_m_ret .ne. gate_spacing_m) THEN
        write(6,805) gate_spacing_m_ret,gate_spacing_m
  805   format(' Error, returned gate spacing is different from',/,
     :         ' parameter value',2F12.2)
        istatus = 0
        GO TO 998
      END IF

      call GET_AZIMUTHS_DEG ( i_tilt,
     :                  max_rays,
     :                  az_array,
     :                  istatus )
      IF( istatus .ne. 1 )  GO TO 998

      DO 50 i = 1,N_GATES
        slant_ranges_m(i) = float(i) * gate_spacing_m_ret
   50 CONTINUE

      write(6,815)
  815 format(' READ_DATA_88D > Loading arrays')

      DO 100 i_ray = 1,N_RAYS
        call ld_ray(
     :              i_ray,
     :              n_gates,
     :              r_missing_data,
     :              i_missing_data,
     :              velocity(1,i_ray),
     :              reflect(1,i_ray),
     :              istatus)
        IF( istatus .ne. 1 )  GO TO 998
  100 CONTINUE

      RETURN

  998 write(6,825)
  825 format(' Bad Status in Readdata')
      RETURN
      END
