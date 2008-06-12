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
     :               namelist_parms,
     :               gate_spacing_m_ret,   ! Output
     :               i_scan_mode,
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
      integer i_tilt
      real vel_thr_rtau
      real r_missing_data
c
c     Output variables
c
      real gate_spacing_m_ret
      integer num_sweeps
      real elevation_deg
      integer n_rays
      integer n_gates
      real slant_ranges_m(max_gates)
      real velocity(max_gates,MAX_RAY_TILT)
      real reflect(max_gates,MAX_RAY_TILT)
      real az_array(MAX_RAY_TILT)
      real v_nyquist_tilt
      integer istatus
c
c     Misc local variables
c
      integer i,i_scan_mode,i_ray
      real v_nyquist_ray(MAX_RAY_TILT)
c
c     Get housekeeping info for radar scan.
c
      call get_scandata(
     :          i_tilt,
     :          r_missing_data,
     :          MAX_RAY_TILT,
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
     :                  MAX_RAY_TILT,
     :                  az_array,
     :                  istatus )
      IF( istatus .ne. 1 )  GO TO 998

      DO 50 i = 1,N_GATES
        slant_ranges_m(i) = float(i) * gate_spacing_m_ret
   50 CONTINUE

      call lgate_lut_gen(gsp_ref_m_cmn,gsp_vel_m_cmn
     1                  ,n_ref_gates_cmn,n_vel_gates_cmn) 

      write(6,815)
  815 format(' READ_DATA_88D > Loading arrays')

      DO 100 i_ray = 1,N_RAYS
        call ld_ray(
     :              i_ray,
     :              n_gates,
     :              r_missing_data,
!    :              i_missing_data,
     :              velocity(1,i_ray),
     :              reflect(1,i_ray),
     :              istatus)
        IF( istatus .ne. 1 )  GO TO 998
  100 CONTINUE

      if(namelist_parms%l_line_ref_qc)then
          write(6,*)' read_data_88d: calling rayqckz...'
!         iscale = 2
!         miss = nint(b_missing_data)
!         call rayqckz(MAX_REF_GATES,N_RAYS,iscale
!    1                ,miss,reflect,az_array)
      endif

!     This has been disabled since it apparently uses up the memory when
!     we run on the IBM
      if(.false.)then
!     if(namelist_parms%l_unfold)then
          write(6,*)' read_data_88d: v_nyquist_tilt = ',v_nyquist_tilt       
          if(v_nyquist_tilt .ne. r_missing_data .and. 
     1       v_nyquist_tilt .ne. 0.                   )then
              write(6,*)' read_data_88d: calling unfold...'       
!             call unfold(MAX_VEL_GATES,N_RAYS,velocity
!    1                   ,az_array,v_nyquist_tilt,r_missing_data)       
              v_nyquist_ray = r_missing_data ! prevent further dealiasing
              v_nyquist_tilt = r_missing_data ! prevent further dealiasing
          else
              write(6,*)' Nyquist is missing - skipping unfold call'
          endif
      endif

      RETURN

  998 write(6,825)
  825 format(' Bad Status in Readdata')
      RETURN
      END
