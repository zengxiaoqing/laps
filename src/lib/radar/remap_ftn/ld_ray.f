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
      subroutine LD_RAY(i_ray,                  ! Input
     :                  ngates_remap,           ! Input
     :                  r_missing_data,         ! Input
!    :                  i_missing_data,         ! Input
     :                  velocity,               ! Output
     :                  reflect,                ! Output
     :                  istatus)                ! Output
c
c     PURPOSE:
c       Transform a ray in the common block to output arrays that correspond
c       to the lookup tables. In the output arrays, the reflectivity and
c       velocity have been remapped to the same gate spacing.
c
c     Changed buffer storage so vel, ref stored as integer*2.
c     This program converts reflectivity and velocity to real
c     Keith Brewster, CAPS
c
      implicit none
c
c     Input variables
c
      integer i_ray
      integer ngates_remap
      real r_missing_data
c
c     Output variables
c
      real velocity(ngates_remap)
      real reflect(ngates_remap)
      integer istatus
c
c     Include files
c
      include 'remap_dims.inc'
      include 'remap_buffer.cmn'
      include 'remap_constants.dat'
c
c     Misc. Local Variables
c
      integer i,igate_88d
      real ratio_ref, ratio_vel
c
c     Make sure the array dimensions are consistent
c
      IF (ngates_remap .ne. 1840) THEN
        write(6,805)
  805   format
     :(' LD_RAY is hard wired for 1840 gates, error in LD_RAY')
        write(6,810) ngates_remap
  810   format(' ngates_remap = ',I12)
        istatus = 0
        RETURN
      END IF

      IF (n_vel_gates_cmn .ne. 920 .and. i_ray .eq. 1) THEN
        write(6,*)' NOTE: LD_RAY using non-standard # of vel gates'
      END IF

!     IF (n_vel_gates_cmn .ne. MAX_VEL_GATES) THEN
!       write(6,*)
!    1       ' ERROR: LD_RAY is expecting n_vel_gates_cmn=MAX_VEL_GATES'       
!       write(6,820) n_vel_gates_cmn,MAX_VEL_GATES
! 820   format(' n_vel_gates_cmn/MAX_VEL_GATES = ',2I12)
!       istatus = 0
!       RETURN
!     END IF

      IF (n_ref_gates_cmn .ne. 460 .and. i_ray .eq. 1) THEN
        write(6,*)' NOTE: LD_RAY using non-standard # of ref gates'
      END IF

!     IF (n_ref_gates_cmn .ne. MAX_REF_GATES) THEN
!       write(6,*)
!    1       ' ERROR: LD_RAY is expecting n_ref_gates_cmn=MAX_REF_GATES'       
!       write(6,830) n_ref_gates_cmn,MAX_REF_GATES
! 830   format(' n_ref_gates_cmn/MAX_REF_GATES = ',2I12)
!       istatus = 0
!       RETURN
!     END IF
c
c     Fill velocity and reflectivity ray.
c      Output gates are 250 m spacing,
c             1840 are used for both reflecivity and velocity.
c
c      For a typical NEXRAD radar....
c     
c      Input gates are as follows:
c             460 gates, 1000m spacing for reflectivity,
c             920 gates,  250m for velocity.
c      Reflectivity is remapped by gate replication (factor of 4).
c      Velocity is remapped by a 1 to 1 mapping for the 920 gates,
c        output gates 921-1840 are filled by the missing data value.
c
c      Note that other input gate spacings may not be fully supported yet.
c

c     Velocity

      ratio_vel = gsp_vel_m_cmn / GATE_SPACING_M
      velocity = r_missing_data           ! Initialize array all at once

      DO i = 1,n_vel_gates_cmn

        IF (b_vel_cmn(i,i_ray) .ne. b_missing_cmn) THEN
          velocity(i) = 0.5*float(b_vel_cmn(i,i_ray))
!       ELSE 
!         velocity(i) = r_missing_data
        END IF

      ENDDO

c     Reflectivity

      ratio_ref = gsp_ref_m_cmn / GATE_SPACING_M
      reflect = r_missing_data            ! Initialize array all at once

      DO igate_88d = 1,n_ref_gates_cmn

        i = nint(float(igate_88d) * ratio_ref)

        if(i .le. ngates_remap)then
            IF (b_ref_cmn(igate_88d,i_ray) .ne. b_missing_cmn) THEN
              reflect(i) = 0.5*float(b_ref_cmn(igate_88d,i_ray))
!           ELSE
!             reflect(i) = r_missing_data
            END IF

        endif ! i

      ENDDO
c
      if(i_ray .eq. 1)then
          write(6,*)' LD_RAY: ratio_ref = ',ratio_ref
     1                     ,' ratio_vel = ',ratio_vel
      endif

      istatus = 1
      RETURN
      END
