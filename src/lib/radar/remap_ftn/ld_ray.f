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
      subroutine LD_RAY(i_ray,                  ! Input
     :                  ngates_remap,           ! Input
     :                  r_missing_data,         ! Input
     :                  i_missing_data,         ! Input
     :                  velocity,               ! Output
     :                  reflect,                ! Output
     :                  istatus)                ! Output
c
c     PURPOSE:
c       Transform a ray common block 88D data to output arrays.
c
c     Changed buffer storage so vel, ref stored as integer*2.
c     This program converts reflectivity and velocity to real*4
c     Keith Brewster, CAPS
c
      implicit none
c
c     Input variables
c
      integer*4 i_ray
      integer*4 ngates_remap
      real*4 r_missing_data
      integer*4 i_missing_data
c
c     Output variables
c
      real*4 velocity(ngates_remap)
      real*4 reflect(ngates_remap)
      integer*4 istatus
c
c     Include files
c
      include 'remap_dims.inc'
      include 'remap_buffer.cmn'
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

      IF (n_vel_gates_cmn .ne. 920) THEN
        write(6,815)
  815   format
     :(' LD_RAY is hard wired for 920 vel gates, error in LD_RAY')
        write(6,820) n_vel_gates_cmn
  820   format(' n_vel_gates_cmn = ',I12)
        istatus = 0
        RETURN
      END IF

      IF (n_ref_gates_cmn .ne. 460) THEN
        write(6,825)
  825   format
     :(' LD_RAY is hard wired for 460 ref gates, error in LD_RAY')
        write(6,830) n_ref_gates_cmn
  830   format(' n_ref_gates_cmn = ',I12)
        istatus = 0
        RETURN
      END IF
c
c     Fill velocity and reflectivity ray.
c      Output gates are 250 m spacing,
c             1840 are used for both reflecivity and velocity.
c      Input gates are 460 gates, 1000m spacing for reflectivity,
c                      920 gates,  250m for velocity.
c      Reflectivity is remapped by gate replication (factor of 4).
c      Velocity is remapped by a 1 to 1 mapping for the 920 gates,
c        output gates 921-1840 are filled by the missing data value.
c

c     Velocity

      ratio_vel = gsp_vel_m_cmn / 250.

      DO i = 1,n_vel_gates_cmn

        IF (b_vel_cmn(i,i_ray) .ne. b_missing_cmn) THEN
          velocity(i) = 0.5*float(b_vel_cmn(i,i_ray))
        ELSE 
          velocity(i) = r_missing_data
        END IF

      ENDDO

c     Reflectivity

      ratio_ref = gsp_ref_m_cmn / 250.

      DO igate_88d = 1,n_ref_gates_cmn

        i = nint(float(igate_88d) * ratio_ref)

        if(i .le. ngates_remap)then
            IF (b_ref_cmn(igate_88d,i_ray) .ne. b_missing_cmn) THEN
              reflect(i) = 0.5*float(b_ref_cmn(igate_88d,i_ray))
            ELSE
              reflect(i) = r_missing_data
            END IF

        endif ! i

      ENDDO
c
      istatus = 1
      RETURN
      END
