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
      subroutine GET_AZIMUTHS_DEG ( i_tilt,    ! Input
     :              Max_rays,        ! Input # rays dimensioned for remapper
     :              az_array,        ! Output
     :              istatus )        ! Output
c
c     PURPOSE:
c       Fill azimuth array with azimuths from the common buffer.
c
c     AUTHOR:
c       Steve Albers, FSL
c     MODIFICATIONS:
c       Clean-up and improve efficiency Keith Brwester, CAPS, June, 1994
c

      implicit none
c
c     Input variables
c
      integer*4 i_tilt
      integer*4 Max_rays
c
c     Output variables
c
      real*4    az_array(Max_rays)
      integer*4 istatus
c
c     Include files
c
      include 'remap_buffer.cmn'
c
c     Misc internal variables
c
      integer i
c
      IF( n_rays_cmn .gt. Max_rays ) THEN
        write(6,805)
 805    format(' n_rays_cmn exceeds Max_rays in GET_AZIMUTHS_DEG')
        write(6,810) n_rays_cmn,Max_rays
 810    format(2I12)
        istatus = 0
        RETURN
      END IF
c
c     Transfer common data to output array
c
      DO 100 i = 1,n_rays_cmn
        az_array(i) = azim_cmn(i)
 100  CONTINUE
      istatus = 1
c
      RETURN
      END
