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
      integer i_tilt
      integer Max_rays
c
c     Output variables
c
      real    az_array(Max_rays)
      integer istatus
c
c     Include files
c
      include 'remap_dims.inc'
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
