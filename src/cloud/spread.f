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

      subroutine spread(a,imax,jmax,kmax,i,j,k,b)
c this subroutine currently inserts the observation into the
c three dimensional array at one point location.
      dimension a(imax,jmax,kmax)

      if(.true.)return

      a(i,j,k) = b
c     if(b .ne. 1.00)write(6,101)i,j,k,b
101   format(' Spread ',3i3,f8.2)

      return
      end

      subroutine spread2(a_array,b_array,i_array,j_array,n,max_array,kma
     1x
     1                          ,i,j,k,a,b)

c     This subroutine inserts the cloud sounding into the analysis arrays
c     at one point location.

      real a_array(max_array,kmax)
      real b_array(max_array,kmax)
      integer i_array(max_array)
      integer j_array(max_array)

      a_array(n,k) = a
      b_array(n,k) = b
      i_array(n) = i
      j_array(n) = j

      return
      end
