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

      real*4 a_array(max_array,kmax)
      real*4 b_array(max_array,kmax)
      integer*4 i_array(max_array)
      integer*4 j_array(max_array)

      a_array(n,k) = a
      b_array(n,k) = b
      i_array(n) = i
      j_array(n) = j

      return
      end
