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


      subroutine adjust_heights(temp_3d,heights_3d,ht_500_fg
     1                         ,ni,nj,nk,k_ref,istatus)

!     This routine changes the reference for the hydrostatic integration
!     of temperature into height to the model 500mb ht field.
!     The heights normally are referenced using unreduced absolute pressure
!     at the surface (laps var = 'PS').

!     This is currently being used for experimental purposes only.

!     Dec 1994    Steve Albers

      real temp_3d(ni,nj,nk)    ! Input
      real heights_3d(ni,nj,nk) ! Input/Output
      real ht_500_fg(ni,nj)     ! Local

      write(6,*)' Subroutine adjust_heights'

      resid_max = 0.

      do i = 1,ni
      do j = 1,nj

!         Calculate residual = model 500 ht - first guess height
          residual = ht_500_fg(i,j) - heights_3d(i,j,k_ref)

!         Apply residual to column of height values.
          do k = 1,nk
              heights_3d(i,j,k) = heights_3d(i,j,k) + residual
          enddo ! k

          resid_max = max(resid_max,abs(residual))

      enddo ! j
      enddo ! i

      write(6,*)'         Maximum height adjustment (m) = ',resid_max

      istatus = 1
      return
      end
