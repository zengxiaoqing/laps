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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine ghbry (i4time,p_3d,pb,sfc_t,lt1dat,htby,
     1     ii,jj,kk,istatus)
    
c     routine to return boundary layer top in pressure units
     
      implicit none

c     input parameters
      
      integer i4time
      integer ii,jj,kk,istatus
      real p_3d (ii,jj,kk)
      real pb (ii,jj)           ! surface station pressure (hPa)
      real sfc_t(ii,jj)         ! surface station temperature (k)
      real lt1dat(ii,jj,kk)
      real htby (ii,jj)
      
      
c     dynamic allocation 
      
      real n2(ii,jj,kk)         ! brunt-vaisala frequency
      real theta (ii,jj,kk)     ! potential temperature (k)
      
c     internal variables

      integer i,j,k
      real x1,x2,y1,y2          ! brackets for zero brunt-vaisala freq. interp
      real g                    ! gravity
      real r                    ! gas constant of dry air
      real kkk                  ! (r/cp) where cp is the specific heat
      real pdepth (ii,jj)

      save g,r,kkk              ! inserted save to be safe, probably not needed

      data g / 980.665/         !  cm/s**2  (cgs units)
      data r / 2.8704e6/        !  erg/g/k  (cgs units)
      data kkk / 0.2857/        ! (dimensionless)

c     note: in the computation of n2 the units of pressure cancel since
c     the pressure occurs both in the numerator and denominator....therefore
c     mb are sufficient units for pressure and the relationship
c     1  mb = 1000 cgs units does not have to be applied.
      
c     ------------begin exe last revised 10/26/99 db
      
c     put surface pressure in as first guess of boundary level top
c     matrix algebra assignment

      htby = pb 

c     alternate method of boundary layer analysis (~ from RUC), 
c     look for the level at which the upper level theta exceeds the surface
c     theta.  This does not use virtual temperature at this time to reduce
c     changes to code.

c     compute theta for all levels, reference pressure (p_naught) is pb or 
c     surface pressure at that location.

      do k = 1,kk
         do j = 1, jj
            do i = 1,ii

               theta(i,j,k) = lt1dat(i,j,k) *
     1              ( pb(i,j) / p_3d(i,j,k) )**kkk

            enddo
         enddo
      enddo

c     now seek out the level at which the surface temperature is greater than 
c     theta from the bottom up

      do i = 1,ii
         do j = 1,jj
            do k = 1,kk

c     note that 3K are added here in this test to the sfc temperature 
c     to give some flexibility to the logic due to roundoff error and
c     surface temperature measurement uncertainty
               if( theta (i,j,k) .ge. sfc_t(i,j)+3.0 ) then
                  htby(i,j) = p_3d(i,j,k)
                  htby(i,j) = min (htby(i,j), pb(i,j) )
c     overwrite this entry with interpolated (linear density) value
c     note that the level (k) is at this time the upper bound of the
c     interpolated computation, the interpolation begins to test whether
c     the hight of the boundary is equal to the surface pressure.  if
c     this is the case, then things will be left alone
                  if(pb(i,j) .eq. htby(i,j) .or. k .eq. 1) then ! continue
                     continue
                  else! compute the interpolated value
                     call interp (sfc_t(i,j),theta(i,j,k-1),
     1                    theta(i,j,k),log(p_3d(i,j,k-1)),
     1                    log(p_3d(i,j,k)),htby(i,j))
                     htby(i,j) = exp (htby(i,j))
                  endif
                  go to 432     !skip out of this loop
               endif
            enddo               ! k
 432        continue            ! height has been found here.
         enddo                  ! j
      enddo                     ! i



c     SKIP THIS SECTION OF OLD CODE
c     section is for original boundary layer analysis (currently skipped over)
c      
c      do j = 1,jj
c         do i = 1,ii
c            do k = 2,kk-1  !note that n2(,,k) end points are not filled.
c               
c               n2(i,j,k) = (g/lt1dat(i,j,k) )**2 * (p_3d(i,j,k)/r) *
c     1              (
c     1              kkk*lt1dat(i,j,k)/p_3d(i,j,k)
c     1              -  ( lt1dat(i,j,k+1)-lt1dat(i,j,k-1) )
c     1              / ( p_3d(i,j,k+1) - p_3d(i,j,k-1) )
c     1              )
c               
c               
c            enddo               ! k
c         enddo
c      enddo
c 
cc     now decide where the "height of the boundary layer" is
c
c      do j = 1,jj
c         do i = 1,ii
c
c            y2 = 2.             ! number greater than zero for test below
c
cc     search upward for n2(k+1) being negative, then interpolate there
c
c            do k = 2,kk-2       ! avoid using n2(,,k) endpoints, (not filled)
c                  x1 = p_3d(i,j,k)
c                  y1 = n2(i,j,k)
c
c               if (pb(i,j) .ge. p_3d(i,j,k+1) 
c     1              .and. n2(i,j,k+1) .lt. 0.0) then
c
c                  x2 = p_3d(i,j,k+1)
c                  y2 = n2(i,j,k+1) 
cc     bail out of loop here to not affect regions above 1st inversion
c                  go to 111
c
c               endif
c               
c            enddo
cc     divert code here to not fall into section 111
c            go to 112
c 111        continue
c            
c            if (y2.lt.0.0)   then !actually found negative level
cc     this test prevents going to the top of the column w/o inversion
cc     should never happen, but this is safeguard.
c
cc     sun machines seem to have a problem when y2 and y1 are very close.
cc     to avoid this situation, such close numbers indicate that the 
cc     top of the boundary is very close to x1 so we assign this
cc     here,  if the difference is large enough, then we use the interp
cc     routine.
c
cc               write (6,*) 'TEMPP ', abs(y2-y1)
c
c               if (abs(y2-y1) .le. 1.e-6) then
c                  htby(i,j) = x1
c               else
c               
cc     interpolate in height space
c                  call interp( 0.,y1,y2,log(x1),log(x2),htby(i,j) )
c                  htby(i,j) = exp(htby(i,j))
c
cc     double safeguard on making sure htby is not below ground level.
c                  htby(i,j) = min (htby(i,j),pb(i,j))
c               endif
c
cc     check for runaway adjustment, assign to base value
c               if (htby(i,j).lt.pb(i,j)-400. ) then ! regard as too high
c                  write (6,*) 'i,j,htby(i,j), gt 400mb depth, adjust',
c     1                 ' to surface pressure - 400 mb', i,j,htby(i,j)
c                  htby(i,j) = pb(i,j)-400.
c               endif
c
c            endif
c               
c 112        continue
c            
c         enddo
c         
c      enddo
c      
c     write out the pbl pressure top
C     END OF SKIPPED SECTION


      pdepth = pb-htby
c     put in check for below ground pressures (assign as sfc pressure)
      where (pdepth < 0.0) 
         htby = pb
         pdepth = pb-htby
      end where

      call check_nan2 (htby,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'NaN values in var:htby routine:ghbry.f'
         return
      endif
      
      call gen_bl_file (i4time,htby,ii,jj,istatus)
      
      if (istatus.eq.0) print*, 'Error in gen_bl_file routine'
      
      return
      end
