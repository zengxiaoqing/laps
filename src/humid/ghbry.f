cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps
cdis
cdis    this software and its documentation are in the public domain and
cdis    are furnished "as is."  the united states government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  they assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  if significant modifications or enhancements
cdis    are made to this software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
        subroutine ghbry (i4time,plevel,pb,lt1dat,htby,
     1                      ii,jj,kk,istatus)


c       $log: ghbry.for,v $
c revision 1.2  1995/09/13  21:35:21  birk
c added disclaimer to files
c
c revision 1.1  1994/04/25  15:12:05  birk
c initial revision
c

c       routine to return boundary layer depth in meters
c       currently this is fixed at 100 meters.  logically, this
c       should be determined by some type of model


        implicit none

c        include 'parmtrs.inc'

c input parameters

      integer i4time
      integer ii,jj,kk,istatus
      real plevel (kk)
      real pb (ii,jj)  ! surface station pressure
      real lt1dat(ii,jj,kk)
      real htby (ii,jj)

c dynamic allocation 

        real*4 n2(ii,jj,kk) ! brunt-vaisala frequency

c internal variables

        real*4 x1,x2,y1,y2 ! brackets for zero brunt-vaisala freq. interp
        real*4 g ! gravity
        real*4 r ! gas constant of dry air
        real*4 kkk ! (r/cp) where cp is the specific heat
        real*4 y1min

        data g / 980.665/   !  cm/s**2  (cgs units)
        data r / 2.8704e6/  !  erg/g/k  (cgs units)
        data kkk / 0.2857/    ! (dimensionless)

        logical*4 test  ! to test whether we are above the surface or not

c       note: in the computation of n2 the units of pressure cancel since
c       the pressure occurs both in the numerator and denominator....therefore
c       mb are sufficient units for pressure and the relationship
c     1  mb = 1000 cgs units does not have to be applied.

        integer*4 i,j,k

        do j = 1,jj
        do i = 1,ii

        htby(i,j) = pb(i,j)  !put surface pressure in as first guess of
c               boundary level

        enddo
        enddo


c       istatus = 1
c       if(1.eq.1) then
c       return
c       endif

c       code here currently under test march 13, 1992  dan b.
c       section is for better boundary layer analysis

        do j = 1,jj
        do i = 1,ii
        do k = 2,kk-1

        n2(i,j,k) = (g/lt1dat(i,j,k) )**2 * (plevel(k)/r) *
     1  (
     1  kkk*lt1dat(i,j,k)/plevel(k)
     1  -  ( lt1dat(i,j,k+1)-lt1dat(i,j,k-1) )
     1  / ( plevel(k+1) - plevel(k-1) )
     1  )


        enddo ! k
        enddo
        enddo




c        now decide where the "height of the boundary layer" is


        y1min = 0.0
        do j = 1,jj
        do i = 1,ii

        test = .false.
        do k = 1,kk

        if (pb(i,j) .ge. plevel(k) .and. n2(i,j,k) .lt. 0.0) then
        test = .true.
        x1 = plevel(k)
        y1 = n2(i,j,k)
        y1min = min(y1,y1min)
        endif

        if (test .and. pb(i,j) .ge. plevel(k) .and. n2(i,j,k) .ge. 0.0) 
     1then
        test = .false.
        x2 = plevel(k)
        y2 = n2(i,j,k)
c       bail out of loop here to not affect regions above 1st inversion
        go to 111
        endif


        enddo
c       divert code here to not fall into section 111
        go to 112
111     continue


        if (y1.ne.y2)   then

c       interpolate in height space
                call interp( 0.,y1,y2,log(x1),log(x2),htby(i,j) )
                htby(i,j) = exp(htby(i,j))


        endif

112     continue

        enddo

c       print*, j,y1min,y1
c       assign y1 zero for monitoring purposes only when using above type
c       statement
        y1 = 0.0
        enddo

c       write out the pbl pressure top

        call gen_bl_file (i4time,htby,ii,jj,istatus)

        if (istatus.eq.0) print*, 'error in gen_bl_file routine'



        return
        end
