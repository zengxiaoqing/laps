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
     1     ii,jj,kk,istatus)
    
c     routine to return boundary layer depth in meters
c     currently this is fixed at 100 meters.  logically, this
c     should be determined by some type of model
     
      implicit none

c     input parameters
      
      integer i4time
      integer ii,jj,kk,istatus
      real plevel (kk)
      real pb (ii,jj)           ! surface station pressure
      real lt1dat(ii,jj,kk)
      real htby (ii,jj)
      
c     dynamic allocation 
      
      real n2(ii,jj,kk)         ! brunt-vaisala frequency
      
c     internal variables

      integer i,j,k
      real x1,x2,y1,y2          ! brackets for zero brunt-vaisala freq. interp
      real g                    ! gravity
      real r                    ! gas constant of dry air
      real kkk                  ! (r/cp) where cp is the specific heat
      save g,r,kkk              ! inserted save to be safe, probably not needed

      data g / 980.665/         !  cm/s**2  (cgs units)
      data r / 2.8704e6/        !  erg/g/k  (cgs units)
      data kkk / 0.2857/        ! (dimensionless)

c     note: in the computation of n2 the units of pressure cancel since
c     the pressure occurs both in the numerator and denominator....therefore
c     mb are sufficient units for pressure and the relationship
c     1  mb = 1000 cgs units does not have to be applied.
      

c     ------------begin exe last revised 10/26/99 db
      
      do j = 1,jj
         do i = 1,ii
            
            htby(i,j) = pb(i,j) !put surface pressure in as first guess of
c     boundary level
            
         enddo
      enddo

c     section is for better boundary layer analysis
      
      do j = 1,jj
         do i = 1,ii
            do k = 2,kk-1
               
               n2(i,j,k) = (g/lt1dat(i,j,k) )**2 * (plevel(k)/r) *
     1              (
     1              kkk*lt1dat(i,j,k)/plevel(k)
     1              -  ( lt1dat(i,j,k+1)-lt1dat(i,j,k-1) )
     1              / ( plevel(k+1) - plevel(k-1) )
     1              )
               
               
            enddo               ! k
         enddo
      enddo
 
c     now decide where the "height of the boundary layer" is

      do j = 1,jj
         do i = 1,ii

            y2 = 2.             ! number greater than zero for test below

c     search upward for n2(k+1) being negative, then interpolate there

            do k = 1,kk-1
                  x1 = plevel(k)
                  y1 = n2(i,j,k)

               if (pb(i,j) .ge. plevel(k+1) 
     1              .and. n2(i,j,k+1) .lt. 0.0) then

                  x2 = plevel(k+1)
                  y2 = n2(i,j,k+1)
c     bail out of loop here to not affect regions above 1st inversion
                  go to 111

               endif
               
            enddo
c     divert code here to not fall into section 111
            go to 112
 111        continue
            
            if (y2.lt.0.0)   then !actually found negative level
c     this test prevents going to the top of the column w/o inversion
c     should never happen, but this is safeguard.
               
c     interpolate in height space
               call interp( 0.,y1,y2,log(x1),log(x2),htby(i,j) )
               htby(i,j) = exp(htby(i,j))

c     double safeguard on making sure htby is not below ground level.
               htby(i,j) = min (htby(i,j),pb(i,j))

            endif
            
 112        continue
            
         enddo

      enddo
      
c     write out the pbl pressure top
      
      call gen_bl_file (i4time,htby,ii,jj,istatus)
      
      if (istatus.eq.0) print*, 'Error in gen_bl_file routine'
      
      return
      end
