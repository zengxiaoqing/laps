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
      subroutine tiros (
     1  sh,                 ! specific humidity g/g
     1  lat,lon,            ! lat and longitude (deg)
     1  i4time,             ! i4time of run (seconds)
     1  p,                  ! pressure hpa (laps vert grid)
     1  cloud,              ! cloud array
     1  t,                  ! lt1 (laps 3d temps)
     1  ntiros,             ! satellite number
     1  ii,jj,kk            ! grid dimensions
     1  )

c   The module tiros.f is the current interface to use
c   TOVS radiances in LAPS through variational methods. 
c   This algorithm is now using the ECMWF RATAOV forward
c   model for its radiance computations and the
c   computation of the Jacobian matrix used in the
c   solution of the 1-D variational adjustment of the
c   atmospheric profile.
c   
c   Since the RTATOV software is proprietary code, it is
c   not offered with LAPS.  Instead a signed agreement
c   is required with ECMWF to establish the rights to
c   use the code.  Once this is in hand, FSL will
c   provide and interface to RTATOV.
c   
c   Author:  Dan Birkenheuer
c   Date:    April 1997



      implicit none

      include 'lapsparms.for'

c  parameter list variables
      integer ii,jj,kk
      real    sh(ii,jj,kk)
      real    lat(ii,jj),lon(ii,jj)
      integer i4time
      real    t(ii,jj,kk),p(kk)
      real    cloud(ii,jj,kk)
      integer ntiros


c internal variables
      integer istatus
      integer i,j,k


c climate model variables
      integer*4 julian_day
      real standard_press(40)
      real tempertur_guess(40)
      real mixratio_guess(40)


c missing data flag
      real rmd

c RTATOV variables
      integer nlev
      integer nchan
      parameter (nlev = 67)    ! background profile dim
      parameter (nchan = 40)   ! channels computed
      real ulad (nlev)         ! background profile
      real radiance(nchan)     ! computed radiances
      real b_temps (nchan)     ! computed brightness temps
      real xktbb (nlev,nchan)  ! adjoint of Jacobian matrix
      real inv_k (nlev,nchan)  ! inverse multiplier matrix
      real mx (nlev)           ! modified profile
      real ym (nchan)          ! radiance data (btemp K)

c dynamic dependent variables
      real lt1(nz_l,1)  ! one location temp
      real lq3(nz_l,1)  ! one location q
      real psfc (nx_l, ny_l)
      real t_surf (nx_l, ny_l)
      real td_surf (nx_l, ny_l)


c  cloud variables
      real cld(nx_l,ny_l)





c      data    model_p/.1,.2,.5,1.,1.5,2.,3.,4.,5.,7.,10.,15.,
c     1 20.,25.,30.,
c     1 50.,60.,70.,85.,100.,115.,135.,150.,200.,250.,300.,350.,400.,
c     1 430.,475.,500.,570.,620.,670.,700.,780.,850.,920.,950.,1000./




      character*9 filename
      character*9 grid_name

c     this is a place holder routine for the real routine

      write(6,*) 'Tiros code activated but not present'
      write(6,*) 'Use of RTATOV requires license from ECMWF'
      write(6,*) 'Stopping lq3 run, invalid moisture_switch.txt file'

      stop

      return
      end

