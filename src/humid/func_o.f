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
      real function func (x)


c   This routine interfaces GOES 8/10 satellite broadcast network data (and
c   local GVAR data) to the LAPS moisture analysis.  In 1999, this routine
c   was modified from an earlier version that used the University of
c   Wisconsin -- Madison's forward model to a new model developed at
c   NESDIS.  OPTRAN (optical transmittance) forward model was developed by
c   Thomas Kleespies (NESDIS) and questions about this model should be
c   directed to him.  Forecast Systems Laboratory does not in any way
c   guarantee the validity of OPTRAN and distributes this software on an
c   as-is basis.  MOREOVER, FSL HAS PERMISSION TO DISTRIBUTE OPTRAN AS PART
c   OF LAPS TO FEDERAL AGENCIES.  NON-FEDERAL ENTITIES NEED TO INQUIRE WITH
c   NESDIS TO ESTABLISH THEIR RIGHTS AND OBLIGATIONS WITH REGARD TO OPTRAN.
c   
c   The version of OPTRAN with which this software is used, has been
c   modified by FSL to include both sounder and imager channels for a
c   particular satellite in one call to the routine.  Thus a user only need
c   to setup OPTRAN for a particular satellite.  After doing such, either
c   the imager or sounding instrument can be used with the software without
c   further recompilation.

      implicit none
      save
      include 'Constants.inc'
      
c     parameter variables
      
      real x(3)
      
c     optran specific arrays for powell function calling
      
      real radiance_ob (Nchan)
      integer cost_kk
      real cost_p(Nlevel)
      real cost_t_l(Nlevel)
      real cost_mr_l(Nlevel)
      real cost_tskin
      real cost_psfc
      integer cost_julian_day
      real cost_lat
      real cost_theta
      integer cost_isnd

      common /cost_optran/radiance_ob, cost_kk, cost_p, cost_t_l,
     1     cost_mr_l, cost_tskin, cost_psfc, cost_julian_day, cost_lat,
     1     cost_theta, cost_isnd

c     local analogs to common block variables for input to parameters

      integer kk
      real p(Nlevel)
      real t_l(Nlevel)
      real mr_l(Nlevel)
      real tskin
      real psfc
      integer julian_day
      real lat
      real theta

c     local variables

      integer i,j
      integer kan(18)
      real w(Nlevel)
      real tbest(nchan)
      logical first_time
      data first_time /.true./
      integer lvl500, lvl700, lvl100
      real var_weights(7)       ! weights computing in func
      
c     externals
      
      real plango
      
c     code
      
      if (first_time) then
         first_time  = .false.
         
c     set up for sounder if needed instead of imager
         
         if (cost_isnd.eq.1) then ! sounder radiances used
            
            kan(1) = 10         ! 7.4 like ch3
            kan(2) = 8          ! 11.02 like ch4
            kan(3) = 7          ! 12.02 like ch5
            kan(4) = 11
            kan(5) = 16
            kan(6) = 6
            kan(7) = 12
           
            var_weights(1) = .0022
            var_weights(2) = 0.
            var_weights(3) = 0.
            var_weights(4) = .0034
            var_weights(5) = 0.
            var_weights(6) = 0.
            var_weights(7) = 0.0036
            
         else                   ! imager radiances used

            var_weights(1) = .0022
            var_weights(2) = 0.
            var_weights(3) = 0.
            
         endif
         
c     determine layer levels 
         
         do i = cost_kk,1, -1
            if (cost_p(i) .le. 100.) lvl100=i
            if (cost_p(i) .le. 700.) lvl700=i
            if (cost_p(i) .le. 500.) lvl500=i
         enddo
         
      endif                     ! first_time
      
c     fill local variables from common block for parameter list in ofm
      
      do i = 1, nlevel
         p(i) = cost_p(i)
         t_l(i) = cost_t_l(i)
         mr_l(i) = cost_mr_l(i)
      enddo
      kk = cost_kk
      tskin = cost_tskin
      psfc = cost_psfc
      julian_day = cost_julian_day
      lat = cost_lat
      theta = cost_theta
      
c     modify mixing ratio per predifined pressure layers.
      
      do i = 1,cost_kk
         
         if(i.lt. lvl700) then  ! sfc to 780
            mr_l(i) = abs(x(1)) * cost_mr_l(i)
         elseif (i.ge.lvl700 .and. i.le. lvl500) then ! 700 to 500
            mr_l(i) = abs(x(2)) * cost_mr_l(i)
         elseif (i.gt.lvl500 .and. i.le. lvl100) then ! between 475 and 100

c     the corresponding change must also be made in GOES_SBN.f where
c     this information is applied.
            mr_l(i) = 
     1           ((abs(x(3))-1.)*(p(i)/500.) + 1.) *
     1           cost_mr_l(i)
         else
            mr_l(i) =  cost_mr_l(i)
         endif
         
      enddo
      
c     perform forward model computation for radiance
c     here is the main call to optran in this whole module, the real
c     workhorse part of the code.
      
      call ofm ( kk, p, t_l, 
     1     mr_l, tskin, psfc,
     1     julian_day, lat, theta, tbest) 
      
c     compute cost function
      
      func = 0.0
      
      if (cost_isnd.eq.1) then  ! SOUNDER radiances used
         
         do j = 1,7             ! radiance _ob(1-7) is sounder btemp
            func = func + var_weights(j)*( radiance_ob(j) -
     1           tbest(j) )**2/2.
         enddo
         
      else                      ! IMAGER situation (only 3 channels)

         do j = 1,3             ! radiance _ob(1-3) is imager btemp
            func = func + var_weights(j)*( radiance_ob(j) -
     1           tbest(j+7) )**2/2.
         enddo

      endif

c     stability cost is identical for both imager and sounder

      do j = 1,3
         func = func + .1 * ((x(j) - 1.)**2 )
      enddo

      return
      end
