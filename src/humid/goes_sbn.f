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
      subroutine goes_sbn (
     1     sh,                  ! specific humidity g/g
     1     lat,lon,             ! lat and longitude (deg)
     1     i4time,              !i4time of run (seconds)
     1     p,                   !pressure hpa (laps vert grid)
     1     cloud,               !cloud array
     1     t,                   ! lt1 (laps 3d temps)
     1     ngoes,               ! goes satellite number
     1     isnd,                ! sounder switch
     1     sat_skip,            ! normally 1 for full resolution
     1     ii,jj,kk             ! grid dimensions
     1     )

 
c   By inclusion of the goes_sbn data into the laps moisture analysis, an
c   improvement in upper level moisture (above 500 mb) can be anticipated to be
c   about 70%.  Current research is pursuing using the satellite data in other
c   levels and other variables such as temperature. 
c
c   19 October, 1999, Tuesday
c   
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

c   
      implicit none
      save
      include 'Constants.inc'
      include 'grid_fname.cmn'

c     include 'lapsparms.for'

c     parameter list variables

      integer ii,jj,kk
      real sh(ii,jj,kk)
      real lat(ii,jj),lon(ii,jj)
      integer i4time
      real t(ii,jj,kk),p(kk)
      real cloud(ii,jj,kk)
      integer ngoes
      integer isnd
      integer sat_skip


c internal variables

      integer istatus
      integer i4time_sat
      integer i,j,k,k2
      real local_model_p(40)
      real dummy
      data local_model_p/.1,.2,.5,1.,1.5,2.,3.,4.,5.,7.,10.,15.,
     1     20.,25.,30.,
     1     50.,60.,70.,85.,100.,115.,135.,150.,200.,250.,300.,350.,400.,
     1     430.,475.,500.,570.,620.,670.,700.,780.,850.,920.,950.,1000./
      integer k500, k700

c climate model variables
      integer*4 julian_day
      real
     1     standard_press(40),
     1     tempertur_guess(40),
     1     mixratio_guess(40)
      real rmd
      integer n_snd_ch
      parameter (n_snd_ch = 22)
      integer kanch(7)
      data kanch /10,8,7,11,16,6,12/


c dynamic dependent variables

      real ch3(ii,jj),ch4(ii,jj),ch5(ii,jj)
      real mr(ii,jj,kk)
      real t_l(kk,ii,jj), mr_l (kk,ii,jj)

      real model_t(40,ii,jj), model_mr(40,ii,jj)


c forward model variarles

c new optran variables
       real tbest(n_snd_ch)


c old gimtau.f variables
      real radiance(ii,jj,18),tskin(ii,jj),psfc(ii,jj),
     1  theta(ii,jj),
     1  ozo(40),gimrad,tau(40)
      real emiss
      integer kan,lsfc(ii,jj)
      real model_p(40)
      real t_fm(40),w_fm(40),ozo_fm(40)
      common/atmos/model_p,t_fm,w_fm,ozo_fm
      real btemp(ii,jj,18),britgo,plango
      real zenith               ! function call
      real pi, d2r

c       powell specific arrays
      real x(3)
      real xi(3,3)
      real ftol,fret
      integer iter(ii,jj)
      real func         ! function typing for cost function
      external func

c   optran specific arrays for powell function calling

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
      real bias_correction ! function
 

      common /cost_optran/radiance_ob, cost_kk, cost_p, cost_t_l,
     1    cost_mr_l, cost_tskin, cost_psfc, cost_julian_day, cost_lat,
     1    cost_theta, cost_isnd

      integer goes_number

      common /sat_id/ goes_number

c  analysis of the factor field
      integer pn
      real points(ii*jj,3)
      real data_anal(ii,jj)
      real ave,adev,sdev,var,skew,curt
      real upper_limit, lower_limit

c  cloud variables
      real cld(ii,jj)

c  moisture modified field
      real factor(ii,jj), factor2(ii,jj)

c  get latest filename
      character*256 path

c  laplace solver variables
      integer mask(ii,jj)

c  misc variables
      integer failures
      character*4 blank

      real rads (ii,jj,n_snd_ch)

      character*9 filename1,  filename
      character*9 grid_name

      integer len

c     check sat_skip for zero, if zero skip routine

      if (sat_skip.le.0) then
         write (6,*) 'sat_skip parameter <= 0, skipping sat entirely'
         return
      endif
      
c     assign sounder/imager parameter for powell method
      
      cost_isnd = isnd
      
c     assign goes number for common block to make avail where needed for OPTRAN
      
      goes_number = ngoes
      
      
c     assign pressure to global array
      
      do i = 1,40
         model_p(i) = local_model_p (i)
      enddo                     !i
      
c     constants
      
      call get_r_missing_data(rmd, istatus)
      
      pi = acos(-1.0)
      d2r = pi/180.
      blank = '  '
      
c     set laps grid

      call get_laps_config(grid_fnam_common,istatus)
 
c      grid_name = 'nest7grid'
c      call get_laps_config(grid_name,istatus)
      
      do j = 1,jj
         do i = 1,ii
            mask(i,j) = 0
         enddo
      enddo
      
c     get satellite IMAGE radiance data for the laps grid
      
      if (isnd .eq. 0) then     ! seek imager data
         
         write(6,*) 'Attemping moisture analysis with imager'
         
         call get_directory('lvd',path,len)
         
c     install new changes for revised satellite path
         
         
         if (ngoes.eq.8) then
            path = path(1:len)//'goes08/'
            len = len + 7
         elseif (ngoes.eq.9) then
            path = path(1:len)//'goes09/'
            len = len + 7
         elseif (ngoes.eq.10) then
            path = path(1:len)//'goes10/'
            len = len + 7
         endif
         
         
         call get_latest_file (path,i4time,filename1,istatus)
         
         if (istatus.ne.1) return
         
         write (6,*) 'Attempting: ', filename1
c     convert filename to i4time_sat
         call i4time_fname_lp (filename1,i4time_sat,istatus)
         write (6,*) 'Getting satellite radainces (lvd)'
         call read_lvd_3_4_5 (path,i4time_sat,ch3,ch4,ch5,
     1        ii,jj,kk,ngoes,istatus)
         
         if (istatus.ne.1) then
            write(6,*) 'error getting satellite data'
            write(6,*) 'aborting goes_sbn module'
            return
         endif
                  
         write(6,*) ' '
         write(6,*) ' '
         write(6,*) 'Using LVD data from: ', filename1
         write(6,*) ' '
         write(6,*) ' '
         
      endif                     ! get IMAGER data only
      
c     acquire sounder data

      if(isnd.eq.1) then ! get SOUNDER data only

       call rsr (i4time, rads, ii,jj,18,ngoes, istatus)
       if (istatus .ne. 1) then
          write (6,*) 'error obtaining sounder radiances'
          return
       endif

      endif ! only get SOUNDER data

c --------- at this point, the existance of satellite data is established,
c ---------  should be cost effective to continue.


c   set up time for regular laps interval
c   generate filename from 14time for julian day extraction later

        call make_fnam_lp (i4time, filename, istatus)

c   get laps surface temperature
        print*, 'getting surface temperature (lsx)'
        call glst(i4time,tskin,ii,jj,istatus)

        if(istatus.ne.1) then

           write(6,*) ' '
           write(6,*) ' '
           write(6,*) 'Failed to get LSX temp data for forward model'
           write(6,*) ' '
           write(6,*) ' '
           return
 
        endif


c   get laps surface pressure

        print*, 'getting surface pressure (lsx)'
        call glsp(i4time,psfc,ii,jj,istatus)

        if(istatus.ne.1) then

           write(6,*) ' '
           write(6,*) ' '
           write(6,*) 'Failed getting LSX pres for forward model'
           write(6,*) ' '
           write(6,*) ' '
           return

        endif

c     convert pressure to hpa
        do j = 1,jj
           do i = 1,ii
              psfc(i,j) = psfc(i,j)/100.
           enddo
        enddo

c     setup cloud test (cloud array passed in)

        do j = 1,jj
           do i = 1,ii
              cld(i,j) = 0.0
              do k = 1,kk
                 cld(i,j) = max(cld(i,j),cloud(i,j,k))
              enddo
              if(cld(i,j).gt.1.) cld(i,j) = 1.0
              if(cld(i,j).le.0.1) cld(i,j) = 0.0
           enddo
        enddo

        write (6,*) 'Running GOES',ngoes,' forward model OPTRAN vsn'

        do j = 1,jj
           do i = 1,ii
              do k = kk,1,-1

                 if(cloud(i,j,k).ge.1.0) then ! assume cloud top

                    if(p(k).lt.psfc(i,j)) then ! above ground level, assign
                       psfc(i,j) = p(k)
                       tskin(i,j) = t(i,j,k)
                       cld(i,j) = cloud(i,j,k)

                    else

                       print*, 'cloud below ground'

                    endif
                    go to 55
                 endif
              enddo

 55           continue

           enddo
        enddo

c       modify sounding to convert sh to mr and model organization
c       assign 0.0 moisture where there is missing data.

        do i = 1,ii
           do j = 1,jj
              do k = 1,kk

                 if(sh(i,j,k) .ne. -1.e30) then
                    call sh2mr (sh(i,j,k), mr(i,j,k) )
                    mr_l(k,i,j) = mr(i,j,k)
                 else
                    mr_l(k,i,j) = 0.0
                 endif
                 t_l (k,i,j) = t(i,j,k)

              enddo
           enddo
        enddo

        read (filename(3:5),22) julian_day
 22     format (i3)

c     prepare to use forward model functions
c     here use goes 8 for reference (goes 10 not avail)

        call pfcgim (8)

        if(isnd .eq.1) then     ! use sounder data for ch3, ch4, ch5
           do j = 1, jj
              do i = 1, ii
                 if(rads(i,j,10).eq.rmd) then
                    ch3(i,j) = rmd
                 else
                    ch3(i,j) = bias_correction (britgo(rads(i,j,10),10),
     1                   ngoes, 1, 10)
                 endif
                 if(rads(i,j,8).eq.rmd) then
                    ch4(i,j) = rmd
                 else
                    ch4(i,j) = bias_correction (britgo(rads(i,j,8),8),
     1                   ngoes, 1, 8)
                 endif
                 if(rads(i,j,7).eq.rmd) then
                    ch5(i,j) = rmd
                 else
                    ch5(i,j) = bias_correction (britgo(rads(i,j,7),7),
     1                   ngoes, 1, 7)
                 endif
              enddo
           enddo
        endif                   ! sounder used

c  do for each gridpoint

        do j = 1,jj,sat_skip
           do i = 1,ii,sat_skip

c     compute zenith angle for model

              theta(i,j) = zenith(lat(i,j)*d2r,
     1             lon(i,j)*d2r,0.*d2r,-75.*d2r)

c     insert call for OPTRAN for initial comparison with gimtau.f
c     note that optran is configured to return both sounder and imager
c     channels used in this algorithm.

              call ofm ( kk, p, t_l(1,i,j), 
     1             mr_l(1,i,j), tskin(i,j), psfc(i,j),
     1             julian_day, lat(i,j),theta(i,j), tbest) 

              if(isnd.eq.0) then ! IMAGER computation

                 do kan = 1,3

                    btemp(i,j,kan) = tbest (kan+7)

                 enddo          !kan

              endif             ! end IMAGER computation

              if(isnd.eq.1) then ! SOUNDER computation

                 do kan = 1,7

                    btemp(i,j,kan) = tbest(kan)

                 enddo          ! kan

              endif             ! end SOUNDER computation

           enddo                ! j
        enddo                   ! i

c     Execute powell method correction of layer humidity in clear areas

      failures = 0

      do j = 1,jj
         do i = 1,ii

            factor (i,j) = rmd
            factor2(i,j) = rmd

         enddo ! i
      enddo ! j



      do j = 1,jj,sat_skip
         do i = 1,ii,sat_skip

            if (i .eq. 1 .and. j.eq.1) then !first time set
               x(1) = 1.0
               x(2) = 1.5
               x(3) = 0.8
            endif

            do k   = 1,3
               x(k) = 1.0
            enddo

            if (ch3(i,j).eq.rmd) then
               print*, 'missing data in channel 3 abort', i,j

            elseif (ch4(i,j).eq.rmd) then
               print*, 'missing data in channel 4 abort', i,j

            elseif (ch5(i,j).eq.rmd) then
               print*, 'missing data in channel 5 abort', i,j

            else

               if (isnd.eq.1) then
                  do k = 4,7
                     if (rads(i,j,k) .eq. rmd) then
                        print*, 'missing data in sounder channel ',
     1                       k,' index ',i,j
                        go to 145
                     endif
                  enddo
               endif

               continue

c               if( cld(i,j) .eq. 0 .or. cld(i,j).ge.1.) then ! clear
               if( (cld(i,j) .eq. 0 .or. cld(i,j).ge.1.)  
     1              .and.
     1              abs(ch4(i,j)-btemp(i,j,2)).le.1.) then !clear

c   print out the "clear" radiances for 6.7 micron only
c   and compare these to the forward model radiances

                  write(6,32) ' Observed=',ch3(i,j),' Modeled='
     1                 ,btemp(i,j,1),' Diff=',(ch3(i,j)-btemp(i,j,1))
 32               format(1x,a10,f8.3,a9,f8.3,a6,f8.3)
                  write(6,*) ch4(i,j),btemp(i,j,2)
                  write(6,*) ch5(i,j), btemp(i,j,3)

                  do k = 1,3
                  do k2 = 1,3
                    xi(k,k2) = 0.0
                     if(k.eq.k2)  xi(k,k) = -.0001
                  enddo
                  enddo

                  if(isnd.eq.0) then ! USE AS IMAGER DATA, btemps
                     radiance_ob(1) = ch3(i,j)
                     radiance_ob(2) = ch4(i,j)
                     radiance_ob(3) = ch5(i,j)
                  endif

                  if(isnd.eq.1) then ! USE AS SOUNDER DATA, btemps
                     radiance_ob(1) = ch3(i,j)
                     radiance_ob(2) = ch4(i,j)
                     radiance_ob(3) = ch5(i,j)
                     radiance_ob(4) = bias_correction(britgo(
     1                    rads(i,j,kanch(4)),kanch(4)),ngoes,1,kanch(4))
                     radiance_ob(5) = bias_correction(britgo(
     1                    rads(i,j,kanch(5)),kanch(5)),ngoes,1,kanch(5))
                     radiance_ob(6) = bias_correction(britgo(
     1                    rads(i,j,kanch(6)),kanch(6)),ngoes,1,kanch(6))
                     radiance_ob(7) = bias_correction(britgo(
     1                    rads(i,j,kanch(7)),kanch(7)),ngoes,1,kanch(7))
                  endif

c fill powell common block with profile data for optran

                  do k = 1, kk
                     cost_p(k) = p(k)
                     cost_t_l(k) = t_l(k,i,j)
                     cost_mr_l(k) = mr_l (k,i,j)
                  enddo
                  cost_kk = kk
                  cost_tskin = tskin (i,j)
                  cost_psfc = psfc (i,j)
                  cost_julian_day = julian_day
                  cost_lat = lat (i,j)
                  cost_theta = theta (i,j)

                  if(cld(i,j).eq.0.) then
c     don't match low atmosphere (use func, not func3)

                     call powell (x,xi,3,3,ftol,iter(i,j),fret,func)

c     else !clouds... don't match low atmosphere
c     call powell (x,xi,3,3,ftol,iter(i,j),fret,func)
                  endif

                  write(6,33) abs(x(1)), abs(x(2)),abs(x(3)),
     1                 i,j,fret,iter(i,j)
 33               format(3(f7.2,2x),i3,i3,1x,f7.2,i3)



                  if (cld(i,j) .eq. 0. .and. iter(i,j) .lt. 50
     1                 .and. abs(abs(x(1))-1.) .lt. .1 .and.
     1                 iter(i,j) .gt. 1 
     1                 .and.  abs(x(3)).ne.0.0
     1                 .and.  abs(x(2)).ne.0.0) then
                     factor(i,j)  = abs(x(3))
                     factor2(i,j) = abs(x(2))
                  elseif (cld(i,j).gt.0.)then
                     write(6,*) '  .... coordinate rejected, cloudy'
                  else
                     write(6,*) i,j, '  .... coordinate rejected', 
     1                    abs(x(1)),iter(i,j), cld(i,j)
                     
                     failures = failures + 1
                     
                  endif
                  
                  write(6,*) blank
                  
               endif            !end of powell function
               
            endif               ! end of missing data flag test
            
 145        continue            !(placed here for sounder missing data flag test)

         enddo
      enddo
      
      write(6,*) failures,' failures occurred due to layer confusion' 
      write(6,*) '...non-convergence, or clouds'
     
c     modify original lq3 file with new factors for comparison tests.
c     modify lq3 only in clear areas as defined by lc3.
      
c     analyze top level adjustments.
      
      pn = 0
      
      do j = 1,jj
         do i = 1,ii
            mask(i,j) = 0
            data_anal(i,j) = 1.
            if (factor(i,j).ne.rmd ) then
               pn = pn+1
               points(pn,1) = factor(i,j)
               points(pn,2) = i
               points(pn,3) = j
               mask(i,j) = 1
               data_anal(i,j) = factor(i,j)
            endif
         enddo
      enddo
      
c     derive field statistics to determine outliers
      call moment_b (points(1,1),pn,ave,adev,sdev,var,skew,
     1     curt,istatus)
      upper_limit = ave + 3.*sdev
      lower_limit = ave - 3.*sdev
      write (6,*) 
      write (6,*) 
      write (6,*) 'Classify acceptable data' 
      write (6,*) 
      write (6,*) 'acceptable range', lower_limit, upper_limit
      
      do i = 1,pn
         if(points(i,1) .lt. upper_limit .and.
     1        points(i,1) .gt. lower_limit) then
            data_anal(int(points(i,2)),int(points(i,3))) = points (i,1)
            write(6,*) points(i,1), 'assigned'
         else
            write(6,*) points(i,1), 'rejected ******************'
            points(i,2) = 0     ! flag for bad point in prep grid
         endif
      enddo
      
      if (pn.ne.0) then
         
         call prep_grid(ii,jj,data_anal,ii*jj,points,pn,istatus)
         if(istatus.eq.1) then
            call slv_laplc (data_anal,mask,ii,jj)
            call smooth_grid2 (ii,jj,data_anal,1)
            call two_d_stats(ii,jj,data_anal,rmd)
            
         else
            write(6,*) 'not enough data to process, skipping slv_lapc'
         endif
         
      else
         write(6,*) 
     1        'pn = 0,no acceptable data to analyze for adjustment'
         return
         
      endif
      
      
c     find k500 (level at or above 500)
      
      do k = 1, kk
         if (p(k) .le. 500.)then
            k500 = k
            go to 475
         endif
      enddo
 475  continue
      
c     find k700 (level at or above 700)
      do k = 1, kk
         if (p(k) .le. 700.)then
            k700 = k
            go to 476
         endif
      enddo
 476  continue
      
      
c     modify lq3 field  top level
      
      do j = 1,jj
         do i = 1,ii
            do k = k500+1, kk   !between 475 and 100 mb
               
               sh(i,j,k) = sh(i,j,k) * 
     1              ((abs(data_anal(i,j))-1.) * 
     1              (p(k)/500.) +1.)

            enddo
         enddo
      enddo
      
c     analyze second level adjustments.
      
      pn = 0
      
      do j = 1,jj
         do i = 1,ii
            mask(i,j) = 0
            data_anal(i,j) = 1.
            if (factor2(i,j).ne.rmd ) then
               pn = pn+1
               points(pn,1) = factor2(i,j)
               points(pn,2) = i
               points(pn,3) = j
               mask(i,j) = 1
               data_anal(i,j) = factor2(i,j)
            endif
         enddo
      enddo
      
c     derive field statistics to determine outliers
      call moment_b (points(1,1),pn,ave,adev,sdev,var,skew,
     1     curt,istatus)
      upper_limit = ave + 3.*sdev
      lower_limit = ave - 3.*sdev
      write (6,*) 
      write (6,*) 
      write (6,*) 'Classify acceptable data' 
      write (6,*) 
      write (6,*) 'acceptable range', lower_limit, upper_limit
      
      do i = 1,pn
         if(points(i,1) .lt. upper_limit .and.
     1        points(i,1) .gt. lower_limit) then
            data_anal(int(points(i,2)),int(points(i,3))) = points (i,1)
            write(6,*) points(i,1), 'assigned'
         else
            write(6,*) points(i,1), 'rejected ******************'
            points(i,2) = 0     ! flag for bad point in prep grid
         endif
      enddo
      
      if (pn.ne.0) then
         
         call prep_grid(ii,jj,data_anal,ii*jj,points,pn,istatus)
         if(istatus.eq.1) then
            call slv_laplc (data_anal,mask,ii,jj)
            call smooth_grid2 (ii,jj,data_anal,1)
            call two_d_stats(ii,jj,data_anal,rmd)
            
         else
            write(6,*) 'not enough data to process, skipping slv_lapc'
         endif
         
      else
         write(6,*) 
     1        'pn = 0,no acceptable data to analyze for adjustment'
         return
         
      endif
      
c     modify lq3 field  second level
      
      do j = 1,jj
         do i = 1,ii
            do k = k700,k500    !between 700 and 500 mb
               
               sh(i,j,k) = sh(i,j,k) * data_anal(i,j)
               
            enddo
         enddo
      enddo      
      
      
      
      return
      end
  
