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
     1  sh,                 ! specific humidity g/g
     1  lat,lon,            ! lat and longitude (deg)
     1  i4time,             !i4time of run (seconds)
     1  p,                  !pressure hpa (laps vert grid)
     1  cloud,              !cloud array
     1  t,                  ! lt1 (laps 3d temps)
     1  ngoes,              ! goes satellite number
     1  isnd,               ! sounder switch
     1  ii,jj,kk            ! grid dimensions
     1  )


c   The routine goes_sbn is the link to goes 8/9 satellite broadcast network
c   (sbn) data input to the laps moisture analysis.  This module performs one
c   major task.  reading in the "background" or input specific humidity (sh)
c   data along with the necessary data to run the forward model, and the goes
c   channel radiances from the laps lvd files, it will modify the laps sh field
c   to best agree with the observed satellite radiances in three of the 4 IR
c   goes channels.  The main consideration here was to run this step using GOES
c   image data which are available to the local forecast office through SBN
c   (Satellite Broadcast Network).
c   
c   Navigation of sbn data into laps coordinates is handled by the lvd process.
c   (creator John Smart and Dan Birkenheuer).  This step, using the satellite
c   data in the moisture code was put together by Dan Birkenheuer.  Components
c   of this module are adapted from the university of Wisconsin -- Madison,
c   specifically the goes forward model.  This version of the goes forward
c   model is the latest as of fall 1995.  Coefficients were transformed from
c   ibm specific format to ASCII making them fully portable for a heterogeneous
c   UNIX environment.  The Powell method is used for solving the variational
c   analysis herein.  It may not be the most elegant method, but has proven to
c   work.
c   
c   By inclusion of the goes_sbn data into the laps moisture analysis, an
c   improvement in upper level moisture (above 500 mb) can be anticipated to be
c   about 70%.  Current research is pursuing using the satellite data in other
c   levels and other variables such as temperature. 
c   
c             ========== BEGIN MAJOR REVISION STATEMENT  ===========
c      
c                 Remarks about the revision made March 18, 1997
c                                 Dan Birkenheuer
c                        CIRA Forecast Systems Laboratory
c   
c   The latest set of changes to this module essentially is a modification to
c   code to enable it to do something for which it was not originally intended.
c   (wecome to the real world of research code)  Hence, the code is not exactly
c   suited to this mod.  Here is a bit of history.
c   
c   Initially module GOES_SBN was designed to utilize AWIPS satellite broadcast
c   network (SBN) GOES imagery and compare imagery radiances with radiances
c   derived using the forward model and through veriational methods, scale the
c   upper level moisture (500-100hpa in LAPS) so the modeled radiances match
c   those measured.  The three radiances used were imager channels 3,4,5 or the
c   6.7, 11., and 12 micron IR channels.
c   
c   The next chapter in this story comes about when it was desired to utilize
c   the same code, but run it using sounder instrument radiances instead of
c   imager radiances.  This is where we are at the moment.  It is hoped that at
c   some point down the road, we can go one step further and use MORE sounder
c   radiances to the problem.  However, that would probably mean a much larger
c   scale change to this code than these minor adjustments.  I envision that
c   would mean working up an entirely new module.  Plus a bit of research which
c   is not included in any budgets at the moment.  So this is what we have.
c   
c   To adapt to sounder data, a new flag was added - ISND - which is now read
c   from moisture.txt in the LAPS static area (the fourth parameter to this
c   ever-growing file).  When ISND is 0 the module should act as it used to
c   processing imager data from possibly a SBN source, (GVAR works too).  If
c   ISND is 1 however, the code will seek out sounder data and substitute
c   sounder channels 10, 8, and 7 for imager channels 3, 4, and 5 respectively.
c   It will then automatically relay this info to the FUNC routine used by the
c   variational adjustment method causing it to use the forward model
c   configured for the sounder.  The module should be sounder capable with this
c   arrangement.
c   
c      ================   END MAJOR REVISION STATEMENT  ===================


c $log: goes_sbn.for,v $
c revision 1.2  1996/08/30  20:46:15  birk
c modified to run under laps
c
c revision 1.1  1996/08/12  19:44:59  birk
c initial revision
c


        implicit none

c        include 'lapsparms.for'

c       parameter list variables

      integer ii,jj,kk
      real sh(ii,jj,kk)
      real lat(ii,jj),lon(ii,jj)
      integer i4time
      real t(ii,jj,kk),p(kk)
      real cloud(ii,jj,kk)
      integer ngoes
      integer isnd


c internal variables

      integer istatus
      integer i4time_sat
      integer i,j,k,k2
      real local_model_p(40)

c climate model variables
      integer*4 julian_day
      real
     1     standard_press(40),
     1     tempertur_guess(40),
     1     mixratio_guess(40)
      real rmd
      integer n_snd_ch
      parameter (n_snd_ch = 18)
      integer kanch(3)
      data kanch /10,8,7/


c dynamic dependent variables

      real ch3(ii,jj),ch4(ii,jj),ch5(ii,jj)
      real mr(ii,jj,kk)
      real t_l(kk,ii,jj), mr_l (kk,ii,jj)

      real model_t(40,ii,jj), model_mr(40,ii,jj)


c forward model variarles
      real radiance(ii,jj,3),tskin(ii,jj),psfc(ii,jj),
     1  theta(ii,jj),
     1  ozo(40),gimrad,tau(40)
      real emiss
      integer kan,lsfc(ii,jj)
      real model_p(40)
      real t_fm(40),w_fm(40),ozo_fm(40)
      common/atmos/model_p,t_fm,w_fm,ozo_fm
      real btemp(ii,jj,3),britgo,plango
      real zenith               ! function call
      real pi, d2r

c       powell specific arrays
      real x(3)
      real xi(3,3)
      real ftol,fret
      integer iter(ii,jj)
      integer ngoes_cost,isnd_cost
      real radiance_ob(3),p_cost(40),t_cost(40),ozo_cost(40),
     1  tskin_cost,psfc_cost,
     1        theta_cost,w_cost(40)
      integer lsfc_cost
      common/cost_var/radiance_ob, p_cost, t_cost,ozo_cost,
     1        tskin_cost, lsfc_cost,psfc_cost,theta_cost,w_cost,
     1        ngoes_cost,isnd_cost
      real func, func3          ! function typing for cost function
      external func,func3

c  analysis of the factor field
      integer pn
      real points(3,ii*jj)
      real data_anal(ii,jj)

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
cdline        integer*4 mlevel(kk)
cdline        character*125 commentline
cdline        character*64 dummy_string
cdline        real research_output (ii,jj,kk)
cdline        real s_btemp(ii,jj,18)  !sounder b_temp
cdline        real s_radiance(ii,jj,18)  ! sounder radiance
cdline        real w_model (39)  ! forward model weighting function
cdline        real p_dm (39)  ! derivate pressures
      real rads (ii,jj,n_snd_ch)

        data local_model_p/.1,.2,.5,1.,1.5,2.,3.,4.,5.,7.,10.,15.,
     1  20.,25.,30.,
     1  50.,60.,70.,85.,100.,115.,135.,150.,200.,250.,300.,350.,400.,
     1  430.,475.,500.,570.,620.,670.,700.,780.,850.,920.,950.,1000./




        character*9 filename1,  filename
        character*9 grid_name

        integer len

c assign satellite number for func routine in powell

           ngoes_cost = ngoes

c assign pressure to global array

           do i = 1,40
              model_p(i) = local_model_p (i)
           enddo !i

c pass imager/sounder information to routine func for powell

           isnd_cost = isnd

c       constants


        call get_r_missing_data(rmd, istatus)

        pi = acos(-1.0)
        d2r = pi/180.
        blank = '  '

cdline              do i = 1,kk
cdline              mlevel(i) = p(i)
cdline              enddo




c       set laps grid                               !test removal later
        grid_name = 'nest7grid'
        call get_laps_config(grid_name,istatus)

        do j = 1,jj
        do i = 1,ii
        mask(i,j) = 0
        enddo
        enddo

c       get satellite IMAGE radiance data for the laps grid

      if (isnd .eq. 0) then ! seek imager data

         call get_directory('lvd',path,len)

         call get_latest_file (path,i4time,filename1,istatus)

        if (istatus.ne.1) return

        write (6,*) 'Attempting: ', filename1
c       convert filename to i4time_sat
        call i4time_fname_lp (filename1,i4time_sat,istatus)
        write (6,*) 'Getting satellite radainces (lvd)'
        call read_lvd_3_4_5 (i4time_sat,ch3,ch4,ch5,
     1      ii,jj,kk,ngoes,istatus)

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

      endif  ! get IMAGER data only

c  acquire sounder data (experimental)

      if(isnd.eq.1) then ! get SOUNDER data only, still experimental


       call rsr (i4time, rads, ii,jj,n_snd_ch,ngoes, istatus)
       if (istatus .ne. 1) then
          write (6,*) 'error obtaining sounder radiances'
          return
       endif



      endif ! only get SOUNDER data , experimental commment out.


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

c  convert pressure to hpa
        do j = 1,jj
        do i = 1,ii
        psfc(i,j) = psfc(i,j)/100.
        enddo
        enddo


c  setup cloud test (cloud array passed in)

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


        write (6,*) 'Running GOES',ngoes,' forward model'


c       now with surface temp, pressure and cloud info we can proceed
c       to process a "radiance surface".  that is, this surface sits
c       at the bottom of the clear column above.  below this we cannot
c       model in the ir.  above this we can use the information to analyze
c       in laps.


c       assume that we have a clear column

        do j = 1,jj
        do i = 1,ii

        call locate(model_p,40,psfc(i,j),lsfc(i,j))

        enddo
        enddo

c       we now have tskin, psfc, lsfc defined for clear column
c       now modify for cloud top.


        do j = 1,jj
           do i = 1,ii

c       print*, ' '

              do k = kk,1,-1

                 if(cloud(i,j,k).ge.1.0) then ! assume cloud top

                    if(p(k).lt.psfc(i,j)) then ! above ground level, assign
                       psfc(i,j) = p(k)
                       tskin(i,j) = t(i,j,k)
                       cld(i,j) = cloud(i,j,k)

                       call locate(model_p,40,psfc(i,j),lsfc(i,j))

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

c  interpolate to forward model coordinates

c  we have laps variables in t_l, and mr_l organized with pressure
c  p (laps p, given), desire in model_p space for the forward model
c

        do k = 20,40            ! for each desired pressure

           call locate(p,kk,model_p(k),k2)

           do i = 1,ii
              do j = 1,jj


                 if (p(k2).eq.model_p(k)) then
                    model_t(k,i,j) = t_l(k2,i,j)
                    model_mr(k,i,j) = mr_l(k2,i,j) * 1000. ! g/kg
                 else
                    call interp (model_p(k),p(k2),p(k2+1),
     1                   t_l(k2,i,j),t_l(k2+1,i,j),model_t(k,i,j))
                    call interp (model_p(k),p(k2),p(k2+1),
     1                   mr_l(k2,i,j)*1000.,mr_l(k2+1,i,j)*1000.,
     1                   model_mr(k,i,j)) ! g/kg
                 endif




              enddo
           enddo
        enddo

c  now we need to fill in the climotology in the upper levels of the
c   sounding in model_p space.



        read (filename(3:5),22) julian_day
 22     format (i3)

        call climate_sm (lat(1,1),julian_day,standard_press,
     1       tempertur_guess,mixratio_guess,istatus)

        do j = 1,jj
           do i = 1,ii
              do k = 1,19


                 model_t(k,i,j) = tempertur_guess(k)
                 model_mr(k,i,j) = mixratio_guess(k) ! g/kg


              enddo
           enddo
        enddo


c  prepare to use forward model functions

        call pfcgim (ngoes)


c now in experimental code, we have both sounder and imager data
c we now test and use the appropriate one.  default d_line code.
c if d_lines are not activated, then we will use IMAGE data at all times

        if(isnd .eq.1) then     ! use sounder data for ch3, ch4, ch5
           do j = 1, jj
              do i = 1, ii
                 if(rads(i,j,10).eq.rmd) then
                    ch3(i,j) = rmd
                 else
                    ch3(i,j) = britgo(rads(i,j,10),10)
                 endif
                 if(rads(i,j,8).eq.rmd) then
                    ch4(i,j) = rmd
                 else
                    ch4(i,j) = britgo(rads(i,j,8),8)
                 endif
                 if(rads(i,j,7).eq.rmd) then
                    ch5(i,j) = rmd
                 else
                    ch5(i,j) = britgo(rads(i,j,7),7)
                 endif
              enddo
           enddo
        endif                   ! sounder used



c  call forward model with these profiles and output radiances

        call pfcgim (ngoes)

c  do for each gridpoint



        do j = 1,jj
           do i = 1,ii


              theta(i,j) = zenith(lat(i,j)*d2r,
     1             lon(i,j)*d2r,0.*d2r,-75.*d2r)
              emiss = .99




c perform forward model computation for radiance

              if(isnd.eq.0) then ! IMAGER computation
                 do kan = 1,3

                    call taugim(model_t(1,i,j),model_mr(1,i,j),ozo,
     1                   theta(i,j),ngoes,kan+22,tau)
                    radiance(i,j,kan) = gimrad(tau,model_t(1,i,j),
     1                   tskin(i,j),
     1                   kan+22,lsfc(i,j),psfc(i,j),emiss)
                    btemp(i,j,kan) = britgo(radiance(i,j,kan),kan+22)

                 enddo          ! kan
              endif             ! IMAGER computation




              if(isnd.eq.1) then ! SOUNDER computation
                 do kan = 1,3

                    call taugim(model_t(1,i,j),model_mr(1,i,j),ozo,
     1                   theta(i,j),ngoes,kanch(kan),tau)
                    radiance(i,j,kan) = gimrad(tau,model_t(1,i,j),
     1                   tskin(i,j),
     1                   kanch(kan),lsfc(i,j),psfc(i,j),emiss)
                    btemp(i,j,kan) = britgo(radiance(i,j,kan),
     1                   kanch(kan))

 

                 enddo          ! kan
              endif             ! SOUNDER computation

          


cdline          do kan = 1,n_snd_ch !sounder channels for research

cdline        call taugim(model_t(1,i,j),model_mr(1,i,j),ozo,
cdline    1                  theta(i,j),ngoes,kan,tau)
cdline        call weight_func (tau,model_p,40,w_model,p_dm,39)
cdline        s_radiance(i,j,kan) = gimrad(tau,model_t(1,i,j),tskin(i,j),
cdline    1                              kan,lsfc(i,j),psfc(i,j),emiss)
cdline        s_btemp(i,j,kan) = britgo(s_radiance(i,j,kan),kan)

cdline          enddo ! Kan for sounder


      enddo
      enddo




c  generate table of clear sounder btemps, computed and observed

cdline       open (34,file='sounder.out',form='formatted')

cdline       do j = 1, jj
cdline       do i = 1, ii

cdline        do kan = 1,n_snd_ch
cdline        if(rads(i,j,kan).gt.0.0 .and. rads(i,j,kan).lt.500.) then
cdline        rads(i,j,kan) = britgo(rads(i,j,kan),kan)
cdline        endif
cdline        enddo ! kan


cdline          write(34,77) i,j,(s_btemp(i,j,kan), rads(i,j,kan) 
cdline    1                      ,kan=1,n_snd_ch)
cdline77        format (i2,1x,i2,1x,36(f7.3,1x))

cd          write(34,*) i,j,ch4(i,j),rads(i,j,8)

cdline       enddo
cdline       enddo

cdline       close(34)




c   attempt powell method correction of layer humidity in clear areas
c   only for starters.


      do j = 1,jj
         do i = 1,ii

            factor(i,j) = rmd
            factor2(i,j) = rmd



            if (ch3(i,j).eq.rmd) then
               print*, 'missing data in channel 3 abort', i,j

            elseif (ch4(i,j).eq.rmd) then
               print*, 'missing data in channel 4 abort', i,j

            elseif (ch5(i,j).eq.rmd) then
               print*, 'missing data in channel 5 abort', i,j

            else
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

c  this is a diagnostic output for plotting with psiplot (fort24.dat)

c                  write (24, *) ch4(i,j), btemp(i,j,2)
c                  write (23, *) ch3(i,j), btemp(i,j,1)
c                  write (25, *) ch5(i,j), btemp(i,j,3)



                  do k = 1,3
                     xi(k,k) = -.0001
                     x(k) = 1.0
                  enddo


                  if(isnd.eq.0) then ! USE AS IMAGER DATA
                     radiance_ob(1) = plango(ch3(i,j),23)
                     radiance_ob(2) = plango(ch4(i,j),24)
                     radiance_ob(3) = plango(ch5(i,j),25)
                  endif

                  if(isnd.eq.1) then ! USE AS SOUNDER DATA
                     radiance_ob(1) = plango(ch3(i,j),10)
                     radiance_ob(2) = plango(ch4(i,j),8)
                     radiance_ob(3) = plango(ch5(i,j),7)
                  endif


                  do k = 1,40
                     w_cost(k) = model_mr(k,i,j)
                     t_cost(k) = model_t (k,i,j)
                     ozo_cost(k) = ozo(k)
                     p_cost(k) = model_p(k)
                  enddo

                  tskin_cost = tskin(i,j)
                  theta_cost = theta(i,j)
                  lsfc_cost = lsfc(i,j)
                  psfc_cost = psfc(i,j)



                  if(cld(i,j).eq.0.) then
c don't match low atmosphere (use func, not func3)

                     call powell (x,xi,3,3,ftol,iter(i,j),fret,func)

c               else !clouds... don't match low atmosphere
c               call powell (x,xi,3,3,ftol,iter(i,j),fret,func)
                  endif


                  write(6,33) abs(x(1)), abs(x(2)),abs(x(3)),
     1                 i,j,psfc(i,j),iter(i,j)
 33               format(3(f7.2,2x),i3,i3,1x,f7.2,i3)



                  if (cld(i,j) .eq. 0. .and. iter(i,j) .lt. 50
     1                 .and. abs(abs(x(2))-1.) .lt. .05 ) then
                     factor(i,j) = abs(x(3))
                     factor2(i,j) = abs(x(2))
                  else
                     write(6,*) i,j, '  .... coordinate rejected', 
     1                    abs(x(2)),iter(i,j), cld(i,j)

                     failures = failures + 1

                  endif


                  write(6,*) blank

               endif            !end of powell function



            endif               ! end of missing data flag test

         enddo
      enddo

       write(6,*) failures,' failures occurred due to layer confusion' 
       write(6,*) '...non-convergence, or clouds'




c  modify original lq3 file with new factors for comparison tests.
c modify lq3 only in clear areas as defined by lc3.


c  analyze top level adjustments.

       pn = 0

       do j = 1,jj
          do i = 1,ii
             if (factor(i,j).ne.rmd ) then
                pn = pn+1
                points(1,pn) = factor(i,j)
                points(2,pn) = i
                points(3,pn) = j
                mask(i,j) = 1
                data_anal(i,j) = factor(i,j)
             endif
          enddo
       enddo

c  d-lines that follow are for  RESEARCH purposes!

cdline        commentline = 'Early version without satellite input'
cdline        path='/data/mdlg/birk/lq3_before/'
cdline         do k = 1,kk
cdline         do j = 1,jj
cdline         do i = 1,ii
cdline             if(sh(i,j,k) .le. 0.0 ) then
cdline                 research_output(i,j,k) = rmd
cdline             else
cdline                 research_output(i,j,k) = sh(i,j,k)
cdline             endif
cdline         enddo
cdline         enddo
cdline         enddo
cdline        call writef_sp
cdline    1          (i4time,commentline,mlevel,path,research_output,
cdline    1            ii,jj,kk,istatus)


cdline        call opngks
cdline        call plotfield(sh(1,1,17),ii,jj)

      if (pn.ne.0) then

         call prep_grid(ii,jj,data_anal,points,pn)
         call slv_laplc (data_anal,mask,ii,jj)
         call smooth_grid2 (ii,jj,data_anal,1)
cdline        call plotfield (data_anal,ii,jj)



      else
         write(6,*) 
     1        'pn = 0,no acceptable data to analyze for adjustment'
         return

      endif

c       modify lq3 field  top level

      do j = 1,jj
         do i = 1,ii
            do k = 14,21        !between 475 and 100 mb

               sh(i,j,k) = sh(i,j,k) * data_anal(i,j)


            enddo
         enddo
      enddo



cdline        call plotfield(sh(1,1,17),ii,jj)
cdline        call clsgks

cdline        dummy_string = 'mv gmeta /data/mdlg/birk/gmeta/'//filename

cdline        call system (dummy_string)
cdline        commentline = 'Final version with satellite input'
cdline        path='/data/mdlg/birk/lq3_after/'
cdline         do k = 1,kk
cdline         do j = 1,jj
cdline         do i = 1,ii
cdline             if(sh(i,j,k) .le. 0.0 ) then
cdline                 research_output(i,j,k) = rmd
cdline             else
cdline                 research_output(i,j,k) = sh(i,j,k)
cdline             endif
cdline         enddo
cdline         enddo
cdline         enddo
cdline        call writef_sp
cdline    1     (i4time,commentline,mlevel,path,research_output,
cdline    1       ii,jj,kk, istatus)



      return
      end

