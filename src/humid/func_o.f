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
cdis
cdis
cdis
cdis
cdis
cdis
c FORTRAN90 constructs
      real function func(x)

c     This module now includes more than just satllite data.  It is the 
c     minimization area for sounder radiance, GVAP, GPS, and cloud.

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
      include 'constants_optran.inc'
      
c     parameter variables
      
      real x(3)
      
c     optran specific arrays for powell function calling
      
      real btemp_ob (Nchan)
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
      integer cost_rad_istatus
      real cost_sec_za          ! secant zenith angle for O90
      real cost_sfc_emis        ! surface emissivity for O90
      real cost_sfc_refl        ! surface reflectance for O90
      real cost_sec_solar       ! solar secant angle for O90


c     background covariance common block
      common /cost_background/ background_covar,cost_covar,covar_s
      real background_covar (3,93,65)
      real cost_covar(3)
      integer level7, level5, covar_s
      real covar_sum
      integer covar_count

c     common from lq3driver main program, direct insert from namelist
c     common block for cloud weighting .. direct insert into func_o.f
      common /func_o_insert/ cloud_weight
      real cloud_weight      

c     radiomenter direct insert of weight form driver and namelist
      common /radiometer/ radio_wt
      real radio_wt

c     optran common 

      common /cost_optran/btemp_ob, cost_kk, cost_p, cost_t_l,
     1     cost_mr_l, cost_tskin, cost_psfc, cost_julian_day, cost_lat,
     1     cost_theta, cost_isnd, cost_rad_istatus, cost_sec_za,
     1     cost_sfc_emis, cost_sfc_refl,cost_sec_solar

c     gvap common

      common /cost_gvap/cost_w1_x,cost_w2_x,cost_w3_x,
     1     cost_w1_y,cost_w2_y,cost_w3_y,
     1     cost_gvap_p,cost_weight,
     1     cost_gvap_istatus,cost_data,cost_kstart,cost_qs,
     1     cost_ps, cost_p1d, cost_mdf,cost_grad_x,cost_grad_y
      real cost_w1_x,cost_w2_x,cost_w3_x,cost_gvap_p,cost_weight,
     1     cost_w1_y,cost_w2_y,cost_w3_y
      integer cost_gvap_istatus
      real cost_data(500)
      real cost_grad_x(500)
      real cost_grad_y(500)
      integer cost_kstart
      real cost_qs
      real cost_ps
      real cost_p1d(500)
      real cost_mdf

c     cloud common

      common /cost_cloud/cost_cloud,cost_cld,cost_cloud_istatus,
     1     cost_sat,cost_qadjust
      integer cost_cloud_istatus
      real cost_cloud(500)
      real cost_cld
      real cost_sat(500)
      real :: cost_qadjust(500)
      real cloud_temp           ! temp used in subroutine call

c     gps common

      common /cost_gps/cost_gps_data, cost_gps_weight,cost_gps_istatus
      integer cost_gps_istatus
      real cost_gps_data
      real cost_gps_weight

c     SND common block
      common /cost_snd/cost_snd_data, cost_snd_wt,cost_snd_istatus
      real cost_snd_data(500)   !mixing ratio
      real cost_snd_wt(500)
      integer cost_snd_istatus

c     Science common block
      common/cost_science/cost_comment_switch
      integer cost_comment_switch

c     Radiance display common block
      common /cost_display/ display_btemps
      real, dimension (1000) :: display_btemps

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
      real sec_za
      real sfc_emis
      real sfc_refl
      real sec_solar
      real lpw1,lpw2,lpw3

c     local monitor variables

      real max_func_rad, min_func_rad
      real max_func_back, min_func_back
      real max_func_gvap, min_func_gvap
      real max_func_gvap1
      real max_func_gvap2
      real max_func_gvap3
      real max_func_cloud, min_func_cloud
      real max_func_gps, min_func_gps
      real max_func_snd, min_func_snd

c     lcal variables

      integer i,j,k
      integer kan(18)
      real w(Nlevel)
      real tbest(nchan)
      real radest(nchan)
      logical first_time,first_gvap
      save first_time, first_gvap
      data first_time /.true./
      data first_gvap /.true./
      integer lvl500, lvl700, lvl100
      save lvl500, lvl700, lvl100
      real p1,p2,p3             !pressure tops for gvap layers
      real cloud_thresh         ! percent cloud for action (0.6)
      real GT  ! cloud functions
      real ipw                  !integrated water for GPS minimization
      real cloud_integral
      
c     externals
      
      real plango

c     code

      cloud_thresh = 0.6
      tbest = 0.0 ! initialize tbest array to zero each call (no carryover)
      cloud_integral = 0.0

c     simulate current variational code... remove when new cloud
c     routine is implemented this just tests good cloud i/o read
c     it has nothing to do with clouds at this current location

c      if(cost_cloud_istatus.ne.1 .or. cost_cld .ne. 0.0) then ! cloud 
c         func = 0.0
c         return                 ! ignore function 
c      endif

      func = 0.0                ! Default is start at minima
      max_func_rad = 0.0
      max_func_snd = 0.0

c     define G parameter

      GT = 1.0

c     constrain x to positive
      do i =1,3
         x(i) = abs(x(i))
      enddo
      
      if (first_time) then
         first_time  = .false.

c     initialize min variables

         min_func_rad = 1.e38
         min_func_back = 1.e38
         min_func_gvap = 1.e38
         min_func_cloud = 1.e38
         min_func_gps = 1.e38
         min_func_snd = 1.e38

c     set up for sounder if needed instead of imager
         
         if (cost_isnd == 1) then ! sounder radiances used
            
            kan(1) = 10         ! 7.4 like ch3
            kan(2) = 8          ! 11.02 like ch4
            kan(3) = 7          ! 12.02 like ch5
            kan(4) = 11
            kan(5) = 16
            kan(6) = 6
            kan(7) = 12


            
         else                   ! imager radiances used
            
         endif
         
c     determine layer levels 
         
         do i = cost_kk,1, -1
            if (cost_p(i) < 100.) lvl100=i
            if (cost_p(i) < 700.) lvl700=i
            if (cost_p(i) <= 500.) lvl500=i
         enddo
         
      endif                     ! first_time
      
c     fill local variables from common block for parameter list in ofm
c     and now snd

c     modify mixing ratio per predifined pressure layers.
         
      do i = 1,cost_kk
         
         if(i < lvl700) then  ! sfc to 780
            mr_l(i) = abs(x(1)) * cost_mr_l(i)
         elseif (i >= lvl700 .and. i.le. lvl500) then ! 700 to 500
            mr_l(i) = abs(x(2)) * cost_mr_l(i)
         elseif (i > lvl500 .and. i.le. lvl100) then ! between 475 and 100
            
c     the corresponding change must also be made in variational.f where
c     this information is applied.
            
            mr_l(i) = abs(x(3))*cost_mr_l(i)
c     1              ((abs(x(3))-1.)*(p(i)/500.) + 1.) *
c     1              cost_mr_l(i)
         else
            mr_l(i) =  cost_mr_l(i)
         endif
         
      enddo

C----------------END PRELIMINARY COMPUTATIONS AND PREP -----------




      
C     SATELLITE RADIANCE SECTION
c     note that radiance and tbest are brightness temperatures (not
c     radiance as the code might lend you to believe

      if (cost_rad_istatus .eq. 1) then
         
         do i = 1, nlevel
            p(i) = cost_p(i)
            t_l(i) = cost_t_l(i)
c            mr_l(i) = cost_mr_l(i) ! insert to remove effect of radiance term
         enddo
         kk = cost_kk
         tskin = cost_tskin
         psfc = cost_psfc
         julian_day = cost_julian_day
         lat = cost_lat
         theta = cost_theta
         sec_za = cost_sec_za
         sfc_emis = cost_sfc_emis
         sfc_refl = cost_sfc_refl
         sec_solar = cost_sec_solar
         
         
c     perform forward model computation for radiance
c     here is the main call to optran in this whole module, the real
c     time-consuming part of the code.
         
         call ofm ( kk, p, t_l, 
     1        mr_l, tskin, psfc,
     1        julian_day, lat, theta, tbest, radest,
     1        sec_za,sfc_emis,
     1        sfc_refl,
     1        sec_solar) 
         
c     compute cost function

         if (cost_isnd == 1) then ! SOUNDER radiances used
            if(cost_isnd == 1 .and. cost_cld >= cloud_thresh) then ! report conflict
               GT = 0.25 ! reduce influence of this term by 3/4 due to 
c     conflict with cloud analysis
            endif
            
            do j = 1,7          ! btemp_ob(1-7) is sounder btemp
               func = func + ( btemp_ob(j) -
     1              tbest(kan(j)) )**2/2.
            enddo
            func = func /0.02 !brightnessT variance (=1.5^2)

         else                   ! IMAGER situation (only 3 channels)
            
            do j = 1,3          ! radiance _ob(1-3) is imager btemp
               func = func + ( btemp_ob(j) -
     1              tbest(j+7) )**2/2.
            enddo
            func = func /0.02 !brightnessT variance (=1.5^2)
            
         endif

         func = func * GT ! importance reduced by cloud influence

         max_func_rad = func
         if (max_func_rad .ne. 0.0) then
            min_func_rad = min(min_func_rad,max_func_rad)
         endif
            
         
c     stability cost is identical for both imager and sounder
         
      endif                     ! cost_rad_istatus


c     fill display_btemps section for later display
      display_btemps(1:7) = tbest (1:7)


c
c     END SATELLITE RADIANCE SECTION
















C     BACKGROUND SECTION

c     background weighting, in effect even if radiance data are not present.

      if(covar_s .eq.0) then !old method, skip nam covar
         max_func_back = 0.0

         do j = 1,3
            max_func_back =   ((x(j) - 1.)**2) + max_func_back
         enddo

         max_func_back = max_func_back/(0.005**2) ! background error small

         func = func + max_func_back

         if (max_func_back .ne. 0.0) then
            min_func_back = min(min_func_back,max_func_back)
         endif
      endif

c     OKYEON SECTION FOR 3 TESTS

      if (covar_s.ne.0) then ! divide term by proper covar

         max_func_back = 0.0

         do j = 1,3
            max_func_back=((x(j) - 1.)**2)/cost_covar(j)+max_func_back
         enddo

         max_func_back = max_func_back/(0.005**2) ! background error small

         func = func + max_func_back

         if (max_func_back .ne. 0.0) then
            min_func_back = min(min_func_back,max_func_back)
         endif

      endif

c      write(6,*) 'func 1, ',func
c     END BACKGROUND SECTION














      
c     GVAP SECTION -- UNITS mm

      if (cost_gvap_istatus ==  1) then

c     test for weight of measurement
         if(cost_weight == cost_mdf) then !skip this step
            continue            ! skip this iteration
         else                   ! process gvap
            if (first_gvap) then
               first_gvap = .false.
               write(6,*) 'TEMP GVAP accepted'
            endif
            
c     integrate q for gvap layers
c     determine sigma level pressure analogs

            if (cost_gvap_p == cost_mdf) then
               cost_gvap_p  = cost_ps ! substitute for missing gvap press
            endif
            
            call sigma_to_p (0.1, cost_gvap_p, 0.9, p1)
            call sigma_to_p (0.1, cost_gvap_p, 0.7, p2)
            call sigma_to_p (0.1, cost_gvap_p, 0.3, p3)

c     now broken into gradiants in x and y
c     first the x grad

            call int_layerpw(x,cost_grad_x,cost_kstart,
     1           cost_qs,cost_gvap_p,cost_p1d,p1,p2,p3,lpw1,lpw2,lpw3,
     1           cost_kk,cost_mdf)
            
            if (p1 <= 300.0) then
               write(6,*)'TEMM ', x, p1,p2,p3,lpw1,lpw2,lpw3,
     1              cost_w1_x,cost_w2_x,cost_w3_x
            endif
            
            if (lpw2 == cost_mdf) then
               i = i
            endif

            max_func_gvap1 = 0.0
            max_func_gvap2 = 0.0
            max_func_gvap3 = 0.0
            
            if (abs(lpw1) .lt.  1000.) then
               max_func_gvap1 =  
     1              (lpw1-cost_w1_x)**2*cost_weight
            endif
            if (abs(lpw2) .lt.  1000.) then
               max_func_gvap2 =
     1              (lpw2-cost_w2_x)**2*cost_weight
            endif
            if (abs(lpw3) .lt. 1000.) then
               max_func_gvap3 =  
     1              (lpw3-cost_w3_x)**2*cost_weight
            endif
          
            max_func_gvap = (max_func_gvap1+max_func_gvap2
     1           +max_func_gvap3)

c     now repeat for the y gradient

            call int_layerpw(x,cost_grad_y,cost_kstart,
     1           cost_qs,cost_gvap_p,cost_p1d,p1,p2,p3,lpw1,lpw2,lpw3,
     1           cost_kk,cost_mdf)
            
            if (p1 <= 300.0) then
               write(6,*)'TEMM ', x, p1,p2,p3,lpw1,lpw2,lpw3,
     1              cost_w1_y,cost_w2_y,cost_w3_y
            endif
            
            if (lpw2 == cost_mdf) then
               i = i
            endif

            if (abs(lpw1) .lt.1000.) then
               max_func_gvap1 =  
     1              (lpw1-cost_w1_y)**2*cost_weight
            endif
            if (abs(lpw2) .lt.1000.) then
               max_func_gvap2 =
     1              (lpw2-cost_w2_y)**2*cost_weight
            endif
            if (abs(lpw3) .lt.1000.) then
               max_func_gvap3 =  
     1              (lpw3-cost_w3_y)**2*cost_weight
            endif
          
            max_func_gvap = (max_func_gvap1+max_func_gvap2
     1           +max_func_gvap3) + max_func_gvap
            max_func_gvap = max_func_gvap * 10000. ! analytic immportance

c     max_func_gvap is in mm (above) now divide by mm error to make
c     dimensionless
            max_func_gvap = max_func_gvap / ((3.27)**2) ! from SFOV worst case 
            func = func + max_func_gvap
            if (max_func_gvap .ne. 0.0) then
               min_func_gvap = min(min_func_gvap,max_func_gvap)
            endif
            
            
            
         endif                  !weight function test
      endif                     !data present test
C     END GVAP SECTION





















      
c    CLOUD SECTION -- UNITS (none, just a fraction 0->1)
      
      if (cost_cloud_istatus == 1) then ! cloud data present
         max_func_cloud = 0.0
         do k = 1,cost_kk
            cloud_integral = cloud_integral + cost_cloud(k)
            if (cost_data(k) /= cost_mdf .and. cost_data(k) > 0.0) then
               if(cost_cloud(k) >= cloud_thresh) then
                  cloud_temp = cost_data(k)
                  call cloud_sat (cost_cloud(k),cost_qadjust(k),
     1                 cloud_temp)! cloud temp is "cloud forming q"
                  if(k < lvl700 ) then ! sfc to 700
                     max_func_cloud = max_func_cloud + 
     1                    (cost_data(k)*x(1) - cloud_temp)**2 
                  elseif (k < lvl500) then ! 700-500
                     max_func_cloud = max_func_cloud + 
     1                    (cost_data(k)*x(2) - cloud_temp)**2
                  elseif (k < lvl100) then ! 500-100
                     max_func_cloud = max_func_cloud + 
     1                    (cost_data(k)*x(3) - cloud_temp)**2  
                  endif         ! layer test (1-3)
               endif            ! if cloudy check (cloudy enough?)
            endif               ! mdf check and bad value check
         enddo                  ! enddo k level

c here the weight of the cloud function is hard coded as 0.5  this will 
c now be modified to become a namelist parameter to help in the improve
c ment of the cloud analysis for Steve Albers needs. 3/8/12   DB
         max_func_cloud = max_func_cloud * cloud_weight
         func = func + max_func_cloud

c the following IF is suspect as a coding bug, it is deemed to have
c little or no significance on the processing
         if (max_func_cloud .ne. 0.0) then
            min_func_cloud = min(min_func_cloud,max_func_cloud)
         endif


      endif                     ! cloud data present
      
C     END CLOUD SECTION














      
c     GPS SECTION  !  UNITS cm
      
      if (cost_gps_istatus == 1) then

         call int_ipw (x,cost_p1d,cost_data,cost_kstart,
     1        cost_qs,cost_ps,ipw,cost_mdf,cost_kk)

         max_func_gps = 0.0
         
         max_func_gps = (cost_gps_data-ipw)**2*cost_gps_weight

         max_func_gps = max_func_gps/(0.03**2)  ! variance in cm**2

         func = func + max_func_gps
         if (max_func_gps .ne. 0.0) then
            min_func_gps = min(min_func_gps,max_func_gps)
         endif

      else
         continue
      endif
C     END GPS SECTION



      









c     RAOB SECTION (SND)

c     Error term computation revised 4/30/04 db 
c     now uses a 5% approximation in RH converted to approximate 
c     SH units

      if (cost_snd_istatus == 1) then
         
         do i = 1, cost_kk      ! all laps levels
            
            if (mr_l(i).ne.0.0.and.cost_snd_data(i).ne.cost_mdf) then
               max_func_snd = max_func_snd +
     1              (mr_l(i)-cost_snd_data(i))**2 * 
     1              cost_snd_wt(i) !weighted squared difference
     1              * radio_wt
     1              / (0.05*mr_l(i))**2
            endif
            
         enddo
         
      endif
      
      func = func + max_func_snd 
      if (max_func_snd .ne. 0.0) then
         min_func_snd = min(min_func_snd,max_func_snd)
      endif

C     END RAOB SECTION

      if (func .ne. func ) then ! Nan
         write (6,*) 'func is a nan in func_o.f'
         stop
      endif



















C     BOOKEEPING/MONITOR SECTION

c     print test output

c      write (24,*) x,
c     1     max_func_back/func,
c     1     max_func_rad/func,
c     1     max_func_gvap1/func,
c     1     max_func_gvap2/func,
c     1     max_func_gvap3/func,
c     1     max_func_cloud/func,
c     1     max_func_gps/func,func
C     END BOOKEEPING/MONITOR SECTION


      
      return
      end
