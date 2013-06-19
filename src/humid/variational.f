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
c FORTRAN 90 CONSTRUCTS USED OVER F77 CODE

      subroutine variational (
     1     sh,                  ! specific humidity g/g
     1     sfc_data,            ! struct type lbsi surface data
     1     i4time,              ! i4time of run (seconds)
     1     p_3d,                ! pressure hpa (laps vert grid)
     1     cloud,               ! cloud array
     1     istatus_cloud,       ! clouds istatus
     1     sat,                 ! saturated specific hum.
     1     qadjust,             ! moisture needed for cloud formation
     1     t,                   ! lt1 (laps 3d temps)
     1     mdf,
     1     qs,kstart,
     1     ngoes,               ! goes satellite number
     1     isnd,                ! sounder switch
     1     sat_skip,            ! normally 1 for full resolution
     1     gw1,gw2,gw3,
     1     gww1,gww2,gww3,
     1     gvap_p,istatus_gvap,
     1     gps_data,
     1     gps_w,
     1     istatus_gps,
     1     q_snd,
     1     weight_snd,
     1     raob_switch,
     1     ii,jj,kk,             ! grid dimensions
     1     print_switch,
     1     covar_switch,path2covar)

 
c   By inclusion of the goes_sbn data into the laps moisture analysis, an
c   improvement in upper level moisture (above 500 mb) can be anticipated to be
c   about 70%.  Current research is pursuing using the satellite data in other
c   levels and other variables such as temperature. 
c
c   This routine interfaces GOES 8/10/11 satellite broadcast network data (and
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
c


      use module_sfc_structure
      implicit none
      include 'constants_optran.inc'
      include 'grid_fname.cmn'

c     include 'lapsparms.for'

c     parameter list variables


      integer ::  ii,jj,kk,print_switch,covar_switch
      type (lbsi), dimension (ii,jj) :: sfc_data
      real :: sh(ii,jj,kk)
      integer :: i4time
      real :: sat (ii,jj,kk)
      real ::  t(ii,jj,kk),p_3d(ii,jj,kk),mdf
      real :: cloud(ii,jj,kk)
      real :: qs(ii,jj)
      integer :: istatus_cloud
      real :: gw1(ii,jj),gww1(ii,jj)
      real :: gw2(ii,jj),gww2(ii,jj)
      real :: gw3(ii,jj),gww3(ii,jj)
      real :: gvap_p (ii,jj)
      integer :: istatus_gvap
      integer :: ngoes
      integer :: isnd
      integer :: sat_skip
      integer :: kstart(ii,jj)
      real :: gps_data(ii,jj)
      real :: gps_w (ii,jj)
      integer :: istatus_gps
      real :: q_snd(ii,jj,kk)
      real :: weight_snd(ii,jj,kk)
      integer :: raob_switch
      real, dimension(18)  :: testestrad,testbest
      character *256 cpath,tpath,path2covar



c internal variables

c      real covar (19,19,93,65),covar_temp(19,19,93,65)
      real, dimension (:,:,:,:), allocatable :: covar
      real, dimension (:,:,:,:), allocatable :: covar_temp
      integer l
      integer :: istatus
      integer :: i4time_sat
      integer :: i,j,k,k2,ijk
      real, dimension (40) :: local_model_p = (/.1,.2,.5,1.,1.5,2.,
     1     3.,4.,5.,7.,10.,15.,
     1     20.,25.,30.,
     1     50.,60.,70.,85.,100.,115.,135.,150.,200.,250.,300.,350.,400.,
     1     430.,475.,500.,570.,620.,670.,700.,780.,850.,920.,
     1     950.,1000./)
      real :: dummy
      real :: qadjust(ii,jj,kk)
      integer :: k500, k700
      integer :: sat_index


c climate model variables
      integer :: julian_day
      real :: standard_press(40),
     1     tempertur_guess(40),
     1     mixratio_guess(40)
      real :: rmd
      integer, parameter :: n_snd_ch = 22
      integer, dimension (7):: kanch =(/10,8,7,11,16,6,12/)
      integer :: restore_cost_rad_istatus


c dynamic dependent variables

      real :: ch3(ii,jj),ch4(ii,jj),ch5(ii,jj)
c      real :: mr(ii,jj,kk)
      real, dimension (:,:,:), allocatable:: mr
c      real :: t_l(kk,ii,jj), mr_l (kk,ii,jj),p_l(kk,ii,jj)
      real, dimension (:,:,:), allocatable :: t_l,mr_l,p_l

      real :: model_t(40,ii,jj), model_mr(40,ii,jj)

c     forward model variarles
      
c     new optran variables
      
      real :: tbest(n_snd_ch), radest(n_snd_ch)
      
c     old gimtau.f variables
      
      real ::  radiance(ii,jj,18),tskin(ii,jj),psfc(ii,jj),
     1     theta(ii,jj),
     1     ozo(40),gimrad,tau(40)
      real :: emiss
      integer :: kan,lsfc(ii,jj)
      real, dimension (40) :: model_p
      real :: t_fm(40),w_fm(40),ozo_fm(40)
      common/atmos/model_p,t_fm,w_fm,ozo_fm
      real :: btemp(ii,jj,18),plango
      real :: zenith               ! function call
      real :: pi, d2r
      real, external :: britgo
      
c     powell specific arrays
      real, dimension (3) :: x
      real, dimension (3,3) :: xi
      real :: ftol,fret
      integer, dimension (ii,jj) :: iter
      real, external :: func    ! function typing for cost function

      
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
      integer goes_good
      real cost_sec_za
      real cost_sfc_emis
      real cost_sfc_refl
      real cost_sec_solar
      
      real bias_correction      ! function

c     background covariance common block
      common /cost_background/ background_covar,cost_covar,covar_s
      real background_covar (3,93,65)
      real cost_covar (3)
      integer level7, level5, covar_s
      real covar_sum
      integer covar_count
      
c     optran common block
      
      common /cost_optran/btemp_ob, cost_kk, cost_p, cost_t_l,
     1     cost_mr_l, cost_tskin, cost_psfc, cost_julian_day, cost_lat,
     1     cost_theta, cost_isnd, cost_rad_istatus, cost_sec_za,
     1     cost_sfc_emis, cost_sfc_refl,cost_sec_solar
      
c     gvap common block
      
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
      
c     cloud common block
      
      common /cost_cloud/cost_cloud,cost_cld,cost_cloud_istatus,
     1     cost_sat,cost_qadjust
      integer cost_cloud_istatus
      real cost_cloud(500)
      real cost_cld,cost_sat(500)
      real :: cost_qadjust(500)
      
      integer goes_number
      
      common /sat_id/ goes_number

c     common block for gps
      common /cost_gps/cost_gps_data, cost_gps_weight,cost_gps_istatus
      integer cost_gps_istatus
      real cost_gps_data
      real cost_gps_weight

c     SND common block
      common /cost_snd/cost_snd_data, cost_snd_wt,cost_snd_istatus
      real cost_snd_data(500)
      real cost_snd_wt(500)
      integer cost_snd_istatus

c     Science common block
      common/cost_science/cost_comment_switch
      integer cost_comment_switch

c     Display common block
      common /cost_display/ display_btemps
      real, dimension(1000) :: display_btemps
      real, dimension(Nchan,ii,jj) :: field_display_btemps
      
c     analysis of the factor field
      integer pn
      real points(ii*jj,3)
      real data_anal(ii,jj)
      real ave,adev,sdev,var,skew,curt
      real upper_limit, lower_limit
      
c     cloud variables
      real cld(ii,jj)
      
c     moisture modified field
      real factor(ii,jj), factor2(ii,jj), factor3(ii,jj)
      
c     get latest filename
      character*256 path
      
c     laplace solver variables
      integer mask(ii,jj)
      
c     misc variables
      integer failures
      character*4 blank
      
      real rads (ii,jj,n_snd_ch)
      
      character*9 filename1,  filename
      character*9 grid_name

      integer len

c     code *******************************************************************

c     allocate temp arrays
      allocate (mr (ii,jj,kk))
      allocate (t_l (kk,ii,jj))
      allocate (mr_l (kk,ii,jj))
      allocate (p_l (kk,ii,jj))


c     read covariance data ----set for nam model ---------------------------------------

      if (covar_switch .ne. 0 ) then
c     allocate temp covar arrays
      allocate (covar(19,19,93,65),covar_temp(19,19,93,65))
         level7 = 7  !Okyeon needs to determine 700 level
         level5 = 11  !need to determine 500 level
         covar_s = covar_switch
         call s_len(path2covar,len)
         write (6,*) path2covar(1:len),len
         open(26, file=path2covar(1:len),form='unformatted',
     1        access='sequential',status='old')
         read (26) ((((covar(k,l,i,j),k = 1,19),l=1,19),i=1,93),j=1,65)
         close (26)
c         write (6,*)((((covar(k,l,i,j),k = 1,19),l=1,19),i=1,93),j=1,65)

C     INVERT COVAR since covar is direct from NAM and NOT laps.

c     use f90 construct to make arrays identical
         covar_temp = covar

         do k=1,19
            do l = 1,19
               do i = 1,93
                  do j = 1,65
                     covar (20-k,20-l,i,j) = covar_temp (k,l,i,j)
                  enddo
               enddo
            enddo
         enddo
                     

         if (covar_switch .eq. 1) then ! constant diagonal terms

            covar_sum = 0.0
            covar_count = 0

            do k = 1,19
               do i = 1,93
                  do j = 1,65
                     covar_sum = covar(k,k,i,j) +covar_sum
                     covar_count = covar_count +1
                  enddo
               enddo
            enddo
            covar_sum = covar_sum /float(covar_count) !covar_sum is avg

            do k = 1,3
               do i = 1, 93
                  do j = 1,65
                     background_covar(k,i,j) = covar_sum ! avg covar
                  enddo
               enddo
            enddo



         endif
         if (covar_switch.eq.2) then ! constant diagonal terms at each i,j



            do j = 1,65
               do i = 1,93
                  covar_sum = 0.0
                  covar_count = 0
                  do k = 1,19
                     covar_sum = covar(k,k,i,j) +covar_sum
                     covar_count = covar_count +1
                  enddo
                  covar_sum = covar_sum/float (covar_count)
                  do k = 1,19
                     background_covar(k,i,j) = covar_sum
                  enddo
               enddo
            enddo
           

         endif
         if (covar_switch.eq. 3) then! 3-layer diagnoal avg at each i,j

            do j = 1,65
               do i = 1,93
                  covar_sum = 0.0
                  covar_count = 0
                  do k = 1,level7
                     covar_sum = covar(k,k,i,j) +covar_sum
                     covar_count = covar_count +1
                  enddo
                  covar_sum = covar_sum/float (covar_count)
                  do k = 1,level7
                     background_covar(k,i,j) = covar_sum
                  enddo

                  covar_sum = 0.0
                  covar_count = 0
                  do k = level7+1,level5
                     covar_sum = covar(k,k,i,j) +covar_sum
                     covar_count = covar_count +1
                  enddo
                  covar_sum = covar_sum/float (covar_count)
                  do k = level7+1,level5
                     background_covar(k,i,j) = covar_sum
                  enddo

                  covar_sum = 0.0
                  covar_count = 0
                  do k = level5+1,19
                     covar_sum = covar(k,k,i,j) +covar_sum
                     covar_count = covar_count +1
                  enddo
                  covar_sum = covar_sum/float (covar_count)
                  do k = level5+1,19
                     background_covar(k,i,j) = covar_sum
                  enddo




               enddo
            enddo
           


         endif

c     deallocate covar temp arrays
      deallocate (covar)
      deallocate (covar_temp)

      endif

c     end read covariance data -----------------------------------------------
      
      write(6,*) 'new version 3/7/05 uses satellite gradients'
      write (6,*)

      cost_rad_istatus = 1      ! assume good
      goes_good = 1             ! assume good


c     code change to handle situation where direct radiance assimilation
c     is not desired.  Note here that ngoes comes into this routine 
c     assigned as zero (0).  This test takes this situation and assigns
c     cost_rad_istatus and goes_good to zero to insure that it not be used
c     however, it was discovered that this code still needs to compute 
c     synthetic radiances and to avoid a divide by zero, the FAKE assignment
c     of ngoes to goes 12 is made here.  This fake assignment should be 
c     goes 9 for asia if this is desired to run in that mode.  
      if (ngoes == 0) then
         cost_rad_istatus = 0
         goes_good = 0
         ngoes = 12 ! fake a satellite but don't use it in variational
c                     adjustment
      endif


      if (istatus_gps .eq. 1) write (6,*) 'GPS usage is a GO'

c     check sat_skip for zero, if zero skip routine
      
      if (sat_skip.le.0) then
         write (6,*) 'sat_skip parameter <=0, skipping POWELL entirely'
         return
      endif
      
c     assign sounder/imager parameter for powell method
      
      cost_isnd = isnd
      
c     assign goes number for common block to make avail where needed for OPTRAN
      
      goes_number = ngoes
      if (goes_number == 8) sat_index = 1 ! goes east
      if (goes_number == 10) sat_index = 2 ! goes west
      if (goes_number == 9) sat_index = 3 ! pacific goes
      if (goes_number == 12) sat_index = 1 ! new goes east
      if (goes_number == 11) sat_index = 2 ! new goes west
      
      
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
      
c     grid_fnam_common = 'nest7grid'
      call get_laps_config(grid_fnam_common,istatus)
      
c     grid_name = 'nest7grid'
c     call get_laps_config(grid_name,istatus)
      
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
         elseif (ngoes.eq.12) then
            path = path(1:len)//'goes12/'
         elseif (ngoes.eq.11) then
            path = path(1:len)//'goes11/'
            len = len + 7
         endif
         
         
         call get_latest_file (path,i4time,filename1,istatus)
         
         if (istatus.ne.1) then
            write(6,*) 'variational.f::
     1           Failure on call to get_latest_file'
            return
         endif
         
         write (6,*) 'Attempting: ', filename1
c     convert filename to i4time_sat
         call i4time_fname_lp (filename1,i4time_sat,istatus)
         write (6,*) 'Getting satellite radiances (lvd)'
         call read_lvd_3_4_5 (path,i4time_sat,ch3,ch4,ch5,
     1        ii,jj,kk,ngoes,istatus)
         
         if (istatus.ne.1) then
            write(6,*) 'error getting satellite data'
            write(6,*) 'aborting goes_sbn module'
            goes_good = 0
         endif
         
         do j = 1,jj
            do i = 1,ii
               if(ch3(i,j).le.0.0 .or. ch4(i,j).le.0.0 .or.
     1              ch5(i,j).le.0.0) then
                  istatus = 0
c                  write(6,*) 'Zero value radiance discovered'
c                  write(6,*) 'IMAGER DATA'
c                  write(6,*) '   aborting satellite moisture'
c                  write(6,*) '   data untrustworthy'
                  goes_good  = 0
               endif
            enddo
         enddo    
                
         
         write(6,*) ' '
         write(6,*) ' '
         write(6,*) 'Using LVD data from: ', filename1
         write(6,*) ' '
         write(6,*) ' '
         
      endif                     ! get IMAGER data only
      
c     acquire sounder data
      
      if(isnd.eq.1) then        ! get SOUNDER data only
         
         call rsr (i4time, rads, ii,jj,18,ngoes, istatus)
         if (istatus .ne. 1) then
            write (6,*) 'error obtaining sounder radiances'
            goes_good = 0
            write(6,*) 'goes_good is ',goes_good
         endif
         
      endif                     ! only get SOUNDER data
      
c     --------- at this point, the existance of satellite data is established,
c     ---------  should be cost effective to continue.
      
      
c     set up time for regular laps interval
c     generate filename from 14time for julian day extraction later
      
      call make_fnam_lp (i4time, filename, istatus)
      
c     fill tskin with laps sfc temp from structure

      tskin = sfc_data%sfc_temp

c     get laps surface pressure now replaced with value from structure sfc_data
      
      print*, 'Transferring surface pressure from structure'

      psfc = sfc_data%sfc_pres

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
      
      write (6,*) 'Running GOES',ngoes,' forward model OPTRAN90 vsn'
      write (6,*) 'Modify tskin and psfc for cloud top if present'
      
      do j = 1,jj
         do i = 1,ii
            do k = kk,1,-1
               
               if(cloud(i,j,k).ge.1.0) then ! assume cloud top
                  
                  if(p_3d(i,j,k).lt.psfc(i,j)) then ! above ground level
                     psfc(i,j) = p_3d(i,j,k)
                     tskin(i,j) = t(i,j,k)
                     cld(i,j) = cloud(i,j,k)
                     
                  else
                     print*, 'cloud below ground'
                  endif
                  go to 55
               endif
            enddo
            
 55         continue
            
         enddo
      enddo
      
c     modify sounding to convert sh to mr and model organization
c     assign 0.0 moisture where there is missing data.
      
      do i = 1,ii
         do j = 1,jj
            do k = 1,kk
               
               if(sh(i,j,k) .ne. -1.e30) then
                  call sh2mr (sh(i,j,k), mr(i,j,k),istatus )
                  if(istatus .eq.1) then
                     mr_l(k,i,j) = mr(i,j,k)
                  else
                     Write(6,*) 'trap sh2ml, i,j,k,mr ',i,j,k,mr(i,j,k)
                     write(6,*) 'assigning mr to zero'
                     mr_l(k,i,j) = 0.0
                  endif
               else
                  mr_l(k,i,j) = 0.0
               endif
               t_l (k,i,j) = t(i,j,k)
               p_l (k,i,j) = p_3d(i,j,k)
               
            enddo
         enddo
      enddo
      
      read (filename(3:5),22) julian_day
 22   format (i3)
      
c     prepare to use forward model functions
c     here use goes 8 for reference (goes 10 not avail)
      
      call pfcgim (8)
      
      if(isnd .eq.1 .and. goes_good .eq.1) then ! use sounder data for ch3, ch4, ch5
         do j = 1, jj
            do i = 1, ii
               if(rads(i,j,10).eq.rmd .or.
     1              rads(i,j,10).le. 0.0) then
                  ch3(i,j) = rmd
                  if (rads(i,j,10).le.0.0) then
                     write(6,*) 'Zero in ch10 ',rads(i,j,10), i,j
                  endif
               else
                  ch3(i,j) = bias_correction (britgo(rads(i,j,10),10),
     1                 ngoes, 1, 10)
               endif
               if(rads(i,j,8).eq.rmd .or.
     1              rads(i,j,8).le.0.0) then
                  ch4(i,j) = rmd
                  if (rads(i,j,8).le.0.0) then
                     write(6,*) 'Zero in ch8 ',rads(i,j,8), i,j
                  endif
               else
                  ch4(i,j) = bias_correction (britgo(rads(i,j,8),8),
     1                 ngoes, 1, 8)
               endif
               if(rads(i,j,7).eq.rmd .or.
     1              rads(i,j,7).le.0.0) then
                  ch5(i,j) = rmd
                  if (rads(i,j,7).le.0.0) then
                     write(6,*) 'Zero in ch7 ',rads(i,j,7), i,j
                  endif
               else
                  ch5(i,j) = bias_correction (britgo(rads(i,j,7),7),
     1                 ngoes, 1, 7)
               endif
            enddo
         enddo
      endif                     ! sounder used
      
c     do for each gridpoint

c      write(29,*) 'i j diff ',(' ob',i,' est',i,i=1,18)

      do j = 1,jj,sat_skip
         do i = 1,ii,sat_skip
            
c     retrieve zenith angle for model from sfc_data structure
            
            theta(i,j) = 1./ sfc_data(i,j)%secza(sat_index) ! (cos of zen)
            theta(i,j) = acos(theta(i,j))/d2r
            
            if(abs(theta(i,j)) .ge. 70.) then
               ch3(i,j) = rmd   ! designed to through out processing
c     at the location where there is no possiblity of running the forward
c     model
               go to 864        !skip ofm computation here
            endif
            
c     insert call for OPTRAN for initial comparison with gimtau.f
c     note that optran is configured to return both sounder and imager
c     channels used in this algorithm.
            
            call ofm ( kk, p_l(1,i,j), t_l(1,i,j), 
     1           mr_l(1,i,j), sfc_data(i,j)%sfc_temp, 
     1           psfc(i,j),
     1           julian_day, sfc_data(i,j)%lat, theta(i,j), tbest,
     1           radest,
     1           sfc_data(i,j)%secza(sat_index),
     1           sfc_data(i,j)%sfc_emiss(1),
     1           sfc_data(i,j)%sfc_refl(1),
     1           sfc_data(i,j)%secsola
     1           )
            
            if(isnd.eq.0) then  ! IMAGER computation
               
               do kan = 1,3
                  
                  btemp(i,j,kan) = tbest (kan+1) ! optran90 2,3,4 (chan 3,4,5)
                  
               enddo            !kan
               
            endif               ! end IMAGER computation
            
            if(isnd.eq.1) then  ! SOUNDER computation
               
               do kan = 1,7
                  
                  btemp(i,j,kan) = tbest(kanch(kan))
                  
               enddo            ! kan
               
            endif               ! end SOUNDER computation
 864        continue
c            if (rads(i,j,1) .ne. mdf .and.
c     1           abs(sfc_data(i,j)%sfc_temp - tbest(8)).lt. 5.0) then 
c               write (29, *)i,j,
c     1              (rads(i,j,8)-radest(8)),
c     1              (rads(i,j,kan),radest(kan),kan=1,18)

c            endif
         enddo                  ! j
      enddo                     ! i

      close (29)

 865  continue
      
c     Execute powell method correction of layer humidity in clear areas

c     at this point in the code, both the forward model and measured radiances
c     are in arrays.  They can be printed out at this point as optional files

c     compare rads with new variable of modeled radiances for compare.
c     channels

      
      failures = 0
      

      factor  = rmd
      factor2 = rmd
      factor3 = rmd
            

      restore_cost_rad_istatus = cost_rad_istatus
      
      do j = 1,jj,sat_skip
         do i = 1,ii,sat_skip
            
            cost_rad_istatus = restore_cost_rad_istatus
            
            if(goes_good .eq. 0) cost_rad_istatus = 0
            
            do k   = 1,3
               x(k) = 1.0
            enddo
            
            if (cost_rad_istatus.ne.0) then ! reduce work 
            
               if (ch3(i,j).eq.rmd) then
                  print*, 'missing data in channel 3 noted', i,j
                  cost_rad_istatus = 0
                  
               elseif (ch4(i,j).eq.rmd) then
                  print*, 'missing data in channel 4 noted', i,j
                  cost_rad_istatus = 0
                  
               elseif (ch5(i,j).eq.rmd) then
                  print*, 'missing data in channel 5 noted', i,j
                  cost_rad_istatus = 0
                  
               endif ! end test of channel good 
            endif
            
            if (isnd.eq.1 .and. cost_rad_istatus.ne.0) then
               do k = 4,7
                  if (rads(i,j,kanch(k)) .eq. rmd) then
                     print*, 'missing data in sounder channel ',
     1                    kanch(k),' index ',i,j
                     cost_rad_istatus = 0 ! fail using radiance
                  endif
                  if (rads(i,j,kanch(k)) .le. 0.0) then
                     write(6,*) 'negative radiance,r,chan,i,j'
                     write(6,*) rads(i,j,kanch(k)), kanch(k),i,j
                     cost_rad_istatus = 0 ! fail using radiance
                  endif
               enddo
            endif
            
            continue
            
c     new check for viable radiance location
            
            if(cost_rad_istatus.eq.1 .and. 
     1           abs(ch4(i,j)-btemp(i,j,2)).lt.5.) then
               if(print_switch .eq. 1) then
                  write(6,*) 'Radiance passed clear chan test',
     1                 abs(ch4(i,j)-btemp(i,j,2))
               endif
               continue         ! pass clear channel test
            elseif (cost_rad_istatus .eq. 1) then
               if (print_switch .eq.1) then
                  write (6,*)'Radiance failed clear chan test',
     1                 abs(ch4(i,j)-btemp(i,j,2)), cost_rad_istatus
               endif
               cost_rad_istatus = 0 ! failing clear channel test
            endif
            
c     radiance data quality known, continue with normal run
            
            if (cost_rad_istatus.eq.1) then
               if (print_switch .eq. 1) then
                  write(6,32) ' Observed=',ch3(i,j),' Modeled='
     1                 ,btemp(i,j,1),' Diff=',(ch3(i,j)-btemp(i,j,1))
 32               format(1x,a10,f8.3,a9,f8.3,a6,f8.3)
                  write(6,*) ch4(i,j),btemp(i,j,2)
                  write(6,*) ch5(i,j), btemp(i,j,3)
               endif
            else
               if(print_switch .eq. 1 )then
                  write(6,*)'Radiance data not included in variational'
               endif
            endif
            
c     initialize cost function vector for scaling output
            
            do k = 1,3
               do k2 = 1,3
                  xi(k,k2) = 0.0
                  if(k.eq.k2)  xi(k,k) = 1.0
               enddo
            enddo
            
c     copy imager data into radiance arrays
            
            if(isnd.eq.0 .and. cost_rad_istatus.eq.1) then ! USE 
                                !AS IMAGER DATA, btemps
               btemp_ob(1) = ch3(i,j)
               btemp_ob(2) = ch4(i,j)
               btemp_ob(3) = ch5(i,j)
            endif
            
c     copy sounder data into radiance arrays
            
            if(isnd.eq.1 .and. cost_rad_istatus.eq.1) then ! USE
                                ! AS SOUNDER DATA, btemps
               btemp_ob(1) = ch3(i,j)
               btemp_ob(2) = ch4(i,j)
               btemp_ob(3) = ch5(i,j)
               btemp_ob(4) = bias_correction(britgo(
     1              rads(i,j,kanch(4)),kanch(4)),ngoes,1,kanch(4))
               btemp_ob(5) = bias_correction(britgo(
     1              rads(i,j,kanch(5)),kanch(5)),ngoes,1,kanch(5))
               btemp_ob(6) = bias_correction(britgo(
     1              rads(i,j,kanch(6)),kanch(6)),ngoes,1,kanch(6))
               btemp_ob(7) = bias_correction(britgo(
     1              rads(i,j,kanch(7)),kanch(7)),ngoes,1,kanch(7))
               
c     check for bad data in btemp_ob
               
               do k = 1,7
                  if (btemp_ob(k) .le.0.0 ) then
                     write(6,*) 'bad btemp_ob', btemp_ob(k),
     1                    kanch(k), filename1,' aborting'
                     istatus = 0
                     return
                  endif
               enddo
               
            endif
            
c     fill powell common block with profile data for routine variational
c     this code executed for all types of data

c     OKYEON COVARIANCE TESTING

            if (covar_switch.ne.0)then
               do k = 1,3
                  cost_covar(k) = background_covar(k,i,j)
               enddo
            endif
        
            
c     fill cost function for background atmosphere
            
            do k = 1, kk
               cost_p(k) = p_l(k,i,j)
               cost_t_l(k) = t_l(k,i,j) 
               cost_mr_l(k) = mr_l (k,i,j)
               cost_p1d(k) = p_3d(i,j,k)
               if (i.ne.1 .and. j.ne.1 .and. i.ne.ii .and. j.ne.jj)then
                  cost_grad_x(k) = (sh(i+1,j,k)-sh(i-1,j,k))*0.5
                  cost_grad_y(k) = (sh(i,j+1,k)-sh(i,j-1,k))*0.5
               else
                  cost_grad_x(k) = 0.0
                  cost_grad_y(k) = 0.0
               endif
               cost_data(k) = sh(i,j,k)
            enddo

            cost_kk = kk
            cost_tskin = tskin(i,j)
            cost_psfc = psfc(i,j)
            cost_julian_day = julian_day
            cost_lat = sfc_data(i,j)%lat
            cost_theta = theta (i,j)
            cost_sec_za = sfc_data(i,j)%secza(sat_index)
            cost_sfc_emis = sfc_data(i,j)%sfc_emiss(1)
            cost_sfc_refl = sfc_data(i,j)%sfc_refl(1)
            cost_sec_solar = sfc_data(i,j)%secsola

c     cost function data for gvap (now include gradients in x and y instead of actual values.
            
            cost_gvap_istatus = istatus_gvap
            if (i.gt.1 .and. j.gt.1 .and. i.lt.ii .and. j.lt.jj)then             
               cost_w1_x = (gw1(i+1,j)-gw1(i-1,j))*0.5
               cost_w2_x = (gw2(i+1,j)-gw2(i-1,j))*0.5
               cost_w3_x = (gw3(i+1,j)-gw3(i-1,j))*0.5
               cost_w1_y = (gw1(i,j+1)-gw1(i,j-1))*0.5
               cost_w2_y = (gw2(i,j+1)-gw2(i,j-1))*0.5
               cost_w3_y = (gw3(i,j+1)-gw3(i,j-1))*0.5
               IF (gw1(i+1,j)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw2(i+1,j)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw3(i+1,j)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw1(i-1,j)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw2(i-1,j)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw3(i-1,j)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw1(i,j+1)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw2(i,j+1)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw3(i,j+1)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw1(i,j-1)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw2(i,j-1)==MDF) cost_gvap_istatus = 0 ! skip
               IF (gw3(i,j-1)==MDF) cost_gvap_istatus = 0 ! skip
            else
               cost_w1_x = 0.0
               cost_w2_x = 0.0
               cost_w3_x = 0.0
               cost_w1_y = 0.0
               cost_w2_y = 0.0
               cost_w3_y = 0.0
               cost_gvap_istatus = 0  ! skip this boundary part
            endif

            cost_weight = gww1(i,j)
            cost_gvap_p = gvap_p(i,j)
            cost_kstart = kstart (i,j)
            cost_qs = qs (i,j)
            cost_ps = sfc_data(i,j)%sfc_pres
            cost_mdf = mdf
            
c     cost function data for cloud analysis
            
            cost_cloud_istatus = istatus_cloud
            cost_cld = cld(i,j)
            do k = 1,kk
               cost_cloud(k) = cloud(i,j,k)
               cost_sat(k) = sat(i,j,k)
               cost_qadjust(k) = qadjust (i,j,k)
            enddo

c     cost function data for gps

            cost_gps_data = gps_data(i,j)
            cost_gps_weight = gps_w (i,j)
            cost_gps_istatus = istatus_gps

c     cost function data for SND, fill with mixing ratio

            do k = 1,kk

               if(q_snd(i,j,k).ne.rmd) then
                  call sh2mr(q_snd(i,j,k)*1.e-3, ! change to g/g for call
     1                 cost_snd_data(k),istatus)
                  if(istatus.ne.1) cost_snd_data(k) = rmd
               else
                  cost_snd_data(k) = rmd
               endif

               cost_snd_wt(k) = weight_snd(i,j,k)
            enddo
            cost_snd_istatus = raob_switch
            
c     executed variational search

            ftol = 0.05

            cost_comment_switch = 1 ! turn on for one print of gps&cloud

            call powell (x,xi,3,3,ftol,iter(i,j),fret,func,print_switch)
            
c     check output of variational search for fret of 0.0 that
c     indicates no convergence (fret is result of func)  we assume
c     that the func will never be non-zero in real search.

            field_display_btemps (1:7,i,j) = display_btemps(1:7)
            
            if (fret.eq.0.0) then ! assume that func set to no
                                !convergence
               if (print_switch .eq. 1) then
                  write(6,*) ' FRET = 0, assume no converge, assign 1.0'
               endif
             
               do ijk = 1,3
                  x(ijk) = 1.0
               enddo
            endif
            
c     write out solution details
            if (iter(i,j) .eq. 0) then
               if (print_switch .eq. 1) then
                  write(6,*) 'No iterations..',i,j,cost_gps_istatus, 
     1                 cost_gvap_istatus, cost_cloud_istatus, fret
               else
                  write(6,33) 'TEML ',abs(x(1)), abs(x(2)),abs(x(3)),
     1                 i,j,fret,iter(i,j)
 33               format(a5,3(e11.6,2x),i3,1x,i3,1x,e11.6,i3)
               endif
            endif
               
c     criterion to accept the result is based on variational 
c     performance, the cloud test is in the process of being dropped.
            
c     if (cld(i,j) .eq. 0. .and. iter(i,j) .lt. 5
            if ( iter(i,j) .lt. 5
c     1           .and. abs(abs(x(1))-1.) .lt. .1
     1           .and.
     1           iter(i,j) .gt. 0 
     1           .and.  abs(x(3)).ne.0.0
     1           .and.  abs(x(2)).ne.0.0
     1           .and.  abs(x(1)).ne.0.0) then
               factor(i,j)  = abs(x(3))
               factor2(i,j) = abs(x(2))
               factor3(i,j) = abs(x(1)) ! low level
c     elseif (cld(i,j).gt.0.)then
c     write(6,*) '  .... coordinate rejected, cloudy'
            else
               if (print_switch .eq. 1) then
                  write(6,*) i,j, '  .... coordinate rejected', 
c     1                 abs(x(1)),iter(i,j), cld(i,j)
     1                 abs(x(1)),iter(i,j)
               endif
               
               failures = failures + 1
               
            endif
            if (print_switch .eq. 1) then
               write(6,*) blank
            endif
         enddo
      enddo
      
      if (print_switch .eq. 1) then
         write(6,*) failures,' failures from variational processing' 
         write(6,*) '...non-convergence '
      endif

c     ************* Section on processing resulting scaling fields  ****
      
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
      write (6,*) 'acceptable range', lower_limit, upper_limit,
     1     (lower_limit+upper_limit)/2., 'top'
      
      do i = 1,pn
         if(points(i,1) .lt. upper_limit .and.
     1        points(i,1) .gt. lower_limit) then
            data_anal(int(points(i,2)),int(points(i,3))) = points (i,1)
            if (print_switch .eq. 1) then
               write(6,*) points(i,1), 'assigned'
            endif
         else
            if (print_switch .eq. 1) then
               write(6,*) points(i,1), 'rejected ******************'
            endif
            points(i,2) = 0     ! flag for bad point in prep grid
         endif
      enddo
      
      if (pn.ne.0) then
         
         call prep_grid(ii,jj,data_anal,ii*jj,points,pn,istatus)
         if(istatus.eq.1) then
            call slv_laplc (data_anal,mask,ii,jj)
            call smooth_grid2 (ii,jj,data_anal,1)
            call two_d_stats(ii,jj,data_anal,rmd)
            if (print_switch .eq. 1) then
               write(6,*) 'TEMP processed slv_lapc 1'
            endif
            
         else
            if (print_switch .eq. 1) then
               write(6,*) 'TEMP not enough data, skipping slv_lapc'
            endif
            data_anal = 1.0
         endif
         
      else
         if (print_switch .eq. 1) then
            write(6,*) 
     1      'TEMP pn = 0,no acceptable data to analyze for adjustment'
         endif
         data_anal = 1.0
         return
         
      endif

      do j = 1,jj
         do i = 1,ii
      
      
c     find k500 (level at or above 500)
            
            do k = 1, kk
               if (p_3d(i,j,k) .le. 500.)then
                  k500 = k
                  go to 475
               endif
            enddo
 475        continue
            
c     find k700 (level at or above 700)
            do k = 1, kk
               if (p_3d(i,j,k) .le. 700.)then
                  k700 = k
                  go to 476
               endif
            enddo
 476        continue
            
            
c     modify lq3 field  top level
            
            
            do k = k500+1, kk   !between 475 and 100 mb
               
               sh(i,j,k) = sh(i,j,k) * data_anal(i,j)
c     1              ((abs(data_anal(i,j))-1.) * 
c     1              (p_3d(i,j,k)/500.) +1.)
               
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
      write (6,*) 'acceptable range', lower_limit, upper_limit,
     1     (lower_limit+upper_limit)/2., 'mid'
      
      do i = 1,pn
         if(points(i,1) .lt. upper_limit .and.
     1        points(i,1) .gt. lower_limit) then
            data_anal(int(points(i,2)),int(points(i,3))) = points (i,1)
            if (print_switch .eq. 1) then
               write(6,*) points(i,1), 'assigned'
            endif
         else
            if (print_switch .eq. 1) then
               write(6,*) points(i,1), 'rejected ******************'
            endif
            points(i,2) = 0     ! flag for bad point in prep grid
         endif
      enddo
      
      if (pn.ne.0) then
         
         call prep_grid(ii,jj,data_anal,ii*jj,points,pn,istatus)
         if(istatus.eq.1) then
            call slv_laplc (data_anal,mask,ii,jj)
            call smooth_grid2 (ii,jj,data_anal,1)
            call two_d_stats(ii,jj,data_anal,rmd)
            if (print_switch .eq. 1) then
               write (6,*) 'TEMP processed slv_lapc 2'
            endif
            
         else
            if (print_switch .eq. 1) then
               write(6,*) 'TEMP not enough data, skipping slv_lapc'
            endif
            data_anal = 1.0
         endif
         
      else
         if (print_switch .eq. 1) then
            write(6,*) 
     1        'TEMP pn = 0,no acceptable data to analyze for adjustment'
         endif
         data_anal = 1.0
         return
         
      endif
      
c     modify lq3 field  second level
      
      do j = 1,jj
         do i = 1,ii
            do k = k700+1,k500    !between 700 and 500 mb
               
               sh(i,j,k) = sh(i,j,k) * data_anal(i,j)
               
            enddo
         enddo
      enddo      
      
c     modify third (lowest level)

      
      pn = 0
      
      do j = 1,jj
         do i = 1,ii
            mask(i,j) = 0
            data_anal(i,j) = 1.
            if (factor3(i,j).ne.rmd ) then
               pn = pn+1
               points(pn,1) = factor3(i,j)
               points(pn,2) = i
               points(pn,3) = j
               mask(i,j) = 1
               data_anal(i,j) = factor3(i,j)
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
      write (6,*) 'acceptable range', lower_limit, upper_limit,
     1     (lower_limit+upper_limit)/2., 'low'
      
      do i = 1,pn
         if(points(i,1) .lt. upper_limit .and.
     1        points(i,1) .gt. lower_limit) then
            data_anal(int(points(i,2)),int(points(i,3))) = points (i,1)
            if (print_switch .eq. 1) then
               write(6,*) points(i,1), 'assigned'
            endif
         else
            if (print_switch .eq. 1) then
               write(6,*) points(i,1), 'rejected ******************'
            endif
            points(i,2) = 0     ! flag for bad point in prep grid
         endif
      enddo
      
      if (pn.ne.0) then
         
         call prep_grid(ii,jj,data_anal,ii*jj,points,pn,istatus)
         if(istatus.eq.1) then
            call slv_laplc (data_anal,mask,ii,jj)
            call smooth_grid2 (ii,jj,data_anal,1)
            call two_d_stats(ii,jj,data_anal,rmd)
            write (6,*) 'TEMP processed slv_lapc 2'
            
         else
            write(6,*) 'TEMP not enough data, skipping slv_lapc'
            data_anal = 1.0
         endif
         
      else
         write(6,*) 
     1        'TEMP pn = 0,no acceptable data to analyze for adjustment'
         data_anal = 1.0        ! assign entire array 1.0
         return
         
      endif
      
c     modify lq3 field low  level
      
      do j = 1,jj
         do i = 1,ii
            do k = 1,k700    !between sfc and 700 mb
               if(sh(i,j,k) .ne. mdf)
     1              sh(i,j,k) = sh(i,j,k) * data_anal(i,j)
               
            enddo
         enddo
      enddo  

c     fill field_display_btemps with missing data flag where not zero

      where (field_display_btemps < 1.0e-7 )
         field_display_btemps = mdf
      endwhere

c     all done using optran90 need to deallocate the coefficient arrays

      call optran_deallocate (istatus)

c     deallocate temp arrays
      deallocate (mr)
      deallocate (t_l)
      deallocate (mr_l)
      deallocate (p_l)
    

      return
      end
  
