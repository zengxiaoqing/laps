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

       subroutine laps_anl(uobs,vobs
     1     ,obs_point,max_obs,nobs_point                               ! I
     1     ,n_radars,istat_radar_vel                                   ! I
     1     ,vr_obs_unfltrd,vr_nyq,v_nyquist_in,idx_radar_a             ! I
!    1     ,upass1,vpass1                                              ! O
     1     ,n_var                                                      ! I
     1     ,uanl,vanl                                                  ! O
     1     ,wt_p,weight_bkg_const,rms_thresh_wind                      ! I/L
     1     ,max_radars                                                 ! I
     1     ,n_radarobs_tot_unfltrd,rlat_radar,rlon_radar,rheight_radar ! I
     1     ,thresh_2_radarobs_lvl_unfltrd                              ! I
     1     ,thresh_4_radarobs_lvl_unfltrd                              ! I
     1     ,u_laps_bkg,v_laps_bkg                                      ! I/L
     1     ,imax,jmax,kmax,lat,lon                                     ! I
     1     ,i4time,grid_spacing_m                                      ! I
     1     ,r_missing_data                                             ! I
!    1     ,i_3d                                                       ! I
     1     ,l_derived_output,l_grid_north,l_3pass,l_correct_unfolding  ! I
!    1     ,n_iter_wind_in
     1     ,weight_cdw,weight_sfc,weight_pirep,weight_prof,weight_radar! I
     1     ,istatus)

!     This routine uses the inputted wind data and actually does the analysis

!     1992          Steve Albers
!     1992 Apr      Steve Albers      Subroutine to make u/v radar obs
!     1992 Apr      Steve Albers      Looping for multiple radars
!     1992 Aug      LincLb/S. Albers  More efficient, combining u+v analyses
!     1992 Aug      LincLb/S. Albers  Barnes_multivariate changed to be more
!                                     parallelizable
!     1992 Oct      S. Albers         Add intermediate pass with multi-Doppler
!                                     obs only
!     1994 May      S. Albers         Add Local wt_p_spread array so that
!                                     wt_p array is not modified by analysis
!     1994 Nov 01   S. Albers         Remove lapsprms.for include, add arguments
!     1995 Aug      S. Albers         Some U/V arrays now equivalenced
!     1995 Sep      K. Brewster       Added nyquist as input grid array

!********************ARGUMENT LIST********************************************

      integer*4 max_obs
!     parameter (max_obs = 40000)       
      include 'barnesob.inc'
      type (barnesob) obs_point(max_obs)      ! Full Wind Obs  - Non-radar data
      type (barnesob) obs_point_qced(max_obs) ! QC'd Obs       - Non-radar data
      type (barnesob) obs_radar(max_obs)      ! Full Wind Obs  - Radar data
      type (barnesob) obs_barnes(max_obs)     ! Full Wind Obs  - All Data

      integer n_var                                                ! Input
      integer*4 imax,jmax,kmax        ! 3D array dimensions        ! Input

!     3D arrays of u/v observations, all data sources, one datum per gridpoint.
      real*4 uobs(imax,jmax,kmax),vobs(imax,jmax,kmax)             ! Input

      integer*4 n_radars   ! Actual number of radars having data     Input

!     First pass analyzed winds
!     real*4 upass1(imax,jmax,kmax), vpass1(imax,jmax,kmax)        ! Output

!     Final pass analyzed winds
      real*4 uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)             ! Output
      real*4 varbuff(imax,jmax,kmax,n_var)                         ! Equiv Abv

!     3D array of observation weights, depends on data type
!     The choices are outlined below
      real*4    wt_p(imax,jmax,kmax)                               ! Input
      real*4    wt_p_radar(imax,jmax,kmax)                         ! Local

!     Model background field
      real*4 u_laps_bkg(imax,jmax,kmax),v_laps_bkg(imax,jmax,kmax) ! Input

!     Arrays of lat and lon for each gridpoint
      real*4 lat(imax,jmax),lon(imax,jmax)                         ! Input

      integer*4 i4time                                             ! Input

      real*4 grid_spacing_m                                        ! Input
      real*4 r_missing_data                  ! missing data value    Input

!-----Radar Data -----------------------------------------------------------

      integer*4 max_radars            ! max possible # of radars         Input

!     4D Radial velocity array (all radars)
      real*4 vr_obs_unfltrd(imax,jmax,kmax,max_radars)                 ! Input
      real*4 vr_nyq(imax,jmax,kmax,max_radars)                         ! Input
      integer*4 idx_radar_a(max_radars)                                ! Input

!     Nyquist velocity (if known and constant) for each radar
      real*4 v_nyquist_in(max_radars)                                  ! Input

!     Location of each radar
      real*4 rlat_radar(max_radars),rlon_radar(max_radars)             ! Input
     1                     ,rheight_radar(max_radars)

      integer*4 thresh_2_radarobs_lvl_unfltrd                          ! Input
     1         ,thresh_4_radarobs_lvl_unfltrd

!     # of radar obs before filtering for each radar (modified by QC)
      integer*4 n_radarobs_tot_unfltrd(max_radars)                     ! Input/Modified

!--------------------------------------------------------------------------------

!     real*4 uobs_diff(imax,jmax,kmax),vobs_diff(imax,jmax,kmax)       ! Local
      real*4, allocatable, dimension(:,:,:) :: uobs_diff               ! Local
      real*4, allocatable, dimension(:,:,:) :: vobs_diff               ! Local
      real*4, allocatable, dimension(:,:,:) :: upass1                  ! Local
      real*4, allocatable, dimension(:,:,:) :: vpass1                  ! Local
      real*4, allocatable, dimension(:,:,:) :: pres_3d                 ! Local

      real*4 varobs_diff_spread(imax,jmax,kmax,n_var)                  ! Local

      integer*4 n_obs_lvl(kmax)                                        ! Local
      logical  l_analyze(kmax) ! This depends on presence of radar obs ! Local
      logical  l_derived_output ! Flag for producing derived output    ! Input
      logical  l_grid_north     ! Flag for grid north or true north    ! Input
      logical  l_3pass          ! Flag for doing 3 pass analysis       ! Input
      logical  l_correct_unfolding ! Flag for dealiasing               ! Input
      logical  l_point_struct

      real*4   rms_thresh                                              ! Input
      real*4   weight_bkg_const                                        ! Input

!     These are the weights of the various data types (filling the 3D array)
      real*4 weight_cdw,weight_sfc,weight_pirep,weight_prof
     1      ,weight_radar                                              ! Input

      integer*4 istatus         ! (1 is good)                          ! Output

      integer*4  n_fnorm_dum

      character*3 c3_string

!****************END DECLARATIONS *********************************************

      write(6,*)' Subroutine laps_anl...'

      l_point_struct = .true. 

csms$serial(default=ignore)  begin              

!     Compare background to obs
      call compare_wind(
     1            u_laps_bkg,v_laps_bkg,' FG ',
     1            istat_radar_vel,max_radars,vr_obs_unfltrd,n_radars,
     1            rlat_radar,rlon_radar,rheight_radar,
     1            lat,lon,
     1            imax,jmax,kmax,r_missing_data,
     1            obs_point,max_obs,nobs_point,l_point_struct,
     1            weight_pirep,weight_prof,weight_sfc,weight_cdw,
     1            uobs,vobs,wt_p,istatus)
csms$serial end

      rms_thresh_norm = rms_thresh_wind          

      n_iter_wind = 1

      do iter = 1,n_iter_wind

csms$serial(<wt_p_radar, varobs_diff_spread, 
csms$>       rms_thresh , out>:default=ignore)  begin
      if(.true.)then ! Experimental
          if(iter .ge. 2)then
              write(6,*)
     1        ' Replacing background with the analyzed winds, iter=',ite
     1r
              do i = 1,imax
              do j = 1,jmax
              do k = 1,kmax
                  u_laps_bkg(i,j,k) = uanl(i,j,k)
                  v_laps_bkg(i,j,k) = vanl(i,j,k)
              enddo ! k
              enddo ! j
              enddo ! i
          endif
      endif

!     Subtract the background from the non-radar obs, then apply QC thresholds
!     and spread the obs vertically.

      write(6,91)
91    format(1x,' Subtracting the background from the obs'
     1         ,' then spreading the obs vertically.'
     1         /'       i    j    k  kk   udf   vdf     '
     1         ,'uob   vob     ubg   vbg vcdf  wt')

      n_qc_acrft_bad = 0
      n_qc_cdw_bad = 0
      n_qc_sfc_bad = 0
      n_qc_prof_bad = 0
      n_qc_total_bad = 0

      n_qc_acrft_good = 0
      n_qc_cdw_good = 0
      n_qc_sfc_good = 0
      n_qc_prof_good = 0
      n_qc_total_good = 0

      iwrite = 1

      qc_thresh = 30. ! Threshold speed for throwing out the ob

      allocate( pres_3d(imax,jmax,kmax), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate pres_3d'
          stop
      endif

      call get_pres_3d(i4time,imax,jmax,kmax,pres_3d,istatus)
      if(istatus .ne. 1)return

      call get_rep_pres_intvl(pres_3d,imax,jmax,kmax,rep_pres_intvl
     1                       ,istatus)

      deallocate(pres_3d)

      if(.not. l_point_struct)then

       allocate( uobs_diff(imax,jmax,kmax), STAT=istat_alloc )
       if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate uobs_diff'
          stop
       endif

       allocate( vobs_diff(imax,jmax,kmax), STAT=istat_alloc )
       if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate vobs_diff'
          stop
       endif

       do j=1,jmax
       do i=1,imax
        do k = 1,kmax
          if(wt_p(i,j,k) .ne. r_missing_data)then
            if(uobs(i,j,k) .ne. r_missing_data      .and.
     1         vobs(i,j,k) .ne. r_missing_data              )then

              speed_bkg  = sqrt(u_laps_bkg(i,j,k)**2
     1                        + v_laps_bkg(i,j,k)**2)

              uobs_diff(i,j,k) = uobs(i,j,k) - u_laps_bkg(i,j,k)
              vobs_diff(i,j,k) = vobs(i,j,k) - v_laps_bkg(i,j,k)
              speed_diff = sqrt(uobs_diff(i,j,k)**2
     1                        + vobs_diff(i,j,k)**2)

!             Apply QC check to the OB against the background analysis
              if(
!                Make sure we actually have a real reference background
     1           (speed_bkg .gt. 0.) .and.

!                General QC check
     1           (speed_diff .gt. qc_thresh

!              Stricter QC check for pireps
     1                       .OR. 
     1         (speed_diff .gt. 10. .and. wt_p(i,j,k) .eq. weight_pirep)

!              Stricter QC check for Cloud Drift Winds
     1                       .OR. 
     1         (speed_diff .gt. 10. .and. wt_p(i,j,k) .eq. weight_cdw)

!              Stricter QC check for profilers
     1                       .OR. 
     1         (speed_diff .gt. 22. .and. wt_p(i,j,k) .eq. weight_prof)
     1                                                                 )

     1                                                          )then

                ! Throw out the ob
                  if(wt_p(i,j,k) .eq. weight_pirep)then
                      n_qc_acrft_bad = n_qc_acrft_bad + 1
                      write(6,101,err=199)i,j,k
     1                  ,uobs_diff(i,j,k)
     1                  ,vobs_diff(i,j,k)
     1                  ,uobs(i,j,k)
     1                  ,vobs(i,j,k)
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,wt_p(i,j,k)
101                   format(' Prp QCed out - ',2i5,i4,1x,3(2x,2f5.0)
     1                                         ,f5.0,f5.2)

                  elseif(wt_p(i,j,k) .eq. weight_cdw)then
                      n_qc_cdw_bad = n_qc_cdw_bad + 1
                      write(6,111,err=199)i,j,k
     1                  ,uobs_diff(i,j,k)
     1                  ,vobs_diff(i,j,k)
     1                  ,uobs(i,j,k)
     1                  ,vobs(i,j,k)
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,wt_p(i,j,k)
111                  format(' Cdw QCed out - ',2i5,i4,1x,3(2x,2f5.0)
     1                                        ,f5.0,f5.2)

                  elseif(wt_p(i,j,k) .eq. weight_sfc)then
                      n_qc_sfc_bad = n_qc_sfc_bad + 1
                      write(6,121,err=199)i,j,k
     1                  ,uobs_diff(i,j,k)
     1                  ,vobs_diff(i,j,k)
     1                  ,uobs(i,j,k)
     1                  ,vobs(i,j,k)
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,wt_p(i,j,k)
121                  format(' Sfc QCed out - ',2i5,i4,1x,3(2x,2f5.0)
     1                                        ,f5.0,f5.2)

                  elseif(wt_p(i,j,k) .eq. weight_prof)then
                      n_qc_prof_bad = n_qc_prof_bad + 1
                      write(6,131,err=199)i,j,k
     1                  ,uobs_diff(i,j,k)
     1                  ,vobs_diff(i,j,k)
     1                  ,uobs(i,j,k)
     1                  ,vobs(i,j,k)
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,wt_p(i,j,k)
131                  format(' Prf QCed out - ',2i5,i4,1x,3(2x,2f5.0)
     1                                        ,f5.0,f5.2)

                  endif ! Type of OB to write out

!                 Set the difference OB to missing (original ob left alone)
199               uobs_diff(i,j,k) = r_missing_data
                  vobs_diff(i,j,k) = r_missing_data
                  wt_p(i,j,k) = r_missing_data

                  n_qc_total_bad = n_qc_total_bad + 1

              else ! write out the good OB
                  if(wt_p(i,j,k) .eq. weight_pirep)then
                      n_qc_acrft_good = n_qc_acrft_good + 1
                      c3_string = 'Prp'
                  endif

                  if(wt_p(i,j,k) .eq. weight_cdw)then
                      n_qc_cdw_good  = n_qc_cdw_good + 1
                      c3_string = 'Cdw'
                  endif

                  if(wt_p(i,j,k) .eq. weight_sfc)then
                      n_qc_sfc_good   = n_qc_sfc_good + 1
                      c3_string = 'Sfc'
                  endif

                  if(wt_p(i,j,k) .eq. weight_prof)then
                      n_qc_prof_good  = n_qc_prof_good + 1
                      c3_string = 'Prf'
                  endif

                  n_qc_total_good = n_qc_total_good + 1

                  if(n_qc_total_good .le. 500 .OR. j .eq. (j/10)*10)then
                      iwrite = 1
                  else
                      iwrite = 0
                  endif

                  if(iwrite .eq. 1)then
                      write(6,201,err=202)c3_string,i,j,k
     1                  ,uobs_diff(i,j,k)
     1                  ,vobs_diff(i,j,k)
     1                  ,uobs(i,j,k)
     1                  ,vobs(i,j,k)
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,wt_p(i,j,k)
201                   format(1x,a3,2i5,i4,4x,f6.1,f6.1,2(2x,2f6.1)
     1                         ,f5.1,f5.2)
202                   continue
                  endif

              endif

              if(speed_bkg .gt. 200.)then
                  write(6,*)
     1                ' Error: first guess winds > 200. m/s detected'
                  istatus = 0
                  return
              endif

            else
              write(6,*)' ERROR in laps_anl: MISSING DATA'
              istatus = 0
              return

            endif ! U and V .ne. MISSING


          else  ! wt_p .eq. MISSING; set these to missing just in case
            uobs_diff(i,j,k) = r_missing_data
            vobs_diff(i,j,k) = r_missing_data

          endif ! wt_p .ne. MISSING

        enddo ! k

       enddo ! i
       enddo ! j

       varobs_diff_spread(:,:,:,1) = uobs_diff
       varobs_diff_spread(:,:,:,2) = vobs_diff

       deallocate(uobs_diff)
       deallocate(vobs_diff)

      else ! Perform QC of data structure info (l_point_struct = .true)
          do i_ob = 1,nobs_point
              i = obs_point(i_ob)%i                       
              j = obs_point(i_ob)%j                       
              k = obs_point(i_ob)%k                       

              u = obs_point(i_ob)%valuef(1)
              v = obs_point(i_ob)%valuef(2)

              speed_bkg  = sqrt(u_laps_bkg(i,j,k)**2
     1                        + v_laps_bkg(i,j,k)**2)

              u_diff = u - u_laps_bkg(i,j,k)
              v_diff = v - v_laps_bkg(i,j,k)
              speed_diff = sqrt(u_diff**2 + v_diff**2)

!             speed_thresh = max(10.,0.2 * speed_bkg)

!             Apply QC check to the OB against the background analysis
              if(
!                Make sure we actually have a real reference background
     1           (speed_bkg .gt. 0.) .and.

!                General QC check
     1           (speed_diff .gt. qc_thresh

!              Stricter QC check for pireps
     1                       .OR. 
     1         (speed_diff .gt. 10. .and. 
     1                              obs_point(i_ob)%file .eq. 'pin')        

!              Stricter QC check for Cloud Drift Winds
     1                       .OR. 
     1         (speed_diff .gt. 10. .and. 
     1                              obs_point(i_ob)%file .eq. 'cdw')        

!              Stricter QC check for profilers
     1                       .OR. 
     1         (speed_diff .gt. 22. .and. 
     1                              obs_point(i_ob)%file .eq. 'pro')        
     1                                                                 )

     1                                                          )then

!                 Throw out the ob
                  if(obs_point(i_ob)%file .eq. 'pin')then
                      n_qc_acrft_bad = n_qc_acrft_bad + 1
                  elseif(obs_point(i_ob)%file .eq. 'cdw')then
                      n_qc_cdw_bad = n_qc_cdw_bad + 1
                  elseif(obs_point(i_ob)%file .eq. 'lso')then
                      n_qc_sfc_bad = n_qc_sfc_bad + 1
                  elseif(obs_point(i_ob)%file .eq. 'pro')then
                      n_qc_prof_bad = n_qc_prof_bad + 1
                  endif

                  write(6,231,err=232)obs_point(i_ob)%type(1:5)
     1                  ,i,j,k
     1                  ,u_diff
     1                  ,v_diff
     1                  ,u
     1                  ,v
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,obs_point(i_ob)%weight
231               format(a5,' QCed out - ',2i5,i4,1x,3(2x,2f5.0)
     1                                        ,f5.0,f5.2)
232               continue

                  n_qc_total_bad = n_qc_total_bad + 1

              else ! keep and write out the good OB
                  if(obs_point(i_ob)%file .eq. 'pin')then
                      n_qc_acrft_good = n_qc_acrft_good + 1
                      c3_string = 'Prp'
                  endif

                  if(obs_point(i_ob)%file .eq. 'cdw')then
                      n_qc_cdw_good  = n_qc_cdw_good + 1
                      c3_string = 'Cdw'
                  endif

                  if(obs_point(i_ob)%file .eq. 'lso')then
                      n_qc_sfc_good   = n_qc_sfc_good + 1
                      c3_string = 'Sfc'
                  endif

                  if(obs_point(i_ob)%file .eq. 'pro')then
                      n_qc_prof_good  = n_qc_prof_good + 1
                      c3_string = 'Prf'
                  endif

                  n_qc_total_good = n_qc_total_good + 1

!                 Assign data structure element (using difference ob)
                  obs_point_qced(n_qc_total_good) = obs_point(i_ob)
                  obs_point_qced(n_qc_total_good)%value(1) = u_diff
                  obs_point_qced(n_qc_total_good)%value(2) = v_diff

                  if(n_qc_total_good .le. 500 .OR. 
     1               n_qc_total_good .eq. (n_qc_total_good/10)*10)then
                      iwrite = 1
                  else
                      iwrite = 0
                  endif

                  if(iwrite .eq. 1)then
                      write(6,201,err=302)c3_string,i,j,k
     1                  ,u_diff
     1                  ,v_diff
     1                  ,u
     1                  ,v
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,obs_point_qced(n_qc_total_good)%weight
302                   continue
                  endif

              endif ! passed the QC test

          enddo ! i_ob

          ncnt_total = n_qc_total_good

      endif ! .not. l_point_struct

      write(6,*)
      write(6,*)' QC info for non-radar data (after remapping to grid)'
      write(6,601)n_qc_acrft_good,n_qc_acrft_bad
     1           ,pct_rejected(n_qc_acrft_good,n_qc_acrft_bad)
 601  format(' # of Aircraft   GOOD/BAD QC = ',2i6,7x
     1      ,'% rejected = ',f6.1)

      write(6,602)n_qc_cdw_good,n_qc_cdw_bad
     1           ,pct_rejected(n_qc_cdw_good,n_qc_cdw_bad)
 602  format(' # of CDWs       GOOD/BAD QC = ',2i6,7x
     1      ,'% rejected = ',f6.1)

      write(6,603)n_qc_sfc_good,n_qc_sfc_bad
     1           ,pct_rejected(n_qc_sfc_good,n_qc_sfc_bad)
 603  format(' # of SFC        GOOD/BAD QC = ',2i6,7x
     1      ,'% rejected = ',f6.1)

      write(6,604)n_qc_prof_good,n_qc_prof_bad
     1           ,pct_rejected(n_qc_prof_good,n_qc_prof_bad)
 604  format(' # of PROFs      GOOD/BAD QC = ',2i6,7x
     1      ,'% rejected = ',f6.1)

      write(6,605)n_qc_total_good,n_qc_total_bad
     1           ,pct_rejected(n_qc_total_good,n_qc_total_bad)
 605  format(/' # of Non-Radar  GOOD/BAD QC = ',2i6,7x
     1       ,'% rejected = ',f6.1)

      I4_elapsed = ishow_timer()

csms$serial end

      call get_inst_err2(r_missing_data                               ! I
     1                  ,obs_point_qced,max_obs,n_qc_total_good       ! I
     1                  ,rms_thresh_norm                              ! I
     1                  ,rms_inst,rms_thresh)                         ! O

      call barnes_multivariate(varbuff                                ! O
     1        ,n_var,ncnt_total,obs_point_qced                        ! I
     1        ,imax,jmax,kmax,grid_spacing_m,rep_pres_intvl           ! I
     1        ,varobs_diff_spread                                     ! O (aerr)
     1        ,wt_p,fnorm_dum,n_fnorm_dum                             ! I
     1        ,l_analyze_dum,.false.,rms_thresh,weight_bkg_const      ! I
     1        ,topo_dum,rland_frac_dum,1,1                            ! I
     1        ,n_obs_lvl,istatus)                                     ! O
      if(istatus .ne. 1)return

      wt_p_radar = wt_p

csms$print_mode(async) begin
csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 1 processor=',me

csms$serial(default=ignore)  begin              

      write(6,*)' Allocating upass1,vpass1'

      allocate( upass1(imax,jmax,kmax), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate upass1'
     1             ,istat_alloc,imax,jmax,kmax
          stop
      endif
!     call maxminavIJK(upass1,imax,jmax,kmax)

      allocate( vpass1(imax,jmax,kmax), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate vpass1'
     1             ,istat_alloc,imax,jmax,kmax
          stop
      endif
!     call maxminavIJK(vpass1,imax,jmax,kmax)

      call move_3d(varbuff(1,1,1,1),upass1,imax,jmax,kmax)
      call move_3d(varbuff(1,1,1,2),vpass1,imax,jmax,kmax)

      I4_elapsed = ishow_timer()

!     Perform radar QC by differencing radial velocities and first pass analysis
      do i_radar = 1,n_radars

          v_nyquist_global = v_nyquist_in(i_radar)
          
          write(6,*)' Radar QC for radar #/v_nyq/lat/lon ',i_radar
     1             ,v_nyquist_global
     1             ,rlat_radar(i_radar),rlon_radar(i_radar)
          call qc_radar_obs(
     1           imax,jmax,kmax                             ! Input
     1          ,r_missing_data                             ! Input
     1          ,vr_obs_unfltrd(1,1,1,i_radar)              ! Input/Output
     1          ,vr_nyq(1,1,1,i_radar)                      ! Input
     1          ,n_radarobs_tot_unfltrd(i_radar)            ! Input/Output
     1          ,lat,lon                                    ! Input
     1          ,rlat_radar(i_radar),rlon_radar(i_radar)    ! Input
     1          ,rheight_radar(i_radar)                     ! Input
     1          ,upass1,vpass1  ! 1st pass anal             ! Input
     1          ,u_laps_bkg,v_laps_bkg                      ! Input
     1          ,v_nyquist_global                           ! Input
     1          ,l_correct_unfolding,l_grid_north           ! Input
     1          ,istatus                                    ! Input/Output
     1                                                          )
      enddo ! i_radar

!     Perform analysis with radar data added in
      do k = 1,kmax
          l_analyze(k) = .false.
      enddo ! k

!     Fill 'obs_barnes' at this point in case there is no radar data
      obs_barnes = obs_point_qced
      ncnt_total = n_qc_total_good

csms$serial end

      if(n_radars .le. 1 .or. .not. l_3pass)then ! Single Doppler (or no radar) Option

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 2 processor=',me

csms$serial(<wt_p_radar , varobs_diff_spread,  
csms$>       icount_radar_total, out>:default=ignore)  begin

          mode = 1 ! All radar obs (in this case single Doppler)

!         Take the data from all the radars and add the derived radar obs into
!         uobs_diff_spread and vobs_diff_spread (varobs_diff_spread)
          if(l_point_struct)then
              wt_p_radar = r_missing_data                 ! Initialize
              varobs_diff_spread = r_missing_data         ! Initialize
          endif

          call insert_derived_radar_obs(
     1         mode                                       ! Input
     1        ,n_radars,max_radars,idx_radar_a            ! Input
     1        ,imax,jmax,kmax                             ! Input
     1        ,r_missing_data                             ! Input
     1        ,vr_obs_unfltrd                             ! Input
     1        ,thresh_2_radarobs_lvl_unfltrd              ! Input
     1        ,thresh_4_radarobs_lvl_unfltrd              ! Input
     1        ,i4time                                     ! Input
     1        ,lat,lon                                    ! Input
     1        ,rlat_radar,rlon_radar                      ! Input
     1        ,rheight_radar                              ! Input
     1        ,upass1,vpass1                              ! Input
     1        ,u_laps_bkg,v_laps_bkg                      ! Input
     1        ,weight_radar                               ! Input
     1        ,l_derived_output,l_grid_north              ! Input
     1        ,wt_p_radar                                 ! Input/Output
     1        ,varobs_diff_spread(1,1,1,1),varobs_diff_spread(1,1,1,2) ! I/O
     1        ,l_analyze,icount_radar_total               ! Output
     1        ,n_radarobs_tot_unfltrd                     ! Input
     1        ,istatus                                    ! Input/Output
     1                                                          )

csms$serial end

          if(icount_radar_total .gt. 0)then 

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 3 processor=',me

csms$serial(<rms_thresh, out>:default=ignore)  begin              
              I4_elapsed = ishow_timer()

              write(6,*)' Calling barnes with single radar obs added'       

csms$serial end

              if(l_point_struct)then
                  call arrays_to_barnesobs(imax,jmax,kmax             ! I
     1                              ,r_missing_data                   ! I
     1                              ,varobs_diff_spread,wt_p_radar    ! I
     1                              ,n_var,max_obs,obs_radar          ! I/O
     1                              ,ncnt_radar,weight_radar_total    ! O
     1                              ,istatus)                         ! O

!                 Combine radar (obs_radar) and non-radar (obs_point_qced) 
!                 data structures into new structure (obs_barnes)
                  obs_barnes = obs_point_qced
                  ncnt_total = n_qc_total_good
                  do i = 1,ncnt_radar
                      ncnt_total = ncnt_total + 1
                      obs_barnes(ncnt_total) = obs_radar(i)
                      obs_barnes(ncnt_total)%type = 'radar'
                  enddo ! i
              endif

              call get_inst_err2(r_missing_data                       ! I
     1                  ,obs_barnes,max_obs,ncnt_total                ! I
     1                  ,rms_thresh_norm                              ! I
     1                  ,rms_inst,rms_thresh)                         ! O

              call barnes_multivariate(varbuff                           ! O
     1           ,n_var,ncnt_total,obs_barnes                            ! I
     1           ,imax,jmax,kmax,grid_spacing_m,rep_pres_intvl           ! I
     1           ,varobs_diff_spread                                  ! O (aerr)
     1           ,wt_p_radar,fnorm_dum,n_fnorm_dum                       ! I
     1           ,l_analyze_dum,.false.,rms_thresh,weight_bkg_const      ! I
     1           ,topo_dum,rland_frac_dum,1,1                            ! I
     1           ,n_obs_lvl,istatus)                                     ! O

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 4 processor=',me

csms$serial(default=ignore)  begin              

              call move_3d(varbuff(1,1,1,1),uanl,imax,jmax,kmax)
              call move_3d(varbuff(1,1,1,2),vanl,imax,jmax,kmax)

              if(istatus .ne. 1)return

              I4_elapsed = ishow_timer()

csms$serial end

          endif ! There is any radar data

      else ! n_radars .gt. 1

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 5 processor=',me

csms$serial(<wt_p_radar, varobs_diff_spread, icount_radar_total, 
csms$>                    out>:default=ignore) begin       

          mode = 2 ! Only multi-Doppler obs

!         Take the data from all the radars and add the derived radar obs into
!         uobs_diff_spread and vobs_diff_spread (varobs_diff_spread)
          if(l_point_struct)then
              wt_p_radar = r_missing_data                 ! Initialize
              varobs_diff_spread = r_missing_data         ! Initialize
          endif

          call insert_derived_radar_obs(
     1         mode                                       ! Input
     1        ,n_radars,max_radars,idx_radar_a            ! Input
     1        ,imax,jmax,kmax                             ! Input
     1        ,r_missing_data                             ! Input
     1        ,vr_obs_unfltrd                             ! Input
     1        ,thresh_2_radarobs_lvl_unfltrd              ! Input
     1        ,thresh_4_radarobs_lvl_unfltrd              ! Input
     1        ,i4time                                     ! Input
     1        ,lat,lon                                    ! Input
     1        ,rlat_radar,rlon_radar                      ! Input
     1        ,rheight_radar                              ! Input
     1        ,upass1,vpass1                              ! Input
     1        ,u_laps_bkg,v_laps_bkg                      ! Input
     1        ,weight_radar                               ! Input
     1        ,l_derived_output,l_grid_north              ! Input
     1        ,wt_p_radar                                 ! Input/Output
     1        ,varobs_diff_spread(1,1,1,1),varobs_diff_spread(1,1,1,2) ! I/O
     1        ,l_analyze,icount_radar_total               ! Output
     1        ,n_radarobs_tot_unfltrd                     ! Input
     1        ,istatus                                    ! Input/Output
     1                                                          )

csms$serial end

          if(icount_radar_total .gt. 0)then 

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 6 processor=',me

csms$serial(<rms_thresh, out>:default=ignore)  begin              

              I4_elapsed = ishow_timer()

              write(6,401)icount_radar_total
 401          format(1x,' Analyzing with ',i5
     1                 ,' multi-doppler grid points')     

csms$serial end

              if(l_point_struct)then
                  call arrays_to_barnesobs(imax,jmax,kmax             ! I
     1                              ,r_missing_data                   ! I
     1                              ,varobs_diff_spread,wt_p_radar    ! I
     1                              ,n_var,max_obs,obs_radar          ! I/O
     1                              ,ncnt_radar,weight_radar_total    ! O
     1                              ,istatus)                         ! O

!                 Combine radar (obs_radar) and non-radar (obs_point_qced) 
!                 data structures into new structure (obs_barnes)
                  obs_barnes = obs_point_qced
                  ncnt_total = n_qc_total_good
                  do i = 1,ncnt_radar
                      ncnt_total = ncnt_total + 1
                      obs_barnes(ncnt_total) = obs_radar(i)
                      obs_barnes(ncnt_total)%type = 'radar'
                  enddo ! i

              endif

              call get_inst_err2(r_missing_data                         ! I
     1                  ,obs_barnes,max_obs,ncnt_total                  ! I
     1                  ,rms_thresh_norm                                ! I
     1                  ,rms_inst,rms_thresh)                           ! O

              call barnes_multivariate(varbuff                          ! O
     1          ,n_var,ncnt_total                                       ! I
     1          ,obs_barnes,imax,jmax,kmax,grid_spacing_m,rep_pres_intvl! I   
     1          ,varobs_diff_spread                                     ! O (aerr)
     1          ,wt_p_radar,fnorm_dum,n_fnorm_dum                       ! I
     1          ,l_analyze_dum,.false.,rms_thresh,weight_bkg_const      ! I
     1          ,topo_dum,rland_frac_dum,1,1                            ! I
     1          ,n_obs_lvl,istatus)                                     ! O

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 7 processor=',me

csms$serial(default=ignore)  begin              

              call move_3d(varbuff(1,1,1,1),uanl,imax,jmax,kmax)
              call move_3d(varbuff(1,1,1,2),vanl,imax,jmax,kmax)

              if(istatus .ne. 1)return

csms$serial end

          endif

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 8 processor=',me

csms$serial(<wt_p_radar , varobs_diff_spread, rms_thresh, out>
csms$>                                     :default=ignore)  begin

!         Make sure each level of uanl and vanl is initialized in the event it
!         was not analyzed.
          do k = 1,kmax
              if(.not. l_analyze(k)
     1               .or. icount_radar_total .eq. 0)then
                   write(6,412)k
412                format(' No multi-radar obs at lvl',i3,
     1                          ' Insert 1st pass into uanl,vanl')
                   do j=1,jmax
                   do i=1,imax
                       uanl(i,j,k) = upass1(i,j,k)
                       vanl(i,j,k) = vpass1(i,j,k)

                   enddo ! i
                   enddo ! j
              endif

          enddo ! k

          I4_elapsed = ishow_timer()

          mode = 1 ! Single and multi-Doppler obs

!         Take the data from all the radars and add the derived radar obs into
!         uobs_diff_spread and vobs_diff_spread
          if(l_point_struct)then
              wt_p_radar = r_missing_data                 ! Initialize
              varobs_diff_spread = r_missing_data         ! Initialize
          endif

          call insert_derived_radar_obs(
     1         mode                                       ! Input
     1        ,n_radars,max_radars,idx_radar_a            ! Input
     1        ,imax,jmax,kmax                             ! Input
     1        ,r_missing_data                             ! Input
     1        ,vr_obs_unfltrd                             ! Input
     1        ,thresh_2_radarobs_lvl_unfltrd              ! Input
     1        ,thresh_4_radarobs_lvl_unfltrd              ! Input
     1        ,i4time                                     ! Input
     1        ,lat,lon                                    ! Input
     1        ,rlat_radar,rlon_radar                      ! Input
     1        ,rheight_radar                              ! Input
     1        ,uanl,vanl                                  ! Input
     1        ,u_laps_bkg,v_laps_bkg                      ! Input
     1        ,weight_radar                               ! Input
     1        ,l_derived_output,l_grid_north              ! Input
     1        ,wt_p_radar                                 ! Input/Output
     1        ,varobs_diff_spread(1,1,1,1),varobs_diff_spread(1,1,1,2)  ! I/O
     1        ,l_analyze,icount_radar_total               ! Output
     1        ,n_radarobs_tot_unfltrd                     ! Input
     1        ,istatus                                    ! Input/Output
     1                                                          )

          I4_elapsed = ishow_timer()

          write(6,*)' Calling barnes with single+multi radar obs added'       

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 9 processor=',me
csms$serial end
csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 10 processor=',me

          if(l_point_struct)then
              call arrays_to_barnesobs(imax,jmax,kmax             ! I
     1                              ,r_missing_data                   ! I
     1                              ,varobs_diff_spread,wt_p_radar    ! I
     1                              ,n_var,max_obs,obs_radar          ! I/O
     1                              ,ncnt_radar,weight_radar_total    ! O
     1                              ,istatus)                         ! O

!             Combine radar (obs_radar) and non-radar (obs_point_qced) 
!             data structures into new structure (obs_barnes)
              obs_barnes = obs_point_qced
              ncnt_total = n_qc_total_good

              do i = 1,ncnt_radar
                  ncnt_total = ncnt_total + 1
                  obs_barnes(ncnt_total) = obs_radar(i)
                  obs_barnes(ncnt_total)%type = 'radar'
              enddo ! i

          endif

          call get_inst_err2(r_missing_data                           ! I
     1                  ,obs_barnes,max_obs,ncnt_total                ! I
     1                  ,rms_thresh_norm                              ! I
     1                  ,rms_inst,rms_thresh)                         ! O

          call barnes_multivariate(varbuff                            ! O
     1       ,n_var,ncnt_total,obs_barnes                             ! I
     1       ,imax,jmax,kmax                                          ! I
     1       ,grid_spacing_m,rep_pres_intvl                           ! I
     1       ,varobs_diff_spread                                      ! O (aerr)
     1       ,wt_p_radar,fnorm_dum,n_fnorm_dum                        ! I
     1       ,l_analyze_dum,.false.,rms_thresh,weight_bkg_const       ! I
     1       ,topo_dum,rland_frac_dum,1,1                             ! I
     1       ,n_obs_lvl,istatus)                                      ! O

csms$insert      call nnt_me(me)
csms$insert      print *, 'got to 11 processor=',me

csms$serial(default=ignore)  begin              

          call move_3d(varbuff(1,1,1,1),uanl,imax,jmax,kmax)
          call move_3d(varbuff(1,1,1,2),vanl,imax,jmax,kmax)

          if(istatus .ne. 1)return

          I4_elapsed = ishow_timer()

csms$serial end

      endif ! n_radars

csms$serial(default=ignore)  begin              

      write(6,*)' Adding analyzed differences to the background '
     1         ,'to reconstruct full analyses'

      do k=1,kmax ! Add back differences for first pass

          do j=1,jmax
          do i=1,imax
              if(upass1(i,j,k) .ne. r_missing_data)then
                  upass1(i,j,k) = upass1(i,j,k) + u_laps_bkg(i,j,k)
                  vpass1(i,j,k) = vpass1(i,j,k) + v_laps_bkg(i,j,k)
              else
                  write(6,*)
     1            ' ERROR: Missing data value(s) detected in first'
     1           ,' pass at lvl',k
                  istatus = 0
                  return
              endif
          enddo ! i
          enddo ! j


          if(l_analyze(k) .OR. (.true. .and. icount_radar_total .gt. 0) ! l_3d
     1                         )then ! This depends on the presence of radar obs
              write(6,511)k
511           format(' Use 2nd Pass for lvl',i3)

              do j=1,jmax
              do i=1,imax
                  uanl(i,j,k) = uanl(i,j,k) + u_laps_bkg(i,j,k)
                  vanl(i,j,k) = vanl(i,j,k) + v_laps_bkg(i,j,k)

              enddo ! i
              enddo ! j

          else
              write(6,512)k
512           format(' Use 1st Pass for lvl',i3)
              do j=1,jmax
              do i=1,imax
                  uanl(i,j,k) = upass1(i,j,k)
                  vanl(i,j,k) = vpass1(i,j,k)

              enddo ! i
              enddo ! j

          endif
      enddo ! k

      if(iter .eq. n_iter_wind)then
!         Compare 1st pass analysis to obs
          call compare_wind(
     1            upass1,vpass1,'PS1 ',
     1            istat_radar_vel,max_radars,vr_obs_unfltrd,n_radars,
     1            rlat_radar,rlon_radar,rheight_radar,
     1            lat,lon,
     1            imax,jmax,kmax,r_missing_data,
     1            obs_barnes,max_obs,ncnt_total,l_point_struct,
     1            weight_pirep,weight_prof,weight_sfc,weight_cdw,
     1            uobs,vobs,wt_p,istatus)

      endif ! Last iteration

      write(6,*)' Deallocate upass1, vpass1'
      deallocate(upass1)
      deallocate(vpass1)

csms$serial end

      enddo ! n_iter_wind

csms$serial(default=ignore)  begin              

!     Compare final analysis to obs
      call compare_wind(
     1            uanl,vanl,'LAPS',
     1            istat_radar_vel,max_radars,vr_obs_unfltrd,n_radars,
     1            rlat_radar,rlon_radar,rheight_radar,
     1            lat,lon,
     1            imax,jmax,kmax,r_missing_data,
     1            obs_barnes,max_obs,ncnt_total,l_point_struct,
     1            weight_pirep,weight_prof,weight_sfc,weight_cdw,
     1            uobs,vobs,wt_p,istatus)

      istatus = 1

      write(6,*)' End of subroutine laps_anl'

csms$serial end

      return
      end


      subroutine get_inst_err(imax,jmax,kmax,r_missing_data        ! I
     1                       ,wt_p,rms_thresh_norm                 ! I
     1                       ,rms_inst,rms_thresh)                 ! O

      real*4    wt_p(imax,jmax,kmax)                               ! Input

csms$ignore begin

      write(6,*)
      write(6,*)' subroutine get_inst_err...'

      n_obs_total = 0
      wt_p_inv_total = 0.

      do i = 1,imax
      do j = 1,jmax
      do k = 1,kmax
          if(wt_p(i,j,k) .ne. r_missing_data)then
              n_obs_total = n_obs_total + 1
              wt_p_inv_total = wt_p_inv_total + 1.0 / wt_p(i,j,k)
          endif
      enddo ! k
      enddo ! j
      enddo ! i

      if(n_obs_total .gt. 0)then
          wt_p_inv_ave = wt_p_inv_total / float(n_obs_total)
          rms_inst = sqrt(wt_p_inv_ave)
      else
          wt_p_inv_ave = 0.
          rms_inst = 0.
      endif

      rms_thresh = rms_inst * rms_thresh_norm

      write(6,*)' n_obs_total = ',n_obs_total
      write(6,*)' wt_p_inv_total,wt_p_inv_ave = '
     1           ,wt_p_inv_total,wt_p_inv_ave
      write(6,*)' rms_inst, rms_thresh = ',rms_inst,rms_thresh
csms$ignore end

      return
      end


      subroutine get_inst_err2(r_missing_data                       ! I
     1                        ,obs_barnes,max_obs,nobs_barnes       ! I
     1                        ,rms_thresh_norm                      ! I
     1                        ,rms_inst,rms_thresh)                 ! O

      integer*4 max_obs
      include 'barnesob.inc'
      type (barnesob) obs_barnes(max_obs)                           

csms$ignore begin

      write(6,*)
      write(6,*)' subroutine get_inst_err2...'

      n_obs_total = 0
      wt_p_inv_total = 0.

      do i = 1,nobs_barnes
          n_obs_total = n_obs_total + 1
          wt_p_inv_total = wt_p_inv_total + 1.0 / obs_barnes(i)%weight       
      enddo ! i

      if(n_obs_total .gt. 0)then
          wt_p_inv_ave = wt_p_inv_total / float(n_obs_total)
          rms_inst = sqrt(wt_p_inv_ave)
      else
          wt_p_inv_ave = 0.
          rms_inst = 0.
      endif

      rms_thresh = rms_inst * rms_thresh_norm

      write(6,*)' n_obs_total = ',n_obs_total
      write(6,*)' wt_p_inv_total,wt_p_inv_ave = '
     1           ,wt_p_inv_total,wt_p_inv_ave
      write(6,*)' rms_inst, rms_thresh = ',rms_inst,rms_thresh
csms$ignore end

      return
      end

