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

       subroutine laps_anl(uobs,vobs,n_radars
     1     ,istat_radar_vel                                      ! Input
     1     ,vr_obs_unfltrd,vr_nyq,v_nyquist_in
!    1     ,upass1,vpass1                                        ! Output
     1     ,n_var,n_fnorm                                        ! Input
     1     ,uanl,vanl                                            ! Output
     1     ,wt_p,weight_bkg_const,rms_thresh_wind                ! Input/Local
     1     ,max_radars
     1     ,n_radarobs_tot_unfltrd,rlat_radar,rlon_radar,rheight_radar
     1     ,u_laps_bkg,v_laps_bkg                                ! Input/Local
     1     ,imax,jmax,kmax,lat,lon
     1     ,i4time,grid_spacing_m
     1     ,r_missing_data
     1     ,i_3d                                                 ! Input
     1     ,l_derived_output,l_grid_north,l_3pass,l_correct_unfolding
     1     ,n_iter_wind_in
     1     ,weight_cdw,weight_sfc,weight_pirep,weight_prof,weight_radar
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
      parameter (max_obs = 40000)       
      include 'barnesob.inc'
      type (barnesob) obs_barnes(max_obs)                           

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
      real*4    wt_p_spread(imax,jmax,kmax)                        ! Local

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

!     Nyquist velocity (if known and constant) for each radar
      real*4 v_nyquist_in(max_radars)                                  ! Input

!     Location of each radar
      real*4 rlat_radar(max_radars),rlon_radar(max_radars)             ! Input
     1                     ,rheight_radar(max_radars)

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
      logical  l_3d,l_not_struct

      real*4   rms_thresh                                              ! Input
      real*4   weight_bkg_const                                        ! Input

!     These are the weights of the various data types (filling the 3D array)
      real*4 weight_cdw,weight_sfc,weight_pirep,weight_prof
     1      ,weight_radar ! Input

      integer*4 istatus         ! (1 is good)                          ! Output

      integer*4  n_fnorm

      character*3 c3_string

      dimension fnorm(0:n_fnorm)

!****************END DECLARATIONS *********************************************

      write(6,*)' Subroutine laps_anl...'

csms$serial(default=ignore)  begin              

!     Compare background to obs
      call compare_wind(
     1            u_laps_bkg,v_laps_bkg,' FG ',
     1            istat_radar_vel,max_radars,vr_obs_unfltrd,n_radars,
     1            rlat_radar,rlon_radar,rheight_radar,
     1            lat,lon,
     1            imax,jmax,kmax,r_missing_data,
     1            weight_pirep,weight_prof,weight_sfc,weight_cdw,
     1            uobs,vobs,wt_p,istatus)
csms$serial end

      l_3d = .true.
      l_not_struct = .false.

      rms_thresh_norm = rms_thresh_wind          ! Not used if l_3d = .false.

      n_iter_wind = 1

      do iter = 1,n_iter_wind

csms$serial(<varobs_diff_spread, wt_p_spread,
csms$>       fnorm, l_analyze, rms_thresh , out>:default=ignore)  begin
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

      n_qc_pirep_bad = 0
      n_qc_cdw_bad = 0
      n_qc_sfc_bad = 0
      n_qc_prof_bad = 0
      n_qc_total_bad = 0

      n_qc_pirep_good = 0
      n_qc_cdw_good = 0
      n_qc_sfc_good = 0
      n_qc_prof_good = 0
      n_qc_total_good = 0

      iwrite = 1

      qc_thresh = 30. ! Threshold speed for throwing out the ob

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

      allocate( pres_3d(imax,jmax,kmax), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate pres_3d'
          stop
      endif

      call get_pres_3d(i4time,imax,jmax,kmax,pres_3d,istatus)
      if(istatus .ne. 1)return

      call get_rep_pres_intvl(pres_3d,imax,jmax,kmax,rep_pres_intvl
     1                       ,istatus)

      do j=1,jmax
      do i=1,imax
        do k = 1,kmax
!         wt_p_spread(i,j,k) = wt_p(i,j,k)        ! Initialize the local array

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

!       1       .OR. (abs(k-16) .le. 5 .and. (i .ne. 50 .or. j .ne. 14)) )
     1                                                          )then

                ! Throw out the ob
                  if(wt_p(i,j,k) .eq. weight_pirep)then
                      n_qc_pirep_bad = n_qc_pirep_bad + 1
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
                      n_qc_pirep_good = n_qc_pirep_good + 1
                      c3_string = 'Prp'
                      n_good_thistype = n_qc_pirep_good
                  endif

                  if(wt_p(i,j,k) .eq. weight_cdw)then
                      n_qc_cdw_good  = n_qc_cdw_good + 1
                      c3_string = 'Cdw'
                      n_good_thistype = n_qc_cdw_good
                  endif

                  if(wt_p(i,j,k) .eq. weight_sfc)then
                      n_qc_sfc_good   = n_qc_sfc_good + 1
                      c3_string = 'Sfc'
                      n_good_thistype = n_qc_sfc_good
                  endif

                  if(wt_p(i,j,k) .eq. weight_prof)then
                      n_qc_prof_good  = n_qc_prof_good + 1
                      c3_string = 'Prf'
                      n_good_thistype = n_qc_prof_good
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

!       Spread the difference ob vertically
        call spread_vert(uobs_diff,vobs_diff,l_3d,iwrite
     1          ,varobs_diff_spread(1,1,1,1),varobs_diff_spread(1,1,1,2)
     1          ,wt_p,wt_p_spread,pres_3d,i,j,imax,jmax,kmax,istatus)

        if(istatus .ne. 1)return

      enddo ! i
      enddo ! j

      deallocate(uobs_diff)
      deallocate(vobs_diff)
      deallocate(pres_3d)

      write(6,*)
      write(6,*)' QC info for non-radar data (after remapping to grid)'
      write(6,601)n_qc_pirep_good,n_qc_pirep_bad
     1           ,pct_rejected(n_qc_pirep_good,n_qc_pirep_bad)
 601  format(' # of PIREPs     GOOD/BAD QC = ',2i6,7x
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

!     Perform 1st pass analysis of non-radar difference obs
      do k = 1,kmax
          l_analyze(k) = .true.
      enddo ! k

      call get_inst_err(imax,jmax,kmax,r_missing_data
     1        ,wt_p_spread,rms_thresh_norm,rms_inst,rms_thresh)

csms$serial end

      call arrays_to_barnesobs      (imax,jmax,kmax                   ! I
     1                              ,r_missing_data                   ! I
     1                              ,varobs_diff_spread,wt_p_spread   ! I
     1                              ,n_var,max_obs,obs_barnes         ! I/O
     1                              ,ncnt_total,weight_total          ! O
     1                              ,istatus)                         ! O

      call barnes_multivariate(varbuff,n_var,ncnt_total,obs_barnes
!     call barnes_multivariate(varbuff,n_var,max_obs,obs_barnes
     1        ,imax,jmax,kmax,grid_spacing_m,rep_pres_intvl
     1        ,varobs_diff_spread
     1        ,wt_p_spread,fnorm,n_fnorm
     1        ,l_analyze,l_not_struct,rms_thresh,weight_bkg_const
     1        ,n_obs_lvl,istatus)
      if(istatus .ne. 1)return
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
          if(v_nyquist_in(i_radar) .ne. r_missing_data
     1       .and. v_nyquist_in(i_radar).lt.200.
     1       .and. v_nyquist_in(i_radar).gt.1.0 ) then
              v_nyquist_2_in = 2. * v_nyquist_in(i_radar)
              unfolding_thresh_in = 0.7 * v_nyquist_in(i_radar)
          else
              v_nyquist_2_in = r_missing_data
          endif
          write(6,*)' Radar QC for radar #/v_nyq*2 '
     1             ,i_radar,v_nyquist_2_in
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
     1          ,v_nyquist_2_in,unfolding_thresh_in         ! Input
     1          ,l_correct_unfolding,l_grid_north           ! Input
     1          ,istatus                                    ! Input/Output
     1                                                          )
      enddo

!     Perform analysis with radar data added in
      do k = 1,kmax
          l_analyze(k) = .false.
      enddo ! k

      if(n_radars .le. 1 .or. .not. l_3pass)then ! Single Doppler (or no radar) Option

          mode = 1 ! All radar obs (in this case single Doppler)

!         Take the data from all the radars and add the derived radar obs into
!         uobs_diff_spread and vobs_diff_spread (varobs_diff_spread)
          call insert_derived_radar_obs(
     1         mode                                       ! Input
     1        ,n_radars,max_radars                        ! Input
     1        ,imax,jmax,kmax                             ! Input
     1        ,r_missing_data                             ! Input
     1        ,vr_obs_unfltrd                             ! Input
     1        ,i4time                                     ! Input
     1        ,lat,lon                                    ! Input
     1        ,rlat_radar,rlon_radar                      ! Input
     1        ,rheight_radar                              ! Input
     1        ,upass1,vpass1                              ! Input
     1        ,u_laps_bkg,v_laps_bkg                      ! Input
     1        ,weight_radar                               ! Input
     1        ,l_derived_output,l_grid_north              ! Input
     1        ,wt_p_spread                                ! Input/Output
     1        ,varobs_diff_spread(1,1,1,1),varobs_diff_spread(1,1,1,2) ! I/O
     1        ,l_analyze,icount_radar_total               ! Output
     1        ,n_radarobs_tot_unfltrd                     ! Input
     1        ,istatus                                    ! Input/Output
     1                                                          )
          if(icount_radar_total .gt. 0 .or. .not. .true.)then ! l_3d

csms$insert      print *, 'Error Parallelization not done for this code section'
csms$insert      stop
              I4_elapsed = ishow_timer()

              write(6,*)' Calling barnes with modified radar obs added'

              call move_3d(uanl,varbuff(1,1,1,1),imax,jmax,kmax)
              call move_3d(vanl,varbuff(1,1,1,2),imax,jmax,kmax)

              call get_inst_err(imax,jmax,kmax,r_missing_data
     1            ,wt_p_spread,rms_thresh_norm,rms_inst,rms_thresh)

              call arrays_to_barnesobs(imax,jmax,kmax                 ! I
     1                              ,r_missing_data                   ! I
     1                              ,varobs_diff_spread,wt_p_spread   ! I
     1                              ,n_var,max_obs,obs_barnes         ! I/O
     1                              ,ncnt_total,weight_total          ! O
     1                              ,istatus)                         ! O

              call barnes_multivariate
     1                             (varbuff,n_var,ncnt_total,obs_barnes
!             call barnes_multivariate(varbuff,n_var,max_obs,obs_barnes       
     1           ,imax,jmax,kmax,grid_spacing_m,rep_pres_intvl
     1           ,varobs_diff_spread
     1           ,wt_p_spread,fnorm,n_fnorm
     1           ,l_analyze,l_not_struct,rms_thresh,weight_bkg_const
     1           ,n_obs_lvl,istatus)

              call move_3d(varbuff(1,1,1,1),uanl,imax,jmax,kmax)
              call move_3d(varbuff(1,1,1,2),vanl,imax,jmax,kmax)

              if(istatus .ne. 1)return

              I4_elapsed = ishow_timer()

          endif ! There is any radar data

      else ! n_radars .gt. 1

csms$insert      print *, 'Error Parallelization not done for this code section'
csms$insert      stop
          mode = 2 ! Only multi-Doppler obs

!         Take the data from all the radars and add the derived radar obs into
!         uobs_diff_spread and vobs_diff_spread (varobs_diff_spread)
          call insert_derived_radar_obs(
     1         mode                                       ! Input
     1        ,n_radars,max_radars                        ! Input
     1        ,imax,jmax,kmax                             ! Input
     1        ,r_missing_data                             ! Input
     1        ,vr_obs_unfltrd                             ! Input
     1        ,i4time                                     ! Input
     1        ,lat,lon                                    ! Input
     1        ,rlat_radar,rlon_radar                      ! Input
     1        ,rheight_radar                              ! Input
     1        ,upass1,vpass1                              ! Input
     1        ,u_laps_bkg,v_laps_bkg                      ! Input
     1        ,weight_radar                               ! Input
     1        ,l_derived_output,l_grid_north              ! Input
     1        ,wt_p_spread                                ! Input/Output
     1        ,varobs_diff_spread(1,1,1,1),varobs_diff_spread(1,1,1,2) ! I/O
     1        ,l_analyze,icount_radar_total               ! Output
     1        ,n_radarobs_tot_unfltrd                     ! Input
     1        ,istatus                                    ! Input/Output
     1                                                          )

          if(icount_radar_total .gt. 0 .or. .not. .true.)then ! l_3d

          I4_elapsed = ishow_timer()

          write(6,*)' Calling barnes with only multi-doppler obs '
     1          ,'creating an intermediate analysis'

          call move_3d(uanl,varbuff(1,1,1,1),imax,jmax,kmax)
          call move_3d(vanl,varbuff(1,1,1,2),imax,jmax,kmax)

          call get_inst_err(imax,jmax,kmax,r_missing_data
     1        ,wt_p_spread,rms_thresh_norm,rms_inst,rms_thresh)

          call arrays_to_barnesobs  (imax,jmax,kmax                   ! I
     1                              ,r_missing_data                   ! I
     1                              ,varobs_diff_spread,wt_p_spread   ! I
     1                              ,n_var,max_obs,obs_barnes         ! I/O
     1                              ,ncnt_total,weight_total          ! O
     1                              ,istatus)                         ! O

          call barnes_multivariate(varbuff,n_var,ncnt_total,obs_barnes       
!         call barnes_multivariate(varbuff,n_var,max_obs,obs_barnes
     1       ,imax,jmax,kmax,grid_spacing_m,rep_pres_intvl
     1       ,varobs_diff_spread
     1       ,wt_p_spread,fnorm,n_fnorm
     1       ,l_analyze,l_not_struct,rms_thresh,weight_bkg_const
     1       ,n_obs_lvl,istatus)

          call move_3d(varbuff(1,1,1,1),uanl,imax,jmax,kmax)
          call move_3d(varbuff(1,1,1,2),vanl,imax,jmax,kmax)

          if(istatus .ne. 1)return

          endif

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
          call insert_derived_radar_obs(
     1         mode                                       ! Input
     1        ,n_radars,max_radars                        ! Input
     1        ,imax,jmax,kmax                             ! Input
     1        ,r_missing_data                             ! Input
     1        ,vr_obs_unfltrd                             ! Input
     1        ,i4time                                     ! Input
     1        ,lat,lon                                    ! Input
     1        ,rlat_radar,rlon_radar                      ! Input
     1        ,rheight_radar                              ! Input
     1        ,uanl,vanl                                  ! Input
     1        ,u_laps_bkg,v_laps_bkg                      ! Input
     1        ,weight_radar                               ! Input
     1        ,l_derived_output,l_grid_north              ! Input
     1        ,wt_p_spread                                ! Input/Output
     1        ,varobs_diff_spread(1,1,1,1),varobs_diff_spread(1,1,1,2)  ! I/O
     1        ,l_analyze,icount_radar_total               ! Output
     1        ,n_radarobs_tot_unfltrd                     ! Input
     1        ,istatus                                    ! Input/Output
     1                                                          )

          I4_elapsed = ishow_timer()

          write(6,*)' Calling barnes with modified radar obs added'

          call move_3d(uanl,varbuff(1,1,1,1),imax,jmax,kmax)
          call move_3d(vanl,varbuff(1,1,1,2),imax,jmax,kmax)

          call get_inst_err(imax,jmax,kmax,r_missing_data
     1        ,wt_p_spread,rms_thresh_norm,rms_inst,rms_thresh)

          call arrays_to_barnesobs  (imax,jmax,kmax                   ! I
     1                              ,r_missing_data                   ! I
     1                              ,varobs_diff_spread,wt_p_spread   ! I
     1                              ,n_var,max_obs,obs_barnes         ! I/O
     1                              ,ncnt_total,weight_total          ! O
     1                              ,istatus)                         ! O

          call barnes_multivariate(varbuff,n_var,ncnt_total,obs_barnes
!         call barnes_multivariate(varbuff,n_var,max_obs,obs_barnes
     1       ,imax,jmax,kmax
     1       ,grid_spacing_m,rep_pres_intvl
     1       ,varobs_diff_spread
     1       ,wt_p_spread,fnorm,n_fnorm
     1       ,l_analyze,l_not_struct,rms_thresh,weight_bkg_const
     1       ,n_obs_lvl,istatus)

          call move_3d(varbuff(1,1,1,1),uanl,imax,jmax,kmax)
          call move_3d(varbuff(1,1,1,2),vanl,imax,jmax,kmax)

          if(istatus .ne. 1)return

          I4_elapsed = ishow_timer()

      endif ! n_radars

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
!         write(6,*)' Calling comparisons'
!         call comparisons(
!    1            upass1,vpass1,istat_radar_vel,max_radars,
!    1            vr_obs_unfltrd,
!    1            rlat_radar,rlon_radar,rheight_radar,
!    1            lat,lon,
!    1            uanl,vanl,u_laps_bkg,v_laps_bkg,
!    1            istat_bal,
!    1            imax,jmax,kmax,r_missing_data,
!    1            weight_pirep,weight_prof,weight_sfc,weight_cdw,
!    1            uobs,vobs,wt_p,
!    1            n_radars)

!         Compare 1st pass analysis to obs
          call compare_wind(
     1            upass1,vpass1,'PS1 ',
     1            istat_radar_vel,max_radars,vr_obs_unfltrd,n_radars,
     1            rlat_radar,rlon_radar,rheight_radar,
     1            lat,lon,
     1            imax,jmax,kmax,r_missing_data,
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
     1            weight_pirep,weight_prof,weight_sfc,weight_cdw,
     1            uobs,vobs,wt_p,istatus)

      istatus = 1

      write(6,*)' End of subroutine laps_anl'

csms$serial end

      return
      end


      subroutine insert_derived_radar_obs(
     1   mode                                       ! Input
     1  ,n_radars,max_radars                        ! Input
     1  ,imax,jmax,kmax                             ! Input
     1  ,r_missing_data                             ! Input
     1  ,vr_obs_unfltrd                             ! Input
     1  ,i4time                                     ! Input
     1  ,lat,lon                                    ! Input
     1  ,rlat_radar,rlon_radar                      ! Input
     1  ,rheight_radar                              ! Input
     1  ,upass1,vpass1                              ! Input
     1  ,u_laps_bkg,v_laps_bkg                      ! Input
     1  ,weight_radar                               ! Input
     1  ,l_derived_output,l_grid_north              ! Input
     1  ,wt_p_spread                                ! Input/Output
     1  ,uobs_diff_spread,vobs_diff_spread          ! Input/Output
     1  ,l_analyze,icount_radar_total               ! Output
     1  ,n_radarobs_tot_unfltrd                     ! Input
     1  ,istatus                                    ! Input/Output
     1                                                          )

      real*4   vr_obs_unfltrd(imax,jmax,kmax,max_radars)
      real*4   rlat_radar(max_radars),rlon_radar(max_radars)
      real*4   rheight_radar(max_radars)
      real*4   n_radarobs_tot_unfltrd(max_radars)
      real*4   lat(imax,jmax),lon(imax,jmax)
      real*4   upass1(imax,jmax,kmax),vpass1(imax,jmax,kmax)
      real*4   u_laps_bkg(imax,jmax,kmax),v_laps_bkg(imax,jmax,kmax)
      real*4   uobs_diff_spread(imax,jmax,kmax)
     1        ,vobs_diff_spread(imax,jmax,kmax)
      real*4   wt_p_spread(imax,jmax,kmax)

      real*4   vr_obs_fltrd(imax,jmax,kmax)                          ! Local
      real*4   upass1_buf(imax,jmax,kmax)                            ! Local
      real*4   vpass1_buf(imax,jmax,kmax)                            ! Local
      logical  l_good_multi_doppler_ob(imax,jmax,kmax)               ! Local
      logical  l_analyze(kmax),l_derived_output,l_grid_north

csms$ignore begin
      write(6,*)' Entering insert_derived_radar_obs, mode =',mode

!     Initialize l_good_multi_doppler_ob if mode = 2
      if(mode .eq. 2)then
          do k = 1,kmax
          do j = 1,jmax
          do i = 1,imax
              l_good_multi_doppler_ob(i,j,k) = .false.
          enddo ! i
          enddo ! j
          enddo ! k
      endif

!     This routine takes the data from all the radars and adds the derived
!     radar obs into uobs_diff_spread and vobs_diff_spread

!     Loop through the radars
      do i_radar = 1,n_radars

          write(6,*)' Generate derived radar "vector" obs, radar # ',i_r
     1adar

          if(i_radar .eq. 1)then ! Use original 1st pass arrays
              call make_derived_radar_obs(
     1   imax,jmax,kmax                             ! Input
     1  ,mode                                       ! Input
     1  ,r_missing_data                             ! Input
     1  ,i_radar                                    ! Input
     1  ,vr_obs_unfltrd(1,1,1,i_radar)              ! Input
     1  ,i4time                                     ! Input
     1  ,lat,lon                                    ! Input
     1  ,rlat_radar(i_radar),rlon_radar(i_radar)    ! Input
     1  ,rheight_radar(i_radar)                     ! Input
     1  ,upass1,vpass1  ! 1st pass anal             ! Input
     1  ,u_laps_bkg,v_laps_bkg                      ! Input
     1  ,weight_radar                               ! Input
     1  ,l_derived_output,l_grid_north              ! Input
     1  ,wt_p_spread                                ! Input/Output
     1  ,uobs_diff_spread,vobs_diff_spread          ! Input/Output
     1  ,n_radarobs_tot_unfltrd(i_radar)            ! Input
     1  ,vr_obs_fltrd                               ! Local
     1  ,l_good_multi_doppler_ob                    ! Input/Output
     1  ,istatus                                    ! Input/Output
     1                                                          )
          else ! Use buffer 1st pass arrays
              call make_derived_radar_obs(
     1   imax,jmax,kmax                             ! Input
     1  ,mode                                       ! Input
     1  ,r_missing_data                             ! Input
     1  ,i_radar                                    ! Input
     1  ,vr_obs_unfltrd(1,1,1,i_radar)              ! Input
     1  ,i4time                                     ! Input
     1  ,lat,lon                                    ! Input
     1  ,rlat_radar(i_radar),rlon_radar(i_radar)    ! Input
     1  ,rheight_radar(i_radar)                     ! Input
     1  ,upass1_buf,vpass1_buf ! 1st pass anal + derived radar ! Input
     1  ,u_laps_bkg,v_laps_bkg                      ! Input
     1  ,weight_radar                               ! Input
     1  ,l_derived_output,l_grid_north              ! Input
     1  ,wt_p_spread                                ! Input/Output
     1  ,uobs_diff_spread,vobs_diff_spread          ! Input/Output
     1  ,n_radarobs_tot_unfltrd(i_radar)            ! Input
     1  ,vr_obs_fltrd                               ! Local
     1  ,l_good_multi_doppler_ob                    ! Input/Output
     1  ,istatus                                    ! Input/Output
     1                                                          )
          endif

        ! Should we insert the derived radar obs into the 1st pass buffer arrays?
          if(i_radar .lt. n_radars)then

              if(i_radar .eq. 1)then ! Initialize the 1st pass buffer arrays

                  do k = 1,kmax
                  do j = 1,jmax
                  do i = 1,imax
                      upass1_buf(i,j,k) = upass1(i,j,k)
                      vpass1_buf(i,j,k) = vpass1(i,j,k)
                  enddo ! i
                  enddo ! j
                  enddo ! k

              endif

            ! Insert the derived radar obs into the 1st pass buffer arrays
              do k = 1,kmax
              do j = 1,jmax
              do i = 1,imax
                  if(wt_p_spread(i,j,k) .eq. weight_radar)then
                      upass1_buf(i,j,k) = uobs_diff_spread(i,j,k)
                      vpass1_buf(i,j,k) = vobs_diff_spread(i,j,k)
                  endif
              enddo ! i
              enddo ! j
              enddo ! k
          endif

      enddo ! i_radar

      icount_radar_total = 0

!     Use only multiple Doppler obs if mode = 2, use all obs if mode = 1
      do k = 1,kmax
          icount_good_lvl = 0
          if(mode .eq. 1)then ! Use all Doppler obs
              do j = 1,jmax
              do i = 1,imax
                  if(wt_p_spread(i,j,k) .eq. weight_radar)then ! Good single or multi ob
                      icount_good_lvl = icount_good_lvl + 1
                  endif
              enddo ! i
              enddo ! j
          elseif(mode .eq. 2)then ! Throw out single Doppler obs
              do j = 1,jmax
              do i = 1,imax
                  if(wt_p_spread(i,j,k) .eq. weight_radar)then ! Good single or multi ob
                      if(.not. l_good_multi_doppler_ob(i,j,k))then
                          uobs_diff_spread(i,j,k) = r_missing_data
                          vobs_diff_spread(i,j,k) = r_missing_data
                          wt_p_spread(i,j,k) = r_missing_data
                      else
                          icount_good_lvl = icount_good_lvl + 1
                      endif
                  endif
              enddo ! i
              enddo ! j
          endif ! mode

          if(icount_good_lvl .gt. 0)l_analyze(k) = .true.

          if(mode .eq. 1)then ! Use all Doppler obs
              write(6,504)k,icount_good_lvl,l_analyze(k)
504           format(' LVL',i3,' # sngl+multi = ',i6,l2)
          elseif(mode .eq. 2)then ! Use only multi Doppler obs
              write(6,505)k,icount_good_lvl,l_analyze(k)
505           format(' LVL',i3,' # multi = ',i6,l2)
          endif

          icount_radar_total = icount_radar_total + icount_good_lvl

      enddo ! k

      write(6,*)
     1     ' Finished insert_derived_radar_obs, icount_radar_total = '
     1    ,icount_radar_total

csms$ignore end
      return
      end


      subroutine make_derived_radar_obs(
     1   imax,jmax,kmax                             ! Input
     1  ,mode                                       ! Input
     1  ,r_missing_data                             ! Input
     1  ,i_radar                                    ! Input
     1  ,vr_obs_unfltrd                             ! Input
     1  ,i4time                                     ! Input
     1  ,lat,lon                                    ! Input
     1  ,rlat_radar,rlon_radar,rheight_radar        ! Input
     1  ,upass1,vpass1                              ! Input
     1  ,u_laps_bkg,v_laps_bkg                      ! Input
     1  ,weight_radar                               ! Input
     1  ,l_derived_output,l_grid_north              ! Input
     1  ,wt_p_spread                                ! Input/Output
     1  ,uobs_diff_spread,vobs_diff_spread          ! Input/Output
     1  ,n_radarobs_tot_unfltrd                     ! Input
     1  ,vr_obs_fltrd                               ! Local
     1  ,l_good_multi_doppler_ob                    ! Input/Output
     1  ,istatus                                    ! Input/Output
     1                                                          )

      real*4   vr_obs_unfltrd(imax,jmax,kmax)
      real*4   lat(imax,jmax),lon(imax,jmax)
      real*4   upass1(imax,jmax,kmax),vpass1(imax,jmax,kmax)
      real*4   u_laps_bkg(imax,jmax,kmax),v_laps_bkg(imax,jmax,kmax)
      real*4   uobs_diff_spread(imax,jmax,kmax)
     1        ,vobs_diff_spread(imax,jmax,kmax)
      real*4   wt_p_spread(imax,jmax,kmax)
      real*4   vr_obs_fltrd(imax,jmax,kmax)

      logical  l_good_multi_doppler_ob(imax,jmax,kmax),l_derived_output
      logical  l_grid_north

      character*31 ext

csms$ignore begin
      write(6,*)' Filtering radar obs into superobs',rlat_radar,rlon_rad
     1ar
      n_radarobs_tot_fltrd = 0
      i_radar_reject = 0

      write(6,*)' LVL  # Obs  Intvl  # FLTR'

      do k=1,kmax

!       Count number of unfiltered obs after rejecting obs having non-radar data
        n_radarobs_lvl_unfltrd = 0
        do j=1,jmax
        do i=1,imax
            if(vr_obs_unfltrd(i,j,k) .ne. r_missing_data)then
                if(wt_p_spread(i,j,k) .ne. weight_radar .and.
     1     wt_p_spread(i,j,k) .ne. r_missing_data)then ! Non-radar ob
                    vr_obs_unfltrd(i,j,k) = r_missing_data
                    i_radar_reject = i_radar_reject + 1
                else
                    n_radarobs_lvl_unfltrd = n_radarobs_lvl_unfltrd + 1
                endif
            endif
        enddo ! i
        enddo ! j

!       Filter to make the filtered ob array more sparse
        call filter_radar_obs(
     1                  n_radarobs_lvl_unfltrd,       ! Input
     1                  imax,jmax,                    ! Input
     1                  vr_obs_unfltrd(1,1,k),        ! Input
     1                  r_missing_data,               ! Input
     1                  vr_obs_fltrd(1,1,k),          ! Input/Output
     1                  intvl_rad)                    ! Output

!       Count number of filtered obs
        n_radarobs_lvl_fltrd = 0
        do j=1,jmax
        do i=1,imax
          if(vr_obs_fltrd(i,j,k) .ne. r_missing_data)then
            n_radarobs_lvl_fltrd = n_radarobs_lvl_fltrd + 1
          endif
        enddo ! i
        enddo ! j

        n_radarobs_tot_fltrd = n_radarobs_tot_fltrd 
     1                       + n_radarobs_lvl_fltrd

        write(6,501)k,n_radarobs_lvl_unfltrd,intvl_rad
     1             ,n_radarobs_lvl_fltrd
501     format(1x,i3,i6,i7,i9)

      enddo ! k

      write(6,*)' # Radar Obs Rejected due to other data = ',i_radar_rej
     1ect
      write(6,502)n_radarobs_tot_unfltrd,n_radarobs_tot_fltrd
502   format(1x,' # Radar Obs TOTAL UNFILTERED / FILTERED = ',2i7)


c  convert radar obs into u & v by using tangential component of first pass
      write(6,*)' Generating derived radar obs, opening d00 file, i4time
     1 = '
     1                  ,i4time
      write(6,*)'  i   j   k    df    vr    fgr'

      if(l_derived_output)then
          if(i_radar .le. 99)then
              write(ext,531)i_radar
 531          format('d',i2.2)
          else
              ext = 'dxx'
          endif
          call open_lapsprd_file(61,i4time,ext,istatus)
      endif

      height_grid = 0. ! This approximation won't hurt the azimuth

      do k=1,kmax
        icount_output = 0

        do j=1,jmax
        do i=1,imax
          if(vr_obs_fltrd(i,j,k) .ne. r_missing_data)then
            call latlon_to_radar(lat(i,j),lon(i,j),height_grid
     1              ,azimuth,slant_range,elev
     1                  ,rlat_radar,rlon_radar,rheight_radar)


            if(abs(upass1(i,j,k)) .ge. 1e6 .or.
     1         abs(vpass1(i,j,k)) .ge. 1e6)then
                ierr_count = ierr_count + 1
                if(ierr_count .lt. 100)write(6,*)
     1      ' Error in upass1,vpass1',i,j,k,upass1(i,j,k),vpass1(i,j,k)
                istatus = 0
                return

            else ! valid 1st pass winds
                if(l_grid_north)then

                    call uvgrid_to_radar(
     1                       upass1(i,j,k) + u_laps_bkg(i,j,k),
     1                       vpass1(i,j,k) + v_laps_bkg(i,j,k),
     1                       t_radar,
     1                       r_radar,
     1                       azimuth,
     1                       lat(i,j),
     1                       lon(i,j) )

                    call radar_to_uvgrid(t_radar,
     1                       vr_obs_fltrd(i,j,k),
     1                       u_wind,
     1                       v_wind,
     1                       azimuth,
     1                       lat(i,j),
     1                       lon(i,j) )

                else ! we are using true north winds

                    call uvtrue_to_radar(
     1                       upass1(i,j,k) + u_laps_bkg(i,j,k),
     1                       vpass1(i,j,k) + v_laps_bkg(i,j,k),
     1                       t_radar,
     1                       r_radar,
     1                       azimuth)

                    call radar_to_uvtrue(t_radar,
     1                       vr_obs_fltrd(i,j,k),
     1                       u_wind,
     1                       v_wind,
     1                       azimuth)

                endif ! l_grid_north

                call uv_to_disp(u_wind,
     1                          v_wind,
     1                          di_wind,
     1                          speed)


!               Compare radar radial velocity to 1st pass analysis
                diff_radial = vr_obs_fltrd(i,j,k) - r_radar

                icount_output = icount_output + 1
                if(icount_output .le. 3)write(6,310)i,j,k
     1                                 ,diff_radial
     1                                 ,vr_obs_fltrd(i,j,k),r_radar
310             format(3i4,3f6.1)

!               if(k .eq. 13 .and. i .eq. 29 .and. j .eq. 23)then
!                   write(6,*)' Magic Pt: t_radar',t_radar
!               endif


                if(.true.)then

!                 Set flag for good_multi_doppler_ob
!                 This checks for multi Doppler when mode = 2. Note that the
!                 geometry of the multi-Doppler measurements is not yet considered
                  if(mode .eq. 2)then
                      if(wt_p_spread(i,j,k) .eq. weight_radar)then
                          l_good_multi_doppler_ob(i,j,k) = .true.
                      endif
                  endif

!                 Subtract background from radar ob to get difference radar ob
                  uobs_diff_spread(i,j,k) = u_wind - u_laps_bkg(i,j,k)
                  vobs_diff_spread(i,j,k) = v_wind - v_laps_bkg(i,j,k)

!                 write(6,320)i,j,k,di_wind,speed,u_wind,v_wind
!320              format(1x,3i2,2f6.1,2f6.1)

                  if(l_derived_output)then
                      write(61,321)i-1,j-1,k-1,di_wind,speed
321                   format(1x,3i4,2f6.1,2f6.1)
                  endif

                  wt_p_spread(i,j,k) = weight_radar

                endif

            endif ! Error condition

          endif ! Not missing data
        enddo ! i
        enddo ! j
      enddo ! k

      if(l_derived_output)then
          close(61)
      endif

csms$ignore end
      return
      end



      subroutine filter_radar_obs(
     1                  n_radarobs_lvl_unfltrd,       ! Input
     1                  imax,jmax,                    ! Input
     1                  vr_obs_unfltrd,               ! Input
     1                  r_missing_data,               ! Input
     1                  vr_obs_fltrd,                 ! Input/Output
     1                  intvl_rad)                    ! Output

!       Filter to make the filtered ob array more sparse

        implicit none

!       Threshold number of radar obs on a given level
        integer*4 thresh_2_radarobs_lvl_unfltrd
     1           ,thresh_4_radarobs_lvl_unfltrd
        parameter (thresh_2_radarobs_lvl_unfltrd = 300)
        parameter (thresh_4_radarobs_lvl_unfltrd = 600)

        integer*4 n_radarobs_lvl_unfltrd, intvl_rad, imax, jmax
        real*4 vr_obs_unfltrd(imax,jmax)
        real*4 vr_obs_fltrd(imax,jmax)
        real*4 r_missing_data

        logical l_found_one, l_imax_odd, l_jmax_odd
        integer i,j,ii,jj

csms$ignore begin
        l_imax_odd = 1 .eq. mod(imax,2)
        l_jmax_odd = 1 .eq. mod(jmax,2)

!       Test against threshold for number of radar obs on a given level
!       This logic is supposed to select a sparse subset of the radial
!       velocities yet retain grid boxes that are isolated
        if(n_radarobs_lvl_unfltrd .gt. thresh_4_radarobs_lvl_unfltrd
     1                                                         )then
           ! Keep only every fourth ob.  Keep one ob out of every `quad'.
           intvl_rad = 4
           do j=1,jmax-1,2
           do i=1,imax-1,2
              l_found_one = .false.
              do jj = j,j+1
              do ii = i,i+1
                 if  ( l_found_one ) then
                    vr_obs_fltrd(ii,jj) = r_missing_data
                 elseif ( vr_obs_unfltrd(ii,jj) .ne. r_missing_data
     1                                                            )then
                    vr_obs_fltrd(ii,jj) = vr_obs_unfltrd(ii,jj)
                    l_found_one = .true.
                 else
                    vr_obs_fltrd(ii,jj) = r_missing_data
                 endif
              enddo ! ii
              enddo ! jj
           enddo ! i

           if ( l_imax_odd ) then
              vr_obs_fltrd(imax,j) = r_missing_data
              vr_obs_fltrd(imax,j+1) = r_missing_data
           endif

           enddo ! j

           if ( l_imax_odd ) then
              do i=1,imax
                 vr_obs_fltrd(i,jmax) = r_missing_data
              enddo
           endif

        elseif(n_radarobs_lvl_unfltrd .gt.
     1         thresh_2_radarobs_lvl_unfltrd)then
           ! Keep every other ob.  Keep one ob out of every `pair'.
           intvl_rad = 2
           do j=1,jmax
           do i=1,imax
             vr_obs_fltrd(i,j) = r_missing_data
             if((i+j) .eq. ((i+j)/2)*2 .and. i .lt. imax)then ! Sum is even
               l_found_one = .false.
               do ii = i,i+1
                 if (vr_obs_unfltrd(ii,j) .ne. r_missing_data) then
                   if(.not. l_found_one)then
                      vr_obs_fltrd(ii,j) = vr_obs_unfltrd(ii,j)
                      l_found_one = .true.
                   endif
                 endif
               enddo ! ii
             endif ! On checkerboard
           enddo ! i
           enddo ! j

        else
           ! Keep every ob
           intvl_rad = 1
           do j=1,jmax
           do i=1,imax
              vr_obs_fltrd(i,j) = vr_obs_unfltrd(i,j)
           enddo
           enddo
      endif

csms$ignore end
      return
      end


      subroutine spread_vert(uobs_in,vobs_in,l_3d,iwrite,uobs_out
     1                      ,vobs_out,wt_p,weights,pres_3d,i,j
     1                      ,imax,jmax,kmax,istatus)

!     Modified 7/94 S. Albers to allow better handling of variable
!     vertical resolution of pressure
!     3/97 S. Albers: Obtain parameters from runtime file

      include 'windparms.inc'

      real*4 uobs_in(imax,jmax,kmax)
      real*4 uobs_out(imax,jmax,kmax)
      real*4 vobs_in(imax,jmax,kmax)
      real*4 vobs_out(imax,jmax,kmax)
      real*4 wt_p(imax,jmax,kmax)                      ! Input
      real*4 weights(imax,jmax,kmax)                   ! Output (spread)
      real*4 pres_3d(imax,jmax,kmax)                   ! Input

      integer*4 vert_rad_pirep
      integer*4 vert_rad_sfc
      integer*4 vert_rad_cdw
      integer*4 vert_rad_prof

      logical l_3d

!     Vertical radius of influence for each data source (pascals)
      real*4 r0_vert_pirep,r0_vert_cdw,r0_vert_sfc,r0_vert_prof
      parameter (r0_vert_pirep = 2500.)
      parameter (r0_vert_cdw  = 2500.)
      parameter (r0_vert_sfc   = 2500.)
      parameter (r0_vert_prof  = 2500.)

csms$ignore begin
      call get_r_missing_data(r_missing_data, istatus)
      if(istatus .ne. 1)return

      vert_rad_pirep = 0
      vert_rad_sfc   = 0
      vert_rad_cdw   = 0
      vert_rad_prof  = 0

      pres_range_pirep = 0.
      pres_range_cdw   = 0.
      pres_range_sfc   = 0.
      pres_range_prof  = 0.

!     Initialize the obs_out columns
      do k = 1,kmax
          uobs_out(i,j,k) = uobs_in(i,j,k)
          vobs_out(i,j,k) = vobs_in(i,j,k)
          weights(i,j,k) = wt_p(i,j,k)
      enddo ! k

      do k = 1,kmax
          if(.true.)then ! l_3d
              kklow = k
              kkhigh = k
          else
              kklow = 1
              kkhigh = kmax
          endif

          if(weights(i,j,k) .eq. weight_pirep)then ! Spread this pirep vertically
              do kk = kklow,kkhigh
                  if(.true.)then ! l_3d
                      if(k .ne. kk)stop
                      dist_pa = 0.
                  else
                      dist_pa = abs(pres_3d(i,j,kk) - pres_3d(i,j,k))
                  endif

                  if(weights(i,j,kk) .eq. r_missing_data   .and.
     1               dist_pa         .le. pres_range_pirep       )then
                      weights(i,j,kk) = weight_pirep
     1                  * exp(-(dist_pa/r0_vert_pirep))
                      uobs_out(i,j,kk) = uobs_in(i,j,k)
                      vobs_out(i,j,kk) = vobs_in(i,j,k)
                      if(iwrite .eq. 1)write(6,101)i,j,k,kk
     1                                            ,uobs_out(i,j,kk)
     1                                            ,vobs_out(i,j,kk)
     1                                            ,weights(i,j,kk)
101                   format(' Prp',2i5,2i4,2f6.1,f8.5)
                  endif
              enddo
          endif

          if(weights(i,j,k) .eq. weight_cdw)then ! Spread this meso vertically
              do kk = kklow,kkhigh
                  if(.true.)then ! l_3d
                      if(k .ne. kk)stop
                      dist_pa = 0.
                  else
                      dist_pa = abs(pres_3d(i,j,kk) - pres_3d(i,j,k))
                  endif

                  if(weights(i,j,kk) .eq. r_missing_data   .and.
     1               dist_pa         .le. pres_range_cdw         )then
                      weights(i,j,kk) = weight_cdw
     1                          * exp(-(dist_pa/r0_vert_cdw))
                      uobs_out(i,j,kk) = uobs_in(i,j,k)
                      vobs_out(i,j,kk) = vobs_in(i,j,k)
                      if(iwrite .eq. 1)write(6,201)i,j,k,kk
     1                                            ,uobs_out(i,j,kk)
     1                                            ,vobs_out(i,j,kk)
     1                                            ,weights(i,j,kk)
201                   format(' Cdw',2i5,2i4,2f6.1,f8.5)
                  endif
              enddo
          endif

          if(weights(i,j,k) .eq. weight_sfc)then ! Spread this Sfc vertically
              do kk = kklow,kkhigh
                  if(.true.)then ! l_3d
                      if(k .ne. kk)stop
                      dist_pa = 0.
                  else
                      dist_pa = abs(pres_3d(i,j,kk) - pres_3d(i,j,k))
                  endif

                  if(weights(i,j,kk) .eq. r_missing_data   .and.
     1               dist_pa         .le. pres_range_sfc         )then
                      weights(i,j,kk) = weight_sfc
     1                          * exp(-(dist_pa/r0_vert_sfc))
                      uobs_out(i,j,kk) = uobs_in(i,j,k)
                      vobs_out(i,j,kk) = vobs_in(i,j,k)
                      if(iwrite .eq. 1)write(6,301)i,j,k,kk
     1                                            ,uobs_out(i,j,kk)
     1                                            ,vobs_out(i,j,kk)
     1                                            ,weights(i,j,kk)
301                   format(' Sfc',2i5,2i4,2f6.1,f8.5)
                  endif
              enddo
          endif

          if(weights(i,j,k) .eq. weight_prof)then ! Spread this profiler vertically

!             Spread on high side
              kp1 = min(kmax,k+1)
              if(weights(i,j,kp1) .eq. r_missing_data)then
              do kk = kklow,kkhigh
                  if(.true.)then ! l_3d
                      if(k .ne. kk)stop
                      dist_pa = 0.
                  else
                      dist_pa = abs(pres_3d(i,j,kk) - pres_3d(i,j,k))
                  endif

                  if(weights(i,j,kk) .eq. r_missing_data   .and.
     1               dist_pa         .le. pres_range_prof  .and.
     1               kk              .gt. k                      )then
                      weights(i,j,kk) = weight_prof
     1                          * exp(-(dist_pa/r0_vert_prof))
                      uobs_out(i,j,kk) = uobs_in(i,j,k)
                      vobs_out(i,j,kk) = vobs_in(i,j,k)
                      if(iwrite .eq. 1)write(6,401)i,j,k,kk
     1                                            ,uobs_out(i,j,kk)       
     1                                            ,vobs_out(i,j,kk)
     1                                            ,weights(i,j,kk)
401                   format(' Prf',2i5,2i4,2f6.1,f8.5)
                  endif
                enddo ! kk
              endif

!             Spread on low side
              km1 = max(1,k-1)
              if(weights(i,j,km1) .eq. r_missing_data)then
              do kk = kklow,kkhigh
                  if(.true.)then ! l_3d
                      if(k .ne. kk)stop
                      dist_pa = 0.
                  else
                      dist_pa = abs(pres_3d(i,j,kk) - pres_3d(i,j,k))
                  endif

                  if(weights(i,j,kk) .eq. r_missing_data   .and.
     1               dist_pa         .le. pres_range_prof  .and.
     1               kk              .lt. k                      )then
                      weights(i,j,kk) = weight_prof
     1                          * exp(-(dist_pa/r0_vert_prof))
                      uobs_out(i,j,kk) = uobs_in(i,j,k)
                      vobs_out(i,j,kk) = vobs_in(i,j,k)
                      if(iwrite .eq. 1)write(6,401)i,j,k,kk
     1                                            ,uobs_out(i,j,kk)
     1                                            ,vobs_out(i,j,kk)
     1                                            ,weights(i,j,kk)
                  endif
                enddo ! kk
              endif
          endif
      enddo ! k
csms$ignore end

      return
      end


      subroutine get_inst_err(imax,jmax,kmax,r_missing_data
     1           ,wt_p_spread,rms_thresh_norm,rms_inst,rms_thresh)

      real*4    wt_p_spread(imax,jmax,kmax)                        ! Input

csms$ignore begin

      write(6,*)
      write(6,*)' subroutine get_inst_err...'

      n_obs_total = 0
      wt_p_inv_total = 0.

      do i = 1,imax
      do j = 1,jmax
      do k = 1,kmax
          if(wt_p_spread(i,j,k) .ne. r_missing_data)then
              n_obs_total = n_obs_total + 1
              wt_p_inv_total = wt_p_inv_total + 1.0 / wt_p_spread(i,j,k)
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

