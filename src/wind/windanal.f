       subroutine laps_anl(uobs,vobs,n_radars,
     1      vr_obs_unfltrd,vr_nyq,v_nyquist_in
     1     ,upass1,vpass1                                        ! Output
     1     ,n_var,n_fnorm                                        ! Input
     1     ,uanl,vanl,varanl                                     ! Output
     1     ,wt_p,wt_p_spread,weight_bkg_const                    ! Input/Local
     1     ,max_radars
     1     ,n_radarobs_tot_unfltrd,rlat_radar,rlon_radar,rheight_radar
     1     ,u_laps_bkg,v_laps_bkg                                ! Input/Output
     1     ,imax,jmax,kmax,lat,lon
     1     ,i4time,grid_spacing_m
     1     ,r_missing_data
     1     ,l_good_multi_doppler_ob,l_analyze
     1     ,l_derived_output,l_grid_north,l_3pass,l_correct_unfolding
     1     ,n_iter_wind
     1     ,weight_meso,weight_sao,weight_pirep,weight_prof,weight_radar     
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

      integer*4 imax,jmax,kmax        ! 3D array dimensions        ! Input

!     3D arrays of u/v observations, all data sources, one datum per gridpoint.
      real*4 uobs(imax,jmax,kmax),vobs(imax,jmax,kmax)             ! Input

      integer*4 n_radars   ! Actual number of radars having data     Input

!     First pass analyzed winds
      real*4 upass1(imax,jmax,kmax), vpass1(imax,jmax,kmax)        ! Output
      real*4 varpass1(imax,jmax,kmax,n_var)                        ! Equiv Abv

!     Final pass analyzed winds
      real*4 uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)             ! Output
      real*4 varanl(imax,jmax,kmax,n_var)                          ! Equiv Abv

!     3D array of observation weights, depends on data type
!     The choices are outlined below
      real*4    wt_p(imax,jmax,kmax)                               ! Input
      real*4    wt_p_spread(imax,jmax,kmax)                        ! Local

!     Model background field
      real*4 u_laps_bkg(imax,jmax,kmax),v_laps_bkg(imax,jmax,kmax) ! Input

!     Arrays of lat and lon for each gridpoint
      real*4 lat(imax,jmax),lon(imax,jmax)                         ! Input

      integer*4 i4time                                             ! Input

      real*4 r_missing_data                  ! missing data value    Input

!-----Radar Data -----------------------------------------------------------

      integer*4 max_radars            ! max possible # of radars         Input

!     4D Radial velocity array (all radars)
      real*4 vr_obs_unfltrd(imax,jmax,kmax,max_radars)                 ! Input
      real*4 vr_nyq(imax,jmax,kmax,max_radars)                         ! Input

!     Nyquist velocity (if known and constant) for each radar
      real*4 v_nyquist_in(max_radars)                                  ! Input

!     3D Radial velocity for a given radar after filtering
      real*4 vr_obs_fltrd(imax,jmax,kmax)                              ! Local

!     Location of each radar
      real*4 rlat_radar(max_radars),rlon_radar(max_radars)             ! Input
     1                     ,rheight_radar(max_radars)

!     # of radar obs before filtering for each radar (modified by QC)
      integer*4 n_radarobs_tot_unfltrd(max_radars)                     ! Input/Modified

!--------------------------------------------------------------------------------

      real*4    r0_array_out(imax,jmax)                                ! Local
      real*4    density_array_in(imax,jmax)                            ! Local

      real*4 uobs_diff(imax,jmax,kmax),vobs_diff(imax,jmax,kmax)       ! Local
      real*4 uobs_diff_spread(imax,jmax,kmax)                          ! Local
      real*4 vobs_diff_spread(imax,jmax,kmax)                          ! Local
      real*4 varobs_diff_spread(imax,jmax,kmax,n_var)                  ! Local

      integer*4 n_obs_lvl(kmax)                                        ! Local
      real*4 upass1_buf(imax,jmax,kmax)                                ! Local
      real*4 vpass1_buf(imax,jmax,kmax)                                ! Local
      logical l_good_multi_doppler_ob(imax,jmax,kmax)                  ! Local
      logical  l_analyze(kmax) ! This depends on presence of radar obs ! Local
      logical  l_derived_output ! Flag for producing derived output    ! Input
      logical  l_grid_north     ! Flag for grid north or true north    ! Input
      logical  l_3pass          ! Flag for doing 3 pass analysis       ! Input
      logical  l_correct_unfolding ! Flag for dealiasing               ! Input

!     These are the weights of the various data types (filling the 3D array)
      real*4 weight_meso,weight_sao,weight_pirep,weight_prof,weight_rada
     1r ! Input

      integer*4 istatus         ! (1 is good)                          ! Output

!****************END ARGUMENT LIST *******************************************

      integer*4  n_fnorm

      dimension fnorm(0:n_fnorm)

      do iter = 1,n_iter_wind

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
91    format(1x,' Subtracting the background from the obs, then spreadin
     1g the'
     1  ,' obs vertically.'
     1  /'      i   j   k  kk  udf  vdf       uob  vob     ubg   vbg  sp
     1  wt')

      qc_thresh = 30. ! Threshold speed for throwing out the ob

      do j=1,jmax
      do i=1,imax
        do k = 1,kmax
!         wt_p_spread(i,j,k) = wt_p(i,j,k)                         ! Initialize the local array

          if(wt_p(i,j,k) .ne. r_missing_data)then
            if(uobs(i,j,k) .ne. r_missing_data          .and.
     1       vobs(i,j,k) .ne. r_missing_data                        )the
     1n

              speed_bkg  = sqrt(u_laps_bkg(i,j,k)**2 + v_laps_bkg(i,j,k)
     1**2)

              uobs_diff(i,j,k) = uobs(i,j,k) - u_laps_bkg(i,j,k)
              vobs_diff(i,j,k) = vobs(i,j,k) - v_laps_bkg(i,j,k)
              speed_diff = sqrt(uobs_diff(i,j,k)**2  + vobs_diff(i,j,k)*
     1*2)

!             Apply QC check to the OB against the background analysis
              if(
!                    Make sure we actually have a real reference background
     1    (speed_bkg .gt. 0.) .and.

!                    General QC check
     1    (speed_diff .gt. qc_thresh

!                    Stricter QC check for pireps
     1      .or. (speed_diff .gt. 10. .and. wt_p(i,j,k) .eq. weight_pire
     1p)

!                    Stricter QC check for profilers
     1      .or. (speed_diff .gt. 15. .and. wt_p(i,j,k) .eq. weight_prof
     1)
     1                                                               )

!       1       .OR. (abs(k-16) .le. 5 .and. (i .ne. 50 .or. j .ne. 14)) )
     1                                                          )then

                ! Throw out the ob
                  if(wt_p(i,j,k) .eq. weight_pirep)then
                      write(6,101,err=199)i,j,k
     1                  ,uobs_diff(i,j,k)
     1                  ,vobs_diff(i,j,k)
     1                  ,uobs(i,j,k)
     1                  ,vobs(i,j,k)
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,wt_p(i,j,k)
101                   format(' Prp QCed out - ',2i4,i3,3(2x,2f5.0),f5.0,
     1f5.2)

                  elseif(wt_p(i,j,k) .eq. weight_meso)then
                      write(6,111,err=199)i,j,k
     1                  ,uobs_diff(i,j,k)
     1                  ,vobs_diff(i,j,k)
     1                  ,uobs(i,j,k)
     1                  ,vobs(i,j,k)
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,wt_p(i,j,k)
111                  format(' Mso QCed out - ',2i4,i3,3(2x,2f5.0),f5.0,f
     15.2)

                  elseif(wt_p(i,j,k) .eq. weight_sao)then
                      write(6,121,err=199)i,j,k
     1                  ,uobs_diff(i,j,k)
     1                  ,vobs_diff(i,j,k)
     1                  ,uobs(i,j,k)
     1                  ,vobs(i,j,k)
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,wt_p(i,j,k)
121                  format(' Sao QCed out - ',2i4,i3,3(2x,2f5.0),f5.0,f
     15.2)

                  elseif(wt_p(i,j,k) .eq. weight_prof)then
                      write(6,131,err=199)i,j,k
     1                  ,uobs_diff(i,j,k)
     1                  ,vobs_diff(i,j,k)
     1                  ,uobs(i,j,k)
     1                  ,vobs(i,j,k)
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,wt_p(i,j,k)
131                  format(' Prf QCed out - ',2i4,i3,3(2x,2f5.0),f5.0,f
     15.2)

                  endif ! Type of OB to write out

!                 Set the difference OB to missing (original ob left alone)
199               uobs_diff(i,j,k) = r_missing_data
                  vobs_diff(i,j,k) = r_missing_data
                  wt_p(i,j,k) = r_missing_data

              else ! write out obs
                  write(6,201,err=202)i,j,k
     1                  ,uobs_diff(i,j,k)
     1                  ,vobs_diff(i,j,k)
     1                  ,uobs(i,j,k)
     1                  ,vobs(i,j,k)
     1                  ,u_laps_bkg(i,j,k)
     1                  ,v_laps_bkg(i,j,k)
     1                  ,speed_diff
     1                  ,wt_p(i,j,k)
201               format(' OB ',2i4,i3,3x,f6.1,f6.1,2(2x,2f6.1),f5.1,f5.
     12)
202               continue

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


          else ! Set these to missing just in case
            uobs_diff(i,j,k) = r_missing_data
            vobs_diff(i,j,k) = r_missing_data

          endif ! wt_p .ne. MISSING

        enddo ! k

!       Spread the difference ob vertically
        call spread_vert(uobs_diff,vobs_diff,uobs_diff_spread,vobs_diff_
     1spread
     1          ,wt_p,wt_p_spread,i,j,imax,jmax,kmax,istatus)

        if(istatus .ne. 1)return

      enddo ! i
      enddo ! j

!     Initialize fnorm array used in barnes_multivariate
      write(6,*)' Creating fnorm LUT'
      r0_norm = 10. ! When the distance = r0_norm, the fnorm is effectively 1
      r0_norm_sq = r0_norm**2
      exp_offset = 70.
      expm80 = exp(-80.)
      do iii = 0,n_fnorm
          dist_norm_sq = (float(iii)/r0_norm_sq)
          arg = dist_norm_sq - exp_offset
          if(arg .le. 80.)then
              fnorm(iii) = exp(-arg)
          else
              fnorm(iii) = expm80
          endif
      enddo

      I4_elapsed = ishow_timer()

!     Perform 1st pass analysis of non-radar difference obs
      do k = 1,kmax
          l_analyze(k) = .true.
      enddo ! k

      call move_3d(upass1,varpass1(1,1,1,1),imax,jmax,kmax)
      call move_3d(vpass1,varpass1(1,1,1,2),imax,jmax,kmax)
      call move_3d(uobs_diff_spread,varobs_diff_spread(1,1,1,1)
     1                             ,imax,jmax,kmax)
      call move_3d(vobs_diff_spread,varobs_diff_spread(1,1,1,2)
     1                             ,imax,jmax,kmax)

      call barnes_multivariate(varpass1,n_var
     1        ,imax,jmax,kmax,grid_spacing_m
     1        ,varobs_diff_spread
     1        ,wt_p_spread,fnorm,n_fnorm
     1        ,l_analyze,weight_bkg_const,r0_array_out
     1        ,n_obs_lvl,istatus)
      if(istatus .ne. 1)return

      call move_3d(varpass1(1,1,1,1),upass1,imax,jmax,kmax)
      call move_3d(varpass1(1,1,1,2),vpass1,imax,jmax,kmax)
      call move_3d(varobs_diff_spread(1,1,1,1),uobs_diff_spread
     1                             ,imax,jmax,kmax)
      call move_3d(varobs_diff_spread(1,1,1,2),vobs_diff_spread
     1                             ,imax,jmax,kmax)

      I4_elapsed = ishow_timer()

      do k = 1,kmax
          if(n_obs_lvl(k) .eq. 0)then
              write(6,311)k
311           format(1x,' No obs at lvl',i3,' Using Zero array as 1st Pa
     1ss')
              do j=1,jmax
              do i=1,imax
                  upass1(i,j,k) = 0.
                  vpass1(i,j,k) = 0.
              enddo ! i
              enddo ! j
          endif

      enddo ! k

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
!         uobs_diff_spread and vobs_diff_spread
          call insert_derived_radar_obs(
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
     1  ,l_analyze                                  ! Output
     1  ,n_radarobs_tot_unfltrd                     ! Input
     1  ,vr_obs_fltrd                               ! Local
     1  ,upass1_buf,vpass1_buf                      ! Local
     1  ,l_good_multi_doppler_ob                    ! Local
     1  ,istatus                                    ! Input/Output
     1                                                          )

          I4_elapsed = ishow_timer()

          write(6,*)' Calling barnes with modified radar obs added'

          call move_3d(uanl,varanl(1,1,1,1),imax,jmax,kmax)
          call move_3d(vanl,varanl(1,1,1,2),imax,jmax,kmax)

          call move_3d(uobs_diff_spread,varobs_diff_spread(1,1,1,1)
     1                             ,imax,jmax,kmax)
          call move_3d(vobs_diff_spread,varobs_diff_spread(1,1,1,2)
     1                             ,imax,jmax,kmax)

          call barnes_multivariate(varanl,n_var
     1       ,imax,jmax,kmax,grid_spacing_m
     1       ,varobs_diff_spread
     1       ,wt_p_spread,fnorm,n_fnorm
     1       ,l_analyze,weight_bkg_const,r0_array_out
     1       ,n_obs_lvl,istatus)

          call move_3d(varanl(1,1,1,1),uanl,imax,jmax,kmax)
          call move_3d(varanl(1,1,1,2),vanl,imax,jmax,kmax)
          call move_3d(varobs_diff_spread(1,1,1,1),uobs_diff_spread
     1                             ,imax,jmax,kmax)
          call move_3d(varobs_diff_spread(1,1,1,2),vobs_diff_spread
     1                             ,imax,jmax,kmax)

          if(istatus .ne. 1)return

          I4_elapsed = ishow_timer()

      else ! n_radars .ne. 1

          mode = 2 ! Only multi-Doppler obs

!         Take the data from all the radars and add the derived radar obs into
!         uobs_diff_spread and vobs_diff_spread
          call insert_derived_radar_obs(
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
     1  ,l_analyze                                  ! Output
     1  ,n_radarobs_tot_unfltrd                     ! Input
     1  ,vr_obs_fltrd                               ! Local
     1  ,upass1_buf,vpass1_buf                      ! Local
     1  ,l_good_multi_doppler_ob                    ! Local
     1  ,istatus                                    ! Input/Output
     1                                                          )

          I4_elapsed = ishow_timer()

          write(6,*)' Calling barnes with only multi-doppler obs '
     1          ,'creating an intermediate analysis'

          call move_3d(uanl,varanl(1,1,1,1),imax,jmax,kmax)
          call move_3d(vanl,varanl(1,1,1,2),imax,jmax,kmax)

          call move_3d(uobs_diff_spread,varobs_diff_spread(1,1,1,1)
     1                             ,imax,jmax,kmax)
          call move_3d(vobs_diff_spread,varobs_diff_spread(1,1,1,2)
     1                             ,imax,jmax,kmax)

          call barnes_multivariate(varanl,n_var
     1       ,imax,jmax,kmax,grid_spacing_m
     1       ,varobs_diff_spread
     1       ,wt_p_spread,fnorm,n_fnorm
     1       ,l_analyze,weight_bkg_const,r0_array_out
     1       ,n_obs_lvl,istatus)

          call move_3d(varanl(1,1,1,1),uanl,imax,jmax,kmax)
          call move_3d(varanl(1,1,1,2),vanl,imax,jmax,kmax)
          call move_3d(varobs_diff_spread(1,1,1,1),uobs_diff_spread
     1                             ,imax,jmax,kmax)
          call move_3d(varobs_diff_spread(1,1,1,2),vobs_diff_spread
     1                             ,imax,jmax,kmax)

          if(istatus .ne. 1)return

!         Make sure each level of uanl and vanl is initialized in the event it
!         was not analyzed.
          do k = 1,kmax
              if(n_obs_lvl(k) .eq. 0)then
                  write(6,411)k
411               format(1x,' No obs at lvl',i3,
     1                  ' Insert Zero array into uanl,vanl')
                  do j=1,jmax
                  do i=1,imax
                      uanl(i,j,k) = 0.
                      vanl(i,j,k) = 0.
                  enddo ! i
                  enddo ! j
              elseif(.not. l_analyze(k))then
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
     1   mode                                       ! Input
     1  ,n_radars,max_radars                        ! Input
     1  ,imax,jmax,kmax                             ! Input
     1  ,r_missing_data                             ! Input
     1  ,vr_obs_unfltrd                             ! Input
     1  ,i4time                                     ! Input
     1  ,lat,lon                                    ! Input
     1  ,rlat_radar,rlon_radar                      ! Input
     1  ,rheight_radar                              ! Input
     1  ,uanl,vanl                                  ! Input
     1  ,u_laps_bkg,v_laps_bkg                      ! Input
     1  ,weight_radar                               ! Input
     1  ,l_derived_output,l_grid_north              ! Input
     1  ,wt_p_spread                                ! Input/Output
     1  ,uobs_diff_spread,vobs_diff_spread          ! Input/Output
     1  ,l_analyze                                  ! Output
     1  ,n_radarobs_tot_unfltrd                     ! Input
     1  ,vr_obs_fltrd                               ! Local
     1  ,upass1_buf,vpass1_buf                      ! Local
     1  ,l_good_multi_doppler_ob                    ! Local
     1  ,istatus                                    ! Input/Output
     1                                                          )

          I4_elapsed = ishow_timer()

          write(6,*)' Calling barnes with modified radar obs added'

          call move_3d(uanl,varanl(1,1,1,1),imax,jmax,kmax)
          call move_3d(vanl,varanl(1,1,1,2),imax,jmax,kmax)

          call move_3d(uobs_diff_spread,varobs_diff_spread(1,1,1,1)
     1                             ,imax,jmax,kmax)
          call move_3d(vobs_diff_spread,varobs_diff_spread(1,1,1,2)
     1                             ,imax,jmax,kmax)

          call barnes_multivariate(varanl,n_var,imax,jmax,kmax
     1       ,grid_spacing_m
     1       ,varobs_diff_spread
     1       ,wt_p_spread,fnorm,n_fnorm
     1       ,l_analyze,weight_bkg_const,r0_array_out
     1       ,n_obs_lvl,istatus)

          call move_3d(varanl(1,1,1,1),uanl,imax,jmax,kmax)
          call move_3d(varanl(1,1,1,2),vanl,imax,jmax,kmax)
          call move_3d(varobs_diff_spread(1,1,1,1),uobs_diff_spread
     1                             ,imax,jmax,kmax)
          call move_3d(varobs_diff_spread(1,1,1,2),vobs_diff_spread
     1                             ,imax,jmax,kmax)

          if(istatus .ne. 1)return

            I4_elapsed = ishow_timer()


      endif ! n_radars

      write(6,*)' Adding analyzed differences to the background to recon
     1struct'
     1  ,' full analyses'

      do k=1,kmax ! Add back differences for first pass

          do j=1,jmax
          do i=1,imax
              if(upass1(i,j,k) .ne. r_missing_data)then
                  upass1(i,j,k) = upass1(i,j,k) + u_laps_bkg(i,j,k)
                  vpass1(i,j,k) = vpass1(i,j,k) + v_laps_bkg(i,j,k)
              else
                  write(6,*)' ERROR: Missing data value(s) detected in f
     1irst'
     1      ,' pass at lvl',k
                  istatus = 0
                  return
              endif
          enddo ! i
          enddo ! j


          if(l_analyze(k))then ! This depends on the presence of radar obs
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

      enddo ! n_iter_wind

      istatus = 1
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
     1  ,l_analyze                                  ! Output
     1  ,n_radarobs_tot_unfltrd                     ! Input
     1  ,vr_obs_fltrd                               ! Local
     1  ,upass1_buf,vpass1_buf                      ! Local
     1  ,l_good_multi_doppler_ob                    ! Local
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
      real*4   vr_obs_fltrd(imax,jmax,kmax)

      real*4 upass1_buf(imax,jmax,kmax)
      real*4 vpass1_buf(imax,jmax,kmax)

      logical  l_good_multi_doppler_ob(imax,jmax,kmax)
      logical  l_analyze(kmax),l_derived_output,l_grid_north

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

      enddo ! k

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

      character*13 filename13
      character*31 ext

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

        n_radarobs_tot_fltrd = n_radarobs_tot_fltrd + n_radarobs_lvl_flt
     1rd

        write(6,501)k,n_radarobs_lvl_unfltrd,intvl_rad,n_radarobs_lvl_fl
     1trd
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
          ext = 'd20'
          if(i_radar .eq. 1)ext = 'd01'
          if(i_radar .eq. 2)ext = 'd02'
          if(i_radar .eq. 3)ext = 'd03'
          if(i_radar .eq. 4)ext = 'd04'
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


            if(abs(upass1(i,j,k)) .ge. 1e6 .or. abs(vpass1(i,j,k)) .ge. 
     11e6)then
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
     1                       lon(i,j) )

                    call radar_to_uvgrid(t_radar,
     1                       vr_obs_fltrd(i,j,k),
     1                       u_wind,
     1                       v_wind,
     1                       azimuth,
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
                if(icount_output .le. 3)
     1    write(6,310)i,j,k,diff_radial,vr_obs_fltrd(i,j,k),r_radar
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
                 elseif ( vr_obs_unfltrd(ii,jj) .ne. r_missing_data ) th
     1en
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

      return
      end



      subroutine spread_vert(uobs_in,vobs_in,uobs_out,vobs_out,
     1      wt_p,weights,i,j,imax,jmax,kmax,istatus)

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

      real*4 PRESSURE_INTERVAL_L

      integer*4 vert_rad_pirep
      integer*4 vert_rad_sao
      integer*4 vert_rad_meso
      integer*4 vert_rad_prof

!     Vertical radius of influence for each data source (pascals)
      real*4 r0_vert_pirep,r0_vert_meso,r0_vert_sao,r0_vert_prof
      parameter (r0_vert_pirep = 2500.)
      parameter (r0_vert_meso  = 2500.)
      parameter (r0_vert_sao   = 2500.)
      parameter (r0_vert_prof  = 2500.)

      integer*4 vert_rad_pirep_s,vert_rad_meso_s,vert_rad_sao_s
     1         ,vert_rad_prof_s

      call get_r_missing_data(r_missing_data, istatus)
      if(istatus .ne. 1)return

      call get_vert_rads       (vert_rad_pirep,
     1                          vert_rad_sao,
     1                          vert_rad_meso,
     1                          vert_rad_prof,
     1                          istatus)
      if(istatus .ne. 1)return

      call get_pressure_interval(PRESSURE_INTERVAL_L,istatus)
      if(istatus .ne. 1)return

!     Initialize the obs_out columns
      do k = 1,kmax
          uobs_out(i,j,k) = uobs_in(i,j,k)
          vobs_out(i,j,k) = vobs_in(i,j,k)
          weights(i,j,k) = wt_p(i,j,k)
      enddo ! k

      iscale = nint(5000. / PRESSURE_INTERVAL_L) ! Affects # of grid points looped in the vertical
      vert_rad_pirep_s = vert_rad_pirep * iscale
      vert_rad_meso_s  = vert_rad_meso  * iscale
      vert_rad_sao_s   = vert_rad_sao   * iscale
      vert_rad_prof_s  = vert_rad_prof  * iscale

      do k = 1,kmax
          if(weights(i,j,k) .eq. weight_pirep)then ! Spread this pirep vertically
              do kk = max(1,k-vert_rad_pirep_s)
     1               ,min(kmax,k+vert_rad_pirep_s)
                  dist_pa = abs(k - kk) * PRESSURE_INTERVAL_L
                  if(weights(i,j,kk) .eq. r_missing_data)then
                      weights(i,j,kk) = weight_pirep
     1                  * exp(-(dist_pa/r0_vert_pirep))
                      uobs_out(i,j,kk) = uobs_in(i,j,k)
                      vobs_out(i,j,kk) = vobs_in(i,j,k)
                      write(6,101)i,j,k,kk,uobs_out(i,j,kk),vobs_out(i,j
     1,kk)
     1                                                  ,weights(i,j,kk)
101                   format(' Prp',2i4,2i3,2f6.1,f8.5)
                  endif
              enddo
          endif

          if(weights(i,j,k) .eq. weight_meso)then ! Spread this meso vertically
              do kk = max(1,k-vert_rad_meso_s)
     1               ,min(kmax,k+vert_rad_meso_s)
                  dist_pa = abs(k - kk) * PRESSURE_INTERVAL_L
                  if(weights(i,j,kk) .eq. r_missing_data)then
                      weights(i,j,kk) = weight_meso
     1                          * exp(-(dist_pa/r0_vert_meso))
                      uobs_out(i,j,kk) = uobs_in(i,j,k)
                      vobs_out(i,j,kk) = vobs_in(i,j,k)
                      write(6,201)i,j,k,kk,uobs_out(i,j,kk),vobs_out(i,j
     1,kk)
     1                                                  ,weights(i,j,kk)
201                   format(' Mso',2i4,2i3,2f6.1,f8.5)
                  endif
              enddo
          endif

          if(weights(i,j,k) .eq. weight_sao)then ! Spread this sao vertically
              do kk = max(1,k-vert_rad_sao_s)
     1               ,min(kmax,k+vert_rad_sao_s)
                  dist_pa = abs(k - kk) * PRESSURE_INTERVAL_L
                  if(weights(i,j,kk) .eq. r_missing_data)then
                      weights(i,j,kk) = weight_sao
     1                          * exp(-(dist_pa/r0_vert_sao))
                      uobs_out(i,j,kk) = uobs_in(i,j,k)
                      vobs_out(i,j,kk) = vobs_in(i,j,k)
                      write(6,301)i,j,k,kk,uobs_out(i,j,kk),vobs_out(i,j
     1,kk)
     1                                                  ,weights(i,j,kk)
301                   format(' Sao',2i4,2i3,2f6.1,f8.5)
                  endif
              enddo
          endif

          if(weights(i,j,k) .eq. weight_prof)then ! Spread this profiler vertically

!             Spread on high side
              kp1 = min(kmax,k+1)
              if(weights(i,j,kp1) .eq. r_missing_data)then

                do kk = kp1,min(kmax,k+vert_rad_prof_s)
                  dist_pa = abs(k - kk) * PRESSURE_INTERVAL_L
                  if(weights(i,j,kk) .eq. r_missing_data)then
                      weights(i,j,kk) = weight_prof
     1                          * exp(-(dist_pa/r0_vert_prof))
                      uobs_out(i,j,kk) = uobs_in(i,j,k)
                      vobs_out(i,j,kk) = vobs_in(i,j,k)
                      write(6,401)i,j,k,kk,uobs_out(i,j,kk),vobs_out(i,j
     1,kk)
     1                                                  ,weights(i,j,kk)
401                   format(' Prf',2i4,2i3,2f6.1,f8.5)
                  endif
                enddo ! kk

              endif

!             Spread on low side
              km1 = max(1,k-1)
              if(weights(i,j,km1) .eq. r_missing_data)then

                do kk = km1,max(1,k-vert_rad_prof_s),-1
                  dist_pa = abs(k - kk) * PRESSURE_INTERVAL_L
                  if(weights(i,j,kk) .eq. r_missing_data)then
                      weights(i,j,kk) = weight_prof
     1                          * exp(-(dist_pa/r0_vert_prof))
                      uobs_out(i,j,kk) = uobs_in(i,j,k)
                      vobs_out(i,j,kk) = vobs_in(i,j,k)
                      write(6,401)i,j,k,kk,uobs_out(i,j,kk),vobs_out(i,j
     1,kk)
     1                                                  ,weights(i,j,kk)
                  endif
                enddo ! kk

              endif
          endif
      enddo ! k

      return
      end

      subroutine get_pressure_interval(pressure_interval,istatus)

      include 'lapsparms.cmn' ! PRESSURE_INTERVAL_L

!     This routine accesses the pressure_interval variable from the
!     .parms file via the common block. Note the variable name in the
!     argument list is different in the calling routine

      if(iflag_lapsparms_cmn .ne. 1)then
          write(6,*)' ERROR, get_laps_config not called'
          istatus = 0
          return
!         stop
      endif

      pressure_interval = PRESSURE_INTERVAL_L

      return
      end

