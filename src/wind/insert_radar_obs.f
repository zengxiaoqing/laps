
      subroutine insert_derived_radar_obs(
     1   mode                                       ! Input
     1  ,n_radars,max_radars,idx_radar_a            ! Input
     1  ,imax,jmax,kmax                             ! Input
     1  ,r_missing_data                             ! Input
     1  ,vr_obs_unfltrd                             ! Input
     1  ,thresh_2_radarobs_lvl_unfltrd              ! Input
     1  ,thresh_4_radarobs_lvl_unfltrd              ! Input
     1  ,i4time                                     ! Input
     1  ,lat,lon                                    ! Input
     1  ,rlat_radar,rlon_radar                      ! Input
     1  ,rheight_radar                              ! Input
     1  ,upass1,vpass1                              ! Input
     1  ,u_laps_bkg,v_laps_bkg                      ! Input
     1  ,weight_radar                               ! Input
     1  ,l_derived_output,l_grid_north              ! Input
     1  ,wt_p_radar                                 ! Input/Output
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
      real*4   wt_p_radar(imax,jmax,kmax)

      real*4   vr_obs_fltrd(imax,jmax,kmax)                          ! Local
      real*4   upass1_buf(imax,jmax,kmax)                            ! Local
      real*4   vpass1_buf(imax,jmax,kmax)                            ! Local

      integer*4 idx_radar_a(max_radars)
      integer*4 thresh_2_radarobs_lvl_unfltrd
     1         ,thresh_4_radarobs_lvl_unfltrd

      logical  l_good_multi_doppler_ob(imax,jmax,kmax)               ! Local
      logical  l_analyze(kmax),l_derived_output,l_grid_north
      logical  l_multi_doppler_new

      data l_multi_doppler_new /.false./ ! Flag for new CWB routine

csms$ignore begin
      write(6,*)' Entering insert_derived_radar_obs, mode =',mode

!     Initialize l_good_multi_doppler_ob if mode = 2
      if(mode .eq. 2)then
          l_good_multi_doppler_ob = .false.
      endif

!     This routine takes the data from all the radars and adds the derived
!     radar obs into uobs_diff_spread and vobs_diff_spread

      if(.not. l_multi_doppler_new)then

!       Loop through the radars
        do i_radar = 1,n_radars

          write(6,*)' Generate derived radar "vector" obs, radar # '
     1             ,i_radar

          if(i_radar .eq. 1)then ! Use original 1st pass arrays
              call make_derived_radar_obs(
     1   imax,jmax,kmax                             ! Input
     1  ,mode                                       ! Input
     1  ,r_missing_data                             ! Input
     1  ,i_radar,idx_radar_a(i_radar)               ! Input
     1  ,vr_obs_unfltrd(1,1,1,i_radar)              ! Input
     1  ,thresh_2_radarobs_lvl_unfltrd              ! Input
     1  ,thresh_4_radarobs_lvl_unfltrd              ! Input
     1  ,i4time                                     ! Input
     1  ,lat,lon                                    ! Input
     1  ,rlat_radar(i_radar),rlon_radar(i_radar)    ! Input
     1  ,rheight_radar(i_radar)                     ! Input
     1  ,upass1,vpass1  ! 1st pass anal             ! Input
     1  ,u_laps_bkg,v_laps_bkg                      ! Input
     1  ,weight_radar                               ! Input
     1  ,l_derived_output,l_grid_north              ! Input
     1  ,wt_p_radar                                 ! Input/Output
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
     1  ,i_radar,idx_radar_a(i_radar)               ! Input
     1  ,vr_obs_unfltrd(1,1,1,i_radar)              ! Input
     1  ,thresh_2_radarobs_lvl_unfltrd              ! Input
     1  ,thresh_4_radarobs_lvl_unfltrd              ! Input
     1  ,i4time                                     ! Input
     1  ,lat,lon                                    ! Input
     1  ,rlat_radar(i_radar),rlon_radar(i_radar)    ! Input
     1  ,rheight_radar(i_radar)                     ! Input
     1  ,upass1_buf,vpass1_buf ! 1st pass anal + derived radar ! Input
     1  ,u_laps_bkg,v_laps_bkg                      ! Input
     1  ,weight_radar                               ! Input
     1  ,l_derived_output,l_grid_north              ! Input
     1  ,wt_p_radar                                 ! Input/Output
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
                  if(wt_p_radar(i,j,k) .eq. weight_radar)then
                      upass1_buf(i,j,k) = uobs_diff_spread(i,j,k)
                      vpass1_buf(i,j,k) = vobs_diff_spread(i,j,k)
                  endif
              enddo ! i
              enddo ! j
              enddo ! k
          endif

        enddo ! i_radar

      endif ! l_multi_doppler_new

      icount_radar_total = 0

!     Use only multiple Doppler obs if mode = 2, use all obs if mode = 1
      do k = 1,kmax
          icount_good_lvl = 0
          if(mode .eq. 1)then ! Use all Doppler obs
              do j = 1,jmax
              do i = 1,imax
                  if(wt_p_radar(i,j,k) .eq. weight_radar)then ! Good single or multi ob
                      icount_good_lvl = icount_good_lvl + 1
                  endif
              enddo ! i
              enddo ! j
          elseif(mode .eq. 2)then ! Throw out single Doppler obs
              do j = 1,jmax
              do i = 1,imax
                  if(wt_p_radar(i,j,k) .eq. weight_radar)then ! Good single or multi ob
                      if(.not. l_good_multi_doppler_ob(i,j,k))then
                          uobs_diff_spread(i,j,k) = r_missing_data
                          vobs_diff_spread(i,j,k) = r_missing_data
                          wt_p_radar(i,j,k) = r_missing_data
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
     1  ,i_radar,idx_radar                          ! Input
     1  ,vr_obs_unfltrd                             ! Input
     1  ,thresh_2_radarobs_lvl_unfltrd              ! Input
     1  ,thresh_4_radarobs_lvl_unfltrd              ! Input
     1  ,i4time                                     ! Input
     1  ,lat,lon                                    ! Input
     1  ,rlat_radar,rlon_radar,rheight_radar        ! Input
     1  ,upass1,vpass1                              ! Input
     1  ,u_laps_bkg,v_laps_bkg                      ! Input
     1  ,weight_radar                               ! Input
     1  ,l_derived_output,l_grid_north              ! Input
     1  ,wt_p_radar                                 ! Input/Output
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
      real*4   wt_p_radar(imax,jmax,kmax)
      real*4   vr_obs_fltrd(imax,jmax,kmax)

      logical  l_good_multi_doppler_ob(imax,jmax,kmax),l_derived_output
      logical  l_grid_north, l_multi_doppler_new

      integer*4 thresh_2_radarobs_lvl_unfltrd
     1         ,thresh_4_radarobs_lvl_unfltrd

      character*31 ext

      data l_multi_doppler_new /.false./

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
                if(wt_p_radar(i,j,k) .ne. weight_radar .and.
     1     wt_p_radar(i,j,k) .ne. r_missing_data)then ! Non-radar ob
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
     1                  thresh_2_radarobs_lvl_unfltrd,! Input
     1                  thresh_4_radarobs_lvl_unfltrd,! Input
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
      write(6,*)
     1   ' Generating derived radar obs, opening dxx file, i4time = '       
     1                  ,i4time
      write(6,*)'  i   j   k    df    vr    fgr   vt'

      if(l_derived_output)then
          if(i_radar .le. 99)then
              write(ext,531)idx_radar ! i_radar
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

                if(icount_output .eq. (icount_output/50)*50)then
                    write(6,310)i,j,k
     1                         ,diff_radial
     1                         ,vr_obs_fltrd(i,j,k),r_radar
     1                         ,t_radar
                endif
                icount_output = icount_output + 1
310             format(3i4,4f6.1)

!               if(k .eq. 13 .and. i .eq. 29 .and. j .eq. 23)then
!                   write(6,*)' Magic Pt: t_radar',t_radar
!               endif


                if(.true.)then

!                 Set flag for good_multi_doppler_ob
!                 This checks for multi Doppler when mode = 2. Note that the
!                 geometry of the multi-Doppler measurements is not yet considered
                  if(mode .eq. 2)then
                      if(wt_p_radar(i,j,k) .eq. weight_radar)then
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

                  wt_p_radar(i,j,k) = weight_radar

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
     1                  thresh_2_radarobs_lvl_unfltrd,! Input
     1                  thresh_4_radarobs_lvl_unfltrd,! Input
     1                  r_missing_data,               ! Input
     1                  vr_obs_fltrd,                 ! Input/Output
     1                  intvl_rad)                    ! Output

!       Filter to make the filtered ob array more sparse

        implicit none

!       Threshold number of radar obs on a given level
        integer*4 thresh_2_radarobs_lvl_unfltrd
     1           ,thresh_4_radarobs_lvl_unfltrd
!       parameter (thresh_2_radarobs_lvl_unfltrd = 300)
!       parameter (thresh_4_radarobs_lvl_unfltrd = 600)

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

