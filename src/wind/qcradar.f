
      subroutine qc_radar_obs(
     1           imax,jmax,kmax                             ! Input
     1          ,r_missing_data                             ! Input
     1          ,vr_obs                                     ! Input/Output
     1          ,vr_nyq                                     ! Input
     1          ,n_radarobs_tot_unfltrd                     ! Input
     1          ,lat,lon                                    ! Input
     1          ,rlat_radar,rlon_radar                      ! Input
     1          ,rheight_radar                              ! Input
     1          ,upass1,vpass1  ! 1st pass anal             ! Input
     1          ,u_laps_bkg,v_laps_bkg                      ! Input
     1          ,v_nyquist_2,unfolding_thresh               ! Input
     1          ,l_correct_unfolding,l_grid_north           ! Input
     1          ,istatus                                    ! Input/Output
     1                                                          )

      real*4   vr_obs(imax,jmax,kmax)
      real*4   vr_nyq(imax,jmax,kmax)
      real*4   lat(imax,jmax),lon(imax,jmax)
      real*4   upass1(imax,jmax,kmax),vpass1(imax,jmax,kmax)
      real*4   u_laps_bkg(imax,jmax,kmax),v_laps_bkg(imax,jmax,kmax)

      logical l_correct_unfolding,l_grid_north

      write(6,*)' LVL  # Obs  Intvl  # QC'

      write(6,*)' Also applying QC and dealiasing'
      write(6,*)'  i   j   k    df    vr    fgr'

      height_grid = 0. ! This approximation won't hurt the azimuth

      icount_unfld = 0
      icount_good_qc = 0
      icount_bad_qc = 0

      do k=1,kmax
        height_k = height_of_level(k)
        do j=1,jmax
        do i=1,imax
          if(vr_obs(i,j,k) .ne. r_missing_data)then
            call latlon_to_radar(lat(i,j),lon(i,j),height_k
     1              ,azimuth,slant_range,elev
     1                  ,rlat_radar,rlon_radar,rheight_radar)


            if(abs(upass1(i,j,k)) .ge. 1e6 .or. 
     1         abs(vpass1(i,j,k)) .ge. 1e6)then
                ierr_count = ierr_count + 1
                if(ierr_count .lt. 100)write(6,*)
     1      ' Error in upass1,vpass1',i,j,k,upass1(i,j,k),vpass1(i,j,k)
                istatus = 0
                return

            else
                if(l_grid_north)then
                    call uvgrid_to_radar(
     1                       upass1(i,j,k) + u_laps_bkg(i,j,k),
     1                       vpass1(i,j,k) + v_laps_bkg(i,j,k),
     1                       t_pass1,
     1                       r_pass1,
     1                       azimuth,
     1                       lat(i,j),
     1                       lon(i,j))

                else
                    call uvtrue_to_radar(
     1                       upass1(i,j,k) + u_laps_bkg(i,j,k),
     1                       vpass1(i,j,k) + v_laps_bkg(i,j,k),
     1                       t_pass1,
     1                       r_pass1,
     1                       azimuth)

                endif ! l_grid_north

                call radar_to_uvtrue(t_pass1,
     1                       vr_obs(i,j,k),
     1                       u_true,
     1                       v_true,
     1                       azimuth)

                call uv_to_disp(u_true,
     1                          v_true,
     1                          di_true,
     1                          speed)


!               Compare radar radial velocity to 1st pass analysis
                diff_radial = r_pass1 - vr_obs(i,j,k)

!               Dealias the radar data if appropriate
!               Check for folded radar data (up to +/- 1 Nyquist velocity)
!               Test if residual is more than about 1.3 V Nyquist and
!                                   less than about 2.7 V Nyquist
!
                if(v_nyquist_2 .ne. r_missing_data)then ! Use global Nyquist vel
                  if(abs(abs(diff_radial)-v_nyquist_2) 
     1                                        .lt. unfolding_thresh)then       
c                   call latlon_to_radar(lat(i,j),lon(i,j),height_k
c    1                  ,azimuth,slant_range,elev
c    1                  ,rlat_radar,rlon_radar,rheight_radar)

                    velold = vr_obs(i,j,k)

!                   Adjust the velocity value
                    if(L_correct_unfolding)then
                        r_nyquist_number = nint(diff_radial/v_nyquist_2)
                        diff_radial = diff_radial
     1                                  - r_nyquist_number * v_nyquist_2
                        vr_obs(i,j,k) = vr_obs(i,j,k)
     1                                  + r_nyquist_number * v_nyquist_2
                        icount_unfld=icount_unfld+1

                        if(icount_unfld .le. 50)then
                            write(6,102)i,j,k
     1                                  ,r_pass1
     1                                  ,velold
     1                                  ,vr_obs(i,j,k)
     1                                  ,diff_radial
     1                                  ,nint(azimuth)
     1                                  ,nint(slant_range/1000.)
     1                                  ,elev
                        endif ! icount_unfld < 50
                    endif ! L_correct_unfolding

102                 format(1x,'Folding at'
     1            ,i3,i3,i3,' vp,vr,vrnw',3f6.1,' df',f6.1
     1            ,' azran=',i3,'/',i3,' el=',f4.1)

                  endif ! Data appears folded

                else ! Invalid global Nyquist velocity, try gridpoint Nyquist
                  if (vr_nyq(i,j,k) .gt. 0.) THEN  ! Valid gridpoint Nyquist vel
                    v_nyq_2=2.*vr_nyq(i,j,k)
                    uf_thresh = 0.7 * vr_nyq(i,j,k)
                    if(abs(abs(diff_radial)-v_nyq_2).lt.uf_thresh)then
c                     call latlon_to_radar(lat(i,j),lon(i,j),height_k
c    1                  ,azimuth,slant_range,elev
c    1                  ,rlat_radar,rlon_radar,rheight_radar)

                      velold = vr_obs(i,j,k)

!                     Adjust the velocity value
                      if(L_correct_unfolding)then
                        r_nyquist_number = nint(diff_radial/v_nyq_2)
                        diff_radial = diff_radial
     1                                  - r_nyquist_number * v_nyq_2
                        vr_obs(i,j,k) = vr_obs(i,j,k)
     1                                  + r_nyquist_number * v_nyq_2
                        icount_unfld=icount_unfld+1

                        if(icount_unfld .le. 50)then
                            write(6,102)i,j,k
     1                                  ,r_pass1
     1                                  ,velold
     1                                  ,vr_obs(i,j,k)
     1                                  ,diff_radial
     1                                  ,nint(azimuth)
     1                                  ,nint(slant_range/1000.)
     1                                  ,elev
                        endif ! icount_unfld < 50
                      endif ! L_correct_unfolding
                    endif ! Less than unfolding thresh
                  endif ! Valid gridpoint Nyquist vel
                endif ! We have a valid global nyquist velocity for the radar

!               QC check for radar (current threshold is 12 m/s)
                if(abs(diff_radial) .lt. 12.)then
!                 write(6,320)i,j,k,di_true,speed,u_true,v_true
!320              format(1x,3i2,2f6.1,2f6.1)

!                 write(61,321)lat(i,j),lon(i,j),k,di_true,speed
!321              format(1x,f6.3,f8.3,i2,2f6.1,2f6.1)

                  icount_good_qc = icount_good_qc + 1

                else ! Radar Ob is QC'ed out
                  if(icount_bad_qc .lt. 50)
     1       write(6,311)i,j,k,diff_radial,vr_obs(i,j,k),r_pass1
311               format(3i4,3f6.1,' Radar OB QCed out')
                  vr_obs(i,j,k) = r_missing_data
                  icount_bad_qc = icount_bad_qc + 1

                endif

            endif ! Error condition

          endif ! Not missing data
        enddo ! i
        enddo ! j
      enddo ! k

      write(6,*)' # of Input (unfiltered) Radar obs = '
     1                                      ,n_radarobs_tot_unfltrd

      n_radarobs_tot_unfltrd = icount_good_qc

      write(6,*)' # of Velocities Unfolded = ',icount_unfld
      write(6,*)' # of RADAR GOOD/BAD QC = ',icount_good_qc
     1                                      ,icount_bad_qc
      write(6,*)' # of Output (unfiltered) Radar obs = '
     1                                      ,n_radarobs_tot_unfltrd       

      return
      end

