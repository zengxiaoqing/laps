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

      subroutine qc_radar_obs(
     1           imax,jmax,kmax                             ! Input
     1          ,r_missing_data                             ! Input
     1          ,nx_r,ny_r,ioffset,joffset                  ! Input
     1          ,vr_obs                                     ! Input/Output
     1          ,vr_nyq                                     ! Input
     1          ,n_radarobs_tot_unfltrd                     ! Input
     1          ,lat,lon                                    ! Input
     1          ,rlat_radar,rlon_radar                      ! Input
     1          ,rheight_radar                              ! Input
     1          ,upass1,vpass1  ! 1st pass anal             ! Input
     1          ,u_laps_bkg,v_laps_bkg                      ! Input
     1          ,v_nyquist_global                           ! Input
     1          ,l_correct_unfolding,l_grid_north           ! Input
     1          ,istatus                                    ! Input/Output
     1                                                          )

      real   vr_obs(nx_r,ny_r,kmax)
      real   vr_nyq(nx_r,ny_r,kmax)
      real   lat(imax,jmax),lon(imax,jmax)
      real   upass1(imax,jmax,kmax),vpass1(imax,jmax,kmax)
      real   u_laps_bkg(imax,jmax,kmax),v_laps_bkg(imax,jmax,kmax)

      logical l_correct_unfolding,l_grid_north

      parameter (unf_nyq_frac = 0.7)

      write(6,*)' Global Nyquist variable = ',v_nyquist_global

      write(6,*)' LVL  # Obs  Intvl  # QC'

      write(6,*)' Also applying QC and dealiasing'
      write(6,*)'  i   j   k   df    vr    fgr'

      height_grid = 0. ! This approximation won't hurt the azimuth

      icount_unfld = 0
      icount_good_qc = 0
      icount_bad_qc = 0

      if(v_nyquist_global .ne. r_missing_data
     1   .and. v_nyquist_global .lt. 200.
     1   .and. v_nyquist_global .gt. 1.0 ) then
          v_nyquist_2 = 2. * v_nyquist_global
          uf_thresh = v_nyquist_global * unf_nyq_frac 
          v_nyq = v_nyquist_global
      else
          v_nyquist_2 = r_missing_data
          v_nyq = r_missing_data
      endif


      do k=1,kmax
        height_k = height_of_level(k)
        do jo=1,ny_r
        do io=1,nx_r

         i = io + ioffset
         j = jo + joffset

         if(i .ge. 1 .and. i .le. imax .and. 
     1      j .ge. 1 .and. j .le. jmax       )then 

          if(vr_obs(io,jo,k) .ne. r_missing_data)then
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
     1                       vr_obs(io,jo,k),
     1                       u_true,
     1                       v_true,
     1                       azimuth)

                call uv_to_disp(u_true,
     1                          v_true,
     1                          di_true,
     1                          speed)


!               Compare radar radial velocity to 1st pass analysis
                diff_radial = r_pass1 - vr_obs(io,jo,k)

!               Dealias the radar data if appropriate
!               Check for folded radar data (up to +/- 1 Nyquist velocity)
!               Test if residual is more than about 1.3 V Nyquist and
!                                   less than about 2.7 V Nyquist
!
                if(v_nyquist_2 .ne. r_missing_data)then ! Use global Nyquist vel
                    v_nyq_2 = v_nyquist_2
                    v_nyq = v_nyquist_global
                elseif(vr_nyq(io,jo,k) .ne. r_missing_data)then
                    v_nyq_2 = 2. * vr_nyq(io,jo,k)
                    v_nyq = vr_nyq(io,jo,k)
                else
                    v_nyq_2 = r_missing_data
                    v_nyq = r_missing_data
                endif

                if(v_nyq .ne. r_missing_data)then
                  if (v_nyq .gt. 0.) THEN  ! Valid Nyquist vel
                    uf_thresh = unf_nyq_frac * v_nyq

!                   Test if we could benefit from a single de-aliasing
                    if(abs(abs(diff_radial)-v_nyq_2).lt.uf_thresh)then
c                     call latlon_to_radar(lat(i,j),lon(i,j),height_k
c    1                  ,azimuth,slant_range,elev
c    1                  ,rlat_radar,rlon_radar,rheight_radar)

                      velold = vr_obs(io,jo,k)

!                     Adjust the velocity value
                      if(L_correct_unfolding)then
                        r_nyquist_number = nint(diff_radial/v_nyq_2)
                        diff_radial = diff_radial
     1                                  - r_nyquist_number * v_nyq_2
                        vr_obs(io,jo,k) = vr_obs(io,jo,k)
     1                                  + r_nyquist_number * v_nyq_2
                        icount_unfld=icount_unfld+1

                        if(icount_unfld .le. 50)then
                            write(6,102)i,j,k
     1                                  ,r_pass1
     1                                  ,velold
     1                                  ,vr_obs(io,jo,k)
     1                                  ,diff_radial
     1                                  ,nint(azimuth)
     1                                  ,nint(slant_range/1000.)
     1                                  ,elev
102                         format(1x,'Folding at'
     1                        ,i3,i3,i3,' vp,vr,vrnw',3f6.1,' df',f6.1       
     1                        ,' azran=',i3,'/',i3,' el=',f4.1)

                        endif ! icount_unfld < 50
                      endif ! L_correct_unfolding
                    endif ! Less than unfolding thresh
                  endif ! Nyquist vel > 0
                endif ! non-missing nyquist velocity 

!               QC check for radar (current threshold is 12 m/s)
                if(abs(diff_radial) .lt. 12.)then
!                 write(6,320)i,j,k,di_true,speed,u_true,v_true
!320              format(1x,3i2,2f6.1,2f6.1)

!                 write(61,321)lat(i,j),lon(i,j),k,di_true,speed
!321              format(1x,f6.3,f8.3,i2,2f6.1,2f6.1)

                  icount_good_qc = icount_good_qc + 1

                else ! Radar Ob is QC'ed out
                  ! Temporarily comment out the following two lines for writing to save
                  ! output time and file size. Steve or Yuanfu will add a parameter controling what to write
                  ! Yuanfu July 2013
!                  if(icount_bad_qc .lt. 50)
!     1       write(6,311)i,j,k,diff_radial,vr_obs(io,jo,k),r_pass1
311               format(3i4,3f6.1,' Radar OB QCed out')
                  vr_obs(io,jo,k) = r_missing_data
                  icount_bad_qc = icount_bad_qc + 1

                endif

            endif ! Error condition

          endif ! Not missing data
         endif ! i/j within domain
        enddo ! io
        enddo ! jo
      enddo ! k

      write(6,*)' # of Input (unfiltered) Radar obs = '
     1                                      ,n_radarobs_tot_unfltrd

      n_radarobs_tot_unfltrd = icount_good_qc

      write(6,*)' # of Velocities Unfolded = ',icount_unfld
      write(6,604)icount_good_qc,icount_bad_qc
     1           ,pct_rejected(icount_good_qc,icount_bad_qc)
 604  format(' # of RADAR      GOOD/BAD QC = ',2i6,7x
     1      ,'% rejected = ',f6.1)       

      write(6,*)' # of Output (unfiltered) Radar obs = '
     1                                      ,n_radarobs_tot_unfltrd       

      istatus = 1
      return
      end

