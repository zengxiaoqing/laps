cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 

        subroutine insert_tsnd(i4time               ! Input
     1               ,lat,lon                       ! Input
     1               ,heights_3d                    ! Input
     1               ,sh_3d                         ! Input
     1               ,pres_3d                       ! Input
     1               ,temp_3d                       ! Input/Output
     1               ,ilaps_cycle_time              ! Input
     1               ,l_use_raob                    ! Input
     1               ,weight_bkg_const              ! Input
     1               ,i4time_raob_window            ! Input
     1               ,ni,nj,nk                      ! Input
     1               ,grid_spacing_m                ! Input
     1               ,istatus)                      ! Output

!       Nov. 1992               Steve Albers
!       So far, this routine uses only one RASS or combined RASS sounding.
!       The top part of this routine evaluates the biases and does more QC
!       The bottom part inserts the RASS

!       Oct. 1993               Steve Albers
!       Code upgraded to do a Barnes analysis of multiple RASS

!       Dec. 1995               Keith Brewster
!       Changed to call read_tsnd and process temperature data from
!       rawinsondes as well as RASS.

!       Dec 11 1995             Steve Albers
!       Misc bug fixes and variable name changes

!       Jun 16 1997             Ken Dritz
!       Changed NZ_L_MAX to nk.

!       Jun 16 1997             Ken Dritz
!       Added call to get_r_missing_data.  Pass r_missing_data to
!       read_tsnd and analyze_tsnd.

        integer*4 max_snd_grid
        parameter (max_snd_grid = 500)              ! Total number of profiles

        real*4 lat(ni,nj),lon(ni,nj)
        real*4 temp_3d(ni,nj,nk)
        real*4 sh_3d(ni,nj,nk)
        real*4 pres_3d(ni,nj,nk)
        real*4 heights_3d(ni,nj,nk)

        real*4 bias_tsnd(max_snd_grid,nk),bias_htlow(max_snd_grid)
        real*4 wt_tsnd(max_snd_grid,nk)

        real*4 lat_tsnd(max_snd_grid),lon_tsnd(max_snd_grid)
        integer*4 igrid_tsnd(max_snd_grid),jgrid_tsnd(max_snd_grid)
        real*4 tsnd (max_snd_grid,nk) ! Vertically interpolated TSND temp
        character*5 c5_name(max_snd_grid) 
        character*4 c4_obstype(max_snd_grid) 

        logical l_qc,l_flag_vv,l_good_tsnd(max_snd_grid),l_use_raob

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting r_missing_data'
           return
        endif

        do i_pr = 1,max_snd_grid
            l_good_tsnd(i_pr) = .false.
            do k = 1,nk
                bias_tsnd(i_pr,k) = r_missing_data
                wt_tsnd(i_pr,k) = r_missing_data
            enddo ! k
        enddo

!       Read in TSND and Temperature sonde data

        max_snd_obs = max_snd_grid                            ! rass + raobs

        call read_tsnd(i4time,heights_3d,                     ! Input
     1                   temp_3d,sh_3d,pres_3d,               ! Input
     1                   lat_tsnd,lon_tsnd,                   ! Output
     1                   lat,lon,                             ! Input
     1                   max_snd_grid,max_snd_obs,            ! Input
     1                   tsnd,                                ! Output
     1                   c5_name,c4_obstype,                  ! Output
     1                   l_use_raob,                          ! Input
     1                   i4time_raob_window,                  ! Input
!    1                   t_maps_inc,                          ! Input
     1                   bias_htlow,                          ! Output
     1                   n_rass,n_snde,n_tsnd,                ! Output
     1                   ilaps_cycle_time,                    ! Input
     1                   ni,nj,nk,                            ! Input
     1                   r_missing_data,                      ! Input
     1                   istatus)                             ! Output

        if(istatus .ne. 1)then
            write(6,*)' bad istatus returned from read_tsnd'
            return
        endif

        n_good_tsnd = 0

        do i_tsnd = 1,n_tsnd

          iwrite = 1

          call latlon_to_rlapsgrid(lat_tsnd(i_tsnd),lon_tsnd(i_tsnd)
     1                          ,lat,lon,ni,nj,ri,rj,istatus)

          igrid_tsnd(i_tsnd) = nint(ri)
          jgrid_tsnd(i_tsnd) = nint(rj)


          if(istatus .eq. 1)then

!           Find Temperature bias for each level
            write(6,*)
            write(6,*)' Temperature bias, sounding # ',i_tsnd,'  '
     1                ,c5_name(i_tsnd),'  ',c4_obstype(i_tsnd)       
            if(iwrite .eq. 1)write(6,*)
     1      '   k     Tobs        sh      tamb      tlaps      bias'

            l_qc = .false.
!           l_flag_vv = .true.
            l_flag_vv = .false.

            do k = 1,nk
              if(tsnd(i_tsnd,k) .ne. r_missing_data)then

                IF(c4_obstype(i_tsnd) .eq. 'RASS') THEN
!                   Convert from virtual temperature to temperature
                    tvir = tsnd(i_tsnd,k)
                    sh = sh_3d(igrid_tsnd(i_tsnd),jgrid_tsnd(i_tsnd),k)       
                    p_pa = 
     1                 pres_3d(igrid_tsnd(i_tsnd),jgrid_tsnd(i_tsnd),k)    
                    tamb = devirt_sh(tvir,sh,p_pa)
                ELSE
                    sh = 0.
                    tamb = tsnd(i_tsnd,k)
                END IF

                bias_tsnd(i_tsnd,k) =  tamb -
     1                  temp_3d(igrid_tsnd(i_tsnd),jgrid_tsnd(i_tsnd),k)
                wt_tsnd(i_tsnd,k) = 1.0

                if(iwrite .eq. 1)write(6,1)k,tsnd(i_tsnd,k),sh,tamb,
     1                  temp_3d(igrid_tsnd(i_tsnd),jgrid_tsnd(i_tsnd),k)
     1                                ,bias_tsnd(i_tsnd,k)
1               format(i4,f10.1,f10.4,3f10.1)
                if(abs(bias_tsnd(i_tsnd,k)) .gt. 10.)then
                   l_qc = .true.
                   if(iwrite .eq. 1)write(6,*)
     1                 ' ABS(Temp - FIRST GUESS) > 10., Temp NOT USED'       
                endif

!               This should discriminate the vertical velocity data (inactive)
                IF(c4_obstype(i_tsnd) .eq. 'RASS') THEN
                    if(abs(tamb-267.7) .gt. 3.0)then
                        l_flag_vv = .false.
                    endif
                ENDIF

              endif ! Valid data for this TSND at this level
            enddo ! k

!           Find range of levels with a bias value
            k_highest = 0
            k_lowest = nk+1
            do k = 1,nk
                if(bias_tsnd(i_tsnd,k) .ne. r_missing_data)then
                    k_highest = max(k_highest,k)
                    k_lowest  = min(k_lowest,k)
                endif
            enddo ! k

!           Fill in lower parts of Bias field
            if(k_lowest .ge. 2 .and. k_lowest .le. nk)then
                if(c4_obstype(i_tsnd) .ne. 'SAT') THEN
                    do k = 1,k_lowest-1
                        bias_tsnd(i_tsnd,k) = bias_htlow(i_tsnd)
                        wt_tsnd(i_tsnd,k) = 1.0
                    enddo ! k
                endif ! c4_obstype
            endif

!           Add ramp above top of bias data
            if(k_highest .ge. 1 .and. k_highest .le. nk-1)then
                do k = k_highest+1,k_highest+1
                    bias_tsnd(i_tsnd,k) = 
     1              bias_tsnd(i_tsnd,k_highest) * 1.0
                    wt_tsnd(i_tsnd,k) = .04
                enddo ! k
            endif

            write(6,*)' Vertically blended bias field, old/new temps # '
     1                ,i_tsnd
            do k = 1,nk
                if(iwrite .eq. 1)write(6,11,err=12)k,bias_tsnd(i_tsnd,k)       
     1                 ,temp_3d(igrid_tsnd(i_tsnd),jgrid_tsnd(i_tsnd),k)
     1                 ,temp_3d(igrid_tsnd(i_tsnd),jgrid_tsnd(i_tsnd),k)    
     1                + bias_tsnd(i_tsnd,k),wt_tsnd(i_tsnd,k)
11              format(1x,i4,f7.1,3f8.1)
12              continue
            enddo ! k


!           Apply TSND bias to 3D array
            if((.not. l_qc) .and. (.not. l_flag_vv))then
                write(6,*)' Applying the TSND bias corrections'
                l_good_tsnd(i_tsnd) = .true.
                n_good_tsnd = n_good_tsnd + 1
            else
                write(6,*)' Not using TSND data due to QC problems'
                write(6,*)' l_qc =      ',l_qc
                write(6,*)' l_flag_vv = ',l_flag_vv
                l_good_tsnd(i_tsnd) = .false.
                do k = 1,nk
                    bias_tsnd(i_tsnd,k) = r_missing_data
                enddo ! k
            endif

          else  ! failed istatus test
            write(6,*)' TSND station # ',i_tsnd,' is outside of LAPS'  
     1               ,' grid (or missing)'

          endif ! istatus test

        enddo ! i_tsnd

        call analyze_tsnd(n_tsnd,n_good_tsnd,ni,nj,nk,l_good_tsnd
     1      ,weight_bkg_const                                  ! Input
     1      ,grid_spacing_m,max_snd_grid                       ! Input
     1      ,r_missing_data                                    ! Input
     1      ,wt_tsnd,igrid_tsnd,jgrid_tsnd,bias_tsnd,temp_3d,istatus)       

        if(istatus .ne. 1)then
            write(6,*)' Bad istatus returned from analyze_tsnd'
            return
        endif

        istatus = 1
        return
        end


        subroutine analyze_tsnd(n_tsnd,n_good_tsnd,ni,nj,nk,l_good_tsnd     
     1      ,weight_bkg_const                                  ! Input
     1      ,grid_spacing_m,max_snd_grid                       ! Input
     1      ,r_missing_data                                    ! Input
     1      ,wt_tsnd,igrid_tsnd,jgrid_tsnd,bias_tsnd,temp_3d,istatus)       

!       Original Version        Steve Albers

!       Jun 16 1997             Ken Dritz
!       Changed NZ_L_MAX to nk.

!       Jun 16 1997             Ken Dritz
!       Added r_missing_data as dummy argument.  Pass r_missing_data to
!       barnes_univariate_shell.

        real*4 temp_3d(ni,nj,nk)

!       These arrays are passed in
        real*4 bias_3d(ni,nj,nk)
        real*4 wt_3d(ni,nj,nk)
        real*4 bias_obs_3d(ni,nj,nk)
        real*4 bias_tsnd(max_snd_grid,nk)
        real*4 wt_tsnd(max_snd_grid,nk)

        integer*4 igrid_tsnd(max_snd_grid),jgrid_tsnd(max_snd_grid)

        logical l_good_tsnd(max_snd_grid),l_analyze(nk)

        write(6,*)
        write(6,*)' Subroutine analyze_tsnd'

        if(.false.)then ! Insert a single TSND
          do i_tsnd = 1,n_tsnd
            write(6,*)' Look at TSND # ',i_tsnd
            if(l_good_tsnd(i_tsnd))then
                write(6,*)' Potentially adding in the biases'
     1          ,igrid_tsnd(i_tsnd),jgrid_tsnd(i_tsnd)
                if(i_tsnd .eq. 2)then
                    write(6,*)' Actually adding in the biases'
                    do k = 1,nk
                        if(bias_tsnd(i_tsnd,k) .ne. r_missing_data)then
                            do j = 1,nj
                            do i = 1,ni
                                temp_3d(i,j,k) = temp_3d(i,j,k) 
     1                                         + bias_tsnd(i_tsnd,k)
                            enddo ! i
                            enddo ! j
                        endif ! Bias is missing
                    enddo ! k
                endif
            endif
          enddo

        else ! Use Barnes analysis routine (multi-TSND)

          if(n_good_tsnd .gt. 0)then
             write(6,*)' # of TSND stations passing QC = ',n_good_tsnd

           ! This is in effect a single pass Barnes with a spatially varying
           ! radius of influence to account for clustering of data

             call get_fnorm_max(ni,nj,r0_norm,r0_value_min,fnorm_max)
             n_fnorm = int(fnorm_max) + 1

             write(6,*)' Calling new barnes_univariate_shell routine'

             call barnes_univariate_shell(ni,nj,nk           ! Inputs
     1               ,r_missing_data                         ! Input
     1               ,grid_spacing_m                         ! Input
     1               ,max_snd_grid                           ! Input
     1               ,l_good_tsnd,n_tsnd                     ! Inputs
     1               ,bias_tsnd                              ! Input
     1               ,bias_3d                                ! Output
     1               ,bias_obs_3d                            ! Dummy
     1               ,l_analyze                              ! Output
     1               ,wt_3d                                  ! Dummy
     1               ,wt_tsnd,igrid_tsnd,jgrid_tsnd          ! Inputs
     1               ,weight_bkg_const                       ! Input
     1               ,n_fnorm                                ! Input
     1               ,istatus)                               ! Output

             if(istatus .ne. 1)then
                 write(6,*)' Bad status ret fm barnes_univariate'
                 return
             endif

             write(6,*)' Adding back in the biases'
             do k = 1,nk
               if(l_analyze(k))then
                 do j = 1,nj
                 do i = 1,ni
                   if(bias_3d(i,j,k) .ne. r_missing_data)then
                     temp_3d(i,j,k) = temp_3d(i,j,k) + bias_3d(i,j,k)       
                   endif
                 enddo ! i
                 enddo ! j
               endif ! l_analyze(k)
             enddo ! k

          else
             write(6,*)' No Good TSND data, Barnes skipped'
             istatus = 1

          endif ! N GOOD TSNDes > 0

        endif ! Do Barnes

        return
        end


        subroutine barnes_univariate_shell(ni,nj,nk              ! Inputs
     1                   ,r_missing_data                         ! Input
     1                   ,grid_spacing_m                         ! Input
     1                   ,max_snd_grid                           ! Input
     1                   ,l_good_tsnd,n_tsnd                     ! Inputs
     1                   ,bias_tsnd                              ! Input
     1                   ,bias_3d                                ! Output
     1                   ,bias_obs_3d                            ! Input
     1                   ,l_analyze                              ! Output
     1                   ,wt_3d                                  ! Local
     1                   ,wt_tsnd,igrid_tsnd,jgrid_tsnd          ! Inputs
     1                   ,weight_bkg_const                       ! Input
     1                   ,n_fnorm                                ! Input
     1                   ,istatus)                               ! Output

!       Jun 16 1997             Ken Dritz
!       Changed NZ_L_MAX to nk.

!       Jun 16 1997             Ken Dritz
!       Added r_missing_data as dummy argument.

        logical l_good_tsnd(max_snd_grid)
        real*4 bias_tsnd(max_snd_grid,nk)
        integer*4 igrid_tsnd(max_snd_grid),jgrid_tsnd(max_snd_grid)
        real*4 wt_tsnd(max_snd_grid,nk)

        real*4 bias_obs_3d(ni,nj,nk)
        real*4 bias_3d(ni,nj,nk)
        real*4 wt_3d(ni,nj,nk)
        integer*4 n_obs_lvl(nk)                                ! Local

        logical l_analyze(nk)

        integer*4  n_fnorm

        dimension fnorm(0:n_fnorm)

        do k = 1,nk

            l_analyze(k) = .false.

            do i = 1,ni
            do j = 1,nj

                wt_3d(i,j,k) = 0.
                bias_obs_3d(i,j,k) = r_missing_data

            enddo ! j
            enddo ! i

        enddo ! k


        do i_tsnd = 1,n_tsnd

          do k = 1,nk

            if(l_good_tsnd(i_tsnd) .and. 
     1         bias_tsnd(i_tsnd,k) .ne. r_missing_data)then

!               Count obs and determine i,j of obs.

                i = igrid_tsnd(i_tsnd)
                j = jgrid_tsnd(i_tsnd)

                bias_obs_3d(i,j,k) = bias_tsnd(i_tsnd,k)
                wt_3d(i,j,k) = wt_tsnd(i_tsnd,k)

                l_analyze(k) = .true.

                write(6,71)i,j,k,bias_obs_3d(i,j,k),wt_3d(i,j,k)
71              format(1x,3i3,2e11.4)

            endif ! We have an good tsnd

          enddo ! k

        enddo ! Loop through TSND

!        ncnt_total = ncnt_total + ncnt

!       Initialize fnorm array used in barnes_new
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

        n_var = 1

        call barnes_multivariate(
     1                      bias_3d                         ! Outputs
     1                     ,n_var                           ! Input
     1                     ,ni,nj,nk,grid_spacing_m         ! Inputs
     1                     ,bias_obs_3d,wt_3d,fnorm,n_fnorm ! Inputs
     1                     ,l_analyze                       ! Input
     1                     ,weight_bkg_const                ! Input
     1                     ,n_obs_lvl,istatus)              ! Outputs

        return
        end

