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
        subroutine rd_acars_t(i4time,heights_3d
     1  ,N_ACARS,n_acars_obs,ext_in
     1  ,u_maps_inc,v_maps_inc,ni,nj,nk
     1  ,lat,lon
     1  ,acars_i,acars_j,acars_ht,acars_temp
     1  ,temp_obs,max_obs,n_obs                         ! temp data structure
     1                                                  ,istatus)

!       1998        Steve Albers  Called Sequentially for acars temps, 
!                                 then satellite sounding temps.

!******************************************************************************

        include 'tempobs.inc'

!       LAPS Grid Dimensions
        real*4 lat(ni,nj)
        real*4 lon(ni,nj)

!       Acars

        integer acars_i(N_ACARS)  ! X acars coordinates
        integer acars_j(N_ACARS)  ! Y acars coordinates
        integer acars_ht(N_ACARS) ! HT acars
        real    acars_temp(N_ACARS)  ! u acars component

!******************************************************************************

        real*4 heights_3d(ni,nj,nk)
        real*4 u_maps_inc(ni,nj,nk)
        real*4 v_maps_inc(ni,nj,nk)

        character*9 asc9_tim_acars
        character ext*31, ext_in*3

        logical l_eof

        n_good_acars = 0
        n_bad_acars = 0

!       Open input intermediate data file
        lun_in = 31
        call open_lapsprd_file_read(lun_in,i4time,ext_in,istatus)
        if(istatus .ne. 1)go to 999

        lun_tmg = 32
        ext = 'tmg'

!       Open output intermediate graphics file
        call open_lapsprd_file_append(lun_tmg,i4time,ext,istatus)       
        if(istatus .ne. 1)go to 888

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' Error getting laps_cycle_time'
            return
        endif

        call get_tempob_time_window('ACARS',i4_window_acars,istatus)
        if(istatus .eq. 1)then
            write(6,*)' i4_window_acars = ',i4_window_acars
        else
            write(6,*)' Error getting i4_window_acars'
            return
        endif

        write(6,*)
        write(6,*)'             Reading ACARS Obs: ',ext_in
        write(6,*)
     1  '   n   i  j  k   temp    azi    ran '

10      i_qc = 1

!       if(ext_in .eq. 'pin')then
            call read_laps_acars_temp(lun_in,xlat,xlon,elev,temp
     1                                          ,asc9_tim_acars,l_eof)
            if(elev .eq. 0.)i_qc = 0
!       else
!           call read_laps_ssd_temp(lun_in,xlat,xlon,pres,temp
!    1                                          ,asc9_tim_acars,l_eof)
!       endif

        if(l_eof)goto900

        call cv_asc_i4time(asc9_tim_acars,i4time_acars)

        if(abs(i4time_acars - i4time) .le. i4_window_acars)then ! in time window

            rcycles = float(i4time - i4time_acars) 
     1              / float(ilaps_cycle_time)

!           Climo QC check
            if(temp .lt. 500. .and. i_qc .eq. 1)then

                call latlon_to_rlapsgrid(xlat,xlon,lat,lon,ni,nj
     1                                  ,ri,rj,istatus)
                i_grid = nint(ri)
                j_grid = nint(rj)

                if(i_grid .ge.  1 .and. j_grid .ge. 1 .and.
     1             i_grid .le. ni .and. j_grid .le. nj)then

!                   ACARS is in horizontal domain

                    if(ext_in .eq. 'pin')then
!                       Assume ACARS elev is geometric height MSL
                        rk = height_to_zcoord2(elev,heights_3d
     1                       ,ni,nj,nk,i_grid,j_grid,istatus_rk)
                        if(istatus_rk .ne. 1)then
                            write(6,*)' WARNING: rejecting ACARS ',
     1                      'apparently above top of domain ',elev
                        endif

                    else ! ssd
                        istatus = 0
                        return

!                       rk = zcoord_of_pressure(pres)
!                       istatus_rk = 1

                    endif

                    k_grid = nint(rk)

                    if(istatus_rk .eq. 1
     1             .and. k_grid .le. nk
     1             .and. k_grid .ge. 1    )then ! ACARS is in vertical domain

                        n_acars_obs = n_acars_obs + 1

                        if(n_acars_obs .gt. N_ACARS)then
                           write(6,*)' Warning: Too many acarss, '
     1                              ,'limit is ',N_ACARS
                           istatus = 0
                           return
                        endif

                        acars_i(n_acars_obs) = i_grid
                        acars_j(n_acars_obs) = j_grid

                        acars_ht(n_acars_obs) = elev

                        t_diff = 0.
!                       call get_time_term(t_diff)

!                       call interp_tobs_to_laps(
!    1                             elev,temp,                            ! I
!    1                             t_diff,temp_bkg_3d,                   ! I
!    1                             t_interp,                             ! O
!    1                             1,iwrite,level,.true.,                ! I
!    1                             1,                                    ! I
!    1                             lat_pr,lon_pr,i_grid,j_grid,          ! I
!    1                             ni,nj,nk,                             ! I
!    1                             1,1,r_missing_data,                   ! I
!    1                             heights_3d)                           ! I

                        write(lun_tmg,*)ri-1.,rj-1.,rk-1.,t_interp

                        write(6,101)xlat,xlon,dd,ff,rk
     1                             ,t_buff,acars_temp(n_acars_obs)
101                     format(2f8.2,2f8.1,f8.1,3f8.2)

!                       Calculate observation bias
                        bias = t_interp - 
     1                         temp_bkg_3d(i_grid,j_grid,k_grid)

!                       QC check of bias
                        if(abs(bias) .le. 10.)then
                            n_good_acars = n_good_acars + 1            

!                           Insert ob into data structure

                        else
                            n_bad_acars = n_bad_acars + 1            
          
                        endif


                    endif ! In vertical bounds


                    write(6,20)n_acars_obs,
     1                 acars_i(n_acars_obs),
     1                 acars_j(n_acars_obs),
     1                 acars_ht(n_acars_obs),
     1                 acars_temp(n_acars_obs)
20                  format(i4,1x,3i3,2f7.1,2x,2f7.1,2x,f7.1)

                else
                    write(6,*)' Out of horizontal bounds',i_grid,j_grid        

                endif ! In horizontal bounds
            endif ! Good data

        else
            write(6,*)' Out of temporal bounds'
     1                              ,abs(i4time_acars - i4time)

        endif ! In temporal bounds

100     goto10

900     write(6,*)' End of ACARS ',ext_in,' file'

        write(6,*)' # of ACARS read in = ',n_acars_obs
        write(6,*)' # of ACARS passing gross climo check= ',n_acars_obs       
        write(6,*)' # of ACARS passing QC check= ',n_good_acars
        write(6,*)' # of ACARS failing QC check= ',n_bad_acars

        close(lun_in)
        close(lun_tmg)

        istatus = 1

        return

999     write(6,*)' No acars data present'
        istatus = 1
        return


888     write(6,*)' Open error for TMG file'
        istatus = 0
        close(lun_tmg)
        return

        end


        subroutine read_laps_ssd_temp(lun,xlat,xlon,pres,temp
     1                                          ,asc9_tim_acars,l_eof)

        real*4 pres ! pa
        real*4 temp ! degrees K  (99999. is missing)

        character*9 asc9_tim_acars

        logical l_eof

        l_eof = .false.

100     read(lun,895,err=100,end=900)xlat,xlon,pres,dd,ff,asc9_tim_acars       
895     FORMAT(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)

        return

 900    l_eof = .true.

        return
        end

        subroutine read_laps_acars_temp(lun,xlat,xlon,elev,temp
     1                                          ,asc9_tim_acars,l_eof)

        real*4 elev ! meters
        real*4 temp ! degrees K  (99999. is missing)

        character*9 asc9_tim_acars,asc9_tim_rcvd
        character*80 string

        logical l_eof

        dd = 99999.
        ff = 99999.

        l_eof = .false.

5       read(lun,101,end=900,err=5)string(1:6)
101     format(a6)

        if(string(2:5) .eq. 'Time')then
!           a9time = string(30:39)
            read(lun,151)asc9_tim_acars,asc9_tim_rcvd
151         format(1x,a9,2x,a9)
            write(6,151)asc9_tim_acars,asc9_tim_rcvd
        endif

        if(string(2:4) .eq. 'Lat')then
            read(lun,201)xlat,xlon,elev
201         format(2(f8.3,2x), f6.0,2i5)
        endif

        if(string(2:5) .eq. 'Temp')then
            read(lun,202)temp
 202        format (1x, i3,7x, f6.1)
 220        format (' ', i3, ' deg @ ', f6.1, ' m/s')
            write(6,220)temp
            temp = temp
            return
        endif

500     goto5

900     l_eof = .true.

        return

        end

