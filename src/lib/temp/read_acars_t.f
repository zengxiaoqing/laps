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
        subroutine rd_acars_t(i4time,heights_3d,temp_bkg_3d         ! I
     1                       ,pres_3d                               ! I
     1                       ,MAX_ACARS                             ! I
     1                       ,n_good_acars                          ! O
     1                       ,ext_in                                ! I
!    1                       ,u_maps_inc,v_maps_inc                 ! I
     1                       ,ni,nj,nk                              ! I
     1                       ,lat,lon,r_missing_data                ! I
     1                       ,temp_obs,max_obs,n_obs                ! I/O
     1                       ,istatus)                              ! O

!       1999        Steve Albers, FSL      Called for acars temps

!******************************************************************************

        use mem_namelist, ONLY: iwrite_output

        include 'tempobs.inc'

!       LAPS Grid Dimensions
        real lat(ni,nj)
        real lon(ni,nj)

!       Acars

        integer acars_i(MAX_ACARS)     ! I acars gridpoint
        integer acars_j(MAX_ACARS)     ! J acars gridpoint
        integer acars_ht(MAX_ACARS)    ! HT acars
        real    acars_temp(MAX_ACARS)  ! acars temp

!******************************************************************************

        real heights_3d(ni,nj,nk)
        real pres_3d(ni,nj,nk)
        real temp_bkg_3d(ni,nj,nk)

        real u_maps_inc(ni,nj,nk)
        real v_maps_inc(ni,nj,nk)

        character*9 asc9_tim_acars
        character ext*31, ext_in*3
        character*8 c8_acars_type

        logical l_eof, l_geoalt

        write(6,*)
        write(6,*)' Subroutine rd_acars_t...'

        n_acars_read = 0
        n_acars_obs = 0
        n_good_acars = 0
        n_bad_acars = 0

!       Open input intermediate data file
        lun_in = 31
        call open_lapsprd_file_read(lun_in,i4time,ext_in,istatus)
        if(istatus .ne. 1)go to 999

        lun_tmg = 32
        ext = 'tmg'
        close(lun_tmg)

!       Open output intermediate graphics file
        if(iwrite_output .ge. 1)then
            call open_lapsprd_file_append(lun_tmg,i4time,ext,istatus)       
            if(istatus .ne. 1)go to 888
        endif

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
     1  '   n   i  j  k   T-Raw  T-grid  bias  '

10      i_qc = 1

        if(n_acars_read .le. 200)then
            iwrite = 1
            write(6,*)
        else
            iwrite = 0
        endif

        call read_acars_ob(lun_in,'temp',xlat,xlon,elev_in,temp_ob,arg2       
     1                                  ,asc9_tim_acars,iwrite
     1                                  ,l_geoalt,l_eof)

        if(l_eof)goto900

        if(elev_in .eq. 0.)i_qc = 0

        n_acars_read = n_acars_read + 1

        call cv_asc_i4time(asc9_tim_acars,i4time_acars)

        if(abs(i4time_acars - i4time) .le. i4_window_acars)then ! in time window

            rcycles = float(i4time - i4time_acars) 
     1              / float(ilaps_cycle_time)

!           Climo QC check
            if(temp_ob .lt. 500. .and. i_qc .eq. 1)then

                call latlon_to_rlapsgrid(xlat,xlon,lat,lon,ni,nj
     1                                  ,ri,rj,istatus)
                i_grid = nint(ri)
                j_grid = nint(rj)

                if(i_grid .ge.  1 .and. j_grid .ge. 1 .and.
     1             i_grid .le. ni .and. j_grid .le. nj)then

!                   ACARS is in horizontal domain

                    if(ext_in .eq. 'pin')then

                        if(l_geoalt)then ! ACARS elev is geometric altitude MSL
                            elev_geo = elev_in

                        else !             ACARS elev is pressure altitude MSL
                            elev_std = elev_in

                            if(abs(elev_std) .lt. 90000.)then ! Within flag value
                                pres_mb = ztopsa(elev_std)
                                pres_pa = pres_mb * 100.
                                call pressure_to_height(pres_pa
     1                                             ,heights_3d       
     1                                             ,ni,nj,nk
     1                                             ,i_grid,j_grid
     1                                             ,elev_geo
     1                                             ,istatus_rk)      
                            else
                                istatus_rk = 0
                            endif

                        endif

                        if(istatus_rk .eq. 1)then
                            rk = height_to_zcoord2(elev_geo,heights_3d       
     1                           ,ni,nj,nk,i_grid,j_grid,istatus_rk)       
                        endif

                        if(istatus_rk .ne. 1)then
                            write(6,*)' WARNING: rejecting ACARS ',
     1                      'elevation questionable ',elev_std,elev_geo       
                        endif

                    else 
                        istatus = 0
                        return

                    endif

                    k_grid = nint(rk)

                    if(istatus_rk .eq. 1
     1             .and. k_grid .le. nk
     1             .and. k_grid .ge. 1    )then ! ACARS is in vertical domain

                        n_acars_obs = n_acars_obs + 1

                        if(n_acars_obs .gt. MAX_ACARS)then
                           write(6,*)' Warning: Too many acarss, '
     1                              ,'limit is ',MAX_ACARS
                           istatus = 0
                           return
                        endif

                        acars_i(n_acars_obs) = i_grid
                        acars_j(n_acars_obs) = j_grid

                        acars_ht(n_acars_obs) = elev_geo

                        t_diff = 0.
!                       call get_time_term(t_diff)

                        pres_ob = r_missing_data

                        call interp_tobs_to_laps(
     1                             elev_geo,temp_ob,                     ! I
     1                             pres_ob,                              ! I
     1                             t_diff,temp_bkg_3d,                   ! I
     1                             t_interp,                             ! O
     1                             1,iwrite,k_grid,.true.,               ! I
     1                             1,                                    ! I
     1                             lat_pr,lon_pr,i_grid,j_grid,          ! I
     1                             ni,nj,nk,                             ! I
     1                             1,1,r_missing_data,                   ! I
     1                             pres_3d,                              ! I
     1                             heights_3d)                           ! I

                        if(t_interp .eq. r_missing_data)then
                            write(6,*)' ERROR: t_interp = ',t_interp
                            istatus = 0
                            return
                        endif

                        if(l_geoalt)then ! elev is geometric altitude MSL
                            c8_acars_type = 'WISDOM'
                        else
                            c8_acars_type = 'ACARS'
                        endif

                        if(iwrite_output .ge. 1)then
                            write(lun_tmg,*)ri,rj,k_grid
     1                                     ,t_interp,c8_acars_type
                        endif

!                       Calculate observation bias
                        bias = t_interp - 
     1                         temp_bkg_3d(i_grid,j_grid,k_grid)

!                       QC check of bias
                        if(abs(bias) .le. 10.)then
                            n_good_acars = n_good_acars + 1            
                            n_obs = n_obs + 1

                            if(n_obs .gt. max_obs)then
                                write(6,*)
     1                        ' Error - too many obs in data structure'       
                                write(6,*)
     1                        ' Increase max_obs parameter from',max_obs     
                                istatus = 0
                                return
                            endif

                            call get_temp_obstype(c8_acars_type
     1                                           ,itype_acars,1)

!                           Insert ob into data structure
                            temp_obs(n_obs,i_ri) = i_grid
                            temp_obs(n_obs,i_rj) = j_grid
                            temp_obs(n_obs,i_rk) = rk
                            temp_obs(n_obs,i_ob_raw) = temp_ob
                            temp_obs(n_obs,i_i) = i_grid
                            temp_obs(n_obs,i_j) = j_grid
                            temp_obs(n_obs,i_k) = k_grid
                            temp_obs(n_obs,i_ob_grid) = t_interp
                            temp_obs(n_obs,i_wt) = 1.0
                            temp_obs(n_obs,i_bias) = bias
                            temp_obs(n_obs,i_inst_err) = 1.0
                            temp_obs(n_obs,i_obstype) = itype_acars
 

                        else
                            n_bad_acars = n_bad_acars + 1            
          
                        endif


                        if(iwrite .eq. 1)write(6,20,err=21)
     1                                    n_acars_read,n_acars_obs       
     1                                   ,i_grid,j_grid,k_grid       
     1                                   ,rk ! ,rk_pspace
     1                                   ,temp_ob,t_interp,bias
20                      format(2i5,1x,3i4,2x,f8.3,2x,3f7.1)
21                      continue

                    else
                        write(6,*)' Note: Out of vertical Bounds'
     1                             ,n_acars_read       

                    endif ! In vertical bounds

                else
                    if(n_acars_read .le. 100 .OR.
     1                 n_acars_read .eq. (n_acars_read/100)*100)then
                        write(6,*)' Out of horizontal bounds'
     1                            ,n_acars_read,i_grid,j_grid        
                    endif

                endif ! In horizontal bounds
            endif ! Good data

        else
            write(6,*)' Out of temporal bounds',n_acars_read
     1                              ,abs(i4time_acars - i4time)

        endif ! In temporal bounds

100     goto10

900     write(6,*)' End of ACARS ',ext_in,' file'

        write(6,*)' # of ACARS read in = ',n_acars_read
        write(6,*)' # of ACARS passing bounds checks = ',n_acars_obs      
        write(6,*)' # of ACARS passing QC check = ',n_good_acars
        write(6,*)' # of ACARS failing QC check = ',n_bad_acars
        write(6,*)' % of ACARS failing QC check = ',
     1                     pct_rejected(n_good_acars,n_bad_acars)

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


