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

        subroutine read_tdsnd(i4time_sys,heights_3d,temp_bkg_3d,  ! Input
     1                   pres_3d,                                 ! Input
     1                   lat_tdsnd,lon_tdsnd,                     ! Output
     1                   lat,lon,                                 ! Input
     1                   max_snd,max_snd_levels,                  ! Input
     1                   ob_pr_td,inst_err_tdsnd,                 ! Output
     1                   c5_name,c8_sndtype,                      ! Output
     1                   l_read_raob,l_3d,                        ! Input
     1                   i4_window_raob_file,                     ! Input
     1                   bias_htlow,                              ! Output
     1                   n_tdsnd,                                 ! Output
     1                   ilaps_cycle_time,                        ! Input
     1                   imax,jmax,kmax,                          ! Input
     1                   r_missing_data,                          ! Input
     1                   istatus)                                 ! Output

        character*256 path_to_gps

!       Output arrays
        real lat_tdsnd(max_snd)
        real lon_tdsnd(max_snd)
        real bias_htlow(max_snd)
        real ob_pr_td(max_snd,kmax) ! Vertically interpolated TSND temp
        real inst_err_tdsnd(max_snd)
        character*5 c5_name(max_snd) 
        character*8 c8_sndtype(max_snd) 

!       Local arrays
        integer num_pr(max_snd)
        integer nlevels(max_snd),nlevels_good(max_snd)
        real ob_pr_ht_obs(max_snd,max_snd_levels)
        real ob_pr_pr_obs(max_snd,max_snd_levels)
        real ob_pr_td_obs(max_snd,max_snd_levels)
        real elev_tdsnd(max_snd)

        character*9 a9time

        real heights_3d(imax,jmax,kmax)
        real temp_bkg_3d(imax,jmax,kmax)
        real sh_3d(imax,jmax,kmax)
        real pres_3d(imax,jmax,kmax)
        real lat(imax,jmax)
        real lon(imax,jmax)

!       These two arrays (not used yet) serve for incrementing the out of
!       date rass obs according to the model rates of change.
!       real t_maps_inc(imax,jmax,kmax)

        character ext*31
        character ext_uc*31
        character*255 c_filespec

        logical l_read_raob,l_3d,l_string_contains,l_exist

!       GPS obs
        integer gps_n, gps_indomain
        parameter (gps_n = 10000)
        real gps_tpw(gps_n)
        real gps_wet(gps_n)
        real gps_error(gps_n)
        real gps_xy(2,gps_n)
        real gps_elv(gps_n)
        real gps_tim(gps_n)

        write(6,*)
     1  ' Subroutine read_tdsnd -- reads SND and GPS to make HMG'

        lun_hmg = 32

        iwrite_output = 1

        n_rass = 0
        n_snde = 0
        n_tdsnd = 0

        do i_tdsnd = 1,max_snd
            nlevels(i_tdsnd) = 0
            nlevels_good(i_tdsnd) = 0
        enddo

        do i_tdsnd = 1,max_snd ! Initialize some output arrays
            inst_err_tdsnd(i_tdsnd) = r_missing_data
            do level = 1,kmax
                ob_pr_td(i_tdsnd,level) = r_missing_data
            enddo
        enddo

        ext = 'hmg'

        call lapsprd_file_exist(i4time_sys,ext,l_exist,istatus)
        if(l_exist)then
            write(6,*)' HMG file already exists, exiting read_tdsnd'
            istatus = 1
            return
        endif

        if(iwrite_output .ge. 1)then
            call open_lapsprd_file(lun_hmg,i4time_sys,ext,istatus)
            write(6,*)' Could not open file for writing: '
     1               ,lun_hmg,i4time_sys,ext
            if(istatus .ne. 1)return
        endif

        n_rass = 0
c
c       Process sounding data
c
        if(.not. l_read_raob)then
            write(6,*)' Skipping read of SND data'
            istatus = 1
            goto 900
        endif

        i4time_snd = i4time_sys
        lag_time = 0 ! sounding files are time stamped hourly
        rcycles = float(i4time_sys - i4time_snd + lag_time)
     1                                  / float(ilaps_cycle_time)

! ***   Read in SND data  ***************************************

        ext = 'snd'
        call upcase(ext,ext_uc)
        call get_tempob_time_window(ext_uc(1:3),i4_window_ob,istatus)
        if(istatus .ne. 1)return
        call get_filespec(ext,2,c_filespec,istatus)
        call get_file_time(c_filespec,i4time_sys,i4time_nearest)

        i4time_diff = abs(i4time_sys - i4time_nearest)
        if(i4time_diff .le. i4_window_raob_file)then
          write(6,*)' Nearest SND file is within time window'
     1                ,i4time_diff,i4_window_raob_file
        else
          write(6,*)' Nearest SND file is outside time window'
     1                ,i4time_diff,i4_window_raob_file
          go to 890
        endif

        i4time_snd = i4time_nearest

        call open_lapsprd_file_read(12,i4time_snd,ext,istatus)
        if(istatus .ne. 1)go to 890

        do i_tdsnd = n_rass+1,max_snd

            if(i_tdsnd .le. 200)then
                iwrite = 1
            elseif(i_tdsnd .le. 1000 .and. 
     1             i_tdsnd .eq. (i_tdsnd/10)*10)then
                iwrite = 1
            elseif(i_tdsnd .eq. (i_tdsnd/100)*100)then
                iwrite = 1
            else
                iwrite = 0
            endif

640         continue

            read(12,801,err=706,end=800)
     1      ista,nlevels_in,lat_tdsnd(i_tdsnd),lon_tdsnd(i_tdsnd)
     1          ,elev_tdsnd(i_tdsnd) 
     1          ,c5_name(i_tdsnd),a9time,c8_sndtype(i_tdsnd)
801         format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)

!           Determine if SND is in the LAPS domain
706         call latlon_to_rlapsgrid(lat_tdsnd(i_tdsnd)
     1                              ,lon_tdsnd(i_tdsnd)          
     1                              ,lat,lon
     1                              ,imax,jmax,ri,rj,istatus)

            i_ob = nint(ri)
            j_ob = nint(rj)

            if(iwrite .eq. 1)write(6,707,err=708)i_tdsnd,ista,nlevels_in       
     1                 ,lat_tdsnd(i_tdsnd),lon_tdsnd(i_tdsnd)
     1                 ,elev_tdsnd(i_tdsnd),i_ob,j_ob,c5_name(i_tdsnd)
     1                 ,a9time,c8_sndtype(i_tdsnd)
707         format(/' SND #',i5,i6,i5,2f8.2,e10.3,2i4,1x,a5,3x,a9,1x,a8)       

            if(nlevels_in .gt. max_snd_levels)then
                write(6,*)' ERROR: too many levels in SND file '       
     1                   ,i_tdsnd,nlevels_in,max_snd_levels
                istatus = 0
                return
            endif

            if(l_string_contains(c8_sndtype(i_tdsnd),'SAT',istatus))then       
                inst_err_tdsnd(i_tdsnd) = 5.0
            else
                inst_err_tdsnd(i_tdsnd) = 1.0
            endif                

708         do level = 1,nlevels_in

                read(12,*,err=640)ht_in,pr_in,t_in,td_in,dd_in,ff_in

                i_qc = 1

!               Test this by deliberately setting ht_in to missing
!               ht_in = r_missing_data

!               Determine whether we need to supply our own height 
!                                                (only pres given)
                if(ht_in .eq. r_missing_data .and. 
     1             pr_in .ne. r_missing_data                      )then       

                    if(i_ob .ge. 1 .and. i_ob .le. imax .and.
     1                 j_ob .ge. 1 .and. j_ob .le. jmax      )then
                        pr_in_pa = pr_in * 100.
                        call pressure_to_height(pr_in_pa,heights_3d
     1                     ,imax,jmax,kmax,i_ob,j_ob,ht_buff,istatus)
                        if(istatus .ne. 1)goto710
                        ht_in = ht_buff
                        if(iwrite .eq. 1)write(6,*)
     1                      ' Pressure was given, ht was derived:'       
     1                      ,pr_in,ht_in
                    endif
                endif

                if(l_string_contains(c8_sndtype(i_tdsnd),'RADIO'
     1                                                 ,istatus))then        
                    ht_agl = ht_in - elev_tdsnd(i_tdsnd)
                    if(ht_agl .gt. 3000.)then
                        i_qc = 0 ! Reject radiometer temps more than 3000m agl
                        write(6,*)' rejecting upper level radiometer'        
     1                           ,level,ht_agl
                    endif                    
                endif

710             if( abs(td_in)        .lt. 99.
     1             .and.  abs(ht_in) .lt. 1e6
     1             .and.  i_qc       .eq. 1
     1             .and.  level      .le. max_snd_levels )then
                    nlevels_good(i_tdsnd) = nlevels_good(i_tdsnd) + 1       

                    ob_pr_ht_obs(i_tdsnd,nlevels_good(i_tdsnd)) = ht_in
                    ob_pr_pr_obs(i_tdsnd,nlevels_good(i_tdsnd)) = pr_in       
                    ob_pr_td_obs(i_tdsnd,nlevels_good(i_tdsnd)) 
     1                  =  td_in + 273.15

c                   if(iwrite .eq. 1)write(6,611,err=312)ista,i_tdsnd       
c       1                ,ob_pr_ht_obs(i_tdsnd,nlevels_good(i_tdsnd))
c       1                ,ob_pr_td_obs(i_tdsnd,nlevels_good(i_tdsnd))
c611                format(1x,i6,i4,5f8.1)

                endif

612             continue
            enddo ! level

            call cv_asc_i4time(a9time,i4time_ob)
            if(abs(i4time_ob-i4time_sys) .gt. i4_window_ob)then
                if(iwrite .eq. 1)write(6,*)' Out of time bounds:'
     1                                          ,i4time_ob-i4time_sys       
     1                                          ,i4_window_ob
                nlevels_good(i_tdsnd)=0 ! This effectively throws out the sounding

            else ! In time bounds
                if(i_ob .ge. 1 .and. i_ob .le. imax .and.
     1             j_ob .ge. 1 .and. j_ob .le. jmax )then ! within domain

                    if(nlevels_in .ge. 2)then ! Upper Lvl profile 
                                              ! + sfc temps present

                        if(iwrite .eq. 1)write(6,*)
     1                   '  In Bounds - Vertically Interpolating the '       
     1                  ,'Sonde'

                    elseif(nlevels_in .eq. 1)then
                        if(iwrite .eq. 1)write(6,*)
     1                      '  Single level ',nlevels_in

                        if(.not. l_3d)then
                            write(6,*)
     1                      ' ERROR, l_3d is false for 1 level sounding'       
                            nlevels_good(i_tdsnd)=0 ! effectively throws out Sonde
                        endif

                    else ! less than 1 level
                        if(iwrite .eq. 1)write(6,*)
     1                      '  Less than 1 level ',nlevels_in

                        nlevels_good(i_tdsnd)=0 ! effectively throws out Sonde

                    endif

                else ! outside domain

                    if(iwrite .eq. 1)write(6,*)
     1                      '  Out of Bounds ',nlevels_in
                    nlevels_good(i_tdsnd)=0 ! effectively throws out Sonde

                endif ! within domain

            endif

            if(nlevels_good(i_tdsnd) .gt. 0)then

!          ***  Interpolate the sonde observations to the LAPS grid levels  ****

                do level = 1,kmax

                    t_diff = 0. ! t_maps_inc(i_ob,j_ob,level) * rcycles

                    call interp_tobs_to_laps(ob_pr_ht_obs,ob_pr_td_obs, ! I
     1                         ob_pr_pr_obs,                            ! I
     1                         t_diff,temp_bkg_3d,                      ! I
     1                         ob_pr_td(i_tdsnd,level),                 ! O
     1                         i_tdsnd,iwrite,level,l_3d,               ! I
     1                         nlevels_good,                            ! I
     1                         lat_tdsnd,lon_tdsnd,i_ob,j_ob,           ! I
     1                         imax,jmax,kmax,                          ! I
     1                         max_snd,max_snd_levels,                  ! I
     1                         r_missing_data,                          ! I
     1                         pres_3d,                                 ! I
     1                         heights_3d)                              ! I

c                   write(6,711,err=712)ista,i_tdsnd,level
c       1                ,ob_pr_td(i_tdsnd,level)
c       1                ,temp_bkg_3d(i_ob,j_ob,level)
c       1                ,heights_3d(i_ob,j_ob,level)
c       1                ,t_diff
711                 format(1x,i6,2i4,f7.1,1x,f7.1,f8.0,f6.1)

                    if(iwrite_output .ge. 1)then
712                     if(c8_sndtype(i_tdsnd)(1:3) .ne. 'GPS'  .and.
     1                     c8_sndtype(i_tdsnd)(1:4) .ne. 'RAOB' .and.
     1                     c8_sndtype(i_tdsnd)(1:4) .ne. 'RADI'  )then
                            write(lun_hmg,*)ri,rj,level       
     1                      ,ob_pr_td(i_tdsnd,level),c8_sndtype(i_tdsnd) 
     1                      ,' NONAME'
                        else
                            write(lun_hmg,*)ri,rj,level       
     1                      ,ob_pr_td(i_tdsnd,level),c8_sndtype(i_tdsnd) 
     1                      ,' ',c5_name(i_tdsnd)      
                        endif
                    endif
                enddo ! level

            endif ! # levels > 0

        enddo  ! i_tdsnd
        write(6,*)' ERROR: Used all space in temperature arrays'
        write(6,*)' while reading SND.  Check max_snd: '
     1            ,max_snd
        istatus = 0

800     continue ! Exit out of loop when file is done
        n_snde = i_tdsnd - 1 - n_rass 
        close(12)

! ***   Read in GPSdata  ***************************************
        call get_gps_path(path_to_gps,istatus)

        write(6,*)' call read_gps_obs'
        i4beg = i4time_sys - 960
        i4end = i4time_sys + 840
        itry = 1
810     istatus = 0

        call read_gps_obs (lun_hmg, path_to_gps, i4beg, i4end,
     1     imax, jmax, lat, lon, bad_sfc,
     1     gps_tpw, gps_wet, gps_error, gps_xy, gps_elv, 
     1     gps_tim, gps_indomain, gps_n, istatus)
        if(istatus .eq. 1)then
           write(6,*)' Success reading GPS obs, #obs = ',gps_indomain
        else
           if(itry .le. 1)then
              write(6,*)' Failure reading GPS obs - try again...'
              itry = itry + 1
              i4beg = i4time_sys - 0  
              i4end = i4time_sys + 900
              goto 810
           else
              write(6,*)' Failure reading GPS obs'
           endif
        endif

        if(.false.)then
          n_gps = 0
          level = 0
          do i = 1,gps_num
            call latlon_to_rlapsgrid(gps_xy(1,i),gps_xy(2,i),lat,lon
     1                                            ,imax,jmax,ri,rj)
            write(6,*)i,gps_xy(:,i)
            if(ri .ge. 1. .and. ri .le. float(imax) .AND.
     1         rj .ge. 1. .and. rj .le. float(jmax)      )then
                n_gps = n_gps + 1
                write(6      ,*)ri,rj,level       
     1             ,gps_tpw(i),'GPS     '
                write(lun_hmg,*)ri,rj,level       
     1             ,gps_tpw(i),'GPS     '
            endif
          enddo ! i
        endif 

        istatus = 1
        goto 900

890     write(6,*)' Warning: could not open current SND file'
        istatus = 1

900     n_tdsnd = n_rass + n_snde
 
        write(6,*) ' Read ',n_snde,' SND sounding(s)'
        write(6,*) ' Read ',n_tdsnd,' Total LRS+SND sounding(s)'
        write(6,*) ' Read ',gps_indomain,' GPS obs in domain'

        close(lun_hmg)

        RETURN
        end
