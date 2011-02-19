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

        subroutine read_tsnd(i4time_sys,heights_3d,temp_bkg_3d,   ! Input
     1                   sh_3d,pres_3d,                           ! Input
     1                   lat_tsnd,lon_tsnd,                       ! Output
     1                   lat,lon,                                 ! Input
     1                   max_snd,max_snd_levels,                  ! Input
     1                   ob_pr_t,inst_err_tsnd,                   ! Output
     1                   c5_name,c8_sndtype,                      ! Output
     1                   l_read_raob,l_3d,                        ! Input
     1                   i4_window_raob_file,                     ! Input
!    1                   t_maps_inc,                              ! Input
     1                   bias_htlow,                              ! Output
     1                   n_rass,n_snde,n_tsnd,                    ! Output
     1                   ilaps_cycle_time,                        ! Input
     1                   imax,jmax,kmax,                          ! Input
     1                   r_missing_data,                          ! Input
     1                   istatus)                                 ! Output

!       1992     Steve Albers   Read RASS data from lrs files
!       1994     Steve Albers   Withold RASS surface ob
!       1994     Steve Albers   Use SH instead of RH
c       1995     Keith Brewster, CAPS, Added reading of .snd files for
c                                soundings
c       1996     Steve Albers   Read nearest LRS file, even if its time does
c                               not exactly match the LAPS analysis time. 
!       1997     Ken Dritz      Change NZ_L_MAX to kmax, making ob_pr_t an
!                               automatic array (resizability change).
!       1997     Ken Dritz      Add r_missing_data as a dummy argument.
!       1998 Jan Steve Albers   General cleanup including error messages.
!                               Improved the handling of observation times. 
!                               The a9time is now being read in for both rass
!                               and raobs. Time thresholding was introduced 
!                               for rass.
!       1998 Feb Steve Albers   Added feature to calculate the height from
!                               the pressure if the height is missing.

        use mem_namelist, ONLY: iwrite_output,radiometer_ht_temp

        real surface_rass_buffer
        parameter (surface_rass_buffer = 30.)


!       Output arrays
        real lat_tsnd(max_snd)
        real lon_tsnd(max_snd)
        real bias_htlow(max_snd)
        real ob_pr_t (max_snd,kmax) ! Vertically interpolated TSND temp
        real inst_err_tsnd(max_snd)
        character*5 c5_name(max_snd) 
        character*8 c8_sndtype(max_snd) 

!       Local arrays
        integer num_pr(max_snd)
        integer nlevels(max_snd),nlevels_good(max_snd)
        real ob_pr_ht_obs(max_snd,max_snd_levels)
        real ob_pr_pr_obs(max_snd,max_snd_levels)
        real ob_pr_t_obs(max_snd,max_snd_levels)
        real elev_tsnd(max_snd)

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

        logical l_read_raob,l_3d,l_string_contains

!       Initialize

        write(6,*)' Subroutine read_tsnd -- reads LRS and SND'

        n_rass = 0
        n_snde = 0
        n_tsnd = 0

        do i_tsnd = 1,max_snd
            nlevels(i_tsnd) = 0
            nlevels_good(i_tsnd) = 0
        enddo

        do i_tsnd = 1,max_snd ! Initialize some output arrays
            inst_err_tsnd(i_tsnd) = r_missing_data
            do level = 1,kmax
                ob_pr_t(i_tsnd,level) = r_missing_data
            enddo
        enddo

        call get_tempob_time_window('LRS',i4_window_ob,istatus)
        if(istatus .ne. 1)return

        i4_window_rass_file = 3600

        ext = 'tmg'
        if(iwrite_output .ge. 1)then
            call open_lapsprd_file(32,i4time_sys,ext,istatus)
            if(istatus .ne. 1)return
        endif

! ***   Read in rass data from nearest filetime ******************************

        ext = 'lrs'
        call upcase(ext,ext_uc)
        call get_filespec(ext,2,c_filespec,istatus)
        call get_file_time(c_filespec,i4time_sys,i4time_file)

        lag_time = 0 ! Middle of rass hourly sampling period
        i4time_rass_offset = i4time_sys - (i4time_file + lag_time)
        rcycles = float(i4time_rass_offset) / float(ilaps_cycle_time)     

        write(6,*)' i4time_rass_offset/rcycles = '
     1             ,i4time_rass_offset,rcycles

        if(abs(i4time_rass_offset) .gt. i4_window_rass_file)then        
            write(6,*)' RASS file/offset is > 60 minutes from LAPS time'    
            write(6,*)' Skipping the use of RASS'
            go to 590
        endif

        call open_lapsprd_file_read(12,i4time_file,ext,istatus)
        if(istatus .ne. 1)go to 590

400     do i_tsnd = 1,max_snd

340         iwrite=1

            read(12,401,err=406,end=500)
     1      ista,nlevels_in,lat_tsnd(i_tsnd),lon_tsnd(i_tsnd)
     1                     ,elev_tsnd(i_tsnd)       
     1                     ,c5_name(i_tsnd),a9time,c8_sndtype(i_tsnd)       
401         format(i12,i12,f11.0,f15.0,f15.0,5x,a5,3x,a9,1x,a8)

            if(nlevels_in .gt. max_snd_levels)then
                write(6,*)' ERROR: too many levels in file ',ext_uc(1:3)       
     1                   ,i_tsnd,nlevels_in,max_snd_levels
                istatus = 0
                return
            endif

!           Determine if rass is in the LAPS domain
406         call latlon_to_rlapsgrid(lat_tsnd(i_tsnd),lon_tsnd(i_tsnd)       
     1                              ,lat,lon    
     1                              ,imax,jmax,ri,rj,istatus)

            i_ob = nint(ri)
            j_ob = nint(rj)

            write(6,407,err=408)ext_uc(1:3),i_tsnd,ista,nlevels_in       
     1                 ,lat_tsnd(i_tsnd),lon_tsnd(i_tsnd)
     1                 ,elev_tsnd(i_tsnd),i_ob,j_ob,c5_name(i_tsnd)
     1                 ,a9time,c8_sndtype(i_tsnd)       
407         format(/' ',a3,' #',i5,i6,i5,2f8.2,e10.3,2i4,1x,a5,3x,a9
     1                       ,1x,a8)

            if(l_string_contains(c8_sndtype(i_tsnd),'SAT',istatus))then       
                inst_err_tsnd(i_tsnd) = 5.0
            else
                inst_err_tsnd(i_tsnd) = 1.0
            endif                

408         do level = 1,nlevels_in

                read(12,*,err=340)ht_in,t_in,i_qc

                if(       i_qc .eq. 1
     1             .and.  ht_in .gt. 
     1                    elev_tsnd(i_tsnd) + surface_rass_buffer
     1             .and.  level .le. max_snd_levels )then
                    nlevels_good(i_tsnd) = nlevels_good(i_tsnd) + 1       

                    ob_pr_ht_obs(i_tsnd,nlevels_good(i_tsnd)) = ht_in       
                    ob_pr_t_obs(i_tsnd,nlevels_good(i_tsnd)) =  t_in
                    ob_pr_pr_obs(i_tsnd,nlevels_good(i_tsnd)) = 
     1                                                 r_missing_data

c                   write(6,311,err=312)ista,i_tsnd
c       1                ,ob_pr_ht_obs(i_tsnd,nlevels_good(i_tsnd))
c       1                ,ob_pr_t_obs(i_tsnd,nlevels_good(i_tsnd))
c311                format(1x,i6,i4,5f8.1)

                endif

312             continue
            enddo ! level

            if(i_ob .ge. 1 .and. i_ob .le. imax .and.
     1         j_ob .ge. 1 .and. j_ob .le. jmax .and.
     1         nlevels_in .ge. 2  ! Upper Lvl profile + sfc temps present
     1                                                  )then
                write(6,*)'  In Bounds - Vertically Interpolating the '
     1                   ,c8_sndtype(i_tsnd)
            else
                if(nlevels_in .lt. 2)then
                    write(6,*)'  Less than 2 levels ',nlevels_in
                else
                    write(6,*)'  Out of Bounds ',nlevels_in
                endif
                nlevels_good(i_tsnd)=0 ! This effectively throws out the Rass
            endif

            call cv_asc_i4time(a9time,i4time_ob)
            if(abs(i4time_ob-i4time_sys) .gt. i4_window_ob)then
                write(6,*)' Out of time bounds:',i4time_ob-i4time_sys
     1                                          ,i4_window_ob
                nlevels_good(i_tsnd)=0 ! This effectively throws out the Rass    
            endif

            if(nlevels_good(i_tsnd) .gt. 0)then

!          ***  Interpolate the rass observations to the LAPS grid levels  ****

                do level = 1,kmax

                    t_diff = 0. ! t_maps_inc(i_ob,j_ob,level) * rcycles

                    call interp_tobs_to_laps(ob_pr_ht_obs,ob_pr_t_obs, ! I
     1                          ob_pr_pr_obs,                          ! I
     1                          t_diff,temp_bkg_3d,                    ! I
     1                          ob_pr_t(i_tsnd,level),                 ! O
     1                          i_tsnd,iwrite,                         ! I
     1                          level,l_3d,                            ! I
     1                          nlevels_good,                          ! I
     1                          lat_tsnd,lon_tsnd,i_ob,j_ob,           ! I
     1                          imax,jmax,kmax,                        ! I
     1                          max_snd,max_snd_levels,                ! I
     1                          r_missing_data,                        ! I
     1                          pres_3d,                               ! I
     1                          heights_3d)                            ! I

c                   write(6,411,err=412)ista,i_tsnd,level
c       1                ,ob_pr_t(i_tsnd,level)
c       1                ,temp_bkg_3d(i_ob,j_ob,level)
c       1                ,heights_3d(i_ob,j_ob,level)
c       1                ,t_diff
411                 format(1x,i6,2i4,f7.1,1x,f7.1,f8.0,f6.1)

                    if(iwrite_output .ge. 1)then
412                     write(32,*)ri,rj,level       
     1                        ,ob_pr_t(i_tsnd,level),c8_sndtype(i_tsnd)       
                    endif
                enddo ! level

            endif ! # levels > 0 (good rass)

        enddo  ! i_tsnd
        write(6,*)' WARNING: Used all space in temperature arrays'
        write(6,*)' while reading RASS.  Check max_snd: '
     1            ,max_snd

500     continue ! Exit out of loop when file is done

        n_rass = i_tsnd - 1
        close(12)


        goto600
590     write(6,*)' Warning, could not open current file ',ext_uc(1:3)       
        
600     CONTINUE
        write(6,*) ' Read ',n_rass,' ',ext_uc(1:3)
     1            ,' temperature sounding(s)'
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

        do i_tsnd = n_rass+1,max_snd

            if(i_tsnd .le. 200)then
                iwrite = 1
            elseif(i_t_snd .le. 1000 .and. 
     1             i_t_snd .eq. (i_tsnd/10)*10)then
                iwrite = 1
            elseif(i_tsnd .eq. (i_tsnd/100)*100)then
                iwrite = 1
            else
                iwrite = 0
            endif

640         continue

            read(12,801,err=706,end=800)
     1      ista,nlevels_in,lat_tsnd(i_tsnd),lon_tsnd(i_tsnd)
     1          ,elev_tsnd(i_tsnd) 
     1          ,c5_name(i_tsnd),a9time,c8_sndtype(i_tsnd)
801         format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)

!           Determine if SND is in the LAPS domain
706         call latlon_to_rlapsgrid(lat_tsnd(i_tsnd),lon_tsnd(i_tsnd)       
     1                              ,lat,lon
     1                              ,imax,jmax,ri,rj,istatus)

            i_ob = nint(ri)
            j_ob = nint(rj)

            if(iwrite .eq. 1)write(6,707,err=708)i_tsnd,ista,nlevels_in       
     1                 ,lat_tsnd(i_tsnd),lon_tsnd(i_tsnd)
     1                 ,elev_tsnd(i_tsnd),i_ob,j_ob,c5_name(i_tsnd)
     1                 ,a9time,c8_sndtype(i_tsnd)
707         format(/' SND #',i5,i6,i5,2f8.2,e10.3,2i4,1x,a5,3x,a9,1x,a8)       

            if(nlevels_in .gt. max_snd_levels)then
                write(6,*)' ERROR: too many levels in SND file '       
     1                   ,i_tsnd,nlevels_in,max_snd_levels
                istatus = 0
                return
            endif

            if(l_string_contains(c8_sndtype(i_tsnd),'SAT',istatus))then       
                inst_err_tsnd(i_tsnd) = 5.0
            else
                inst_err_tsnd(i_tsnd) = 1.0
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

                if(l_string_contains(c8_sndtype(i_tsnd),'RADIO'
     1                                                 ,istatus))then        
                    ht_agl = ht_in - elev_tsnd(i_tsnd)
                    if(ht_agl .gt. radiometer_ht_temp)then
                        i_qc = 0 ! Reject radiometer temps more than 1000m agl
                        write(6,*)' rejecting upper level radiometer'        
     1                           ,level,ht_agl
                    endif                    
                endif

710             if( abs(t_in)        .lt. 99.
     1             .and.  abs(ht_in) .lt. 1e6
     1             .and.  i_qc       .eq. 1
     1             .and.  level      .le. max_snd_levels )then
                    nlevels_good(i_tsnd) = nlevels_good(i_tsnd) + 1       

                    ob_pr_ht_obs(i_tsnd,nlevels_good(i_tsnd)) = ht_in
                    ob_pr_pr_obs(i_tsnd,nlevels_good(i_tsnd)) = pr_in       
                    ob_pr_t_obs(i_tsnd,nlevels_good(i_tsnd)) 
     1                  =  t_in + 273.15

c                   if(iwrite .eq. 1)write(6,611,err=312)ista,i_tsnd       
c       1                ,ob_pr_ht_obs(i_tsnd,nlevels_good(i_tsnd))
c       1                ,ob_pr_t_obs(i_tsnd,nlevels_good(i_tsnd))
c611                format(1x,i6,i4,5f8.1)

                endif

612             continue
            enddo ! level

            call cv_asc_i4time(a9time,i4time_ob)
            if(abs(i4time_ob-i4time_sys) .gt. i4_window_ob)then
                if(iwrite .eq. 1)write(6,*)' Out of time bounds:'
     1                                          ,i4time_ob-i4time_sys       
     1                                          ,i4_window_ob
                nlevels_good(i_tsnd)=0 ! This effectively throws out the sounding

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
                            nlevels_good(i_tsnd)=0 ! effectively throws out Sonde
                        endif

                    else ! less than 1 level
                        if(iwrite .eq. 1)write(6,*)
     1                      '  Less than 1 level ',nlevels_in

                        nlevels_good(i_tsnd)=0 ! effectively throws out Sonde

                    endif

                else ! outside domain

                    if(iwrite .eq. 1)write(6,*)
     1                      '  Out of Bounds ',nlevels_in
                    nlevels_good(i_tsnd)=0 ! effectively throws out Sonde

                endif ! within domain

            endif

            if(nlevels_good(i_tsnd) .gt. 0)then

!          ***  Interpolate the sonde observations to the LAPS grid levels  ****

                do level = 1,kmax

                    t_diff = 0. ! t_maps_inc(i_ob,j_ob,level) * rcycles

                    call interp_tobs_to_laps(ob_pr_ht_obs,ob_pr_t_obs, ! I
     1                         ob_pr_pr_obs,                           ! I
     1                         t_diff,temp_bkg_3d,                     ! I
     1                         ob_pr_t(i_tsnd,level),                  ! O
     1                         i_tsnd,iwrite,                          ! I
     1                         level,l_3d,                             ! I
     1                         nlevels_good,                           ! I
     1                         lat_tsnd,lon_tsnd,i_ob,j_ob,            ! I
     1                         imax,jmax,kmax,                         ! I
     1                         max_snd,max_snd_levels,                 ! I
     1                         r_missing_data,                         ! I
     1                         pres_3d,                                ! I
     1                         heights_3d)                             ! I

c                   write(6,711,err=712)ista,i_tsnd,level
c       1                ,ob_pr_t(i_tsnd,level)
c       1                ,temp_bkg_3d(i_ob,j_ob,level)
c       1                ,heights_3d(i_ob,j_ob,level)
c       1                ,t_diff
711                 format(1x,i6,2i4,f7.1,1x,f7.1,f8.0,f6.1)

                    if(iwrite_output .ge. 1)then
712                     write(32,*)ri,rj,level       
     1                        ,ob_pr_t(i_tsnd,level),c8_sndtype(i_tsnd)       
                    endif
                enddo ! level

            endif ! # levels > 0

        enddo  ! i_tsnd
        write(6,*)' ERROR: Used all space in temperature arrays'
        write(6,*)' while reading SND.  Check max_snd: '
     1            ,max_snd
        istatus = 0

800     continue ! Exit out of loop when file is done
        n_snde = i_tsnd - 1 - n_rass 
        close(12)
        istatus = 1
        goto 900

890     write(6,*)' Warning: could not open current SND file'
        istatus = 1

900     n_tsnd = n_rass + n_snde
 
        write(6,*) ' Read ',n_snde,' SND sounding(s)'
        write(6,*) ' Read ',n_tsnd,' Total LRS+SND sounding(s)'

        RETURN
        end
