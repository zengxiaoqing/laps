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

        subroutine read_tsnd(i4time_sys,heights_3d,temp_3d,sh_3d,
     1                   lat_pr,lon_pr,
     1                   lat,lon,
     1                   ob_pr_t,
     1                   c5_name,
     1                   i4_window_raob_file,
!    1                   t_maps_inc,
     1                   bias_htlow,
     1                   n_rass,n_snde,n_tsnd,
     1                   ilaps_cycle_time,
     1                   imax,jmax,kmax,
     1                   r_missing_data
     1                                  )

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

!       Note that max_snd needs to be the same throughout insert_tsnd.f
!                                                     and read_tsnd.f
        integer*4 max_snd,max_snd_levels
        parameter (max_snd = 30)
        parameter (max_snd_levels = 128)

        real*4 surface_rass_buffer
        parameter (surface_rass_buffer = 30.)


!       Declarations
        integer nlevels(max_snd),nlevels_good(max_snd)
        real lat_pr(max_snd)
        real lon_pr(max_snd)
        real elev_pr(max_snd)
        integer num_pr(max_snd)

        real bias_htlow(max_snd)
        real ob_pr_t (max_snd,kmax) ! Vertically interpolated RASS temp
        real ob_pr_ht_obs(max_snd,max_snd_levels)
        real ob_pr_t_obs(max_snd,max_snd_levels)
        character*5 c5_name(max_snd) 
        character*9 a9time

        real*4 heights_3d(imax,jmax,kmax)
        real*4 temp_3d(imax,jmax,kmax)
        real*4 sh_3d(imax,jmax,kmax)
        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)

!       These two arrays (not used yet) serve for incrementing the out of
!       date rass obs according to the model rates of change.
!       real*4 t_maps_inc(imax,jmax,kmax)

        character*13 filename13
        character ext*31
        character*255 c_filespec

!       Initialize

        write(6,*)' Subroutine read_rass -- reads RASS and Sondes'

        n_rass = 0
        n_snde = 0
        n_tsnd = 0

        do i_pr = 1,max_snd
            nlevels(i_pr) = 0
            nlevels_good(i_pr) = 0
        enddo

        do i_pr = 1,max_snd
            do level = 1,kmax
                ob_pr_t(i_pr,level)  = r_missing_data
            enddo
        enddo

        i4_window_rass_ob   = ilaps_cycle_time
        i4_window_rass_file = 3600

! ***   Read in rass data from nearest filetime ******************************

        ext = 'lrs'
        call get_filespec(ext,2,c_filespec,istatus)
        call get_file_time(c_filespec,i4time_sys,i4time_file)

        lag_time = 0 ! Middle of rass hourly sampling period
        i4time_rass_offset = i4time_sys - (i4time_file + lag_time)
        rcycles = float(i4time_rass_offset) / float(ilaps_cycle_time)     

        write(6,*)' i4time_rass_offset/rcycles = '
     1             ,i4time_rass_offset,rcycles

        if(i4time_rass_offset .gt. i4_window_rass_file)then 
            write(6,*)' RASS file is > 60 minutes from LAPS time'       
            write(6,*)' Skipping the use of RASS'
            go to 590
        endif

        call open_lapsprd_file(12,i4time_file,ext,istatus)
        if(istatus .ne. 1)go to 590

400     do i_pr = 1,max_snd

340         continue

            read(12,401,err=406,end=500)
     1     ista,nlevels_in,lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr)
     1                                       ,c5_name(i_pr),a9time
401         format(i12,i12,f11.0,f15.0,f15.0,5x,a5,3x,a9)

!           Determine if rass is in the LAPS domain
406         call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr),lat,lon,i
     1max,jmax
     1                          ,ri,rj,istatus)

            i_ob = nint(ri)
            j_ob = nint(rj)

            write(6,407,err=408)i_pr,ista,nlevels_in
     1                 ,lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr)
     1                 ,i_ob,j_ob,c5_name(i_pr),a9time
407         format(/' RASS #',i3,i6,i5,2f8.2,e10.3,2i4,1x,a5,3x,a9)

408         do level = 1,nlevels_in

                read(12,*,err=340)ht_in,t_in,i_qc

                if(       i_qc .eq. 1
     1             .and.  ht_in .gt. elev_pr(i_pr) + surface_rass_buffer
     1             .and.  level .le. max_snd_levels )then
                    nlevels_good(i_pr) = nlevels_good(i_pr) + 1
                    ob_pr_ht_obs(i_pr,nlevels_good(i_pr)) = ht_in
                    ob_pr_t_obs(i_pr,nlevels_good(i_pr)) =  t_in

c                   write(6,311,err=312)ista,i_pr
c       1                ,ob_pr_ht_obs(i_pr,nlevels_good(i_pr))
c       1                ,ob_pr_t_obs(i_pr,nlevels_good(i_pr))
c311                format(1x,i6,i4,5f8.1)

                endif

312             continue
            enddo ! level

            if(i_ob .ge. 1 .and. i_ob .le. imax .and.
     1         j_ob .ge. 1 .and. j_ob .le. jmax .and.
     1         nlevels_in .ge. 2  ! Upper Lvl profile + sfc temps present
     1                                                  )then
                write(6,*)'  In Bounds - Vertically Interpolating the Ra
     1ss'
            else
                write(6,*)'  Out of Bounds or < 2 levels',nlevels_in
                nlevels_good(i_pr)=0 ! This effectively throws out the Rass
            endif

            call cv_asc_i4time(a9time,i4time_ob)
            if(abs(i4time_ob-i4time_sys) .gt. i4_window_rass_ob)then
                write(6,*)' Out of time bounds:',i4time_ob-i4time_sys
     1                                          ,i4_window_rass_ob
                nlevels_good(i_pr)=0 ! This effectively throws out the Rass    
            endif

            if(nlevels_good(i_pr) .gt. 0)then

!          ***  Interpolate the rass observations to the LAPS grid levels  ****

                do level = 1,kmax

                    t_diff = 0. ! t_maps_inc(i_ob,j_ob,level) * rcycles

                    call interp_rass_to_laps(ob_pr_ht_obs,ob_pr_t_obs,
     1                          t_diff,
     1                          ob_pr_t(i_pr,level),
     1                          i_pr,
     1                          level,
     1                          nlevels_good,
     1                          lat_pr,lon_pr,i_ob,j_ob,
     1                          imax,jmax,kmax,
     1                          max_snd,max_snd_levels,r_missing_data,      
     1                          heights_3d)

c                   write(6,411,err=412)ista,i_pr,level
c       1                ,ob_pr_t(i_pr,level)
c       1                ,temp_3d(i_ob,j_ob,level)
c       1                ,heights_3d(i_ob,j_ob,level)
c       1                ,t_diff
411                 format(1x,i6,2i4,f7.1,1x,f7.1,f8.0,f6.1)

412                 continue
                enddo ! level


!          ***  Interpolate LAPS to the LOWEST RASS level  ****

                do level = 1,1

                    t_diff = 0. ! t_maps_inc(i_ob,j_ob,level) * rcycles

                    call interp_laps_to_rass(ob_pr_ht_obs,ob_pr_t_obs,
     1                         t_interp_laps,p_interp_laps,
     1                         i_pr,
     1                         level,
     1                         nlevels_good,
     1                         lat_pr,lon_pr,i_ob,j_ob,
     1                         imax,jmax,kmax,
     1                         max_snd,max_snd_levels,r_missing_data,     
     1                         temp_3d,heights_3d)

                    call interp_laps_to_rass(ob_pr_ht_obs,ob_pr_t_obs,
     1                         sh_interp_laps,p_interp_laps,
     1                         i_pr,
     1                         level,
     1                         nlevels_good,
     1                         lat_pr,lon_pr,i_ob,j_ob,
     1                         imax,jmax,kmax,
     1                         max_snd,max_snd_levels,r_missing_data,
     1                         sh_3d,heights_3d)

                    tvir = ob_pr_t_obs(i_pr,level) + t_diff
                    p_pa = p_interp_laps
                    tamb = devirt_sh(tvir,sh_interp_laps,p_pa)

                    bias_htlow(i_pr) = tamb - t_interp_laps
!                   bias_htlow(i_pr) = tvir - t_interp_laps

                    write(6,*)' TRASS-VIR= ',tvir
                    write(6,*)' SH       = ',sh_interp_laps
                    write(6,*)' TRASS-AMB= ',tamb
                    write(6,*)' TLAPS    = ',t_interp_laps
                    write(6,*)' Low bias = ',bias_htlow(i_pr)

!d                   write(6,511,err=512)ista,i_pr,level
!d      1                ,ob_pr_t(i_pr,level)
!d      1                ,t_diff
511                 format(1x,i6,2i4,f8.0,1x,f7.1,f6.1)

512                 continue
                enddo ! level

            endif ! # levels > 0 (good rass)

        enddo  ! i_pr
        write(6,*)' WARNING: Used all space in temperature arrays'
        write(6,*)' while reading RASS.  Check max_snd: ',max_snd

500     continue ! Exit out of loop when file is done

        n_rass = i_pr - 1
        close(12)


        goto600
590     write(6,*)' Error opening LRS (or equivalent) file'
        
600     CONTINUE
        write(6,*) ' Read ',n_rass,' RASS temperature sounding(s)'
c
c       Process sounding data
c
        i4time_snd = i4time_sys
        lag_time = 0 ! sounding files are time stamped hourly
        rcycles = float(i4time_sys - i4time_snd + lag_time)
     1                                  / float(ilaps_cycle_time)

! ***   Read in Sonde data  ***************************************

        ext = 'snd'
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

        call open_lapsprd_file(12,i4time_snd,ext,istatus)
        if(istatus .ne. 1)go to 890

        do i_pr = n_rass+1,max_snd

640         continue

            read(12,801,err=706,end=800)
     1      ista,nlevels_in,lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr)
     1                                        ,c5_name(i_pr),a9time
801         format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9)

!           Determine if sonde is in the LAPS domain
706         call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr),lat,lon
     1                              ,imax,jmax,ri,rj,istatus)

            i_ob = nint(ri)
            j_ob = nint(rj)

            write(6,707,err=708)i_pr,ista,nlevels_in
     1                 ,lat_pr(i_pr),lon_pr(i_pr)
     1                 ,elev_pr(i_pr),i_ob,j_ob,c5_name(i_pr),a9time
707         format(/' Sonde #',i3,i6,i5,2f8.2,e10.3,2i4,1x,a5,3x,a9)

708         do level = 1,nlevels_in

                read(12,*,err=640)ht_in,pr_in,t_in,td_in,dd_in,ff_in

                if( abs(t_in)        .lt. 99.
     1             .and.  abs(ht_in) .lt. 1e6
     1             .and.  level      .le. max_snd_levels )then
                    nlevels_good(i_pr) = nlevels_good(i_pr) + 1
                    ob_pr_ht_obs(i_pr,nlevels_good(i_pr)) = ht_in
                    ob_pr_t_obs(i_pr,nlevels_good(i_pr)) =  t_in + 273.15

c                   write(6,611,err=312)ista,i_pr
c       1                ,ob_pr_ht_obs(i_pr,nlevels_good(i_pr))
c       1                ,ob_pr_t_obs(i_pr,nlevels_good(i_pr))
c611                format(1x,i6,i4,5f8.1)

                endif

612             continue
            enddo ! level

            if(i_ob .ge. 1 .and. i_ob .le. imax .and.
     1         j_ob .ge. 1 .and. j_ob .le. jmax .and.
     1         nlevels_in .ge. 2  ! Upper Lvl profile + sfc temps present
     1                                                  )then
                write(6,*)'  In Bounds - Vertically Interpolating the So
     1nde'
            else
                write(6,*)'  Out of Bounds or < 2 levels ',nlevels_in
                nlevels_good(i_pr)=0 ! This effectively throws out the Sonde
            endif

            if(nlevels_good(i_pr) .gt. 0)then

!          ***  Interpolate the sonde observations to the LAPS grid levels  ****

                do level = 1,kmax

                    t_diff = 0. ! t_maps_inc(i_ob,j_ob,level) * rcycles

                    call interp_rass_to_laps(ob_pr_ht_obs,ob_pr_t_obs,
     1                         t_diff,
     1                         ob_pr_t(i_pr,level),
     1                         i_pr,
     1                         level,
     1                         nlevels_good,
     1                         lat_pr,lon_pr,i_ob,j_ob,
     1                         imax,jmax,kmax,
     1                         max_snd,max_snd_levels,r_missing_data,
     1                         heights_3d)

c                   write(6,711,err=712)ista,i_pr,level
c       1                ,ob_pr_t(i_pr,level)
c       1                ,temp_3d(i_ob,j_ob,level)
c       1                ,heights_3d(i_ob,j_ob,level)
c       1                ,t_diff
711                 format(1x,i6,2i4,f7.1,1x,f7.1,f8.0,f6.1)

712                 continue
                enddo ! level


!          ***  Interpolate LAPS to the LOWEST Sonde level  ****

                level=1

                t_diff = 0. ! t_maps_inc(i_ob,j_ob,level) * rcycles

                call interp_laps_to_rass(ob_pr_ht_obs,ob_pr_t_obs,
     1                      t_interp_laps,p_interp_laps,
     1                      i_pr,
     1                      level,
     1                      nlevels_good,
     1                      lat_pr,lon_pr,i_ob,j_ob,
     1                      imax,jmax,kmax,
     1                      max_snd,max_snd_levels,r_missing_data,
     1                      temp_3d,heights_3d)

                tamb = ob_pr_t_obs(i_pr,level) + t_diff

                bias_htlow(i_pr) = tamb - t_interp_laps

                write(6,*)' TSonde-AMB= ',tamb
                write(6,*)' TLAPS    = ',t_interp_laps
                write(6,*)' Low bias = ',bias_htlow(i_pr)

!d                   write(6,811,err=812)ista,i_pr,level
!d      1                ,ob_pr_t(i_pr,level)
!d      1                ,t_diff
811                 format(1x,i6,2i4,f8.0,1x,f7.1,f6.1)

812                 continue

            endif ! # levels > 0

        enddo  ! i_pr
        write(6,*)' ERROR: Used all space in temperature arrays'
        write(6,*)' while reading sondes.  Check max_snd: ',max_snd

800     continue ! Exit out of loop when file is done
        n_snde = i_pr - 1 - n_rass 
        close(12)
        goto 900

890     write(6,*)' Warning: could not open current SND file'

900     n_tsnd = n_rass + n_snde
 
        write(6,*) ' Read ',n_snde,' temperature sonde(s)'
        write(6,*) ' Read ',n_tsnd,' Total RASS+Sonde sounding(s)'

        RETURN
        end
