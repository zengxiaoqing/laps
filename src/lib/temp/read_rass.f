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

        subroutine read_rass(i4time,heights_3d,temp_3d,sh_3d,
     1                   lat_pr,lon_pr,
     1                   lat,lon,
     1                   ob_pr_t,
!    1                   t_maps_inc,
     1                   bias_htlow,
     1                   n_rass,
     1                   ilaps_cycle_time,
     1                   imax,jmax,kmax
     1                                  )

!       1992 Steve Albers
!       1994 Steve Albers   Withold RASS surface ob
!       1994 Steve Albers   Use SH instead of RH

        include 'lapsparms.inc'

!       Note that max_rs needs to be the same throughout insert_rass.f
!                                                    and read_rass.f
        integer*4 max_rs,max_rs_levels
        parameter (max_rs = 10)
        parameter (max_rs_levels = 128)

        real*4 surface_rass_buffer
        parameter (surface_rass_buffer = 30.)


!       Declarations
        integer nlevels(max_rs),nlevels_good(max_rs)
        real lat_pr(max_rs)
        real lon_pr(max_rs)
        real elev_pr(max_rs)
        integer num_pr(max_rs)

        real bias_htlow(max_rs)
        real ob_pr_t (max_rs,kmax) ! Local, Vertically interpolated RASS temp
        real ob_pr_ht_obs(max_rs,max_rs_levels)
        real ob_pr_t_obs(max_rs,max_rs_levels)

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
        character*5 c5_name

!       Initialize

        write(6,*)' Subroutine read_rass'

        n_rass = 0

        do i_pr = 1,max_rs
            nlevels(i_pr) = 0
            nlevels_good(i_pr) = 0
        enddo

        do i_pr = 1,max_rs
            do level = 1,kmax

                ob_pr_t(i_pr,level)  = r_missing_data

            enddo
        enddo

        i4time_rass = i4time
        lag_time = 0 ! Middle of rass hourly sampling period
        rcycles = float(i4time - i4time_rass + lag_time)
     1                                  / float(ilaps_cycle_time)

! ***   Read in rass data    ***************************************

        ext = 'lrs'
        call open_lapsprd_file(12,i4time_rass,ext,istatus)
        if(istatus .ne. 1)go to 390
        goto400

390     write(6,*)' Trying for RASS data from previous cycle'
        i4time_rass = i4time_rass - ilaps_cycle_time

        call open_lapsprd_file(12,i4time_rass,ext,istatus)
        if(istatus .ne. 1)go to 890

400     do i_pr = 1,max_rs

            read(12,401,err=406,end=500)
     1     ista,nlevels_in,lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr)
     1                                                  ,c5_name
401         format(i12,i12,f11.0,f15.0,f15.0,5x,a5)

!           Determine if rass is in the LAPS domain
406         call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr),lat,lon,i
     1max,jmax
     1                          ,ri,rj,istatus)

            i_ob = nint(ri)
            j_ob = nint(rj)

            write(6,407,err=408)i_pr,ista,nlevels_in,lat_pr(i_pr),lon_pr
     1(i_pr)
     1                  ,elev_pr(i_pr),i_ob,j_ob,c5_name
407         format(/' RASS #',i3,i6,i5,2f8.2,e10.3,2i4,1x,a5)

408         do level = 1,nlevels_in

                read(12,*,err=312)ht_in,t_in,i_qc

                if(       i_qc .eq. 1
     1             .and.  ht_in .gt. elev_pr(i_pr) + surface_rass_buffer
     1                                                          )then
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
     1       j_ob .ge. 1 .and. j_ob .le. jmax .and.
     1    nlevels_in .ge. 2  ! Upper Lvl profile + sfc temps present
     1                                                  )then
                write(6,*)'  In Bounds - Vertically Interpolating the Ra
     1ss'
            else
                write(6,*)'  Out of Bounds or < 2 levels',nlevels_in
                nlevels_good(i_pr)=0 ! This effectively throws out the Rass
            endif

            if(nlevels_good(i_pr) .gt. 0)then

!          ***  Interpolate the rass observations to the LAPS grid levels  ****

                do level = 1,kmax

                    t_diff = 0. ! t_maps_inc(i_ob,j_ob,level) * rcycles

                    call interp_rass_to_laps(ob_pr_ht_obs,ob_pr_t_obs,
     1                             t_diff,
     1                             ob_pr_t(i_pr,level),
     1                             i_pr,
     1                             level,
     1                             nlevels_good,
     1                             lat_pr,lon_pr,i_ob,j_ob,
     1                             imax,jmax,kmax,
     1                             max_rs,max_rs_levels,r_missing_data,
     1                             heights_3d)

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
     1                             t_interp_laps,p_interp_laps,
     1                             i_pr,
     1                             level,
     1                             nlevels_good,
     1                             lat_pr,lon_pr,i_ob,j_ob,
     1                             imax,jmax,kmax,
     1                             max_rs,max_rs_levels,r_missing_data,
     1                             temp_3d,heights_3d)

                    call interp_laps_to_rass(ob_pr_ht_obs,ob_pr_t_obs,
     1                             sh_interp_laps,p_interp_laps,
     1                             i_pr,
     1                             level,
     1                             nlevels_good,
     1                             lat_pr,lon_pr,i_ob,j_ob,
     1                             imax,jmax,kmax,
     1                             max_rs,max_rs_levels,r_missing_data,
     1                             sh_3d,heights_3d)

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

            endif ! # levels > 0

        enddo  ! i_pr

        n_rass = i_pr
        close(12)

        goto900

500     continue ! Exit out of loop when file is done

        n_rass = i_pr - 1
        close(12)

        if(n_rass .ge. max_rs)then
            write(6,*)' WARNING: MAX # of RASSs (parameter max_rs) allow
     1ed = '
     1                                          ,max_rs
        endif

        goto900

880     write(6,*)' Error opening PRG (or equivalent) file'
        return

890     write(6,*)' Error opening LRS (or equivalent) file'
        return

900     return

        end
