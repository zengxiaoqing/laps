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
        subroutine plot_obs(k_level,l_ask_questions,asc9_tim,i_radar
     1    ,imax,jmax,kmax,grid_ra_ref,grid_ra_vel,lat,lon,topo,mode)

!       Steve A         Nov  1989       Original Version
!       Steve A         Nov  1991       Adjustable Dimensions

        include 'trigd.inc'

        real*4 grid_ra_ref(imax,jmax,kmax)
        real*4 grid_ra_vel(imax,jmax,kmax)
        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)
        real*4 topo(imax,jmax)

        character*150 directory
        character*31 ext
        character*10  units_2d
        character*125 comment_2d
        character*13 filename13
        character*3 var_2d
        character*1 c1_plottype

        character*13 fileradar
        character*24 atime_24
        character*9 asc9_tim,asc9_tim_radar
        logical l_ask_questions

        data init/0/

        Real*4
     1  azimuth_deg,
     1  X,
     1  Y,
     1  spd_kt,
     1  SPEED_ms,
     1  DIR,
     1  DU,
     1  PROJROT,
     1  range_km,
     1  pix_per_km,
     1  mspkt

        data
     1  PROJROT/0./,
     1  mspkt/.518/

        Character
     1  c_obs_type*1,c_map*1,c_mode*1,c_anl*1,c_radial*1,
     1  STRING

        common /plotobs/ c_obs_type,c_map,c_mode,c_anl,c_radial

        integer*4               vert_rad_pirep,
     1                          vert_rad_sao,
     1                          vert_rad_meso,
     1                          vert_rad_prof

!       *** set up constants

        call get_r_missing_data(r_missing_data, istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in plot_obs'
        endif

        call get_vert_rads (    vert_rad_pirep,
     1                          vert_rad_sao,
     1                          vert_rad_meso,
     1                          vert_rad_prof,
     1                          istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error in plot_obs'
        endif

!       mode (1: call from 'wd'/'ob' options, 2: call from 'rd'/'vi')

        grid_scale = 3.1
        pix_per_km = 1.65

        x0 = 170.
        y0 = 245.
        pix_per_km = 1.37
        dubase = 2.4

        i_low_offset  = -nint(pix_per_km * 1.51 * grid_scale)
        i_high_offset = i_low_offset 
     1                + nint(pix_per_km * 3.02 * grid_scale)
        dusmall = dubase * pix_per_km/1.65

        size_factor = float(max(imax,jmax)) / 300.

        size_prof = 3.  * size_factor
        size_pirep = 3. * size_factor
        size_maps = 2.  * size_factor
        size_vad = 2.   * size_factor
        size_anl = 0.7  * size_factor
        size_suw = 1.   * size_factor
        size_radar = 1. * size_factor
        size_meso = 2.  * size_factor

        LUN_IN = 5

        call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)

        write(6,*)' Subroutine plot_obs at ',asc9_tim,': Level '
     1           ,k_level

        call setusv_dum(2hIN,201)

!       Plot Radar Obs   *********************************************************************

        write(6,95)
95      format(' Want Radial Velocities     [y,n,a,<RET>=y] ? ',$)
        if(l_ask_questions)read(lun_in,211)c_radial

        if(c_radial(1:1) .eq. 'n')goto205

        retheight = height_of_level(k_level)

        c1_plottype = 'y'

        if(mode .eq. 1)then
            call cv_asc_i4time(asc9_tim,i4time_needed)

!            write(6,110)
!110         format(' Radar #',30x,'? ',$)
!            if(l_ask_questions)read(lun_in,*)i_radar

            if(i_radar .le. 9)then
                write(ext,151)i_radar
151             format('v0',i1)
            else
                write(ext,152)i_radar
152             format('v',i2)
            endif

            var_2d = 'VEL'

            write(6,*)' mode 1, calling get_laps_3dgrid'

            call get_laps_3dgrid(i4time_needed,200000000,i4time_found
     1          ,imax,jmax,kmax,ext,var_2d
     1          ,units_2d,comment_2d,grid_ra_vel,istatus)

            if(istatus.ne.1)goto221

        endif

        do k_grid = 1,kmax

          if(k_grid .gt. k_level .and. k_level .gt. 0)goto221

          do j_grid = 1,jmax
          do i_grid = 1,imax

            k_sfc = nint(height_to_zcoord(topo(i_grid,j_grid),istatus))

            if(k_grid .eq. k_level .or.
     1           k_level .eq. 0 .and. k_sfc .eq. k_grid)then

              if(grid_ra_vel(i_grid,j_grid,k_grid)
     1                          .ne. r_missing_data)then
                call plot_vr(i_grid,j_grid
     1              ,grid_ra_vel(i_grid,j_grid,k_grid),imax,jmax
     1              ,c1_plottype,alat_radar,alon_radar)
              endif

            endif

          enddo ! i_grid
          enddo ! j_grid

        enddo ! k_grid

221     continue

        if(mode .eq. 2)return

!******  Derived Radar Obs ********************************

205     call cv_asc_i4time(asc9_tim,i4time)

        write(6,210)
210     format(' Derived radar obs  [r, Ret = None]',30x,'? ',$)    
        if(l_ask_questions)read(lun_in,211)c_obs_type
211     format(a1)

        if(c_obs_type(1:1) .ne. ' ')then
        if(c_obs_type(1:1) .eq. 'r')then

        write(6,*)
        write(6,*)' Derived Radar Obs'

        call setusv_dum(2hIN,7)

        lun = 61

        if(i_radar .le. 9)then
            write(ext,251)i_radar
 251        format('d0',i1)
        else
            write(ext,252)i_radar
 252        format('d',i2)
        endif

        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=1300)

1211    read(61,*,end=1300)ri,rj,k,dir,speed_ms
!421    format(1x,f6.3,f8.3,i2,2f6.1)
        ri = ri + 1.
        rj = rj + 1.
        k_ob = k + 1.

        k_sfc = nint(height_to_zcoord(topo(nint(ri),nint(rj)),istatus))

        if(k_ob .eq. k_level .or.
     1           k_level .eq. 0 .and. k_sfc .eq. k_ob)then

!               write(6,421)alat,alon,k,dir,speed_ms
!               call latlon_ram_laps(alat,alon,x,y,init,'p')
!               write(6,323)alat,alon,nint(x),nint(y)
1322             format(' old ',2f8.3,2i5)
1323             format(' new ',2f8.3,10x,2i5)

                spd_kt = speed_ms / mspkt

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_radar)

        else if(k_ob .gt. k_level)then
                goto1300

        endif ! k_ob .eq. k_level

        goto1211

1300    continue
        endif


!       Plot SUW obs   ****************************************************************************************************************

        if(c_obs_type(1:1) .eq. 's')then

222     call setusv_dum(2hIN,201)

        lun = 32
        ext = 'suw'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=240)

230     read(32,*,end=240)azimuth_deg,range_km,elev_deg,dir,speed_ms

        range_m = range_km * 1000.

        range_km = range_km * cosd(elev_deg)

        call radar_to_latlon(alat,alon,retheight
     1                      ,azimuth_deg,range_m,elev_deg
     1                  ,rlat_radar,rlon_radar,rheight_radar)

          k = nint(height_to_zcoord(retheight,istatus))

            if ( k .eq. k_level) then ! plot even range_km rings

              if ( range_km .gt. 20. )  then

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = SPEED_ms  / mspkt

c               write(6,112)elev_deg,k,range_km,azimuth_deg,dir,spd_kt
112             format(1x,f6.1,i4,2f7.0,4x,2f7.0,i4)

                call plot_windob(dir,spd_kt,ri,rj,lat,lon
     1                          ,imax,jmax,size_suw)
                endif


              end if ! Valid Wind

        goto230
240     continue
        endif
        endif


        dularge = dusmall * 2.

!       Plot VAD wind  *************************************************
        write(6,*)
        write(6,*)' Plotting VAD Wind'

        call setusv_dum(2hIN,229)

        lun = 21
        ext = 'vad'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=30)

        do while (.true.)
            read(21,*,end=30)el,dir,speed_ms

            retheight = el * 1000. + height_radar

            k = nint(height_to_zcoord(retheight,istatus))

            write(6,*)el,retheight,k

            if(speed_ms .lt. 90. .and. k .eq. k_level)then

                alat = lat_radar
                alon = lon_radar

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = SPEED_ms  / mspkt

                write(6,111)alat,alon,max(dir,-99.),spd_kt
111             format(1x,2f8.1,4x,f7.0,f7.0,i4,f8.3)

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_vad)


            endif
        enddo

 30     continue

!       Plot Mesonet winds  ***********************************************
        if(.false.)then

        write(6,*)' Mesonet Data'

        call setusv_dum(2hIN,14)

        lun = 32
        ext = 'msg'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=50)

!35     read(32,*,end=40)alat,alon,dir,speed_ms,retheight
35      read(32,*,end=40)ri,rj,rk,dir,speed_ms
        ri = ri + 1.
        rj = rj + 1.
        rk = rk + 1.

        k = nint(rk)

        if(k_level .gt. 0)then
            wt_vert = weight_vertical(retheight,k_level,istatus)
        else
            wt_vert = 0.
        endif

        if(abs(k - k_level) .le. vert_rad_meso .or. k_level .eq. 0)then

                write(6,111)ri,rj,max(dir,-99.)
     1                  ,speed_ms,k,wt_vert

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = SPEED_ms  / mspkt

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax,size
     1_meso)

        endif ! k .eq. k_level

        goto35

40      continue

        close(32)

        endif ! .false.

!       Plot Sfc/METAR winds  ***********************************************
50      write(6,*)' Sfc/METAR Data'

        call setusv_dum(2hIN,14)

        lun = 32
        ext = 'sag'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=811)

55      read(32,*,end=60)ri,rj,rk,dir,speed_ms
        ri = ri + 1.
        rj = rj + 1.
        rk = rk + 1.

        k = nint(rk)

        if(k_level .gt. 0)then
            wt_vert = weight_vertical(retheight,k_level,istatus)
        else
            wt_vert = 0.
        endif

        if(abs(k - k_level) .le. vert_rad_sao .or. k_level .eq. 0)then

                write(6,111)ri,rj,max(dir,-99.)
     1                     ,speed_ms,k,wt_vert

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = SPEED_ms  / mspkt

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax,size
     1_meso)

        endif ! k .eq. k_level

        goto55

60      continue

        close(32)

        call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)

!       Plot Profiler winds  ***********************************************
811     write(6,*)
        write(6,*)' Profiler Winds'

        call setusv_dum(2hIN,17)


        lun = 32
        ext = 'prg'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=911)

        do while (.true.)
            read(32,*,end=41)ri,rj,rk,dir,speed_ms
            ri = ri + 1.
            rj = rj + 1.
            rk = rk + 1.

            k = nint(rk)
!           write(6,*)k,alat,alon,retheight,dir,speed_ms
!           k_sfc = nint(height_to_zcoord(topo(i,j),istatus))
            k_sfc = 2

            if(k .eq. k_level
     1  .or. k_level .eq. 0 .and. k .eq. k_sfc
     1  .and. dir .ne. r_missing_data
     1  .and. speed_ms .ne. r_missing_data
     1                                                  )then

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = SPEED_ms  / mspkt

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax,size
     1_prof)
                write(6,111,err=121)alat,alon,dir,spd_kt
121             continue

            endif ! k .eq. k_level

        enddo

41      continue

        close(32)


!       Plot Pirep winds  ***********************************************
911     write(6,*)
        write(6,*)' ACARS Winds'

        lun = 32
        ext = 'pig'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=1011)

        do while (.true.)
            read(32,*,end=51)ri,rj,rk,dir,speed_ms
            ri = ri + 1.
            rj = rj + 1.
            rk = rk + 1.

            k = nint(rk)

            if(abs(k - k_level) .le. vert_rad_pirep)then

                if(k .eq. k_level)then
                    call setusv_dum(2hIN,12)
                else
                    call setusv_dum(2hIN,15)
                endif

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = SPEED_ms  / mspkt

                write(6,111)ri,rj,max(dir,-99.),spd_kt
                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax,size
     1_pirep)

            endif ! k .eq. k_level

        enddo

51      continue

        close(32)



!       Plot MODEL winds  ***********************************************
921     write(6,*)
        write(6,*)' MODEL Winds'

        call setusv_dum(2hIN,203)


        lun = 32
        ext = 'mpg'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=1011)

        do while (.true.)

            read(32,*,end=922)
     1  k_grid,dum,dum,i_grid,j_grid,u_grid,v_grid

            alat = lat(i_grid,j_grid)
            alon = lon(i_grid,j_grid)

            if(k_grid .eq. k_level)then

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')
                call uvgrid_to_disptrue(u_grid,
     1                          v_grid,
     1                          dir,
     1                          speed_ms,
     1                          alon)

                spd_kt = SPEED_ms / mspkt

!               write(6,111)alat,alon,max(dir,-99.),spd_kt
                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax,size
     1_maps)

            endif ! k .eq. k_level

        enddo ! i

922     continue

        close(32)

1011    continue

        return
        end


        subroutine plot_temp_obs(k_level,i4time,imax,jmax,kmax
     1                          ,r_missing_data,lat,lon,topo)

        character*3 ext
        character*150 directory
        character*13 filename13

!       Plot Temperature Obs  ***********************************************

!       size_temp = 8. * float(max(imax,jmax)) / 300.
        size_temp = 3.33

        write(6,*)
        write(6,*)' Plot Temperature Obs, size_temp = ',size_temp

        icol_in = 17
        call setusv_dum(2hIN,17)

        lun = 32
        ext = 'tmg'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=41)

        td = r_missing_data
        p = r_missing_data
        dir = r_missing_data
        spd_kt = r_missing_data
        gust = r_missing_data

        do while (.true.)
            read(32,*,end=41)ri,rj,rk,t_k
            ri = ri + 1.
            rj = rj + 1.
            rk = rk + 1.

            k = nint(rk)
            k_sfc = 2

            if(k .eq. k_level                      .or.
     1         k_level .eq. 0 .and. k .eq. k_sfc   .and.
     1         t_k .ne. r_missing_data                    )then

                t_c = t_k - 273.15

!               spd_kt = SPEED_ms  / mspkt

                iflag = 2

                call plot_mesoob(dir,spd_kt,gust,t_c,td,p,ri,rj
     1                          ,lat,lon,imax,jmax,size_temp,icol_in
     1                          ,iflag)


                write(6,111,err=121)ri,rj,t_c
111             format(1x,3f8.1)
121             continue

            endif ! k .eq. k_level

        enddo

41      continue

        close(32)


        return
        end
