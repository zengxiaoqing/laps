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
        subroutine plot_obs(k_level,l_ask_questions,asc9_tim
     1    ,i_radar_start,i_radar_end
     1    ,imax,jmax,kmax,n_plotted,grid_ra_ref,grid_ra_vel,lat,lon
     1    ,topo,mode)

!       Steve A         Nov  1989       Original Version
!       Steve A         Nov  1991       Adjustable Dimensions

        include 'trigd.inc'

        real*4 grid_ra_ref(imax,jmax,kmax)
        real*4 grid_ra_vel(imax,jmax,kmax)
        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)
        real*4 topo(imax,jmax)

        integer*4 n_plotted(imax,jmax)

        character*150 directory
        character*31 ext,c3_obsext
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

        write(6,*)' Skipping call to ht_of_level to get retheight'       
!       retheight = ht_of_level(k_level)

        c1_plottype = 'y'

        do i_radar = i_radar_start,i_radar_end

          write(6,*)' mode / radar = ',mode,i_radar

          if(mode .eq. 1)then
            call cv_asc_i4time(asc9_tim,i4time_needed)

            if(i_radar .le. 9)then
                write(ext,151)i_radar
151             format('v0',i1)
            else
                write(ext,152)i_radar
152             format('v',i2)
            endif

            var_2d = 'VEL'

            write(6,*)' calling get_laps_3dgrid'

            call get_laps_3dgrid(i4time_needed,200000000,i4time_found
     1          ,imax,jmax,kmax,ext,var_2d
     1          ,units_2d,comment_2d,grid_ra_vel,istatus)

            if(istatus.ne.1)goto221

            write(6,*)comment_2d

          endif

          do k_grid = 1,kmax

            if(k_grid .gt. k_level .and. k_level .gt. 0)goto221

            do j_grid = 1,jmax
            do i_grid = 1,imax

              if(k_level .eq. 0)then
                k_sfc = 
     1              nint(height_to_zcoord(topo(i_grid,j_grid),istatus))     
              else
                k_sfc = -99
              endif

              if(k_grid  .eq. k_level                   .OR.
     1           k_level .eq. 0 .and. k_sfc .eq. k_grid      )then

                if(grid_ra_vel(i_grid,j_grid,k_grid)
     1                          .ne. r_missing_data)then

                  n_plotted(i_grid,j_grid) = n_plotted(i_grid,j_grid)+1      

                  call plot_vr(i_grid,j_grid
     1              ,grid_ra_vel(i_grid,j_grid,k_grid),imax,jmax
     1              ,c1_plottype,n_plotted(i_grid,j_grid))

                endif

              endif ! at the right level to plot

            enddo ! i_grid
            enddo ! j_grid

          enddo ! k_grid

221       continue

        enddo ! i_radar

        if(mode .eq. 2)return

!******  Derived Radar Obs ********************************

205     call cv_asc_i4time(asc9_tim,i4time)

        write(6,210)
210     format(' Derived radar obs  [r, Ret = None]',30x,'? ',$)    
        if(l_ask_questions)read(lun_in,211)c_obs_type
211     format(a1)

!       if(c_obs_type(1:1) .ne. ' ')then
        if(c_obs_type(1:1) .eq. 'r')then

        write(6,*)
        write(6,*)' Derived Radar Obs'

        call setusv_dum(2hIN,7)

        lun = 61

        do i_radar = i_radar_start,i_radar_end

          if(i_radar .le. 9)then
            write(ext,251)i_radar
 251        format('d0',i1)
          else
            write(ext,252)i_radar
 252        format('d',i2)
          endif

          call get_directory(ext,directory,len_dir)
          open(lun
     1        ,file=directory(1:len_dir)//filename13(i4time,ext(1:3))      
     1        ,status='old',err=1300)

1211      read(61,*,end=1300)ri,rj,k,dir,speed_ms
!421      format(1x,f6.3,f8.3,i2,2f6.1)
          ri = ri + 1.
          rj = rj + 1.
          k_ob = k + 1.

          k_sfc = 
     1    nint(height_to_zcoord(topo(nint(ri),nint(rj)),istatus))

          if(k_ob .eq. k_level .or.
     1       k_level .eq. 0 .and. k_sfc .eq. k_ob)then

!               write(6,421)alat,alon,k,dir,speed_ms
1322            format(' old ',2f8.3,2i5)
1323            format(' new ',2f8.3,10x,2i5)

                spd_kt = speed_ms / mspkt

                write(6,*)nint(ri),nint(rj),dir,spd_kt

!               Note the 'false' passed as DXX winds are grid north
                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_radar,'false')

          elseif(k_ob .gt. k_level)then
                goto1300

          endif ! k_ob .eq. k_level

          goto1211

1300      continue

        enddo ! i_radar

        endif


!       Plot SUW obs   ****************************************************************************************************************

        if(c_obs_type(1:1) .eq. 's')then

222         call setusv_dum(2hIN,201)

            lun = 32
            ext = 'suw'
            call get_directory(ext,directory,len_dir)
            open(lun
     1          ,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1          ,status='old',err=240)

230         read(32,*,end=240)azimuth_deg,range_km,elev_deg,dir,speed_ms

            range_m = range_km * 1000.

            range_km = range_km * cosd(elev_deg)

            call radar_to_latlon(alat,alon,retheight
     1                      ,azimuth_deg,range_m,elev_deg
     1                      ,rlat_radar,rlon_radar,rheight_radar)

            k = nint(height_to_zcoord(retheight,istatus))

            if ( k .eq. k_level) then ! plot even range_km rings

              if ( range_km .gt. 20. )  then

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = SPEED_ms  / mspkt

c               write(6,112)elev_deg,k,range_km,azimuth_deg,dir,spd_kt
112             format(1x,f6.1,i4,2f7.0,4x,2f7.0,i4)

                call plot_windob(dir,spd_kt,ri,rj,lat,lon
     1                          ,imax,jmax,size_suw,'true')
              endif

            end if ! Valid Wind

        goto230
240     continue

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

            retheight = el * 1000. + rheight_radar

            k = nint(height_to_zcoord(retheight,istatus))

            write(6,*)el,retheight,k

            if(speed_ms .lt. 90. .and. k .eq. k_level)then

                alat = rlat_radar
                alon = rlon_radar

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = SPEED_ms  / mspkt

                write(6,111)alat,alon,max(dir,-99.),spd_kt
111             format(1x,2f8.1,4x,f7.0,f7.0,i4,f8.3)

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_vad,'true')


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

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_meso,'true')

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

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_meso,'true')

        endif ! k .eq. k_level

        goto55

60      continue

        close(32)

        call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)

!       Plot Profile winds  ***********************************************
811     write(6,*)
        write(6,*)
     1          ' Profile Winds (e.g. Profiler, Tower, Raob, Dropsonde)'       

        call setusv_dum(2hIN,17)

        lun = 32
        ext = 'prg'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=911)

        do while (.true.)
            read(32,*,end=41)ri,rj,rk,dir,speed_ms
            ri = ri ! + 1.
            rj = rj ! + 1.
            rk = rk ! + 1.

            k = nint(rk)
!           write(6,*)k,alat,alon,retheight,dir,speed_ms
!           k_sfc = nint(height_to_zcoord(topo(i,j),istatus))
            k_sfc = 2

            if(k .eq. k_level
     1    .or. k_level .eq. 0 .and. k .eq. k_sfc
     1   .and. dir .ne. r_missing_data
     1   .and. speed_ms .ne. r_missing_data
     1                                                  )then

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = SPEED_ms  / mspkt

                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_prof,'true')
                write(6,111,err=121)alat,alon,dir,spd_kt
121             continue

            endif ! k .eq. k_level

        enddo

41      continue

        close(32)

!       Plot Pirep winds  ***********************************************
911     write(6,*)
        write(6,*)' ACARS/CDW Winds'

        lun = 32
        ext = 'pig'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=1011)

        do while (.true.)
            read(32,*,end=51)ri,rj,rk,dir,speed_ms,c3_obsext
            ri = ri + 1.
            rj = rj + 1.
            rk = rk + 1.

            k = nint(rk)

            if(abs(k - k_level) .le. vert_rad_pirep)then

                if(k .eq. k_level)then
                    if(c3_obsext .eq. 'pin')then
                        call setusv_dum(2hIN,12)
                    else
                        call setusv_dum(2hIN,8)
                    endif
                else
                    if(c3_obsext .eq. 'pin')then
                        call setusv_dum(2hIN,15)
                    else
                        call setusv_dum(2hIN,34)
                    endif
                endif

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')

                spd_kt = SPEED_ms  / mspkt

                write(6,921)ri,rj,rk,max(dir,-99.),spd_kt,c3_obsext
921             format(1x,3f8.1,4x,f7.0,f7.0,2x,a3)
                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_pirep,'true')

            endif ! k .eq. k_level

        enddo

51      continue

        close(32)



!       Plot MODEL winds  ***********************************************
        write(6,*)
        write(6,*)' MODEL Winds'

        call setusv_dum(2hIN,203)


        lun = 32
        ext = 'mpg'
        call get_directory(ext,directory,len_dir)
        open(lun,file=directory(1:len_dir)//filename13(i4time,ext(1:3))
     1  ,status='old',err=1011)

        do while (.true.)

            read(32,*,end=922)
     1          k_grid,dum,dum,i_grid,j_grid,u_grid,v_grid

            alat = lat(i_grid,j_grid)
            alon = lon(i_grid,j_grid)

            if(k_grid .eq. k_level)then

!               call latlon_ram(alat,alon,x,y,x0,y0,pix_per_km)
!               call latlon_ram_laps(alat,alon,x,y,init,'p')
                call uvgrid_to_disptrue(u_grid,
     1                          v_grid,
     1                          dir,
     1                          speed_ms,
     1                          alat,
     1                          alon)

                spd_kt = SPEED_ms / mspkt

!               write(6,111)alat,alon,max(dir,-99.),spd_kt
                call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                          ,size_maps,'true')

            endif ! k .eq. k_level

        enddo ! i

922     continue

        close(32)

1011    continue

        return
        end


        subroutine plot_temp_obs(k_level,i4time,imax,jmax,kmax
     1                          ,r_missing_data,lat,lon,topo,zoom
     1                          ,plot_parms)

        include 'lapsplot.inc'

        character*3 ext
        character*8 c8_obstype
        character*150 directory
        character*150 filename
        character*13 filename13

!       Plot Temperature Obs  ***********************************************

!       size_temp = 8. * float(max(imax,jmax)) / 300.
        size_temp = 1.1 ! 3.33

        write(6,*)
        write(6,*)' Plot Temperature Obs, size_temp = ',size_temp

        lun = 32
        ext = 'tmg'
        call get_directory(ext,directory,len_dir)

        filename = directory(1:len_dir)//filename13(i4time,ext(1:3))
        call s_len(filename,len_file)

        open(lun,file=filename(1:len_file),status='old',err=31)
        go to 32

 31     write(6,*)' Could not open ',filename(1:len_file)
        go to 41

 32     td = r_missing_data
        p = r_missing_data
        dir = r_missing_data
        spd_kt = r_missing_data
        gust = r_missing_data

        nobs_temp = 0

        do while (.true.) ! Count the temperature obs
            read(32,*,end=41,err=36)ri,rj,rk,t_k,c8_obstype
 36         continue

            ri = ri + 1.
            rj = rj + 1.
            rk = rk + 1.

            k = nint(rk)
            k_sfc = 2

            if( k .eq. k_level  
     1                      .OR.
     1         (k_level .eq. 0 .and. k .eq. k_sfc)  )then

              if(t_k .ne. r_missing_data)then

                nobs_temp = nobs_temp + 1
 
              endif ! t_k .ne. r_missing_data

            endif ! k .eq. k_level

        enddo

41      continue

        rewind(32)

        write(6,*)' Number of temperature obs = ',nobs_temp

        do while (.true.) ! Plot the temp obs
            read(32,*,end=141,err=150)ri,rj,rk,t_k,c8_obstype
150         continue

            ri = ri + 1.
            rj = rj + 1.
            rk = rk + 1.

            k = nint(rk)
            k_sfc = 2

            if( k .eq. k_level  
     1                      .OR.
     1         (k_level .eq. 0 .and. k .eq. k_sfc)  )then

              if(t_k .ne. r_missing_data)then

                t_c = t_k - 273.15

!               spd_kt = SPEED_ms  / mspkt

                iflag = 3

                if(c8_obstype(1:2) .eq. 'RA')then       ! RASS, RAOB
                    icol_in = 12
                elseif(c8_obstype(1:2) .eq. 'SA')then   ! SATSND
                    icol_in = 17
                else                                    ! ACARS
                    icol_in = 3
                endif

                call setusv_dum(2hIN,icol_in)

                iflag_cv = 0

                call plot_mesoob(dir,spd_kt,gust,t_c,td,p,ri,rj
     1                          ,lat,lon,imax,jmax,size_temp
     1                          ,zoom,nobs_temp
     1                          ,icol_in,du_loc,wx
     1                          ,iflag,iflag_cv,plot_parms)


                write(6,111,err=121)ri,rj,t_c,c8_obstype
111             format(1x,3f8.1,1x,a8)
121             continue

              endif ! t_k .ne. r_missing_data

            endif ! k .eq. k_level

        enddo

141     continue

        close(32)


        return
        end
