
        subroutine plot_allsky(i4time_ref,lun,NX_L,NY_L,NZ_L
     1                          ,minalt,maxalt,minazi,maxazi
     1                          ,r_missing_data,laps_cycle_time,maxstns
     1                          ,i_overlay,plot_parms,namelist_parms
     1                          ,l_polar,l_cyl)       

        use mem_namelist, ONLY: max_snd_grid, max_snd_levels
     1                        , grid_spacing_m, aod

        include 'lapsplot.inc'

        addlogs(x,y) = log10(10.**x + 10.**y)

        real pres_3d(NX_L,NY_L,NZ_L)
        real field_3d(NX_L,NY_L,NZ_L)
        real heights_3d(NX_L,NY_L,NZ_L)
        real clwc_3d(NX_L,NY_L,NZ_L)
        real cice_3d(NX_L,NY_L,NZ_L)
        real rain_3d(NX_L,NY_L,NZ_L)
        real snow_3d(NX_L,NY_L,NZ_L)
        real aod_3d(NX_L,NY_L,NZ_L)
!       real rh_3d(NX_L,NY_L,NZ_L)

        parameter (nc = 3)

        real pres_2d(NX_L,NY_L)
        real t_2d(NX_L,NY_L)
        real td_2d(NX_L,NY_L)
        real u_2d(NX_L,NY_L)
        real v_2d(NX_L,NY_L)
        real pw_2d(NX_L,NY_L)
        real cape_2d(NX_L,NY_L)
        real lil_2d(NX_L,NY_L)
        real lic_2d(NX_L,NY_L)
        real swi_2d(NX_L,NY_L)
        real static_albedo(NX_L,NY_L)
        real land_use(NX_L,NY_L)
        real topo_albedo_2d(nc,NX_L,NY_L)
        real lat(NX_L,NY_L)
        real lon(NX_L,NY_L)
        real topo(NX_L,NY_L)
        real sol_alt_2d(NX_L,NY_L)
        real sol_azi_2d(NX_L,NY_L)
        real moon_alt_2d(NX_L,NY_L)
        real moon_azi_2d(NX_L,NY_L)
        real moon_mag,moon_mag_thr

        real temp_vert(max_snd_levels)
        real ht_vert(max_snd_levels)
        real u_vert(max_snd_levels)
        real v_vert(max_snd_levels)
        real rh_vert(max_snd_levels)
        real sh_vert(max_snd_levels) ! dimensionless
        real td_vert(max_snd_levels)
        real lwc_vert(max_snd_levels)
        real ice_vert(max_snd_levels)
        real rain_vert(max_snd_levels)
        real snow_vert(max_snd_levels)
        real pice_vert(max_snd_levels)

        real pres_1d(max_snd_levels)
        real logp_1d(max_snd_levels), logp_bottom, logp_top, logp
     1                              , logp_sfc
        real lil_sfc, lic_sfc, lil_cpt, lic_cpt
  
        real k_to_c, make_td, make_ssh

        character*1 c_prodtype, c_plotobs
        character*3 var_2d
        character*150  directory
        character*31  ext
        character*10  units_2d
        character*125 comment_2d
        character*40 c_model
        character*9 a9time
        character*5 fcst_hhmm
        character*3 c3_string
        character*4 c4_string
        character*40 c_label
        character*16 c16_latlon
        character*11 c_pw
        character*11 c_cape
        character*20 c20_x, c20_y
        character*255 new_dataroot
        logical l_latlon, l_parse, l_plotobs
        logical l_idl /.false./
        logical l_cyl               
        logical l_polar               

        integer i_overlay

        include 'icolors.inc'

!       Sounding observation declarations
        real ob_pr_ht_obs(max_snd_grid,max_snd_levels)
        real ob_pr_pr_obs(max_snd_grid,max_snd_levels)
        real ob_pr_u_obs(max_snd_grid,max_snd_levels) 
        real ob_pr_v_obs(max_snd_grid,max_snd_levels)
        real ob_pr_t_obs(max_snd_grid,max_snd_levels)
        real ob_pr_td_obs(max_snd_grid,max_snd_levels)
        real lat_pr(max_snd_grid)
        real lon_pr(max_snd_grid)
        real elev_pr(max_snd_grid)
        integer nlevels_obs_pr(max_snd_grid)
        integer i4time_ob_pr(max_snd_grid)
        character*5 c5_name, c5_name_a(max_snd_grid), c5_name_min
        character*8 obstype(max_snd_grid)

        real r_cloud_3d(minalt:maxalt,minazi:maxazi) ! cloud opacity
        real cloud_od(minalt:maxalt,minazi:maxazi)   ! cloud optical depth
        real blog_v_roll(minalt:maxalt,minazi:maxazi)
        real blog_moon_roll(minalt:maxalt,minazi:maxazi)
        real glow_stars(minalt:maxalt,minazi:maxazi)
        real elong_roll(minalt:maxalt,minazi:maxazi)
        real airmass_2_cloud_3d(minalt:maxalt,minazi:maxazi)
        real airmass_2_topo_3d(minalt:maxalt,minazi:maxazi)
        real topo_swi(minalt:maxalt,minazi:maxazi)
        real topo_albedo(nc,minalt:maxalt,minazi:maxazi)
        real aod_2_cloud(minalt:maxalt,minazi:maxazi)
        real aod_2_topo(minalt:maxalt,minazi:maxazi)
        real r_cloud_trans(minalt:maxalt,minazi:maxazi)  ! sun to cloud transmissivity (direct+fwd scat)
        real cloud_rad_c(nc,minalt:maxalt,minazi:maxazi) ! sun to cloud transmissivity (direct+fwd scat) * solar color/int
        real clear_rad_c(nc,minalt:maxalt,minazi:maxazi) ! integrated fraction of air illuminated by the sun along line of sight
        real alt_a_roll(minalt:maxalt,minazi:maxazi)
        real azi_a_roll(minalt:maxalt,minazi:maxazi)

        parameter (ni_polar = 511)
        parameter (nj_polar = 511)
        real r_cloud_3d_polar(ni_polar,nj_polar)
        real blog_v_roll_polar(ni_polar,nj_polar)
        real alt_a_polar(ni_polar,nj_polar)
        real azi_a_polar(ni_polar,nj_polar)
        real elong_a_polar(ni_polar,nj_polar)

        real sky_rgb_polar(0:2,ni_polar,nj_polar)
        integer isky_rgb_polar(0:2,ni_polar,nj_polar)
        real sky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi)
        integer isky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi)

        data ilun /0/
        character*3 clun
 
        common /image/ n_image

        skewt(t_c,logp) = t_c - (logp - logp_bottom) * 32.

        if(maxalt .eq. 180)then
            alt_scale = 0.5 
        else
            alt_scale = 1.0 
        endif

        if(maxazi .eq. 720)then
            azi_scale = 0.5 
        else
            azi_scale = 1.0
        endif

        if(l_polar .eqv. .true.)then
            mode_polar = 2
        else
            mode_polar = 0
        endif

        nsmooth = plot_parms%obs_size
        if(nsmooth .ne. 3)then
            nsmooth = 1
        endif

        I4_elapsed = ishow_timer()

        write(6,*)
        write(6,*)' subroutine plot_allsky: nsmooth/aod is ',nsmooth,aod
        write(6,*)' l_cyl/l_polar = ',l_cyl,l_polar
        write(6,*)' minalt/maxalt = ',minalt,maxalt
        write(6,*)' minazi/maxazi = ',minazi,maxazi
        write(6,*)' alt_scale/azi_scale = ',alt_scale,azi_scale

        od_atm_a = aod

        itd = 2 ! dashed dewpoint lines

        n_image = 0

        ext = 'static'

!       Get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

        var_2d='LAT'
        call read_static_grid(NX_L,NY_L,var_2d,lat,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-lat'
            return
        endif

        var_2d='LON'
        call read_static_grid(NX_L,NY_L,var_2d,lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-lon'
            return
        endif

        var_2d='AVG'
        call read_static_grid(NX_L,NY_L,var_2d,topo,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-topo'
            return
        endif

        write(6,*)
        write(6,*)' Input number of all-sky locations...'
        read(5,*)nloc

 80     write(6,*)
        write(6,*)' Input x grid point (or latitude) for all-sky...'
        read(5,*)c20_x

        call s_len(c20_x,lenx)
        l_latlon = l_parse(c20_x(1:lenx),'l')

        if(l_latlon)then ! x value was flagged as latitude with "l" at the end 
            read(c20_x(1:lenx-1),*)soundlat

            write(6,*)' Input longitude for allsky plot...'       
            read(5,*)c20_y
            call s_len(c20_y,leny)
            if(l_parse(c20_y(1:leny),'l'))then
                read(c20_y(1:leny-1),*)soundlon
            else
                read(c20_y(1:leny),*)soundlon
            endif

            call latlon_to_rlapsgrid(soundlat,soundlon,lat,lon
     1                              ,NX_L,NY_L,xsound,ysound,istatus)       

            if(istatus .ne. 1)then
                write(6,*)' Station is outside domain - try again...'
                return
            endif

            if(xsound .lt. 1. .or. xsound .gt. float(NX_L) .OR.
     1         ysound .lt. 1. .or. ysound .gt. float(NY_L)    )then
                write(6,*)' Station is outside domain - try again...'
                return
            endif

        else
            read(c20_x(1:lenx),*)xsound
            write(6,*)' Input y grid point for allsky plot...'
            read(5,*)ysound
            if(xsound .lt. 1.0 .OR. ysound .lt. 1.0)then ! scale domain
                write(6,*)' Values less than 1.0, scale to domain'
                xsound = nint(xsound * (NX_L-1)) + 1
                ysound = nint(ysound * (NY_L-1)) + 1
                soundlat = lat(nint(xsound),nint(ysound))
                soundlon = lon(nint(xsound),nint(ysound))
            endif
        endif

        write(6,*)' soundlat/soundlon ',soundlat,soundlon
        write(6,*)' xsound/ysound ',xsound,ysound

 40     continue
!40     write(6,*)' Enter c_plotobs'
!       read(5,*)c_plotobs
!       if(c_plotobs .eq. '1')then
!           l_plotobs = .true.
!       elseif(c_plotobs .eq. '0')then
            l_plotobs = .false.
!       else
!           write(6,*)' Unknown c_plotobs, will quit ',c_plotobs
!           go to 900
!       endif

        if(.true.)then ! force config with new dataroot
            write(6,*)' Enter new dataroot:'
            read(5,17)new_dataroot
 17         format(a)
            call s_len(new_dataroot,lenroot)

            if(new_dataroot(1:1) .eq. 'q')then
                write(6,*)' Unknown dataroot, will quit'
                go to 900
            else
                write(6,*)' new dataroot is ',new_dataroot(1:lenroot) 
            endif

            call force_get_laps_config(new_dataroot(1:lenroot),istatus)
            if(istatus .ne. 1)then
                write(6,*)' Bad status returned from force_laps_config'
                return
            else
                write(6,*)' Forced config to ',new_dataroot(1:lenroot)
            endif
        endif

        tlow_c = -30.
        thigh_c = +50.

!       Get 3-D pressure field
        call get_pres_3d(i4_valid,NX_L,NY_L,NZ_L,pres_3d,istatus)
        if(istatus .ne. 1)go to 900

        isound = nint(xsound)
        jsound = nint(ysound)

        rlat = lat(isound,jsound)
        rlon = lon(isound,jsound)

        write(c16_latlon,101)rlat,rlon
 101    format(2f8.2)

        do iz = 1,NZ_L
            pres_1d(iz) = pres_3d(isound,jsound,iz)
            logp_1d(iz) = log(pres_1d(iz))
        enddo ! iz

        if(l_plotobs .eqv. .false.)then

          n_lvls_snd = NZ_L

          write(6,*)' i4time_ref = ',i4time_ref

!         Read appropriate 3-D fields
50        call input_product_info(i4time_ref            ! I
     1                         ,laps_cycle_time         ! I
     1                         ,3                       ! I
     1                         ,c_prodtype              ! O
     1                         ,ext                     ! O
     1                         ,directory               ! O
     1                         ,a9time                  ! O
     1                         ,fcst_hhmm               ! O
     1                         ,i4_initial              ! O
     1                         ,i4_valid                ! O
     1                         ,istatus)                ! O

          write(6,*)' a9time is ',a9time

        endif

200     continue

        if(l_plotobs .eqv. .false.)then

!         Read Height
          if(c_prodtype .eq. 'A')then
            iflag_temp = 2 ! Returns Height

            write(6,*)' Calling get_laps_3dgrid for heights'
            var_2d = 'HT'
            ext = 'lt1'
            call get_laps_3dgrid
     1          (i4time_ref,10800,i4time_nearest,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,heights_3d,istatus)

!           call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
!    1                      ,NX_L,NY_L,NZ_L,heights_3d,istatus)
            if(istatus .ne. 1)then
                write(6,*)' Error reading LAPS Heights'
                return
            endif

          elseif(c_prodtype .eq. 'N')then
            call get_directory('balance',directory,len_dir)
            ext = 'lt1'
            directory = directory(1:len_dir)//ext(1:3)

            var_2d = 'HT'

            call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_nearest       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,NX_L,NY_L,NZ_L,heights_3d,istatus)       

          elseif(c_prodtype .eq. 'B' .or. c_prodtype .eq. 'F')then ! Bkg or Fcst
            var_2d = 'HT'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,heights_3d
     1                              ,istatus)
            if(istatus .ne. 1)goto300

          else
            write(6,*)' Unknown choice, will quit'
            go to 900

          endif

          goto400

!         Read RH/SH 
 300      istat_td = 0
          if(c_prodtype .eq. 'A')then ! Read RH
            var_2d = 'RHL'
            ext = 'lh3'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_rh)
            if(istat_rh .ne. 1)goto1000

          elseif(c_prodtype .eq. 'N')then ! Read RH
            call get_directory('balance',directory,len_dir)
            ext = 'lh3'
            directory = directory(1:len_dir)//ext(1:3)

            var_2d = 'RHL'

            call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_nearest       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,NX_L,NY_L,NZ_L,field_3d,istat_rh)       
            if(istat_rh .ne. 1)goto1000

          elseif(c_prodtype .eq. 'B' .or. c_prodtype .eq. 'F')then ! Bkg or Fcst
            var_2d = 'SH'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_sh)
            if(istat_sh .ne. 1)goto1000

          else
            write(6,*)' Sorry, RH/SH not yet supported for prodtype: '
     1               ,c_prodtype
            istat_rh = 0
            goto1000

          endif

          if(c_prodtype .eq. 'A' .or. c_prodtype .eq. 'N')then
            if(nsmooth .gt. 1)then
                call smooth_box_3d(field_3d,NX_L,NY_L,NZ_L,nsmooth)
            endif
            call interp_3d(field_3d,rh_vert,xsound,xsound,ysound,ysound,       
     1                     NX_L,NY_L,NZ_L,1,NZ_L,r_missing_data)

            do iz = 1,NZ_L
                t_c         = k_to_c(temp_vert(iz)) 
                td_vert(iz) = DWPT(t_c,rh_vert(iz))
            enddo ! iz

            istat_td = 1

          else
            if(nsmooth .gt. 1)then
                call smooth_box_3d(field_3d,NX_L,NY_L,NZ_L,nsmooth)
            endif
            call interp_3d(field_3d,sh_vert,xsound,xsound,ysound,ysound,       
     1                     NX_L,NY_L,NZ_L,1,NZ_L,r_missing_data)

            t_ref = -199.
            do iz = 1,NZ_L
                t_c         = k_to_c(temp_vert(iz)) 
                p_mb        = pres_1d(iz)/100.
                q_gkg       = sh_vert(iz) * 1000.
                td_vert(iz) = make_td(p_mb,t_c,q_gkg,t_ref)
                rh_vert(iz) = humidity(t_c,td_vert(iz))
            enddo ! iz

            istat_td = 1

          endif

400       continue

!         Read Cloud Liquid
          istat_lwc = 0
          if(c_prodtype .eq. 'A')then ! Read Cloud Liquid
            var_2d = 'LWC'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_ref,10800,i4time_lwc,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,clwc_3d,istat_lwc)
            call make_fnam_lp(i4time_lwc,a9time,istatus)
            i4time_solar = i4time_lwc
          elseif(c_prodtype .eq. 'F')then 
            var_2d = 'LWC'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,clwc_3d
     1                              ,istat_lwc)
            i4time_solar = i4_valid
          endif

          if(istat_lwc .ne. 1)then
              write(6,*)' Error reading LWC field in plot_allsky'
              return
          endif

!         Calculate solar position for 2D array of grid points
          do i = 1,NX_L
          do j = 1,NY_L
            call solar_position(soundlat,soundlon,i4time_solar
     1                         ,sol_alt_2d(i,j),solar_dec,solar_ha)
            call equ_to_altaz_d(solar_dec,solar_ha,soundlat                 
     1                         ,altdum,sol_azi_2d(i,j))               
            if(sol_azi_2d(i,j) .lt. 0.)sol_azi_2d(i,j) = 
     1                                 sol_azi_2d(i,j) + 360.
          enddo ! j
          enddo ! i

          write(6,*)' solar alt/az (2d array)',sol_alt_2d(isound,jsound)       
     1                                        ,sol_azi_2d(isound,jsound)       

!         Calculate solar position for all-sky point
          call solar_position(soundlat,soundlon,i4time_solar,solar_alt     
     1                                    ,solar_dec,solar_ha)
          call equ_to_altaz_d(solar_dec,solar_ha,soundlat                 
     1                                ,altdum,solar_az)               
          if(solar_az .lt. 0.)solar_az = solar_az + 360.

          write(6,*)' solar alt/az (all-sky)',solar_alt,solar_az

          write(6,*)' call sun_moon at grid point ',isound,jsound
          call sun_moon(i4time_solar,lat,lon,NX_L,NY_L,isound,jsound
     1                 ,alm,azm,elgms,moon_mag)
          write(6,24)alm,azm,elgms,moon_mag 
24        format('  alt/az/elg/mag = ',4f8.2)

!         alm = -90.          ! Test for disabling
!         moon_mag = -4.0    ! Test for disabling
          moon_mag_thr = -6.0

          moon_alt_2d = alm
          moon_azi_2d = azm

          if(istat_lwc .eq. 1)then
            call interp_3d(field_3d,lwc_vert,xsound,xsound
     1                    ,ysound,ysound,NX_L,NY_L,NZ_L,1,NZ_L
     1                    ,r_missing_data)
          else
            lwc_vert = -999.
          endif

!       Read Cloud Ice
          istat_ice = 0
          if(c_prodtype .eq. 'A')then ! Read Cloud Ice
            var_2d = 'ICE'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_lwc,0,i4time_nearest,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,cice_3d,istat_ice)
          elseif(c_prodtype .eq. 'F')then 
            var_2d = 'ICE'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,cice_3d
     1                              ,istat_ice)
!           if(istat_ice .ne. 1)goto1000
          endif

          if(istat_ice .eq. 1)then
            call interp_3d(field_3d,ice_vert,xsound,xsound
     1                    ,ysound,ysound,NX_L,NY_L,NZ_L,1,NZ_L
     1                    ,r_missing_data)
          else
            write(6,*)' Error reading ICE field in plot_allsky'
            return
          endif
 
!         goto500

!         Read Precipitating Rain
          istat_rain = 0
          if(c_prodtype .eq. 'A')then ! Read Precipitating Rain
            var_2d = 'RAI'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_lwc,0,i4time_nearest,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,rain_3d,istat_rain)
          elseif(c_prodtype .eq. 'F')then 
            var_2d = 'RAI'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,rain_3d
     1                              ,istat_rain)
!           if(istat_rain .ne. 1)goto1000
          endif

          if(istat_rain .eq. 1)then
            call interp_3d(field_3d,rain_vert,xsound,xsound
     1                    ,ysound,ysound,NX_L,NY_L,NZ_L,1,NZ_L
     1                    ,r_missing_data)
          else
            write(6,*)' Error reading RAI field in plot_allsky'
            return
          endif

!         Read Precipitating Snow
          istat_snow = 0
          if(c_prodtype .eq. 'A')then ! Read Precipitating Snow
            var_2d = 'SNO'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_lwc,0,i4time_nearest,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,snow_3d,istat_snow)
          elseif(c_prodtype .eq. 'F')then 
            var_2d = 'SNO'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,snow_3d
     1                              ,istat_snow)
!           if(istat_snow .ne. 1)goto1000
          endif

          if(istat_snow .eq. 1)then
            call interp_3d(field_3d,snow_vert,xsound,xsound
     1                    ,ysound,ysound,NX_L,NY_L,NZ_L,1,NZ_L
     1                    ,r_missing_data)
          else
            write(6,*)' Error reading SNO field in plot_allsky'
            return
          endif

          goto500

!         Read Precipitating Ice
          istat_pice = 0
          if(c_prodtype .eq. 'A')then ! Read Precipitating Ice
            var_2d = 'PIC'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_lwc,0,i4time_nearest,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_pice)
          elseif(c_prodtype .eq. 'F')then 
            var_2d = 'PIC'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istat_pice)
!           if(istat_pice .ne. 1)goto1000
          endif

          if(istat_pice .eq. 1)then
            call interp_3d(field_3d,pice_vert,xsound,xsound
     1                    ,ysound,ysound,NX_L,NY_L,NZ_L,1,NZ_L
     1                    ,r_missing_data)
          else
            write(6,*)' Error reading PIC field in plot_allsky'
            return
          endif

500       continue

          ext = 'lcv'

!         Read in swi data
          var_2d = 'SWI'
          call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,swi_2d,0,istat_sfc)
          if(istat_sfc .ne. 1)then
              write(6,*)' Error reading LCV/SWI field in plot_allsky'      
              return
          endif

          call get_static_field_interp('albedo',i4time_lwc,NX_L,NY_L
     1                                ,static_albedo,istat_sfc)
          if(istat_sfc .ne. 1)then
              write(6,*)' Error reading albedo field in plot_allsky'      
              return
          endif

          if(.true.)then ! initial setting of albedo
            do j = 1,NY_L
            do i = 1,NX_L
              if(static_albedo(i,j) .ne. 0.08)then ! brown land
                topo_albedo_2d(1,i,j) = static_albedo(i,j)
                topo_albedo_2d(2,i,j) = static_albedo(i,j)
                topo_albedo_2d(3,i,j) = 0.01
              else ! blue lakes
                topo_albedo_2d(1,i,j) = 0.04
                topo_albedo_2d(2,i,j) = 0.06
                topo_albedo_2d(3,i,j) = 0.12
              endif
            enddo ! i
            enddo ! j

          endif

          if(.true.)then
!           Consider additional albedo info based on land use 
!           Eventually this can include snow cover
            var_2d='USE'
            call read_static_grid(NX_L,NY_L,var_2d,land_use,istatus)
            if(istatus .ne. 1)then
               print*,' Warning: could not read static-landuse'
               return
            endif

            write(6,*)' Set 3-color albedo based on land use'
            call land_albedo(land_use,NX_L,NY_L,topo_albedo_2d)
 
          endif

          goto600

!       Read in sfc data (pressure, temp, dewpoint, u, v, tpw, cape)
          if(c_prodtype .eq. 'A')then ! Read LSX
            ext = 'lsx'

            var_2d = 'PS'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,pres_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'T'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,t_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'TD'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,td_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'U'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,u_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'V'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,v_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            ext = 'lh4'
            var_2d = 'TPW'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,pw_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            ext = 'lst'
            var_2d = 'PBE'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,cape_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            ext = 'lil'
            var_2d = 'LIL'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,lil_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'LIC'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,lic_2d,0,istat_lic)
            if(istat_lic .ne. 1)lic_2d = r_missing_data

          elseif(c_prodtype .eq. 'B' .or. c_prodtype .eq. 'F')then ! Bkg or Fcst
            write(6,*)' Look for Bkg/Fcst sfc fields'

            call s_len(directory,len_dir)

            write(6,*)' Orig directory = ',directory
            write(6,*)' ext = ',ext

            if(c_prodtype .eq. 'B')then
                directory = directory(1:len_dir)//'../lgb'
            else ! Fcst
                ext = 'fsf'
                call s_len(directory,len_dir)
                do i = 1,len_dir-2
                    if(directory(i:i+2) .eq. 'fua')then
                        directory(i:i+2) = ext(1:3)
                        write(6,*)'Substituted directory string'
                    endif
                enddo ! i
            endif

            write(6,*)' New directory = ',directory
            write(6,*)' ext = ',ext

            var_2d = 'PSF'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,NX_L,NY_L
     1                              ,pres_2d
     1                              ,istat_sfc)
!           if(istat_sfc .ne. 1)goto1000

            var_2d = 'TSF'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,NX_L,NY_L
     1                              ,t_2d
     1                              ,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'DSF'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,NX_L,NY_L
     1                              ,td_2d
     1                              ,istat_sfc)
            if(istat_sfc .ne. 1)goto1000

            var_2d = 'TPW'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,NX_L,NY_L
     1                              ,pw_2d
     1                              ,istat_sfc)
            if(istat_sfc .ne. 1)pw_2d = r_missing_data

            var_2d = 'PBE'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,NX_L,NY_L
     1                              ,cape_2d
     1                              ,istat_sfc)
            if(istat_sfc .ne. 1)cape_2d = r_missing_data

            var_2d = 'LIL'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,NX_L,NY_L
     1                              ,lil_2d
     1                              ,istat_sfc)
            if(istat_sfc .ne. 1)lil_2d = r_missing_data

            var_2d = 'LIC'
            call get_lapsdata_2d(i4_initial,i4_valid
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d
     1                              ,NX_L,NY_L
     1                              ,lic_2d
     1                              ,istat_lic)
            if(istat_lic .ne. 1)lic_2d = r_missing_data

          else
            istat_sfc = 0
            go to 1000

          endif

!         Convert sfc variables
          p_sfc_pa = pres_2d(isound,jsound)
          logp_sfc = log(p_sfc_pa)
          t_sfc_k  = t_2d(isound,jsound)
          td_sfc_k = td_2d(isound,jsound)
          pw_sfc   = pw_2d(isound,jsound)
          cape_sfc = cape_2d(isound,jsound)
          lil_sfc  = lil_2d(isound,jsound)
          lic_sfc  = lic_2d(isound,jsound)
600       topo_sfc = topo(isound,jsound)
        endif ! l_plotobs is FALSE

        write(6,*)' a9time is ',a9time

!       Get Atmospheric Optical Depth (3D field)
        call get_aod_3d(pres_3d,heights_3d,topo,NX_L,NY_L,NZ_L
     1                 ,aod_3d)

        do iloc = 1,nloc

          kstart = 0 ! 0 means sfc, otherwise level of start

!         Get line of sight from isound/jsound
          call get_cloud_rays(i4time,clwc_3d,cice_3d,heights_3d
     1                     ,rain_3d,snow_3d
     1                     ,pres_3d,aod_3d,topo_sfc
     1                     ,topo,swi_2d,topo_swi
     1                     ,topo_albedo_2d,topo_albedo
     1                     ,aod_2_cloud,aod_2_topo
     1                     ,r_cloud_3d,cloud_od
     1                     ,r_cloud_trans,cloud_rad_c
     1                     ,clear_rad_c
     1                     ,airmass_2_cloud_3d,airmass_2_topo_3d
     1                     ,NX_L,NY_L,NZ_L,isound,jsound,kstart
     1                     ,view_alt,view_az,sol_alt_2d,sol_azi_2d
     1                     ,moon_alt_2d,moon_azi_2d
     1                     ,moon_mag,moon_mag_thr
     1                     ,minalt,maxalt,minazi,maxazi
     1                     ,alt_scale,azi_scale
     1                     ,grid_spacing_m,r_missing_data)

          write(6,*)' Return from get_cloud_rays...',a9time

 900      continue

1000      continue

          ilun = ilun + 1
          write(clun,14)ilun
14        format(i3.3)

!         Get alt_a_roll and azi_a_roll arrays (needs to be generalized)?
          do i = minalt,maxalt
            call get_val(i,minalt,alt_scale,altobj)
            alt_a_roll(i,:) = altobj
          enddo 
          do j = minazi,maxazi
            call get_val(j,minazi,azi_scale,aziobj)
            azi_a_roll(:,j) = aziobj
          enddo

          blog_moon_roll = 0.
          if(moon_mag .lt. moon_mag_thr .AND.
     1     alm      .gt. 0.                 )then
            write(6,*)' Moon glow being calculated: ',alm,azm
            call get_glow_obj(i4time,alt_a_roll,azi_a_roll
     1                       ,minalt,maxalt,minazi,maxazi 
     1                       ,rlat,rlon,alt_scale,azi_scale
     1                       ,alm,azm,moon_mag
     1                       ,blog_moon_roll)

            write(6,*)' range of blog_moon_roll is',
     1          minval(blog_moon_roll),maxval(blog_moon_roll)
          endif

!         if(solar_alt .ge. 0.)then
          if(.true.)then
            I4_elapsed = ishow_timer()
            write(6,*)' call get_skyglow_cyl for sun or moon...'

            if(solar_alt .ge. -16.)then
                write(6,*)' Sun is significant, alt is:',solar_alt
                call skyglow_cyl(solar_alt,solar_az,blog_v_roll
     1                          ,elong_roll,od_atm_a
     1                          ,minalt,maxalt,minazi,maxazi
     1                          ,alt_scale,azi_scale)
                I4_elapsed = ishow_timer()

                if(moon_mag .lt. moon_mag_thr .AND.
     1             alm      .gt. 0.                 )then
                    if(solar_alt .ge. 0.)then ! valid only in daytime
                        write(6,*)' Moonglow is being added to daylight'       
                        blog_v_roll = 
     1                              addlogs(blog_v_roll,blog_moon_roll)      
                    endif
                endif
            endif ! solar_alt .ge. -16.

            if(solar_alt .lt. -16. .AND. moon_mag .lt. moon_mag_thr
     1                             .AND. alm .gt. 0.         )then
                write(6,*)' Moon skyglow significant: mag ',moon_mag
                call skyglow_cyl(alm,azm,blog_v_roll,elong_roll,od_atm_a
     1                          ,minalt,maxalt,minazi,maxazi
     1                          ,alt_scale,azi_scale)
                blog_v_roll = blog_v_roll + (-26.7 - moon_mag) * 0.4
                write(6,*)' Range of blog_v_roll for moon is',
     1                    minval(blog_v_roll),maxval(blog_v_roll)
                I4_elapsed = ishow_timer()
            endif

          else
            blog_v_roll = 8.0

          endif

!         Reproject Polar Cloud Plot
          lunsky = 60 
          write(lunsky,*)rmaglim_v
          call cyl_to_polar(r_cloud_3d,r_cloud_3d_polar,minalt,maxalt
     1                   ,maxazi,alt_scale,azi_scale
     1                   ,alt_a_polar,azi_a_polar
     1                   ,ni_polar,nj_polar)
!         write(6,*)' cyl slice at 40alt ',r_cloud_3d(40,:)
!         write(6,*)' polar slice at 256 ',r_cloud_3d_polar(256,:)

!         Reproject Skyglow Field
          call cyl_to_polar(blog_v_roll,blog_v_roll_polar
     1                               ,minalt,maxalt,maxazi
     1                               ,alt_scale,azi_scale
     1                               ,alt_a_polar,azi_a_polar
     1                               ,ni_polar,nj_polar)

!         Write time label
          open(53,file='label.'//clun,status='unknown')
          write(53,*)a9time
          close(53)

!         Write lat/lon label
          open(54,file='label2.txt',status='unknown')
          write(54,54)soundlat,soundlon
          close(54)
 54       format(2f8.2)

          if(mode_polar .eq. 1)then
            continue
          endif ! l_polar

!         if((l_cyl .eqv. .true.) .OR. (mode_polar .eq. 2))then              
          if(.true.)then

!           Get all sky for cyl   
            ni_cyl = maxalt - minalt + 1
            nj_cyl = maxazi - minazi + 1

            I4_elapsed = ishow_timer()

            if(solar_alt .lt. 0.)then
              write(6,*)' call get_starglow with cyl data'
              call get_starglow(i4time_solar,alt_a_roll,azi_a_roll       ! I
     1                     ,minalt,maxalt,minazi,maxazi                  ! I
     1                     ,rlat,rlon,alt_scale,azi_scale                ! I
     1                     ,glow_stars)                                  ! O

              write(6,*)' range of glow_stars (before) is',
     1             minval(glow_stars),maxval(glow_stars)

              write(6,*)' Moonglow is being added to starlight'       
              glow_stars = addlogs(glow_stars,blog_moon_roll)

              write(6,*)' range of glow_stars (after) is',
     1             minval(glow_stars),maxval(glow_stars)
            endif

            I4_elapsed = ishow_timer()

            write(6,*)' call get_sky_rgb with cyl data'
            call get_sky_rgb(r_cloud_3d              ! cloud opacity
     1                    ,cloud_od                  ! cloud optical depth
     1                    ,r_cloud_trans             ! cloud solar transmittance
     1                    ,cloud_rad_c               ! cloud solar transmittance / color
     1                    ,clear_rad_c               ! clear sky illumination by sun     
     1                    ,blog_v_roll               ! skyglow
     1                    ,glow_stars                ! starglow
     1                    ,od_atm_a
     1                    ,airmass_2_cloud_3d      
     1                    ,airmass_2_topo_3d      
     1                    ,topo_swi,topo_albedo
     1                    ,aod_2_cloud,aod_2_topo
     1                    ,alt_a_roll,azi_a_roll       
     1                    ,elong_roll    
     1                    ,ni_cyl,nj_cyl  
     1                    ,solar_alt,solar_az        ! sun alt/az
     1                    ,alm,azm,moon_mag          ! moon alt/az/mag
     1                    ,sky_rgb_cyl)   

            if(l_cyl .eqv. .true.)then
!             Write all sky for cyl
              isky_rgb_cyl = sky_rgb_cyl   
              npts = 3*(maxalt-minalt+1)*(maxazi-minazi+1)
              write(6,*)' Write all sky cyl text file ',npts
              open(55,file='allsky_rgb_cyl.'//clun,status='unknown')
              write(55,*)isky_rgb_cyl           
              close(55)
            endif

            if(l_polar .eqv. .true.)then
!             Reproject sky_rgb array from cyl to polar    
              do iaz = 0,maxazi,20
               write(6,*)'iaz,cyl(maxalt/2,iaz)',iaz
     1                       ,sky_rgb_cyl(1,maxalt/2,iaz)
              enddo ! iaz

              write(6,*)' Call cyl_to_polar with sky rgb data'

              do ic = 0,nc-1
                call cyl_to_polar(sky_rgb_cyl(ic,:,:)
     1                           ,sky_rgb_polar(ic,:,:)
     1                           ,minalt,maxalt,maxazi
     1                           ,alt_scale,azi_scale
     1                           ,alt_a_polar,azi_a_polar
     1                           ,ni_polar,nj_polar)
              enddo ! ic

!             Write all sky for polar
              isky_rgb_polar = sky_rgb_polar
              npts = 3*ni_polar*nj_polar
              write(6,*)' Write all sky polar text file'
     1                  ,isky_rgb_polar(:,255,255),npts
              open(54,file='allsky_rgb_polar.'//clun,status='unknown')
              write(54,*)isky_rgb_polar
              close(54)
            endif

          endif ! mode_polar = 0 or 2

          write(6,*)' End of plot_allsky for iloc...',iloc
          write(6,*)

          I4_elapsed = ishow_timer()

        enddo ! iloc

        write(6,*)
        write(6,*)' End of plot_allsky...'
        write(6,*)

        return
        end
