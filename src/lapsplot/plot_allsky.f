
        subroutine plot_allsky(i4time_ref,lun,NX_L,NY_L,NZ_L
     1                          ,ni_polar,nj_polar,ipolar_sizeparm
     1                          ,density     
     1                          ,iplo,iphi,jplo,jphi
     1                          ,r_missing_data,laps_cycle_time
     1                          ,l_polar,l_cyl)       

        use mem_namelist, ONLY: max_snd_grid, max_snd_levels
     1                        , grid_spacing_m, aod, aero_scaleht
     1                        , earth_radius
        use ppm
        use mem_allsky

        include 'trigd.inc'

        addlogs(x,y) = log10(10.**x + 10.**y)
        horz_depf(htmsl,erad) = acosd(erad/(erad+htmsl))

!       real pres_3d(NX_L,NY_L,NZ_L)
        real field_3d(NX_L,NY_L,NZ_L)
!       real heights_3d(NX_L,NY_L,NZ_L)
!       real clwc_3d(NX_L,NY_L,NZ_L)
!       real cice_3d(NX_L,NY_L,NZ_L)
!       real rain_3d(NX_L,NY_L,NZ_L)
!       real snow_3d(NX_L,NY_L,NZ_L)
!       real aod_3d(NX_L,NY_L,NZ_L)

        parameter (nc = 3)
!       real transm_3d(NX_L,NY_L,NZ_L)
!       real transm_4d(NX_L,NY_L,NZ_L,nc) 

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
!       real ghi_2d(NX_L,NY_L)
!       real dhi_2d(NX_L,NY_L)
        real static_albedo(NX_L,NY_L)
        real land_use(NX_L,NY_L)
        real land_frac(NX_L,NY_L)
        real snow_cover(NX_L,NY_L)
        real topo_albedo_2d(nc,NX_L,NY_L)
        real albedo_bm(nc,NX_L,NY_L)
        integer ialbedo_bm(nc,NX_L,NY_L)
        real albedo_usgs(nc,NX_L,NY_L)
        real lat(NX_L,NY_L)
        real lon(NX_L,NY_L)
        real topo(NX_L,NY_L)
        real dx(NX_L,NY_L)
        real dy(NX_L,NY_L)
        real alt_norm(NX_L,NY_L)     ! Solar Alt w.r.t. terrain normal
        real sol_alt_2d(NX_L,NY_L)
        real sol_azi_2d(NX_L,NY_L)
        real eobsc(NX_L,NY_L)        ! array of 'eobsl' values
        real moon_alt_2d(NX_L,NY_L)
        real moon_azi_2d(NX_L,NY_L)
        real moon_mag,moon_mag_thr

        real pres_1d(NZ_L)

        real lil_sfc, lic_sfc, lil_cpt, lic_cpt
  
        real k_to_c, make_td, make_ssh

        character*1 c_prodtype, c_plotobs
        character*3 var_2d
        character*150  directory, filename
        character*31  ext
        character*10  units_2d
        character*125 comment_2d
        character*40 c_model
        character*9 a9time
        character*24 a24time
        character*5 fcst_hhmm
        character*3 c3_string
        character*4 c4_string
        character*40 c_label
        character*11 c_pw
        character*11 c_cape
        character*20 c20_x, c20_y
        character*255 new_dataroot
        logical l_latlon, l_parse, l_plotobs, l_solar_eclipse
        logical l_idl /.false./
        logical l_cyl, l_polar, l_water_world 
        logical l_binary /.false./
        logical l_terrain_following /.false./
        logical l_require_all_fields ! requiring all LAPS fields to run
        logical l_test_cloud /.false./

        include 'icolors.inc'

!       Sounding observation declarations
        real lat_pr(max_snd_grid)
        real lon_pr(max_snd_grid)
        real elev_pr(max_snd_grid)
        integer nlevels_obs_pr(max_snd_grid)
        character*5 c5_name, c5_name_a(max_snd_grid), c5_name_min
        character*8 obstype(max_snd_grid)

        parameter (nsp = 4)

        real, allocatable, dimension(:,:) :: alt_a_roll
        real, allocatable, dimension(:,:) :: azi_a_roll
        real, allocatable, dimension(:,:) :: cloud_od
        integer, allocatable, dimension(:,:) :: camera_cloud_mask
        real*8, allocatable, dimension(:,:) :: dist_2_topo

        real alt_a_polar(iplo:iphi,jplo:jphi)
        real azi_a_polar(iplo:iphi,jplo:jphi)
        real elong_a_polar(iplo:iphi,jplo:jphi)

        real, allocatable, dimension(:,:,:) :: sky_rgb_polar
        integer, allocatable, dimension(:,:,:) :: isky_rgb_polar
!       real sky_rgb_polar(0:2,ni_polar,nj_polar)
!       integer isky_rgb_polar(0:2,ni_polar,nj_polar)

        real, allocatable, dimension(:,:,:) :: sky_rgb_cyl
        integer, allocatable, dimension(:,:,:) :: isky_rgb_cyl
        integer mode_cloud_mask /1/ ! ignore the mask

        integer maxloc
        parameter (maxloc = 1000)

        real xsound(maxloc),ysound(maxloc)
        real soundlat(maxloc),soundlon(maxloc)
        real htagl(maxloc),corr1_a(maxloc)

        data ilun /0/
        character*3 clun
 
        common /image/ n_image

        rpd = 3.141592653589/180.

        call alloc_allsky(NX_L,NY_L,NZ_L,nc,istatus)

        if(l_polar .eqv. .true.)then
            mode_polar = 2
        else
            mode_polar = 0
        endif

        nsmooth = 1

        I4_elapsed = ishow_timer()

        write(6,*)
        write(6,*)' subroutine plot_allsky: nsmooth/aod is ',nsmooth,aod
        write(6,*)' l_cyl/l_polar = ',l_cyl,l_polar
        write(6,*)' ipolar_sizeparm = ',ipolar_sizeparm
        write(6,*)' density = ',density

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
        read(lun,*)nloc

        write(6,*)' Number of all-sky locations from file is ',nloc

        do iloc = 1,nloc

 80       write(6,*)
          write(6,*)' Input x grid point (or latitude) for all-sky...'
          read(5,*)c20_x

          call s_len(c20_x,lenx)
          l_latlon = l_parse(c20_x(1:lenx),'l')

          if(l_latlon)then ! x value was flagged as latitude with "l" at the end 
            read(c20_x(1:lenx-1),*)soundlat(iloc)

            write(6,*)' Input longitude for allsky plot...'       
            read(5,*)c20_y
            call s_len(c20_y,leny)
            if(l_parse(c20_y(1:leny),'l'))then
                read(c20_y(1:leny-1),*)soundlon(iloc)
            else
                read(c20_y(1:leny),*)soundlon(iloc)
            endif

            call latlon_to_rlapsgrid(soundlat(iloc),soundlon(iloc)
     1                              ,lat,lon,NX_L,NY_L
     1                              ,xsound(iloc),ysound(iloc),istatus)

            if(istatus .ne. 1)then
                write(6,*)' Station is outside domain - try again...'
                return
            endif

            if(xsound(iloc) .lt. 1. .or. xsound(iloc) .gt. float(NX_L)
     1                              .or.
     1         ysound(iloc) .lt. 1. .or. ysound(iloc) .gt. float(NY_L)    
     1                                                             )then
                write(6,*)' Station is outside domain - try again...'
                return
            endif

          else
            read(c20_x(1:lenx),*)xsound(iloc)
            write(6,*)' Input y grid point for allsky plot...'
            read(5,*)ysound
            if(xsound(iloc) .lt. 1.0 .OR. ysound(iloc) .lt. 1.0)then ! scale domain
                write(6,*)' Values less than 1.0, scale to domain'
                xsound(iloc) = nint(xsound(iloc) * (NX_L-1)) + 1
                ysound(iloc) = nint(ysound(iloc) * (NY_L-1)) + 1
                soundlat(iloc) = 
     1              lat(nint(xsound(iloc)),nint(ysound(iloc)))
                soundlon(iloc) = 
     1              lon(nint(xsound(iloc)),nint(ysound(iloc)))
            endif
          endif

          read(5,*)htagl(iloc)

          write(6,*)' soundlat/soundlon ',soundlat(iloc),soundlon(iloc)
          write(6,*)' xsound/ysound ',xsound(iloc),ysound(iloc)
          write(6,*)' htagl ',htagl(iloc)

        enddo ! iloc

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
            read(5,15)new_dataroot
 15         format(a)
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

        write(6,*)' time diff is ',i4time_ref - i4time_now_gg()

        l_water_world = .false.
        if(l_parse(directory,'fim'))then
          l_require_all_fields = .false.
        elseif(l_parse(directory,'rams'))then
          l_require_all_fields = .false.
        elseif(i4time_ref - i4time_now_gg() .lt. 1e6)then  ! present/past
          l_require_all_fields = .true.
        elseif(i4time_ref - i4time_now_gg() .gt. 75e6)then ! >2.5y future 
          l_require_all_fields = .false.                   ! water only
          l_test_cloud = .false.
          l_water_world = .true.
        elseif(i4time_ref - i4time_now_gg() .gt. 45e6)then ! >1.5y future
          l_require_all_fields = .false.                   ! test clouds
          l_test_cloud = .true.
        else                                               ! >10d future
          l_require_all_fields = .false.
        endif

        if(l_require_all_fields .eqv. .true.)then

          n_lvls_snd = NZ_L

          write(6,*)' requiring fields at i4time_ref = ',i4time_ref

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

          write(6,*)' a9time (product info) is ',a9time

          i4_wdw = 108000

!         Read Height
          if(c_prodtype .eq. 'A')then
            iflag_temp = 2 ! Returns Height

            write(6,*)' Calling get_laps_3dgrid for heights'
            var_2d = 'HT'
            ext = 'lt1'
            call get_laps_3dgrid
     1          (i4time_ref,i4_wdw,i4time_nearest,NX_L,NY_L,NZ_L       
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

!         Read RH/SH (removed this section)
300       continue 

400       continue

          if(.not. l_test_cloud)then

!          Read Cloud Liquid
           istat_lwc = 0
           if(c_prodtype .eq. 'A')then ! Read Cloud Liquid
            var_2d = 'LWC'
            ext = 'lwc'
            call get_laps_3dgrid
     1          (i4time_ref,i4_wdw,i4time_lwc,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,clwc_3d,istat_lwc)
!           call make_fnam_lp(i4time_lwc,a9time,istatus)
!           call cv_i4tim_asc_lp(i4time_lwc,a24time,istatus)
            i4time_data  = i4time_lwc
           elseif(c_prodtype .eq. 'F')then 
            var_2d = 'LWC'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,clwc_3d
     1                              ,istat_lwc)
            call cv_i4tim_asc_lp(i4_valid,a24time,istatus)
            i4time_lwc  = i4_valid
            i4time_data = i4_valid
           endif

           if(istat_lwc .ne. 1)then
              write(6,*)' Error reading LWC field in plot_allsky'
              return
           endif

           if(istat_lwc .eq. 1)then
            continue
           else
            continue
           endif

!          Read Cloud Ice
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
            continue
           else
            write(6,*)' Error reading ICE field in plot_allsky'
            return
           endif
 
!          goto500

!          Read Precipitating Rain
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
            continue
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
            continue
           else
            write(6,*)' Error reading SNO field in plot_allsky'
            return
           endif

           i4wdw_sfc = 0

           goto500

!          Read Precipitating Ice
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
            continue
           else
            write(6,*)' Error reading PIC field in plot_allsky'
            return
           endif

          else ! l_test_cloud
           write(6,*)' generate idealized cloud fields'
           clwc_3d(:,:,14) = .000
           i4time_lwc = i4time_ref
           i4wdw_sfc = 100000

          endif ! l_test_cloud is F

500       continue

          if(c_prodtype .eq. 'A')then 
!           Read in swi data
            ext = 'lcv'
            var_2d = 'SWI'
            call get_laps_2dgrid(i4time_lwc,i4wdw_sfc,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,swi_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)then
              write(6,*)' Error reading LCV/SWI field in plot_allsky'      
              return
            endif

!           Read in snow cover data
            ext = 'lm2'
            var_2d = 'SC'
            call get_laps_2dgrid(i4time_lwc,i4wdw_sfc,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,snow_cover,0,istat_sfc)
            if(istat_sfc .ne. 1 .and. istat_sfc .ne. -1)then
                write(6,*)' Error reading LM2/SC field in plot_allsky'      
                return
            endif
!           snow_cover = 0. ! testing
            write(6,*)' range of snow_cover is',
     1                minval(snow_cover),maxval(snow_cover)

!           Read in pw data
            ext = 'lh4'
            var_2d = 'TPW'
            call get_laps_2dgrid(i4time_lwc,i4wdw_sfc,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,pw_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)then
                write(6,*)' Error reading LH4/PW field in plot_allsky'      
                return
            endif
            write(6,*)' range of pw_2d is',
     1                minval(pw_2d),maxval(pw_2d)

          else
            snow_cover = r_missing_data
            pw_2d = r_missing_data
          endif

          write(6,*)
     1          ' calling get_static_field_interp for static_albedo'
          call get_static_field_interp('albedo',i4time_lwc,NX_L,NY_L
     1                                ,static_albedo,istat_sfc)
          if(istat_sfc .ne. 1)then
              write(6,*)' Error in get_static_field_interp'      
              return
          else
              write(6,*)' success from get_static_field_interp'
          endif

          if(.true.)then ! initial setting of albedo (just a placeholder)
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

          I4_elapsed = ishow_timer()

          if(.true.)then ! experimental
            i4time_solar = i4time_ref
          else
            i4time_solar = i4time_data
          endif

!         Calculate solar position for 2D array of grid points
          write(6,*)' call solar_position for 2D array'
          do i = 1,NX_L
          do j = 1,NY_L
            call solar_position(lat(i,j),lon(i,j),i4time_solar
     1                         ,sol_alt_2d(i,j),solar_dec,solar_ha)
            call equ_to_altaz_d(solar_dec,solar_ha,lat(i,j)
     1                         ,altdum,sol_azi_2d(i,j))               
            if(sol_azi_2d(i,j) .lt. 0.)sol_azi_2d(i,j) = 
     1                                 sol_azi_2d(i,j) + 360.
          enddo ! j
          enddo ! i

        else ! l_require_all_fields = F
          if(l_test_cloud)then
            write(6,*)' generate idealized cloud fields'
            clwc_3d(:,:,14) = .000 ! .001
          endif

          if(l_parse(directory,'fim'))then
            filename = '/scratch/staging/fab/albers/fimlarge.nc'
            write(6,*)' Looking for FIM LWC data in ',trim(filename)
            call get_fim_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,pres_3d
     +                   ,clwc_3d
     +                   ,lun_out
     +                   ,istatus)
            istatus_ht = 0
            write(6,*)' returned from get_fim_data ',istatus
            
          elseif(l_parse(directory,'rams'))then
            write(6,*)' Looking for RAMS LWC data in ',trim(directory)
            filename = '/home/fab/albers/muri/rams_micro_v3.nc'
            call get_rams_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,pres_3d
     +                   ,heights_3d
     +                   ,clwc_3d
     +                   ,cice_3d
     +                   ,rain_3d
     +                   ,snow_3d
     +                   ,lun_out
     +                   ,istatus)
            istatus_ht = 1
            write(6,*)' returned from get_rams_data'
          else
            istatus_ht = 0
          endif

          i4time_solar = i4time_ref

!         Calculate solar position for 2D array of grid points
          do i = 1,NX_L
          do j = 1,NY_L
            call solar_position(lat(i,j),lon(i,j),i4time_solar
     1                         ,sol_alt_2d(i,j),solar_dec,solar_ha)
            call equ_to_altaz_d(solar_dec,solar_ha,lat(i,j)
     1                         ,altdum,sol_azi_2d(i,j))               
            if(sol_azi_2d(i,j) .lt. 0.)sol_azi_2d(i,j) = 
     1                                 sol_azi_2d(i,j) + 360.
            if(i .eq. 71 .and. j .eq. 275)then
              write(6,*)'i/j/lat/lon/solar_dec/solar_ha/sol_alt_2d(i,j)'
     1                  ,i,j,lat(i,j),lon(i,j),solar_dec,solar_ha
     1                  ,sol_alt_2d(i,j)
            endif
   
            swi_2d(i,j) = 1300. * sind(max(sol_alt_2d(i,j),0.))

          enddo ! j
          enddo ! i

          read(lun,*) ! advance through input data
          write(6,*)' Running without LAPS cloud and other current data'
          snow_cover = 0. ! r_missing_data

!         Use standard atmosphere for heights (uniform pressure grid)
          if(istatus_ht .eq. 0)then
            write(6,*)' Getting 1D pressure levels'
            call get_pres_1d(i4time_ref,NZ_L,pres_1d,istatus)
            if(istatus .ne. 1)then
              write(6,*)' error getting 1d pressures'
              return
            endif
            write(6,*)' Getting 3D heights from 1D pressure levels'
            do k = 1,NZ_L
              heights_3d(:,:,k) = psatoz(pres_1d(k)/100.)
            enddo ! k
          endif

        endif ! l_require_all_fields is TRUE

        call make_fnam_lp(i4time_solar,a9time,istatus)
        call cv_i4tim_asc_lp(i4time_solar,a24time,istatus)

        I4_elapsed = ishow_timer()

!       Consider additional albedo info based on land use 
!       This also includes snow cover
        var_2d='USE'
        write(6,*)' calling read_static_grid for land use'
        call read_static_grid(NX_L,NY_L,var_2d,land_use,istatus)
        if(istatus .ne. 1)then
           print*,' Warning: could not read static-landuse'
           return
        else
           write(6,*)' Successful return from read_static_grid'
        endif

        var_2d='LDF'
        write(6,*)' calling read_static_grid for land frac'
        call read_static_grid(NX_L,NY_L,var_2d,land_frac,istatus)
        if(istatus .ne. 1)then
           print*,' Warning: could not read static-landfrac'
           return
        else
           write(6,*)' Successful return from read_static_grid'
        endif

        I4_elapsed = ishow_timer()

        where(topo(:,:) .ge. 3200.); land_use(:,:) = 19.; end where

        call land_albedo(land_use,NX_L,NY_L,albedo_usgs)
        write(6,*)' debias albedo_usgs by a factor of 1.25'
        albedo_usgs = albedo_usgs / 1.25

        call land_albedo_bm(lat,lon,NX_L,NY_L,i4time_solar
     1                     ,albedo_bm,istat_bm)

        I4_elapsed = ishow_timer()
 
        if(istat_bm .eq. 1)then
            write(6,*)' Use 3-color albedo based on Blue Marble'

!           Test a flop if needed
!           do i = 1,NX_L
!           do j = 1,NY_L
!               topo_albedo_2d(:,NX_L+1-i,NY_L+1-j) 
!    1             = albedo_bm(:,i,j)
!               topo_albedo_2d(:,i,j) 
!    1             = albedo_bm(:,i,j)
!           enddo ! j
!           enddo ! i

            topo_albedo_2d = albedo_bm

            write(6,*)
     1          ' Row of multi-spectral albedo (through domain center)'
            jrow = NY_L/2
            do i = 1,NX_L
                if(i .eq. (i/5)*5 .OR. abs(i-NX_L/2) .lt. 20)then
                    write(6,16)i,lat(i,jrow),lon(i,jrow)
     1                       ,land_frac(i,jrow),topo_albedo_2d(:,i,jrow)
16                  format(i5,3f9.3,2x,3f9.3)                 
                endif
            enddo ! i
            alb_min = minval(topo_albedo_2d(2,:,jrow))
            alb_max = maxval(topo_albedo_2d(2,:,jrow))
            write(6,*)' Min/Max in row is ',alb_min,alb_max

            alb_max = maxval(topo_albedo_2d(2,:,:))
            do i = 1,NX_L
            do j = 1,NY_L
                if(topo_albedo_2d(2,i,j) .eq. alb_max)then
                  if(alb_max .le. 1.0)then
                    write(6,*)' Max albedo at ',i,j
     1                        ,topo_albedo_2d(:,i,j)
                  endif
                endif
            enddo ! j
            enddo ! i

            if(.false.)then
                write(6,*)' Write test blue marble albedo image'
                call get_directory('static',directory,len_dir)
                ialbedo_bm(:,:,:) = nint(topo_albedo_2d(:,:,:)*255.)
                call writeppm3Matrix(
     1               ialbedo_bm(1,:,:),ialbedo_bm(2,:,:)
     1              ,ialbedo_bm(3,:,:),trim(directory)//'/ialbedo_bm')
            endif
        else
            write(6,*)' Use 3-color albedo based on land use'
            topo_albedo_2d = albedo_usgs
        endif

        call compare_land_albedo(land_use,NX_L,NY_L,albedo_usgs
     1                          ,albedo_bm,static_albedo)
        
        write(6,*)' Row of albedo_bm / snow / snowalb / albedo / topo'
        jrow = NY_L/2
        do i = 1,NX_L
          do j = 1,NY_L
            if(snow_cover(i,j) .ne. r_missing_data)then
              if(topo(i,j) .gt. 1900. .and. topo(i,j) .le. 3500.)then
                snowalb = snow_cover(i,j)**2. * 0.7
              else
                snowalb = snow_cover(i,j) * 0.7
              endif
              do ic = 1,3
!               topo_albedo_2d(ic,i,j) = 
!    1            max(topo_albedo_2d(ic,i,j),snowalb)
                topo_albedo_2d(ic,i,j) = snow_cover(i,j) * snowalb + 
     1                (1.0-snow_cover(i,j)) * topo_albedo_2d(ic,i,j)
              enddo ! ic
            endif
            if(j .eq. jrow)then
              if(i .eq. (i/5)*5 .OR. abs(i-NX_L/2) .lt. 20)then
                write(6,18)i,albedo_bm(2,i,jrow),snow_cover(i,jrow)
     1                ,snowalb,topo_albedo_2d(2,i,jrow),topo(i,jrow)
     1                ,lat(i,jrow),lon(i,jrow),nint(land_use(i,jrow))
18              format(i5,2x,4f9.3,f9.0,4x,2f9.2,i4)
              endif
            endif
          enddo ! j
        enddo ! i 

        I4_elapsed = ishow_timer()

        call get_grid_spacing_array(lat,lon,NX_L,NY_L,dx,dy)

        I4_elapsed = ishow_timer()

        call solar_normal(NX_L,NY_L,topo,dx,dy,lat,lon ! I
     1                   ,sol_alt_2d,sol_azi_2d        ! I
     1                   ,alt_norm)                    ! O

        I4_elapsed = ishow_timer()

        write(6,*)' a9time (solar) is ',a9time

!       Determine aod_ref as aerosol optical depth
        pw_ref = pw_2d(NX_L/2,NY_L/2)
        if(pw_ref .eq. r_missing_data)then
           write(6,*)' pw_ref is missing, use default value'
           aod_ref = .07
        elseif(aod .lt. 0.)then
           write(6,*)' aod is < 0, use as scale factor with pw_ref'
           aod_ref = pw_ref * (-aod)
        else
           write(6,*)' aod is >= 0, use directly for aod_ref'
           aod_ref = aod
        endif

        write(6,19,err=20)pw_ref,aod,aod_ref
19      format(' pw_ref,aod,aod_ref = ',3f10.3)
20      continue

        if(l_parse(directory,'rams'))then
           mode_aero_cld = 3
           write(6,*)' RAMS run: set mode_aero_cld = ',mode_aero_cld

!          Zero out hydrometeors (e.g. from RAMS)
           if(.true.)then
              write(6,*)' RAMS run: zero out hydrometeors'
              clwc_3d = 0.
              cice_3d = 0.
              rain_3d = 0.
              snow_3d = 0.
           endif
        endif

!       Get Aersol Extinction Coefficient (3D field)
        call get_aod_3d(pres_3d,heights_3d,topo,NX_L,NY_L,NZ_L,aod_ref
     1                 ,aod_3d)

        I4_elapsed = ishow_timer()

        rlat_last = -999.
        rlon_last = -999.
        i4time_last = 0.
      
        do iloc = 1,nloc
          i_obs = nint(xsound(iloc))
          j_obs = nint(ysound(iloc))

          ri_obs = xsound(iloc)
          rj_obs = ysound(iloc)

          rlat = lat(i_obs,j_obs)
          rlon = lon(i_obs,j_obs)

          rlat = soundlat(iloc)
          rlon = soundlon(iloc)

          write(6,*)
          write(6,*)' iloc = ',iloc,rlat,rlon

          write(6,*)' Enter minalt,maxalt (e.g. 0,90 0,180)'
          read(lun,*)minalt,maxalt           

          write(6,*)' Enter minazi,maxazi'                     
          read(lun,*)minazi,maxazi           

          write(6,*)' Enter alt_scale,azi_scale'               
          read(lun,*)alt_scale,azi_scale     

          write(6,*)' minalt/maxalt = ',minalt,maxalt
          write(6,*)' minazi/maxazi = ',minazi,maxazi
          write(6,*)' alt_scale/azi_scale = ',alt_scale,azi_scale

          if(minazi .eq. maxazi)then
            write(6,*)' Error minazi = maxazi'
            stop
          endif

          allocate(alt_a_roll(minalt:maxalt,minazi:maxazi))
          allocate(azi_a_roll(minalt:maxalt,minazi:maxazi))
          allocate(cloud_od(minalt:maxalt,minazi:maxazi))
          allocate(camera_cloud_mask(minalt:maxalt,minazi:maxazi))
          allocate(dist_2_topo(minalt:maxalt,minazi:maxazi))
          allocate(sky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi))
          allocate(isky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi))

          if(.false.)then
              topo_sfc = topo(i_obs,j_obs)
          else
              call bilinear_laps(ri_obs,rj_obs,NX_L,NY_L,topo,topo_sfc)
          endif
 
          write(6,*)' ri/rj/topo_sfc ',ri_obs,rj_obs,topo_sfc
          write(6,22)topo_albedo_2d(:,i_obs,j_obs)
22        format('  albedo RGB of observer ',3f9.3)

          write(6,*)' snow cover of observer: ',snow_cover(i_obs,j_obs)  

          write(6,*)' solar alt/az (2d array)',sol_alt_2d(i_obs,j_obs)       
     1                                        ,sol_azi_2d(i_obs,j_obs) 

!         Calculate solar position for all-sky point
          if(.true.)then ! test this again to allow fractional gridpoints?
              call solar_position(soundlat(iloc),soundlon(iloc)
     1                           ,i4time_solar,solar_alt     
     1                           ,solar_dec,solar_ha)
              call equ_to_altaz_d(solar_dec,solar_ha,soundlat(iloc)
     1                           ,altdum,solar_az)               
              if(solar_az .lt. 0.)solar_az = solar_az + 360.
              solar_lat = solar_dec
              solar_lon = soundlon(iloc) - solar_ha
          else ! ensure consistency between both solar positions
              solar_alt = sol_alt_2d(i_obs,j_obs)
              solar_az = sol_azi_2d(i_obs,j_obs)              
          endif

          write(6,*)' solar alt/az (observer)',solar_alt,solar_az

          hdist_loc = sqrt((rlat-rlat_last)**2 + (rlon-rlon_last)**2)
!         if(solar_alt .gt. 4.0)then
!             hdist_loc_thr = 0.3
!         else
!             hdist_loc_thr = 0.0
!         endif

          thr1 = ((grid_spacing_m / 10000.) / rpd) * sind(solar_alt) 
          thr2 = solar_alt * 0.1
          thr3 = solar_alt - 3.0
          hdist_loc_thr = max(min(thr1,thr2,thr3),0.0)

          if(hdist_loc .gt. hdist_loc_thr .OR. 
     1       i4time_solar .ne. i4time_last    )then
            newloc = 1
            rlat_last = rlat; rlon_last = rlon
            i4time_last = i4time_solar
          else
            newloc = 0
          endif

          write(6,*)' i4time_last/i4time_solar = '
     1               ,i4time_last,i4time_solar
          write(6,23)hdist_loc,hdist_loc_thr,newloc
23        format('  hdist_loc/hdist_loc_thr/newloc = ',2f9.3,i3)

          write(6,*)' call sun_moon at observer sfc grid point ',i_obs
     1                                                          ,j_obs
          idebug = 2
          call sun_moon(i4time_solar,lat,lon,NX_L,NY_L,i_obs,j_obs   ! I
     1                 ,alm,azm                                      ! O
     1                 ,idebug,0.,earth_radius                       ! I
     1                 ,elgms,moon_mag,rmn                           ! O 
     1                 ,geo_dec,geo_ra,geo_sublon,geo_dist           ! O
     1                 ,emag,eobsf,eobsl)                            ! O

          if(elgms .lt. 1.4)then
            write(6,24)emag,eobsf,eobsl
 24         format(' NOTE: Solar Eclipse Conditions: mag/obsc = '
     1            ,f9.6,2f11.6)
            l_solar_eclipse = .true.
          elseif(elgms .lt. 0.6)then
            write(6,*)' NOTE: Possible Solar Eclipse Conditions'
            l_solar_eclipse = .true.
          else
            l_solar_eclipse = .false.
          endif

          write(6,*)' call sun_moon at observer sfc grid point ',i_obs
     1                                                          ,j_obs
          call sun_moon(i4time_solar,lat,lon,NX_L,NY_L,i_obs,j_obs   ! I
     1                 ,alm,azm                                      ! O
     1                 ,idebug,htagl(iloc),earth_radius              ! I
     1                 ,elgms,moon_mag,rmn                           ! O 
     1                 ,geo_dec,geo_ra,geo_sublon,geo_dist           ! O
     1                 ,emag,eobsf,eobsl)                            ! O
          write(6,25)alm,azm,elgms,moon_mag,rmn
 25       format('  alt/az/elg/mnmag/rmn = ',2f8.2,f9.4,f8.2,f9.6)

!         Consider passing 'topo_flag' into 'sun_moon' to consider either
!         solar or lunar eclipses
!         http://www.jgisen.de/eclipse

!         'emag' is solar eclipse magnitude for the observer
!         'eobsf' is observer obscuration
!         'eobsl' includes limb darkening for observer obscuration

          if(maxval(sol_alt_2d) .gt. 0. .and. 
     1       l_solar_eclipse .eqv. .true.)then
            write(6,*)' Calculate gridded sfc eclipse obscuration'
            elgmin = 9999.
            idebug = 0
            do i = 1,NX_L
            do j = 1,NY_L
              call sun_moon(i4time_solar,lat,lon,NX_L,NY_L,i,j         ! I
     1                     ,almgrd,azmgrd                              ! O
     1                     ,idebug,0.,earth_radius                     ! I
     1                     ,elggrd,grdmoon_mag,rmn                     ! O
     1                     ,geo_dec,geo_ra,geo_sublon,geo_dist         ! O
     1                     ,solar_eclipse_magnitude,eobsf,eobsc(i,j))  ! O

              if(elggrd .lt. elgmin .and. sol_alt_2d(i,j) .gt. 0.)then
                elgmin = elggrd
                emagmax = solar_eclipse_magnitude
                eobsmax = eobsc(i,j)
                imxecl = i
                jmxecl = j
              endif

            enddo ! j
            enddo ! i

            write(6,*)' range of eobsc is',minval(eobsc),maxval(eobsc)
            write(6,*)' best elg/mag/eobs ',elgmin,emagmax,eobsmax,
     1                ' at',i,j,lat(imxecl,jmxecl),lon(imxecl,jmxecl)

          else
            eobsc = 0.
          endif

!         l_solar_eclipse = .false. ! test

!         alm = -90.          ! Test for disabling
!         moon_mag = -4.0    ! Test for disabling
          moon_mag_thr = -6.0

          moon_alt_2d = alm
          moon_azi_2d = azm

!         Get alt_a_roll and azi_a_roll arrays
          do i = minalt,maxalt
            call get_val(i,minalt,alt_scale,altobj)
            alt_a_roll(i,:) = altobj
          enddo 
          do j = minazi,maxazi
            call get_val(j,minazi,azi_scale,aziobj)
            azi_a_roll(:,j) = aziobj
          enddo

          write(6,*)' alt range is ',alt_a_roll(minalt,minazi)
     1                              ,alt_a_roll(maxalt,minazi)
          write(6,*)' azi range is ',azi_a_roll(minalt,minazi)
     1                              ,azi_a_roll(minalt,maxazi)

          ilun = ilun + 1
          write(clun,34)ilun
34        format(i3.3)

          allocate(aod_ill_opac(minalt:maxalt,minazi:maxazi))
          allocate(aod_ill_opac_potl(minalt:maxalt,minazi:maxazi))

          exposure = density

          if(l_water_world)then ! water world experimental simulation
            land_frac = 0.
            snow_cover = 0.
            topo = 0.
            topo_sfc = 0.
            topo_albedo_2d(1,:,:) = .003
            topo_albedo_2d(2,:,:) = .007
            topo_albedo_2d(3,:,:) = .028
          endif

          if(.true.)then
            write(6,*)' call calc_allsky at i4time_solar:',i4time_solar
            call calc_allsky(i4time_solar,exposure ! ,clwc_3d,cice_3d
!    1                     ,heights_3d                              ! I
!    1                     ,rain_3d,snow_3d                         ! I
!    1                     ,pres_3d,aod_3d                          ! I
     1                     ,topo_sfc,topo,swi_2d                    ! I
     1                     ,topo_albedo_2d,land_frac,snow_cover     ! I
     1                     ,htagl(iloc)                             ! I
     1                     ,aod_ref                                 ! I
     1                     ,NX_L,NY_L,NZ_L,newloc                   ! I
     1                     ,ri_obs,rj_obs                           ! I
     1                     ,alt_a_roll,azi_a_roll                   ! I
     1                     ,sol_alt_2d,sol_azi_2d                   ! I
     1                     ,solar_alt,solar_az                      ! I
     1                     ,solar_lat,solar_lon                     ! I
     1                     ,alt_norm                                ! I
     1                     ,moon_alt_2d,moon_azi_2d,alm,azm         ! I
     1                     ,moon_mag,moon_mag_thr,elgms             ! I
     1                     ,l_solar_eclipse,eobsc,emag              ! I
     1                     ,rlat,rlon,lat,lon                       ! I
     1                     ,minalt,maxalt,minazi,maxazi,nc,nsp      ! I
     1                     ,ni_cyl,nj_cyl                           ! O
     1                     ,alt_scale,azi_scale                     ! I
     1                     ,grid_spacing_m,r_missing_data           ! I
     1                     ,l_binary,l_terrain_following            ! I
     1                     ,mode_cloud_mask,camera_cloud_mask       ! I
     1                     ,cloud_od,dist_2_topo                    ! O
     1                     ,sky_rgb_cyl,istatus)                    ! O
            if(istatus .ne. 1)then
              write(6,*)' Error istatus returned from calc_allsky'
              return
            endif

          else
            continue

          endif ! call calc_allsky

          deallocate(aod_ill_opac)
          deallocate(aod_ill_opac_potl)
          deallocate(alt_a_roll)
          deallocate(azi_a_roll)
          deallocate(cloud_od)
          deallocate(camera_cloud_mask)
          deallocate(dist_2_topo)

          write(6,*)' end of subroutine call block - write labels'

!           Write time label
            open(53,file='label.'//clun,status='unknown')
            write(53,*)a9time
            write(53,*)a24time(1:17)
            close(53)


!           Write lat/lon and other info for label
            open(54,file='label2.'//clun,status='unknown')
            write(54,54)soundlat(iloc),soundlon(iloc),
     1                  minalt,maxalt,minazi,maxazi,ni_cyl,nj_cyl,
     1                  solar_alt,solar_az,alt_scale,azi_scale,
     1                  ni_polar,nj_polar
 54         format(2f8.2/6i8/2f8.2,2f7.2/2i6)

            if(ghi_sim .eq. r_missing_data)then
                write(54,*)
            elseif(ghi_sim .gt. 100.)then
                write(54,55)nint(ghi_sim)
 55             format(i4)
            elseif(ghi_sim .gt. 1.)then
                write(54,56)ghi_sim
 56             format(f8.1)
            else
                write(54,57)ghi_sim
 57             format(f8.4)
            endif
                
            close(54)

            if(l_cyl .eqv. .true.)then
!             Write all sky for cyl
              isky_rgb_cyl = sky_rgb_cyl   
              npts = 3*(maxalt-minalt+1)*(maxazi-minazi+1)
!             write(6,*)' Write all sky cyl text file ',npts
!             open(55,file='allsky_rgb_cyl.'//clun,status='unknown')
!             write(55,*)isky_rgb_cyl           
!             close(55)
              write(6,*)' Write all sky cyl ppm file '
              call writeppm3Matrix(
     1                   isky_rgb_cyl(0,:,:),isky_rgb_cyl(1,:,:)
     1                  ,isky_rgb_cyl(2,:,:),'allsky_rgb_cyl_'//clun)
              do iaz = minazi,maxazi,40
               write(6,*)'iaz,cyl(maxalt/1,iaz)',iaz
     1                       ,isky_rgb_cyl(:,maxalt/1,iaz)
              enddo ! iaz
            endif

            if(l_polar .eqv. .true.)then

              allocate(sky_rgb_polar(0:2,iplo:iphi,jplo:jphi))
              allocate(isky_rgb_polar(0:2,iplo:iphi,jplo:jphi))

!             Reproject sky_rgb array from cyl to polar    
              do iaz = minazi,maxazi,20
               write(6,*)'iaz,cyl((maxalt+minalt)/2,iaz)',iaz
     1                       ,sky_rgb_cyl(1,(maxalt+minalt)/2,iaz)
              enddo ! iaz

              write(6,*)' Call cyl_to_polar with sky rgb data'

              if(htagl(iloc) .le. 18000. .and. 
     1                  (maxalt*2 + minalt) .gt. 0)then
                polat = +90. ! +/-90 for zenith or nadir at center of plot
              else
                polat = -90.
              endif

              horz_dep = horz_depf(htagl(iloc),earth_radius)

!             if(ipolar_sizeparm .ge. 3)then
              if(htagl(iloc) .gt. earth_radius*2.5)then
                pomag = htagl(iloc) / (earth_radius*0.7)
              elseif(htagl(iloc) .gt. earth_radius*1.1)then
                pomag = 3.0
              elseif(htagl(iloc) .gt. earth_radius*0.75)then
                pomag = 2.5
              elseif(htagl(iloc) .gt. earth_radius*0.5)then
                pomag = 2.0
              elseif(polat .eq. -90.)then ! looking down
                write(6,*)' horz_dep = ',horz_dep
                pomag = (90. / (90. - horz_dep)) * 0.96
                pomag = max(pomag,1.0)
              else                        ! looking up
                pomag = 1.0
              endif

              write(6,*)'htrat/pomag',htagl(iloc)/earth_radius,pomag

              do ic = 0,nc-1
                call cyl_to_polar(sky_rgb_cyl(ic,:,:)
     1                           ,sky_rgb_polar(ic,:,:)
     1                           ,minalt,maxalt,minazi,maxazi
     1                           ,alt_scale,azi_scale,polat,pomag
     1                           ,alt_a_polar,azi_a_polar
     1                           ,iplo,iphi,jplo,jphi
     1                           ,ni_polar,nj_polar)
              enddo ! ic

!             Write all sky for polar
              where(sky_rgb_polar .eq. r_missing_data)
                  sky_rgb_polar = 0.
                  isky_rgb_polar = 0
              endwhere
              isky_rgb_polar = sky_rgb_polar
              write(6,*)' max polar 1 is ',maxval(sky_rgb_polar)
     1                                    ,maxval(isky_rgb_polar)
              write(6,*)' min polar 1 is ',minval(sky_rgb_polar)
     1                                    ,minval(isky_rgb_polar)

              if(minval(isky_rgb_polar) .lt. 0)then
                write(6,*)' WARNING: maxval isky_rgb_polar < 0'
              
                where(isky_rgb_polar .lt. 0)
                      isky_rgb_polar = 0
                endwhere
                write(6,*)' max polar 2 is ',maxval(sky_rgb_polar)
     1                                      ,maxval(isky_rgb_polar)
                write(6,*)' min polar 2 is ',minval(sky_rgb_polar)
     1                                      ,minval(isky_rgb_polar)
              endif

!             write(6,*)' ipolar array',
!    1            isky_rgb_polar(2,ni_polar/2,1:nj_polar)
              nip_crop = iphi-iplo+1
              njp_crop = jphi-jplo+1
              npts = 3*nip_crop*njp_crop
!             write(6,*)' Write all sky polar text file'
!    1                  ,isky_rgb_polar(:,255,255),npts
!             open(54,file='allsky_rgb_polar.'//clun,status='unknown')
!             write(54,*)isky_rgb_polar
!             close(54)
              write(6,*)' Write all sky polar ppm file '
              call writeppm3Matrix(
     1                  isky_rgb_polar(0,:,:),isky_rgb_polar(1,:,:)
     1                 ,isky_rgb_polar(2,:,:),'allsky_rgb_polar_'//clun)

              deallocate(sky_rgb_polar)
              deallocate(isky_rgb_polar)

            endif ! l_polar

!         endif ! mode_polar = 0 or 2

 900      continue

1000      continue

          write(6,*)' End of plot_allsky for iloc...',iloc
          write(6,*)

          I4_elapsed = ishow_timer()

          deallocate(sky_rgb_cyl)
          deallocate(isky_rgb_cyl)

        enddo ! iloc

        call dealloc_allsky

        write(6,*)
        write(6,*)' End of plot_allsky...'
        write(6,*)

        return
        end
