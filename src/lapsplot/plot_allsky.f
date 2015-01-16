
        subroutine plot_allsky(i4time_ref,lun,NX_L,NY_L,NZ_L
!    1                          ,minalt,maxalt,minazi,maxazi
     1                          ,r_missing_data,laps_cycle_time,maxstns
     1                          ,i_overlay,plot_parms,namelist_parms
     1                          ,l_polar,l_cyl)       

        use mem_namelist, ONLY: max_snd_grid, max_snd_levels
     1                        , grid_spacing_m, aod, aero_scaleht
        use ppm

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
        real transm_3d(NX_L,NY_L,NZ_L)
        real transm_4d(NX_L,NY_L,NZ_L,nc) 

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
        real snow_cover(NX_L,NY_L)
        real topo_albedo_2d(nc,NX_L,NY_L)
        real albedo_bm(nc,NX_L,NY_L)
        integer ialbedo_bm(nc,NX_L,NY_L)
        real lat(NX_L,NY_L)
        real lon(NX_L,NY_L)
        real topo(NX_L,NY_L)
        real dx(NX_L,NY_L)
        real dy(NX_L,NY_L)
        real alt_norm(NX_L,NY_L)
        real sol_alt_2d(NX_L,NY_L)
        real sol_azi_2d(NX_L,NY_L)
        real eobsc(NX_L,NY_L)
        real moon_alt_2d(NX_L,NY_L)
        real moon_azi_2d(NX_L,NY_L)
        real moon_mag,moon_mag_thr

        real pres_1d(NZ_L)

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
        logical l_cyl               
        logical l_polar 
        logical l_allsky_clouds /.true./ ! requiring cloud data to run

        integer i_overlay

        include 'icolors.inc'

!       Sounding observation declarations
        real lat_pr(max_snd_grid)
        real lon_pr(max_snd_grid)
        real elev_pr(max_snd_grid)
        integer nlevels_obs_pr(max_snd_grid)
        character*5 c5_name, c5_name_a(max_snd_grid), c5_name_min
        character*8 obstype(max_snd_grid)

        parameter (nsp = 4)
        real, allocatable, dimension(:,:) :: r_cloud_3d ! cloud opacity        
        real, allocatable, dimension(:,:) :: cloud_od   ! cloud optical depth  
        real, allocatable, dimension(:,:,:) :: cloud_od_sp ! cloud species tau
        real, allocatable, dimension(:,:) :: blog_v_roll                         
        real, allocatable, dimension(:,:) :: blog_moon_roll                         
        real, allocatable, dimension(:,:) :: blog_sun_roll                         
        real, allocatable, dimension(:,:,:) :: glow_stars                        
        real, allocatable, dimension(:,:) :: elong_roll                        
        real, allocatable, dimension(:,:) :: airmass_2_cloud_3d
        real, allocatable, dimension(:,:) :: airmass_2_topo_3d                        
        real, allocatable, dimension(:,:) :: topo_swi                        
        real, allocatable, dimension(:,:,:) :: topo_albedo
        real, allocatable, dimension(:,:,:) :: ghic
        real, allocatable, dimension(:,:) :: aod_2_cloud
        real, allocatable, dimension(:,:) :: aod_2_topo
        real, allocatable, dimension(:,:) :: dist_2_topo
        real, allocatable, dimension(:,:) :: aod_ill
        real, allocatable, dimension(:,:) :: aod_ill_dir
        real, allocatable, dimension(:,:) :: aod_tot
        real, allocatable, dimension(:,:) :: r_cloud_trans ! sun to cloud transmissivity (direct+fwd scat)
        real, allocatable, dimension(:,:,:) :: cloud_rad_c ! sun to cloud transmissivity (direct+fwd scat) * solar color/int
        real, allocatable, dimension(:,:) :: cloud_rad_w   ! sun to cloud transmissivity (direct+fwd scat) * rad
        real, allocatable, dimension(:,:,:) :: clear_rad_c ! clear sky illumination
        real, allocatable, dimension(:,:,:) :: clear_radf_c! integrated fraction of air illuminated by the sun along line of sight
        real, allocatable, dimension(:,:) :: alt_a_roll
        real, allocatable, dimension(:,:) :: azi_a_roll

        parameter (ni_polar = 511)
        parameter (nj_polar = 511)
        real r_cloud_3d_polar(ni_polar,nj_polar)
        real blog_v_roll_polar(ni_polar,nj_polar)
        real alt_a_polar(ni_polar,nj_polar)
        real azi_a_polar(ni_polar,nj_polar)
        real elong_a_polar(ni_polar,nj_polar)

        real sky_rgb_polar(0:2,ni_polar,nj_polar)
        integer isky_rgb_polar(0:2,ni_polar,nj_polar)

        real, allocatable, dimension(:,:,:) :: sky_rgb_cyl
        integer, allocatable, dimension(:,:,:) :: isky_rgb_cyl

        integer maxloc
        parameter (maxloc = 10)

        real xsound(maxloc),ysound(maxloc)
        real soundlat(maxloc),soundlon(maxloc)
        real htagl(maxloc)

        data ilun /0/
        character*3 clun
 
        common /image/ n_image

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

        enddo ! nloc

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

        if(l_allsky_clouds .eqv. .true.)then

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

            istat_td = 1

          else
            if(nsmooth .gt. 1)then
                call smooth_box_3d(field_3d,NX_L,NY_L,NZ_L,nsmooth)
            endif

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
!           call make_fnam_lp(i4time_lwc,a9time,istatus)
!           call cv_i4tim_asc_lp(i4time_lwc,a24time,istatus)
            i4time_solar = i4time_lwc
          elseif(c_prodtype .eq. 'F')then 
            var_2d = 'LWC'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,clwc_3d
     1                              ,istat_lwc)
            call cv_i4tim_asc_lp(i4_valid,a24time,istatus)
            i4time_solar = i4_valid
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
            continue
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
            continue
          else
            write(6,*)' Error reading PIC field in plot_allsky'
            return
          endif

500       continue

          if(c_prodtype .eq. 'A')then 
!           Read in swi data
            ext = 'lcv'
            var_2d = 'SWI'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,swi_2d,0,istat_sfc)
            if(istat_sfc .ne. 1)then
              write(6,*)' Error reading LCV/SWI field in plot_allsky'      
              return
            endif

!           Read in snow cover data
            ext = 'lm2'
            var_2d = 'SC'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
     1                      ,ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                      ,snow_cover,0,istat_sfc)
            if(istat_sfc .ne. 1)then
                write(6,*)' Error reading LM2/SC field in plot_allsky'      
                return
            endif
            write(6,*)' range of snow_cover is',
     1                minval(snow_cover),maxval(snow_cover)

!           Read in pw data
            ext = 'lh4'
            var_2d = 'TPW'
            call get_laps_2dgrid(i4time_lwc,0,i4time_nearest
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

          I4_elapsed = ishow_timer()

!         Calculate solar position for 2D array of grid points
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

        else
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
          enddo ! j
          enddo ! i

          swi_2d = 1300. * sind(max(sol_alt_2d,0.))

          read(lun,*) ! advance through input data
          write(6,*)' Running without cloud and other current data'
          snow_cover = r_missing_data

!         Use standard atmosphere for heights (uniform pressure grid)
          call get_pres_1d(i4time_ref,NZ_L,pres_1d,istatus)
          if(istatus .ne. 1)then
            write(6,*)' error getting 1d pressures'
            return
          endif
          do k = 1,NZ_L
            heights_3d(:,:,k) = psatoz(pres_1d(k)/100.)
          enddo ! k

        endif ! l_allsky_clouds is TRUE

        call make_fnam_lp(i4time_solar,a9time,istatus)
        call cv_i4tim_asc_lp(i4time_solar,a24time,istatus)

!       Consider additional albedo info based on land use 
!       This also includes snow cover
        var_2d='USE'
        call read_static_grid(NX_L,NY_L,var_2d,land_use,istatus)
        if(istatus .ne. 1)then
           print*,' Warning: could not read static-landuse'
           return
        endif

        where(topo(:,:) .ge. 3200.); land_use(:,:) = 19.; end where

        call land_albedo_bm(lat,lon,NX_L,NY_L,albedo_bm,istatus)
 
        if(istatus .eq. 1)then
            write(6,*)' Set 3-color albedo based on Blue Marble'

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

            write(6,*)' Row of green albedo'
            jrow = NY_L/2
            do i = 1,NX_L
                if(i .eq. (i/5)*5 .OR. abs(i-NX_L/2) .lt. 20)then
                    write(6,*)i,topo_albedo_2d(2,i,jrow)
                endif
            enddo ! i
            alb_min = minval(topo_albedo_2d(2,:,jrow))
            alb_max = maxval(topo_albedo_2d(2,:,jrow))
            write(6,*)' Min/Max in row is ',alb_min,alb_max

            alb_max = maxval(topo_albedo_2d(2,:,:))
            do i = 1,NX_L
            do j = 1,NY_L
                if(topo_albedo_2d(2,i,j) .eq. alb_max)then
                    write(6,*)' Max albedo at ',i,j
     1                        ,topo_albedo_2d(:,i,j)
                endif
            enddo ! j
            enddo ! i
!           ialbedo_bm(:,:,:) = nint(topo_albedo_2d(:,:,:)*255.)
!           call writeppm3Matrix(
!    1               ialbedo_bm(1,:,:),ialbedo_bm(2,:,:)
!    1              ,ialbedo_bm(3,:,:),'ialbedo_bm')
        else
            write(6,*)' Set 3-color albedo based on land use'
            call land_albedo(land_use,NX_L,NY_L,topo_albedo_2d)
        endif
        
        write(6,*)' Row of albedobm / snow / albedo / topo'
        jrow = NY_L/2
        do i = 1,NX_L
          do j = 1,NY_L
            if(snow_cover(i,j) .ne. r_missing_data)then
              if(topo(i,j) .gt. 2000. .and. topo(i,j) .le. 3500.)then
                snowalb = min(snow_cover(i,j),0.2)
              else
                snowalb = snow_cover(i,j)
              endif
              do ic = 1,3
                topo_albedo_2d(ic,i,j) = 
     1            max(topo_albedo_2d(ic,i,j),snowalb)
              enddo ! ic
            endif
          enddo ! j
          if(i .eq. (i/5)*5 .OR. abs(i-NX_L/2) .lt. 20)then
            write(6,51)i,albedo_bm(2,i,jrow),snow_cover(i,jrow)
     1               ,topo_albedo_2d(2,i,jrow),topo(i,jrow)
 51         format(i5,3f9.3,f9.0)
          endif
        enddo ! i 

        call get_grid_spacing_array(lat,lon,NX_L,NY_L,dx,dy)

        I4_elapsed = ishow_timer()

        call solar_normal(NX_L,NY_L,topo,dx,dy,lat,lon ! I
     1                   ,sol_alt_2d,sol_azi_2d        ! I
     1                   ,alt_norm)                    ! O

        I4_elapsed = ishow_timer()

        write(6,*)' a9time is ',a9time

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

        write(6,91,err=92)pw_ref,aod,aod_ref
91      format(' pw_ref,aod,aod_ref = ',3f10.3)
92      continue

!       Get Atmospheric Optical Depth (3D field)
        call get_aod_3d(pres_3d,heights_3d,topo,NX_L,NY_L,NZ_L,aod_ref
     1                 ,aod_3d)

        I4_elapsed = ishow_timer()
      
        do iloc = 1,nloc
          write(6,*)
          write(6,*)' iloc = ',iloc

          write(6,*)' Enter minalt,maxalt (e.g. 0,90 0,180)'
          read(lun,*)minalt,maxalt           

          write(6,*)' Enter minazi,maxazi'                     
          read(lun,*)minazi,maxazi           
!         minazi = 0
!         maxazi = maxalt * 4

          write(6,*)' Enter alt_scale,azi_scale'               
          read(lun,*)alt_scale,azi_scale     
!         alt_scale = 90. / float(maxalt)
!         azi_scale = alt_scale

          write(6,*)' minalt/maxalt = ',minalt,maxalt
          write(6,*)' minazi/maxazi = ',minazi,maxazi
          write(6,*)' alt_scale/azi_scale = ',alt_scale,azi_scale

          allocate(r_cloud_3d(minalt:maxalt,minazi:maxazi))
          allocate(cloud_od(minalt:maxalt,minazi:maxazi))
          allocate(cloud_od_sp(minalt:maxalt,minazi:maxazi,nsp))
          allocate(blog_v_roll(minalt:maxalt,minazi:maxazi))
          allocate(blog_moon_roll(minalt:maxalt,minazi:maxazi))
          allocate(blog_sun_roll(minalt:maxalt,minazi:maxazi))
          allocate(glow_stars(nc,minalt:maxalt,minazi:maxazi))
          allocate(elong_roll(minalt:maxalt,minazi:maxazi))
          allocate(airmass_2_cloud_3d(minalt:maxalt,minazi:maxazi))
          allocate(airmass_2_topo_3d(minalt:maxalt,minazi:maxazi))
          allocate(topo_swi(minalt:maxalt,minazi:maxazi))
          allocate(topo_albedo(nc,minalt:maxalt,minazi:maxazi))
          allocate(ghic(nc,minalt:maxalt,minazi:maxazi))
          allocate(aod_2_cloud(minalt:maxalt,minazi:maxazi))
          allocate(aod_2_topo(minalt:maxalt,minazi:maxazi))
          allocate(dist_2_topo(minalt:maxalt,minazi:maxazi))
          allocate(aod_ill(minalt:maxalt,minazi:maxazi))
          allocate(aod_ill_dir(minalt:maxalt,minazi:maxazi))
          allocate(aod_tot(minalt:maxalt,minazi:maxazi))
          allocate(r_cloud_trans(minalt:maxalt,minazi:maxazi))
          allocate(cloud_rad_c(nc,minalt:maxalt,minazi:maxazi))
          allocate(cloud_rad_w(minalt:maxalt,minazi:maxazi))
          allocate(clear_rad_c(nc,minalt:maxalt,minazi:maxazi))
          allocate(clear_radf_c(nc,minalt:maxalt,minazi:maxazi))
          allocate(alt_a_roll(minalt:maxalt,minazi:maxazi))
          allocate(azi_a_roll(minalt:maxalt,minazi:maxazi))
          allocate(sky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi))
          allocate(isky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi))

          isound = nint(xsound(iloc))
          jsound = nint(ysound(iloc))

          topo_sfc = topo(isound,jsound)
 
          write(6,*)' i/j/topo_sfc ',isound,jsound,topo_sfc
          write(6,*)' albedo RGB',topo_albedo_2d(:,isound,jsound)

          rlat = lat(isound,jsound)
          rlon = lon(isound,jsound)

          write(6,*)' solar alt/az (2d array)',sol_alt_2d(isound,jsound)       
     1                                        ,sol_azi_2d(isound,jsound) 

!         Calculate solar position for all-sky point
          if(.false.)then
              call solar_position(soundlat(iloc),soundlon(iloc)
     1                           ,i4time_solar,solar_alt     
     1                           ,solar_dec,solar_ha)
              call equ_to_altaz_d(solar_dec,solar_ha,soundlat(iloc)
     1                           ,altdum,solar_az)               
              if(solar_az .lt. 0.)solar_az = solar_az + 360.
          else ! ensure consistency between both solar positions
              solar_alt = sol_alt_2d(isound,jsound)
              solar_az = sol_azi_2d(isound,jsound)              
          endif

          write(6,*)' solar alt/az (all-sky)',solar_alt,solar_az

          write(6,*)' call sun_moon at grid point ',isound,jsound
          idebug = 2
          call sun_moon(i4time_solar,lat,lon,NX_L,NY_L,isound,jsound
     1                 ,alm,azm,idebug,elgms,moon_mag 
     1                 ,solar_eclipse_magnitude,eobsf,eobsl)
          write(6,24)alm,azm,elgms,moon_mag 
24        format('  alt/az/elg/mnmag = ',2f8.2,f9.4,f8.2)

!         Consider passing 'topo_flag' into 'sun_moon' to consider either
!         solar or lunar eclipses
!         http://www.jgisen.de/eclipse
          if(elgms .lt. 0.5)then
            write(6,25)solar_eclipse_magnitude,eobsf,eobsl
25          format(' NOTE: Solar Eclipse Conditions: mag/obsc = '
     1            ,f9.6,2f9.6)
            l_solar_eclipse = .true.
          elseif(elgms .lt. 0.6)then
            write(6,*)' NOTE: Possible Solar Eclipse Conditions'
            l_solar_eclipse = .true.
          else
            l_solar_eclipse = .false.
          endif

          if(maxval(sol_alt_2d) .gt. 0. .and. 
     1       l_solar_eclipse .eqv. .true.)then
            write(6,*)' Calculate gridded eclipse obscuration'
            elgmin = 9999.
            idebug = 0
            do i = 1,NX_L
            do j = 1,NY_L
              call sun_moon(i4time_solar,lat,lon,NX_L,NY_L,i,j
     1                     ,almgrd,azmgrd,idebug,elggrd,grdmoon_mag
     1                     ,solar_eclipse_magnitude,eobsf,eobsc(i,j))

              if(elggrd .lt. elgmin)then
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

          twi_0 = 0.

!         Get line of sight from isound/jsound
          call get_cloud_rays(i4time_solar,clwc_3d,cice_3d
     1                     ,heights_3d                           ! I
     1                     ,rain_3d,snow_3d                      ! I
     1                     ,pres_3d,aod_3d,topo_sfc,topo,swi_2d  ! I
     1                     ,topo_albedo_2d                       ! I
     1                     ,topo_swi,topo_albedo,ghic            ! O
!    1                     ,ghi_2d,dhi_2d                        ! O
     1                     ,aod_vrt,aod_2_cloud,aod_2_topo       ! O
     1                     ,dist_2_topo                          ! O
     1                     ,aod_ill,aod_ill_dir                  ! O
     1                     ,aod_tot,transm_obs                   ! O
     1                     ,transm_3d,transm_4d                  ! O
     1                     ,r_cloud_3d,cloud_od,cloud_od_sp      ! O
     1                     ,r_cloud_trans,cloud_rad_c,cloud_rad_w! O
     1                     ,clear_rad_c,clear_radf_c,patm        ! O
     1                     ,airmass_2_cloud_3d,airmass_2_topo_3d ! O
     1                     ,htmsl                                ! O
     1                     ,htagl(iloc)                          ! I
     1                     ,aod_ref                              ! I
     1                     ,NX_L,NY_L,NZ_L,isound,jsound         ! I
     1                     ,alt_a_roll,azi_a_roll                ! I
     1                     ,sol_alt_2d,sol_azi_2d                ! I
     1                     ,alt_norm                             ! I
     1                     ,moon_alt_2d,moon_azi_2d              ! I
     1                     ,moon_mag,moon_mag_thr,twi_0          ! I
     1                     ,l_solar_eclipse,rlat,rlon,lat,lon    ! I
     1                     ,minalt,maxalt,minazi,maxazi          ! I
     1                     ,alt_scale,azi_scale                  ! I
     1                     ,grid_spacing_m,r_missing_data)       ! I

          write(6,*)' Return from get_cloud_rays: ',a9time
     1             ,' aod_vrt is ',aod_vrt

 900      continue

1000      continue

          ilun = ilun + 1
          write(clun,14)ilun
14        format(i3.3)

!         Moon glow in cylindrical coordinates (add color info)?                   
          blog_moon_roll = 0.
!         if(moon_mag .lt. moon_mag_thr .AND.
          if(.true.                     .AND.
     1       alm      .gt. 0.           .AND.
     1       l_solar_eclipse .eqv. .false.    )then
            write(6,*)' Moon glow being calculated: ',alm,azm
            diam_deg = 0.5
            call get_glow_obj(i4time,alt_a_roll,azi_a_roll
     1                       ,minalt,maxalt,minazi,maxazi 
     1                       ,alt_scale,azi_scale
     1                       ,alm,azm,moon_mag,diam_deg
     1                       ,blog_moon_roll)

            write(6,*)' range of blog_moon_roll is',
     1          minval(blog_moon_roll),maxval(blog_moon_roll)
          endif

!         Sun glow in cylindrical coordinates, treated as round?
          blog_sun_roll = 0
          if(l_solar_eclipse .eqv. .true.)then
              if(eobsl .ge. 1.0)then
                  s_mag = -12.5
                  diam_deg = 8.0    
              else
                  s_mag = -26.74 - (log10(1.0-eobsl))*2.5
                  if(s_mag .lt. -12.5)then
                      diam_deg = 0.5    
                  else ! show corona even outside totality
                      diam_deg = 8.0    
                  endif
              endif
          else
              s_mag = -26.74
              diam_deg = 0.5
          endif
          write(6,*)' Sun glow being calculated: '
     1                 ,solar_alt,solar_az,s_mag
          call get_glow_obj(i4time,alt_a_roll,azi_a_roll
     1                     ,minalt,maxalt,minazi,maxazi 
     1                     ,alt_scale,azi_scale
     1                     ,solar_alt,solar_az,s_mag,diam_deg
     1                     ,blog_sun_roll)
          write(6,*)' range of blog_sun_roll is',
     1        minval(blog_sun_roll),maxval(blog_sun_roll),diam_deg

!         if(solar_alt .ge. 0.)then
          if(.true.)then
            I4_elapsed = ishow_timer()
!           write(6,*)' call get_skyglow_cyl for sun or moon...'

!           if(solar_alt .ge. -16.)then
            if(.false.)then
                write(6,*)' Sun is significant, alt is:',solar_alt
                call skyglow_cyl(solar_alt,solar_az,blog_v_roll  ! IO
     1                          ,elong_roll,aod_vrt              ! OI
     1                          ,minalt,maxalt,minazi,maxazi     ! I
     1                          ,alt_scale,azi_scale)            ! I
                I4_elapsed = ishow_timer()

                if(moon_mag .lt. moon_mag_thr .AND.
     1             alm      .gt. 0.                 )then
                    if(solar_alt .ge. 0.)then ! valid only in daytime
                        write(6,*)' Moonglow is being added to daylight'       
                        do j = minazi,maxazi
                        do i = minalt,maxalt
                           blog_v_roll(i,j) = 
     1                     addlogs(blog_v_roll(i,j),blog_moon_roll(i,j))
                        enddo ! i
                        enddo ! j
                    endif
                endif
            endif ! solar_alt .ge. -16.

!           if(solar_alt .lt. -16. .AND. moon_mag .lt. moon_mag_thr
!    1                             .AND. alm .gt. 0.         )then
            if(.false.)then
                write(6,*)' Moon skyglow significant: mag ',moon_mag
                call skyglow_cyl(alm,azm,blog_v_roll,elong_roll,aod_vrt 
     1                          ,minalt,maxalt,minazi,maxazi
     1                          ,alt_scale,azi_scale)
                blog_v_roll = blog_v_roll + (-26.74 - moon_mag) * 0.4
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
!         call cyl_to_polar(r_cloud_3d,r_cloud_3d_polar,minalt,maxalt
!    1                   ,maxazi,alt_scale,azi_scale
!    1                   ,alt_a_polar,azi_a_polar
!    1                   ,ni_polar,nj_polar)
!         write(6,*)' cyl slice at 40alt ',r_cloud_3d(40,:)
!         write(6,*)' polar slice at 256 ',r_cloud_3d_polar(256,:)

!         Reproject Skyglow Field
!         call cyl_to_polar(blog_v_roll,blog_v_roll_polar
!    1                               ,minalt,maxalt,maxazi
!    1                               ,alt_scale,azi_scale
!    1                               ,alt_a_polar,azi_a_polar
!    1                               ,ni_polar,nj_polar)

!         Write time label
          open(53,file='label.'//clun,status='unknown')
          write(53,*)a9time
          write(53,*)a24time(1:17)
          close(53)

          if(.true.)then

!           Get all sky for cyl   
            ni_cyl = maxalt - minalt + 1
            nj_cyl = maxazi - minazi + 1

!           Write lat/lon and other info for label
            open(54,file='label2.'//clun,status='unknown')
            write(54,54)soundlat(iloc),soundlon(iloc),
     1                  minalt,maxalt,minazi,maxazi,ni_cyl,nj_cyl,
     1                  solar_alt,solar_az,alt_scale,azi_scale
            close(54)
 54         format(2f8.2/6i8/2f8.2,2f7.2)

            I4_elapsed = ishow_timer()

            if(solar_alt .lt. 0. .OR. l_solar_eclipse .eqv. .true.)then
              write(6,*)' call get_starglow with cyl data'
              call get_starglow(i4time_solar,alt_a_roll,azi_a_roll       ! I
     1                     ,minalt,maxalt,minazi,maxazi                  ! I
     1                     ,rlat,rlon,alt_scale,azi_scale                ! I
     1                     ,glow_stars)                                  ! O

              write(6,*)' range of glow_stars (before) is',
     1             minval(glow_stars(2,:,:)),maxval(glow_stars(2,:,:))

              write(6,*)' range of moonglow is',
     1             minval(blog_moon_roll),maxval(blog_moon_roll)

              write(6,*)' Moonglow is being added to starlight'       
              do j = minazi,maxazi
              do i = minalt,maxalt
                do ic = 1,nc
                  glow_stars(ic,i,j) = 
     1              addlogs(glow_stars(ic,i,j),blog_moon_roll(i,j))
                enddo ! ic
              enddo ! i 
              enddo ! j 

              write(6,*)' range of glow_stars (after) is',
     1             minval(glow_stars(2,:,:)),maxval(glow_stars(2,:,:))
            endif

            I4_elapsed = ishow_timer()

            if(solar_alt .gt. -3.)then
                call get_idx(solar_alt,minalt,alt_scale,ialt_sun)
                call get_idx(solar_az ,minazi,azi_scale,jazi_sun)
                ialt_sun = ialt_sun - minalt + 1
                jazi_sun = jazi_sun - minazi + 1
                write(6,*)' solar_alt,minalt,alt_scale = '
     1                     ,solar_alt,minalt,alt_scale
                write(6,*)' ialt_sun,jazi_sun = ',ialt_sun,jazi_sun
            endif

            write(6,*)' call get_sky_rgb with cyl data'
            call get_sky_rgb(r_cloud_3d      ! cloud opacity
     1                    ,cloud_od          ! cloud optical depth
     1                    ,cloud_od_sp,nsp   ! cloud species optical depth
     1                    ,r_cloud_trans     ! cloud solar transmittance
     1                    ,cloud_rad_c       ! cloud solar transmittance / color
     1                    ,cloud_rad_w       ! cloud solar transmittance * rad
     1                    ,clear_rad_c       ! clear sky illumination by sun     
     1                    ,l_solar_eclipse,i4time_solar,rlat,rlon,eobsl
     1                    ,clear_radf_c      ! clear sky frac illumination by sun     
     1                    ,patm,htmsl
     1                    ,blog_v_roll       ! skyglow
     1                    ,blog_sun_roll     ! sunglow
     1                    ,blog_moon_roll    ! moonglow
     1                    ,glow_stars        ! starglow
     1                    ,aod_vrt 
     1                    ,transm_obs        ! observer illumination
     1                    ,ialt_sun,jazi_sun ! sun location
     1                    ,airmass_2_cloud_3d      
     1                    ,airmass_2_topo_3d      
     1                    ,topo_swi,topo_albedo
     1                    ,topo_albedo_2d(2,isound,jsound)
     1                    ,aod_2_cloud,aod_2_topo,aod_ill,aod_ill_dir
     1                    ,dist_2_topo
     1                    ,alt_a_roll,azi_a_roll ! I   
     1                    ,elong_roll    
     1                    ,ni_cyl,nj_cyl  
     1                    ,solar_alt,solar_az,twi_0 
     1                    ,alm,azm,moon_mag  ! moon alt/az/mag
     1                    ,sky_rgb_cyl)   

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
               write(6,*)'iaz,cyl(maxalt/2,iaz)',iaz
     1                       ,isky_rgb_cyl(:,maxalt/2,iaz)
              enddo ! iaz
            endif

            if(l_polar .eqv. .true.)then
!             Reproject sky_rgb array from cyl to polar    
              do iaz = minazi,maxazi,20
               write(6,*)'iaz,cyl(maxalt/2,iaz)',iaz
     1                       ,sky_rgb_cyl(1,maxalt/2,iaz)
              enddo ! iaz

              write(6,*)' Call cyl_to_polar with sky rgb data'

              do ic = 0,nc-1
                call cyl_to_polar(sky_rgb_cyl(ic,:,:)
     1                           ,sky_rgb_polar(ic,:,:)
     1                           ,minalt,maxalt,minazi,maxazi
     1                           ,alt_scale,azi_scale
     1                           ,alt_a_polar,azi_a_polar
     1                           ,ni_polar,nj_polar)
              enddo ! ic

!             Write all sky for polar
              isky_rgb_polar = sky_rgb_polar
              where(sky_rgb_polar .eq. r_missing_data)
     1              isky_rgb_polar = 0
!             write(6,*)' ipolar array',
!    1            isky_rgb_polar(2,ni_polar/2,1:nj_polar)
              npts = 3*ni_polar*nj_polar
!             write(6,*)' Write all sky polar text file'
!    1                  ,isky_rgb_polar(:,255,255),npts
!             open(54,file='allsky_rgb_polar.'//clun,status='unknown')
!             write(54,*)isky_rgb_polar
!             close(54)
              write(6,*)' Write all sky polar ppm file '
              call writeppm3Matrix(
     1                  isky_rgb_polar(0,:,:),isky_rgb_polar(1,:,:)
     1                 ,isky_rgb_polar(2,:,:),'allsky_rgb_polar_'//clun)
            endif

          endif ! mode_polar = 0 or 2

          write(6,*)' End of plot_allsky for iloc...',iloc
          write(6,*)

          I4_elapsed = ishow_timer()

          deallocate(r_cloud_3d)
          deallocate(cloud_od)
          deallocate(cloud_od_sp)
          deallocate(blog_v_roll)
          deallocate(blog_moon_roll)
          deallocate(blog_sun_roll)
          deallocate(glow_stars)
          deallocate(elong_roll)
          deallocate(airmass_2_cloud_3d)
          deallocate(airmass_2_topo_3d)
          deallocate(topo_swi)
          deallocate(topo_albedo)
          deallocate(ghic)
          deallocate(aod_2_cloud)
          deallocate(aod_2_topo)
          deallocate(dist_2_topo)
          deallocate(aod_ill)
          deallocate(aod_ill_dir)
          deallocate(r_cloud_trans)
          deallocate(cloud_rad_c)
          deallocate(cloud_rad_w)
          deallocate(clear_rad_c)
          deallocate(clear_radf_c)
          deallocate(alt_a_roll)
          deallocate(azi_a_roll)
          deallocate(sky_rgb_cyl)
          deallocate(isky_rgb_cyl)

        enddo ! iloc

        write(6,*)
        write(6,*)' End of plot_allsky...'
        write(6,*)

        return
        end
