
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
        use wrf_lga

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
        real snow_albedo_max(NX_L,NY_L)
        real snow_depth(NX_L,NY_L)
        real seaice(NX_L,NY_L)
        real topo_albedo_2d(nc,NX_L,NY_L)
        real albedo_bm(nc,NX_L,NY_L)
        real bm_counts(nc,NX_L,NY_L)
        integer ialbedo_bm(nc,NX_L,NY_L)
        real albedo_usgs(nc,NX_L,NY_L)
        real lat(NX_L,NY_L)
        real lon(NX_L,NY_L)
        real*8 dlat(NX_L,NY_L)
        real*8 dlon(NX_L,NY_L)
        real topo(NX_L,NY_L)
        real dum_2d(NX_L,NY_L)
        real dx(NX_L,NY_L)
        real dy(NX_L,NY_L)
        real alt_norm(NX_L,NY_L)     ! Solar Alt w.r.t. terrain normal
        real sol_alt_2d(NX_L,NY_L)
        real sol_azi_2d(NX_L,NY_L)
        real eobsc(NX_L,NY_L)        ! array of 'eobsl' values
        real moon_alt_2d(NX_L,NY_L)
        real moon_azi_2d(NX_L,NY_L)
        real moon_mag,moon_mag_thr

        parameter (mxopt = 10)
        real*8 a_vec(mxopt),a_last(mxopt),depth_this_run,f_merit
        real*8 dstep(mxopt),dstep_gran(mxopt)

        real pres_1d(NZ_L)

        real lil_sfc, lic_sfc, lil_cpt, lic_cpt
  
        real k_to_c, make_td, make_ssh, mfpath

        character*1 c_prodtype, c_plotobs, c1_lat, c1_lon
        character*3 var_2d
        character*150  directory, filename, wrfout_full, ramsout_full
        character*31  ext
        character*10  units_2d
        character*125 comment_2d
        character*40 c_model /''/
        character*9 a9time
        character*24 a24time
        character*5 fcst_hhmm
        character*3 c3_string
        character*4 c4_string
        character*40 c_label
        character*11 c_pw
        character*11 c_cape
        character*20 c20_x, c20_y
        character*255 new_dataroot, filename_ppm
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
        real, allocatable, dimension(:,:,:) :: camera_rgbf
        real*8, allocatable, dimension(:,:) :: dist_2_topo

        real alt_a_polar(iplo:iphi,jplo:jphi)
        real azi_a_polar(iplo:iphi,jplo:jphi)
        real elong_a_polar(iplo:iphi,jplo:jphi)

        integer maxloc
        parameter (maxloc = 1000)
        integer minalt_a(maxloc),maxalt_a(maxloc)
        integer minazi_a(maxloc),maxazi_a(maxloc)
        real alt_scale_a(maxloc),azi_scale_a(maxloc)

        real, allocatable, dimension(:,:,:) :: sky_rgb_polar
        integer, allocatable, dimension(:,:,:) :: isky_rgb_polar
!       real sky_rgb_polar(0:2,ni_polar,nj_polar)
!       integer isky_rgb_polar(0:2,ni_polar,nj_polar)

        real, allocatable, dimension(:,:,:) :: sky_rgb_cyl
        integer, allocatable, dimension(:,:,:) :: isky_rgb_cyl
        real correlation(nc,maxloc)
        integer mode_cloud_mask /4/ ! ignore the mask

        real*8 xsound(maxloc),ysound(maxloc)
        real*8 soundlat(maxloc),soundlon(maxloc)
        real fsoundlat(maxloc),fsoundlon(maxloc)
        real htagl(maxloc),exposure_a(maxloc)

        data ilun /0/
        character*3 clun
        character*10 clun_loop
 
        data i_aero_synplume /0/
        data i_aero_1d /1/

        common /image/ n_image

        rpd = 3.141592653589/180.

!       Initialize
        snow_cover = r_missing_data        
        snow_albedo_max = r_missing_data        
        snow_depth = r_missing_data        
        seaice = r_missing_data

        mil = 1
        mih = NX_L
        mjl = 1
        mjh = NY_L

        call alloc_allsky(NX_L,NY_L,NZ_L,nc,istatus)

        if(l_polar .eqv. .true.)then
            mode_polar = 2
        else
            mode_polar = 0
        endif

        nsmooth = 1

        I4_elapsed = ishow_timer()

        call GETENV('C_MODEL',c_model)

        aod_tmp  = getenv_real('AOD',r_missing_data)
        if(aod_tmp .ne. r_missing_data)aod = aod_tmp

        write(6,*)
        write(6,*)' subroutine plot_allsky: nsmooth/aod is ',nsmooth,aod
        write(6,*)' l_cyl/l_polar = ',l_cyl,l_polar
        write(6,*)' ipolar_sizeparm = ',ipolar_sizeparm
        write(6,*)' density = ',density
        write(6,*)' c_model = ',c_model

        if(trim(c_model) .eq. 'optimize')then
          mode_cloud_mask = 5         
        endif

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
        dlat = dble(lat)

        var_2d='LON'
        call read_static_grid(NX_L,NY_L,var_2d,lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-lon'
            return
        endif
        dlon = dble(lon)

        var_2d='AVG'
        call read_static_grid(NX_L,NY_L,var_2d,topo,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-topo'
            return
        endif

!       Topo offset (based on topo_1s data over Colorado)
        if(grid_spacing_m .le. 30.)then
          ioff = -1  ! move the terrain eastward this amount
          joff = +13 ! move the terrain northward this amount

          ioff = +0  ! move the terrain eastward this amount
          joff = +12 ! move the terrain northward this amount

          do idum = 1,NX_L
          do jdum = 1,NY_L
            i = min(max(idum-ioff,1),NX_L)
            j = min(max(jdum-joff,1),NY_L)
            dum_2d(idum,jdum) = topo(i,j)
          enddo ! j
          enddo ! i
          topo(:,:) = dum_2d(:,:)

!         Find Summits diagnostic near domain center
          ibox = 13
          ibox2 = ibox/2
          imid = NX_L/2
          jmid = NY_L/2
          do is = imid-300,imid+300
          do js = jmid-300,jmid+300
            topomin = 99999.                
            topomax = 0.
            do iis = is-ibox2,is+ibox2
            do jjs = js-ibox2,js+ibox2
              if(is .ne. iis .or. js .ne. jjs)then
                topomax = max(topomax,topo(iis,jjs))
                topomin = min(topomin,topo(iis,jjs))
              endif
            enddo ! jjs
            enddo ! iis
            if(topo(is,js) .gt. topomax .and.
     1         topo(is,js) .gt. topomin+150. .and.
     1         topo(is,js) .gt. 2430.              )then
              write(6,14)is,js,lat(is,js),lon(is,js),topo(is,js),topomax
     1                                                          ,topomin
14            format(' Peak detected at',2i5,2f11.5,3f9.2)           
            endif
          enddo ! js
          enddo ! is

        endif

!       Generate default 'snow_albedo_max' field
        do i = 1,NX_L
          do j = 1,NY_L
            if(topo(i,j) .gt. 1900. .and. topo(i,j) .le. 3500.)then
              snow_albedo_max(i,j) = 0.4
            else
              snow_albedo_max(i,j) = 0.7
            endif
          enddo ! j
        enddo ! i

        write(6,*)' line 290 aod is ',aod

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
            fsoundlat(iloc) = sngl(soundlat(iloc))

            write(6,*)' Input longitude for allsky plot...'       
            read(5,*)c20_y
            call s_len(c20_y,leny)
            if(l_parse(c20_y(1:leny),'l'))then
                read(c20_y(1:leny-1),*)soundlon(iloc)
            else
                read(c20_y(1:leny),*)soundlon(iloc)
            endif
            fsoundlon(iloc) = sngl(soundlon(iloc))

            write(6,*)' observer lat/lon read in =',soundlat(iloc)
     1                                             ,soundlon(iloc)


            write(6,*)' line 330 aod is ',aod

            call latlon_db_rlapsgrid(soundlat(iloc),soundlon(iloc)
     1                              ,dlat,dlon,NX_L,NY_L
     1                              ,xsound(iloc),ysound(iloc),istatus)
            if(istatus .ne. 1)then
!             if(trim(c_model) .ne. 'hrrr_smoke')then
              if(len(c_model) .eq. 0)then
                write(6,*)' Station is outside domain - try again...'
                return
              else
                write(6,*)' WARNING: Station is outside domain...'
              endif
            endif

            if(xsound(iloc) .lt. 1d0 .or. xsound(iloc) .gt. dble(NX_L)
     1                              .or.
     1         ysound(iloc) .lt. 1d0 .or. ysound(iloc) .gt. dble(NY_L)    
     1                                                             )then
!             if(trim(c_model) .ne. 'hrrr_smoke')then
              if(len(c_model) .eq. 0)then
                write(6,*)' Station is outside domain - try again...'
                return
              else
                write(6,*)' WARNING: Station is outside domain...'
              endif
            endif

          else
            read(c20_x(1:lenx),*)xsound(iloc)
            write(6,*)' Input y grid point for allsky plot...'
            read(5,*)ysound
            if(xsound(iloc) .lt. 1d0 .OR. ysound(iloc) .lt. 1d0)then ! scale domain
                write(6,*)' Values less than 1.0, scale to domain'
                xsound(iloc) = nint(xsound(iloc) * (NX_L-1)) + 1
                ysound(iloc) = nint(ysound(iloc) * (NY_L-1)) + 1
                soundlat(iloc) = 
     1              lat(nint(xsound(iloc)),nint(ysound(iloc)))
                soundlon(iloc) = 
     1              lon(nint(xsound(iloc)),nint(ysound(iloc)))
            endif
          endif

          read(5,*)htagl(iloc),exposure_a(iloc)

          write(6,*)' soundlat/soundlon ',soundlat(iloc),soundlon(iloc)      
          write(6,*)' xsound/ysound ',xsound(iloc),ysound(iloc)
          write(6,*)' htagl/exposure_a ',htagl(iloc),exposure_a(iloc)

        enddo ! iloc

        write(6,*)' line 370 aod is ',aod

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

            aod_tmp  = getenv_real('AOD',r_missing_data)
            if(aod_tmp .ne. r_missing_data)aod = aod_tmp
        endif

        tlow_c = -30.
        thigh_c = +50.

        write(6,*)' line 410 aod is ',aod

!       Get 3-D pressure field
        call get_pres_3d(i4_valid,NX_L,NY_L,NZ_L,pres_3d,istatus)
        if(istatus .ne. 1)go to 900

        write(6,*)' time diff is ',i4time_ref - i4time_now_gg()

        l_water_world = .false.
        if(l_parse(directory,'fim'))then                   ! special routines
          l_require_all_fields = .false.
        elseif(l_parse(directory,'rams'))then
          l_require_all_fields = .false.
        elseif(l_parse(directory,'navgem'))then
          l_require_all_fields = .false.
        elseif(c_model(1:5) .eq. 'cloud')then  ! e.g. cloud_13_13_00001_00002
          write(6,*)' parse c_model for idealized cloud tests'
          l_require_all_fields = .false.                   ! 
          l_test_cloud = .true.
          read(c_model(7:8),*)lvl1 
          read(c_model(10:11),*)lvl2
          read(c_model(13:17),*)clwc_ideal
          read(c_model(19:23),*)cice_ideal
          clwc_ideal = clwc_ideal * .00001
          cice_ideal = cice_ideal * .00001
        elseif(c_model(1:9) .eq. 'aerocloud')then  
          l_require_all_fields = .false.                   ! 
          mode_aero_cld = 3
          i_aero_synplume = 2
!         read(c_model(11:13),*)aero_synfactor
        elseif(c_model(1:4) .eq. 'aero')then               ! e.g. aeroloop
          l_require_all_fields = .true.                    ! 
!         mode_aero_cld = 3
!         i_aero_synplume = 2
!         read(c_model(11:13),*)aero_synfactor
        elseif(trim(c_model) .eq. 'optimize')then          
          l_require_all_fields = .true.
        elseif(i4time_ref - i4time_now_gg() .gt. 20e6)then ! >0.6y future
          l_require_all_fields = .false.                   ! test clouds
          l_test_cloud = .true.
          lvl1 = 13
          lvl2 = 13
          clwc_ideal = .000
          cice_ideal = .000
        elseif(trim(c_model) .ne. '')then                  ! e.g. HRRR-AK
          l_require_all_fields = .false.
        elseif(grid_spacing_m .le. 30.)then
          l_require_all_fields = .false.                   ! water only
        elseif(i4time_ref - i4time_now_gg() .lt. 1e6)then  ! present/past
          l_require_all_fields = .true.
        elseif(i4time_ref - i4time_now_gg() .gt. 75e6)then ! >2.5y future 
          l_require_all_fields = .false.                   ! water only
          l_test_cloud = .false.
          l_water_world = .true.
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
           write(6,*)' generate idealized cloud fields - 1'
           clwc_3d(:,:,14) = .001 ! 14/775mb .001
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
          snow_cover = 0. ! initialize

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
            call GETENV('RAMSOUT_FULL',ramsout_full)
            filename = trim(ramsout_full)
!           filename = '/home/fab/albers/muri/rams_micro_v3.nc'
            call get_rams_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L,NZ_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lat,lon
     +                   ,pres_3d
     +                   ,heights_3d
     +                   ,clwc_3d
     +                   ,cice_3d
     +                   ,rain_3d
     +                   ,snow_3d
     +                   ,aod_3d
     +                   ,lun_out
     +                   ,istatus)
            write(6,*)' returned from get_rams_data ',istatus
            if(istatus .ne. 1)then
                stop
            endif

            istatus_ht = 1

            mode_aero_cld = 3
            aod = 0.

            aod_3d = max(aod_3d,1e-9)
            if(.false.)then ! zero out hydrometeors
              clwc_3d = max(clwc_3d * 1e-3,0.)
              cice_3d = max(cice_3d * 1e-3,0.)
              rain_3d = max(rain_3d * 1e-3,0.)
              snow_3d = max(snow_3d * 1e-3,0.)
            endif

            i_aero_1d = 0 ! retain the 3D aerosols read in
            write(6,*)' RAMS run: set mode_aero_cld = ',mode_aero_cld
            write(6,*)' aod_3d range is ',minval(aod_3d),maxval(aod_3d)
            write(6,*)' clwc_3d range is ',minval(clwc_3d)
     1                                    ,maxval(clwc_3d)
            write(6,*)' cice_3d range is ',minval(cice_3d)
     1                                    ,maxval(cice_3d)
            write(6,*)' rain_3d range is ',minval(rain_3d)
     1                                    ,maxval(rain_3d)
            write(6,*)' snow_3d range is ',minval(snow_3d)
     1                                    ,maxval(snow_3d)

!           Zero out hydrometeors
            if(.false.)then
              write(6,*)' RAMS run: zero out hydrometeors'
              clwc_3d = 0.
              cice_3d = 0.
              rain_3d = 0.
              snow_3d = 0.
            endif

          elseif(l_parse(directory,'navgem'))then
            write(6,*)' Looking for NAVGEM data in ',trim(directory)
            filename = '/home/fab/albers/muri/navgem/file.nc'
            call get_navgem_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L,NZ_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,pres_3d
     +                   ,clwc_3d
     +                   ,cice_3d
     +                   ,rain_3d
     +                   ,snow_3d
     +                   ,seaice,snow_depth
     +                   ,lun_out
     +                   ,istatus)
            istatus_ht = 0 ! we can perhaps try the height_AGL 
            write(6,*)' returned from get_navgem_data'

            mode_aero_cld = 3
            write(6,*)' NAVGEM run: mode_aero_cld = ',mode_aero_cld

            i_aero_1d = 0 ! retain the 3D aerosols read in from NAVGEM
            i_aero_synplume = 0

!           Zero out hydrometeors
            if(.false.)then
              write(6,*)' NAVGEM run: zero out hydrometeors'
              clwc_3d = 0.
              cice_3d = 0.
              rain_3d = 0.
              snow_3d = 0.
            endif

          elseif(trim(c_model) .eq. 'hrrr_smoke' .or.
     +           trim(c_model) .eq. 'wrf_chem' .or. 
     +           trim(c_model) .eq. 'hrrr_ak'         )then

            write(6,*)' Getting 1D pressure levels'
            call get_pres_1d(i4time_sys,NZ_L,pres_1d,istatus)

            if(trim(c_model) .eq. 'hrrr_smoke')then
              mode_aero_cld = 3
              aod = 0.
              i_aero_1d = 0 ! retain the 3D aerosols read in

!             Speed up in 'get_cloud_rad'
              mil = 1
              mih = NX_L/2
              mjl = NY_L/2
              mjh = NY_L
              itype_aod = 2 ! 1 is TAOD5503D, 2 is tracer_1a              
            elseif(trim(c_model) .eq. 'hrrr_ak')then
              mode_aero_cld = 3
              aod = 0.
              i_aero_1d = 0 ! retain the 3D aerosols read in
              itype_aod = 2 ! 1 is TAOD5503D, 2 is tracer_1a              
            elseif(trim(c_model) .eq. 'wrf_chem')then ! MURI
              mode_aero_cld = 3
              aod = 0.
              i_aero_1d = 0 ! retain the 3D aerosols read in
              itype_aod = 3 ! 1 is TAOD5503D, 2 is tracer_1a, 3 is MURI
            endif

!           if(trim(c_model) .eq. 'hrrr_ak')then
              call GETENV('WRFOUT_FULL',wrfout_full)
              filename = trim(wrfout_full)
!           else
!             filename = 'file.nc'
!           endif


            write(6,*)' Looking for WRF SWIM data in ',trim(filename) ! trim(directory)
            call wrf2swim(filename,i4time_sys,itype_aod,NX_L,NY_L,NZ_L
     +           ,lat,lon    
     +           ,pres_1d,land_frac,snow_cover,snow_albedo_max
     +           ,istatus)
            if(istatus .ne. 1)then
               write(6,*)' returning: wrf2swim istatus is ',istatus
               return
            endif
            istatus_ht = 0 ! unless we are reading height
            write(6,*)' returned from wrf2swim for ',trim(c_model)

            if(trim(c_model) .eq. 'hrrr_smoke')then
              write(6,*)'   aod max, mean free path, vsby at each level'
              write(6,*)'     m^-1          m                m'
              do ka = 1,NZ_L
                aodmax = maxval(aod_3d(:,:,ka))
                do ia = 1,NX_L
                do ja = 1,NY_L
                  if(aod_3d(ia,ja,ka) .eq. aodmax)then
                    imaxa = ia
                    jmaxa = ja 
                  endif
                enddo ! ja
                enddo ! ia
                if(aodmax .gt. 0.)then
                  mfpath = 1. / aodmax
                  vismax =  mfpath * (-log(.02))
                else
                  mfpath = 0.
                  vismax = 0.
                endif
                write(6,151)ka,pres_1d(ka),aodmax,mfpath,vismax
     1                        ,lat(imaxa,jmaxa),lon(imaxa,jmaxa)
151             format(i3,f9.0,e13.5,2x,2f12.1,2f9.3)            
              enddo
            endif

            if(trim(c_model) .eq. 'hrrr_smoke' .and. .false.)then
              write(6,*)' Zero out hydrometeors'
              clwc_3d = 0.
              cice_3d = 0.
              rain_3d = 0.
              snow_3d = 0.
            endif

          else
            istatus_ht = 0

          endif

          i4time_solar = i4time_ref

          I4_elapsed = ishow_timer()

!         Calculate solar position for 2D array of grid points
          write(6,*)' Getting 2D solar position'
          if(.true.)then
           call get_solaltaz_2d(lat,lon,i4time_solar,NX_L,NY_L
     1                         ,sol_alt_2d,sol_azi_2d)
           icen = NX_L/2
           jcen = NY_L/2
           write(6,*)' solar alt/az (center)',sol_alt_2d(icen,jcen)
     1                                       ,sol_azi_2d(icen,jcen)

           swi_2d(:,:) = 1300. * sin(max(sol_alt_2d(:,:)*rpd,0.))
           write(6,*)' solar swi (center)',swi_2d(icen,jcen)

          else
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
          endif

          I4_elapsed = ishow_timer()

          read(lun,*) ! advance through input data
          write(6,*)' Running without LAPS cloud and other current data'

          if(l_test_cloud)then ! recently active section
            write(6,*)' generate idealized cloud fields - 2',clwc_ideal
     1                                                      ,cice_ideal
            clwc_3d(:,:,:) = 0.
            cice_3d(:,:,:) = 0.
            clwc_3d(:,:,lvl1:lvl2) = clwc_ideal ! 13/800mb .001
            cice_3d(:,:,lvl1:lvl2) = cice_ideal ! 13/800mb .001
          endif

!         Use standard atmosphere for heights (uniform pressure grid)
          if(istatus_ht .eq. 0)then
            write(6,*)' Getting 1D pressure levels'
            call get_pres_1d(i4time_ref,NZ_L,pres_1d,istatus)
            if(istatus .ne. 1)then
              write(6,*)' error getting 1d pressures'
              return
            endif
            I4_elapsed = ishow_timer()
            write(6,*)' Getting 3D heights from 1D pressure levels'
            do k = 1,NZ_L
              heights_3d(:,:,k) = psatoz(pres_1d(k)/100.)
            enddo ! k
            I4_elapsed = ishow_timer()
          endif

        endif ! l_require_all_fields is TRUE

        call make_fnam_lp(i4time_solar,a9time,istatus)
        call cv_i4tim_asc_lp(i4time_solar,a24time,istatus)

        read(a9time(3:5),*)idoy
        r_au = radnorm(idoy)

        write(6,*)' Solar dist (AU) ',r_au

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
     1                     ,albedo_bm,bm_counts,istat_bm)

        I4_elapsed = ishow_timer()
 
        if(istat_bm .eq. 1)then

            if(maxval(albedo_bm) .gt. 1.)then
               write(6,*)' ERROR: land_albedo_bm has invalid values'
               stop
            endif

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
     1        ' Row of spectral albedo / counts (through domain center)'
            jrow = NY_L/2
            do i = 1,NX_L
                if(i .eq. (i/5)*5 .OR. abs(i-NX_L/2) .lt. 20)then
                    write(6,16)i,lat(i,jrow),lon(i,jrow)
     1                    ,land_frac(i,jrow),topo_albedo_2d(:,i,jrow)
     1                    ,bm_counts(:,i,jrow)
16                  format(i5,3f9.3,2x,3f9.3,2x,3f9.3)                 
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
            if(seaice(i,j) .ne. r_missing_data)then ! NAVGEM
              snowalb = 0.8
              do ic = 1,3
                topo_albedo_2d(ic,i,j) = seaice(i,j) * snowalb + 
     1                (1.0-seaice(i,j)) * topo_albedo_2d(ic,i,j)
              enddo ! ic
            elseif(snow_depth(i,j) .ne. r_missing_data)then ! NAVGEM
              snowalb = 0.8
              snowcvr = max(snow_depth(i,j) / 0.03,1.0)
              do ic = 1,3
                topo_albedo_2d(ic,i,j) = snowcvr * snowalb + 
     1                (1.0-snowcvr) * topo_albedo_2d(ic,i,j)
              enddo ! ic
            elseif(snow_cover(i,j) .ne. r_missing_data)then ! Other models
!             snowalb = snow_cover(i,j) * snow_albedo_max(i,j)
              snowalb = snow_albedo_max(i,j)
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
           write(6,*)' aod is >= 0, use directly for aod_ref',aod
           aod_ref = aod
        endif

!       Get Aersol Extinction Coefficient (3D field)
        call get_aod_3d(pres_3d,heights_3d,topo,NX_L,NY_L,NZ_L
     1                 ,aod,aod_ref,i_aero_synplume,i_aero_1d,aod_3d)

        write(6,19,err=20)pw_ref,aod,aod_ref
19      format(' pw_ref,aod,aod_ref = ',3f10.3)
20      continue

        I4_elapsed = ishow_timer()

        rlat_last = -999.
        rlon_last = -999.
        i4time_last = 0.

        if(mode_cloud_mask .eq. 5)then
            nloops = 75
            write(6,*)' Optimize mode: initialize',nloops
            filename_ppm = ''

!           Number and permutation of variables
            nv = 3
            n_id = 2
            n_jd = 1
            n_hm = 3
            n_kd = 4
            n_ao = 5
            
            hm_factor = 1.0d0
            ridisp = 0.d0
            rjdisp = 0.d0
            rkdisp = 0.d0

            a_vec(n_hm) = hm_factor
            a_vec(n_id) = ridisp
            a_vec(n_jd) = rjdisp
            a_vec(n_kd) = rkdisp
            a_vec(n_ao) = aod

            dstep(:) = 0d0 ! initialize
            dstep(n_hm) = 0.2d0 
            dstep(n_id) = 1d0 
            dstep(n_jd) = 1d0 
            dstep(n_kd) = 1d0 
            dstep(n_ao) = 1d-2 

            dstep_gran(n_hm) = 0d0
            dstep_gran(n_id) = 1d0
            dstep_gran(n_jd) = 1d0
            dstep_gran(n_kd) = 1d0
            dstep_gran(n_ao) = 0d0

            init_optmiz = 0
            iexit_optimize = 0

        else ! standard run
            nloops = 1

        endif

        do iloop = 1,nloops ! optimization loop
          if(iloop .gt. 1)then ! optimization case
                write(6,41)iloop,hm_factor,ridisp,rjdisp    
41              format(' Modify fields for optimization ',i4,3e13.5)
                
                clwc_3d(:,:,:) = clwc_3d(:,:,:) * hm_factor
                cice_3d(:,:,:) = cice_3d(:,:,:) * hm_factor

                ido = nint(ridisp)
                jdo = nint(rjdisp)

                idl1 = max(1+ido,1)
                idh1 = min(NX_L+ido,NX_L)

                idl2 = max(idl1-ido,1)
                idh2 = min(idh1-ido,NX_L)
                    
                jdl1 = max(1+jdo,1)
                jdh1 = min(NY_L+jdo,NY_L)

                jdl2 = max(jdl1-jdo,1)
                jdh2 = min(jdh1-jdo,NY_L)

                write(6,42)ido,idl1,idh1,idl2,idh2
42              format(  ' ido,idl1,idh1,idl2,idh2',5i6)             
                write(6,43)jdo,jdl1,jdh1,jdl2,jdh2
43              format(  ' jdo,jdl1,jdh1,jdl2,jdh2',5i6)             

                clwc_3d(idl1:idh1,jdl1:jdh1,:) =
     1          clwc_3d(idl2:idh2,jdl2:jdh2,:)

                cice_3d(idl1:idh1,jdl1:jdh1,:) =
     1          cice_3d(idl2:idh2,jdl2:jdh2,:)

!               Get Aersol Extinction Coefficient (3D field)
                if(n_ao .le. nv)then
                  call get_aod_3d(pres_3d,heights_3d,topo,NX_L,NY_L,NZ_L
     1                   ,aod,aod_ref,i_aero_synplume,i_aero_1d,aod_3d)
                endif
          endif                
      
          do iloc = 1,nloc
            ri_obs = xsound(iloc)
            rj_obs = ysound(iloc)

!           if(trim(c_model) .ne. 'hrrr_smoke')then ! Keep inside model domain
            if(.true.)then
              ri_obs = min(max(ri_obs,1.),float(NX_L))
              rj_obs = min(max(rj_obs,1.),float(NY_L))
!           else
!               ri_obs = 477.
!               rj_obs = 881.             
            endif

            i_obs = nint(ri_obs)
            j_obs = nint(rj_obs)

            rlat = lat(i_obs,j_obs)
            rlon = lon(i_obs,j_obs)

            rlat = fsoundlat(iloc)
            rlon = fsoundlon(iloc)

            write(6,*)
            write(6,*)' iloc = ',iloc,rlat,rlon

            if(iloop .eq. 1)then
              write(6,*)' Enter minalt,maxalt (e.g. 0,90 0,180)'
              read(lun,*)minalt_a(iloc),maxalt_a(iloc)           
            endif
            minalt = minalt_a(iloc)
            maxalt = maxalt_a(iloc)

            if(iloop .eq. 1)then
              write(6,*)' Enter minazi,maxazi'                     
              read(lun,*)minazi_a(iloc),maxazi_a(iloc)           
            endif
            minazi = minazi_a(iloc)
            maxazi = maxazi_a(iloc)

            if(iloop .eq. 1)then
              write(6,*)' Enter alt_scale,azi_scale'               
              read(lun,*)alt_scale_a(iloc),azi_scale_a(iloc)     
            endif
            alt_scale = alt_scale_a(iloc)     
            azi_scale = azi_scale_a(iloc)     

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
            allocate(camera_rgbf(nc,minalt:maxalt,minazi:maxazi))
            allocate(dist_2_topo(minalt:maxalt,minazi:maxazi))
            allocate(sky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi))
            allocate(isky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi))

            if(.false.)then
                topo_sfc = topo(i_obs,j_obs)
            else
                call bilinear_laps(ri_obs,rj_obs,NX_L,NY_L,topo
     1                 ,topo_sfc)
            endif
 
            write(6,*)' observer ri/rj/topo_sfc ',ri_obs,rj_obs
     1                                             ,topo_sfc
            write(6,22)topo_albedo_2d(:,i_obs,j_obs)
22          format('  albedo RGB of observer ',3f9.3)

            write(6,*)' land use / frac of observer: '
     1                 ,land_use(i_obs,j_obs),land_frac(i_obs,j_obs)

            write(6,*)' array of land frac'
            call ascii_map(i_obs,j_obs,NX_L,NY_L,land_frac,10,1)

            write(6,*)' snow cover of observer: '
     1                 ,snow_cover(i_obs,j_obs)  

            write(6,*)' solar alt/az (2d array)'       
     1                ,sol_alt_2d(i_obs,j_obs),sol_azi_2d(i_obs,j_obs) 

!           Calculate solar position for all-sky point
            if(.true.)then ! test this again to allow fractional gridpoints?
                call solar_position(fsoundlat(iloc),fsoundlon(iloc)
     1                             ,i4time_solar,solar_alt     
     1                             ,solar_dec,solar_ha)
                call equ_to_altaz_d(solar_dec,solar_ha,fsoundlat(iloc)
     1                             ,altdum,solar_az)               
                if(solar_az .lt. 0.)solar_az = solar_az + 360.
                  solar_lat = solar_dec
                  solar_lon = soundlon(iloc) - solar_ha
            else ! ensure consistency between both solar positions
                solar_alt = sol_alt_2d(i_obs,j_obs)
                solar_az = sol_azi_2d(i_obs,j_obs)              
            endif

            write(6,*)' solar alt/az (observer)',solar_alt,solar_az

            hdist_loc = sqrt((rlat-rlat_last)**2 + (rlon-rlon_last)**2)
!           if(solar_alt .gt. 4.0)then
!                 hdist_loc_thr = 0.3
!           else
!                 hdist_loc_thr = 0.0
!           endif

            thr1 = ((grid_spacing_m / 10000.) / rpd) * sind(solar_alt) 
            thr2 = solar_alt * 0.1
            thr3 = solar_alt - 3.0
            hdist_loc_thr = max(min(thr1,thr2,thr3),0.0)

            if(hdist_loc .gt. hdist_loc_thr  .OR. 
     1         i4time_solar .ne. i4time_last .OR.   
     1         mode_cloud_mask .eq. 5            )then
                newloc = 1
                rlat_last = rlat; rlon_last = rlon
                i4time_last = i4time_solar
            else
                newloc = 0
            endif

            write(6,*)' i4time_last/i4time_solar = '
     1                 ,i4time_last,i4time_solar
            write(6,23)hdist_loc,hdist_loc_thr,newloc
23          format('  hdist_loc/hdist_loc_thr/newloc = ',2f9.3,i3)

            write(6,*)' call sun_moon at observer sfc grid point '
     1                                                    ,i_obs,j_obs
            idebug = 2
            call sun_moon(i4time_solar,lat,lon,NX_L,NY_L,i_obs,j_obs   ! I
     1                   ,alm,azm                                      ! O
     1                   ,idebug,0.,earth_radius                       ! I
     1                   ,elgms,moon_mag,rmn                           ! O 
     1                   ,geo_dec,geo_ra,geo_sublon,geo_dist           ! O
     1                   ,emag,eobsf,eobsl)                            ! O

            if(elgms .lt. 1.4)then
              write(6,24)emag,eobsf,eobsl
 24           format(' NOTE: Solar Eclipse Conditions: mag/obsc = '
     1              ,f9.6,2f11.6)
              l_solar_eclipse = .true.
            elseif(elgms .lt. 0.6)then
              write(6,*)' NOTE: Possible Solar Eclipse Conditions'
              l_solar_eclipse = .true.
            else
              l_solar_eclipse = .false.
            endif

            write(6,*)' call sun_moon at observer sfc grid point ',i_obs
     1                                                            ,j_obs
            call sun_moon(i4time_solar,lat,lon,NX_L,NY_L,i_obs,j_obs   ! I
     1                   ,alm,azm                                      ! O
     1                   ,idebug,htagl(iloc),earth_radius              ! I
     1                   ,elgms,moon_mag,rmn                           ! O 
     1                   ,geo_dec,geo_ra,geo_sublon,geo_dist           ! O
     1                   ,emag,eobsf,eobsl)                            ! O
            write(6,25)alm,azm,elgms,moon_mag,rmn
 25         format('  alt/az/elg/mnmag/rmn = ',2f8.2,f9.4,f8.2,f9.6)

!           Consider passing 'topo_flag' into 'sun_moon' to consider either
!           solar or lunar eclipses
!           http://www.jgisen.de/eclipse

!           'emag' is solar eclipse magnitude for the observer
!           'eobsf' is observer obscuration
!           'eobsl' includes limb darkening for observer obscuration

            if(maxval(sol_alt_2d) .gt. 0. .and. 
     1       l_solar_eclipse .eqv. .true.)then
              write(6,*)' Calculate gridded sfc eclipse obscuration'
              elgmin = 9999.
              idebug = 0
              do i = 1,NX_L
              do j = 1,NY_L
                call sun_moon(i4time_solar,lat,lon,NX_L,NY_L,i,j         ! I
     1                       ,almgrd,azmgrd                              ! O
     1                       ,idebug,0.,earth_radius                     ! I
     1                       ,elggrd,grdmoon_mag,rmn                     ! O
     1                       ,geo_dec,geo_ra,geo_sublon,geo_dist         ! O
     1                       ,solar_eclipse_magnitude,eobsf,eobsc(i,j))  ! O

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
     1                  ' at',i,j,lat(imxecl,jmxecl),lon(imxecl,jmxecl)

            else
              eobsc = 0.
            endif

!           l_solar_eclipse = .false. ! test

!           alm = -90.          ! Test for disabling
!           moon_mag = -4.0    ! Test for disabling
            moon_mag_thr = -6.0

            moon_alt_2d = alm
            moon_azi_2d = azm

!           Get alt_a_roll and azi_a_roll arrays
            do i = minalt,maxalt
              call get_val(i,minalt,alt_scale,altobj)
              alt_a_roll(i,:) = altobj
              if(altobj .lt. -90.)then
                write(6,*)' ERROR: altobj < -90.',i,altobj
                return
              endif
              if(i .eq. minalt .or.
     1                 (i .eq. maxalt .and. maxalt .gt. 0))then
                if(altobj .ne. nint(altobj))then
                  write(6,*)' ERROR: non-integer altitude bound'
     1                     ,i,alt_scale,altobj
                  return
                endif
              endif
            enddo 

            do j = minazi,maxazi
              call get_val(j,minazi,azi_scale,aziobj)
              azi_a_roll(:,j) = aziobj
              if(j .eq. minazi .or. j .eq. maxazi)then
                if(aziobj .ne. nint(aziobj))then
                  write(6,*)' ERROR: non-integer azimuth bound'
     1                     ,j,azi_scale,aziobj
                  return
                endif
              endif
            enddo

            write(6,*)' alt range is ',alt_a_roll(minalt,minazi)
     1                                ,alt_a_roll(maxalt,minazi)
            write(6,*)' azi range is ',azi_a_roll(minalt,minazi)
     1                                ,azi_a_roll(minalt,maxazi)

            ilun = ilun + 1
            write(clun,34)iloc
34          format(i3.3)

            allocate(aod_ill_opac(minalt:maxalt,minazi:maxazi))
            allocate(aod_ill_opac_potl(minalt:maxalt,minazi:maxazi))

            exposure = exposure_a(iloc) ! density

            if(l_water_world)then ! water world experimental simulation
              land_frac = 0.
              snow_cover = 0.
              topo = 0.
              topo_sfc = 0.
              topo_albedo_2d(1,:,:) = .003
              topo_albedo_2d(2,:,:) = .007
              topo_albedo_2d(3,:,:) = .028
            endif

!           Setup optimization with single site, for multiple sites we can move
!           this up before line 1156 (iloc loop)

            I4_elapsed = ishow_timer()

            write(6,*)' call calc_allsky at i4time_solar:',i4time_solar
     1               ,iloop
            call calc_allsky(i4time_solar,exposure ! ,clwc_3d,cice_3d
!    1                     ,heights_3d                              ! I
!    1                     ,rain_3d,snow_3d                         ! I
!    1                     ,pres_3d,aod_3d                          ! I
     1                     ,topo_sfc,topo,swi_2d,pw_2d              ! I
     1                     ,topo_albedo_2d,land_frac,snow_cover     ! I
     1                     ,htagl(iloc)                             ! I
     1                     ,aod_ref                                 ! I
     1                     ,NX_L,NY_L,NZ_L,newloc                   ! I
     1                     ,ri_obs,rj_obs                           ! I
     1                     ,alt_a_roll,azi_a_roll                   ! I
     1                     ,sol_alt_2d,sol_azi_2d                   ! I
     1                     ,solar_alt,solar_az                      ! I
     1                     ,solar_lat,solar_lon,r_au                ! I
     1                     ,alt_norm                                ! I
     1                     ,moon_alt_2d,moon_azi_2d,alm,azm         ! I
     1                     ,moon_mag,moon_mag_thr,elgms             ! I
     1                     ,l_solar_eclipse,eobsc,emag              ! I
     1                     ,rlat,rlon,lat,lon                       ! I
     1                     ,minalt,maxalt,minazi,maxazi,nsp         ! I
     1                     ,ni_cyl,nj_cyl                           ! O
     1                     ,alt_scale,azi_scale                     ! I
     1                     ,grid_spacing_m,r_missing_data           ! I
     1                     ,l_binary,l_terrain_following            ! I
     1                     ,mode_cloud_mask,camera_cloud_mask       ! I
     1                     ,iloop                                   ! I
     1                     ,cloud_od,dist_2_topo                    ! O
     1                     ,camera_rgbf                             ! O
     1                     ,sky_rgb_cyl,correlation(:,iloc),istatus)! O
            if(istatus .ne. 1)then
              write(6,*)' Error istatus returned from calc_allsky'
              return
            endif

!           For multiple sites we can move this outside 'iloc' loop past
!           line 1725

            if(mode_cloud_mask .eq. 4  .or. mode_cloud_mask .eq. 5)then
              avecorr = sum(correlation(:,iloc))/float(nc)
              write(6,46)correlation(:,iloc),avecorr
46            format('correlation is ',3f9.6,2x,f9.6,' incremental')
              write(6,*)' camera_rgbf checksum = '
     1                  ,sum(min(camera_rgbf,255.))
            else
              correlation = 0.
              avecorr = 0. 
            endif

            if(nloops .gt. 1)then
              write(clun_loop,48)clun,iloop
48            format(a3,'_',i3.3)
            else
              clun_loop = clun
            endif

            if(.true.)then ! write labels and images        
              write(6,*)' end of subroutine call block - write labels'

!             Write time label
              open(53,file='label.'//trim(clun_loop),status='unknown')
              write(53,*)a9time
              write(53,*)a24time(1:17)
              close(53)

              if(soundlat(iloc) .gt. 0.)then
                c1_lat = 'N'
              else
                c1_lat = 'S'
              endif

              if(soundlon(iloc) .gt. 0.)then
                c1_lon = 'E'
              else
                c1_lon = 'W'
              endif

!             Write lat/lon and other info for label
              open(54,file='label2.'//trim(clun_loop),status='unknown')
              write(54,54)soundlat(iloc),c1_lat,
     1                    abs(soundlon(iloc)),c1_lon,
     1                    minalt,maxalt,minazi,maxazi,ni_cyl,nj_cyl,
     1                    solar_alt,solar_az,alt_scale,azi_scale,
     1                    ni_polar,nj_polar,
     1                    htagl(iloc)
 54           format(2(f8.2,a1)/6i8/2f8.2,2f7.2/2i6/f10.0)

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

              write(54,58)correlation(:,iloc),avecorr
 58           format(4f10.3)

              close(54)

              horz_dep = horz_depf(htagl(iloc),earth_radius)
              solalt_limb_true = solar_alt + horz_dep
              write(6,*)' solalt_limb_true is',solalt_limb_true

                if(solalt_limb_true .lt. -10.0)then
                   dither = 2.0
                else
                   dither = 0.0
                endif

              if(l_cyl .eqv. .true.)then
!               Write all sky for cyl

                if(dither .gt. 0.)then

                   I4_elapsed = ishow_timer()

                   write(6,*)' writing dithered image',dither

                   do ialt = minalt,maxalt
                   do iaz = minazi,maxazi
                   do ic = 1,nc
                      call random_number(arand)
                      skydelt = (-1.0 + (2.0 * arand)) * dither
                    isky_rgb_cyl(ic,ialt,iaz) = sky_rgb_cyl(ic,ialt,iaz)
     1                                         + skydelt
                   enddo ! ic

                   isky_rgb_cyl(:,ialt,iaz) =
     1                        min(max(isky_rgb_cyl(:,ialt,iaz),0),255)

                   if(iaz .eq. 3*(minazi+maxazi)/4)then
                         write(6,*)' arand is ',arand
                         write(6,61)ialt,alt_a_roll(ialt,iaz),arand
     1                             ,skydelt,sky_rgb_cyl(2,ialt,iaz)
     1                             ,isky_rgb_cyl(:,ialt,iaz)
61                       format(' dither ',i4,f6.2,2f7.2,f8.2,3i5)
                      endif

                    enddo ! iaz
                    enddo ! ialt 
                else
                    isky_rgb_cyl = sky_rgb_cyl   
                endif                                       

                npts = 3*(maxalt-minalt+1)*(maxazi-minazi+1)
!               write(6,*)' Write all sky cyl text file ',npts
!               open(55,file='allsky_rgb_cyl.'//trim(clun_loop),status='unknown')
!               write(55,*)isky_rgb_cyl           
!               close(55)
                write(6,*)' Write all sky cyl ppm file ',trim(clun_loop)
                call writeppm3Matrix(
     1                isky_rgb_cyl(0,:,:),isky_rgb_cyl(1,:,:)
     1               ,isky_rgb_cyl(2,:,:)
     1               ,'allsky_rgb_cyl_'//trim(clun_loop))     
                do iaz = minazi,maxazi,40
                  write(6,*)'iaz,cyl(maxalt/1,iaz)',iaz
     1                       ,isky_rgb_cyl(:,maxalt/1,iaz)
                enddo ! iaz
              endif

              if(l_polar .eqv. .true.)then

                allocate(sky_rgb_polar(0:2,iplo:iphi,jplo:jphi))
                allocate(isky_rgb_polar(0:2,iplo:iphi,jplo:jphi))

!               Reproject sky_rgb array from cyl to polar    
                do iaz = minazi,maxazi,20
                  write(6,*)'iaz,cyl((maxalt+minalt)/2,iaz)',iaz
     1                       ,sky_rgb_cyl(1,(maxalt+minalt)/2,iaz)
                enddo ! iaz

!               write(6,*)' Call cyl_to_polar with sky rgb data'

                if(htagl(iloc) .le. 20.1 .and. 
     1                  (maxalt*2 + minalt) .gt. 0)then
                  polat = +90. ! +/-90 for zenith or nadir at center of plot
                else
                  polat = -90.
                endif
                polat_tmp = getenv_real('POLAT',r_missing_data)
                if(polat_tmp .ne. r_missing_data)polat = polat_tmp

!               if(ipolar_sizeparm .ge. 3)then
                if(htagl(iloc) .gt. earth_radius*2.5)then
                  pomag = htagl(iloc) / (earth_radius*0.68)
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

                pox = 0.
                poy = 0.

                if(trim(c_model) .eq. 'hrrr_smoke' .and.
     1             htagl(iloc) .ge. 10000e3)then
                  if(abs(rlat) .eq. 0.)then
                    pomag = pomag * 2.4
                    poy = +0.6 ! increase moves down
                    pox = +0.1 ! increase moves right
                  else
                    pomag = pomag * 8.                 
                  endif
                endif

                rotew = getenv_real('ROTEW',r_missing_data) 
 
                if(htagl(iloc) .eq. 400001.)then
                  rotew = +70. ! positive rotation decreases altitude in the south
                  rotz = +135. ! positive rotation moves counterclockwise (left)
                  pomag = pomag * 1.8                 
                  poy = -0.1   ! increase moves down
                elseif(rotew .ne. r_missing_data)then
                  continue
                else
                  rotew = 0.
                  rotz = 0.
                endif

                pomag_tmp = getenv_real('POMAG',r_missing_data)
                if(pomag_tmp .ne. r_missing_data)pomag = pomag_tmp

                pox_tmp = getenv_real('POX',r_missing_data)
                if(pox_tmp .ne. r_missing_data)pox = pox_tmp

                poy_tmp = getenv_real('POY',r_missing_data)
                if(poy_tmp .ne. r_missing_data)poy = poy_tmp

                rotz_tmp  = getenv_real('ROTZ',r_missing_data)
                if(rotz_tmp .ne. r_missing_data)rotz = rotz_tmp

                if(alt_scale .lt. .005)then ! HIGHFLAT regional model
                  pomag = pomag * 5.
                endif                

                write(6,*)'htrat/pomag',htagl(iloc)/earth_radius,pomag

                do ic = 0,nc-1
                  write(6,*)' Call cyl_to_polar with sky rgb data',ic

                  call cyl_to_polar(sky_rgb_cyl(ic,:,:)
     1                             ,sky_rgb_polar(ic,:,:)
     1                             ,minalt,maxalt,minazi,maxazi
     1                             ,alt_scale,azi_scale,polat,pomag
     1                             ,pox,poy,alt_a_polar,azi_a_polar
     1                             ,rotew,rotz
     1                             ,iplo,iphi,jplo,jphi
     1                             ,ni_polar,nj_polar)
                enddo ! ic

!               Write all sky for polar
                where(sky_rgb_polar .eq. r_missing_data)
                  sky_rgb_polar = 0.
                  isky_rgb_polar = 0
                endwhere

                if(dither .gt. 0.)then

                  I4_elapsed = ishow_timer()

                  write(6,*)' writing dithered image',dither

                  do i = iplo,iphi
                  do j = jplo,jphi
                    call random_number(arand)
                    skydelt = (-1.0 + (2.0 * arand)) * dither
                    do ic = 1,nc
                      isky_rgb_polar(ic,i,j) = 
     1                             sky_rgb_polar(ic,i,j) + skydelt
                    enddo ! ic

                    isky_rgb_polar(:,ialt,iaz) =
     1                        min(max(isky_rgb_polar(:,i,j),0),255)

                    if(j .eq. (jplo+jphi)/2)then
                         write(6,61)i,alt_a_roll(i,j),arand
     1                             ,skydelt,sky_rgb_polar(2,i,j)
     1                             ,isky_rgb_polar(:,i,j)
                    endif

                  enddo ! j
                  enddo ! i 
                else
                isky_rgb_polar = sky_rgb_polar
                endif                                       

                write(6,*)' max polar 1 is ',maxval(sky_rgb_polar)
     1                                      ,maxval(isky_rgb_polar)
                write(6,*)' min polar 1 is ',minval(sky_rgb_polar)
     1                                      ,minval(isky_rgb_polar)

                if(minval(isky_rgb_polar) .lt. 0)then
                  write(6,*)' WARNING: maxval isky_rgb_polar < 0'
              
                  where(isky_rgb_polar .lt. 0)
                        isky_rgb_polar = 0
                  endwhere
                  write(6,*)' max polar 2 is ',maxval(sky_rgb_polar)
     1                                        ,maxval(isky_rgb_polar)
                  write(6,*)' min polar 2 is ',minval(sky_rgb_polar)
     1                                        ,minval(isky_rgb_polar)
                endif

!               write(6,*)' ipolar array',
!    1              isky_rgb_polar(2,ni_polar/2,1:nj_polar)
                nip_crop = iphi-iplo+1
                njp_crop = jphi-jplo+1
                npts = 3*nip_crop*njp_crop
!               write(6,*)' Write all sky polar text file'
!    1                    ,isky_rgb_polar(:,255,255),npts
!               open(54,file='allsky_rgb_polar.'//trim(clun_loop),status='unknown')
!               write(54,*)isky_rgb_polar
!               close(54)
                write(6,*)' Write all sky polar ppm file '
     1                   ,trim(clun_loop)
                call writeppm3Matrix(
     1                  isky_rgb_polar(0,:,:),isky_rgb_polar(1,:,:)
     1                 ,isky_rgb_polar(2,:,:)
     1                 ,'allsky_rgb_polar_'//trim(clun_loop))

                deallocate(sky_rgb_polar)
                deallocate(isky_rgb_polar)

              endif ! l_polar

            endif ! write labels and images

!           For multiple sites we can move this outside the 'iloc' loop past
!           line 1725

            deallocate(aod_ill_opac)
            deallocate(aod_ill_opac_potl)
            deallocate(alt_a_roll)
            deallocate(azi_a_roll)
            deallocate(cloud_od)
            deallocate(camera_cloud_mask)
            deallocate(camera_rgbf)
            deallocate(dist_2_topo)

 900        continue
  
1000        continue

            write(6,*)' End of plot_allsky for iloc...',iloc
            write(6,*)

            I4_elapsed = ishow_timer()

            deallocate(sky_rgb_cyl)
            deallocate(isky_rgb_cyl)

          enddo ! iloc

          if(mode_cloud_mask .eq. 4  .or. mode_cloud_mask .eq. 5)then
              avecorr = sum(correlation(:,1:nloc))/float(nc*nloc)
              do iloc = 1,nloc
                write(6,1842)correlation(:,iloc),avecorr
1842            format('correlation is ',3f9.6,2x,f9.6,' allcams')
              enddo ! iloc
          else
              correlation = 0.
              avecorr = 0. 
          endif

          if(mode_cloud_mask .eq. 5)then
              write(6,*)'Update a vector and corresponding parameters'
              a_last(:) = a_vec(:)

              f_merit = 1.0 - avecorr
              f_merit = f_merit + 0.1d0 * (1d0 - a_vec(n_hm))**2
!             f_merit = f_merit + .1d0 * a_vec(n_id)**2
!             f_merit = f_merit + .1d0 * a_vec(n_jd)**2

              if(iexit_optimize .eq. 1)then              
                  write(6,*)' Signal to stop optimization'
                  write(6,1891)a_vec(1:nv),f_merit,iloop
1891              format(' exit optimize: a_vec/f_merit'
     1                  ,3f10.5,5x,f10.6,i4)
                  goto 1900
              endif

              call optimize_wrapper(a_vec,dstep,dstep_gran,nv,f_merit
     1                             ,init_optmiz,iexit_optimize)

              write(6,1892)a_vec(1:nv),f_merit
1892          format(' ret optimize:  a_vec/f_merit',3e13.5,4x,f9.5)

              hm_factor = a_vec(n_hm) / a_last(n_hm)
              ridisp = nint(a_vec(n_id) - a_last(n_id))
              rjdisp = nint(a_vec(n_jd) - a_last(n_jd))
              rkdisp = nint(a_vec(n_kd) - a_last(n_kd))
              aod = a_vec(n_ao)
          endif

        enddo ! iloop for optimize

1900    call dealloc_allsky

        write(6,*)
        write(6,*)' End of plot_allsky...'
        write(6,*)

        return
        end

        subroutine ascii_map(i,j,ni,nj,array,ihw,nd)

        real*4 array(ni,nj)

        imin = max(i-ihw,1)
        imax = min(i+ihw,ni)

        jmin = max(j-ihw,1)
        jmax = min(j+ihw,nj)

        do jr = jmax,jmin,-1
           write(6,1)nint(array(imin:imax,jr))
 1         format(100i1)
        enddo ! jr
        
        return
        end
      
        function getenv_real(c_env,r_missing_data)

        character(*) c_env
        character*20 c_value

        call GETENV(c_env,c_value)

        call s_len(c_value,lenv)
      
        if(lenv .gt. 0)then
            read(c_value,*)getenv_real
        else
            getenv_real = r_missing_data 
        endif
            
        write(6,*)' getenv_real ',c_env,getenv_real

        return
        end
