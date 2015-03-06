
          subroutine calc_allsky(i4time_solar,clwc_3d,cice_3d
     1                     ,heights_3d                           ! I
     1                     ,rain_3d,snow_3d                      ! I
     1                     ,pres_3d,aod_3d,topo_sfc,topo,swi_2d  ! I
     1                     ,topo_albedo_2d                       ! I
     1                     ,htagl                                ! I
     1                     ,aod_ref                              ! I
     1                     ,NX_L,NY_L,NZ_L,isound,jsound,newloc  ! I
     1                     ,alt_a_roll,azi_a_roll                ! I
     1                     ,sol_alt_2d,sol_azi_2d                ! I
     1                     ,alt_norm                             ! I
     1                     ,moon_alt_2d,moon_azi_2d              ! I
     1                     ,moon_mag,moon_mag_thr                ! I
     1                     ,l_solar_eclipse,rlat,rlon,lat,lon    ! I
     1                     ,minalt,maxalt,minazi,maxazi,nc       ! I
     1                     ,alt_scale,azi_scale                  ! I
     1                     ,grid_spacing_m,r_missing_data        ! I
     1                     ,twi_0                                ! I
     1                     ,sky_rgb_cyl)                         ! O

        addlogs(x,y) = log10(10.**x + 10.**y)

!       Input arrays
        real clwc_3d(NX_L,NY_L,NZ_L)      ! Control Variable
        real cice_3d(NX_L,NY_L,NZ_L)      ! Control Variable
        real heights_3d(NX_L,NY_L,NZ_L)   ! Control Variable
        real rain_3d(NX_L,NY_L,NZ_L)      ! Control Variable
        real snow_3d(NX_L,NY_L,NZ_L)      ! Control Variable
        real pres_3d(NX_L,NY_L,NZ_L)      ! Control Variable
        real aod_3d(NX_L,NY_L,NZ_L)       
        real topo(NX_L,NY_L)
        real swi_2d(NX_L,NY_L)
        real topo_albedo_2d(nc,NX_L,NY_L)
        real alt_a_roll(minalt:maxalt,minazi:maxazi)
        real azi_a_roll(minalt:maxalt,minazi:maxazi)
        real sol_alt_2d(NX_L,NY_L)
        real sol_azi_2d(NX_L,NY_L)
        real alt_norm(NX_L,NY_L)
        real moon_alt_2d(NX_L,NY_L)
        real moon_azi_2d(NX_L,NY_L)
        real lat(NX_L,NY_L)
        real lon(NX_L,NY_L)

!       Output arrays
        real sky_rgb_cyl(0:2,minalt:maxalt,minazi:maxazi) ! Observed Variable

!       Local arrays (i.e. outputs from get_cloud_rays)
        real topo_swi(minalt:maxalt,minazi:maxazi)
        real topo_albedo(nc,minalt:maxalt,minazi:maxazi)
        real ghic(nc,minalt:maxalt,minazi:maxazi)
        real aod_2_cloud(minalt:maxalt,minazi:maxazi)
        real aod_2_topo(minalt:maxalt,minazi:maxazi)
        real dist_2_topo(minalt:maxalt,minazi:maxazi)
        real aod_ill(minalt:maxalt,minazi:maxazi)
        real aod_ill_dir(minalt:maxalt,minazi:maxazi)
        real transm_3d(NX_L,NY_L,NZ_L)
        real transm_4d(NX_L,NY_L,NZ_L,nc) 
        real r_cloud_3d(minalt:maxalt,minazi:maxazi)
        real cloud_od(minalt:maxalt,minazi:maxazi)
        real cloud_od_sp(minalt:maxalt,minazi:maxazi,nsp)
        real r_cloud_trans(minalt:maxalt,minazi:maxazi)
        real cloud_rad_c(nc,minalt:maxalt,minazi:maxazi)
        real cloud_rad_w(minalt:maxalt,minazi:maxazi)
        real clear_rad_c(nc,minalt:maxalt,minazi:maxazi)
        real clear_radf_c(nc,minalt:maxalt,minazi:maxazi)
        real airmass_2_cloud_3d(minalt:maxalt,minazi:maxazi)
        real airmass_2_topo_3d(minalt:maxalt,minazi:maxazi)

!       Other Local Arrays
        real blog_moon_roll(minalt:maxalt,minazi:maxazi)
        real blog_sun_roll(minalt:maxalt,minazi:maxazi)
        real glow_stars(nc,minalt:maxalt,minazi:maxazi)

        data ilun /0/
        character*3 clun
        character*24 a24time

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
     1                     ,NX_L,NY_L,NZ_L,isound,jsound,newloc  ! I
     1                     ,alt_a_roll,azi_a_roll                ! I
     1                     ,sol_alt_2d,sol_azi_2d                ! I
     1                     ,alt_norm                             ! I
     1                     ,moon_alt_2d,moon_azi_2d              ! I
     1                     ,moon_mag,moon_mag_thr,twi_0          ! I
     1                     ,l_solar_eclipse,rlat,rlon,lat,lon    ! I
     1                     ,minalt,maxalt,minazi,maxazi          ! I
     1                     ,alt_scale,azi_scale                  ! I
     1                     ,grid_spacing_m,r_missing_data        ! I
     1                     ,sky_rgb_cyl)                         ! O

          write(6,*)' Return from get_cloud_rays: ',a9time
     1             ,' aod_vrt is ',aod_vrt

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
!         if(.true.)then
          I4_elapsed = ishow_timer()
!         write(6,*)' call get_skyglow_cyl for sun or moon...'


!         Reproject Polar Cloud Plot
!         lunsky = 60 
!         write(lunsky,*)rmaglim_v

!         Get all sky for cyl   
          ni_cyl = maxalt - minalt + 1
          nj_cyl = maxazi - minazi + 1

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
     1                    ,aod_vrt,aod_ref 
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

          return
          end
