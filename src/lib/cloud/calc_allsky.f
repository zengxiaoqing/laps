
          subroutine calc_allsky(i4time_solar ! ,clwc_3d,cice_3d
!    1                     ,heights_3d                           ! I
!    1                     ,rain_3d,snow_3d                      ! I
!    1                     ,pres_3d,aod_3d                       ! I
     1                     ,topo_sfc,topo,swi_2d                 ! I
     1                     ,topo_albedo_2d                       ! I
     1                     ,htagl                                ! I
     1                     ,aod_ref                              ! I
     1                     ,NX_L,NY_L,NZ_L,isound,jsound,newloc  ! I
     1                     ,ri_obs,rj_obs                        ! I
     1                     ,alt_a_roll,azi_a_roll                ! I
     1                     ,sol_alt_2d,sol_azi_2d                ! I
     1                     ,solar_alt,solar_az                   ! I
     1                     ,alt_norm                             ! I
     1                     ,moon_alt_2d,moon_azi_2d,alm,azm      ! I
     1                     ,moon_mag,moon_mag_thr                ! I
     1                     ,l_solar_eclipse,rlat,rlon,lat,lon    ! I
     1                     ,minalt,maxalt,minazi,maxazi,nc,nsp   ! I
     1                     ,ni_cyl,nj_cyl                        ! I
     1                     ,alt_scale,azi_scale                  ! I
     1                     ,grid_spacing_m,r_missing_data        ! I
     1                     ,twi_0,l_binary                       ! I
     1                     ,sky_rgb_cyl)                         ! O

        use mem_allsky

        addlogs(x,y) = log10(10.**x + 10.**y)

!       Input arrays
!       real clwc_3d(NX_L,NY_L,NZ_L)      ! Control Variable
!       real cice_3d(NX_L,NY_L,NZ_L)      ! Control Variable
!       real heights_3d(NX_L,NY_L,NZ_L)   ! Control Variable
!       real rain_3d(NX_L,NY_L,NZ_L)      ! Control Variable
!       real snow_3d(NX_L,NY_L,NZ_L)      ! Control Variable
!       real pres_3d(NX_L,NY_L,NZ_L)      ! Control Variable
!       real aod_3d(NX_L,NY_L,NZ_L)       
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

!       Local arrays (e.g. outputs from get_cloud_rays)
!       real transm_3d(NX_L,NY_L,NZ_L)
!       real transm_4d(NX_L,NY_L,NZ_L,nc) 

        real r_cloud_3d(minalt:maxalt,minazi:maxazi)
        real cloud_od(minalt:maxalt,minazi:maxazi)
        real cloud_od_sp(minalt:maxalt,minazi:maxazi,nsp)
        real elong_roll(minalt:maxalt,minazi:maxazi)
        real airmass_2_cloud_3d(minalt:maxalt,minazi:maxazi)
        real airmass_2_topo_3d(minalt:maxalt,minazi:maxazi)
        real topo_swi(minalt:maxalt,minazi:maxazi)
        real topo_albedo(nc,minalt:maxalt,minazi:maxazi)
        real topo_ri(minalt:maxalt,minazi:maxazi)
        real topo_rj(minalt:maxalt,minazi:maxazi)
        real trace_ri(minalt:maxalt,minazi:maxazi)
        real trace_rj(minalt:maxalt,minazi:maxazi)
        real topo_solalt(minalt:maxalt,minazi:maxazi)
        real trace_solalt(minalt:maxalt,minazi:maxazi)
        real gtic(nc,minalt:maxalt,minazi:maxazi)
        real aod_2_cloud(minalt:maxalt,minazi:maxazi)
        real aod_2_topo(minalt:maxalt,minazi:maxazi)
        real dist_2_topo(minalt:maxalt,minazi:maxazi)
        real aod_ill(minalt:maxalt,minazi:maxazi)
        real aod_ill_dir(minalt:maxalt,minazi:maxazi)
        real aod_tot(minalt:maxalt,minazi:maxazi)
        real r_cloud_trans(minalt:maxalt,minazi:maxazi)
        real cloud_rad_c(nc,minalt:maxalt,minazi:maxazi)
        real cloud_rad_w(minalt:maxalt,minazi:maxazi)
        real clear_rad_c(nc,minalt:maxalt,minazi:maxazi)
        real clear_radf_c(nc,minalt:maxalt,minazi:maxazi)

!       Other Local Arrays
        real blog_moon_roll(minalt:maxalt,minazi:maxazi)
        real blog_sun_roll(minalt:maxalt,minazi:maxazi)
        real glow_stars(nc,minalt:maxalt,minazi:maxazi)

        real moon_mag,moon_mag_thr
        logical l_solar_eclipse, l_binary

        write(6,*)' subroutine calc_allsky...'

        write(6,*)' call get_cloud_rays...'

!         Get line of sight from isound/jsound
          call get_cloud_rays(i4time_solar,clwc_3d,cice_3d
     1                     ,heights_3d                           ! I
     1                     ,rain_3d,snow_3d                      ! I
     1                     ,pres_3d,aod_3d,topo_sfc,topo         ! I
     1                     ,topo_albedo_2d                       ! I
     1                     ,topo_swi,topo_albedo,gtic            ! O
     1                     ,topo_ri,topo_rj                      ! O
     1                     ,trace_ri,trace_rj                    ! O
!    1                     ,ghi_2d,dhi_2d                        ! O
     1                     ,aod_vrt,aod_2_cloud,aod_2_topo       ! O
     1                     ,dist_2_topo                          ! O
     1                     ,aod_ill,aod_ill_dir                  ! O
     1                     ,aod_tot,transm_obs,obs_glow_zen      ! O
     1                     ,transm_3d,transm_4d                  ! O
     1                     ,r_cloud_3d,cloud_od,cloud_od_sp      ! O
     1                     ,r_cloud_trans,cloud_rad_c,cloud_rad_w! O
     1                     ,clear_rad_c,clear_radf_c,patm        ! O
     1                     ,airmass_2_cloud_3d,airmass_2_topo_3d ! O
     1                     ,htmsl,horz_dep,twi_0                 ! O
     1                     ,solalt_limb_true                     ! O
     1                     ,htagl                                ! I
     1                     ,aod_ref                              ! I
     1                     ,NX_L,NY_L,NZ_L,isound,jsound,newloc  ! I
     1                     ,ri_obs,rj_obs                        ! I
     1                     ,alt_a_roll,azi_a_roll                ! I
     1                     ,sol_alt_2d,sol_azi_2d                ! I
     1                     ,alt_norm                             ! I
     1                     ,moon_alt_2d,moon_azi_2d              ! I
     1                     ,moon_mag,moon_mag_thr                ! I
     1                     ,l_solar_eclipse,rlat,rlon,lat,lon    ! I
     1                     ,minalt,maxalt,minazi,maxazi          ! I
     1                     ,alt_scale,azi_scale                  ! I
     1                     ,grid_spacing_m,r_missing_data)       ! I

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
     1                       ,alm,azm,moon_mag,diam_deg,horz_dep
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
     1                     ,solar_alt,solar_az,s_mag,diam_deg,horz_dep
     1                     ,blog_sun_roll)
          write(6,*)' range of blog_sun_roll is',
     1        minval(blog_sun_roll),maxval(blog_sun_roll),diam_deg

!         if(solar_alt .ge. 0.)then
!         if(.true.)then
          I4_elapsed = ishow_timer()
!         write(6,*)' call get_skyglow_cyl for sun or moon...'

!         Get all sky for cyl   
          ni_cyl = maxalt - minalt + 1
          nj_cyl = maxazi - minazi + 1

          I4_elapsed = ishow_timer()

          if(solar_alt .lt. 0. .OR. l_solar_eclipse .eqv. .true.)then
              write(6,*)' call get_starglow with cyl data'
              call get_starglow(i4time_solar,alt_a_roll,azi_a_roll       ! I
     1                     ,minalt,maxalt,minazi,maxazi                  ! I
     1                     ,rlat,rlon,alt_scale,azi_scale,horz_dep       ! I
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

          do j = minazi,maxazi
          do i = minalt,maxalt
                trace_solalt(i,j) = solar_alt
                topo_solalt(i,j) = solar_alt

                itrace = nint(trace_ri(i,j))
                jtrace = nint(trace_rj(i,j))
                if(htmsl .gt. 100000. .and. alt_a_roll(i,j) .lt. 0.)then
                    if(itrace .ge. 1 .and. itrace .le. NX_L .and.
     1                 jtrace .ge. 1 .and. jtrace .le. NY_L)then
                        trace_solalt(i,j) = sol_alt_2d(itrace,jtrace)
                    endif
                endif

                itopo = nint(topo_ri(i,j))
                jtopo = nint(topo_rj(i,j))
                if(itopo .ge. 1 .and. itopo .le. NX_L .and.
     1             jtopo .ge. 1 .and. jtopo .le. NY_L)then
                    topo_solalt(i,j) = sol_alt_2d(itopo,jtopo)
                endif

          enddo ! i
          enddo ! j

          if(l_binary .eqv. .false.)then
              write(6,*)' call get_sky_rgb with cyl data'
              if(htagl .eq. 1000e3)then                    ! High alt
                corr1_a = 9.2
              elseif(htagl .eq. 300.)then                  ! BAO
                if(solar_alt .lt. 30.)then
                  corr1_a = 9.2                            ! low sun
                elseif(solar_alt .gt. 60.)then
                  corr1_a = 8.9                            ! high sun
                else
                  corr1_a = 9.2 - (solar_alt-30.)*(0.3/30.)! med sun
                endif
              else
                corr1_a = 9.0
              endif
              if(solar_alt .lt. 0.)corr1_a = 9.26 ! volcanic value

              write(6,*)' corr1 in calc_allsky ',corr1_a

              write(6,*)' call get_sky_rgb with cyl data'
              call get_sky_rgb(r_cloud_3d    ! cloud opacity
     1                    ,cloud_od          ! cloud optical depth
     1                    ,cloud_od_sp,nsp   ! cloud species optical depth
     1                    ,r_cloud_trans     ! cloud solar transmittance
     1                    ,cloud_rad_c       ! cloud solar transmittance / color
     1                    ,cloud_rad_w       ! cloud solar transmittance * rad
     1                    ,clear_rad_c       ! clear sky illumination by sun     
     1                    ,l_solar_eclipse,i4time_solar,rlat,rlon,eobsl
     1                    ,clear_radf_c      ! clear sky frac illumination by sun     
     1                    ,patm,htmsl
!    1                    ,blog_v_roll       ! skyglow
     1                    ,blog_sun_roll     ! sunglow
     1                    ,blog_moon_roll    ! moonglow
     1                    ,glow_stars        ! starglow
     1                    ,aod_vrt,aod_ref 
     1                    ,transm_obs,obs_glow_zen ! observer illumination
     1                    ,ialt_sun,jazi_sun ! sun location
     1                    ,airmass_2_cloud_3d      
     1                    ,airmass_2_topo_3d      
     1                    ,topo_swi,topo_albedo,gtic
     1                    ,topo_albedo_2d(:,isound,jsound)
     1                    ,aod_2_cloud,aod_2_topo,aod_ill,aod_ill_dir
     1                    ,aod_tot
     1                    ,dist_2_topo,topo_solalt,trace_solalt
     1                    ,alt_a_roll,azi_a_roll ! I   
     1                    ,elong_roll    
     1                    ,ni_cyl,nj_cyl,azi_scale  
     1                    ,solar_alt,solar_az
     1                    ,minalt,maxalt,minazi,maxazi                  ! I
     1                    ,twi_0,horz_dep
     1                    ,solalt_limb_true
     1                    ,alm,azm,moon_mag  ! moon alt/az/mag
     1                    ,corr1_a
     1                    ,sky_rgb_cyl)   

          else
              continue ! use cloud_od to drive categorical output

          endif ! l_binary

          return
          end
