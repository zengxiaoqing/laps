         
        subroutine get_cloud_rays(i4time,clwc_3d,cice_3d,heights_3d     ! I
     1                           ,rain_3d,snow_3d                       ! I
     1                           ,pres_3d,topo_sfc,topo_a               ! I
     1                           ,topo_albedo_2d                        ! I
     1                           ,swi_2d                                ! I
     1                           ,topo_gti,topo_albedo                  ! O
     1                           ,gtic,dtic,btic,emic                   ! O
     1                           ,topo_ri,topo_rj                       ! O
     1                           ,trace_ri,trace_rj                     ! O
     1                           ,swi_obs                               ! O
     1                           ,aod_vrt,aod_2_cloud,aod_2_topo        ! O
     1                           ,dist_2_topo                           ! O
     1                           ,aod_ill,aod_ill_dir                   ! O
     1                           ,aod_tot,transm_obs,obs_glow_gnd       ! O
     1                           ,transm_3d,transm_4d                   ! O
     1                           ,r_cloud_3d,cloud_od,cloud_od_sp       ! O
     1                           ,cloud_od_sp_w                         ! O
     1                           ,r_cloud_rad,cloud_rad_c,cloud_rad_w   ! O
     1                           ,clear_rad_c,clear_radf_c,patm         ! O
     1                           ,clear_rad_c_nt                        ! O
     1                           ,airmass_2_cloud_3d,airmass_2_topo_3d  ! O
     1                           ,htstart,horz_dep_d,twi_0              ! O
!    1                           ,solalt_limb_true                      ! O
     1                           ,htagl                                 ! I
!    1                           ,elong                                 ! I
     1                           ,aod,ext_g                             ! I
     1                           ,ni,nj,nk,i,j,newloc,ri_obs,rj_obs     ! I
     1                           ,view_alt,view_az,sol_alt,sol_azi      ! I
     1                           ,alt_norm                              ! I
     1                           ,moon_alt,moon_azi                     ! I
     1                           ,moon_mag,moon_mag_thr                 ! I
     1                           ,l_solar_eclipse,eobsc                 ! I
     1                           ,rlat,rlon,lat,lon                     ! I
     1                           ,minalt,maxalt,minazi,maxazi           ! I
     1                           ,alt_scale,azi_scale                   ! I
     1                           ,l_binary,l_terrain_following          ! I
     1                           ,grid_spacing_m,r_missing_data         ! I
     1                           ,istatus)                              ! O

        use mem_namelist, ONLY: earth_radius,aero_scaleht,redp_lvl
        use mem_allsky, ONLY: aod_3d ! (extinction coefficient)         ! I
        use mem_allsky, ONLY: aod_ill_opac,aod_ill_opac_potl            ! O
        use mem_allsky, ONLY: uprad_4d ! (upward spectral irradiance)   ! L
        use mem_allsky, ONLY: upxrad_3d ! (upward irradiance xcos)      ! L
        use mem_allsky, ONLY: upyrad_3d ! (upward irradiance ycos)      ! L
        use mem_allsky, ONLY: mode_aero_cld
        use cloud_rad

        include 'trigd.inc'
        include 'wa.inc'
        parameter (pi=3.14159265)

!       Statement Functions
        trans(od) = exp(-min(od,80.))
        curvat(hdst,radius) = hdst**2 / (2. * radius)
        resin(x) = sin(x*1.57079)

        real*8 dsdst,dradius_start,daltray,dcurvat2
        dcurvat2(dsdst,dradius_start,daltray) = 
     1      sqrt(dsdst**2 + dradius_start**2 
     1   - (2D0*dsdst*dradius_start*cosd(90D0+daltray))) - dradius_start      
        horz_depf(htmsl,erad) = acosd(erad/(erad+htmsl))
        solocc_f(htmsl,erad,solaltarg) = min(max(
     1          (solaltarg + horz_depf(htmsl,erad))/0.25 ,0.),1.)

!       Airmasses relative to zenith at sea level pressure (true altitude)
        airmass_cosz(cosz) =  
     1    (          1.002432 * cosz**2 + 0.148386  * cosz + 0.0096467)  
     1                                 /  
     1    (cosz**3 + 0.149864 * cosz**2 + 0.0102963 * cosz + .000303978)
        airmassf(z,patm1) = min(patm1*airmass_cosz(cosd(min(z,93.)))  
     1                        ,40.*(1.0+sqrt(max(1.0-patm1,0.)))) 

!       Twilight GHI (http://www.lotek.com/blue_twilight.pdf)
        difftwi(altf) = ! W/m**2                      
     1  4. * exp(0.50 * altf - 0.108 * altf**2 - .0044 * altf**3)

        yinterp(x1,x2,y1,y2,x) = y1 + ((x-x1)/(x2-x1)) * (y2-y1) 

        parameter (efficiency_lum_sun = .136)
        wm2sr_to_nl_550(x) = 2.113e8 * x  
        wm2sr_to_nl_sun(x) = wm2sr_to_nl_550(x) * efficiency_lum_sun

        parameter (rpd = 3.14159 / 180.)

!       Determine shadow regions of the input 3-D cloud array

        real ext_g(nc),twi_trans_c(nc) ! od per airmass, tramsmissivity

        real clwc_3d(ni,nj,nk) ! kg/m**3
        real cice_3d(ni,nj,nk) ! kg/m**3
        real rain_3d(ni,nj,nk) ! kg/m**3
        real snow_3d(ni,nj,nk) ! kg/m**3
        real cond_3d(ni,nj,nk) ! kg/m**3 (effective LWC)
!       real aod_3d(ni,nj,nk)  ! aerosol extinction coefficient
        real bi_coeff(2,2),tri_coeff(2,2,2)
        real heights_3d(ni,nj,nk)    ! MSL
        real transm_3d(ni,nj,nk)     ! O
        real transm_4d(ni,nj,nk,nc)  ! O
        real transm_4d_m(nc)
        real transm_2t(ni,nj)        ! terrain 2D transmission
        real transm_3d_vint(ni,nj,nk)! L
        real uprad_4d_m(nc)
        real sum_aod_rad_opac(nc)
        real slant2_odc(nc)
        real slant2_trans_odc(nc)
        real sum_god(nc)
        real heights_1d(nk)
        real topo_a(ni,nj)
        real grdasp_ll(ni,nj)
        real lat(ni,nj)
        real lon(ni,nj)
        real projrot_2d(ni,nj)
        real topo_albedo_2d(nc,ni,nj)
        real swi_2d(ni,nj)           ! I (use during twilight)
        real pres_3d(ni,nj,nk)
        real pres_1d(nk)
        real view_alt(minalt:maxalt,minazi:maxazi)
        real view_az(minalt:maxalt,minazi:maxazi)

        real sol_alt(ni,nj)
        real sol_azi(ni,nj)
        real alt_norm(ni,nj)             ! Solar Alt w.r.t. terrain normal
        real eobsc(ni,nj)                ! array of 'eobsl' values

        real moon_alt(ni,nj)
        real moon_azi(ni,nj)
        real moon_mag,moon_mag_thr,mean_free_path

        real obj_alt(ni,nj)
        real obj_azi(ni,nj)
        real obj_bri_a(ni,nj)

!       https://en.wikipedia.org/wiki/Irradiance#Irradiance
        real gnd_glow(ni,nj)        ! ground lighting intensity (nL)    
        real gnd_radc(nc,ni,nj)     ! ground spectral radiance (wm2srnm)
        real uprad_3d(nc,ni,nj)     ! layer spectral irradiance                 
        real sfc_glow(ni,nj)        ! pass into get_cloud_rad
        real city_colrat(nc)        /1.5,1.0,0.5/
        real ghi_2d(ni,nj)          ! derived from cloud rad 
        real dhi_2d(ni,nj)          ! diffuse horizontal irradiance
        real dhic_clr(nc)
        real dhic_term_sun(nc)
        real dhic_term_shade(nc)
        real dhi_2d_cld(nc)
        real bni_clr(nc)
        real bhi_clr(nc)
        real absorption(nc)

!       Note that 'ghi_2d' and 'dhi_2d' aren't yet returned for wider use.
!       However they are used to calculate 'topo_gti' and 'gtic' that are 
!       passed back. 'topo_gti' is more legacy and 'gtic' is becoming used.

!       Spectral irradiance W/m2/nm OR normalized to solar spectrum?              
!       For now this is normalized to the solar spectral irradiance
!       http://www.soda-is.com/eng/education/plane_orientations.html
        real ghic_2d(nc,ni,nj)      ! (global horizontal)
        real dhic_2d(nc,ni,nj)      ! (diffuse horizontal)
        real bhic_2d(nc,ni,nj)      ! (direct/beam horizontal) 
        real bnic_2d(nc,ni,nj)      ! (direct/beam normal) 
        real frac_dir_a(ni,nj)      ! debugging
        real transm_tn_a(ni,nj)     ! debugging

        logical l_solar_eclipse, l_radtran /.false./, l_spherical
        logical l_atten_bhd /.true./, l_box, l_latlon_grid, l_binary
        logical l_terrain_following, l_fullres_wdw
        integer idebug_a(minalt:maxalt,minazi:maxazi)

        parameter (nsp = 4)
        real cvr_path_sum_sp(nsp)
        real cvr_path_sum_sp_w(nsp)
        real wt_sp(nsp) ; data wt_sp /1.0,0.5,0.02,0.0714/

        real elong_p(minalt:maxalt,minazi:maxazi)     ! dummy
        real elong(minalt:maxalt,minazi:maxazi)       ! potential future
        real aod_ray_eff(minalt:maxalt,minazi:maxazi) ! zenithal
        real aod_ray_dir(minalt:maxalt,minazi:maxazi) ! zenithal direct (dummy)
        real r_cloud_3d(minalt:maxalt,minazi:maxazi)  ! cloud opacity
        real cloud_od(minalt:maxalt,minazi:maxazi)    ! cloud optical depth
        real cloud_od_sp(minalt:maxalt,minazi:maxazi,nsp)! cloud species tau
        real cloud_od_sp_w(minalt:maxalt,minazi:maxazi,nsp)! cloud species tau (front weighted)
        real r_cloud_rad(minalt:maxalt,minazi:maxazi)    ! sun to cloud transmissivity (direct+fwd scat)
        real cloud_rad_c(nc,minalt:maxalt,minazi:maxazi) ! sun to cloud transmissivity (direct+fwd scat) * solar color/int
        real cloud_rad_w(minalt:maxalt,minazi:maxazi)    ! sun to cloud transmissivity (direct+fwd scat) * trans
        real clear_rad_c(nc,minalt:maxalt,minazi:maxazi) ! clear sky illumination (twilight)
        real clear_radf_c(nc,minalt:maxalt,minazi:maxazi)! integrated 
               ! fraction of gas illuminated by the sun along line of sight
               ! (consider Earth's shadow + clouds, used when sun is below
               !  the horizon), attenuated from consideration behind clouds
        real clear_rad_c_nt(nc,minalt:maxalt,minazi:maxazi)! sky spectral radiance from nlights
        real ag_2d(minalt:maxalt,minazi:maxazi)       ! dummy
        real airmass_2_cloud_3d(minalt:maxalt,minazi:maxazi)
        real airmass_2_topo_3d(minalt:maxalt,minazi:maxazi) ! relative to zenith at std atmos
        real topo_gti(minalt:maxalt,minazi:maxazi) ! global terrain normal irradiance
        real topo_albedo(nc,minalt:maxalt,minazi:maxazi)
        real topo_ri(minalt:maxalt,minazi:maxazi)
        real topo_rj(minalt:maxalt,minazi:maxazi)
        real trace_ri(minalt:maxalt,minazi:maxazi)
        real trace_rj(minalt:maxalt,minazi:maxazi)
        real gtic(nc,minalt:maxalt,minazi:maxazi)  ! spectral terrain GNI
        real dtic(nc,minalt:maxalt,minazi:maxazi)  ! sp terrain diffuse
        real btic(nc,minalt:maxalt,minazi:maxazi)  ! sp beam terrain normal
        real emic(nc,minalt:maxalt,minazi:maxazi)  ! spectral exitance (solar rel)   
        real aod_2_cloud(minalt:maxalt,minazi:maxazi) ! slant path
        real aod_2_topo(minalt:maxalt,minazi:maxazi)  ! slant path
        real*8 dist_2_topo(minalt:maxalt,minazi:maxazi) ! slant dist
        real aod_ill(minalt:maxalt,minazi:maxazi)     ! slant path
        real aod_ill_dir(minalt:maxalt,minazi:maxazi) ! slant path, atten behind clouds
        real aod_tot(minalt:maxalt,minazi:maxazi)     ! slant path
        real sum_odrad_c(nc)
        real sum_odrad_c_last(nc)
        real minalt_deg,maxalt_deg,minazi_deg,maxazi_deg

!       Interpolation info
        integer ialt_intl(minalt:maxalt,minazi:maxazi)
        integer ialt_inth(minalt:maxalt,minazi:maxazi)
        integer jazi_intl(minalt:maxalt,minazi:maxazi)
        integer jazi_inth(minalt:maxalt,minazi:maxazi)
        logical l_pix_trace(minalt:maxalt,minazi:maxazi)

        real*8 dslant1_h,dslant2,dz1_l,dz1_h,ht_l,ht_h

        character*1 cslant
        character var*3,comment*125,ext*31,units*10

        integer icall_rad /0/
        save icall_rad               

        integer icall_uprad /0/
        save icall_uprad               

        real solalt_last /0.0/
        save solalt_last

        do ii = 1,ni
        do jj = 1,nj
           projrot_2d(ii,jj) =
     1     projrot_latlon(lat(ii,jj),lon(ii,jj),istatus)
        enddo ! jj
        enddo ! ii

        crep_thr = 0. ! 0.25
        icd = 1

        kstart = 0 ! 0 means sfc, otherwise level of start

        I4_elapsed = ishow_timer()

        write(6,3)i,j,htagl,ri_obs,rj_obs     
3       format(' Subroutine get_cloud_rays... ',2i5,f8.2,2f9.3)
        
!       Grid north azimuth - true north value
        projrot = projrot_latlon(rlat,rlon,istatus)
        write(6,*)'lat/lon/projrot',rlat,rlon,projrot

!       moon_alt = -10.0
!       moon_azi = 0.

        radius_earth_8_thirds = 6371.e3 * 2.6666666

        pstd = 101325.

        if(lat(1,1) .eq. lat(2,1))then
            l_latlon_grid = .true.
            rimid = float(ni+1)/2.
            write(6,*)' latlon grid with rimid = ',rimid
        else
            l_latlon_grid = .false.
            rimid = float(ni+1)/2.
        endif

!       Initialize
        airmass_2_cloud_3d = 0.
        airmass_2_topo_3d = 0.
!       dist_2_topo = 0.
        topo_gti = 0.
        topo_albedo = 0.0 
        r_cloud_rad = 1.
        cloud_od = 0.
        cloud_od_sp = 0.
        aod_2_cloud = 0.
        aod_2_topo = 0.
        aod_ill = 0.
        aod_ill_dir = 0.
        aod_ill_opac = 0.
        aod_ill_opac_potl = 0.
        aod_tot = 0.
        cloud_rad_w = 1.
        clear_rad_c = 0.
        clear_radf_c = 0.
        topo_ri = 0.
        topo_rj = 0.
        trace_ri = 0.
        trace_rj = 0.
        transm_2t = r_missing_data
        idebug_a = 0

        icloud_tot = 0                          

!       This may need to be conducted at each domain grid point?
!       Special measures may be needed when the observer is in 
!       the night yet can see (from high altitude) terrain in sunlight.
!       We can add the 'horz_dep' to the value of 'sol_alt' for example,
!       or try to use 'solalt_limb_true' calculated below.
        if(((sol_alt(i,j) .lt. -6.0  .AND. moon_alt(i,j) .gt.  0.0) .OR.
     1      (sol_alt(i,j) .lt. -16.0 .AND. moon_alt(i,j) .gt. -6.0))
     1                         .AND. 
     1                moon_mag .le. moon_mag_thr                  )then
            obj_alt = moon_alt
            obj_azi = moon_azi + projrot_2d(i,j)
            obj_bri = 10. ** ((-26.7 - moon_mag)*0.4)
            moon_cond = 1
            write(6,*)' Object is moon, brightness is:',obj_bri
        else
            obj_alt = sol_alt
            obj_azi = sol_azi + projrot_2d(i,j)
            obj_bri = 1.0
            moon_cond = 0
            write(6,*)' Object is sun, brightness is:',obj_bri
        endif

        obj_bri_a(:,:) = obj_bri * (1. - eobsc(:,:))

!       Possibly cloud_rad_c = r_missing_data would help interp
        if(sol_alt(i,j) .gt. 0.)then
            cloud_rad_c = obj_bri
        else
            cloud_rad_c = 0.      
        endif

        write(6,*)' range of clwc,cice is ',maxval(clwc_3d)
     1                                     ,maxval(cice_3d)

        write(6,*)' range of heights is ',minval(heights_3d)
     1                                   ,maxval(heights_3d)

        I4_elapsed = ishow_timer()

        twi_alt = -4.5
      
        if(sol_alt(i,j) .lt. twi_alt)then
          write(6,*)' Call get_sfc_glow'
          call get_sfc_glow(ni,nj,grid_spacing_m,lat,lon
     1                     ,sfc_glow,gnd_glow)
          write(6,*)' Cloud glow for observer (pop nL) is '
     1                     ,sfc_glow(i,j)
          write(6,*)' Grnd glow at observer (pop nL) location is '
     1                     ,gnd_glow(i,j)

!         Update just the gnd_glow with 'get_nlights'
          write(6,*)' Call get_nlights'
          call get_nlights(ni,nj,grid_spacing_m,lat,lon
     1                    ,gnd_glow)

          write(6,*)' range of gnd_glow (nl) is ',minval(gnd_glow)
     1                                           ,maxval(gnd_glow)

          if(minval(gnd_glow) .lt. 0.)then
             write(6,*)' ERROR: Invalid range of gnd_glow'
             write(6,*)' Remove/Regenerate nlights_multispectral.dat?'
             istatus = 0
             return
          endif
          
          write(6,29)i,j,gnd_glow(i,j)
29        format('  Grnd glow (wm2sr) at observer location is ',2i5
     1                                                       ,f10.6)       

          igmin = max( 1,i-120)
          igmax = min(ni,i+120)
          do ig = igmin,igmax
            write(6,30)ig,j,gnd_glow(ig,j),lat(ig,j),lon(ig,j)
30          format(' Gnd glow at',2i5,f11.5,2f9.3)         
          enddo ! ig

          obs_glow_gnd = gnd_glow(i,j)

!         Get upward radiation from ground lights
          if(icall_uprad .eq. 0)then
            I4_elapsed = ishow_timer()
            do ic = 1,nc
!             wm2sr to wm2srnm (= wm2 to wm2nm)
!             convert radiance to spectral radiance
!             (same as solar irradance to solar spectral irradiance)
              call get_fluxsun(wa(ic),1,1,fasun)
              rad_to_sprad = fasun / ghi_zen_toa
              write(6,31)ic,wa(ic),fasun,rad_to_sprad
31            format('ic/wa(um)/fasun/rad_to_sprad',i3,2f9.5,f11.7)
              gnd_radc(ic,:,:) = gnd_glow(:,:) * city_colrat(ic) 
     1                         * rad_to_sprad
            enddo ! ic
            do il = 0,1
              ilevel = il*(nk-1)+1
              write(6,*)' Looping level for call to get_uprad_lyr'
     1                  ,ilevel
              ht = (20000. * il) + 1000.       ! height of aerosol layer
!             use 'aef' in subroutine?
              call get_uprad_lyr(ni,nj,gnd_radc,ht
     1                          ,uprad_4d(:,:,ilevel,:)
     1                          ,upxrad_3d(:,:,ilevel)
     1                          ,upyrad_3d(:,:,ilevel))
              I4_elapsed = ishow_timer()
            enddo ! i
            icall_uprad = 1

            write(6,*)' Vertically interpolate uprad layers'
            do k = 2,nk-1
              frack = float(k-1) / float(nk-1)
              uprad_4d(:,:,k,:) = uprad_4d(:,:,1,:)  * (1.-frack) 
     1                          + uprad_4d(:,:,nk,:) * frack 
              upxrad_3d(:,:,k)  = upxrad_3d(:,:,1)  * (1.-frack) 
     1                          + upxrad_3d(:,:,nk) * frack 
              upyrad_3d(:,:,k)  = upyrad_3d(:,:,1)  * (1.-frack) 
     1                          + upyrad_3d(:,:,nk) * frack 
            enddo ! k
          endif ! icall_uprad

          do k = 1,nk
              write(6,32)k,uprad_4d(i,j,k,2),upxrad_3d(i,j,k)
     1                                      ,upyrad_3d(i,j,k)
32            format(1x,'uprad for observer (wm2nm)',i4,e14.3,2f9.3)
          enddo ! k

!         Temporary conversion from wm2sr to nL
          do ii = 1,ni
          do jj = 1,nj
            gnd_glow(ii,jj) = wm2sr_to_nl_sun(gnd_glow(ii,jj))
          enddo ! jj
          enddo ! ii
          
          write(6,*)' Grnd glow (nL) at observer location is '
     1             ,gnd_glow(i,j)

        else
          write(6,*)' Skip call to get_sfc_glow - solalt is'
     1                                                  ,sol_alt(i,j)      
          uprad_4d(:,:,:,:) = 0.
          obs_glow_gnd = 0.
        endif

        write(6,*)' obs_glow_gnd - (wm2sr) at observer location is'
     1             ,obs_glow_gnd

        I4_elapsed = ishow_timer()

!       Note early twilight will not have 'get_cloud_rad' using a special
!       twilight source as 'get_cloud_rad_faces' is used more explicitly
        grid_size_deg = float(ni) * grid_spacing_m / 110000.
        write(6,*)' grid_size_deg = ',grid_size_deg

!       This is being tested with 3.0 and 4.0 altitude threshold
        if(obj_alt(i,j) .ge. 3.0  .AND. ! as uncorrected for refraction
     1     htagl .lt. 300e3             ! not near geosynchronous
     1                                           )then 
          if(newloc .eq. 1)then
            write(6,*)' call get_cloud_rad (newloc = 1)'
            call get_cloud_rad(obj_alt,obj_azi,sol_alt(i,j),sol_azi(i,j)
     1                    ,clwc_3d,cice_3d
     1                    ,rain_3d,snow_3d,topo_a,lat,lon
     1                    ,heights_3d,transm_3d,transm_4d,i,j,ni,nj,nk
!    1                    ,l_solar_eclipse
     1                    ,twi_alt)  ! I

!           Write out transm_3d field
            var = 'TRN'
            ext = 'trn'
            units = 'none'
            comment = 
     1         '3D downward shortwave from cloud (relative to TOA)'     
            call put_laps_3d(i4time,ext,var,units,comment
     1                      ,transm_4d(:,:,:,1),ni,nj,nk)

          else
            write(6,*)' Skip call to get_cloud_rad (newloc = 0)'
          endif

        else ! called if obj_alt < low threshold
            if( (sol_alt(i,j) - twi_alt) * 
     1          (solalt_last  - twi_alt) .le. 0.)then
                write(6,*)' Crossed twi_alt compared with prior call'
                icall_rad = 0 
            endif

            if(icall_rad .eq. 0)then
              if(htagl .le. 300e3)then
                write(6,*)' call get_cloud_rad_faces...'
                call get_cloud_rad_faces(              
     1            obj_alt,obj_azi,                   ! I
     1            sol_alt(i,j),sol_azi(i,j),         ! I 
     1            clwc_3d,cice_3d,rain_3d,snow_3d,   ! I
     1            topo_a,                            ! I
     1            ni,nj,nk,i,j,                      ! I
     1            heights_3d,                        ! I 
     1            transm_2t,transm_3d,transm_4d)     ! O
              else ! consider experimental version for global domain
                do ii = 1,ni
                do jj = 1,nj
                  if(l_latlon_grid)then
                    grdasp_ll(ii,jj) = min(1.0 / cosd(lat(ii,jj)),20.)
                  else
                    grdasp_ll(ii,jj) = 1.0
                  endif
                enddo ! jj
                enddo ! ii
                write(6,*)' call get_cloud_rad_faces2...'
                htstart = topo_sfc+htagl  ! MSL
                horz_dep_d = horz_depf(max(htstart,0.),earth_radius)
                call get_cloud_rad_faces2(              
     1            obj_alt,obj_azi,horz_dep_d,        ! I
     1            sol_alt(i,j),sol_azi(i,j),         ! I 
     1            clwc_3d,cice_3d,rain_3d,snow_3d,   ! I
     1            topo_a,grdasp_ll,                  ! I
     1            ni,nj,nk,i,j,                      ! I
!    1            ni,nj,nk,ni/2,nj/2,                ! I
     1            heights_3d,                        ! I 
     1            sfc_glow,                          ! I
     1            transm_3d,transm_4d)               ! O
              endif
              icall_rad = 1
            else
              write(6,*)' skip call to get_cloud_rad_faces'
            endif

            I4_elapsed = ishow_timer()

!           Write out transm_3d field
            var = 'TRN'
            ext = 'trn'
            units = 'none'
            comment = 
     1         '3D downward shortwave from cloud (relative to TOA)'     
            call put_laps_3d(i4time,ext,var,units,comment,transm_3d
     1                      ,ni,nj,nk)

        endif

!       Use this version of 'transm_3d' for vertical interpolation
        transm_3d_vint(:,:,:) = transm_3d(:,:,:)
        do ii = 1,ni
        do jj = 1,nj
          do k = 2,nk
            if(heights_3d(ii,jj,k) .gt. topo_a(ii,jj))then
              transm_3d_vint(ii,jj,k-1) = transm_3d(ii,jj,k)
              goto 41
            endif
          enddo ! k
41        continue
        enddo ! jj
        enddo ! ii

        I4_elapsed = ishow_timer()

!       Obtain surface rad from 3D cloud rad fields, etc.
        ghi_zen_sfc = ghi_zen_toa * zen_kt ! from cloud_rad module
        arg2 = asind(10./ghi_zen_sfc)      ! approximately 0.5 degrees
        write(6,*)'ghi_zen_sfc/arg2=',ghi_zen_sfc,arg2

        do ii = 1,ni
        do jj = 1,nj

          bhic_2d(:,ii,jj) = 0. ! initialize
          dhic_2d(:,ii,jj) = 0. ! initialize
          ghic_2d(:,ii,jj) = 0. ! initialize

          ghi_2d(ii,jj) = 1e-10 ! initialize to moderately small value

!         if(ii .eq. (ii/10)*10 .and. jj .eq. nj/2)then
          if(ii .eq. i .and. jj .eq. j)then
             idb_solar = 1
             write(6,*) 
             write(6,*)' obj_alt  observer =',ii,jj,obj_alt(ii,jj)
          else
             idb_solar = 0
          endif

          if(obj_alt(ii,jj) .gt. 0.)then ! sun or moon above horizon
            ecl = 1. - eobsc(ii,jj)
            ecld = max(ecl,.0015)

            do kk = 1,nk
              kkp1 = min(kk+1,nk)

!             Consider first point above terrain instead so we can capture
!             terrain shadowing?
!             if(transm_3d(ii,jj,kk) .gt. 0. .and. 
!    1           transm_3d(ii,jj,kk) .ne. r_missing_data)then
              if(heights_3d(ii,jj,kk) .gt. topo_a(ii,jj))then

                if(sol_alt(i,j) .ge. twi_alt .and.
     1             transm_3d(ii,jj,kkp1) .gt. 0. .and.
     1             transm_3d(ii,jj,kkp1) .ne. r_missing_data)then
                    absorption(:) = 1. - transm_4d(ii,jj,kkp1,:)
     1                                 / transm_3d(ii,jj,kkp1)
                else ! absorbing aerosols at night not yet working
                    absorption(:) = 0.
                endif

                if(transm_2t(ii,jj) .ne. r_missing_data)then
                  transm_tn = transm_2t(ii,jj)
                else
                  transm_tn = transm_3d(ii,jj,kk)
                endif
               
                cloud_albedo = 1. - transm_tn

!               Here we might need only the aerosol absorption and not the
!               gas scattering. Both of these feed into the present calculation
!               of absorption.

                transm_tn_pr1 = transm_tn * (1.-absorption(2))
     1                   / (1. - topo_albedo_2d(2,i,j) * cloud_albedo)
                transm_tn_pr = transm_tn_pr1 ! tn_pr1 or tn (testing)

!               Correct for diffuse radiation at low obj altitude
!               'solalt_eff' is an empirical function to force the effective
!               solar altitude to be 1.5 degrees when the actual is 0.
!               'sb_corr' is an empirical sky brightness correction in 
!               astronomical magnitudes.
!               These empirical forumlae should give a reasonable value for 
!               DHI even when the sun is very near the horizon.
!               'frac_dir' is an estimate of the direct solar radiation now
!               including scattered light in a 5 degree radius of the sun.

!               Note the 1.5 value should be set to 'arg2' (about 0.5) to 
!               match the sun below horizon case. 'transm' should 
!               asymptotically be forced to a value of 1 when the sun is 
!               on the horizon.

!               We might want to calculate the direct independently, then
!               set GHI = Diffuse Horizontal + Direct Horizontal

                objalt_eff = obj_alt(ii,jj) 
     1                     + (1.5 * cosd(obj_alt(ii,jj))**100.)
                ghi_clr = obj_bri * sind(objalt_eff) * ghi_zen_sfc *ecld
                ghi_2d(ii,jj) = transm_tn_pr * ghi_clr    

!               Direct / Diffuse
!               This takes into account clouds though might also consider terrain
                frac_dir1 = max(transm_tn**4.,.0039) ! 5 deg circ of sctrd
                frac_dir = transm_tn + (1.-transm_tn) * frac_dir1

                patm_sfc = ztopsa(topo_a(ii,jj)) / 1013.25
                if(.false.)then
                   sb_corr = 4.0 * (1.0 - (sind(obj_alt(ii,jj))**0.5))
                else
                   sarg = resin(sind(obj_alt(ii,jj)))
                   sb_corr = 3.38 * (1.0 - sarg**0.33)
                endif
                dhi_grn_clr_frac = ext_g(2) * patm_sfc 
     1                           * 10.**(-0.4*sb_corr)
                dhi_2d_clear = ghi_zen_toa * ecld * dhi_grn_clr_frac 
     1                       * obj_bri

                dhi_2d(ii,jj) = ghi_2d(ii,jj) * (1. - frac_dir) 
     1                     * transm_tn_pr
     1                     + dhi_2d_clear * frac_dir
!               dhi_2d(ii,jj) = max(dhi_2d_clear,dhi_2d_est) 
!               dhi_2d(ii,jj) = min(dhi_2d(ii,jj),ghi_2d(ii,jj))

!               Spectral normalized values (under construction)
                airmass_g = airmassf(90.-obj_alt(ii,jj),patm_sfc)
                do ic = 1,nc
                  bni_clr(ic) = obj_bri * ecl 
!    1                        * trans(airmass_g * ext_g(ic))
                enddo ! ic
                bnic_2d(:,ii,jj) = bni_clr(:) * frac_dir
     1                                        * (1.-absorption(:))

                bhi_clr(:) = bni_clr(:) * sind(obj_alt(ii,jj))
     1                                  * (1.-absorption(:))

                bhic_2d(:,ii,jj) = 
     1              bnic_2d(:,ii,jj) * sind(obj_alt(ii,jj))

                colexp = 0.5 ! color based on clearness of sky
                dhic_clr(:) = dhi_grn_clr_frac * (ext_g(:)/.09)**colexp
!               dhic_2d(:,ii,jj) = dhi_2d(ii,jj) / ghi_zen_toa
                dhic_term_sun
     1             = obj_bri * ecld * dhic_clr(:) * frac_dir   ! clr
     1             + bhi_clr(:) * transm_tn_pr * (1.-frac_dir) ! cld
                dhic_term_shade = obj_bri * ecld * dhic_clr(:)
                dhic_2d(:,ii,jj) =    transm_tn  * dhic_term_sun
     1                           +(1.-transm_tn) * dhic_term_shade

!               ghic_clr = ghi_clr / ghi_zen_toa
                ghic_2d(:,ii,jj) = bhic_2d(:,ii,jj) + dhic_2d(:,ii,jj) 

                if(idb_solar .eq. 1)then
                   write(6,*)' 1st point above terrain under observer'
                   write(6,*)' bhi_clr   observer =',bhi_clr  
                   write(6,*)' ghi_clr   observer =',ghi_clr  
!                  write(6,*)' ghic_clr  observer =',ghic_clr  
                   write(6,*)' dhic_clr  observer =',dhic_clr  
                   write(6,*)' frac_dir  observer =',frac_dir
                   write(6,*)' transm_tn observer =',transm_tn
                   write(6,*)' trnalb    observer ='
     1                                   ,topo_albedo_2d(2,i,j)
                   write(6,*)' cldalb    observer =',cloud_albedo
                   write(6,*)' transm_p  observer =',transm_tn_pr1
                   write(6,*)' ecl       observer =',ecl
                   write(6,42)bhic_2d(:,ii,jj)
     1                       ,bhic_2d(2,ii,jj)*ghi_zen_toa
42                 format(' bhic_2d   observer =',f9.5
     1                                           ,' W/m^2 ~=',f10.5)
                   write(6,43)dhic_2d(:,ii,jj)
     1                       ,dhic_2d(2,ii,jj)*ghi_zen_toa
43                 format(' dhic_2d   observer =',f9.5
     1                                           ,' W/m^2 ~=',f10.5)
                   write(6,44)ghic_2d(:,ii,jj)
     1                       ,ghic_2d(2,ii,jj)*ghi_zen_toa
44                 format(' ghic_2d   observer =',f9.5
     1                                           ,' W/m^2 ~=',f10.5)

                endif

                frac_dir_a(ii,jj) = frac_dir
                transm_tn_a(ii,jj) = transm_tn
!               if(jj .eq. 516)then
!                 write(6,*)'ii/jj/kk/transm_tn',ii,jj,kk
!    1                      ,transm_3d(ii,jj,kk:kk+1),transm_tn_a(ii,jj)
!               endif
               
                goto 4 ! we have risen above the terrain
              endif
            enddo ! kk
4           continue

!           Gradual transition from daytime to nighttime equations
            ramp_day = min(obj_alt(ii,jj)/3.,1.) ! ramp from 0-3 deg    
            if(ramp_day .lt. 1.0)then
              diffuse_twi = difftwi(sol_alt(ii,jj)) * obj_bri ! W/m**2
              dhi_2d(ii,jj) = ramp_day * dhi_2d(ii,jj) 
     1                      + (1.-ramp_day) * diffuse_twi
              ghi_2d(ii,jj) = ramp_day * ghi_2d(ii,jj) 
     1                      + (1.-ramp_day) * diffuse_twi
            endif

          endif ! object is above horizon

          if(sol_alt(ii,jj) .lt. 0. .and. sol_alt(ii,jj) .gt. -12.)then

!           compare with compare_analysis_to_rad.f (cloud analysis code)
!           Empirical forumlae for DHI and GHI when the sun is below the horizon.
!           Note they are equal since the direct is zero.
!           We may want to factor in cloud reduction during twilight if
!           it can be obtained via 'swi_2d' or another cloud field.

!           Add solar twilight to DHI and GHI               
            diffuse_twi = difftwi(sol_alt(ii,jj)) ! W/m**2
            if(.true.)then
              dhi_2d(ii,jj) = dhi_2d(ii,jj) + diffuse_twi
              ghi_2d(ii,jj) = ghi_2d(ii,jj) + diffuse_twi
            else ! future possibly
              dhi_2d(ii,jj) = dhi_2d(ii,jj) + swi_2d(ii,jj)
              ghi_2d(ii,jj) = ghi_2d(ii,jj) + swi_2d(ii,jj)
            endif

            if(ii .eq. ni/2 .and. jj .eq. nj/2)then
              write(6,*)' diffuse_twi CTR = '
     1                 ,diffuse_twi,diffuse_twi/ghi_zen_toa
              write(6,*)' swi_2d      CTR = ',swi_2d(ii,jj)
            endif

            if(idb_solar .eq. 1)then
              write(6,*)' dhic_2d  observer 1 =',dhic_2d(:,ii,jj)
              write(6,*)' ghic_2d  observer 1 =',ghic_2d(:,ii,jj)
            endif
            
!           Add solar twilight to spectral normalized values
            dhic_2d(:,ii,jj) = dhic_2d(:,ii,jj) 
     1                       + diffuse_twi / ghi_zen_toa
            ghic_2d(:,ii,jj) = ghic_2d(:,ii,jj) 
     1                       + diffuse_twi / ghi_zen_toa

            if(idb_solar .eq. 1)then
              write(6,*)' dhic_2d  observer 2 =',dhic_2d(:,ii,jj)
              write(6,*)' ghic_2d  observer 2 =',ghic_2d(:,ii,jj)
            endif

          endif

          if(ghi_2d(ii,jj) .lt. 1e-10 .AND. idb_solar .eq. 1)then
            write(6,*)' WARNING, small ghi_2d ',ii,jj,ghi_2d(ii,jj)
     1                                         ,sol_alt(ii,jj)
!           stop
          endif
         
        enddo ! jj
        enddo ! ii
        write(6,*)' range of ghi_2d = ',minval(ghi_2d),maxval(ghi_2d)
        write(6,*)' range of dhi_2d = ',minval(dhi_2d),maxval(dhi_2d)
        write(6,*)' range of dhic_2d(1) = '
     1             ,minval(dhic_2d(1,:,:)),maxval(dhic_2d(1,:,:))
        write(6,*)' range of dhic_2d(2) = '
     1             ,minval(dhic_2d(2,:,:)),maxval(dhic_2d(2,:,:))
        write(6,*)' range of dhic_2d(3) = '
     1             ,minval(dhic_2d(3,:,:)),maxval(dhic_2d(3,:,:))

!       do jj = 1,nj
!       do ii = 1,ni
!           if(gnd_glow(ii,jj) .gt. 0.)then
!               write(6,*)' ii,jj,gndglow',ii,jj,gnd_glow(ii,jj)
!           endif
!       enddo ! i
!       enddo ! j

        I4_elapsed = ishow_timer()

        if(sol_alt(i,j) .lt. twi_alt)then ! Modify moon glow section in green channel
                                      ! Preserve sfc sky glow section in red channel
            do k = 1,nk
              transm_4d(:,:,k,2) = transm_4d(:,:,k,2) * obj_bri_a(:,:) ! correct for sun/moon brightness
            enddo ! k
            write(6,*)
     1    ' transm_4d red is sky brightness, green is moon brightness'
            write(6,*)' Range of transm_4d(green channel) = '
     1            ,minval(transm_4d(:,:,:,2)),maxval(transm_4d(:,:,:,2))
        else
            do ic = 1,nc
             do k = 1,nk
              transm_4d(:,:,k,ic) = transm_4d(:,:,k,ic) * obj_bri_a(:,:) ! correct for sun/moon brightness
             enddo ! k
            enddo ! ic
            write(6,*)' heights_3d column =     ',heights_3d(i,j,:)
            write(6,*)' transm_3d column =      ',transm_3d(i,j,:)
            write(6,*)' transm_3d_vint column = ',transm_3d_vint(i,j,:)
            write(6,*)' transm_4d B column =    ',transm_4d(i,j,:,2)
        endif
        write(6,*)' Range of transm_4d(red channel) = '
     1           ,minval(transm_4d(:,:,:,1)),maxval(transm_4d(:,:,:,1))

        write(6,*)' clwc column = ',clwc_3d(i,j,:)
        write(6,*)' cice column = ',cice_3d(i,j,:)
        write(6,*)' rain column = ',rain_3d(i,j,:)
        write(6,*)' snow column = ',snow_3d(i,j,:)
        write(6,*)' aod_3d column = ',aod_3d(i,j,:)

        I4_elapsed = ishow_timer()

        topo_max_ang = 45.0
        topo_max_ht = maxval(topo_a)

!       MEE values using constants from 'module_cloud_rad.f90'
!       Values are 1.5 / (rho * reff)
!       clwc2alpha = 75. 
        clwc2alpha = 1.5 / (rholiq  * reff_clwc)
        cice2alpha = 1.5 / (rholiq  * reff_cice)
        rain2alpha = 1.5 / (rholiq  * reff_rain)
        snow2alpha = 1.5 / (rhosnow * reff_snow)

        write(6,*)' clwc2alpha (MEE) is ',clwc2alpha

        cond_3d = clwc_3d 
     1          + cice_3d * (cice2alpha/clwc2alpha)
     1          + rain_3d * (rain2alpha/clwc2alpha) 
     1          + snow_3d * (snow2alpha/clwc2alpha)

!       MEE relative to clwc
        wt_sp(1) = clwc2alpha/clwc2alpha
        wt_sp(2) = cice2alpha/clwc2alpha
        wt_sp(3) = rain2alpha/clwc2alpha
        wt_sp(4) = snow2alpha/clwc2alpha

        write(6,*)'wt_sp(1) = ',wt_sp(1)
        write(6,*)'wt_sp(2) = ',wt_sp(2)
        write(6,*)'wt_sp(3) = ',wt_sp(3)
        write(6,*)'wt_sp(4) = ',wt_sp(4)

        ri = ri_obs ! i
        rj = rj_obs ! j

!       kstart = 0 ! 0 means sfc, otherwise level of start

!       Set up starting information for ray trace from observer

        if(.true.)then ! general info
          write(6,*)' i/j/topo_sfc = ',i,j,topo_sfc
          write(6,*)' height column = ',heights_3d(i,j,:)
          do k = 1,nk-1
            if(heights_3d(i,j,k)   .le. topo_sfc .AND.
     1         heights_3d(i,j,k+1) .ge. topo_sfc      )then
                frach = (topo_sfc - heights_3d(i,j,k)) / 
     1                  (heights_3d(i,j,k+1) - heights_3d(i,j,k))
                ksfc = k
                rksfc = float(ksfc) + frach
            endif
          enddo ! k
        endif

        if(htagl .eq. 0.)then ! start from sfc
          htstart = topo_sfc  ! MSL
          write(6,*)' Start at sfc, ksfc = ',ksfc
          rkstart = rksfc
          patm = pres_3d(i,j,ksfc) / 101325.
        else ! start aloft
          htstart = topo_sfc+htagl  ! MSL
          if(htstart .le. heights_3d(i,j,nk))then ! in domain
            do k = 1,nk-1
              if(heights_3d(i,j,k)   .le. htstart .AND.
     1           heights_3d(i,j,k+1) .ge. htstart      )then
                frach = (htstart - heights_3d(i,j,k)) / 
     1                  (heights_3d(i,j,k+1) - heights_3d(i,j,k))
                fracl = 1.0 - frach
                kstart = k
                rkstart = float(kstart) + frach
              endif
            enddo ! k
            if(kstart .eq. 0)then
               write(6,*)' ERROR: kstart = 0',htstart,heights_3d(i,j,:)
               istatus = 0
               return
            endif
            patm = (  pres_3d(i,j,kstart)   * fracl 
     1              + pres_3d(i,j,kstart+1) * frach ) / 101325.
          else ! set up extended k values above the domain
            write(6,*)' starting from above the domain'
            rkstart = float(nk) + (htstart - heights_3d(i,j,nk))/1000. 
            kstart = int(rkstart)
            patm = (ztopsa(htstart)*100.) / 101325.
            horz_dep_tod_d = horz_depf(htstart - heights_3d(i,j,nk)
     1                                ,earth_radius)
            write(6,*)' horizon tod depression / tod ang diam = '
     1                 ,horz_dep_tod_d, 180.-2.*horz_dep_tod_d
          endif
          write(6,*)' Start aloft at k/rk/p = ',kstart,rkstart,patm
        endif

        radius_earth_eff = earth_radius * (1. + (0.33333 * patm))
        radius_start_eff = radius_earth_eff + htstart
        radius_start     = earth_radius     + htstart

!       horz_dep_d = sqrt(2.0 * max(htstart,0.) / earth_radius) * 180./3.14
        horz_dep_d = horz_depf(max(htstart,0.),earth_radius)
        write(6,*)' horizon depression of limb = ',horz_dep_d
        write(6,*)' radius_start = ',radius_start
        if(horz_dep_d .gt. 10.)then
          l_spherical = .true.
        else
          l_spherical = .false.
        endif

!       Always set to true to eliminate small polar artifact at nadir
!       if(grid_spacing_m .ge. 9000. .or. .true.)then
        if(grid_spacing_m .ge. 9000.)then
          l_box = .true.
        else
          l_box = l_spherical
        endif

        if(sol_alt(i,j) .le. 0. .and. sol_alt(i,j) .ge. twi_alt)then
          l_box = .true. ! testing
        endif

        write(6,*)' solalt observer grid / l_box = ',sol_alt(i,j),l_box

!       Used to call skyglow_phys vs skyglow_phys_twi for example
        twi_0 = min(-15.0,-horz_dep_d) ! output
        write(6,*)' twi_0 = ',twi_0

        aod_vrt = aod * exp(-(htstart-redp_lvl)/aero_scaleht)

        write(6,*)' rkstart/htstart/patm = ',rkstart,htstart,patm
        write(6,*)' aero_scaleht = ',aero_scaleht
        write(6,*)' aod/redp_lvl/aod_vrt = ',aod,redp_lvl,aod_vrt
                  
!       Calculate observer mean free paths
        k_obs = nint(rkstart)
        if(k_obs .le. nk)then
          k_vis = k_obs
        else
          k_vis = ksfc
        endif

        if(.true.)then ! general info
          aero_ext_coeff = aod_3d(i,j,k_vis) ! m**-1            

          vsb_const = -log(.02)
          if(aero_ext_coeff .gt. 0.)then
            mean_free_path = 1.0/aero_ext_coeff
            vsb_mi = mean_free_path / 1609.34 * vsb_const
            write(6,5)1.0/aero_ext_coeff,vsb_mi,lat(i,j),lon(i,j),k_vis
 5          format(
     1      ' observer aerosol mean free path / vsb / lat / lon / k '
     1          ,e16.3,'m',2x,f8.2,'mi',5x,2f8.2,i4)
          endif

!         Will need MEE information for hydrometeors
          if(clwc_3d(i,j,k_vis) .gt. 0.)then
            mean_free_path = 1.0/clwc_3d(i,j,k_vis)
            vsb_mi = mean_free_path / 1609.34 * vsb_const
            write(6,*)' observer clwc mean free path / vsb = '
     1                               ,mean_free_path,vsb_mi
          else
            write(6,*)' zero clwc visibility obstruction for observer'
          endif

          if(cice_3d(i,j,k_vis) .gt. 0.)then
            mean_free_path = 1.0/cice_3d(i,j,k_vis)
            vsb_mi = mean_free_path / 1609.34 * vsb_const
            write(6,*)' observer cice mean free path / vsb = '
     1                               ,mean_free_path,vsb_mi
          else
            write(6,*)' zero cice visibility obstruction for observer'
          endif

          if(rain_3d(i,j,k_vis) .gt. 0.)then
            mean_free_path = 1.0/rain_3d(i,j,k_vis)
            vsb_mi = mean_free_path / 1609.34 * vsb_const
            write(6,*)' observer rain mean free path / vsb = '
     1                               ,mean_free_path,vsb_mi
          else
            write(6,*)' zero rain visibility obstruction for observer'
          endif

          if(snow_3d(i,j,k_vis) .gt. 0.)then
            mean_free_path = 1.0/snow_3d(i,j,k_vis)
            vsb_mi = mean_free_path / 1609.34 * vsb_const
            write(6,*)' observer snow mean free path / vsb = '
     1                               ,mean_free_path,vsb_mi
          else
            write(6,*)' zero snow visibility obstruction for observer'
          endif
        endif

        aod_ray_eff = aod_vrt
        aod_ray_dir = aod_vrt

        heights_1d(:) = heights_3d(i,j,:)
        pres_1d(:)    = pres_3d(i,j,:)

        idelt = nint(2. / alt_scale)
        minalt_deg = float(minalt) * alt_scale
        maxalt_deg = float(maxalt) * alt_scale
        minazi_deg = float(minazi) * azi_scale
        maxazi_deg = float(maxazi) * azi_scale

        write(6,*)' range of altitudes is ',minalt_deg,maxalt_deg
        write(6,*)' range of azimuths is  ',minazi_deg,maxazi_deg

!       Optional setup of hi-res window
        if(grid_spacing_m .le. 30.)then
           write(6,*)' Setting up full res window'
           l_fullres_wdw = .true.
           azimin_full = 207.
           azimax_full = 222.
           altmin_full = 0.
           altmax_full = 5.
           jazimin_full = nint(azimin_full / azi_scale)
           jazimax_full = nint(azimax_full / azi_scale)
           ialtmin_full = nint(altmin_full / alt_scale)
           ialtmax_full = nint(altmax_full / alt_scale)
           jazi_delt_outer = nint(1.0 / azi_scale)
        else
           l_fullres_wdw = .false.
           jazi_delt_outer = 1 ! disabled
        endif

        l_pix_trace(:,:) = .false.

!       azid1 = 46. ; azid2 = 226.
        azid1 = 90. ; azid2 = 270.
        if(sol_alt(i,j) .gt. 0.)then
            azid1 = sol_azi(i,j)
            azid2 = mod(azid1+180.,360.)
!           azid2 = azid1
        elseif(moon_alt(i,j) .gt. 0.)then
            azid1 = moon_azi(i,j)
            azid2 = mod(azid1+180.,360.)
        endif
        if(htstart .gt. 100e3)then
            azid1 = int(sol_azi(i,j))
            azid2 = azid1
        endif
        azid1 = 212. ; azid2 = 212.

        write(6,*)'azid1/azid2 = ',azid1,azid2

!       Loop through altitudes as seen from the observer
        do ialt = minalt,maxalt

!        l_process = .false.

         call get_val(ialt,minalt,alt_scale,altray)

!        Fill in pseudo topo outside horizontal domain
!        if(dist_2_topo(ialt,jazi) .eq. 0. .and. 
!    1      view_altitude_deg .lt. -horz_dep_d)then
!           dist_2_topo(ialt,jazi) = htstart / 
!    1                               sind(-view_altitude_deg)
!        endif

         if(altray .lt. -horz_dep_d)then
           call get_topo_info(altray,htstart,earth_radius,0
     1                       ,alt_norm_dum,dist_to_topo)
         else
           dist_to_topo = 0.
         endif

         dist_2_topo(ialt,:) = dist_to_topo ! initialize
         grdasp = (alt_scale / azi_scale) 

         if(grid_spacing_m .lt. 500.)then
           topo_alt_max = 13. ! highest terrain elevation angle
           topo_alt_thresh = max(topo_alt_max,10.)
         else
           topo_alt_thresh = 10.
         endif

         if(altray .ge. 20.)then     ! 2 deg alt, 2 deg azi
             if(ialt .eq. (ialt/idelt)*idelt)then
                 jazi_delt = nint(2. / azi_scale)
             else
                 jazi_delt = maxazi - minazi
             endif
         elseif(altray .ge. topo_alt_thresh)then ! alt_scale, 1 deg azi
             jazi_delt = nint(1. / azi_scale)
         elseif(altray .le. -10. .and. htstart .lt. 100000.)then 
             jazi_delt = nint(1. / azi_scale)
         elseif(altray .le.  -5. .and. htstart .lt. 100000.)then 
             jazi_delt = max(nint(0.4/ azi_scale),1)
             if(jazi_delt .eq. 3)jazi_delt = 2
             if(jazi_delt .eq. 6)jazi_delt = 5
             if(jazi_delt .eq. 7)jazi_delt = 5
             if(jazi_delt .eq. 9)jazi_delt = 8
             if(jazi_delt .gt. 10 .and. jazi_delt .lt. 16)jazi_delt = 10
             if(jazi_delt .gt. 16 .and. jazi_delt .lt. 20)jazi_delt = 16
             if(jazi_delt .gt. 20 .and. jazi_delt .lt. 99)jazi_delt = 20
         elseif(htagl .gt. 30000e3)then  ! near geosynchronous
             grdasp = (alt_scale / azi_scale) 
     1              / cosd(min(abs(altray),89.999))
             jazi_delt_max = nint(8./azi_scale)
             if(altray .lt. 0.)then
               if    (grdasp .ge. 256.0)then
                 jazi_delt = 256
               elseif(grdasp .ge. 128.0)then
                 jazi_delt = 128                       
               elseif(grdasp .ge. 64.0)then
                 jazi_delt = 64                       
               elseif(grdasp .ge. 32.0)then
                 jazi_delt = 32                       
               elseif(grdasp .ge. 16.0)then
                 jazi_delt = 16                       
               elseif(grdasp .ge. 8.0)then
                 jazi_delt = 8                       
               elseif(grdasp .ge. 4.0)then
                 jazi_delt = 4                       
               elseif(grdasp .ge. 2.0)then
                 jazi_delt = 2                       
               else
                 jazi_delt = 1                       
               endif
             else
               jazi_delt = 1                       
             endif
             jazi_delt = min(jazi_delt,jazi_delt_max)
         else                        ! alt_scale, azi_scale
             jazi_delt = 1
         endif

         if(mod((maxazi-minazi),jazi_delt) .ne. 0)then
             write(6,*)' WARNING: jazi_delt has remainder',jazi_delt
         endif

         azi_delt_2 = float(jazi_delt) * azi_scale * 0.5

         call get_htmin(altray,patm,htstart,earth_radius,0,patm2,htmin)

!        This responds to the presence of inner/outer windows
         jazi_delt_eff = max(jazi_delt,jazi_delt_outer)

         if(altray .eq. nint(altray))then
           write(6,*)'altray/htmin/dist_to_topo = '
     1               ,altray,htmin,dist_to_topo
           write(6,52)ialt,altray,jazi_delt,jazi_delt_eff,grdasp
52         format('alt/jazi_delt/jazi_delt_eff,grdasp',i5,f9.3,2i5,f9.3)
         endif

         altray_limb = altray + horz_dep_d
         radius_limb = 90. - horz_dep_d

!        Determine pixels to trace for this altitude
         if(jazi_delt .lt. maxazi-minazi)then ! alt ring is traced
            l_pix_trace(ialt,minazi:maxazi:jazi_delt_eff) = .true.
         endif
         if((l_fullres_wdw .eqv. .true.) .and. ialt .ge. ialtmin_full
     1                   .and. ialt .le. ialtmax_full)then ! fill inner window
            l_pix_trace(ialt,jazimin_full:jazimax_full) = .true.
         endif

         do jazi = minazi,maxazi

          if(l_pix_trace(ialt,jazi))then

           view_azi_deg = float(jazi) * azi_scale
           azigrid = modulo(view_azi_deg + projrot,360.)

           if((abs(view_azi_deg - azid1) .lt. azi_delt_2 .or. 
     1         abs(view_azi_deg - azid2) .lt. azi_delt_2      ) .AND.
     1        (abs(altray) .eq. 12  .or. abs(altray) .eq. 9 .or.
     1         (altray .ge. -5. .and. altray .le. 9.) .or. 
     1         ialt .eq. minalt .or. abs(altray) .eq. 14. .or.
     1         abs(altray) .eq. 16. .or.
     1         altray .eq. -7.5 .or.
     1         abs(altray) .eq. 21.00 .or.
     1         (abs(altray_limb) / radius_limb .le. 0.04) .or.
     1         (altray_limb/radius_limb .lt. .01
     1               .and. ialt .eq. (ialt/10) * 10)      .or.
     1         abs(altray) .eq. 20. .or. abs(altray) .eq. 30. .or.
     1         abs(altray) .eq. 45. .or. abs(altray) .eq. 63.5 .or.
     1         (altray .ge. -75.00 .and. altray .le. -70.00) .or.
     1         abs(altray) .eq. 75.) 
!    1               .AND. altray .eq. nint(altray) 
     1                                                   )then
             idebug = 1
             idebug_a(ialt,jazi) = 1
          else
             idebug = 0
          endif

!         High custom
!         if(ialt .eq. minalt+(maxalt-minalt)*(90-5)/180)then ! -5 degrees alt
          if(.false.)then
             if(view_azi_deg*2. .eq. nint(view_azi_deg*2.) .and.
     1          view_azi_deg .le. 90.                        )then
                idebug = 1
                idebug_a(ialt,jazi) = 1
             else
                idebug = 0
                idebug_a(ialt,jazi) = 0
             endif
          endif

!         Non-verbose (low observer)                    
          if(l_box .eqv. .false.)then
            if(htstart .lt. 25. .or. altray .gt. 0.0)then
              if(mode_aero_cld .lt. 3)then
                idebug = 0 ! ; idebug_a(ialt,jazi) = 0
              endif
            endif
!         else ! extra verbose
!           idebug = 1
          endif

!         Zenith / Nadir
          if(jazi .eq. minazi .and. abs(altray) .eq. 90.)then
              idebug = 1
              idebug_a(ialt,jazi) = 1
          endif

!         Horizon
          if(jazi .eq. minazi .and. altray .eq. 0.)then
              idebug = 1
              idebug_a(ialt,jazi) = 1
          endif

!         Custom
          if(ialt .eq. -108 .and. jazi .eq. 340)then
              idebug = 1
              idebug_a(ialt,jazi) = 1
          endif

!         Extra verbose
          if(idebug .eq. 1)then
            write(6,*)
            write(6,*)'alt/azi = ',ialt,jazi,altray,view_azi_deg
          endif

!         Trace towards sky from each grid point
!         view_altitude_deg = max(altray,0.) ! handle Earth curvature later
          view_altitude_deg = altray

!         Get direction cosines based on azimuth
          xcosg = sind(azigrid)
          ycosg = cosd(azigrid)
          ycost = cosd(view_azi_deg)

!         Initialize variables for this ray
          icloud = 0
          cvr_path_sum = 0.
          cvr_path_sum_sp = 0.
          cvr_path_sum_sp_w = 0.
          ltype_1st = 0
          sum_odrad = 0.
          sum_odrad_c = 0.
          sum_odrad_w = 0.
          sum_clrrad = 0.     ! used for clear_radf_c
          sum_clrrad_pot = 0. ! used for clear_radf_c
          sum_aod = 0.
          sum_aod_ill = 0.
          sum_aod_ill_dir = 0.
          sum_aod_ill_opac = 0.
          sum_aod_ill_opac_potl = 0.
          sum_aod_rad_opac(:) = 0.
          sum_am2cld_num = 0.
          sum_am2cld_den = 0.
          sum_am2cld_atten = 0.
          sum_god(:) = 0.
          sum_xcosup = 0.
          sum_ycosup = 0.
          frac_fntcloud = 1.0
          ray_topo_diff_h = htagl ! 0.
          ray_topo_diff_m = 0.
          ichecksum = 0
          box_max = 0.
          dist_box = 1e30

          if(.true.)then
              r_cloud_3d(ialt,jazi) = 0.
              cloud_od(ialt,jazi) = 0.
              cloud_od_sp(ialt,jazi,:) = 0.
              
              ri1 = ri
              rj1 = rj

              grid_factor = grid_spacing_m / 3000.

              if(view_altitude_deg .lt. -15.)then
                  rkdelt1 = -1.00 * grid_factor
                  rkdelt2 = -1.00 * grid_factor
!             elseif(view_altitude_deg .lt. -4.5)then ! produces red image at -6
!                 rkdelt1 = -0.50 * grid_factor
!                 rkdelt2 = -0.50 * grid_factor
              elseif(view_altitude_deg .lt. 0.0)then
                  rkdelt1 =  0.0
                  rkdelt2 =  0.0
              elseif(view_altitude_deg .le. 0.5)then
                  rkdelt1 = 0.0 ! 0.01 * grid_factor 
                  rkdelt2 = 0.0 ! 0.10 * grid_factor 
              elseif(view_altitude_deg .le. 4.)then
                  rkdelt1 = 0.0 ! 0.10 * grid_factor
                  rkdelt2 = 0.0 ! 0.25 * grid_factor
              elseif(view_altitude_deg .le. 15.)then
                  rkdelt1 = 0.0 ! 0.50 * grid_factor
                  rkdelt2 = 0.0 ! 0.50 * grid_factor
              elseif(view_altitude_deg .le. 45.)then
!                 if(htstart .lt. 15000.)then ! normal
                      rkdelt1 = 1.00 * grid_factor
                      rkdelt2 = 1.00 * grid_factor
!                 else ! high altitude observer
!                     rkdelt1 = 0.00
!                     rkdelt2 = 0.00
!                 endif
              else
                  rkdelt1 = 0.25
                  rkdelt2 = 0.25
              endif

              if(l_box .eqv. .true.)then
                  rkdelt1 =  0.0
                  rkdelt2 =  0.0
              endif

              if(rkstart .gt. float(nk) .and. 
     1           view_altitude_deg .lt. 0.0)then
                  iabove = 1 ! vantage point above domain
                  rkdelt1 =  0.0
                  rkdelt2 =  0.0
                  sind_view = sind(-view_altitude_deg)
                  slant2_optimal = min(250./sind_view,grid_spacing_m)
              else
                  iabove = 0
                  sind_view = sind(abs(view_altitude_deg))
                  slant2_optimal = min(250./sind_view,grid_spacing_m)
!                 slant2_optimal = grid_spacing_m
              endif

!             arg = max(view_altitude_deg,1.0)
!             rkdelt = tand(90. - arg)                     
!             rkdelt = max(min(rkdelt,2.0),0.5)
              
              if(idebug .eq. 1)then
                if(htagl .gt. 8000e3)I4_elapsed = ishow_timer()
                write(6,11)altray,view_azi_deg,azigrid,ialt,jazi
     1                    ,jazi_delt,rkdelt,i,j,slant2_optimal,l_box
11              format(
     1    'Trace the slant path (alt/azi-g/ialt/jazi/jdelt/rkdelt/i/j):'
     1                ,f8.3,2f6.1,i6,2i5,f8.4,2i5,f7.0,l2)                                     
                write(6,12)                                   
12              format('      dz1_l        dz1_h      dxy1_l    dxy1_h',
     1           '  rinew  rjnew   rk    ht_m   topo_m   path     ',
     1           'lwc    ice    rain   snow      slant  cvrpathsum',
     1           '  cloudfrac am2cld  sumclrd smclrdp  am1_h  cld_rd',
     1           ' cld_rd_w  aeroext transm3 aod_sm aod_sm_ill/potl')
              endif

!             Initialize ray
              dxy1_h = 0.
              dz1_h = 0.
              dslant1_h = 0.
              airmass1_h = 0.
              ihit_topo = 0       ! hit topo
              ihit_bounds = 0     ! crossed beyond top of domain
              ioutside_domain = 0 ! outside horizontal domain
              rk_h = rkstart
              rinew_h = i ! for l_spherical = T
              rjnew_h = j ! for l_spherical = T

              rkdelt = rkdelt1

              rk = rkstart          
              ls = 0
              iwrite = 0
              trace_minalt = 1e10
              if(grid_spacing_m .le. 30.)then
                 cvr_path_thr = 2000. 
              else
                 cvr_path_thr = 1.0
              endif

              do while((rk .le. float(nk)-rkdelt .or. iabove .eq. 1)
     1           .AND. ihit_topo .eq. 0
     1           .AND. ihit_bounds .eq. 0
     1           .AND. ioutside_domain .eq. 0
     1           .AND. htmin .lt. 100000.     ! will hit atmosphere
     1           .AND. cvr_path_sum .le. cvr_path_thr) ! Tau < ~75

                ls = ls + 1

                in_h = min(max(nint(rinew_h),1),ni)
                jn_h = min(max(nint(rjnew_h),1),nj)

                if(rkdelt .ne. 0.)then ! trace by pressure levels

                  if(.false.)then ! optimize step size
!                 if(view_altitude_deg .gt. 0.)then ! optimize step size
                    rkdelt_vert = 1.0 - (rk - int(rk))
                    kk_ref = min(int(rk),nk-1)
                    if(.not. l_terrain_following)then
                      delta_grid_height = heights_1d(kk_ref+1)
     1                                  - heights_1d(kk_ref)
                    else
                      delta_grid_height = heights_3d(in_h,jn_h,kk_ref+1)
     1                                  - heights_3d(in_h,jn_h,kk_ref)
                    endif
                    aspect_ratio = delta_grid_height / grid_spacing_m
                    rkdelt_horz = aspect_ratio / tand(view_altitude_deg) 
                    rkdelt = max(min(rkdelt_horz,rkdelt_vert),0.2)
                    if(idebug .eq. 1)then
                      write(6,*)' rk,rkdelt_vert,rkdelt_horz,rkdelt',
     1                            rk,rkdelt_vert,rkdelt_horz,rkdelt
                      write(6,*)' view_altitude_deg,delta_grid_height',
     1                            view_altitude_deg,delta_grid_height
                      write(6,*)' aspect_ratio',
     1                            aspect_ratio
!                     stop
                    endif
                  endif

                  rk = rk + rkdelt

                  rk_l = rk - rkdelt
                  kk_l = int(rk_l)
                  frac_l = rk_l - float(kk_l)

!                 if(kk_l .ge. nk .or. kk_l .le. 0)then
!                   write(6,*)' ERROR: rk/kk_l',rk,kk_l
!                 endif

                  rk_h = rk
                  kk_h = int(rk_h)
                  frac_h = rk_h - float(kk_h)

                  if(rk_h .gt. float(nk) .or. rk_h .lt. 1.)then
                    write(6,*)' ERROR: rk_h/kk_h',rk_h,kk_h,ialt,jazi
                    write(6,*)' ERROR: kk_h is out of bounds '
     1                       ,ls,iabove
     1                       ,rkstart,rk_h,nk,rkdelt,ht_h,ht_m,topo_m
     1                       ,ihit_topo,ihit_bounds,rinew_l,rjnew_l
     1                       ,rjnew_h,jnew_m
     1                       ,dy1_l,dy1_h,rj,ycosg
                    istatus = 0
                    return
                  endif

                  if(.not. l_terrain_following)then
                    ht_l = heights_1d(kk_l) * (1. - frac_l) 
     1                   + heights_1d(kk_l+1) * frac_l

                    pr_l = pres_1d(kk_l) * (1. - frac_l) 
     1                   + pres_1d(kk_l+1) * frac_l

                    if(rk_h .lt. float(nk) .AND. rk_h .ge. 1.0)then
                      ht_h = heights_1d(kk_h) * (1. - frac_h) 
     1                     + heights_1d(kk_h+1) * frac_h

                      pr_h = pres_1d(kk_h) * (1. - frac_h) 
     1                     + pres_1d(kk_h+1) * frac_h

                    elseif(rk_h .eq. float(nk))then
                      ht_h = heights_1d(kk_h)

                      pr_h = pres_1d(kk_h)
                    else
                      write(6,*)' ERROR: rk_h is out of bounds '
     1                       ,ls,iabove
     1                       ,rkstart,rk_h,nk,rkdelt,ht_h,ht_m,topo_m
     1                       ,ihit_topo,rinew_l,rjnew_l,rjnew_h,jnew_m
     1                       ,dy1_l,dy1_h,rj,ycosg
                      istatus = 0
                      return
                    endif

                  else ! l_terrain_following
                    ht_l = heights_3d(in_h,jn_h,kk_l) * (1. - frac_l) 
     1                   + heights_3d(in_h,jn_h,kk_l+1) * frac_l

                    pr_l = pres_3d(in_h,jn_h,kk_l) * (1. - frac_l) 
     1                   + pres_3d(in_h,jn_h,kk_l+1) * frac_l

                    if(rk_h .lt. float(nk) .AND. rk_h .ge. 1.0)then
                      ht_h = heights_3d(in_h,jn_h,kk_h) * (1. - frac_h) 
     1                     + heights_3d(in_h,jn_h,kk_h+1) * frac_h

                      pr_h = pres_3d(in_h,jn_h,kk_h) * (1. - frac_h) 
     1                     + pres_3d(in_h,jn_h,kk_h+1) * frac_h

                    elseif(rk_h .eq. float(nk))then
                      ht_h = heights_3d(in_h,jn_h,kk_h)

                      pr_h = pres_3d(in_h,jn_h,kk_h)
                    else
                      write(6,*)' ERROR: rk_h is out of bounds '
     1                       ,ls,iabove
     1                       ,rkstart,rk_h,nk,rkdelt,ht_h,ht_m,topo_m
     1                       ,ihit_topo,rjnew_l,rjnew_h,jnew_m
     1                       ,dy1_l,dy1_h,rj,ycosg
                      istatus = 0
                      return
                    endif
                  endif ! l_terrain_following

                  dz1_l = ht_l - htstart        
                  dz1_h = ht_h - htstart          
                  dz2   = dz1_h - dz1_l ! layer

                  if(.false.)then
                    dxy1_l = dz1_l * tand(90. - view_altitude_deg)
                    dxy1_h = dz1_h * tand(90. - view_altitude_deg)
                    dxy2   = dz2   * tand(90. - view_altitude_deg)
                  else ! horzdist call (determine dist from height & elev)
                    aterm = 1. / radius_earth_8_thirds

                    if(view_altitude_deg .gt. 89.99)then
                       tanterm = 89.99
                    elseif(view_altitude_deg .lt. -89.99)then
                       tanterm = -89.99
                    else
                       tanterm = view_altitude_deg
                    endif
                    bterm = tand(tanterm)

                    cterm = dz1_l                                            
                    discriminant = 4.*aterm*cterm + bterm**2.

                    if(discriminant .gt. 0.)then     
                      if(view_altitude_deg .gt. 0.)then ! above horizon
                          dxy1_l = ( sqrt(discriminant) - bterm)   
     1                                              / (2.*aterm)
                      else                           ! below horizon
                          dxy1_l = (-sqrt(discriminant) - bterm)   
     1                                              / (2.*aterm)
                      endif
                    else                             ! catch start point
                      dxy1_l = 0.
                    endif

                    cterm = dz1_h                                          
                    discriminant = 4.*aterm*cterm + bterm**2.

                    if(view_altitude_deg .gt. 0.)then ! above horizon
                        dxy1_h = ( sqrt(discriminant) - bterm) 
     1                                           / (2.*aterm)
                    else                           ! below horizon
                        dxy1_h = (-sqrt(discriminant) - bterm) 
     1                                           / (2.*aterm)
                    endif
                 
                    dxy2 = dxy1_h - dxy1_l
                  endif

                  slant2 = sqrt(dxy2**2 + dz2**2) ! layer
                  slant1_l = dslant1_h
                  dslant1_h = dslant1_h + slant2

                  airmassv = (pr_l - pr_h) / pstd
                  airmass2 = airmassv * (slant2 / dz2)

                  if(dxy2 .gt. grid_spacing_m)then
                    cslant = '*'
                  else
                    cslant = ' '
                  endif

                else  ! Trace by slant range (rkdelt = 0. )
                  dz1_l = dz1_h
                  ht_l = htstart + dz1_l

                  if(iabove .eq. 1)then
                    if(ht_l .lt. 10000.)then
                      htarg = 200.
                    else
                      htarg = 500.
                    endif
                    slant2_optimal = min(htarg/sind_view,grid_spacing_m)
                  endif
                  dslant2=max(dble(slant2_optimal),ht_l-30000.D0)
                  if((l_box .eqv. .true.) .and. dist_box .lt. 1e30)then
                    dslant2 = dble(dist_box)
                  endif
                  slant1_l = dslant1_h
                  dslant1_h = dslant1_h + dslant2
                  slant2 = dslant2

                  dxy1_l = dxy1_h
                  dxy1_h = dslant1_h * cosd(altray)

!                 Determine height from distance and elev
                  if(l_spherical .eqv. .false.)then ! more approximate
                    dz1_h = dslant1_h * sind(altray) 
!    1                    + dxy1_h**2 / (2. * radius_earth_eff)
     1                    + curvat(dxy1_h,radius_start_eff)
!                   gc_deg = (dxy1_h/radius_start)/rpd
                  else ! need more accuracy
                    slant1_h = dslant1_h
                    dz1_h = dcurvat2(dslant1_h,dble(radius_start_eff)
     1                                        ,dble(altray))   
                    gc_deg = asind(cosd(altray) 
     1                     * dslant1_h / (radius_start+dz1_h)) 
                  endif

                  rk_l = rk_h
                  ht_h = htstart + dz1_h

                  rk_h = r_missing_data

!                 Pseudo grid values above/below LAPS domain
                  if(.not. l_terrain_following)then
                    if(ht_h .gt. heights_1d(nk) .and. iabove .eq. 1)then
                      rk_h = float(nk) + (ht_h - heights_1d(nk)) / 1000.
                      pr_h = pres_1d(nk)
     1                   * exp(-(ht_h - heights_1d(nk)) / 8000.)
                    elseif(ht_h .lt. heights_1d(1))then
                      rk_h = 1.        + (ht_h - heights_1d(1))  / 1000.
                      pr_h = pres_1d(1)
     1                   * exp(-(ht_h - heights_1d(1))  / 8000.)
                    else ! inside vertical domain
                      do k = 1,nk-1
                        if(heights_1d(k)   .le. ht_h .AND.
     1                     heights_1d(k+1) .ge. ht_h      )then
                            delta_h = heights_1d(k+1) - heights_1d(k)
                            frach = (ht_h - heights_1d(k)) / delta_h
                            rk_h = float(k) + frach
                            pr_h = (pres_1d(k)*(1.-frach)) 
     1                           + (pres_1d(k+1)*frach)
                        endif
                      enddo ! k
                    endif

                  else ! l_terrain_following
                    if(ht_h .gt. heights_3d(in_h,jn_h,nk) .and. 
     1                                         iabove .eq. 1)then
                      rk_h = float(nk) + 
     1                       (ht_h - heights_3d(in_h,jn_h,nk)) / 1000.
                      pr_h = pres_1d(nk)
     1                * exp(-(ht_h - heights_3d(in_h,jn_h,nk)) / 8000.)
                    elseif(ht_h .lt. heights_3d(in_h,jn_h,1))then
                      rk_h = 1.+ (ht_h - heights_3d(in_h,jn_h,1))/1000.       
                      pr_h = pres_3d(in_h,jn_h,1)
     1                   * exp(-(ht_h - heights_3d(in_h,jn_h,1)) /8000.)
                    else ! inside vertical domain
                      do k = 1,nk-1
                        if(heights_3d(in_h,jn_h,k)   .le. ht_h .AND.
     1                     heights_3d(in_h,jn_h,k+1) .ge. ht_h     )then
                            delta_h = heights_3d(in_h,jn_h,k+1) 
     1                              - heights_3d(in_h,jn_h,k)
                            frach = (ht_h - heights_3d(in_h,jn_h,k)) 
     1                            / delta_h
                            rk_h = float(k) + frach
                            pr_h = (pres_3d(in_h,jn_h,k)*(1.-frach)) 
     1                           + (pres_3d(in_h,jn_h,k+1)*frach)
                        endif
                      enddo ! k
                    endif

                  endif ! l_terrain_following

                  rk = rk_h

                  airmass2 = (slant2 / 8000.) * (pr_h / pstd)

                endif ! rkdelt .ne. 0.

!               Extra verbose
                if(altray .le. -60. .and. 
     1              (idebug .eq. 1))then
                  if(ls .le. 15 .or. ls .eq. 100 .or. ls .eq. 1000
     1                          .or. ls .eq. 10000)then
                    write(6,61)ls,rk,ht_h,dslant1_h,slant1_h
61                  format('  ls/rk/ht/dsl1/sl1 =',i6,f9.3,f10.0,2f11.0)
                  elseif(altray .eq. -63.5)then
                    write(6,61)ls,rk,ht_h,dslant1_h,slant1_h
!                   idebug = 1
                  endif
                endif

                rk_m = 0.5 * (rk_l + rk_h)
                ht_m = 0.5 * (ht_l + ht_h)
                k_m = nint(rk_m)

                if((rk_h-rk_l) .gt. 0.)then
                    grid_dh = (ht_h-ht_l) / (rk_h-rk_l)
                    grid_tan = grid_dh / grid_spacing_m
                else
                    grid_tan = 0.25 ! nominal value
                endif

                cslant = ' '

!               ioutside_domain = 0

                if(rk_m    .lt. float(nk) .AND.
     1             rk_m    .gt. 1.0            )then ! in vertical domain

                 if(l_spherical .eqv. .true.)then ! more exact
                  rinew_l = rinew_h
                  rjnew_l = rjnew_h

!                 Obtain latlon from spherical geometry
                  TLat = ASinD(ycost*SinD(gc_deg)*CosD(rLat) 
     1                 + SinD(rLat)*CosD(gc_deg))

                  CosDLon = (CosD(gc_deg) - SinD(rLat)*SinD(TLat)) 
     1                    / (CosD(rLat)*CosD(TLat))
                  If(Abs(CosDLon).gt.1.)CosDLon=Sign(1.,CosDLon)
                  DLon = ACosD(CosDLon)
                  
                  If(view_azi_deg.ge..0.and.view_azi_deg.le.180.)Then ! east
!                 If(azigrid.ge..0.and.azigrid.le.180.)Then ! grid east
                    TLon=rLon+DLon
                  Else ! west
                    TLon=rLon-DLon
                  EndIf

                  if(.not. l_latlon_grid)then ! speed test
                    call latlon_to_rlapsgrid(tlat,tlon,lat,lon,ni,nj
     1                                      ,rinew_h,rjnew_h,istatus)
!                   if(idebug .eq. 1)then
!                      write(6,62)tlat,tlon,gc_deg,rinew_h,rjnew_h,dlon
!62                    format('   check latlon',6f9.3)                    
!                   endif

                    rinew_m = 0.5 * (rinew_l + rinew_h)
                    rjnew_m = 0.5 * (rjnew_l + rjnew_h)

                  else ! assume a lat/lon grid
                    if(tlon .gt. 180.)then
                        tlon1 = tlon - 360.
                    elseif(tlon .lt. -180.)then
                        tlon1 = tlon + 360.
                    else
                        tlon1 = tlon
                    endif
                    rinew_h = yinterp(lon(1,1),lon(ni,1)
     1                               ,1.,float(ni),tlon1)
                    rinew_h = min(max(rinew_h,1.),float(ni))
                    rjnew_h = yinterp(lat(1,1),lat(1,nj)
     1                               ,1.,float(nj),tlat)
                    rjnew_h = min(max(rjnew_h,1.),float(nj))

                    if(abs(rinew_h - rinew_l) .lt. rimid)then
                      rinew_m = 0.5 * (rinew_l + rinew_h)
                      rjnew_m = 0.5 * (rjnew_l + rjnew_h)
                    else
                      rinew_m = 0.5 * (rinew_l + rinew_h)
                      rinew_m = modulo((rinew_m+rimid-1.),float(ni))+0.5    
                      rjnew_m = 0.5 * (rjnew_l + rjnew_h)
                    endif
                  endif

                 else ! more approximate
                  dx1_l = dxy1_l * xcosg
                  dy1_l = dxy1_l * ycosg

                  dx1_h = dxy1_h * xcosg
                  dy1_h = dxy1_h * ycosg

                  ridelt_l = dx1_l / grid_spacing_m
                  rjdelt_l = dy1_l / grid_spacing_m

                  ridelt_h = dx1_h / grid_spacing_m
                  rjdelt_h = dy1_h / grid_spacing_m

                  rinew_l = ri + ridelt_l
                  rinew_h = ri + ridelt_h

                  rjnew_l = rj + rjdelt_l
                  rjnew_h = rj + rjdelt_h

                  rinew_m = 0.5 * (rinew_l + rinew_h)
                  rjnew_m = 0.5 * (rjnew_l + rjnew_h)

                 endif ! l_spherical

!                Can be used to target box boundaries
                 if(l_box .eqv. .true.)then
                   dids = (rinew_h - rinew_l) / slant2
                   djds = (rjnew_h - rjnew_l) / slant2
                   dhds = (ht_h - ht_l) / slant2
                 endif

                 if(ht_m .lt. trace_minalt)then
                   trace_ri(ialt,jazi) = rinew_m
                   trace_rj(ialt,jazi) = rjnew_m
                   trace_minalt = ht_m
                 endif

                 rni = ni
                 rnj = nj

                 airmass1_l = airmass1_h
                 airmass1_h = airmass1_h + airmass2

                 inew_m = nint(rinew_m); jnew_m = nint(rjnew_m)
                 inew_mb = max(min(inew_m,ni),1)

                 if(rinew_h .ge. 1. .and. rinew_h .le. rni .AND. 
     1              rjnew_h .ge. 1. .and. rjnew_h .le. rnj)then ! horz domain

                  i1 = max(min(int(rinew_m),ni-1),1)
                  j1 = max(min(int(rjnew_m),nj-1),1)
                  k1 = max(min(int(rk_m)   ,nk-1),1) 

                  if(i1 .le. 0)then
                   write(6,*)' software error ',i1,rinew_l,rinew_m
     1                                         ,rinew_h,rk_m
                   istatus = 0
                   return
                  endif

!                 Along-ray distance to next box face, along each axis
                  if(l_box .eqv. .true.)then
                    fi = mod(rinew_h,1.)
                    fj = mod(rjnew_h,1.)
                    fk = mod(rk_h,1.)
                    box_epsilon = 0.01

                    if(dids .ne. 0.)then
                      if(dids .lt. 0.)then
                        bi = fi                   
                      else
                        bi = 1.-fi
                      endif
                      if(bi .lt. box_epsilon)then
                        bi = bi + 1.
                      endif
                      disti = bi / abs(dids)
                    else
                      disti = 1e30
                    endif

                    if(djds .ne. 0.)then
                      if(djds .lt. 0.)then
                        bj = fj                   
                      else
                        bj = 1.-fj
                      endif
                      if(bj .lt. box_epsilon)then
                        bj = bj + 1.
                      endif
                      distj = bj / abs(djds)
                    else
                      distj = 1e30
                    endif

                    if(dhds .ne. 0.)then
                      if(dhds .lt. 0.)then          ! downward
                        bk = fk    
                        if(bk .gt. box_epsilon)then ! normal
                          knext = int(rk_h)        
                        else                        ! close to edge
                          knext = int(rk_h) - 1
                        endif
                      elseif(dhds .gt. 0.)then      ! upward
                        bk = 1.-fk
                        if(bk .gt. box_epsilon)then ! normal
                          knext = int(rk_h) + 1
                        else                        ! close to edge
                          knext = int(rk_h) + 2
                        endif
                      endif 
                      if(knext .ge. 1 .and. knext .le. nk)then
                        if(.not. l_terrain_following)then
                          disth = abs((heights_1d(knext) - ht_h) / dhds)
                        else
                          write(6,*)' ERROR in get_cloud_rays'
                          write(6,*)
     1                       ' l_box and l_terrain_following are both T'
                          istatus = 0
                          return
                        endif
                      else
                        disth = 1e30
                      endif
                    else
                      disth = 1e30
                    endif

                    dist_box = min(disti,distj,disth)
                    if(idebug .eq. 1 .and. .false.)then
                       write(6,*)'fi/dids/disti',fi,dids,disti
                       write(6,*)'fj/djds/distj',fj,djds,distj
                       write(6,66)rk_h,knext,fk,dhds,disth
     1                           ,ht_h,ht_l,slant2
66                     format(' rk/knext/fk/dhds/disth/ht/sl'
     1                       ,f9.3,i3,2f9.3,f10.1,3f12.0)
                       write(6,*)'dist_box = ',dist_box
                    endif
                  endif ! l_box

                  fi = rinew_m - i1; i2=i1+1
                  fj = rjnew_m - j1; j2=j1+1
                  fk = rk_m    - k1; k2=k1+1

                  tri_coeff(1,1,1) = (1.-fi) * (1.-fj) * (1.- fk)
                  tri_coeff(2,1,1) = fi      * (1.-fj) * (1.- fk)
                  tri_coeff(1,2,1) = (1.-fi) *     fj  * (1.- fk)
                  tri_coeff(1,1,2) = (1.-fi) * (1.-fj) *      fk
                  tri_coeff(1,2,2) = (1.-fi) *     fj  *      fk
                  tri_coeff(2,1,2) = fi      * (1.-fj) *      fk
                  tri_coeff(2,2,1) = fi      *     fj  * (1.- fk)
                  tri_coeff(2,2,2) = fi      *     fj  *      fk

!                 Update box_max only if crossed into new box
!                 ichecksum_last = ichecksum; ichecksum = i1+j1+k1
!                 if(ichecksum .ne. ichecksum_last)then
                    box_max = maxval(cond_3d(i1:i2,j1:j2,k1:k2))
!                 endif

                  if(cvr_path_sum .le. cvr_path_thr)then ! consider clouds
!                   Cloud present on any of 8 interpolation vertices
                    if(box_max .gt. 0.)then
                      cond_m = sum(tri_coeff(:,:,:) * 
     1                             cond_3d(i1:i2,j1:j2,k1:k2))

                      clwc_m = sum(tri_coeff(:,:,:) * 
     1                             clwc_3d(i1:i2,j1:j2,k1:k2))

                      cice_m = sum(tri_coeff(:,:,:) * 
     1                             cice_3d(i1:i2,j1:j2,k1:k2))

                      rain_m = sum(tri_coeff(:,:,:) * 
     1                             rain_3d(i1:i2,j1:j2,k1:k2))

                      snow_m = sum(tri_coeff(:,:,:) * 
     1                             snow_3d(i1:i2,j1:j2,k1:k2))
                    else
                      cond_m = 0.
                      clwc_m = 0.
                      cice_m = 0.
                      rain_m = 0.
                      snow_m = 0.
                    endif
                  else
                      cond_m = 0.
                      clwc_m = 0.
                      cice_m = 0.
                      rain_m = 0.
                      snow_m = 0.
                  endif

                  if(mode_aero_cld .lt. 3)then
                      aero_ext_coeff = aod_3d(inew_m,jnew_m,k_m)             
                      aod_inc = aero_ext_coeff * slant2
                      cvr_path = cond_m                                 
                  else
                      aero_ext_coeff = sum(tri_coeff(:,:,:) * 
     1                                     aod_3d(i1:i2,j1:j2,k1:k2))
                      aod_inc = aero_ext_coeff * slant2
                      cvr_path = cond_m + aero_ext_coeff / clwc2alpha
                  endif

                  cvr_path_sum_last = cvr_path_sum
                  cvr_path_sum      = cvr_path_sum + cvr_path * slant2
                  cvr_path_sum_sp(1) = 
     1            cvr_path_sum_sp(1) + clwc_m * wt_sp(1) * slant2
                  cvr_path_sum_sp(2) = 
     1            cvr_path_sum_sp(2) + cice_m * wt_sp(2) * slant2
                  cvr_path_sum_sp(3) = 
     1            cvr_path_sum_sp(3) + rain_m * wt_sp(3) * slant2
                  cvr_path_sum_sp(4) = 
     1            cvr_path_sum_sp(4) + snow_m * wt_sp(4) * slant2

!                 if(idebug .eq. 1 .OR. cond_m .gt. 0.)then
                  if(.true.)then

!                   Can this be unfairly reduced when k1 is below the terrain?
                    transm_3d_m = sum(tri_coeff(:,:,:) * 
     1                            transm_3d_vint(i1:i2,j1:j2,k1:k2))

                    if(idebug .eq. 1 .and. .false.)then
                      if(heights_3d(i1,j1,k1) .lt. topo_a(i1,j1))then
                         write(6,*)
     1                   ' Warning: suspect transm_3d_m below terrain'
     1                   ,transm_3d(i1:i2,j1:j2,k1:k2),transm_3d_m
                      endif
                    endif

                    if(.false.)then
                      write(6,*)' ERROR transm_3d_m < 0',transm_3d_m
     1                            ,sum(tri_coeff(:,:,:)),fi,fj,fk
     1                            ,' tri_coeff '
     1                            ,tri_coeff(:,:,:),i1,i2,j1,j2,k1,k2
     1                            ,' transm_3d '
     1                            ,transm_3d(i1:i2,j1:j2,k1:k2)
                      stop
                    else ! prevent negative values from extrapolation
                      transm_3d_m = max(transm_3d_m,0.)
                    endif

                    sum_odrad = sum_odrad + 
     1                 (cvr_path * slant2 * transm_3d_m)       

                    if(.true.)then
                      sum_odrad_w = sum_odrad_w + 
     1                 (cvr_path * slant2 * transm_3d_m**2)       
                    else
                      sum_odrad_w = sum_odrad_w + 
     1                 (cvr_path * slant2 * transm_3d_m 
     1                           * trans(cvr_path_sum/.013))       
                    endif

                    do ic = 1,nc
                      transm_4d_m(ic) = sum(tri_coeff(:,:,:) * 
     1                              transm_4d(i1:i2,j1:j2,k1:k2,ic))
                      if(icall_uprad .gt. 0)then
                        uprad_4d_m(ic) = sum(tri_coeff(:,:,:) * 
     1                              uprad_4d(i1:i2,j1:j2,k1:k2,ic))
                      else
                        uprad_4d_m(ic) = 0.
                      endif
                    enddo ! ic
 
                    if(icall_uprad .gt. 0)then
                      upxrad_3d_m = sum(tri_coeff(:,:,:) *
     1                              upxrad_3d(i1:i2,j1:j2,k1:k2))
                      upyrad_3d_m = sum(tri_coeff(:,:,:) *
     1                              upxrad_3d(i1:i2,j1:j2,k1:k2))

                      sum_xcosup = sum_xcosup + upxrad_3d_m
     1                           * cvr_path * slant2
                      sum_ycosup = sum_ycosup + upyrad_3d_m
     1                           * cvr_path * slant2

                    else
                      upxrad_3d_m = 0.
                      upyrad_3d_m = 0.

                    endif

                    sum_odrad_c(:) = sum_odrad_c(:) + 
     1              (cvr_path * slant2 * transm_4d_m(:))
                  endif ! idebug .eq. 1 .OR. cond_m .gt. 0.

!                 Assess topo height with respect to ray
!                 if(ht_m - htstart .le. 1000.)then ! more accurate topo when low
                  if(ht_m   .le. topo_max_ht .AND. 
     1               altray .le. topo_max_ang     )then 
!                 if(ht_m   .le. topo_max_ht)then 
                    i1 = max(min(int(rinew_m),ni-1),1)
                    fi = rinew_m - i1
                    i2 = i1+1
                    j1 = max(min(int(rjnew_m),nj-1),1)
                    fj = rjnew_m - j1
                    j2 = j1+1

                    bi_coeff(1,1) = (1.-fi) * (1.-fj)
                    bi_coeff(2,1) = fi      * (1.-fj)
                    bi_coeff(1,2) = (1.-fi) *     fj 
                    bi_coeff(2,2) = fi      *     fj 
                    topo_m = sum(bi_coeff(:,:) * topo_a(i1:i2,j1:j2))

                    i1 = min(int(rinew_h),ni-1); fi = rinew_h - i1
                    i2=i1+1
                    j1 = min(int(rjnew_h),nj-1); fj = rjnew_h - j1
                    j2=j1+1

                    bi_coeff(1,1) = (1.-fi) * (1.-fj)
                    bi_coeff(2,1) = fi      * (1.-fj)
                    bi_coeff(1,2) = (1.-fi) *     fj 
                    bi_coeff(2,2) = fi      *     fj 
                    topo_h = sum(bi_coeff(:,:) * topo_a(i1:i2,j1:j2))
                  else
                    topo_m = topo_a(inew_m,jnew_m)
                    topo_h = topo_m
                  endif

!                 Use _h values in addition to _m
!                 ray_topo_diff_m_last = ray_topo_diff_m
                  ray_topo_diff_h_last = ray_topo_diff_h
                  ray_topo_diff_m = ht_m - topo_m
                  ray_topo_diff_h = ht_h - topo_h

!                 Check both mid-point and end-point of ray segment
                  if((topo_m .gt. ht_m .or. topo_h .gt. ht_h) .AND. 
     1                                             ihit_topo .eq. 0)then
                      ihit_topo = 1

                      if(topo_m .gt. ht_m)then ! mid-point hit topo
                          ihit_topo_mid = 1
                          ihit_topo_end = 0
                          frac_step_topo = 0.5 * ray_topo_diff_h_last /
     1                     (ray_topo_diff_h_last-ray_topo_diff_m)
                      else                     ! end-point hit topo
                          ihit_topo_mid = 0
                          ihit_topo_end = 1
                          frac_step_topo = 0.5 * ray_topo_diff_m      /
     1                     (ray_topo_diff_m     -ray_topo_diff_h) + 0.5
                      endif
                  endif ! determine if we are hitting topo

                  sum_aod = sum_aod + aod_inc

!                 Test for hitting tau ~1 (nominal cloud edge)
                  if((cvr_path          .gt. 0.00) .AND.
     1               (cvr_path_sum      .le.  .013     .OR. ! tau ~1
     1                cvr_path_sum_last .eq. 0.  )            )then 
                    aod_2_cloud(ialt,jazi) = sum_aod
                    if(cvr_path_sum_sp(3) .gt. 0.)then
                      ltype_1st = 3 ! first hydrometeors contain rain                     
                    endif
!                   frac_fntcloud = 1.0
                    frac_fntcloud = trans(cvr_path_sum/.013) ! tau
!                   cvr_path_sum_sp_w(:) = cvr_path_sum_sp(:)
                  elseif(cvr_path_sum .gt. .013)then 
!                   frac_fntcloud = 0.0 ! make a more continuous function? 
                    frac_fntcloud = trans(cvr_path_sum/.013) ! tau
                    if(cvr_path_sum_last .le. .013)then ! fractional gridpt
                      frac_cvr_path = (.013-cvr_path_sum_last) / 
     1                                (cvr_path_sum-cvr_path_sum_last)
                      aod_2_cloud(ialt,jazi) = sum_aod - 
     1                                    (1.-frac_cvr_path) * aod_inc
                    endif
                  else
                    frac_fntcloud = 1.0
                  endif

                  if(l_atten_bhd)then ! .true.
!                   if(transm_3d(inew_m,jnew_m,k_m) .gt. 0.
!    1                                             .and. .false.)then
!                     solocc = transm_4d(inew_m,jnew_m,k_m,2) 
!    1                       / transm_3d(inew_m,jnew_m,k_m) 
!                   else
!                     solocc = 1.
!                   endif
                    if(sol_alt(inew_m,jnew_m) .lt. 0.0)then
                       shdw_ht = (-sol_alt(inew_m,jnew_m))**2. * 1000.
                    else
                       shdw_ht = 0.
                    endif

                    solaltarg = sol_alt(inew_m,jnew_m) + 1.
                    if(ht_h .gt. ht_l)then ! ascending ray
                       solocc_l = solocc_f(max(real(ht_l),0.)
     1                                    ,earth_radius,solaltarg)
                       solocc_m = solocc_f(max(real(ht_m),0.)
     1                                    ,earth_radius,solaltarg)
                       solocc_h = solocc_f(max(real(ht_h),0.)
     1                                    ,earth_radius,solaltarg)

!                      horz_dep_l = horz_depf(max(ht_l,0.),earth_radius)
!                      horz_dep_h = horz_depf(max(ht_h,0.),earth_radius)
 
!                      soloccarg = (ht_h - shdw_ht) / (ht_h-ht_l)
!                      solocc = min(max(soloccarg,0.),1.)
                       solocc = (solocc_l + 2. * solocc_m + solocc_h)/4.

                    elseif(ht_m .gt. shdw_ht)then 
                       solocc = 1.
                    else
                       solocc = 0.
                    endif

!                   Approximation to trans4 / trans3
                    if(shdw_ht .gt. 0.)then ! twilight
                      ht_m_eff = max(ht_m,shdw_ht)
                      amapp = 38. * exp(-min(ht_m_eff/8000.,80.))
                      trn4app = trans(amapp * 0.14)
                    elseif(.true.)then      ! daytime
                      patm_arg = exp(-min(ht_m/8000.,80.))
                      amapp = airmassf(90.-solaltarg,patm_arg)
                      trn4app = trans(amapp * 0.14)
                    else
                      trn4app = 1.
                    endif

!                   This weights where we examine illumination fraction
                    clrrad_inc_pot = airmass2
     1                             * frac_fntcloud ! was just in numerator
     1                             * solocc        ! MSL solar occultation
     1                             * trn4app       ! Gas attenuation

!                   Apply shadowing by clouds and terrain (above MSL)
                    clrrad_inc = clrrad_inc_pot 
     1                         * transm_3d(inew_m,jnew_m,k_m) 

                    sum_clrrad_pot = sum_clrrad_pot + clrrad_inc_pot
                    sum_clrrad     = sum_clrrad     + clrrad_inc

                  else
                    sum_clrrad = sum_clrrad 
     1                         + transm_3d(inew_m,jnew_m,k_m) 
     1                         * airmass2                     
                  endif  

                  slant2_od = aero_ext_coeff * slant2
                  slant2_trans_od = slant2_od * trans(sum_aod)

                  do ic = 1,nc
                    slant2_odc(ic) = aero_ext_coeff * slant2 
     1                             + airmass2 * ext_g(ic)
                    slant2_trans_odc(ic) = slant2_odc(ic) 
     1                                   * trans(sum_aod+sum_god(ic))
                  enddo ! ic

!                 Check whether near start of ray and near topo
                  if(rk_m .lt. (rkstart + 1.0) .AND. 
     1               htagl .lt. 1000.)then

!                 if(.false.)then ! near topo
                    sum_aod_ill = sum_aod_ill + aero_ext_coeff * slant2
     1                          * transm_3d(inew_m,jnew_m,int(rk_m)+1)  
                    if(transm_3d(inew_m,jnew_m,int(rk_m)+1) .gt. .1)then
                        transm_3d_dir = 
     1                exp(log(transm_3d(inew_m,jnew_m,int(rk_m)+1))*10.)
                    else
                        transm_3d_dir = 0.
                    endif

                    sum_aod_ill_dir = sum_aod_ill_dir 
     1                         + slant2_od * transm_3d_dir
     1                                     * frac_fntcloud

                    sum_aod_ill_opac = sum_aod_ill_opac 
     1                          + slant2_trans_od
     1                          * transm_3d(inew_m,jnew_m,int(rk_m)+1)  

                    sum_aod_ill_opac_potl = sum_aod_ill_opac_potl 
     1                                    + slant2_trans_od

                    if(icall_uprad .gt. 0)then
                      sum_aod_rad_opac(:) = sum_aod_rad_opac(:)
     1                         + slant2_trans_odc(:)
     1                         * uprad_4d(inew_m,jnew_m,int(rk_m)+1,:)
                    endif

                  else ! typical case far from topo
                    sum_aod_ill = sum_aod_ill + aero_ext_coeff * slant2
     1                          * transm_3d_m

!                   if(sum_aod_ill .lt. 0.)then
!                     write(6,*)' ERROR sum_aod_ill < 0',sum_aod_ill
!    1                            ,aero_ext_coeff,slant2,transm_3d_m
!                     stop
!                   endif

                    if(transm_3d_m .gt. .1)then
                        transm_3d_dir = exp(log(transm_3d_m)*10.)
                    else
                        transm_3d_dir = 0.
                    endif

                    sum_aod_ill_dir = sum_aod_ill_dir 
     1                         + slant2_od * transm_3d_dir
     1                                     * frac_fntcloud

                    sum_aod_ill_opac = sum_aod_ill_opac 
     1                         + slant2_trans_od * transm_3d_m

                    sum_aod_ill_opac_potl = sum_aod_ill_opac_potl
     1                         + slant2_trans_od

!                   Night lights and aerosols and gas irradiance sum
                    sum_aod_rad_opac(:) = sum_aod_rad_opac(:)
     1                         + slant2_trans_odc(:) * uprad_4d_m(:)

                  endif

                  if(l_radtran .eqv. .true.)then
                    dopac = 1.
                    di_a = di_a + dopac * rad * aero_ext_coeff 
     1                          / (aero_ext_coeff + alphabar_g)
                  endif

!                 Determine backscatter threshold for averaging cloud rad
!                 This can be refined by determining the phase of the
!                 cloud so far in the ray traversal
                  if(cvr_path .gt. 0.00)then  
                    frac_liq = (cvr_path_sum_sp(1) + cvr_path_sum_sp(3))
     1                       / cvr_path_sum
!                   if(ltype_1st .eq. 3)then ! 1st hydrometeors have rain
                    if(.true.)then
                        tau_thr = 2. ! 0.5
                    else                   
                        tau_thr = 15. * frac_liq + 7. * (1.-frac_liq)
                    endif
                    bks_thr = tau_thr / clwc2alpha

                    if(cvr_path_sum      .le.  bks_thr .OR. ! tau ~7-15
     1                 cvr_path_sum_last .eq. 0.            )then 
!                     Average value over cloud path (tau ~10)
                      r_cloud_rad(ialt,jazi) = sum_odrad / cvr_path_sum 
                      if(sum_odrad .gt. 0.)then
                          cloud_rad_w(ialt,jazi) = sum_odrad_w 
     1                                           / sum_odrad
                      else
                          cloud_rad_w(ialt,jazi) = 0.
                      endif
                      cloud_rad_c(:,ialt,jazi) 
     1                                = sum_odrad_c(:) / cvr_path_sum 
                      rad_last = r_cloud_rad(ialt,jazi) 
                      rad_last_w = cloud_rad_w(ialt,jazi) 
                      sum_odrad_c_last = cloud_rad_c(:,ialt,jazi)

                      xcos_last = upxrad_3d_m
                      ycos_last = upyrad_3d_m

                      cvr_path_sum_sp_w(:) = cvr_path_sum_sp(:)

                    elseif(cvr_path_sum      .gt.  bks_thr .AND.
     1                     cvr_path_sum_last .le.  bks_thr .AND.
     1                     cvr_path_sum_last .gt.  0.           )then
!                     Find average value right at tau_thr boundary
                      frac_path = (bks_thr - cvr_path_sum_last)
     1                          / (cvr_path_sum - cvr_path_sum_last)

                      rad_new = sum_odrad      / cvr_path_sum
                      if(sum_odrad .gt. 0.)then
                          rad_w_new = sum_odrad_w   / sum_odrad
                      else
                          rad_w_new = 0.
                      endif

                      r_cloud_rad(ialt,jazi) = 
     1                                 rad_last         * (1.-frac_path)
     1                + (sum_odrad      / cvr_path_sum) * (   frac_path)
                      cloud_rad_w(ialt,jazi) = 
     1                                 rad_last_w       * (1.-frac_path)
     1                + (rad_w_new)                     * (   frac_path)
                      cloud_rad_c(:,ialt,jazi) = 
     1                                 sum_odrad_c_last * (1.-frac_path)
     1                + (sum_odrad_c(:) / cvr_path_sum) * (   frac_path)
                      if(idebug .eq. 1)then
                        write(6,71)rad_last,rad_new 
     1                           ,frac_path,r_cloud_rad(ialt,jazi)
     1                           ,cloud_rad_w(ialt,jazi)
71                      format(' last/new/frac/rad/radw',5f9.4)
                      endif

                      cvr_path_sum_sp_w(:)
     1                       = cvr_path_sum_sp_w(:) * (1.-frac_path)
     1                       + cvr_path_sum_sp(:)   *  frac_path

!                     xcosup = ???
!                     ycosup = ???

                    endif ! near tau_thr boundary

                  endif ! cvr_path > 0

!                 Calculated weighted value of airmass to cloud
                  taucloud = clwc2alpha*cvr_path_sum
                  sum_am2cld_atten = sum_am2cld_atten
     1                             + airmass2 * trans(taucloud) 
!                 airmass_2_cloud_3d(ialt,jazi) = sum_am2cld_atten

!                 Determined earlier that we are hitting topo
                  if(ihit_topo .eq. 1)then

!                     Land illumination related to terrain slope
                      if(sol_alt(inew_mb,jnew_m)  .gt. 0. )then  

!                       frac_ghi_dir = 0.9 *
!    1                               sind(sol_alt(inew_mb,jnew_m))**0.5

!                       write(6,*)'ghi',inew_mb,jnew_m
!    1                            ,ghi_2d(inew_mb,jnew_m)
                        if(ghi_2d(inew_mb,jnew_m) .gt. 0.)then
                          frac_ghi_dir = 
     1                   (ghi_2d(inew_mb,jnew_m)-dhi_2d(inew_mb,jnew_m))
     1                    / ghi_2d(inew_mb,jnew_m)
                        else
                          write(6,*)' WARNING: ghi_2d = 0.'
     1                             ,inew_mb,jnew_m
     1                             ,sol_alt(inew_mb,jnew_m)
     1                             ,lat(inew_mb,jnew_m)
     1                             ,lon(inew_mb,jnew_m)
                          frac_ghi_dir = 1.0
                        endif

                        alt_norm_int = sum(bi_coeff(:,:) 
     1                               * alt_norm(i1:i2,j1:j2))
                        if(alt_norm_int .gt. 0. )then
                          solar_corrb = sind(alt_norm_int) 
     1                                / sind(sol_alt(inew_mb,jnew_m))
                          solar_corr = (solar_corrb*frac_ghi_dir) + 
     1                                 1.0 * (1. - frac_ghi_dir) 
                          solar_corr = min(max(solar_corr,0.2),2.0) 
                        else ! land is shadowed (all indirect)
                          solar_corr = 0.0 *     frac_ghi_dir 
     1                               + 1.0 * (1.-frac_ghi_dir)
                          solar_corrb = 0.
                        endif
                      else ! sun is down     
                          solar_corr = 1.0
                          solar_corrb = 0.
                      endif

                      topo_gti(ialt,jazi) = ! instead of swi_2d
     1                  sum(bi_coeff(:,:) * ghi_2d(i1:i2,j1:j2))
     1                                    * solar_corr

!                     Terrain spectral radiance (normalized)
!                     Should gtic = dtic + btic?
                      do ic = 1,nc
                          gtic(ic,ialt,jazi) = 
     1                      sum(bi_coeff(:,:) * ghic_2d(ic,i1:i2,j1:j2))
     1                                        * solar_corr

!                         City lights on the ground (spec exitance - emic)
                          if(sol_alt(inew_mb,jnew_m) .lt. twi_alt)then  

!                           'gnd_glow' is nL, consider w/m2/sr/nm
!                           'emic' is relative to solar isotropic
                            emic(ic,ialt,jazi) = sum(bi_coeff(:,:)    
     1                        * gnd_glow(i1:i2,j1:j2)) * city_colrat(ic)
     1                        / 3e9 ! nl to solar isotropic (i.e. day_int)
                            if(gnd_glow(inew_mb,jnew_m) .gt. 0. .AND. ! .OR. 
     1                                                idebug .eq. 1)then      
                              write(6,82)ialt,jazi,inew_mb,jnew_m
     1                                  ,gnd_glow(inew_mb,jnew_m)
 82                           format(' gnd glow used for emic',4i5
     1                              ,f9.1)
                            endif
                          else
                            emic(ic,ialt,jazi) = 0.                   
                          endif

                          dtic(ic,ialt,jazi) = 
     1                      sum(bi_coeff(:,:) * dhic_2d(ic,i1:i2,j1:j2))

!                         Could be computed from 'bhic' or 'bnic'
                          btic(ic,ialt,jazi) = 
     1                      sum(bi_coeff(:,:) * bhic_2d(ic,i1:i2,j1:j2))
     1                                        * solar_corrb

                          if(.true.)then
                             gtic(ic,ialt,jazi) = btic(ic,ialt,jazi)
     1                                          + dtic(ic,ialt,jazi)
                          endif

                          if(idebug .eq. 1 .and. ic .eq. 2)then
                             if(altray .eq. -90.)then
                                write(6,*)' nadir gtic check (observer)'
                             endif
                             write(6,86)gtic(ic,ialt,jazi)
     1                         ,i1,j1,bnic_2d(ic,i1,j1),bni_clr(ic)
     1                         ,ghic_2d(ic,i1,j1)
     1                         ,bhic_2d(ic,i1,j1)
     1                         ,dhic_2d(ic,i1,j1),frac_dir_a(i1,j1)
     1                         ,transm_tn_a(i1,j1),solar_corr
86                           format(' gtic check ',f9.4,2i5,8f9.5)                          
                             write(6,87)inew_mb,jnew_m
     1                         ,sol_alt(inew_mb,jnew_m)
     1                         ,lat(inew_mb,jnew_m)
     1                         ,lon(inew_mb,jnew_m)
     1                         ,ghi_2d(inew_mb,jnew_m)
     1                         ,dhi_2d(inew_mb,jnew_m)
     1                         ,alt_norm_int
     1                         ,solar_corr
87                           format(' solar_corr ',2i5,7f10.4)
                          endif

                      enddo ! ic                

!                     Topo albedo
                      do ic = 1,nc
                        topo_albedo(ic,ialt,jazi) = sum(bi_coeff(:,:) 
     1                                 * topo_albedo_2d(ic,i1:i2,j1:j2))
                      enddo ! ic
!                     topo_albedo(:,ialt,jazi) = 
!    1                    topo_albedo_2d(:,inew_mb,jnew_m)

                      topo_ri(ialt,jazi) = rinew_m
                      topo_rj(ialt,jazi) = rjnew_m

!                     Test for first segment
                      if(dxy1_l .eq. 0. .AND. frac_step_topo .eq. 0.
     1                                  .AND. htstart .eq. topo_sfc)then
                          airmass_2_topo_3d(ialt,jazi) = 1e-5
                          aod_2_topo(ialt,jazi)        = 1e-4
                          dist_2_topo(ialt,jazi)       = 1d0
                          sum_aod_ill                  = 1e-4 
     1                           * transm_3d(inew_mb,jnew_m,int(rk_m)+1)
                          sum_aod_ill_dir              = 1e-4 
     1                           * transm_3d(inew_mb,jnew_m,int(rk_m)+1)
                          sum_aod_ill_opac             = 1e-4 
     1                           * transm_3d(inew_mb,jnew_m,int(rk_m)+1)
                          sum_aod_ill_opac_potl        = 1e-4 
                      else
                          airmass_2_topo_3d(ialt,jazi)
     1                                = airmass1_l * (1.-frac_step_topo) 
     1                                + airmass1_h * frac_step_topo
                          aod_2_topo(ialt,jazi) = sum_aod
     1                                   - aod_inc * (1.-frac_step_topo)
                          dist_2_topo(ialt,jazi) = dslant1_h 
     1                                  - slant2  * (1d0-frac_step_topo)
!    1                        sqrt(dxy1_h**2 + dz1_h**2)
                          topo_ri(ialt,jazi)
     1                                = rinew_l * (1.-frac_step_topo) 
     1                                + rinew_h * frac_step_topo
                          topo_rj(ialt,jazi)
     1                                = rjnew_l * (1.-frac_step_topo) 
     1                                + rjnew_h * frac_step_topo
                      endif

                      if(idebug .eq. 1)then
                          write(6,91)ht_m,topo_m,frac_step_topo
     1                           ,topo_ri(ialt,jazi),topo_rj(ialt,jazi)
     1                           ,solar_corr  
     1                           ,topo_gti(ialt,jazi)
     1                           ,topo_albedo(:,ialt,jazi)
     1                           ,aod_2_topo(ialt,jazi)
     1                           ,dist_2_topo(ialt,jazi)
     1                           ,sum_aod_ill_opac     
     1                           ,sum_aod_ill_opac_potl 
     1                           ,airmass_2_topo_3d(ialt,jazi)
     1                           ,gnd_glow(inew_mb,jnew_m)
     1                           ,ghi_2d(inew_mb,jnew_m)
     1                           ,dhi_2d(inew_mb,jnew_m)
                      endif
91                    format(' Hit topo',2f8.1,4f8.3,f8.2,4f8.3,f12.0
     1                                  ,4f8.3,2f8.1)
                  endif ! hit topo

                  if(idebug .eq. 1)then
                      cloud_od(ialt,jazi) = clwc2alpha*cvr_path_sum   
                      cloud_od_sp(ialt,jazi,:) 
     1                                = clwc2alpha*cvr_path_sum_sp(:)   
                      r_cloud_3d(ialt,jazi) 
     1                                = 1.-(exp(-cloud_od(ialt,jazi)))
!                     write(6,*)' dslant1_h = ',dslant1_h,slant2
                      write(6,101)dz1_l,dz1_h,dxy1_l,dxy1_h
     1                     ,cslant
     1                     ,rinew_h,rjnew_h,rk,ht_m,topo_m 
     1                     ,cvr_path*1e3
     1                     ,clwc_m*1e3,cice_m*1e3,rain_m*1e3,snow_m*1e3
     1                     ,slant2,cvr_path_sum
     1                     ,r_cloud_3d(ialt,jazi)
     1                     ,sum_am2cld_atten              
     1                     ,sum_clrrad,sum_clrrad_pot,airmass1_h
     1                     ,r_cloud_rad(ialt,jazi)
     1                     ,cloud_rad_w(ialt,jazi)
     1                     ,aero_ext_coeff
     1                     ,transm_3d_m
     1                     ,sum_aod,sum_aod_ill_opac
     1                     ,sum_aod_ill_opac_potl
101                   format(2f13.1,2f10.1,a1,f6.1,f7.1,f6.2,2f8.1
     1                  ,1x,f7.4,2x,4f7.4,f10.1,2f11.4,4f8.3,2f8.4,f10.5
     1                  ,4f7.3)
                  endif ! idebug

                 else ! outside horizontal domain
                  ioutside_domain = 1
                  sum_clrrad     = sum_clrrad     + airmass2  
     1                           * frac_fntcloud
                  sum_clrrad_pot = sum_clrrad_pot + airmass2  
     1                           * frac_fntcloud
!                 sum_aod_ill_dir = sum_aod_ill_dir 
!    1                       + aero_ext_coeff * slant2 * 1.0
!    1                                        * frac_fntcloud
                 endif ! in horizontal domain

                else  ! outside vertical domain
!                ioutside_domain = 1
                 if(.not.(iabove .eq. 1 .and. rk_m .gt. float(nk) 
     1                                  .and. rk_h .lt. rk_m)    )then
                     ihit_bounds = 1
                 endif
                endif ! in vertical domain

                if(outside_domain .eq. 1)then
                  if(idebug .eq. 1 .and. 
     1               (iwrite .le. 10 .or. ls .eq. (ls/40)*40) )then
                      iwrite = iwrite + 1
                      write(6,111)ialt,jazi,ls,rinew_h,rjnew_h,rk,slant2
     1                    ,dxy1_l,dxy1_h,dz1_h,iabove,ihit_bounds
111                   format(' out of bounds ',2i5,i6,4f9.3,3f10.1,2i3)
                      if(iabove .eq. 1)then
                        write(6,112)rk_m,rk_h,sum_clrrad,airmass1_h
112                     format('rk_m/rk_h/sum_clrrad/airmass1_h',4f9.3)
                      endif
                      if(rkdelt .ne. 0.)then
                        write(6,*)' aterm/bterm/cterm/discriminant='
     1                             ,aterm,bterm,cterm,discriminant   
!                       write(6,*)' argd1/argd2=',argd1,argd2
                      endif
                  endif
                endif

                if(rk_h - rkstart .gt. 0.0)then
!                 rkdelt = rkdelt2
!                 rkdelt = 0.5 * (rk_h - rkstart)
                  rkdelt = rkdelt * 1.1             
                  rkdelt = min(rkdelt,rkdelt2)
                endif
              enddo ! while kk/rk

!             Extra verbose
              if(idebug .eq. 1 .and.
     1                (altray .lt. 0. .or. mode_aero_cld .eq. 3) )then
                itrace = nint(trace_ri(ialt,jazi))
                jtrace = nint(trace_rj(ialt,jazi))
                if(itrace .ge. 1 .and. itrace .le. ni .AND.
     1             jtrace .ge. 1 .and. jtrace .le. nj)then
                  trsolalt = sol_alt(itrace,jtrace)
                else
                  trsolalt = -99.
                endif
                write(6,113)ls,rk,ht_h,trace_ri(ialt,jazi)
     1                                ,trace_rj(ialt,jazi),trsolalt
113             format('   ls/rk/ht/slalt = ',i4,f8.3,f12.1,' eor trij'
     1                                       ,3f8.2) 
                write(6,114)ls,rk,rk_h,rk_m,ht_h,ihit_topo,ihit_bounds
     1                     ,ioutside_domain
     1                     ,iabove,gc_deg,dslant1_h,dz1_h,tlat,tlon
     1                     ,cvr_path_sum,idebug
114             format(' end of ray ',
     1      'ls/rk3/ht/ihit-tbo/iab/gc/sl1h/dz1h/latlon/cvr/idebug'
     1              ,i6,3f10.3,f10.1,4i3,f9.3,f13.1,f13.0,2f9.2,f7.2,i3)
              endif        

          endif ! true

!         Store various accumulated sums
          cloud_od(ialt,jazi) = clwc2alpha*cvr_path_sum   
          cloud_od_sp(ialt,jazi,:) = clwc2alpha*cvr_path_sum_sp(:)   
          cloud_od_sp_w(ialt,jazi,:) = clwc2alpha*cvr_path_sum_sp_w(:)   

          if(mode_aero_cld .lt. 3 .or. .true.)then ! prevent dble counting
              r_cloud_3d(ialt,jazi) = 1.-(exp(-cloud_od(ialt,jazi)))
              if(mode_aero_cld .eq. 3 .and. .false.)then ! 
                  frac_aod = sum_aod / cloud_od(ialt,jazi)
              endif
          else
              r_cloud_3d(ialt,jazi)
     1                    = 1.-(exp(-(cloud_od(ialt,jazi)+sum_aod)))
          endif

           airmass_2_cloud_3d(ialt,jazi) = sum_am2cld_atten

!          Account for multiple scattering in aod_ill
           aod_ill(ialt,jazi) = (1.0 - crep_thr) * sum_aod_ill 
     1                               + crep_thr  * sum_aod
!          aod_ill(ialt,jazi) = sum_aod_ill
           if(aod_ill(ialt,jazi) .lt. 0.)then
              write(6,*)' ERROR aod_ill < 0',ialt,jazi,altray
     1                 ,view_azi_deg,sum_aod,sum_aod_ill
              stop
           endif

           aod_ill_dir(ialt,jazi) = min(sum_aod_ill_dir,1e30)
           aod_tot(ialt,jazi) = sum_aod
           aod_ill_opac(ialt,jazi) = sum_aod_ill_opac
           aod_ill_opac_potl(ialt,jazi) = sum_aod_ill_opac_potl
           clear_rad_c_nt(:,ialt,jazi) = sum_aod_rad_opac(:) / (4.*pi)

           if(r_cloud_3d(ialt,jazi) .gt. .5)then
              icloud = 1
           else
              icloud = 0
           endif

           icloud_tot = icloud_tot + icloud

!          l_process(ialt,jazi) = .true.

!          include '../lib/cloud/skyglow_phys.inc'

!          if(sol_alt(i,j) .gt. 0.)then
           if(obj_alt(i,j) .gt. twi_0)then ! experiment

!           Get clear sky daylight brightness ratio at point
!           (i.e. fraction of atmosphere illuminated by the sun)
            if(rkstart .gt. float(nk) .and. sum_clrrad .eq. 0.0)then
              clear_radf_c(:,ialt,jazi) = 1.00 ! ray all above domain
            elseif(rkstart .gt. float(nk-1) .and. 
     1             view_altitude_deg .ge. -horz_dep_d)then
              clear_radf_c(:,ialt,jazi) = 1.00 ! correct inaccuracy
            elseif(sum_clrrad_pot .gt. 0.)then
              clear_radf_c(:,ialt,jazi) = crep_thr ! secondary scattering in cloud shadow
     1                  + (1.0 - crep_thr) * (sum_clrrad/sum_clrrad_pot)
            else ! entire ray in Earth's shadow
              clear_radf_c(:,ialt,jazi) = 1.0
            endif
            if(idebug .eq. 1)then
              write(6,119)rkstart,nk,view_altitude_deg,horz_dep_d
     1             ,cloud_od(ialt,jazi)
     1             ,clear_radf_c(icd,ialt,jazi)
     1             ,aod_ill_opac(ialt,jazi)/aod_ill_opac_potl(ialt,jazi)
     1             ,sum_clrrad,sum_clrrad_pot,airmass1_h,dz1_h
 119          format(
     1          ' rks/nk/alt/dep/cod/radf/aodf/smclrrd/pot/am/dz1_h'
     1              ,f11.2,i3,2f9.4,2f9.3,3f9.3,f13.1,f8.2)
            endif

           endif

!          end include 'skyglow.inc'

          endif ! l_pix_trace is TRUE
         enddo ! jazi

!        I4_elapsed = ishow_timer()

!        Get clear sky twilight brightness in ring of constant altitude
!        Approximate method used only in deep twilight
!        A different condition could be added for very high altitudes
!        by assessing range of solar altitudes around the limb
!        if(sol_alt(i,j) .le. 0.)then
         if(sol_alt(i,j) .le. twi_0 .and. 
     1      sol_alt(i,j) .ge. (-horz_dep_d - 18.0))then ! narrower range
!            write(6,*)' args: ',l_solar_eclipse,i4time,rlat,rlon
             call skyglow_phys_twi(ialt,ialt,1,minazi,maxazi,jazi_delt
     1             ,minalt,maxalt,minazi,maxazi,idebug_a
     1             ,sol_alt(i,j),sol_azi(i,j),view_alt,view_az
     1             ,earth_radius,patm,aod_vrt,aod_ray_eff,aod_ray_dir
     1             ,aero_scaleht,htstart,redp_lvl     ! I
     1             ,aod_ill                           ! I (dummy)
     1             ,l_solar_eclipse,i4time,rlat,rlon  ! I
     1             ,clear_radf_c,horz_dep_d,ag_2d     ! I (ag2d is dummy)
     1             ,clear_rad_c,elong_p             ) ! O (elongp is dummy)
         endif

!        Pixels per degree times 1 and times 2
!        if(jazi_delt .eq. 2 .OR. jazi_delt .eq. 4 .OR. 
!    1      jazi_delt .eq. 5 .OR. jazi_delt .eq. 8 .OR.
!    1      jazi_delt .eq. 10 .OR. jazi_delt .eq. 16 .OR.
!    1      jazi_delt .eq. 20 .OR. jazi_delt .eq. 32 .OR.
!    1      jazi_delt .eq. 64 .OR. jazi_delt .eq. 128 .OR.
!    1      jazi_delt .eq. 256 
!    1                                       )then ! fill missing azimuths

        if(jazi_delt_eff .gt. 1 .and.
     1     jazi_delt_eff .lt. maxazi-minazi  )then ! fill missing azimuths

          do jazi = minazi,maxazi
            call get_interp_parms(minazi,maxazi,jazi_delt_eff,jazi   ! I
     1                           ,fm,fp,jazim,jazip,ir,istatus)      ! O
            if(istatus .ne. 1)then
              write(6,*)' ERROR in jazi call: minazi,maxazi,jazi'
              stop
            endif

            if(ir .ne. 0 .and. (.not. l_pix_trace(ialt,jazi)))then
              r_cloud_3d(ialt,jazi) = 
     1         fm * r_cloud_3d(ialt,jazim) + fp * r_cloud_3d(ialt,jazip)
              cloud_od(ialt,jazi) = 
     1         fm * cloud_od(ialt,jazim)   + fp * cloud_od(ialt,jazip)
              cloud_od_sp(ialt,jazi,:) = 
     1         fm * cloud_od_sp(ialt,jazim,:) 
     1                                  + fp * cloud_od_sp(ialt,jazip,:)
              cloud_od_sp_w(ialt,jazi,:) = 
     1         fm * cloud_od_sp_w(ialt,jazim,:) 
     1                                + fp * cloud_od_sp_w(ialt,jazip,:)
              r_cloud_rad(ialt,jazi) = 
     1         fm * r_cloud_rad(ialt,jazim)+fp * r_cloud_rad(ialt,jazip)      
              cloud_rad_c(:,ialt,jazi) = 
     1         fm * cloud_rad_c(:,ialt,jazim)  
     1                                 + fp *cloud_rad_c(:,ialt,jazip)
              cloud_rad_w(ialt,jazi) = 
     1         fm * cloud_rad_w(ialt,jazim)+fp * cloud_rad_w(ialt,jazip)      

              if(sol_alt(i,j) .le. 0.)then
                clear_rad_c(:,ialt,jazi) = 
     1           fm * clear_rad_c(:,ialt,jazim)   
     1                              + fp * clear_rad_c(:,ialt,jazip)
                clear_rad_c_nt(:,ialt,jazi) = 
     1           fm * clear_rad_c_nt(:,ialt,jazim)   
     1                              + fp * clear_rad_c_nt(:,ialt,jazip)
              endif

              clear_radf_c(:,ialt,jazi) = 
     1         fm * clear_radf_c(:,ialt,jazim)   
     1                                 + fp * clear_radf_c(:,ialt,jazip)  

!             if(altray .eq. -15.0)then
!               write(6,121)jazi,clear_radf_c(:,ialt,jazi)
!    1            ,clear_radf_c(:,ialt,jazim),clear_radf_c(:,ialt,jazip)
121             format('interpa',i5,9f12.0)            
!             endif

              airmass_2_cloud_3d(ialt,jazi) = 
     1                fm * airmass_2_cloud_3d(ialt,jazim) 
     1              + fp * airmass_2_cloud_3d(ialt,jazip)
              airmass_2_topo_3d(ialt,jazi) = 
     1                fm * airmass_2_topo_3d(ialt,jazim) 
     1              + fp * airmass_2_topo_3d(ialt,jazip)
              aod_ill(ialt,jazi) = 
     1                fm * aod_ill(ialt,jazim) 
     1              + fp * aod_ill(ialt,jazip)
              aod_ill_dir(ialt,jazi) = 
     1                fm * aod_ill_dir(ialt,jazim) 
     1              + fp * aod_ill_dir(ialt,jazip)
              aod_ill_opac(ialt,jazi) = 
     1                fm * aod_ill_opac(ialt,jazim) 
     1              + fp * aod_ill_opac(ialt,jazip)
              aod_ill_opac_potl(ialt,jazi) = 
     1                fm * aod_ill_opac_potl(ialt,jazim) 
     1              + fp * aod_ill_opac_potl(ialt,jazip)
              aod_2_cloud(ialt,jazi) = 
     1                fm * aod_2_cloud(ialt,jazim) 
     1              + fp * aod_2_cloud(ialt,jazip)
              aod_2_topo(ialt,jazi) = 
     1                fm * aod_2_topo(ialt,jazim) 
     1              + fp * aod_2_topo(ialt,jazip)
              dist_2_topo(ialt,jazi) = 
     1                fm * dist_2_topo(ialt,jazim) 
     1              + fp * dist_2_topo(ialt,jazip)
              aod_tot(ialt,jazi) = 
     1                fm * aod_tot(ialt,jazim) 
     1              + fp * aod_tot(ialt,jazip)
              topo_gti(ialt,jazi) = 
     1                fm * topo_gti(ialt,jazim) 
     1              + fp * topo_gti(ialt,jazip)
              gtic(:,ialt,jazi) = 
     1                fm * gtic(:,ialt,jazim) 
     1              + fp * gtic(:,ialt,jazip)
              dtic(:,ialt,jazi) = 
     1                fm * dtic(:,ialt,jazim) 
     1              + fp * dtic(:,ialt,jazip)
              btic(:,ialt,jazi) = 
     1                fm * btic(:,ialt,jazim) 
     1              + fp * btic(:,ialt,jazip)
              emic(:,ialt,jazi) = 
     1                fm * emic(:,ialt,jazim) 
     1              + fp * emic(:,ialt,jazip)
              topo_albedo(:,ialt,jazi) = 
     1                fm * topo_albedo(:,ialt,jazim) 
     1              + fp * topo_albedo(:,ialt,jazip)

              if(topo_ri(ialt,jazim) .ne. 0. .and.
     1           topo_ri(ialt,jazip) .ne. 0.      )then 

                if(abs(topo_ri(ialt,jazim) -
     1                 topo_ri(ialt,jazip)) .lt. rimid)then
                  topo_ri(ialt,jazi) = 
     1                fm * topo_ri(ialt,jazim) 
     1              + fp * topo_ri(ialt,jazip)
                else ! cyclic condition
                  if(topo_ri(ialt,jazip) .gt. topo_ri(ialt,jazim))then
                    argm = topo_ri(ialt,jazim) + float(ni)
                    argp = topo_ri(ialt,jazip)
                  else
                    argm = topo_ri(ialt,jazim)
                    argp = topo_ri(ialt,jazip) + float(ni)
                  endif
                  topo_ri(ialt,jazi) =
     1              modulo(fm * argm + fp * argp-0.5, float(ni)) + 0.5
                endif

                topo_rj(ialt,jazi) = 
     1                fm * topo_rj(ialt,jazim) 
     1              + fp * topo_rj(ialt,jazip)
              else
                topo_ri(ialt,jazi) = 0.
                topo_rj(ialt,jazi) = 0.
              endif

              if(abs(trace_ri(ialt,jazim) -
     1               trace_ri(ialt,jazip)) .lt. rimid)then
                trace_ri(ialt,jazi) = 
     1                fm * trace_ri(ialt,jazim) 
     1              + fp * trace_ri(ialt,jazip)
              else ! cyclic condition
                if(trace_ri(ialt,jazip) .gt. trace_ri(ialt,jazim))then
                  argm = trace_ri(ialt,jazim) + float(ni)
                  argp = trace_ri(ialt,jazip)
                else
                  argm = trace_ri(ialt,jazim)
                  argp = trace_ri(ialt,jazip) + float(ni)
                endif
                trace_ri(ialt,jazi) =
     1              modulo(fm * argm + fp * argp-0.5, float(ni)) + 0.5
              endif

              trace_rj(ialt,jazi) = 
     1                fm * trace_rj(ialt,jazim) 
     1              + fp * trace_rj(ialt,jazip)
            else
              if(cloud_rad_c(1,ialt,jazi) .eq. .250 .or. 
     1           cloud_rad_c(1,ialt,jazi) .eq. .500)then
                write(6,201)ialt,jazi,view_alt(ialt,jazi)
     1                               ,view_az(ialt,jazi)          
     1                               ,cloud_rad_c(1,ialt,jazi)
 201            format(' cloud_rad_c before interp:',2i5,2f9.2,f9.3)
              endif

            endif ! ir .ne. 0
          enddo ! jazi interp
         endif ! fill azimuth

        enddo ! ialt

        do ialt = minalt,maxalt
        do jazi = minazi,maxazi
          if(cloud_rad_c(1,ialt,jazi) .eq. .250 .or. 
     1       cloud_rad_c(1,ialt,jazi) .eq. .500)then
            write(6,202)ialt,jazi,view_alt(ialt,jazi)
     1                           ,view_az(ialt,jazi)          
     1                           ,cloud_rad_c(1,ialt,jazi)
 202        format(' cloud_rad_c before alt interp:',2i5,2f9.2,f9.3)
          endif
        enddo ! jazi
        enddo ! ialt

!       We normally fill missing rings above 20 degrees 
        call get_idx(20.,minalt,alt_scale,ialt_min)
        call get_idx(maxalt_deg,minalt,alt_scale,ialt_max)

        I4_elapsed = ishow_timer()

        if(htagl .gt. 30000e3)then ! skip fill for very high altitudes
          ialt_min = minalt+1
          ialt_max = minalt
          write(6,*)' Skip filling missing data in altitude rings'
        else
          write(6,*)' Filling missing data in altitude rings'
          write(6,*)' ialtmin/ialtmax/deg = '
     1               ,ialt_min,ialt_max,maxalt_deg
        endif

        do ialt = ialt_min,ialt_max ! fill in missing alt rings
          call get_interp_parms(minalt,maxalt,idelt,ialt           ! I
     1                         ,fm,fp,ialtm,ialtp,ir,istatus)      ! O
          if(istatus .ne. 1)then
            write(6,*)' ERROR in ialt call: minalt,maxalt,ialt'
            stop
          endif

          if(ir .ne. 0)then
            r_cloud_3d(ialt,:) =
     1       fm * r_cloud_3d(ialtm,:) + fp * r_cloud_3d(ialtp,:)
            cloud_od(ialt,:) =
     1       fm * cloud_od(ialtm,:)   + fp * cloud_od(ialtp,:)
            cloud_od_sp(ialt,:,:) =
     1       fm * cloud_od_sp(ialtm,:,:) + fp * cloud_od_sp(ialtp,:,:)
            cloud_od_sp_w(ialt,:,:) =
     1                                fm * cloud_od_sp_w(ialtm,:,:)
     1                              + fp * cloud_od_sp_w(ialtp,:,:)
            r_cloud_rad(ialt,:) =
     1       fm * r_cloud_rad(ialtm,:)  + fp * r_cloud_rad(ialtp,:)
            cloud_rad_c(:,ialt,:) =
     1       fm * cloud_rad_c(:,ialtm,:) + fp * cloud_rad_c(:,ialtp,:)
            cloud_rad_w(ialt,:) =
     1       fm * cloud_rad_w(ialtm,:)  + fp * cloud_rad_w(ialtp,:)
            clear_rad_c(:,ialt,:) =
     1       fm * clear_rad_c(:,ialtm,:) + fp * clear_rad_c(:,ialtp,:)
            clear_rad_c_nt(:,ialt,:) =
     1                                fm * clear_rad_c_nt(:,ialtm,:) 
     1                              + fp * clear_rad_c_nt(:,ialtp,:)
            clear_radf_c(:,ialt,:) =
     1       fm * clear_radf_c(:,ialtm,:) + fp * clear_radf_c(:,ialtp,:)
            airmass_2_cloud_3d(ialt,:) =
     1           fm * airmass_2_cloud_3d(ialtm,:) 
     1         + fp * airmass_2_cloud_3d(ialtp,:)
            airmass_2_topo_3d(ialt,:) =
     1           fm * airmass_2_topo_3d(ialtm,:) 
     1         + fp * airmass_2_topo_3d(ialtp,:)
            aod_ill(ialt,:) =
     1           fm * aod_ill(ialtm,:) 
     1         + fp * aod_ill(ialtp,:)
            aod_ill_dir(ialt,:) =
     1           fm * aod_ill_dir(ialtm,:) 
     1         + fp * aod_ill_dir(ialtp,:)
            aod_ill_opac(ialt,:) =
     1           fm * aod_ill_opac(ialtm,:) 
     1         + fp * aod_ill_opac(ialtp,:)
            aod_ill_opac_potl(ialt,:) =
     1           fm * aod_ill_opac_potl(ialtm,:) 
     1         + fp * aod_ill_opac_potl(ialtp,:)
            aod_2_cloud(ialt,:) =
     1           fm * aod_2_cloud(ialtm,:) 
     1         + fp * aod_2_cloud(ialtp,:)
            aod_2_topo(ialt,:) =
     1           fm * aod_2_topo(ialtm,:) 
     1         + fp * aod_2_topo(ialtp,:)
            dist_2_topo(ialt,:) =
     1           fm * dist_2_topo(ialtm,:) 
     1         + fp * dist_2_topo(ialtp,:)
            aod_tot(ialt,:) =
     1           fm * aod_tot(ialtm,:) 
     1         + fp * aod_tot(ialtp,:)
            topo_gti(ialt,:) =
     1         fm * topo_gti(ialtm,:) + fp * topo_gti(ialtp,:)
            gtic(:,ialt,:) =
     1         fm * gtic(:,ialtm,:) + fp * gtic(:,ialtp,:)
            dtic(:,ialt,:) =
     1         fm * dtic(:,ialtm,:) + fp * dtic(:,ialtp,:)
            btic(:,ialt,:) =
     1         fm * btic(:,ialtm,:) + fp * btic(:,ialtp,:)
            emic(:,ialt,:) =
     1         fm * emic(:,ialtm,:) + fp * emic(:,ialtp,:)
            topo_albedo(:,ialt,:) =
     1         fm * topo_albedo(:,ialtm,:) + fp * topo_albedo(:,ialtp,:)

!           Check for zero values
            topo_ri(ialt,:) =
     1         fm * topo_ri(ialtm,:) + fp * topo_ri(ialtp,:)
            topo_rj(ialt,:) =
     1         fm * topo_rj(ialtm,:) + fp * topo_rj(ialtp,:)

            trace_ri(ialt,:) =
     1         fm * trace_ri(ialtm,:) + fp * trace_ri(ialtp,:)
            trace_rj(ialt,:) =
     1         fm * trace_rj(ialtm,:) + fp * trace_rj(ialtp,:)
          endif ! ir .ne. 0
        enddo ! ialt

        I4_elapsed = ishow_timer()

        write(6,*)' Number of rays with cloud = ',icloud_tot
     1           ,' out of ',91*361

        write(6,*)' Range of r_cloud_3d = ',minval(r_cloud_3d)
     1                                     ,maxval(r_cloud_3d)

        write(6,*)' Range of cloud_od = ',minval(cloud_od)
     1                                   ,maxval(cloud_od)

        do isp = 1,nsp
          write(6,*)' Range of cloud_od_sp = ',isp
     1                                    ,minval(cloud_od_sp(:,:,isp))
     1                                    ,maxval(cloud_od_sp(:,:,isp))
        enddo ! isp

        write(6,*)' Range of r_cloud_rad 1 = ',minval(r_cloud_rad)
     1                                        ,maxval(r_cloud_rad)

        where(r_cloud_rad .gt. 1e30)r_cloud_rad = 0.

        write(6,*)' Range of r_cloud_rad 2 = ',minval(r_cloud_rad)
     1                                        ,maxval(r_cloud_rad)

        write(6,*)' Range of cloud_rad_c R =',minval(cloud_rad_c(1,:,:))
     1                                       ,maxval(cloud_rad_c(1,:,:))

        write(6,*)' Range of cloud_rad_c G =',minval(cloud_rad_c(2,:,:))
     1                                       ,maxval(cloud_rad_c(2,:,:))

        write(6,*)' Range of cloud_rad_c B =',minval(cloud_rad_c(3,:,:))
     1                                       ,maxval(cloud_rad_c(3,:,:))

        write(6,*)' Range of clear_rad_c 3 =',minval(clear_rad_c(3,:,:))
     1                                       ,maxval(clear_rad_c(3,:,:))

        write(6,*)' Range of clear_rad_c_nt 3 (wm2srnm) ='
     1                                    ,minval(clear_rad_c_nt(3,:,:))
     1                                    ,maxval(clear_rad_c_nt(3,:,:))

!       cloud_rad_w = min(cloud_rad_w,1.)

        write(6,*)' Range of cloud_rad_w = ',minval(cloud_rad_w)
     1                                      ,maxval(cloud_rad_w)

        write(6,*)' Range of clear_radf_c 3 ='
     1                                      ,minval(clear_radf_c(3,:,:))
     1                                      ,maxval(clear_radf_c(3,:,:))

        call get_idx(20.,minalt,alt_scale,ialt_idx)
        if(ialt_idx .gt. minalt .and. ialt_idx .le. maxalt)then
          write(6,*)' Range of clear_radf_c +20. ='
     1                    ,minval(clear_radf_c(3,ialt_idx,:))  
     1                    ,maxval(clear_radf_c(3,ialt_idx,:))
        endif

        call get_idx(-20.,minalt,alt_scale,ialt_idx)
        if(ialt_idx .gt. minalt .and. ialt_idx .le. maxalt)then
          write(6,*)' Range of clear_radf_c -20. ='
     1                    ,minval(clear_radf_c(3,ialt_idx,:))  
     1                    ,maxval(clear_radf_c(3,ialt_idx,:))
        endif

        write(6,*)' Range of topo_gti = ',minval(topo_gti)
     1                                   ,maxval(topo_gti)

        write(6,*)' Range of gtic 2 = ',minval(gtic(2,:,:))
     1                                 ,maxval(gtic(2,:,:))

        write(6,*)' Range of btic 2 = ',minval(btic(2,:,:))
     1                                 ,maxval(btic(2,:,:))

        write(6,*)' Range of dtic 2 = ',minval(dtic(2,:,:))
     1                                 ,maxval(dtic(2,:,:))

        write(6,*)' Range of emic 2 = ',minval(emic(2,:,:))
     1                                 ,maxval(emic(2,:,:))

        write(6,*)' Range of topo_albedo 2 = '
     1                                       ,minval(topo_albedo(2,:,:))
     1                                       ,maxval(topo_albedo(2,:,:))

        write(6,*)' Range of topo_ri = ',minval(topo_ri)
     1                                  ,maxval(topo_ri)

        write(6,*)' Range of topo_rj = ',minval(topo_rj)
     1                                  ,maxval(topo_rj)

        write(6,*)' Range of aod_tot = ',minval(aod_tot)
     1                                  ,maxval(aod_tot)       

        write(6,*)' Range of aod_ill = ',minval(aod_ill)
     1                                  ,maxval(aod_ill)       

        write(6,*)' Range of aod_ill_dir = ',minval(aod_ill_dir)
     1                                      ,maxval(aod_ill_dir)       

        write(6,*)' Range of aod_ill_opac = ',minval(aod_ill_opac)
     1                                       ,maxval(aod_ill_opac)       

        write(6,*)' Range of aod_ill_opac_potl = '
     1        ,minval(aod_ill_opac_potl),maxval(aod_ill_opac_potl)       

        do ialt = minalt,maxalt
        do jazi = minazi,maxazi
          if(cloud_rad_c(1,ialt,jazi) .eq. .250 .or. 
     1       cloud_rad_c(1,ialt,jazi) .eq. .500)then
            write(6,203)ialt,jazi,view_alt(ialt,jazi)
     1                           ,view_az(ialt,jazi)          
     1                           ,cloud_rad_c(1,ialt,jazi)
 203        format(' cloud_rad_c after interp:',2i5,2f9.2,f9.3)
          endif
        enddo ! jazi 
        enddo ! ialt

        if(int(rkstart)+1 .le. nk)then
            transm_obs = transm_3d(i,j,int(rkstart)+1)
        else
            transm_obs = 1
        endif

        swi_obs = ghi_2d(i,j)

!       Terrain shadowing not yet accounted for in this step
        write(6,*)'transm of observer is ',i,j,int(rkstart)+1
     1                                     ,transm_obs
        write(6,302)i,j,ghi_2d(i,j)
302     format(' GHI of observer is        ',2i5,e14.7,' W/m**2')
        write(6,303)i,j,dhi_2d(i,j)
303     format(' Diffuse HI of observer is ',2i5,e14.7,' W/m**2')
        write(6,304)i,j,ghi_2d(i,j)-dhi_2d(i,j)
304     format(' Direct HI of observer is  ',2i5,e14.7,' W/m**2')
        if(sol_alt(i,j) .gt. 0.)then
            write(6,*)'DNI of observer is        ',i,j
     1                ,(ghi_2d(i,j)-dhi_2d(i,j)) / sind(sol_alt(i,j))
        else
            write(6,*)' topo_gti is derived from GHI'
            if(sol_alt(i,j) .gt. -18.)then
                write(6,*)' diffuse twi of observer is '
     1                   ,difftwi(sol_alt(i,j))
            endif
        endif
        write(6,311)bnic_2d(:,i,j)
311     format(' bnic of observer (rgb) is ',3f12.8)
        write(6,312)bhic_2d(:,i,j)
312     format(' bhic of observer (rgb) is ',3f12.8)
        write(6,313)dhic_2d(:,i,j)
313     format(' dhic of observer (rgb) is ',3f12.8)
        write(6,314)ghic_2d(:,i,j)
314     format(' ghic of observer (rgb) is ',3f12.8)

        solalt_last = sol_alt(i,j)

        write(6,*)' Sample of topo_ri, topo_rj, trace_ri, trace_rj'
        write(6,*)
     1    ' ialt   jazi          alt      azi    topo_ri    topo_rj'
     1   ,'    solalt_tp  trace_ri   trace_rj  solalt_tr dist_2_topo'
        do ialt_sample = minalt,min(minalt+100,maxalt)
!          ialt_sample = min(max(-45*4,minalt),maxalt) ! -5 degrees alt
!          ialt_sample = min(max(-15,minalt),maxalt) ! -5 degrees alt
           jazi_sample = minazi

           itopo = nint(topo_ri(ialt_sample,jazi_sample))
           jtopo = nint(topo_rj(ialt_sample,jazi_sample))
           if(itopo .ge. 1 .and. itopo .le. ni .and.
     1        jtopo .ge. 1 .and. jtopo .le. nj       )then
              solalt_topo = sol_alt(itopo,jtopo)
           else
              solalt_topo = 0.
           endif
           itrace = nint(trace_ri(ialt_sample,jazi_sample))
           jtrace = nint(trace_rj(ialt_sample,jazi_sample))
           if(itrace .ge. 1 .and. itrace .le. ni .and.
     1        jtrace .ge. 1 .and. jtrace .le. nj       )then
              solalt_trace = sol_alt(itrace,jtrace)
           else
              solalt_trace = 0.
           endif
           write(6,321)ialt_sample,jazi_sample
     1          ,view_alt(ialt_sample,jazi_sample)
     1          ,view_az(ialt_sample,jazi_sample)
     1          ,topo_ri(ialt_sample,jazi_sample)
     1          ,topo_rj(ialt_sample,jazi_sample)
     1          ,solalt_topo
     1          ,trace_ri(ialt_sample,jazi_sample)
     1          ,trace_rj(ialt_sample,jazi_sample)
     1          ,solalt_trace
     1          ,dist_2_topo(ialt_sample,jazi_sample)
321        format(2i6,' tp/tr',2f9.2,7f11.4)       
        enddo ! jazi_sample
        
        I4_elapsed = ishow_timer()

        write(6,*)' successful return from get_cloud_rays'
 
        istatus = 1
        return
        end

        subroutine get_topo_info(alt,htmsl,earth_radius,idebug
     1                          ,alt_norm,dist_to_topo)

        include 'trigd.inc'

        real alt                  ! I elevation angle
        real htmsl                ! I observer height MSL
        real earth_radius         ! I earth radius (meters)
        real alt_norm             ! O elevation angle rel to sfc normal
        real dist_to_topo         ! O distance to sfc (meters)

!       Altitude relative to surface normal (emission angle)
        if(alt .ne. 0.)then
          slope = tand(90. - abs(alt))
          htrad = (earth_radius+htmsl) / earth_radius
          c = htrad * slope
          call line_ellipse(slope,c,0.,0.,1.,1.,r_missing_data,x1,x2
     1                                                        ,y1,y2)
          gc = atan2d(+y1,-x1) ! great circle dist to ground pt
          if(alt .gt. 0.)then
            alt_norm = alt - gc
          else
            alt_norm = alt + gc
          endif
          distrad = sqrt((htrad+x1)**2 + y1**2)
          dist_to_topo = distrad * earth_radius

          if(idebug .eq. 1)write(6,1)distrad,htrad,x1,y1
 1        format(' distrad,htrad,x1,y1',4f13.8)
        else
          alt_norm = 0.  
        endif

        return
        end
