     
     subroutine get_cloud_rad_faces2(            &
                obj_alt,obj_azi,horz_dep_d_obs,  & ! I
                solalt,solazi,                   & ! I 
                clwc_3d,cice_3d,rain_3d,snow_3d, & ! I
                topo_a,grdasp,                   & ! I
                ni,nj,nk,idb,jdb,                & ! I
                heights_3d,                      & ! I 
                sfc_glow,                        & ! I
                transm_3d,transm_4d)               ! O

     include 'trigd.inc'

!    Calculate 3D radiation field looping through each of the 6 faces in
!    the domain.

!    This version gives better handling of partial illumination of a
!    grid box when the ray is descending below the terrain. A smoother
!    sky illumination results with more consistent 'aodf' values that
!    are responding to the 'transm' values.

     use mem_namelist, ONLY: r_missing_data,earth_radius,grid_spacing_m &
                            ,aod,aero_scaleht,fcterm,redp_lvl
     use mem_allsky, ONLY: aod_3d   ! (extinction coefficient)            ! I
     use mem_allsky, ONLY: uprad_4d ! (upward spectral irradiance)
     use mem_allsky, ONLY: mode_aero_cld
     use cloud_rad ! Cloud Radiation and Microphysics Parameters
     include 'rad.inc' ! e.g. for ext_g, ext_o, o3_du

     trans(od) = exp(-min(od,80.))
     scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! x/scurve range is 0-1

!    For each face, set start/end points to loop in the i,j,k dimension
!    Index of 1 points to the 1st element, 2 points to the last element
!                  T  B  W  E  N  S
     real i1(6)  / 1, 1, 1, 2, 1, 1/            
     real i2(6)  / 2, 2, 1, 2, 2, 2/           
     real j1(6)  / 1, 1, 1, 1, 2, 1/            
     real j2(6)  / 2, 2, 2, 2, 2, 1/           
     real k1(6)  / 2, 1, 1, 1, 1, 1/
     real k2(6)  / 2, 1, 2, 2, 2, 2/
     character*6 cfaces(6) /'Top','Bottom','West','East','North','South'/

     real heights_3d(ni,nj,nk)
     real clwc_3d(ni,nj,nk)  ! kg/m**3
     real cice_3d(ni,nj,nk)  ! kg/m**3
     real rain_3d(ni,nj,nk)  ! kg/m**3
     real snow_3d(ni,nj,nk)  ! kg/m**3
     real b_alpha_3d(ni,nj,nk) ! m**-1          
     real transm_2d(ni,nj)    ! used for terrain surface (the "old" way) 
     real transm_3d(ni,nj,nk) ! direct transmission plus forward scattered
     real transm_4d(ni,nj,nk,nc) ! color information added + limb shadow

     real obj_alt(ni,nj),obj_azi(ni,nj)
     real topo_a(ni,nj)
     real projrot(ni,nj)
     real grdasp(ni,nj)
     real sfc_glow(ni,nj)
     real tbuff(ni,nj)

     real heights_1d(nk)

     real trans_c(nc)
     real bi_coeff(2,2),tri_coeff(2,2,2),b_alpha_cube(2,2,2)
     equivalence (s,scurr) ! needed only during transition from s to scurr

     integer htlutlo,htluthi
     parameter (htlutlo = -1000)
     parameter (htluthi = 30000)
     real htlut(htlutlo:htluthi)

     logical l_same_point

     clwc2alpha = 1.5 / (rholiq  * reff_clwc)
     cice2alpha = 1.5 / (rholiq  * reff_cice)
     rain2alpha = 1.5 / (rholiq  * reff_rain)
     snow2alpha = 1.5 / (rhosnow * reff_snow)
     pice2alpha = 1.5 / (rhograupel * reff_graupel)

     idebug = 0
     angstrom_exp_a = 2.4 - (fcterm * 15.)

     twi_alt = -4.5

!    Initialize
     transm_3d = r_missing_data
     do k = 1,nk
     do i = 1,ni
     do j = 1,nj
       if(topo_a(i,j) .gt. heights_3d(i,j,k))then
         transm_3d(i,j,k) = 0.
       endif
       if(obj_alt(i,j) .lt. twi_alt)then
         transm_3d(i,j,k) = 0.
       endif
     enddo ! j
     enddo ! i
     enddo ! k

     transm_4d = 0.                   

     write(6,*)' transm_3d column 1 = ',transm_3d(idb,jdb,:)

     if(mode_aero_cld .lt. 3)then
       b_alpha_3d = clwc_3d * clwc2alpha * bksct_eff_clwc &
                  + cice_3d * cice2alpha * bksct_eff_cice &
                  + rain_3d * rain2alpha * bksct_eff_rain &
                  + snow_3d * snow2alpha * bksct_eff_snow 
     else
       b_alpha_3d = clwc_3d * clwc2alpha * bksct_eff_clwc &
                  + cice_3d * cice2alpha * bksct_eff_cice &
                  + rain_3d * rain2alpha * bksct_eff_rain &
                  + snow_3d * snow2alpha * bksct_eff_snow &
                  + aod_3d               * bksct_eff_aero
     endif

     write(6,*)' subroutine get_cloud_rad_faces2: solar alt/az ',solalt,solazi

     ntot = 0
     nsteps = 0
     refraction = 0.5 ! initialize

     heights_1d(:) = heights_3d(ni/2,nj/2,:)

     do k = 1,nk-1    
       htlow  = heights_1d(k)
       hthigh = heights_1d(k+1)
       do klut = nint(htlow),nint(hthigh)
         rklut = klut
         if(klut .ge. htlutlo .and. klut .le. htluthi)then
           htlut(klut) = float(k) + (rklut-htlow)/(hthigh-htlow)
         endif
       enddo ! klut
     enddo ! k

     dhmin = heights_1d(2)  - heights_1d(1)
     dhmax = heights_1d(nk) - heights_1d(1)
     write(6,*)' dhmin/dhmax = ',dhmin,dhmax

     faceperim = 0.4
     raysteps =     min(0.5 * grid_spacing_m,1500.)
     raysteps_low = min(grid_spacing_m,3000.)
!    raysteps_low = grid_spacing_m

     refr_mn = 0.0

     objalt_max = maxval(obj_alt)

!    Different criteria might be used at high altitude depending on how far
!    away the limb is, in turn related to horz_dep for the observer
     if(solalt + horz_dep_d_obs .ge. twi_alt)then ! daylight or early twilight

      do if = 1,6
       if(grid_spacing_m .lt. 10000.)then
         facestepij = 0.80 ! 0.85 causes missing points at 15.5 deg (top)
       else ! large domains
         if(if .eq. 1)then
           facestepij = 1.00 ! 0.85 causes missing points at 15.5 deg (top)
         else
           facestepij = 0.80 ! 0.85 causes missing points at 15.5 deg (top)
         endif
       endif

       write(6,*)
       is = 1 + (ni-1)*(i1(if)-1)
       ie = 1 + (ni-1)*(i2(if)-1)
       js = 1 + (nj-1)*(j1(if)-1)
       je = 1 + (nj-1)*(j2(if)-1)
       ks = 1 + (nk-1)*(k1(if)-1)
       ke = 1 + (nk-1)*(k2(if)-1)

       write(6,*)cfaces(if),' face'
       write(6,*)'rangei = ',is,ie  
       write(6,*)'rangej = ',js,je  

       nnew = 0
       nfacesteps = 0

       ris = is-faceperim; rie = ie+faceperim
       rjs = js-faceperim; rje = je+faceperim
       rks = ks          ; rke = ke

!      Exception needed if obs_alt = 90.
       objalt = obj_alt(nint(ris),nint(rjs))
       facesteph = min(max(grid_spacing_m*tand(objalt),dhmin),dhmax) ! * 0.3

       write(6,*)'rangek = ',ks,ke,facesteph

       do rit = ris,rie,facestepij
       do rjt = rjs,rje,facestepij

        if(obj_alt(nint(rit),nint(rjt)) .ge. twi_alt)then

!       do rkt = rks,rke,facestepk
        do htt = heights_1d(ks),heights_1d(ke),facesteph
         rkt = htlut(nint(htt))
         kmarch = nint(rkt)

         it = rit; jt = rjt; kt = rkt
         id = nint(rit); jd = nint(rjt); kd = nint(rkt)

!        kl = max(min(int(rkt),nk-1),1) ; kh = kl + 1
!        fk = rkt - kl
!        htt = heights_1d(kl) * (1.-fk) + heights_1d(kh) * fk

         ihit_terrain = 0
 
         slast = 0
         b_alpha_new = 0.
         btau = 0.

!        We presently do ray marching a constant distance intervals.
!        It may be more efficient to have successive steps march 
!        to the next 3D grid box boundary. We can then interpolate
!        from these "end points" to form the needed integrated values
!        at these grid box boundaries. We'd also want to know the
!        integrated values at the location where the ray comes closest
!        to the center of the grid box that is traverses.

!        When doing bilinear or trilinear interpolation, we're looking at
!        squares/cubes that have grid points lying on the vertices and the
!        'int' operation is used. When assigning the radiance values to the
!        grid, we are considering the nearest grid point using the 'nint'
!        operation. Here the cube is centered on a grid point.

!        Intersections with terrain are considered by finding maxima in the
!        bilinearly interpolated terrain field, relative to the ray height.
!        These maxima are thought to be located along lines connecting two
!        adjacent or diagonally adjacent grid points.

         raysteps_dynamic = min(200. / sind(max(objalt,3.)),raysteps_low)

!        Start ray trace at this point
         if(idebug .eq. 1)write(6,1)if,it,jt,kt
1        format(4i3)
         frac_abv_terrain = 1.0
         do ls = 0,4000 ! max number of ray segments

           if(ls .eq. 0)then ! values at start of trace                
             objalt = obj_alt(id,jd) + refr_mn
             objazi = obj_azi(id,jd)
             dids = -sind(objazi)*cosd(objalt)/grid_spacing_m
             dids = dids * grdasp(id,jd)
             djds = -cosd(objazi)*cosd(objalt)/grid_spacing_m
             dhtds = -sind(objalt)                                 
             dxyds =  cosd(objalt)
             curve = 0.
           else            ! ls > 0, thus after first looping
             slast = scurr ! initialized down below
             htlast = ht
           endif

!          if(idebug .eq. 1)write(6,*)'dids/djds/dhtds = ',dids,djds,dhtds

!          Calculate 'scurr' as the total path traversed so far by the ray.
!          The incremental path length ('ds') will be calculated later on.
!          altthr = 15. - (float(kmarch)/float(nk)) * 10. ! empirical for now
           altthr = 5.
           if(objalt .gt. altthr .AND. if .le. 2)then ! march by height levels
             nksteps = 1
             rkmarch = rkt - float(ls)/float(nksteps)
             if(rkmarch .le. 0.)goto 10 ! going outside the domain
             kmarch = nint(rkmarch) ! kl,kh
             htmarch = heights_1d(max(kmarch,1))
!            scurr = ((htt - htmarch) + curve) / (-dhtds)
             scurr = ((htmarch - htt) - curve) / (+dhtds)

             if(scurr .lt. slast)then
               write(6,*)' ERROR s<slast 1:',s,slast,ls,htt,htmarch,rkmarch
               stop
             endif

           elseif(objalt .lt. 3.0 .and. if .eq. 1)then
             scurr = float(ls) * raysteps_dynamic

           else ! nominally from 3 to 15 degrees
             scurr = float(ls) * raysteps_dynamic

             if(scurr .lt. slast)then
               write(6,*)' ERROR s<slast 2:',s,slast,ls,raysteps_dynamic
               stop
             endif

           endif

           ri = rit + dids*s
           rj = rjt + djds*s

           curve = (dxyds*s)**2 / (2.0*earth_radius)
           ht = htt + dhtds*s + curve

!          This can be used as part of a refraction strategy
!          refk = 0.179
!          ht = htt + dhtds*s + (dxyds*s)**2 / ((2.0/(1.-refk))*earth_radius)

!          if(ht .le. float(htluthi) .and. ht .ge. float(htlutlo))then
           if(ht .le.  1.0*(htluthi) .and. ht .ge.  1.0*(htlutlo))then
             rk = htlut(nint(ht))
!          elseif(ht .ge. htlutlo-500.)then
!            rk = 1.0
           else
             if(idebug .eq. 1)write(6,*)' ht outside lut',ht,heights_1d(1)
             goto 10
           endif

           idlast = id; jdlast = jd; kdlast = kd
           id = nint(ri); jd = nint(rj); kd = nint(rk)

           if(id .eq. idb .and. jd .eq. jdb)then
             idebug = 1
           else
             idebug = 0
           endif

           if(id .lt. 1 .or. id .gt. ni &
         .or. jd .lt. 1 .or. jd .gt. nj &
         .or. kd .lt. 1 .or. kd .gt. nk)then ! outside domain

             if(idebug .eq. 1)write(6,2)s,ri,rj,rk,ht             
2            format('s/ri/rj/rk/ht = ',5f9.2,' outside box')
             goto 10 
           else                              ! inside domain
             if(transm_3d(id,jd,kd) .ne. r_missing_data)then
               if(id .eq. idlast .and. jd .eq. jdlast .and. kd .eq. kdlast)then
                 if(idebug .eq. 1)write(6,3)s,ri,rj,rk,ht,id,jd,kd             
3                format('s/ri/rj/rk/ht = ',5f9.2,' same march',3i5)
               else
                 if(idebug .eq. 1)write(6,4)s,ri,rj,rk,ht,id,jd,kd             
4                format('s/ri/rj/rk/ht = ',5f9.2,' already assigned',3i5)
                 goto 9
               endif
             else ! transm_3d is missing
               if(idebug .eq. 1)write(6,5)s,ri,rj,rk,ht,htt,htmarch,curve,id,jd,kd             
5              format('s/ri/rj/rk/ht/cv = ',8f9.2,' new',3i5)
               nnew = nnew + 1
             endif

!            if(if .eq. 3 .AND. jt .eq. ni/2 .AND. kt .eq. nk/2)then
!              alt_theo2 = -atand((ht-htlast) / (s-slast))
!              write(6,7)id,jd,s,obj_alt(id,jd),alt_theo2
!7             format('key ray',2i5,f9.1,2f9.4)
!            endif

!            Valid trace (even if already assigned)
             illast = il; jllast = jl; kllast = kl

             il = max(min(int(ri),ni-1),1); fi = ri - il; ih=il+1
             jl = max(min(int(rj),nj-1),1); fj = rj - jl; jh=jl+1
             kl = max(min(int(rk),nk-1),1); fk = rk - kl; kh=kl+1

             if(il .eq. illast .and. jl .eq. jllast .and. kl .eq. kllast)then
               l_same_point = .true.
             else
               l_same_point = .false.
             endif

             if(ht - topo_a(id,jd) .le. 1000. .AND. frac_abv_terrain .gt. 0.)then

!              Interpolate to get topography at fractional grid point

!              if(il .le. 0)then
!                write(6,*)' il bounds check',il,id,ri,rj
!                stop
!              endif               

!              if(jl .le. 0)then
!                write(6,*)' jl bounds check',jl,jd,ri,rj
!                stop
!              endif               

               bi_coeff(1,1) = (1.-fi) * (1.-fj)
               bi_coeff(2,1) = fi      * (1.-fj)
               bi_coeff(1,2) = (1.-fi) *     fj 
               bi_coeff(2,2) = fi      *     fj 

               topo_bilin = sum(bi_coeff(:,:) * topo_a(il:ih,jl:jh))

               agl_last = agl
               agl = ht - topo_bilin

               if(ht .lt. topo_bilin)then
                 ihit_terrain = 1
                 if(agl .gt. 0.)then
                   frac_abv_terrain = 1.
                 elseif(agl_last .le. 0.)then
                   frac_abv_terrain = 0.
                 else ! this segment intersects terrain
                   frac_abv_terrain = agl_last / (agl_last - agl)
                 endif
               endif
             endif

             if(frac_abv_terrain .le. 0.)then
!              transm_3d(id,jd,kd) = 0.0                      
               transm_3d(id,jd,kd) = (1. - albedo) ! persistence
             else ! trace (partly) free of terrain
               ds = s - slast
               b_alpha_last = b_alpha_new

               tri_coeff(1,1,1) = (1.-fi) * (1.-fj) * (1.- fk)
               tri_coeff(2,1,1) = fi      * (1.-fj) * (1.- fk)
               tri_coeff(1,2,1) = (1.-fi) *     fj  * (1.- fk)
               tri_coeff(1,1,2) = (1.-fi) * (1.-fj) *      fk
               tri_coeff(1,2,2) = (1.-fi) *     fj  *      fk
               tri_coeff(2,1,2) = fi      * (1.-fj) *      fk
               tri_coeff(2,2,1) = fi      *     fj  * (1.- fk)
               tri_coeff(2,2,2) = fi      *     fj  *      fk

               if(l_same_point .eqv. .false.)then
                 b_alpha_cube(:,:,:) = b_alpha_3d(il:ih,jl:jh,kl:kh)
               endif

               b_alpha_new = sum(tri_coeff(:,:,:)*b_alpha_cube(:,:,:))
               b_alpha_new = max(b_alpha_new,0.) ! clean up far edge extrapolation
               if(ls .gt. 0)then
!                b_alpha_m = 0.5 * (b_alpha_3d(id,jd,kd) + b_alpha_3d(idlast,jdlast,kdlast))
                 b_alpha_m = 0.5 * (b_alpha_new + b_alpha_last)
                 dbtau = ds * b_alpha_m
                 btau = btau + dbtau
                 albedo = btau / (1. + btau)       
               else
                 b_alpha_last = b_alpha_new
                 albedo = 0.
               endif

               if(b_alpha_new .lt. 0.)then
                 write(6,*)' ERROR in b_alpha_new',b_alpha_new,ri,rj,rk,fi,fj,fk
                 write(6,*)'b_alpha_3d',b_alpha_3d(il:ih,jl:jh,kl:kh)
                 stop
               endif
               if(albedo .gt. 1.0 .or. albedo .lt. 0.0)then
                 write(6,*)' ERROR in albedo ',albedo,btau,dbtau,ds,b_alpha_m,b_alpha_new,b_alpha_last,s,slast
                 stop
               endif

               transm_3d(id,jd,kd) = (1. - albedo) ! * frac_abv_terrain              

               if(idebug .eq. 1)then
                   write(6,8)id,jd,kd,transm_3d(id,jd,kd),b_alpha_new,b_alpha_last
8                  format('ijk trn ba',3i5,f10.6,2e12.5)
               endif

               if(frac_abv_terrain .lt. 1.0)then ! 1st terrain intersection
                 transm_2d(id,jd) = (1. - albedo)
               endif
             endif

             nfacesteps = nfacesteps + 1

           endif

9          continue ! pick from from former rays
!          albedo = transm_3d(id,jd,kd)
!          btau = -(albedo) / (albedo-1.)

         enddo ! ls

10       continue

        enddo ! htt
        endif ! daytime / early twilight - trace column
       enddo ! j
       enddo ! i

       write(6,*)' Number new/steps assigned on this face =',nnew,nfacesteps

       ntot = ntot + nnew
       nsteps = nsteps + nfacesteps

       I4_elapsed = ishow_timer()

       write(6,*)' transm_3d column if = ',if,transm_3d(idb,jdb,:)

      enddo ! if

     else  ! solalt < twi_alt
      write(6,*)' solalt < twi_alt, raytrace not needed: ',solalt,horz_dep_d_obs,twi_alt
      I4_elapsed = ishow_timer()

     endif ! solalt

     npts = ni*nj*nk
     fractot = float(ntot)/float(npts)
     write(6,*)' Total is ',ntot, 'Potential pts is',npts,' frac',fractot
     write(6,*)' nsteps is ',nsteps,'frac is',float(nsteps)/float(npts)

     arg1 = minval(transm_3d)
     arg2 = maxval(transm_3d)
     if(arg1 .lt. 0. .OR. arg2 .gt. 1.0)then
       write(6,*)' WARNING: Range of transm_3d = ',arg1,arg2
!      stop
     else
       write(6,*)' range of transm_3d = ',arg1,arg2
     endif

!    where(transm_3d(:,:,:) .eq. r_missing_data)transm_3d(:,:,:) = 1.0

     I4_elapsed = ishow_timer()

     nshadow = 0
     if(solalt + horz_dep_d_obs .ge. twi_alt)then ! daylight or early twilight
       imiss = 0
       do k = 1,nk
!        patm_k = exp(-heights_1d(k)/8000.)
         patm_k = ztopsa(heights_1d(k)) / 1013.
         topo = redp_lvl ! generic topo value
         ht_agl = heights_1d(k) - topo
         patm_o3_msl = patm_o3(heights_1d(k))

!        See http://mintaka.sdsu.edu/GF/explain/atmos_refr/dip.html
         if(ht_agl .gt. 0.)then                               
           horz_dep_d = sqrt(2.0 * ht_agl / earth_radius) * 180./3.14
         else
           horz_dep_d = 0.
         endif

         obj_alt_last = r_missing_data

         tbuff(:,:) = transm_3d(:,:,k)

         do i = 1,ni; do j = 1,nj
           if(i .eq. idb .and. j .eq. jdb .and. k .eq. k)then
             iverbose = 1
           else
             iverbose = 0
           endif

           if(iverbose .eq. 1)then
             write(6,*)' Here iverbose 1 ',i,j,k
           endif

           if(transm_3d(i,j,k) .eq. r_missing_data)then
             imiss = imiss + 1
             if(imiss .le. 10)then
               write(6,11)i,j,k,obj_alt(i,j)
             endif
             if(i .eq. 10 * (i/10) .and. j .eq. 10 * (j/10))then
               write(6,11)i,j,k,obj_alt(i,j)
             endif
11           format(' missing at ',3i5,f9.3)

!            Assign/Interpolate
             if(i .gt. i .and. i .lt. ni .and. j .gt. 1 .and. j .lt. nj)then
               sumval = 0.
               cnt = 0.
               do ii = i-1,i+1   
               do jj = j-1,j+1
                 if(tbuff(ii,jj) .ne. r_missing_data)then
                   sumval = sumval + transm_3d(ii,jj,k)
                   cnt = cnt + 1.
                 endif
               enddo ! jj
               enddo ! ii
               if(cnt .gt. 0.)then
                 transm_3d(i,j,k) = sumval / cnt
               else
                 transm_3d(i,j,k) = 1.0
               endif
             else
               transm_3d(i,j,k) = 1.0
             endif ! inside horizontal domain
           endif ! missing value

           if(iverbose .eq. 1)then
             write(6,*)' Here iverbose 2 ',i,j,k,heights_1d(k),ht_agl,transm_3d(i,j,k)
           endif

           if(transm_3d(i,j,k) .eq. 0.)then
             nshadow = nshadow + 1
!          elseif(transm_3d(i,j,k) .lt. 0.)then
!            write(6,*)' ERROR: transm_3d < 0',i,j,k,transm_3d(i,j,k)
!            stop
!          elseif(transm_3d(i,j,k) .gt. 1.)then
!            write(6,*)' ERROR: transm_3d > 1',i,j,k,transm_3d(i,j,k)
!            stop
           else ! Calculate transm_4d
!            Direct illumination of the cloud is calculated here
!            Indirect illumination is factored in via 'scat_frac'
             obj_alt_thr = .01 ! abs(obj_alt(i,j)) * .00
             if(abs(obj_alt(i,j) - obj_alt_last) .gt. obj_alt_thr .OR. iverbose .eq. 1)then
!              ag = airmassf(cosd(90. - max(obj_alt(i,j),-3.0)),patm_k)
               ag = airmassf(90.-obj_alt(i,j), patm_k)

               if(.true.)then
                 aero_refht = redp_lvl
                 obj_alt_app = obj_alt(i,j) + refraction
                 arght = max(heights_3d(i,j,k),0.)
                 call get_airmass(obj_alt_app,arght              & ! I 
                                 ,patm_k,aero_refht,aero_scaleht & ! I
                                 ,earth_radius,iverbose &          ! I
                                 ,agdum,ao,aa,refr_deg)            ! O
               else
                 aa = 0.
               endif

               obj_alt_last = obj_alt(i,j)
               refraction = refr_deg 
             endif                                                     

             obj_alt_cld = obj_alt(i,j) + horz_dep_d + refraction

             if(iverbose .eq. 1)then
               write(6,*)' Here iverbose 3 ',i,j,k
             endif

!            Estimate solar extinction/reddening by Rayleigh scattering
!            at this cloud altitude
             if(obj_alt_cld .lt. -0.25)then ! (early) twilight cloud lighting
!              twi_int = .1 * 10.**(+obj_alt_cld * 0.4) ! magnitudes per deg
               twi_int = 0.
               rint = twi_int
               gint = twi_int
               bint = twi_int

             else ! object above horizon

               if(iverbose .eq. 1)then
                 write(6,*)' Here iverbose 4 ',i,j,k
               endif

               scat_frac = 1.00
               do ic = 1,nc
                 od_g = ag * ext_g(ic) * scat_frac

                 od_o = (o3_du/300.) * ext_o(ic) * ao * patm_o3_msl

                 ext_a(ic) = (wa(ic)/.55)**(-angstrom_exp_a)
                 od_a = aa * ext_a(ic) * aod

                 trans_c(ic) = trans(od_g + od_o + od_a)

                 if(iverbose .eq. 1)then
                   write(6,21)i,j,k,ic,heights_3d(i,j,k),obj_alt(i,j),obj_alt_cld
21                 format(' ijk/ic/ht/objalt/cld ',2i5,i4,i3,f9.0,2f9.2)
                   write(6,22)ag,agdum,ao,aa,od_g,od_o,od_a,trans_c(ic)
22                 format(' ag/agd/ao/aa/od_g/od_o/od_a/trans',8f9.4)
                 endif
               enddo

!              Fraction of solar disk (approximate)
               if(obj_alt_cld .gt. 0.25)then
                 sol_occ = 1.0
               else
                 occfrac = (obj_alt_cld - (-0.25)) / 0.5             
                 sol_occ = scurve(occfrac)             
               endif

               rint = trans_c(1) * sol_occ
               gint = trans_c(2) * sol_occ       
               bint = trans_c(3) * sol_occ              
             endif  

             transm_4d(i,j,k,1) = transm_3d(i,j,k) * rint
             transm_4d(i,j,k,2) = transm_3d(i,j,k) * gint
             transm_4d(i,j,k,3) = transm_3d(i,j,k) * bint

             if(iverbose .eq. 1)then
               write(6,23)sol_occ,rint,gint,bint,transm_4d(i,j,k,:)
23             format(' solocc/rint/gint/bint/transm_4d',f8.3,2x,3f10.5,2x,3f10.5)
             endif

           endif
         enddo ; enddo ! ij
       enddo ! k
     else ! nighttime: use red channel for sfc lighting
       do k = 1,nk       
         transm_4d(:,:,k,1) = sfc_glow(:,:) * 0.3 ! nominal backsct
       enddo ! k
       transm_3d(:,:,:) = 0.
     endif

     write(6,*)' nshadow = ',nshadow
     write(6,*)' imiss/fracmiss ',imiss,float(imiss)/float(npts)

     if(fractot .lt. 1.0)then
         write(6,*)' WARNING: missing points in get_cloud_rad_faces2',fractot
!        stop
     endif

!    Finalize
     do k = 1,nk
     do i = 1,ni
     do j = 1,nj
       if(topo_a(i,j) .gt. heights_3d(i,j,k))then
         transm_3d(i,j,k) = 0.
       endif
     enddo ! j
     enddo ! i
     enddo ! k

     write(6,*)' transm_3d column 2 = ',transm_3d(idb,jdb,:)

     I4_elapsed = ishow_timer()

     return
     end

!    Notes:
!       above threshold of 15 degrees has banding
!       at 2.2,4 degrees solalt works well - runs slow in top face
!       at 1330UTC -0.45 deg looks OK
!       at 1315UTC -3.19 deg good illumination - some banding
!                  turning off ag for transm_4d still has banding
