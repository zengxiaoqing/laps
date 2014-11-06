     
     subroutine get_cloud_rad_faces(             &
                obj_alt,obj_azi,                 & ! I
                solalt,solazi,                   & ! I 
                clwc_3d,cice_3d,rain_3d,snow_3d, & ! I
                topo_a,                          & ! I
                ni,nj,nk,idb,jdb,                & ! I
                heights_3d,                      & ! I 
                sfc_glow,                        & ! I
                transm_3d,transm_4d)               ! O

     include 'rad.inc'

     trans(od) = exp(-min(od,80.))

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
     real transm_3d(ni,nj,nk)
     real transm_4d(ni,nj,nk,nc)

     real obj_alt(ni,nj),obj_azi(ni,nj)
     real topo_a(ni,nj)
     real sfc_glow(ni,nj)

     real heights_1d(nk)

     real trans_c(nc)
     real bi_coeff(2,2)

     parameter (htlutlo = -1000)
     parameter (htluthi = 30000)
     real htlut(htlutlo:htluthi)

!    Backscattering efficiencies
     real, parameter :: bksct_eff_clwc    = .063
     real, parameter :: bksct_eff_cice    = .14
     real, parameter :: bksct_eff_rain    = .063
     real, parameter :: bksct_eff_snow    = .14
     real, parameter :: bksct_eff_graupel = .30

!    Scattering efficiencies
     real, parameter :: q_clwc    = 2.0
     real, parameter :: q_cice    = 2.0
     real, parameter :: q_rain    = 1.0
     real, parameter :: q_snow    = 1.0
     real, parameter :: q_graupel = 1.0

!    Densities
     real, parameter :: rholiq     =   1e3 ! kilograms per cubic meter
     real, parameter :: rhosnow    = .07e3 ! kilograms per cubic meter
     real, parameter :: rhograupel = .50e3 ! kilograms per cubic meter

!    Effective radii
     real, parameter :: reff_clwc    = .000020 ! m
     real, parameter :: reff_cice    = .000040 ! m
     real, parameter :: reff_rain    = .000750 ! m
     real, parameter :: reff_snow    = .004000 ! m
     real, parameter :: reff_graupel = .010000 ! m


     clwc2alpha = 1.5 / (rholiq  * reff_clwc)
     cice2alpha = 1.5 / (rholiq  * reff_cice)
     rain2alpha = 1.5 / (rholiq  * reff_rain)
     snow2alpha = 1.5 / (rhosnow * reff_snow)
     pice2alpha = 1.5 / (rhograupel * reff_graupel)

     idebug = 0

     r_missing_data = 1e37
     earth_radius = 6371000.

     grid_spacing_m = 500.

     transm_3d = r_missing_data
     transm_4d = 0.                   

     b_alpha_3d = clwc_3d * clwc2alpha * bksct_eff_clwc &
                + cice_3d * cice2alpha * bksct_eff_cice &
                + rain_3d * rain2alpha * bksct_eff_rain &
                + snow_3d * snow2alpha * bksct_eff_snow 

     write(6,*)' subroutine get_cloud_rad_faces: solar alt/az ',solalt,solazi

     ntot = 0
     nsteps = 0

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
     facestepij = 0.80 ! 0.85 causes missing points at 15.5 deg (top face)
     raysteps = 0.5 * grid_spacing_m

     do if = 1,6
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

!      do rkt = rks,rke,facestepk
       do htt = heights_1d(ks),heights_1d(ke),facesteph
         rkt = htlut(nint(htt))

         it = rit; jt = rjt; kt = rkt
         id = nint(rit); jd = nint(rjt); kd = nint(rkt)

!        kl = max(min(int(rkt),nk-1),1) ; kh = kl + 1
!        fk = rkt - kl
!        htt = heights_1d(kl) * (1.-fk) + heights_1d(kh) * fk

         ihit_terrain = 0
 
         slast = 0
         btau = 0.

!        Start ray trace at this point
         if(idebug .eq. 1)write(6,1)if,it,jt,kt
1        format(4i3)
         do ls = 0,4000

           if(ls .eq. 0)then ! values at start of trace                
             objalt = obj_alt(id,jd)
             objazi = obj_azi(id,jd)
             dids = -sind(objazi)*cosd(objalt)/grid_spacing_m
             djds =  cosd(objazi)*cosd(objalt)/grid_spacing_m
             dhtds = -sind(objalt)                                 
             dxyds =  cosd(objalt)
           else
             slast = s
           endif

!          if(idebug .eq. 1)write(6,*)'dids/djds/dhtds = ',dids,djds,dhtds

           if(objalt .gt. 15. .AND. if .le. 2)then ! march by height levels
             nksteps = 2
             rkmarch = rkt - float(ls)/float(nksteps)
             if(rkmarch .le. 0.)goto10 ! going outside the domain
             kmarch = nint(rkmarch) ! kl,kh
             htmarch = heights_1d(max(kmarch,1))
             s = (htt - htmarch) / (-dhtds)

             if(s .lt. slast)then
               write(6,*)' ERROR s<slast 1:',s,slast,ls,htt,htmarch,rkmarch
               stop
             endif

           else
             s = float(ls) * raysteps

             if(s .lt. slast)then
               write(6,*)' ERROR s<slast 2:',s,slast,ls,raysteps
               stop
             endif

           endif

           ri = rit + dids*s
           rj = rjt + djds*s

           ht = htt + dhtds*s + (dxyds*s)**2 / (2.66666*earth_radius)
           if(ht .le. float(htluthi) .and. ht .ge. float(htlutlo))then
             rk = htlut(nint(ht))
!          elseif(ht .ge. htlutlo-500.)then
!            rk = 1.0
           else
             if(idebug .eq. 1)write(6,*)' ht outside lut',ht,heights_1d(1)
             goto10
           endif

           idlast = id; jdlast = jd; kdlast = kd
           id = nint(ri); jd = nint(rj); kd = nint(rk)

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
3                format('s/ri/rj/rk/ht = ',5f9.2,' same march',3i4)
               else
                 if(idebug .eq. 1)write(6,4)s,ri,rj,rk,ht,id,jd,kd             
4                format('s/ri/rj/rk/ht = ',5f9.2,' already assigned',3i4)
                 goto 9
               endif
             else ! transm_3d is missing
               if(idebug .eq. 1)write(6,5)s,ri,rj,rk,ht,id,jd,kd             
5              format('s/ri/rj/rk/ht = ',5f9.2,' new',3i4)
               nnew = nnew + 1
             endif

!            Valid trace (even if already assigned)

             if(ht - topo_a(id,jd) .le. 1000. .AND. ihit_terrain .ne. 1)then

!              Interpolate to get topography at fractional grid point
               il = max(min(int(ri),ni-1),1); fi = ri - il; ih=il+1
               jl = max(min(int(rj),nj-1),1); fj = rj - jl; jh=jl+1

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

               if(ht .lt. topo_bilin)then
                 ihit_terrain = 1
               endif
             endif

             if(ihit_terrain .eq. 1)then
               transm_3d(id,jd,kd) = 0.0                      
!              transm_4d(id,jd,kd,:) = 0.0                      
             else ! trace free of terrain
               ds = s - slast
               db_alpha_m = 0.5 * (b_alpha_3d(id,jd,kd) + b_alpha_3d(idlast,jdlast,kdlast))
               dbtau = ds * db_alpha_m
               btau = btau + dbtau
               albedo = btau / (1. + btau)       
               if(albedo .gt. 1.0 .or. albedo .lt. 0.0)then
                 write(6,*)' ERROR in albedo ',albedo,btau,dbtau,ds,db_alpha_m,s,slast
                 stop
               endif
               transm_3d(id,jd,kd) = 1. - albedo              
             endif

             nfacesteps = nfacesteps + 1

           endif

9          continue

         enddo ! ls

10       continue

       enddo ! k
       enddo ! j
       enddo ! i

       write(6,*)' Number new/steps assigned on this face =',nnew,nfacesteps

       ntot = ntot + nnew
       nsteps = nsteps + nfacesteps

       I4_elapsed = ishow_timer()

     enddo ! if

     npts = ni*nj*nk
     fractot = float(ntot)/float(npts)
     write(6,*)' Total is ',ntot, 'Potential pts is',npts,' frac',fractot
     write(6,*)' nsteps is ',nsteps,'frac is',float(nsteps)/float(npts)

     arg1 = minval(transm_3d)
     arg2 = maxval(transm_3d)
     if(arg1 .lt. 0. .OR. arg2 .gt. 1.0)then
       write(6,*)' ERROR: Range of transm_3d = ',arg1,arg2
       stop
     else
       write(6,*)' range of transm_3d = ',arg1,arg2
     endif

     nshadow = 0
     if(solalt .ge. -4.)then ! daylight or early twilight
       do k = 1,nk
         patm_k = exp(-heights_1d(k)/8000.)
!        patm_k = ztopsa(heights_1d(k)) / 1013.
         topo = 1500. ! generic topo value (possibly redp_lvl?)
         ht_agl = heights_1d(k) - topo

!        See http://mintaka.sdsu.edu/GF/explain/atmos_refr/dip.html
         if(ht_agl .gt. 0.)then                               
           horz_dep_d = sqrt(2.0 * ht_agl / earth_radius) * 180./3.14
         else
           horz_dep_d = 0.
         endif

         obj_alt_last = r_missing_data

         do i = 1,ni; do j = 1,nj
           if(transm_3d(i,j,k) .eq. r_missing_data)then
             if(i .eq. 10 * (i/10) .and. j .eq. 10 * (j/10))then
               write(6,*)' missing at ',i,j,k
             endif
           elseif(transm_3d(i,j,k) .eq. 0.)then
             nshadow = nshadow + 1
!          elseif(transm_3d(i,j,k) .lt. 0.)then
!            write(6,*)' ERROR: transm_3d < 0',i,j,k,transm_3d(i,j,k)
!            stop
!          elseif(transm_3d(i,j,k) .gt. 1.)then
!            write(6,*)' ERROR: transm_3d > 1',i,j,k,transm_3d(i,j,k)
!            stop
           else
!            Calculate transm_4d
             refraction = 0.5 ! typical value near horizon

             obj_alt_cld = obj_alt(i,j) + horz_dep_d + refraction

!            Estimate solar extinction/reddening by Rayleigh scattering
!            at this cloud altitude
             if(obj_alt_cld .lt. 0.)then ! (early) twilight cloud lighting
               twi_int = .1 * 10.**(+obj_alt_cld * 0.4) ! magnitudes per deg
               rint = twi_int
               gint = twi_int
               bint = twi_int

             elseif(obj_alt_cld .ge. 0.)then   ! low daylight sun
!              Direct illumination of the cloud is calculated here
!              Indirect illumination is factored in via 'scat_frac'
               obj_alt_thr = abs(obj_alt(i,j)) * .00
               if(abs(obj_alt(i,j) - obj_alt_last) .gt. obj_alt_thr)then
!                am = airmassf(cosd(90. - max(obj_alt(i,j),-3.0)),patm_k)
                 am = airmassf(90.-obj_alt(i,j), patm_k)
                 obj_alt_last = obj_alt(i,j)
               endif                                                     
               scat_frac = 0.75
               do ic = 1,nc
                 trans_c(ic) = trans(am*ext_g(ic)*scat_frac)
               enddo
               rint = trans_c(1)
               gint = trans_c(2)                 
               bint = trans_c(3)               
             endif  

             transm_4d(i,j,k,1) = transm_3d(i,j,k) * rint
             transm_4d(i,j,k,2) = transm_3d(i,j,k) * gint
             transm_4d(i,j,k,3) = transm_3d(i,j,k) * bint
           endif
         enddo ; enddo ! ij
       enddo ! k
     else ! use red channel for sfc lighting
       do k = 1,nk       
         transm_4d(:,:,k,1) = sfc_glow(:,:) ! /day_int 
       enddo ! k
     endif

     write(6,*)' nshadow = ',nshadow

     if(fractot .lt. 1.0)then
         write(6,*)' ERROR: missing points in get_cloud_rad_faces',fractot
         stop
     endif

     I4_elapsed = ishow_timer()

     return
     end

!    Notes:
!       above threshold of 15 degrees has banding
!       at 2.2,4 degrees solalt works well - runs slow in top face
!       at 1330UTC -0.45 deg looks OK
!       at 1315UTC -3.19 deg good illumination - some banding
!                  turning off am for transm_4d still has banding
