
        subroutine get_airmass(alt,htmsl,patm &          ! I (add htagl?)
                              ,aero_refht,aero_scaleht & ! I
                              ,earth_radius,iverbose &   ! I
                              ,ag,ao,aa)                 ! O

!       Airmasses relative to zenith at sea level pressure (ag - gas)
!       Airmasses relative to zenith at sea level pressure (ao - ozone)
!       Airmasses relative to zenith at aero refht         (aa - aerosol)

        use mem_namelist, ONLY: r_missing_data ! , earth_radius

        include 'trigd.inc'
        include 'rad.inc'

!       Traditional aerosol empirical forumla (valid at H = 1900m)
!       Try using new relationship that includes scale height
        am_aero(z)=1./(COS(Z)+.0123*EXP(-24.5*COS(Z))) ! radians apparent
        am_homo_wiki(z,h,r) = r/h * sqrt(cos(z)**2 + 2.*h/r + (h/r)**2) - (r/h) * cos(z)

        zapp = 90. - alt
        zappi = 90. + alt

        ztrue  = zapp  + refractd_app(alt ,patm)
        ztruei = zappi + refractd_app(abs(alt) ,patm)

        if(iverbose .eq. 1)then
          write(6,11)alt,htmsl,patm,earth_radius
11        format(/' start get_airmass: alt/htmsl/patm/erad' &
                 ,f8.2,f11.1,f9.4,f10.0)
        endif

        call get_htmin(alt,patm,htmsl,earth_radius,iverbose & ! I
                      ,patm2,htmin)                           ! O

!       Gas component for Rayleigh Scattering
        if(alt .lt. -3.)then ! high looking down
          if(htmin .lt. 0.0)then ! hit ground
            patm_gnd = max(ztopsa(aero_refht)/1013.25,patm) ! assumed ground
            ag = airmassf(zappi,(patm_gnd-patm))
            if(iverbose .eq. 1)write(6,21)zappi,patm_gnd,patm,ag
21          format('  high gas   hitting ground',4f9.4)
          else ! grazes atmosphere
            ztrue0 = 90. + refractd_app(0.,patm2)
!           alttrue0 = -refractd_app(0.,patm2)
            ag = (2.*airmassf(ztrue0,patm2)) - airmassf(ztruei,patm)
            if(iverbose .eq. 1)then
              write(6,*)' high gas   looking down at htmin',htmin,ag
            endif
          endif
        else ! standard situation
          ag = airmassf(ztrue,patm)
          if(iverbose .eq. 1)then
            write(6,*)' standard gas   situation',htmsl,alt,ag
          endif
        endif

!       Ozone component
        if(alt .lt. 0.)then ! high looking down
          if(htmin .lt. 0.0)then ! hit ground
            patm_refht = patm_o3(min(aero_refht,htmsl))
            ao = airmasso(zappi,aero_refht)  &
               * (patm_o3(aero_refht)-patm_o3(htmsl))
            ao = max(ao,0.)
            if(iverbose.eq.1)write(6,31)zappi,patm_refht,patm_o3(htmsl),ao
31          format('  high ozone hitting ground',4f9.4)
          else ! grazes atmosphere
            patmo  = patm_o3(htmsl)
            patm2o = patm_o3(htmin)
            ztrue0 = 90. + refractd_app(0.,patm2o)
!           alttrue0 = -refractd_app(0.,patm2o)
!           ao1 = airmasso(ztruei,htmin) * patm2o ! ztruei from htmin
            ao2 = airmasso(ztruei,htmsl) * patmo  ! upward from observer
            ao3 = airmasso(ztrue0,htmin) * patm2o ! 0 deg from htmin
            ao = (ao3-ao2) + ao3
            if(iverbose .eq. 1)then
              write(6,32)htmin,ao2,ao3,ao
32            format('  high ozone looking down at htmin',f12.2,3f9.3)
            endif
          endif
        else
          ao = airmasso(zapp,htmsl) * patm_o3(htmsl)
          if(iverbose .eq. 1)then
            write(6,33)htmsl,hshellf(htmsl),alt,ao
33          format('  standard ozone situation',f12.2,f9.2,f9.2,f9.3)
          endif
        endif

!       Aerosol component

!       Note that near the horizon the aerosol airmass should
!       exceed the gas component by the inverse square root of the
!       respective scale heights.

        patm_aero = exp(-((htmsl-aero_refht) / aero_scaleht))
        if(alt .lt. 0.)then ! high looking down
          if(htmin .lt. 0.0)then ! hit ground
            ZZ = zappi * rpd
!           aa=am_aero(ZZ)                                * (1.0-patm_aero)
            aa=am_homo_wiki(ZZ,aero_scaleht,earth_radius) * (1.0-patm_aero)
            aa = max(aa,0.)
            if(iverbose .eq. 1)then
              write(6,41)zappi,htmin,patm_aero,aa    
41            format('  high aero  hitting ground',4f9.4)           
            endif
          else ! grazes atmosphere
            ZZ = zappi * rpd
            patm2_aero = exp(-((htmin-aero_refht) / aero_scaleht))
!           aa1=am_aero(ZZ)                                * patm2_aero
!           aa2=am_aero(ZZ)                                * patm_aero
            aa2=am_homo_wiki(ZZ,aero_scaleht,earth_radius) * patm_aero

            ZZ = 90. * rpd
!           aa3=am_aero(ZZ)                                * patm2_aero
            aa3=am_homo_wiki(ZZ,aero_scaleht,earth_radius) * patm2_aero

            aa = (aa3-aa2) + aa3
            if(iverbose .eq. 1)then
              write(6,*)' high aero  looking down at htmin',htmin,aa 
            endif
          endif
        else
          ZZ = (min(zapp,90.)) * rpd
!         aa=am_aero(ZZ)                                * patm_aero
          aa=am_homo_wiki(ZZ,aero_scaleht,earth_radius) * patm_aero
          if(iverbose .eq. 1)then
            write(6,*)' standard aero  situation',htmsl,alt,aa       
          endif
        endif

        if(iverbose .eq. 1)write(6,*)' returning from get_airmass'

        return
        end

        subroutine get_htmin(alt,patm,htmsl,earth_radius,iverbose & ! I
                            ,patm2,htmin)                           ! O

!       Find minimum MSL height of a light ray compared to Earth sfc
!       Refraction is taken into account

        include 'trigd.inc'

        parameter (pi=3.14159265)
        parameter (rpd=pi/180.)

        if(alt .lt. 0)then
          patm2 = patm
          niter = 3; iter = 1
          do while (iter .le. niter)
!           rk and htmin could use more accuracy when htmsl is high
            wtm = 0.5 + 0.5 * (1.-patm/patm2)
            wtm2 = (1.-wtm)
            rk = (patm*wtm+patm2*wtm2) * 0.13 ! curvature at midpoint of ray
            erad_eff = earth_radius / (1. - rk)
            htmin_s = htmsl - erad_eff * (alt*rpd)**2 / 2.
            htmin_l = ((erad_eff+htmsl) * cosd(alt)) - erad_eff    
            if(alt .gt. -0.25)then ! small angle forumla
              htmin = htmin_s
            else                 ! large angle formula
              htmin = htmin_l  
            endif
            htmin = max(htmin,-500.)
            patm2 = ztopsa(htmin) / 1013.25
            if(iverbose .eq. 1)then
              write(6,1)iter,htmsl,alt,cosd(alt),htmin_s,htmin_l,htmin,patm,patm2,rk
            endif
1           format('    iter/htmsl/alt/cosd/htmin/patm/patm2/rk = ',i4,f10.1,2f9.3,2f11.0,f11.1,2f10.7,f9.5)
            iter = iter + 1
          enddo ! while
        else
          htmin = htmsl
        endif

        return
        end
