
        parameter (pi = 3.14159265)
        parameter (rpd = pi / 180.)

!       Refraction in degrees using apparent altitude    
        refractd_app(altd,patm) = &
               patm * (.0166666 / tand(altd + 7.31/(altd+4.4))) 

!       Airmasses relative to zenith at sea level pressure (true altitude)
        airmass_cosz(cosz) = &
          (          1.002432 * cosz**2 + 0.148386  * cosz + 0.0096467) &
                                       / &
          (cosz**3 + 0.149864 * cosz**2 + 0.0102963 * cosz + .000303978)
        airmassf(z,patm) = min(patm*airmass_cosz(cosd(min(z,93.))) &
                              ,40.*(1.0+sqrt(max(1.0-patm,0.)))) 

!       Minimum height of ray having negative altitude
        htminf(htmsl,alt,erad) = ((erad+htmsl) * cosd(alt)) - erad

        horz_depf(htmsl,erad) = acosd(erad/(erad+max(htmsl,0.)))

!       Ozone 
        parameter (o3_du = 300.)
        parameter (h_o3 = 22000.)
        parameter (d_o3 = 11000.)

!       Column integrated value per meter
        scalei_o3(htmsl) = exp(-((htmsl-h_o3)/d_o3)**2) / (sqrt(pi) * d_o3)

!       Airmasses estimated for thick shell (both below and above O3 max)
!       Effective scale height above h_o3 varies with 1./d_o3
        relscale_ht_o3(h) = d_o3 / (max(h+2.5*d_o3,h_o3+100.)-h_o3)
        hshellf(h) = max(h_o3-h,d_o3/2.5) * min(1.0,relscale_ht_o3(h))

!       hshellf(h) = max(h_o3-h,d_o3/2.5) * min(1.0,(h_o3/h)**2)
        airmasso(z,htmsl) = 1./SQRT(1.0-(SIND(min(z,90.))/(1.0+hshellf(htmsl)/6371e3))**2)

!       Error function from Wikipedia
        erf_pos(x) = 1.-1./(1.+.278393*x+.230389*x**2+.000972*x**3+.078108*x**4)**4
        erf(x) = sign(erf_pos(abs(x)),x)
!       gaussint(x) = 0.5 * sqrt(pi) * (1. - erf(abs(x)) * x/(abs(x)+1e-6))
        gaussint(x) = 0.5 * sqrt(pi) * (1. - erf(x))

!       patm_o3(ht) = max(min((32000. - ht) / 14000.,1.0),0.)
!       patm_o3(ht) = 1.0 - (gaussint((ht-25000.)/5000.) / 1.77245)
        patm_o3(ht) =        gaussint((ht-h_o3)/d_o3) / 1.77245 

        alpha_o3(ext_o3,htmsl) = (ext_o3 * exp(-((htmsl-h_o3)/d_o3)**2)) / (d_o3 * sqrt(pi))

!       Airmasses for thick shell - direct geometry
        parameter (h1_ha = 13000.)
        parameter (h2_ha = 25000.)
        parameter (h3_ha = 31000.)
        parameter (aod_ha = .015)
        parameter (alpha_ha = aod_ha / ((h2_ha-h1_ha)+0.5*(h3_ha-h2_ha)))
        patm_ah(h,h1,h2) = max(min(1.0,(h2-h)/(h2-h1)),0.)

        s_homo(z,h,re) = &
               re * (sqrt(cosd(z)**2 + 2.*h/re + (h/re)**2) - cosd(z) )
        am_homo_lyr(z,h1,h2,re) = &
               (s_homo(z,h2,re) - s_homo(z,h1,re)) / (h2 - h1)
        patm_a(ht,hr,hs) = exp(-(ht-hr)/hs) ! aerosols in column

!       Rayleigh scattering phase function including secondary scattering
!       This assumes 0.1 effective optical depth / secondary scattering
!                                   (this can be made variable later)
!                                        Primary                               Sec atm    
        rayleigh_pf(theta) = (0.90 * ( (1. + cosd(theta)**2) / (4./3.) ) ) + (0.10 * 1.0) 

        hg(g,pha) = (1.-g**2) / (1. + g**2 - 2.*g*cosd(pha))**1.5

        include 'wa.inc'

!       parameter (nc = 3)
        real ext_g(nc)               ! od per airmass
        data ext_g /.090,.144,.312/  ! 0.14 * (wa/.55)**(-4)              

        real ext_o(nc)               ! od per airmass
        data ext_o /.037,.029,.001/  ! interp from gflash.bas (300DU)

!       real wa(nc)                  ! wavelength (um)
!       data wa    /.615,.546,.450/

        real ext_a(nc)               ! aerosol extinction relative to .55um
!       data ext_a /0.90,1.00,1.10/  ! could change to (wa/.55)**(-1.3)
        real ext_ha(nc)              ! high aerosol ext. relative to .55um

        real angstrom_exp
!       data angstrom_exp /1.3/

!       Brightness conversions
!       https://irmafiles.nps.gov/reference/holding/476525
        s10_to_magsecsq(s10) = 10. - (2.5 * log10(s10/12.96e6))
        s10_to_nl(x) = x * 0.26189                     ! Chad Moore
  
        real magsecsq_to_s10
        magsecsq_to_s10(x) = 10.**((18.89075625-x)*0.4)
