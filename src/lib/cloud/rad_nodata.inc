
        parameter (pi = 3.14159265)
        parameter (rpd = pi / 180.)

        real nl_to_s10

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
!       parameter (o3_du = 300.)
!       parameter (h_o3 = 22000.)
!       parameter (d_o3 = 11000.)

!       Column integrated value per meter
        scalei_o3(htmsl) = exp(-((htmsl-h_o3)/d_o3)**2) / (sqrt(pi) * d_o3)

!       General Airmass routines

!       http://aas.aanda.org/articles/aas/pdf/1998/01/ds1449.pdf (Section 6)
        am_thsh(z,h,r) = 1./ sqrt( 1.-( (r/(r+h))**2 * (sind(z))**2) )

!       http://www.wolframalpha.com/input/?i=integrate+1%2F%28sqrt%281-%28r%2F%28r%2Bh%29%29%5E2+sin%5E2z%29%29+dh+
        am_thsh_indefint(z,h,r) = (2.*h**2 + 4.*h*r + r**2 * cosd(2.*z) + r**2) / &
            (2.*(h+r)*sqrt(1.-r**2*sind(z)**2/(h+r)**2))

        am_thsh_defint(z,h1,h2,r) = am_thsh_indefint(z,h2,r) - am_thsh_indefint(z,h1,r) 

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

!       High altitude aerosols, uniform shell with ramp above
        parameter (h1_ha = 13000.)
        parameter (h2_ha = 25000.)
        parameter (h3_ha = 31000.)
!       parameter (aod_ha = .015)
!       parameter (aod_ha = .004)
!       parameter (alpha_ha = aod_ha / ((h2_ha-h1_ha)+0.5*(h3_ha-h2_ha)))
        patm_ah(h,h1,h2) = max(min(1.0,(h2-h)/(h2-h1)),0.)

!       Airmasses for thick homogenous elevated shell
!       's_homo' is the slant path of a homogenous layer with observer at bottom edge
!       'am_homo_lyr' is airmasses of partial elevated layer (base above observer)
        s_homo(z,h,re) = &
               re * (sqrt(cosd(z)**2 + 2.*h/re + (h/re)**2) - cosd(z) )
        am_homo_lyr(z,h1,h2,re) = &
               (s_homo(z,h2,re) - s_homo(z,h1,re)) / (h2 - h1)

        patm_a(ht,hr,hs) = exp(-(ht-hr)/hs) ! aerosols in column

!       Rayleigh scattering phase function including secondary scattering
!       This assumes 0.1 effective optical depth / secondary scattering
!                                   (this can be made variable later)
!                                            Primary                               Sec atm    
        rayleigh_pf(theta)   =           (1.00 * ( (1. + cosd(theta)**2) / (4./3.) ) ) + (0.00 * 1.0) 
        rayleigh_pf_s(theta) =           (1.00 * ( (1. + cosd(theta)**2) / (4./3.) ) ) + (0.00 * 1.0) 

!       Multiple scattering phase function relates to maximum polarization
        rayl_mxpol(taur) = exp(-min(taur,80.)) ! equivalent to trans()
        rayleigh_pf_m(theta,taur) = rayleigh_pf_s(theta) * rayl_mxpol(taur) + 1.0 * (1.0 - rayl_mxpol(taur))

        hg(g,pha) = (1.-g**2) / (1. + g**2 - 2.*g*cosd(pha))**1.5

!       parameter (nc = 3)
        real ext_g(nc)                  ! od per airmass

!       Multiply Earth values by .006/.038
!       data ext_g /.00142,.00227,.00493/                                       

!       real ext_o(nc)               ! od per airmass
!       data ext_o /.037,.029,.001/  ! interp from gflash.bas (300DU)

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
        nl_to_s10(x) = x / 0.26189                     ! Chad Moore

        real magsecsq_to_s10
        magsecsq_to_s10(x) = 10.**((18.89075625-x)*0.4)

        real nwcm2sr_to_wm2sr
        nwcm2sr_to_wm2sr(x) = x * 1e-5

        parameter (efficiency_lum_sun = .136)
        wm2sr_to_nl_550(x) = 2.113e8 * x  
        wm2sr_to_nl_sun(x) = wm2sr_to_nl_550(x) * efficiency_lum_sun

!       Anisotropic exitance factor for sfc lights (angle above ground looking from lights)
        aef_f(a) = 1. / sind(max(a,+19.47))

!       Microphysics
        reff_cice_f(cice) = .000034 * (max(cice,1e-6)/1.66e-6)**(+0.14)
        reff_clwc_f(clwc) = (3. + 8. * (max(clwc,1e-6)*1.00e3)**(+0.50)) / 1e6

	cice_reff_f(reff) = 1.374408e26 * reff**(50./7.)
	clwc_reff_f(reff) = 1.5625e7 * reff**2 - 93.75 * reff + .000140625

        albedo_to_btau(albedo) = -(albedo) / (albedo-1.)
	
!       Aerosol scattering
        dhg_gterm(pha,g) = 1. / (1. + g**2 - 2.*g*cosd(pha))**1.5
        dhg_fterm(pha,g) = (3. * ((cosd(pha))**2) - 1.) / (2. * (1.+g**2)**1.5)
        dhg(pha,g,f) = (1.-g**2) * (dhg_gterm(pha,g) + f * dhg_fterm(pha,g))
        dhg2(pha,f,p) = (1.-p) * dhg(pha,g1,f) + p * dhg(pha,g2,f) + 0.28/scatter_order * cosd(pha)
