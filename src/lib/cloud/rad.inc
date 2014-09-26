
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

!       Ozone airmasses adapted from Schaefer S&T code
        airmasso(z,htmsl) = 1./SQRT(1.0-(SIND(min(z,90.))/(1.0+(20e3-htmsl)/6378e3))**2)

!       Rayleigh scattering phase function including secondary scattering
!       This assumes 0.1 effective optical depth / secondary scattering
!                                   (this can be made variable later)
!                                        Primary                               Sec atm    
        rayleigh_pf(theta) = (0.90 * ( (1. + cosd(theta)**2) / (4./3.) ) ) + (0.10 * 1.0) 

!       Error function from Wikipedia
        erf(x) = 1.-1./(1.+.278393*x+.230389*x**2+.000972*x**3+.078108*x**4)**4
        gaussint(x) = 0.5 * sqrt(pi) * (1. - erf(abs(x)) * x/(abs(x)+1e-6))

!       patm_o3(ht) = max(min((32000. - ht) / 14000.,1.0),0.)
!       patm_o3(ht) = 1.0 - (gaussint((ht-25000.)/5000.) / 1.77245)
        patm_o3(ht) =        gaussint((ht-25000.)/5000.) / 1.77245 

        hg(g,pha) = (1.-g**2) / (1. + g**2 - 2.*g*cosd(pha))**1.5

        parameter (nc = 3)
        real ext_g(nc)                 ! od per airmass
        data ext_g /.07,.14,.28/       ! refine via Schaeffer

        real ext_o(nc)                 ! od per airmass
        data ext_o /.000,.016,.004/    ! scaled to 1/2 Schaeffer "OZ" values

        real wa(nc)                    ! wavelength (um)
        data wa    /0.68,0.55,0.44/

        real ext_a(nc)                 ! aerosol extinction relative to .55um
        data ext_a /1.00,1.00,1.00/    ! could change to (wa/.55)**(-1.3)
