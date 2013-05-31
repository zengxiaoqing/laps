


       subroutine phase_func_moon(phase_angle_deg,mode,area_rel,area_es,phase_corr_moon)

       include 'trigd.inc'

       Implicit real*8(a-z)

!      Traditional equation
       arg = phase_angle_deg / 100.
       phase_corr_moon = 3.05 * arg - 1.02 * arg**2 + 1.05 * arg**3

!      arg = phase_angle_deg ! Schaefer 1991
!      phase_corr = .026 * arg + 4e-9 * arg**4

!      Opposition effect
       if(phase_angle_deg .lt. 7.0)then
           arg = 0.18 / 7.
           phase_corr_moon = phase_corr_moon + (arg * (phase_angle_deg - 7.))
       endif

!      Supplementary information
       area_rel_7 = (cosd(173.) + 1.) / 2.

       area_rel = (cosd(phase_angle_deg) + 1.) / 2.

!      New formulation matches at 152 degrees / 28 degrees

!      Correct for shadowing effect
       area_rel = area_rel - (area_rel_7 * (1. - area_rel))

       area_rel = max(area_rel,0.d0)

!      Correct Earthshine using Earth's phase function
       arg = (180. - phase_angle_deg) / 100.
       phase_corr_earth = 1.30 * arg - 0.19 * arg**2 + 0.48 * arg**3
       frac_bri_earth = 10**(-phase_corr_earth * 0.4)
       area_es = 9.1201e-5 * frac_bri_earth

       surface_brightness = .073

       frac_full = area_rel * surface_brightness + area_es

       if(mode .eq. 1 .AND. phase_angle_deg .gt. 152.)then ! Modify for Crescent Moon
 
           phase_corr_moon = -log10(frac_full) * 2.5
       
       endif

       return
       end
