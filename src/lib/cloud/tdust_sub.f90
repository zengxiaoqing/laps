
     subroutine tdust(mu0,tau_in,ssa_in,g_in,alb,iverb,t)

!    https://www.swsc-journal.org/articles/swsc/full_html/2015/01/swsc150027/swsc150027.html

     real k,kpterm,mu0

!    t is transmittance  (Output)
!    r is reflectance    (Output)

!    Variable transformation
     g = g_in / (1. + g_in)
     ssa = ssa_in * (1. - g_in**2) / (1. - ssa_in * g_in**2)
     tau = tau_in * (1. - ssa_in * g_in**2)

     if(iverb .eq. 1)write(6,*)' mu0    is ',mu0
     if(iverb .eq. 1)write(6,*)' ssa_in is ',ssa_in
     if(iverb .eq. 1)write(6,*)' ssa    is ',ssa
     if(iverb .eq. 1)write(6,*)' g_in   is ',g_in
     if(iverb .eq. 1)write(6,*)' g      is ',g
     if(iverb .eq. 1)write(6,*)' tau_in is ',tau_in
     if(iverb .eq. 1)write(6,*)' tau    is ',tau
       
     k = sqrt(3. * (1.-ssa) * (1. - g*ssa))       
     pprime = 2./3. * sqrt(3. * (1.-ssa) / (1. - g*ssa))
     alpha = 3./4. * mu0 * ssa * (1. +            g*(1.-ssa)) / (1. - mu0**2 * k**2)
     beta  = 1./2. * mu0 * ssa * (1./mu0 + 3.*mu0*g*(1.-ssa)) / (1. - mu0**2 * k**2)

     c3 = alb + (1.-alb) * alpha - (1.+alb) * beta
     c4 = 1. - alb + pprime * (1. + alb)
     c5 = 1. - alb - pprime * (1. + alb)

     c1num = (1.-pprime) * c3 * exp(-tau/mu0) - (alpha+beta) * c4 * exp( k*tau)
     c1den = (1.+pprime) * c4 * exp(k*tau)    - (1.-pprime)  * c5 * exp(-k*tau)
     c1 = -c1num / c1den

     c2num = (1.+pprime) * c3 * exp(-tau/mu0) - (alpha+beta) * c5 * exp(-k*tau)
     c2den = (1.+pprime) * c4 * exp(k*tau)    - (1.-pprime)  * c5 * exp(-k*tau)
     c2 = c2num / c2den

     t = c1 * exp(-k*tau) * (1. + pprime) + c2 * exp(k*tau) * (1.-pprime) - (alpha+beta-1.) * exp(-tau/mu0)

     if(iverb .eq. 1)write(6,*)' t is ',t

!    Reflectance Section
     b = 0.5 * (1. - g_in)         ! backscatter fraction
     btau = b * tau_in
     r_nonabs = btau / (1. + btau) ! reflectance (non-absorbing case, black surface)
       
     return
     end
