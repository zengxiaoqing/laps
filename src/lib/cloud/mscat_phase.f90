
        subroutine mscat_phase(tau,ssa,g,iverbose,gmean,ssa_eff)

!       Scattering event:
!       1 - trans_tau is scattered at least once          
!       'traus_tau' is in direct beam (scattered zero times)
!       (1 - trans_tau) * trans_tau is scattered only once
!       (1 - trans_tau) * (1 - trans_tau) is scattered at least twice
!       (1 - trans_tau) * (1 - trans_tau) * (trans_tau) is scattered only twice
!       (1 - trans_tau) * (1 - trans_tau) * (1. - trans_tau) is scattered at least thrice

        parameter (n_order=50)

        real frac_order(0:n_order)
        real rel_order(1:n_order)
        real p_order(0:n_order)
        real p_scat_cum(0:n_order)

        trans(od) = exp(-min(od,80.))
        opac(od) = 1.0 - trans(od)

        trans_tau = trans(tau)

        if(iverbose .ge. 1)then
           write(6,*)' tau,ssa,trans(tau) = ',tau,ssa,trans(tau)
           write(6,*)
           write(6,*)' i-1 p_scat_cum p_ord    g(i)'
        endif

        remainder = 1.0

        p_order(0) = 0.
        gsum = 0.
        osum = 0.
        do i = 0,n_order
           p_scat_cum(i) = (1. - trans_tau)**i * ssa**i
           if(i .lt. n_order)then
              p_order(i) = p_scat_cum(i) * trans_tau
           else ! last one is defined as cumulative
              p_order(i) = p_scat_cum(i)
           endif
           if(i .ge. 1)then
              ginc = g**float(i)
              gsum = gsum + g**float(i) * p_order(i)
              osum = osum +               p_order(i)
           else
              ginc = -999.
           endif
           if(iverbose .ge. 1)then
              write(6,11)i,p_scat_cum(i),p_order(i),ginc
           endif
        enddo

11      format(i3,4f9.5)

        ssa_eff = sum(p_order(1:n_order))
        gmean = gsum / osum

        if(iverbose .ge. 1)then
           write(6,*)' sum of p_order/ssa_eff is ',ssa_eff
           write(6,*)' gmean = ', gmean
           write(6,*)
        endif
        
        return
        end

          
