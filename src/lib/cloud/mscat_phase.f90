
        subroutine mscat_phase(tau,ssa,g,iverbose &                     ! I
                              ,rad &                                    ! I
                              ,p_order,rel_order,n_order,gmean,ssa_ret) ! O

!       Estimates fraction of received photons of various scattering orders

!       Scattering event:
!       1 - trans_tau is scattered at least once          
!       'traus_tau' is in direct beam (scattered zero times)
!       (1 - trans_tau) * trans_tau is scattered only once
!       (1 - trans_tau) * (1 - trans_tau) is scattered at least twice
!       (1 - trans_tau) * (1 - trans_tau) * (trans_tau) is scattered only twice
!       (1 - trans_tau) * (1 - trans_tau) * (1. - trans_tau) is scattered at least thrice

!       Presently works for case with low "rad" having light paths traveling through the cloud
!       Can be extended to work for high "rad" with light scattering off the cloud top or edge          

        real frac_order(0:n_order)
        real rel_order(1:n_order)  ! P of received photon with this order
        real p_order(0:n_order)    ! P of emitted photon with this order
        real p_scat_cum(0:n_order) 

        trans(od) = exp(-min(od,80.))
        opac(od) = 1.0 - trans(od)

        trans_tau = trans(tau)

        if(iverbose .ge. 1)then
           write(6,*)' tau,ssa,g,trans(tau) = ',tau,ssa,g,trans(tau)
           write(6,*)'  i  p_scat_cum p_order rel_order'
        endif

        remainder = 1.0

        p_order(0) = 0.
        gsum = 0.
        osum = 0.
        if(rad .lt. 0.8)then ! low rad, mostly transmission
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
           enddo
           ssa_ret = ssa**(tau**2)
        else ! high rad, mostly reflection
           p_order(0) = 0.
           p_order(n_order) = 0.30 * ssa**n_order
           p_order(1) = 1.0 - p_order(n_order)
        endif

11      format(i4,2x,6f9.5)

        ssa_eff = sum(p_order(1:n_order))
        rel_order(1:n_order) = p_order(1:n_order) / ssa_eff
        gmean = gsum / osum

        if(iverbose .ge. 1)then
           do i = 1,n_order
!               if(i .le. 10 .or. i .eq. (i/10)*10)then
                  write(6,11)i,p_scat_cum(i),p_order(i),rel_order(i) ! ,rel_order(i)*ssa_eff
!               endif
           enddo ! i
        endif

        if(iverbose .ge. 1)then
           write(6,*)' sum of p_order/ssa_eff is ',ssa_eff
           write(6,*)' gsum/osum/gmean = ', gsum, osum, gmean
           write(6,*)
        endif
        
        return
        end

          
