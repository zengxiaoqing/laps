     


      subroutine sigma_to_p (ptop, psfc, sigma, p )

c     The purpose of this routine is to take an input sigma level, and its
c     associated pressure top and surface pressure and convert that 
c     to pressure dependent on the reference level pressure.

      real ptop
      real psfc
      real sigma
      real p

      p = sigma*(psfc-ptop) + ptop

      return
      end
