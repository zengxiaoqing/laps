
!     subroutine fireweatherindex(tmk,prs,qvp,ght,ter,sfp,tgk,
!    &                            u10,v10,xlus,fwi,miy,mjx,mkzh)

      subroutine fireweatherindex(t_sfc_k,rh_sfc,p_sfc_mb,u10,v10
     1                           ,miy,mjx
     1                           ,fwi)

      include 'comconst'

      dimension ! tmk(miy,mjx,mkzh), prs(miy,mjx,mkzh), qvp(miy,mjx,mkzh),
!    &          ght(miy,mjx,mkzh),
!    &          ter(miy,mjx), sfp(miy,mjx), tgk(miy,mjx),
     &          u10(miy,mjx), v10(miy,mjx), fwi(miy,mjx),
     $          xlus(miy,mjx)

      real t_sfc_k(miy,mjx)
      real rh_sfc(miy,mjx)
      real p_sfc_mb(miy,mjx)
 
      real t2k, t2f, ht, prs2, q, e, es, rh2, m, n, k_to_f
      logical debug

C      debug = .true.   ! Debug mode
      debug = .false.   ! Production mode

      do 200 j = 1, mjx ! -1
      do 100 i = 1, miy ! -1

        if(.false.)then ! original Seattle code with sigma coordinate inputs

!         if (nint(xlus(i,j)).eq.iwater) then
!          fwi(i,j) = -3.
!          goto 100
!         endif
!         t2k = (tmk(i,j,mkzh) + tgk(i,j)) / 2.
!         t2f = (t2k - celkel) * 1.8 + 32.
!         ht = ght(i,j,mkzh) - ter(i,j)
!         prs2 = (ht-2.)/ht * (sfp(i,j)-prs(i,j,mkzh)) + prs(i,j,mkzh)
!         q = .001 * qvp(i,j,mkzh)  ! g/kg to g/g
!         e = q*prs2/(eps+q)
!         es = ezero * exp( eslcon1*(t2k-celkel)/(t2k-eslcon2) )
!         rh2 = 100.*(e*(prs2-es))/(es*(prs2-e))

        else ! new version to use sfc inputs from LAPS analysis
          rh2 = rh_sfc(i,j)
          prs2 = p_sfc_mb(i,j)
          t2f = k_to_f(t_sfc_k(i,j))

        endif

C---------Use 2.237 to convert m/s to mph
        uuu10 = u10(i,j) * 2.237
        vvv10 = v10(i,j) * 2.237
        if( rh2 .le. 10.5 ) then
           m = 0.03229 + (0.281073 * rh2) - (0.000578 * rh2 * t2f)
        else if( rh2 .gt. 10.5 .and. rh2 .le. 50.5 ) then
           m = 2.22749 + (0.160107 * rh2) - (0.014784 * t2f)
        else if( rh2 .gt. 50.5 .and. rh2 .le. 100 ) then
           m = 21.0606 + (0.005565 * rh2**2) - (0.00035 * rh2 * t2f)
     &        - (0.483199 * rh2)
        else         !   rh2 > 100
           m = 21.0606 + (0.005565 * 100**2) - (0.00035 * 100 * t2f)
     &        - (0.483199 * 100)
           print *, 'fwi calculation: rh2 > 100 ', j, i, rh2, t2f
        endif

        n = 1.0 - 2.0*(m/30.) + 1.5*(m/30.)**2 - 0.5*(m/30.)**3

        fwi(i,j) = (n * sqrt(1.0 + uuu10**2 + vvv10**2)) / 0.3002

!       if( debug .and. mod(i,20) .eq. 0 .and. mod(j,20) .eq. 0 ) then
!          print *, j, i, t2k, t2f, prs2, rh2, fwi(i,j), uuu10, vvv10 
!       endif

 100  continue
 200  continue

      return
      end
