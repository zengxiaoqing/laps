
      subroutine hainesindex(prs,tmk,tdk,haines,miy,mjx,mkzh,prsb,prst)
      include 'comconst'

!     Note that this version differs from the Seattle version to account
!     for flipped vertical coordinates

      dimension tmk(miy,mjx,mkzh), tdk(miy,mjx,mkzh), prs(miy,mjx,mkzh),
     &          haines(miy,mjx)

       do j = 1, mjx ! -1
       do i = 1, miy ! -1
         
         if( prs(i,j,1) .lt. prsb ) then
            haines(i,j) = 0.
         else
            do k = 2, mkzh

              kk  = mkzh+1 - k     ! Account for flipped vertical coordinates
              km1 = kk + 1         ! Account for flipped vertical coordinates

              if( prs(i,j,kk).gt.prst .and. prs(i,j,km1).le.prst ) then       
                 tmkt = tmk(i,j,km1) + (tmk(i,j,kk)-tmk(i,j,km1)) *
     &                  (log(prst)-log(prs(i,j,km1))) /
     &                  (log(prs(i,j,kk))-log(prs(i,j,km1)))
              endif
              if( prs(i,j,kk).gt.prsb .and. prs(i,j,km1).le.prsb ) then
                 tmkb = tmk(i,j,km1) + (tmk(i,j,kk)-tmk(i,j,km1)) *
     &                  (log(prsb)-log(prs(i,j,km1))) /
     &                  (log(prs(i,j,kk))-log(prs(i,j,km1)))
                 tdkb = tdk(i,j,km1) + (tdk(i,j,kk)-tdk(i,j,km1)) *
     &                  (log(prsb)-log(prs(i,j,km1))) /
     &                  (log(prs(i,j,kk))-log(prs(i,j,km1)))
              endif
            enddo

            deltat = tmkb - tmkt
            dpdep  = tmkb - tdkb

            if( nint(prsb) .eq. 700 ) then    ! hainh
               if( deltat .le. 17.5 ) then
                  factor1 = 1.
               else if(deltat .gt. 17.5 .and. deltat .le. 21.5 ) then
                  factor1 = 2.
               else             ! deltat .gt. 21.5
                  factor1 = 3.
               endif
   
               if( dpdep .le. 14.5 ) then
                  factor2 = 1.
               else if( dpdep .gt. 14.5 .and. dpdep .le. 20.5 ) then
                  factor2 = 2.
               else             ! dpdep .gt. 20.5
                  factor2 = 3.
               endif
            endif

            if( nint(prsb) .eq. 850 ) then    ! hainm
               if( deltat .le. 5.5 ) then
                  factor1 = 1.
               else if(deltat .gt. 5.5 .and. deltat .le. 10.5 ) then
                  factor1 = 2.
               else             ! deltat .gt. 10.5
                  factor1 = 3.
               endif
   
               if( dpdep .le. 5.5 ) then
                  factor2 = 1.
               else if( dpdep .gt. 5.5 .and. dpdep .le. 12.5 ) then
                  factor2 = 2.
               else             ! dpdep .gt. 12.5
                  factor2 = 3.
               endif
            endif

            haines(i,j) = factor1 + factor2

C            if( mod(i,10) .eq. 0 .and. mod(j,10) .eq. 0 ) then
C               print *, j, i, tmkb, tmkt, deltat
C               print *, tmkb, tdkb, dpdep
C            endif

         endif

       enddo
       enddo

      return
      end
