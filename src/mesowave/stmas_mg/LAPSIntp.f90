SUBROUTINE LAPSIntp

!==========================================================
!  This routine interpolates gridpoints to observation site
!  and saves the indices and coefficients.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: i,j,ix,iy,it

  DO i=1,numvar
    CALL Grid2Obs(indice(1,1,i),coeffs(1,1,i), &
	rawobs(1,1,i),numobs(i),weight(1,i),numgrd, &
	grdspc,domain)

    ! Compute background values at observation sites:
    DO j=1,numobs(i)

      ! Interpolate background to the obs site:
      bkgobs(j,i) = 0.0
      DO it=3,6,3
        DO iy=2,5,3
          DO ix=1,4,3
            bkgobs(j,i) = bkgobs(j,i) + &
	      bkgrnd(indice(ix,j,i), &
		     indice(iy,j,i), &
		     indice(it,j,i),i)* &
	      coeffs(ix,j,i)*coeffs(iy,j,i)*coeffs(it,j,i)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  ENDDO

END SUBROUTINE LAPSIntp
