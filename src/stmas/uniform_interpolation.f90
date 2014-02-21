!>
!! This is a routine for an interpolation between two uniform grids over the same domain
!!
!! \author Yuanfu Xie
!! \b History: Feb. 2014
!

SUBROUTINE uniform_interpolation(nfrom,nto,vfrom,vto)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nfrom(3),nto(3)
  REAL,    INTENT(IN) :: vfrom(nfrom(1),nfrom(2),nfrom(3))
  REAL,    INTENT(OUT) :: vto(nto(1),nto(2),nto(3))

  ! Local variables:
  INTEGER :: indx(2,3),i,j,k,ii,jj,kk
  REAL    :: coef(2,3)

  DO k=1,nto(3)
    coef(1,3) = FLOAT(nfrom(3)-1)*(k-1)/FLOAT(nto(3)-1)+1.0 ! real position from (1,nfrom)
    indx(1,3) = INT(coef(1,3))                              ! left grid
    coef(2,3) = coef(1,3)-indx(1,3)                         ! distance to left: weight on right
    indx(2,3) = MIN(indx(1,3)+1,nfrom(3))                   ! right grid
    coef(1,3) = 1.0-coef(2,3)                               ! distance to right: weight to left
  DO j=1,nto(2)
    coef(1,2) = FLOAT(nfrom(2)-1)*(j-1)/FLOAT(nto(2)-1)+1.0 ! real position from (1,nfrom)
    indx(1,2) = INT(coef(1,2))                              ! left grid
    coef(2,2) = coef(1,2)-indx(1,2)                         ! distance to left: weight on right
    indx(2,2) = MIN(indx(1,2)+1,nfrom(2))                   ! right grid
    coef(1,2) = 1.0-coef(2,2)                               ! distance to right: weight to left
  DO i=1,nto(1)
    coef(1,1) = FLOAT(nfrom(1)-1)*(i-1)/FLOAT(nto(1)-1)+1.0 ! real position from (1,nfrom)
    indx(1,1) = INT(coef(1,1))                              ! left grid
    coef(2,1) = coef(1,1)-indx(1,1)                         ! distance to left: weight on right
    indx(2,1) = MIN(indx(1,1)+1,nfrom(1))                   ! right grid
    coef(1,1) = 1.0-coef(2,1)                               ! distance to right: weight to left

    ! Interpolation:
    vto(i,j,k) = 0.0
    DO kk=1,2
    DO jj=1,2
    DO ii=1,2
      vto(i,j,k) = vto(i,j,k)+coef(ii,1)*coef(jj,2)*coef(kk,3)* &
                              vfrom(indx(ii,1),indx(jj,2),indx(kk,3))
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  ENDDO

END SUBROUTINE uniform_interpolation
