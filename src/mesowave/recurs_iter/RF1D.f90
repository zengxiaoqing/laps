SUBROUTINE RF1D(u,n,a,np)

!****************************************************
!  This routine applies the recursive filter to u for
!  np times with alpha value of a.
!
!  HISTORY: APR. 2003 by YUANFU XIE.
!****************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n,np
  REAL,    INTENT(IN) :: a
  REAL                :: u(n)

  ! Local variables:
  INTEGER :: i,ip
  REAL    :: one_a,r(n)

  one_a = 1.0-a

  ! Recurisve filter number of np times:
  DO ip=1,np
        
     ! Left to right:
     IF (ip .EQ. 1) THEN
        r(1) = one_a*u(1) 
     ELSE IF (ip == 2 ) THEN
        r(1) = u(1) / ( 1.0 + a )
     ELSE
        r(1) = one_a * ( u(1) - a**3 * u(2) ) / &
             (1.0-a**2)**2
     END IF
     DO i=2,n
        r(i) = a*r(i-1)+one_a*u(i)
     ENDDO
        
     ! Right to left:
     IF (ip .EQ. 1) THEN
        u(n) = r(n)/(1.0+a) 
     ELSE
        u(n) = one_a * ( r(n) - a**3 * r(n-1) ) / &
             (1.0-a**2)**2
     END IF
     DO i=n-1,1,-1
        u(i) = a*u(i+1)+one_a*r(i)
     ENDDO
  ENDDO

END SUBROUTINE RF1D
