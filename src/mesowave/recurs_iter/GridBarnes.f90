SUBROUTINE GridBarnes(f,l,n,b)

!==========================================================
!  This routine applied a Barnes analysis to a gridded two
!  dimensional function as a filter.
!
!  HISTORY: JAN. 2005 by YUANFU XIE.
!==========================================================

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: l(2),n(2)
   REAL,    INTENT(IN) :: f(l(1),l(2))
   REAL,    INTENT(OUT) :: b(l(1),l(2))

   ! Local variables:
   INTEGER :: i,j,ir,ni,nj,ip
   REAL    :: gama,s,w,ws,t(l(1),l(2)),tmp

   gama = 0.2
   s = 2.0e3/gama  ! kapa_0
   ir = 60	   ! neighbors
   
   ! Clean analysis field:
   b = 0.0
   t = 0.0

   ! Barnes analysis:
   DO ip=1,1

      s = s*gama

      ! For every gridpoint:
      DO j=1,n(2)
         DO i=1,n(1)
	    
	    ! For each neighbor within ir:
	    ws = 0.0
	    tmp = 0.0
	    DO nj=-ir,ir,2
	       DO ni=-ir,ir,2
		  w = EXP(-(FLOAT(ni)**2+FLOAT(nj)**2)/s)

		  IF ((i+ni .GE. 1) .AND. (i+ni .LE. n(1)) .AND. &
		      (j+nj .GE. 1) .AND. (j+nj .LE. n(2))) THEN
		     tmp = tmp+(f(i+ni,j+nj)-t(i+ni,j+nj))*w
		     ws = ws + w
		  ENDIF

	       ENDDO
	    ENDDO

	    b(i,j) = b(i,j)+tmp/ws
         ENDDO
      ENDDO

      t = b

   ENDDO

END SUBROUTINE GridBarnes
