SUBROUTINE Iterates(id,bkgd,ldf,nx,ny,ncycles,nvlaps,nfic)

!*************************************************
!  This routine iteratively solves data analysis
!  problem.
!
!  HISTORY: FEB. 2004 by YUANFU XIE.
!*************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: id,nx,ny,ncycles,nvlaps,nfic
  REAL,    INTENT(IN) :: bkgd(nx,ny,ncycles,nvlaps)
  REAL,    INTENT(IN) :: ldf(nx,ny)

  INTEGER :: iter,iobs,i,j,k,no_v
  REAL    :: y0,b(2,3),rms

  s(1:n(1),1:n(2),1:n(3),id) = 0.0

  DO iter=1,nrf(id)

     a(1:n(1),1:n(2),1:n(3),id) = 0.0

     ! QC: bound check:
     IF (((id .EQ. 1) .OR. (id .EQ. 5)) .AND. (iter .EQ. 1)) THEN

        DO iobs=1,nobs

	   IF ((id .EQ. vid(iobs)) .AND. &
	       (idx(1,iobs) .GT. nfic) .AND. &
	       (idx(1,iobs) .LT. n(1)-nfic) .AND. &
	       (idx(2,iobs) .GT. nfic) .AND. &
               (idx(2,iobs) .LT. n(2)-nfic)) THEN

              y0 = 0.0
        
              b(1,1:3) = 1.0-coe(1:3,iobs)
              b(2,1:3) = coe(1:3,iobs)
              DO k=1,2
                 DO j=1,2
                    DO i=1,2
                       Y0 = Y0 + bkgd(idx(1,iobs)+i-1-nfic, &
                                      idx(2,iobs)+j-1-nfic, &
                                      idx(3,iobs)+k-1,vid(iobs))* &
                                 b(i,1)*b(j,2)*b(k,3)
                    ENDDO   
                 ENDDO
              ENDDO
        
              IF (ABS(o(1,iobs)-y0) .GT. 10.0) THEN
	         PRINT*,'Bad QC: ',o(1,iobs),y0,vid(iobs),iobs
	         o(1,iobs) = y0
                 w(iobs) = 0.0
              ENDIF

	   ENDIF

        ENDDO

     ENDIF

     IF (iter .GT. 1) THEN
        CALL Minimize(id)
     ELSE
	a(nfic+1:n(1)-nfic,nfic+1:n(2)-nfic,1:n(3),id) = &
                                 bkgd(1:nx,1:ny,1:n(3),id)
	! Fictitious points:
	DO i=1,nfic
	   a(i,nfic+1:n(2)-nfic,1:n(3),id) = &
                      a(nfic+1,nfic+1:n(2)-nfic,1:n(3),id)
	   a(n(1)-nfic+i,nfic+1:n(2)-nfic,1:n(3),id) = &
                      a(n(1)-nfic,nfic+1:n(2)-nfic,1:n(3),id)
	ENDDO
	DO i=1,nfic
	   a(1:n(1),i,1:n(3),id) = a(1:n(1),nfic+1,1:n(3),id)
	   a(1:n(1),n(2)-nfic+i,1:n(3),id) = &
             a(1:n(1),n(2)-nfic,1:n(3),id)
	ENDDO
     ENDIF

     no_v = 0
     rms = 0.0

     DO iobs=1,nobs

	IF (id .EQ. vid(iobs)) THEN

           y0 = 0.0
        
           b(1,1:3) = 1.0-coe(1:3,iobs)
           b(2,1:3) = coe(1:3,iobs)
           DO k=1,2
              DO j=1,2
                 DO i=1,2
                    Y0 = Y0 + a(idx(1,iobs)+i-1,idx(2,iobs)+j-1, &
                                idx(3,iobs)+k-1,vid(iobs))* &
                              b(i,1)*b(j,2)*b(k,3)
                 ENDDO
              ENDDO
           ENDDO
        
           o(1,iobs) = o(1,iobs)-y0

	   rms = rms + o(1,iobs)*o(1,iobs)
           no_v = no_v + 1

	ENDIF

     ENDDO

     !IF (no_v .NE. 0) THEN
     !   WRITE(11,*) 'RMS: ', SQRT(rms/no_v),iter,al(1:3,id)
     !ELSE
     !   WRITE(11,*) 'RMS: ', SQRT(rms),iter,al(1:3,id)
     !ENDIF

     ! Filter:
     IF (iter .GT. 1)al(1:3,id) = al(1:3,id)*0.8

     ! Accumulate:
     s(1:n(1),1:n(2),1:n(3),id) = s(1:n(1),1:n(2),1:n(3),id)+ &
                                  a(1:n(1),1:n(2),1:n(3),id)

  ENDDO

  ! Land/water weight:
  DO j=1,ny
     DO i=1,nx
  	s(nfic+i,nfic+j,1:n(3),id) = &
		ldf(i,j)*s(nfic+i,nfic+j,1:n(3),id)+ &
                (1.0-ldf(i,j))*bkgd(i,j,1:n(3),id)
     ENDDO
  ENDDO

END SUBROUTINE Iterates
