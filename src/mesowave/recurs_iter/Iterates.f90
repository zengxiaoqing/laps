SUBROUTINE Iterates(id,bkgd,ldf,nx,ny,ds,ncycles,nvlaps,nfic)

!*************************************************
!  This routine iteratively solves data analysis
!  problem.
!
!  HISTORY: FEB. 2004 by YUANFU XIE.
!           Sep. 2004 by YUANFU XIE penalizing div.
!*************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: id,nx,ny,ncycles,nvlaps,nfic
  REAL,    INTENT(IN) :: bkgd(nx,ny,ncycles,nvlaps)
  REAL,    INTENT(IN) :: ldf(nx,ny),ds(3)

  INTEGER :: iter,iobs,i,j,k,no_v,idp,nbqc,numo
  REAL    :: y0,b(2,3),rms,stdv(nvlaps),amx,amn,tmx,tmn

  REAL :: oo(4,mobs)
  oo = o

  ! Unified analysis of velocity:
  idp = id
  IF (id .EQ. 201) idp = id+1

  s(1:n(1),1:n(2),1:n(3),id:idp) = 0.0

  DO iter=1,nrf(id)

     a(1:n(1),1:n(2),1:n(3),id:idp) = 0.0

     ! QC: bound check and compute standard deviations:
     IF (iter .EQ. 1) THEN

	stdv(id) = 0.0
	no_v = 0

	nbqc = 0
	numo = 0
	amx = -1.0e6
	amn = 1.0e6

        DO iobs=1,nobs

	   IF ((id .EQ. vid(iobs)) .AND. &
	       (idx(1,iobs) .GT. nfic) .AND. &
	       (idx(1,iobs) .LT. n(1)-nfic) .AND. &
	       (idx(2,iobs) .GT. nfic) .AND. &
               (idx(2,iobs) .LT. n(2)-nfic)) THEN

              y0 = 0.0

	      numo = numo+1
        
              b(1,1:3) = 1.0-coe(1:3,iobs)
              b(2,1:3) = coe(1:3,iobs)
              DO k=1,2
		 IF ((idx(3,iobs)+k-1 .GE. 1) .AND. &
		     (idx(3,iobs)+k-1 .LE. n(3))) THEN
                 DO j=1,2
                    DO i=1,2
                       Y0 = Y0 + bkgd(idx(1,iobs)+i-1-nfic, &
                                      idx(2,iobs)+j-1-nfic, &
                                      idx(3,iobs)+k-1,vid(iobs))* &
                                 b(i,1)*b(j,2)*b(k,3)
                    ENDDO   
                 ENDDO
		 ENDIF
              ENDDO

	      stdv(id) = stdv(id)+(o(1,iobs)-y0)**2
	      no_v = no_v+1

	      IF (ABS(o(1,iobs)-y0) .GT. qc_cons(id)) THEN
		 PRINT*,'Threshold QC: ',o(1,iobs),y0,id,o(2:4,iobs)
		 o(1,iobs) = y0
		 w(iobs) = 0.0
	      ENDIF

	      ! Bad QC:
              IF ((id .EQ. 1) .OR. (id .EQ. 5)) THEN
                 IF (ABS(o(1,iobs)-y0) .GT. 10.0) THEN
	            !PRINT*,'Bad QC: ',o(1,iobs),y0, &
			! vid(iobs),iobs,o(2:4,iobs)
	            o(1,iobs) = y0
                    w(iobs) = 0.0
		    nbqc = nbqc+1
	         !ELSE
		 !   PRINT*,'God QC: ',o(1,iobs),y0, &
		 !	vid(iobs),iobs,o(2:4,iobs)
                 ENDIF
	      ENDIF

	   ENDIF

        ENDDO

	PRINT*,'Total Bad QC: ',nbqc,id,numo
	if (id .EQ. 4) PRINT*,'MAXimum diff: ',amx

	stdv(id) = SQRT(stdv(id)/no_v)

	!IF (no_v .GT. 0) THEN
        !   PRINT*,'Standard Deviation: ',stdv(id),id
	!ELSE
	!   PRINT*,'Standard Deviation: no observation'
	!ENDIF

	! Standard deviation check: with 4.0*stdv
	nbqc = 0
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
		 IF ((idx(3,iobs)+k-1 .GE. 1) .AND. &
		     (idx(3,iobs)+k-1 .LE. n(3))) THEN
                 DO j=1,2
                    DO i=1,2
                       Y0 = Y0 + bkgd(idx(1,iobs)+i-1-nfic, &
                                      idx(2,iobs)+j-1-nfic, &
                                      idx(3,iobs)+k-1,vid(iobs))* &
                                 b(i,1)*b(j,2)*b(k,3)
                    ENDDO   
                 ENDDO
		 ENDIF
              ENDDO

              IF (ABS(o(1,iobs)-y0) .GT. 3.0*stdv(id)) THEN
	         !PRINT*,'Standard deviation QC: ', &
	    	 !   o(1,iobs),y0,vid(iobs),iobs,stdv(id)
	         nbqc = nbqc+1
	         o(1,iobs) = y0
                 w(iobs) = 0.0
	      !ELSE
	      !  PRINT*,'PASS STD: ',iobs,o(1,iobs),y0,3.0*stdv(id)
              ENDIF

	   ENDIF

        ENDDO

	PRINT*,'Number of STD QCed: ',nbqc,id,stdv(id),numo

     ENDIF

     IF (iter .GT. 1) THEN
        CALL Minimize(id,ds)
     ELSE
	a(nfic+1:n(1)-nfic,nfic+1:n(2)-nfic,1:n(3),id:idp) = &
                                 bkgd(1:nx,1:ny,1:n(3),id:idp)

	! Fictitious points:
	DO i=1,nfic
	   a(i,nfic+1:n(2)-nfic,1:n(3),id:idp) = &
                      a(nfic+1,nfic+1:n(2)-nfic,1:n(3),id:idp)
	   a(n(1)-nfic+i,nfic+1:n(2)-nfic,1:n(3),id:idp) = &
                      a(n(1)-nfic,nfic+1:n(2)-nfic,1:n(3),id:idp)
	ENDDO
	DO i=1,nfic
	   a(1:n(1),i,1:n(3),id:idp) = a(1:n(1),nfic+1,1:n(3),id:idp)
	   a(1:n(1),n(2)-nfic+i,1:n(3),id:idp) = &
             a(1:n(1),n(2)-nfic,1:n(3),id:idp)
	ENDDO
     ENDIF

     no_v = 0
     rms = 0.0

     DO iobs=1,nobs

	IF ((id .EQ. vid(iobs)) .OR. (idp .EQ. vid(iobs))) THEN

           y0 = 0.0
        
           b(1,1:3) = 1.0-coe(1:3,iobs)
           b(2,1:3) = coe(1:3,iobs)
           DO k=1,2
	      IF ((idx(3,iobs)+k-1 .GE. 1) .AND. &
		  (idx(3,iobs)+k-1 .LE. n(3))) THEN
              DO j=1,2
                 DO i=1,2
                    Y0 = Y0 + a(idx(1,iobs)+i-1,idx(2,iobs)+j-1, &
                                idx(3,iobs)+k-1,vid(iobs))* &
                              b(i,1)*b(j,2)*b(k,3)
                 ENDDO
              ENDDO
	      ENDIF
           ENDDO
        
           o(1,iobs) = o(1,iobs)-y0

	   !rms = rms + o(1,iobs)*o(1,iobs)
           !no_v = no_v + 1

	ENDIF

     ENDDO

     !IF (no_v .NE. 0) THEN
     !   WRITE(11,*) 'RMS: ', SQRT(rms/no_v),iter,al(1:3,id)
     !ELSE
     !   WRITE(11,*) 'RMS: ', SQRT(rms),iter,al(1:3,id)
     !ENDIF

     ! Filter:
     IF (iter .GT. 1) THEN
        IF (id .EQ. 6) THEN
	   al(1:3,id) = al(1:3,id)*0.9
        ELSE
           al(1:3,id:idp) = al(1:3,id:idp)*0.9
	ENDIF
     ENDIF

     ! Accumulate:
     s(1:n(1),1:n(2),1:n(3),id:idp) = &
	s(1:n(1),1:n(2),1:n(3),id:idp)+ &
        a(1:n(1),1:n(2),1:n(3),id:idp)

  ENDDO

  ! Temporature problem:
  IF (id .EQ. 1) THEN

     ! Compute the max/min differences from background:
     tmx = MAXVAL(s(nfic+1:nfic+nx,nfic+1:nfic+ny,n(3),1)- &
	       bkgd(1:nx,1:ny,n(3),1))
     tmn = MINVAL(s(nfic+1:nfic+nx,nfic+1:nfic+ny,n(3),1)- &
	       bkgd(1:nx,1:ny,n(3),1))

     PRINT*,'MAX DIFF IN T: ', tmx, tmn

     IF ((tmx .GT. 15.0e0) .OR. (tmn .LT. -15.0e0)) THEN

	OPEN(unit=11,file='badata.dat')
	WRITE(11,*) nx,ny,n,dm(1,3),dm(2,3)
	WRITE(11,*) bkgd(1:n(1),1:n(2),1:n(3),1)
	write(11,*) ldf(1:nx,1:ny)
	WRITE(11,*) nobs
	WRITE(11,*) oo(1:4,1:nobs),vid(1:nobs),w(1:nobs)
	WRITE(11,*) s(1:n(1),1:n(2),1:n(3),1)
	CLOSE(11)
	
        ! DO j=1,ny
        !   DO i=1,nx
	!      s(nfic+i,nfic+j,1:n(3),id) = &
        !        bkgd(i,j,1:n(3),id)
        !   ENDDO
        ! ENDDO

     ENDIF
  ENDIF

  ! Land/water weight:
  IF ((id .NE. 6) .AND. (id .NE. 4)) THEN
     DO j=1,ny
        DO i=1,nx
  	   s(nfic+i,nfic+j,1:n(3),id:idp) = &
		ldf(i,j)*s(nfic+i,nfic+j,1:n(3),id:idp)+ &
                 (1.0-ldf(i,j))*bkgd(i,j,1:n(3),id:idp)
        ENDDO
     ENDDO
  ELSE 
     ! DO j=1,ny
     !    DO i=1,nx
        ! s(nfic+i,nfic+j,1:n(3),id) = bkgd(i,j,1:n(3),id)
     !    ENDDO
     ! ENDDO
  ENDIF

END SUBROUTINE Iterates
