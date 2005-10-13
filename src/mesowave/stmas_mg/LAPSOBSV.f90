SUBROUTINE LAPSOBSV(m)

!==========================================================
!  This routine reads in LSO observation data from LAPS.
!
!  HISTORY:
! 	Creation: YUANFU XIE	8-2005
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: m	! Maximum number of sites

  ! Local variables:
  CHARACTER*24 :: tim		! Obs file time
  CHARACTER :: stn(m)*20	! Station names
  CHARACTER :: prd(m)*11	! Provider names
  CHARACTER :: pwx(m)*25	! Present weather
  CHARACTER :: rtp(m)*6		! Report type
  CHARACTER :: stp(m)*6		! Station type (manual/auto)
  CHARACTER :: amt(m,5)*4	! cloud amount

  INTEGER :: nog,nob		! Number obs over grid/box
  INTEGER :: wid		! WMO id
  INTEGER :: otm(m)		! Observation time

  REAL*4 :: lat(m),lon(m), &	! Lat/Lon
	    elv(m)		! Elevation
  REAL*4 :: tmp(m),tmpea(m), &	! Temperature/expected accuracy
	    dew(m),dewea(m), &	! Dewpoint/EA
	    rhd(m),rhdea(m), &	! Relative Humidity/EA
	    wdi(m),wdiea(m), &	! Wind direction/EA
	    spd(m),spdea(m), &	! Wind speed/EA
	    gdi(m),	     &	! Gust wind direction
	    gsp(m),          &	! Gust wind speed
	    alt(m),altea(m), &	! Altimeter/EA
	    spr(m),prsea(m), &	! Station pressure/EA
	    msp(m), &		! Mean sea level pressure/EA
	    pcc(m),pccea(m), &	! 3-hour pressure change character
	    pch(m),pchea(m), &	! 3-hour pressure change
	    vis(m),visea(m), &	! Visibility/EA
	    sol(m),solea(m), &	! Solar/EA
	    slt(m),sltea(m), &	! Soil/water temperature/EA
	    slm(m),slmea(m), &	! Soil moist/EA
	    pc1(m),pcpea(m), &	! 1-hour precipitation/EA
	    pc3(m),pc6(m), &
	    p24(m), &		! 3,6,24-hour precipitation
	    snw(m),snwea(m), &	! Snow depth/EA
	    cld(m), &		! Number of cloud layers
	    mxt(m),mnt(m)	! 24-hour maximum/minimum temperature
  REAL*4 :: cht(m,5)		! cloud layer heights
  
  INTEGER :: i,j,k,err,iwv
  INTEGER :: hrs,mns,nit	! Time: hours, minutes and mid-night
  REAL :: xyt(3)		! X, Y and T
  REAL :: prs,ALT_2_SFC_PRESS

  ! Read observation data by LAPS time frames:
  numobs = 0
  DO i=1,numtmf

    ! Frame by frame: READ_SURFACE_DATAQC or READ_SURFACE_DATA
    CALL READ_SURFACE_DATA(i4prev(i),tim,nog,nob, &
	otm,wid,stn,prd,pwx,rtp,stp,lat,lon,elv, &
	tmp,dew,rhd,wdi,spd,gdi,gsp,alt,spr,msp,pcc,pch, &
	vis,sol,slt,slm,pc1,pc3,pc6,p24,snw,cld,mxt,mnt, &
	tmpea,dewea,rhdea,wdiea,spdea,altea,prsea,visea, &
	solea,sltea,slmea,pcpea,snwea,amt,cht,mxstts,err)
    IF (err .NE. 1) THEN
      ! LSO data cannot be read in:
      WRITE(*,21) i
    ELSE
      ! Assign LSO data to the corresponding arrays:

      ! Check if time cross midnight: 0 no cross; 1 cross
      nit = 0
      ! Assume time is not lapsed over 12 hours:
      IF (MAXVAL(otm(1:nob))-MINVAL(otm(1:nob)) &
	.GT. 1200.00) nit = 1

      ! Retrieve location and time sequence:
      DO j=1,nob	! Through all obs sites
	! X and Y:
	CALL LATLON_TO_RLAPSGRID(lat(j),lon(j), &
		latgrd,longrd,numgrd(1),numgrd(2), &
		xyt(1),xyt(2),err)

	! T: from LAPS time form: HHMM to seconds
	hrs = otm(j)/100
	mns = otm(j)-hrs*100

	! If cross midnight, set morning hour to 24+:
	IF ((nit .EQ. 1) .AND. (hrs .LE. 11)) hrs = 24+hrs
	! Cross midnight:
	xyt(3) = hrs*3600+mns*60

	! Pass the location/time to obs arrays:
        DO k=1,numvar
	  rawobs(2:4,j+numobs(k),k) = xyt(1:3)
	ENDDO
      ENDDO

      ! Place the observations into right variables:
      DO j=1,numvar
	SELECT CASE (varnam(j))
	CASE ("TEMP")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = tmp(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = tmpea(1:nob)
	CASE ("DEWP")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = dew(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = dewea(1:nob)
	CASE ("VISB")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = vis(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = visea(1:nob)
	CASE ("REDP")
	  DO k=1,nob
	    ! Collect either station pressure or altimeter:
	    IF ((spr(k) .NE. mising) .AND. &
		(spr(k) .NE. badsfc)) THEN
	      prs = spr(k)
	    ELSEIF ((alt(k) .NE. mising) .AND. &
		    (alt(k) .NE. badsfc)) THEN
	      prs = ALT_2_SFC_PRESS(alt(k),elv(k))
	    ELSE
	      prs = badsfc
	    ENDIF

	    ! Convert to reduced pressure:
	    IF (prs .NE. badsfc) THEN
	      CALL REDUCE_P(tmp(k),dew(k),prs,elv(k), &
		lapses(1),lapses(2), &
		rawobs(1,numobs(j)+k,j),0.0,badsfc)
	    ELSE
	      rawobs(1,numobs(j)+k,j) = badsfc
            ENDIF
	  ENDDO
	  weight(1+numobs(j):nob+numobs(j),j) = 1.0 ! altea(1:nob)
	CASE ("WNDU")
	  ! Find the index for v component:
	  iwv = 0
	  DO k=1,numvar
	    IF (varnam(k) .EQ. "WNDV") iwv = k
	  ENDDO
	  IF (iwv .EQ. 0) THEN
	    WRITE(*,24)
	  ENDIF
24 FORMAT('STMAS>LAPSOBS: Warning: no v component wind analysis!')

	  ! Convert wind from direction/speed to U/V:
	  DO k=1,nob
	    IF ((wdi(k) .EQ. mising) .OR. &
	        (wdi(k) .EQ. badsfc) .OR. &
	        (spd(k) .EQ. mising) .OR. &
		(spd(k) .EQ. badsfc)) THEN
	      rawobs(1,numobs(j)+k,j) = badsfc
	    ELSE
	      ! Conversion:
	      CALL DISP_TO_UV(wdi(k),spd(k),xyt(1),xyt(2))
	      CALL UVTRUE_TO_UVGRID(xyt(1),xyt(2), &
		rawobs(1,numobs(j)+k,j), &
		rawobs(1,numobs(j)+k,iwv),lon(k))
	    ENDIF
	  ENDDO
	  weight(1+numobs(j):nob+numobs(j),j) = 1.0
	  weight(1+numobs(j):nob+numobs(j),iwv) = 1.0
	CASE ("WNDV")
	  ! V should be in already. See CASE ("WNDU").
	CASE DEFAULT
	  WRITE(*,22) varnam(j)
	END SELECT
      ENDDO
    ENDIF

    ! Update frm:
    numobs(1:numvar) = numobs(1:numvar)+nob
  ENDDO
21 FORMAT('STMAS>LAPSOBS: Cannot read in LSO data: ',i8)
22 FORMAT('STMAS>LAPSOBS: No such var in LSO data: ',A4)

  ! Check number of obs:
  DO i=1,numvar
    IF (verbal .EQ. 1) WRITE(*,23) varnam(i),numobs(i)
  ENDDO
23 FORMAT('STMAS>LAPSOBSV: NumObs of (raw) ',A4,': ',I8)

  ! Remove invalid data:
  CALL RmvInvld

  ! Remove redundant observations:
  ! CALL RmvDupls

  ! Check the data ranges:
  IF (verbal .EQ. 1) THEN
    DO i=1,numvar
      WRITE(*,25) varnam(i), &
        MINVAL(rawobs(4,1:numobs(i),i)), &
        MAXVAL(rawobs(4,1:numobs(i),i)), &
	(MAXVAL(rawobs(4,1:numobs(i),i))- &
	 MINVAL(rawobs(4,1:numobs(i),i)))/3600.00
      WRITE(*,26) varnam(i), &
		  MINVAL(rawobs(1,1:numobs(i),i)), &
		  MAXVAL(rawobs(1,1:numobs(i),i))
    ENDDO
  ENDIF
25 FORMAT('STMAS>LAPSOBSV: ',A4,' obs time interval: ', &
     2F11.2,/,'STMAS>LAPSOBSV: Time length: ',F4.2,' hours')
26 FORMAT('STMAS>LAPSOBSV: ',A4,' min/max values: ', 2F11.2)

  ! Write out requested obs for testing:
  IF (savdat .EQ. 1) THEN
  OPEN(unit=10,file='STMAS_ob1.dat',form='formatted')
    WRITE(10,*) numobs(saveid),numtmf,domain,grdspc(3)
    WRITE(10,*) rawobs(1:4,1:numobs(saveid),saveid)
    CLOSE(10)
  ENDIF

END SUBROUTINE LAPSOBSV

SUBROUTINE RmvInvld

!==========================================================
!  This routine removes the invalid data (missing/badsfc)
!  data from observations.
!
!  HISTORY:
!	Creation: 8-2005 by YUANFU XIE
!==========================================================

  IMPLICIT NONE

  INTEGER :: i,j,num

  ! For every variables:
  DO i=1,numvar

    num = numobs(i)
    numobs(i) = 0
    ! For every obs:
    DO j=1,num
      IF ((rawobs(1,j,i) .NE. mising) .AND. &
	  (rawobs(1,j,i) .NE. badsfc)) THEN
	! Valid data:
	numobs(i) = numobs(i)+1
	rawobs(1:4,numobs(i),i) = rawobs(1:4,j,i)
	weight(numobs(i),i) = weight(j,i)
      ENDIF
    ENDDO
  ENDDO
 
  ! Check numbers of obs left:
  DO i=1,numvar
    IF (verbal .EQ. 1) WRITE(*,31) varnam(i),numobs(i)
  ENDDO
31 FORMAT('STMAS>LAPS_QCs: NumObs of (Vld) ',A4,': ',I8)

END SUBROUTINE RmvInvld

SUBROUTINE RmvDupls

!==========================================================
!  This routine removes the duplicated observations.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: i,j,k,l,nob
  REAL :: dis

  obsspc = 1.0e10
  DO i=1,numvar

    ! Use weight as flag:
    DO j=1,numobs(i)
      DO k=j+1,numobs(i)
	dis = 0.0
	DO l=2,4
	  dis = dis+(rawobs(l,j,i)-rawobs(l,k,i))* &
		    (rawobs(l,j,i)-rawobs(l,k,i))
	ENDDO

	! Mark those redundants:
	IF (dis .LT. epsiln) THEN
	  weight(k,i) = 0.0
	  IF (verbal .EQ. 1) THEN
	    WRITE(*,11) rawobs(1:4,j,i),rawobs(1:4,k,i), &
		varnam(i),j,k
	  ENDIF
	ELSE
	  ! Minimal observation spacing:
          IF ((ABS(rawobs(2,j,i)-rawobs(2,k,i)) .GT. 0.0) .AND. &
              (ABS(rawobs(2,j,i)-rawobs(2,k,i)) .LT. obsspc(1,i))) &
            obsspc(1,i) = ABS(rawobs(2,j,i)-rawobs(2,k,i))
          IF ((ABS(rawobs(3,j,i)-rawobs(3,k,i)) .GT. 0.0) .AND. &
              (ABS(rawobs(3,j,i)-rawobs(3,k,i)) .LT. obsspc(2,i))) &
            obsspc(2,i) = ABS(rawobs(3,j,i)-rawobs(3,k,i))
          IF ((ABS(rawobs(4,j,i)-rawobs(4,k,i)) .GT. 0.0) .AND. &
              (ABS(rawobs(4,j,i)-rawobs(4,k,i)) .LT. obsspc(3,i))) &
            obsspc(3,i) = ABS(rawobs(4,j,i)-rawobs(4,k,i))
	ENDIF
      ENDDO
11 FORMAT('STMAS>RmvDupls: Redundant data: ',/,4F14.4,/, &
	4F14.4,/,A6,2I8)

    ENDDO
    WRITE(*,12) varnam(i),obsspc(1:3,i)
12 FORMAT('STMAS>RmvDupls: Minimal obs (',A4,') spacing: ',3E12.4)

    ! Remove redundants:
    nob = 0
    DO j=1,numobs(i)
      IF (weight(j,i) .GT. 0.0) THEN
        nob = nob+1
        rawobs(1:4,nob,i) = rawobs(1:4,j,i)
	weight(nob,i) = weight(j,i)
      ENDIF
    ENDDO
    numobs(i) = nob

  ENDDO

END SUBROUTINE RmvDupls
