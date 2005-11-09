SUBROUTINE LAPSBKGD

!==========================================================
!  This routine reads into LGA background fields from LAPS.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	8-2005
!==========================================================

  IMPLICIT NONE

  CHARACTER*31 :: ext	! extension of bkg file used (LAPS)
  INTEGER :: i,j,err	! err: 1 normal; 0 file not found
  INTEGER :: tim	! time of bkg file used (LAPS)
  INTEGER :: iwv	! index for v component of wind

  ! Set LAPS circle time frames:
  DO i=0,numtmf-1
    i4prev(numtmf-i) = i4time-i*lapsdt
  ENDDO

  ! Check LAPS time frames:
  IF (verbal .EQ. 1) THEN
    DO i=1,numtmf
      WRITE(*,11) i,i4prev(i),MOD(i4prev(i),86400)
    ENDDO
  ENDIF
11 FORMAT('STMAS>LAPSBKGD: LAPS Time stamp',I3,':',I11,I7)

  ! Read background fields:
  DO j=1,numvar

    ! Wind:
    IF (varnam(j) .EQ. 'WNDU') THEN
      ! Search index for v component of wind:
      iwv = 0
      DO i=1,numvar
        IF (varnam(i) .EQ. 'WNDV') iwv = i
      ENDDO
      IF (iwv .EQ. 0) THEN
	WRITE(*,12)
	STOP
      ENDIF
12 FORMAT('STMAS>LAPSBKGD: V component of wind is missing!')

      ! Get wind fields:
      DO i=1,numtmf
	CALL GET_BKGWIND_SFC(i4prev(i),ext,tim, &
	  bkgrnd(1,1,i,j),bkgrnd(1,1,i,iwv),lapsdt, &
	  numgrd(1),numgrd(2),err)
	IF (err .EQ. 0) WRITE(*,13) varnam(j),i4prev(i),i,j
      ENDDO

    ! Other fields:
    ELSE IF ((varnam(j) .NE. 'WNDV') .AND. &	! V in with U
	     (varnam(j) .NE. 'CEIL')) THEN	! No ceiling bkg
      DO i=1,numtmf
        CALL GET_BACKGROUND_SFC(i4prev(i),varnam(j),ext,tim, &
	  bkgrnd(1,1,i,j),lapsdt,numgrd(1),numgrd(2),err)
	IF (err .EQ. 0) WRITE(*,13) varnam(j),i4prev(i),i,j
      ENDDO
    ELSE
      IF (needbk(j) .EQ. 0) &
        bkgrnd(1:numgrd(1),1:numgrd(2),1:numtmf,j) = 0.0
    ENDIF
  ENDDO
13 FORMAT('STMAS>LAPSBKGD: Background is not found for: ',A4,i16,2i3)

  ! Save background if requested for debugging:
  IF (savdat .EQ. 1) THEN
    OPEN(unit=10,file='STMAS_bkg.dat',form='formatted')
    WRITE(10,*) numgrd(1:2),numtmf,domain
    WRITE(10,*) bkgrnd(1:numgrd(1),1:numgrd(2),1:numtmf,saveid)
    CLOSE(10)
  ENDIF

END SUBROUTINE LAPSBKGD
