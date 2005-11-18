SUBROUTINE LAPS_QCs

!==========================================================
!  This routine runs quality control over data by threshold
!  values and standard deviation.
!
!  HISTORY:
! 	Creation: YUANFU XIE	6-2005
!==========================================================

  IMPLICIT NONE

  ! Interpolation indices and coefficients:
  CALL LAPSIntp

  ! Optional QCs:
  IF (qc_val .EQ. 1) CALL Thrshold

  ! Save QCed obs:
  CALL CpyQCObs

END SUBROUTINE LAPS_QCs

SUBROUTINE CpyQCObs

!==========================================================
!  This routine copies QCed observation data from rawobs to
!  qc_obs after all QC is done. This routine can be avoid
!  if the qc_obs array keeps the same structure as rawobs.
!
!  HISTORY:
!	Creation: 8-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: i,j

  ! Copy:
  DO i=1,numvar
    DO j=1,numobs(i)
      qc_obs(1:4,j,i) = rawobs(1:4,j,i)

      ! Save innovations:
      IF (needbk(i) .EQ. 1) &
        qc_obs(1,j,i) = qc_obs(1,j,i)-bkgobs(j,i)

    ENDDO
  ENDDO

END SUBROUTINE CpyQCObs

SUBROUTINE Thrshold

!==========================================================
!  This routine does the threshold value QC checks.
!
!  HISTORY:
!	Creation: 8-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: i,j,num

  ! Check:
  DO i=1,numvar
    IF (needbk(i) .EQ. 1) THEN
      num = numobs(i)
      numobs(i) = 0
      DO j=1,num

        ! QC check: avoid bkg = mising with roundoff error:
        IF (ABS(rawobs(1,j,i)-bkgobs(j,i)) .LE. thresh(i)) THEN
	  numobs(i) = numobs(i)+1
	  rawobs(1:4,numobs(i),i) = rawobs(1:4,j,i)
	  weight(numobs(i),i) = weight(j,i)
	  indice(1:6,numobs(i),i) = indice(1:6,j,i)
	  coeffs(1:6,numobs(i),i) = coeffs(1:6,j,i)
	  bkgobs(numobs(i),i) = bkgobs(j,i)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
 
  ! Check numbers of obs left:
  DO i=1,numvar
    IF (verbal .EQ. 1) WRITE(*,31) varnam(i),numobs(i)
  ENDDO
31 FORMAT('STMAS>LAPS_QCs: NumObs of (Vlu) ',A4,': ',I8)

END SUBROUTINE Thrshold
