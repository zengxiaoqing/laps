SUBROUTINE LAPSUnit

!==========================================================
!  This routine converts LSO observation units into a unit
!  consistent with the background.
!
!  HISTORY:
!	Creation: 8-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER :: i

  ! Check all variables:
  DO i=1,numvar

    ! Find necessary conversion:
    SELECT CASE (varnam(i))
    CASE ("TEMP")
      ! Convert to Kelvin from Fahrenheit:
      rawobs(1,1:numobs(i),i) = &
	(rawobs(1,1:numobs(i),i)-32.0)*5.0/9.0+temp_0
    CASE ("DEWP")
      ! Convert to Kelvin from Fahrenheit:
      rawobs(1,1:numobs(i),i) = &
	(rawobs(1,1:numobs(i),i)-32.0)*5.0/9.0+temp_0
    CASE ("WNDU")
      ! Convert to m/s from knots:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*knt2ms
    CASE ("WNDV")
      ! Convert to m/s from knots:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*knt2ms
    CASE ("VISB")
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*mile2m
    CASE ("REDP")
      ! Convert to pascal from mb:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*mb2pas
    END SELECT
  ENDDO

END SUBROUTINE LAPSUnit
