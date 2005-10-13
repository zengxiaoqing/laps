SUBROUTINE LAPSConf

!==========================================================
!  This routine reads in necessary configuration grid data.
!
!  HISTORY: 
!	Creation: 8-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  ! Local variables:
  CHARACTER*9 :: fnm,ext,var
  CHARACTER :: unt*60,com*60
  INTEGER :: i,j,err	! Error indicator

  ! Physical grid points and spacing info:
  ext = 'nest7grid'
  var = 'LAT'
  CALL RD_LAPS_STATIC(dirstc,ext,numgrd(1),numgrd(2),1, &
		      var,unt,com,latgrd,phydxy(1),err)
  IF (err .NE. 1) PRINT*,'STMAS>LAPSInfo: Error getting LAT'
  var = 'LON'
  CALL RD_LAPS_STATIC(dirstc,ext,numgrd(1),numgrd(2),1, &
		      var,unt,com,longrd,phydxy(2),err)
  IF (err .NE. 1) PRINT*,'STMAS>LAPSInfo: Error getting LON'
  var = 'LDF'
  CALL RD_LAPS_STATIC(dirstc,ext,numgrd(1),numgrd(2),1, &
		      var,unt,com,lndfac,grdspc(2),err)
  IF (err .NE. 1) PRINT*,'STMAS>LAPSInfo: Error getting LON'
  ! Removing meaningless land factors:
  DO j=1,numgrd(2)
    DO i=1,numgrd(1)
      IF (lndfac(i,j) .GT. 1.0) lndfac(i,j) = 1.0
      IF (lndfac(i,j) .LT. 0.0) lndfac(i,j) = 0.0
    ENDDO
  ENDDO

  ! time step:
  grdspc(3) = lapsdt

  !*********************
  ! STMAS configuration:
  !*********************

  ! Analysis grid numbers:
  numgrd = numgrd+2*numfic

  ! Analysis domain:
  domain(1,1:2) = 1.0-numfic(1:2)	! X/Y: use grid numbers
  domain(2,1:2) = FLOAT(numgrd(1:2)-numfic(1:2))
  domain(1,3) = MOD(i4time-(numtmf-1)*lapsdt,86400)
  domain(2,3) = domain(1,3)+(numtmf-1)*lapsdt
  grdspc(1:2) = 1.0			! Based the domain setting

END SUBROUTINE LAPSConf
