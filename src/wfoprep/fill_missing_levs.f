SUBROUTINE fill_missing_levs(nx,ny,nz,plevs,data,missval,method)

! Subroutine to fill in missing values in a 3D data array on
! pressure levels.  This routine assumes you have valid data
! on the top most and bottom most levels.

  IMPLICIT NONE
  INTEGER, INTENT(IN)                 :: nx
  INTEGER, INTENT(IN)                 :: ny
  INTEGER, INTENT(IN)                 :: nz
  REAL, INTENT(IN)                    :: plevs(75)  ! Press in Pa
  REAL, INTENT(INOUT)                 :: data(nx,ny,nz)
  REAL, INTENT(IN)                    :: missval
  INTEGER, INTENT(IN)                 :: method

  INTEGER, PARAMETER                  :: METHOD_LIN = 1
  INTEGER, PARAMETER                  :: METHOD_LOG = 2

  INTEGER                             :: k,kbot,ktop
  REAL                                :: weight_top, weight_bot
  LOGICAL                             :: goodlev(nz) 

  goodlev(:) = .false.
  checklevs: DO k = 1,nz
    IF (MAXVAL(data(:,:,k)) .NE. missval) THEN
       goodlev(k) = .true.
    ENDIF
  ENDDO checklevs

  fixlevs: DO k = 1, nz
    IF (goodlev(k)) CYCLE fixlevs
    
    ! Else:

    ! Find first level below that is "good"
    find_kbot: DO kbot = k-1,1,-1
      IF (goodlev(kbot)) EXIT find_kbot
    ENDDO find_kbot

    ! Find first level above that is good
    find_ktop: DO ktop = k+1,nz
      IF (goodlev(ktop)) EXIT find_ktop
    ENDDO find_ktop

    ! Compute interp weights using plevs and method
    IF (method .EQ. METHOD_LIN) THEN   
      ! Linear interp in Pressure
      weight_bot = (plevs(k)-plevs(ktop)) / &
                   (plevs(kbot)-plevs(ktop))
      weight_top = 1. - weight_bot
    ELSE IF (method .EQ. METHOD_LOG) THEN
      ! Linear in ln(P)
      weight_bot = ( ALOG(plevs(k)/plevs(ktop)) ) / &
                   ( ALOG(plevs(kbot)/plevs(ktop)) )
      weight_top = 1. - weight_bot  
    ENDIF
    data(:,:,k) = weight_bot*data(:,:,kbot) + weight_top*data(:,:,ktop)
  ENDDO fixlevs
  RETURN
END SUBROUTINE fill_missing_levs
