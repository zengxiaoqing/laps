!dis
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis
!dis
!dis
!dis
!dis

SUBROUTINE LAPS2Mass(varlaps,imax,jmax,kmax,varmass,eta, &
                     plaps,psurf)

!==========================================================
!  This routine converts the LAPS grid function into GSI
!  grid over mass coordinate.
!
!  HISTORY: MAR. 2006 by YUANFU XIE.
!==========================================================

  ! This routine converts a LAPS variable (varlaps) on LAPS
  ! coordinate into one on a mass coordinate.

  IMPLICIT NONE
      
  INTEGER*4, INTENT(IN) :: imax,jmax,kmax  	! 3D array dimensions
  REAL*4, INTENT(IN) :: varlaps(imax,jmax,kmax)
  REAL*4, INTENT(IN) :: eta(kmax)		! eta=(p-pt)/(ps-pt)
  REAL*4, INTENT(IN) :: plaps(kmax)		! LAPS pressure level
  REAL*4, INTENT(IN) :: psurf(imax,jmax)	! Surface pressure
  REAL*4, INTENT(OUT) :: varmass(imax,jmax,kmax)

  ! Local variables:
  INTEGER :: i,j,k,l,mp
  REAL*4 :: a,b,p

  ! Conversion:
  DO j=1,jmax
    DO i=1,imax
      ! Use pressure decreasing property:
      mp = 1
      DO k=1,kmax
	! Mass pressure:
	p = eta(k)*(psurf(i,j)-plaps(kmax))+plaps(kmax)
	! Search the interval of LAPS pressure levels 
        ! containing the mass level:
	DO l=mp,kmax
	  IF (p .gt. plaps(l)) then
            mp = l
	    a = LOG(p/plaps(l))
	    b = LOG(plaps(l-1)/p)	! Log P interpolation
	    GOTO 10
	  ENDIF
	ENDDO
        ! Reach the top level:
	l = kmax
	a = 0.0
	b = 1.0

10	CONTINUE
	varmass(i,j,k) = (a*varlaps(i,j,l-1)+ &
     			  b*varlaps(i,j,l  ))/(a+b)
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE LAPS2Mass



SUBROUTINE dryairmass(dam,pdam,imax,jmax,kmax,pres_1d, &
     		      heights_3d,p_laps_bkg,t_laps_bkg)

!==========================================================
!  This routine computes the dry air mass and perturbation
!  using temperature field.
!
!  HISTORY: MAR. 2006 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER*4, INTENT(IN) :: imax,jmax,kmax  	! 3D array dimensions
  REAL*4, INTENT(OUT) :: dam(imax,jmax)   ! dry air mass in column
  REAL*4, INTENT(OUT) :: pdam(imax,jmax)  ! perturbation of dam
  REAL*4, INTENT(IN) :: pres_1d(kmax)    ! pressure values@levels
  REAL*4, INTENT(IN) :: heights_3d(imax,jmax,kmax)   ! heights
  REAL*4, INTENT(IN) :: p_laps_bkg(imax,jmax,kmax)   ! p bkg
  REAL*4, INTENT(IN) :: t_laps_bkg(imax,jmax,kmax)   ! t bkg

  ! Local variables:
  INTEGER :: i,j,k

  ! Use pressure equation: dp/p = - g/(RT) dz [pp 20 Holton]
  ! Integration between P_top and P_surface:
  ! int_z=h(surface)^h(top) d ln(p) = - g/R int_h(s)^h(t) dz/T.
  ! Thus: p(h(s)) = p(h(t))*exp(g/R*int_h(s)^h(t) dz/T, where
  ! p(h(t)) = pressure at the highest level: pres_1d(kmax).
  ! Integration is replaced by Riemann sum. 
  DO j=1,jmax
    DO i=1,imax
      ! Riemann sum:
      dam(i,j) = 0.0		! Initial summation
      pdam(i,j) = p_laps_bkg(i,j,1)	! Default value at bottom
      DO k=kmax,2,-1
        ! Summation til height=0
        IF (heights_3d(i,j,k-1) .gt. 0.0) then
	  ! Summation approximates integral:
          dam(i,j) = dam(i,j)+2.0/ &
     		    (t_laps_bkg(i,j,k)+t_laps_bkg(i,j,k-1))* &
     		    (heights_3d(i,j,k)-heights_3d(i,j,k-1))
	ELSE
          ! Interpolate the full pressure between +/- heights:
	  ! Note: the values of these heights are the weights
	  ! so the indices are swapped:
          pdam(i,j) = (p_laps_bkg(i,j,k-1)*heights_3d(i,j,k  )- &
     		       p_laps_bkg(i,j,k  )*heights_3d(i,j,k-1))/ &
     		       (heights_3d(i,j,k)-heights_3d(i,j,k-1))

	  ! Riemann sum over the partial grid:
	  ! Positive height * 0.5*(t_k + interpolated t at z =0):
          dam(i,j) = dam(i,j)+2.0*heights_3d(i,j,k  )/ &
     		(t_laps_bkg(i,j,k)+ &
     		(t_laps_bkg(i,j,k-1)*heights_3d(i,j,k  )- &
     	 	 t_laps_bkg(i,j,k  )*heights_3d(i,j,k-1))/ &
     		(heights_3d(i,j,k)-heights_3d(i,j,k-1)) )
	  GOTO 10
        ENDIF
      ENDDO
      ! Dry air mass: Integral of T -> dry surface pressure:
10    dam(i,j) = exp(9.806/287.0*dam(i,j))*pres_1d(kmax)
    ENDDO
  ENDDO

  ! Perturbation pressure:
  pdam(1:imax,1:jmax) = pdam(1:imax,1:jmax)-dam(1:imax,1:jmax)

END SUBROUTINE dryairmass

