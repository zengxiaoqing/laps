!dis    Forecast Systems Laboratory
!dis    NOAA/OAR/ERL/FSL
!dis    325 Broadway
!dis    Boulder, CO     80303
!dis
!dis    Forecast Research Division
!dis    Local Analysis and Prediction Branch
!dis    LAPS
!dis
!dis    This software and its documentation are in the public domain and
!dis    are furnished "as is."  The United States government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  They assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis     technical support to users.
!dis
!dis    Permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  All modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  If significant modifications or enhancements
!dis    are made to this software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

MODULE Definition

!==========================================================
!  This module defines all parameters and variables needed
!  in STMAS multigrid analysis.
!
!  HISTORY:
!	Creation: YUANFU XIE 8-2005
!==========================================================

!==========================================================
!  Define parameters.
!==========================================================

  INTEGER, PARAMETER :: MAXVAR=15
  INTEGER, PARAMETER :: LSXVAR=21

  !****************
  ! LAPS Constants:
  !****************

  ! Lapse rates: see LAPS mdatlaps.f under sfc:
  ! 1: temperature; 2: Dewpoint
  REAL, PARAMETER :: lapses(2) = (/-0.01167, -0.007/)

  ! Unit conversions:
  REAL, PARAMETER :: temp_0 = 273.16		! Absolute 0
  REAL, PARAMETER :: knt2ms = 0.51444444	! Knot 2 m/s
  REAL, PARAMETER :: mile2m = 1609.0		! Mile to meter
  REAL, PARAMETER :: mb2pas = 100.0		! Mb 2 pascal
  REAL, PARAMETER :: inch2m = 0.0254		! Inch to meter
  REAL, PARAMETER :: gascnt = 287.0		! Gas constant 
						! for dry air
  REAL, PARAMETER :: spheat = 1004.0		! Specific heat 
						! at constant pres

  !****************************
  ! Machine related parameters:
  !****************************

  REAL, PARAMETER :: epsiln = 1.0e-18

!==========================================================
!  This header file defines all necessary variables.
!
!  NOTE: All global variables are six letter long.
!
!  HISTORY: 
!	Creation: YUANFU XIE	6-2005
!       Modified: YUANFU XIE    12-2008 for adding diagnl.
!==========================================================

  !----------------
  ! LAPS variables:
  !----------------
  CHARACTER*4 :: varnam(MAXVAR)
  CHARACTER*256 :: dirstc

  character*8 :: date0,date1
  character*10 :: time0,time1
  character*5 :: zone0,zone1
  integer :: timing0(8),timing1(8) 

  INTEGER :: dirlen

  INTEGER :: lapsdt		! LAPS cycle time (secs)
  INTEGER :: mxstts		! Maximum number obs sites
  INTEGER :: i4time		! LAPS analysis time 
  INTEGER :: i4wndw(2)		! Analysis time window
  INTEGER :: qc_val		! QC flag for threshold check
  INTEGER :: qc_std		! QC flag for standard dev-check
  INTEGER, ALLOCATABLE, DIMENSION(:) &
	  :: i4prev		! LAPStimes

  REAL   :: mising		! Laps bad data flag
  REAL   :: badsfc		! Laps bad sfc flag
  REAL   :: thresh(MAXVAR)	! Threshold for obs agaist bkg
  REAL   :: pnlt_v(MAXVAR)	! Penalty of each variable !added by min-ken hsieh, used in STMASAna
  REAL   :: rdplvl		! Reduced Pressure level definted in LAPS
  REAL, ALLOCATABLE, DIMENSION(:,:) &
	 :: latgrd,longrd,topogr! Grid Lat/Lon/topo
  REAL, ALLOCATABLE, DIMENSION(:,:,:) &
	 :: rawobs		! Raw observations
  REAL, ALLOCATABLE, DIMENSION(:,:) &
	 :: bkgobs		! Background at obs sites

  !-----------------------------------------------------
  ! Interactive variables between LAPS and minimization:
  !-----------------------------------------------------
  INTEGER :: numfic(3)		! Number fictitious points
  INTEGER :: numgrd(3)		! Analysis grid numbers
  INTEGER :: numtmf		! Number of LAPS bkgd frames
  INTEGER :: numvar		! Number of Analysis vars
  INTEGER :: needbk(MAXVAR)	! Add background to incements
  INTEGER :: bounds(MAXVAR)	! Bound constraints
  INTEGER :: radius(MAXVAR)     ! Obs cover grids ! added by min-ken hsieh, used in AddBkgrd
  INTEGER :: lndsea(MAXVAR)     ! Land/sea process! added by min-ken hsieh, used in AddBkgrd
  INTEGER :: slevel(MAXVAR)     ! Starting level of analysis! added by min-ken hsieh, used in STMASAna
  INTEGER :: verbal		! Print message option
  INTEGER :: press_pert		! 1: compute pressure perturbation
  INTEGER :: savdat		! Save background and obs
  INTEGER :: saveid		! Index for saving variable
  INTEGER :: numobs(MAXVAR)	! Maximum number of obs
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) &
	:: indice		! Interpolation indices

  REAL :: grdspc(3)		! Gridspacing in x, y, t
  REAL :: domain(2,3)		! Analysis domain
  REAL :: obsspc(3,MAXVAR)	! Observation spacing
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) &
       :: bkgrnd		! Background fields
  REAL, ALLOCATABLE, DIMENSION(:,:) &
       :: lndfac		! Land factor
  REAL, ALLOCATABLE, DIMENSION(:,:,:) &
       :: qc_obs		! QCed observations
  REAL, ALLOCATABLE, DIMENSION(:,:) &
       :: weight		! Observations weights
  REAL, ALLOCATABLE, DIMENSION(:,:,:) &
       :: coeffs		! Interpolation coefficients
  REAL, ALLOCATABLE, DIMENSION(:,:,:) &
       :: diagnl		! Diagonal array of B for J_b term
  LOGICAL, ALLOCATABLE, DIMENSION(:,:,:,:) &
       :: uncovr		! Uncovered grids !added by min-ken hsieh, used in AddBkgrd and STMASAna

  !------------------------
  ! Minimization variables:
  !------------------------

  INTEGER :: stmasi(1)		! STMAS integer parameters
  INTEGER :: maxitr		! Maximum iterations

  REAL :: phydxy(2)		! Physical spacing for derivatives
  REAL :: stmasr(1)		! STMAS real parameters
  REAL :: penalt		! Penalty parameter
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) &
	:: analys		! Analyzed fields
  !------------------------
  !Verification variables:
  !                             !modified by min-ken hsieh
  !------------------------
  CHARACTER*20,ALLOCATABLE, DIMENSION(:,:) &
	 :: stanam 		!store stn name for each variable

  ! Namelists:
  NAMELIST /STMAS/numfic,numtmf,numgrd,numvar,savdat,saveid,verbal,press_pert
  NAMELIST /STMAS/qc_val,qc_std
  NAMELIST /STMAS/maxitr,stmasi,penalt,stmasr

CONTAINS			! Common utility routines

SUBROUTINE Grid2Obs(indx,coef,obsv,nobs,wght,stna,ngrd,dxyt,domn)

!==========================================================
!  This routine finds the indices and coefficients for an
!  interpolation scheme from a grid function to observation
!  sites.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!       Modification:
!                 25-08-2008 by min-ken hsieh
!                 add parameter stna to map each obs its stn name for STMASVer
!
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3)
  REAL, INTENT(IN) :: dxyt(3),domn(2,3)

  INTEGER, INTENT(OUT) :: indx(6,nobs)
  INTEGER, INTENT(INOUT) :: nobs
  REAL, INTENT(OUT) :: coef(6,nobs)
  REAL, INTENT(INOUT) :: obsv(4,nobs),wght(nobs)
  CHARACTER*20, INTENT(INOUT) :: stna(nobs)		!by min-ken hsieh

  ! Local variables:
  INTEGER :: i,ier
  INTEGER :: nib		! Number of obs in box

  ! Count obs in box:
  nib = 0
  ! Interpolation for each obs:
  DO i=1,nobs
    CALL Intplt3D(obsv(2:4,i),ngrd,dxyt,domn, &
		  indx(1,i),coef(1,i),ier)

    ! Check:
    IF (ier .EQ. 0) THEN
      nib = nib+1

      ! Save the obs and its weight:
      obsv(1:4,nib) = obsv(1:4,i)
      wght(nib) = wght(i)
      stna(nib) = stna(i)

      indx(1:6,nib) = indx(1:6,i)
      coef(1:6,nib) = coef(1:6,i)

    ELSE
	if (verbal .eq. 1) &
         print*,'Grid2Obs: Obs out of the analysi domain ',obsv(2:4,i),i,domn

    ENDIF
  ENDDO

  ! Count inbox obs:
  nobs = nib

END SUBROUTINE Grid2Obs

SUBROUTINE Intplt3D(pstn,ngrd,gspc,domn,indx,coef,ierr)

!==========================================================
!  This routine returns interpolation coefficients and
!  indices of a given location from a grid in 3D space.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3)	! Numbers gridpoint
  REAL,    INTENT(IN) :: pstn(3)	! Interpolate point
  REAL,    INTENT(IN) :: gspc(3)	! Grid spacing
  REAL,    INTENT(IN) :: domn(2,3)	! Domain

  INTEGER, INTENT(OUT) :: indx(6)	! Indices
  INTEGER, INTENT(OUT) :: ierr		! Error flag
  REAL,    INTENT(OUT) :: coef(6)	! Coefficients

  ! Local variables:
  INTEGER :: i

  ierr = 0
  ! Check if it is in box?
  DO i=1,3
     IF (pstn(i) .LT. domn(1,i)) ierr = -i
     IF (pstn(i) .GT. domn(2,i)) ierr = i
  ENDDO

  ! Indices:
  indx(1:3) = (pstn-domn(1,1:3))/gspc

  ! Coefficients:
  coef(1:3) = (pstn-indx(1:3)*gspc-domn(1,1:3))/gspc

  indx(1:3) = indx(1:3)+1

  indx(4:6) = MIN(indx(1:3)+1,ngrd)
  coef(4:6) = coef(1:3)
  coef(1:3) = 1.0-coef(4:6)

END SUBROUTINE Intplt3D

END MODULE Definition
