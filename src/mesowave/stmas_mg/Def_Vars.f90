!==========================================================
!  This header file defines all necessary variables.
!
!  NOTE: All global variables are six letter long.
!
!  HISTORY: 
!	Creation: YUANFU XIE	6-2005
!==========================================================

  !----------------
  ! LAPS variables:
  !----------------
  CHARACTER*4 :: varnam(MAXVAR)
  CHARACTER*256 :: dirstc

  INTEGER :: dirlen

  INTEGER :: lapsdt		! LAPS cycle time (secs)
  INTEGER :: mxstts		! Maximum number obs sites
  INTEGER :: i4time		! LAPS time frames
  INTEGER :: qc_val		! QC flag for threshold check
  INTEGER :: qc_std		! QC flag for standard dev-check
  INTEGER, ALLOCATABLE, DIMENSION(:) &
	  :: i4prev		! LAPStimes

  REAL*4 :: mising		! Laps bad data flag
  REAL   :: badsfc		! Laps bad sfc flag
  REAL   :: thresh(MAXVAR)	! Threshold for obs agaist bkg
  REAL, ALLOCATABLE, DIMENSION(:,:) &
	 :: latgrd,longrd	! Grid Lat/Lon
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
  INTEGER :: verbal		! Print message option
  INTEGER :: savdat		! Save background and obs
  INTEGER :: saveid		! Index for saving variable
  INTEGER :: numobs(MAXVAR)	! Maximum number of obs
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) &
	:: indice		! Interpolation indices

  REAL :: grdspc(3)		! Gridspacing in x, y, t
  REAL :: domain(2,3)		! Analysis domain
  REAL :: obsspc(3,MAXVAR)	! Observation spacing
  REAL*4, ALLOCATABLE, DIMENSION(:,:,:,:) &
       :: bkgrnd		! Background fields
  REAL, ALLOCATABLE, DIMENSION(:,:) &
       :: lndfac		! Land factor
  REAL, ALLOCATABLE, DIMENSION(:,:,:) &
       :: qc_obs		! QCed observations
  REAL, ALLOCATABLE, DIMENSION(:,:) &
       :: weight		! Observations weights
  REAL, ALLOCATABLE, DIMENSION(:,:,:) &
       :: coeffs		! Interpolation coefficients

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

  ! Namelists:
  NAMELIST /STMAS/numfic,numtmf,numgrd,numvar,savdat,saveid,verbal
  NAMELIST /STMAS/qc_val,qc_std
  NAMELIST /STMAS/maxitr,stmasi,penalt,stmasr
