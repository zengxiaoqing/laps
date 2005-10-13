!==========================================================
!  Define parameters.
!==========================================================

  INTEGER, PARAMETER :: MAXVAR=11
  INTEGER, PARAMETER :: LSXVAR=21

  !****************
  ! LAPS Constants:
  !****************

  ! Lapse rates: see LAPS mdatlap.f under sfc:
  ! 1: temperature; 2: Dewpoint
  REAL, PARAMETER :: lapses(2) = (/-0.01167, -0.007/)

  ! Unit conversions:
  REAL, PARAMETER :: temp_0 = 273.16		! Absolute 0
  REAL, PARAMETER :: knt2ms = 0.5277777778	! Knot 2 m/s
  REAL, PARAMETER :: mile2m = 1609.0		! Mile to meter
  REAL, PARAMETER :: mb2pas = 100.0		! Mb 2 pascal
  REAL, PARAMETER :: gascnt = 287.0		! Gas constant 
						! for dry air
  REAL, PARAMETER :: spheat = 1004.0		! Specific heat 
						! at constant pres

  !****************************
  ! Machine related parameters:
  !****************************

  REAL, PARAMETER :: epsiln = 1.0e-18
