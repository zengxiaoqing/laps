!------------------------------------------------------------------------------
!M+
! NAME:
!       parameters
!
! PURPOSE:
!       Module to hold RT model parameter constants
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE parameters
!
! OUTPUTS:
!       Parameters
!       ----------
!       MAX_N_ABSORBERS:             INTEGER parameter defining the maximum
!                                    number of absorbing species.
!
!       MAX_N_PREDICTORS_USED:       INTEGER parameter defining the maximum
!                                    number of predictors used in the absorption
!                                    coefficient calculation.
!
!       MAX_N_STANDARD_PREDICTORS:   INTEGER parameter defining the number of
!                                    standard (i.e. absorber independent) predictors.
!
!       MAX_N_INTEGRATED_PREDICTORS: INTEGER parameter defining the number of
!                                    integrated (i.e. absorber dependent) predictors.
!
!       MAX_N_PREDICTORS:            INTEGER parameter defining the total number
!                                    of predictors for all absorbers,
!                                      MAX_N_STANDARD_PREDICTORS + &
!                                      ( MAX_N_ABSORBERS * MAX_N_INTEGRATED_PREDICTORS )
!
!       MAX_N_LAYERS:                INTEGER parameter defining the maximum
!                                    number of atmospheric layers allowed for
!                                    input.
!
!       MAX_N_ABSORBER_LAYERS:       INTEGER parameter defining the maximum
!                                    number of
!                                    absorber space layers.
!
!       The following are defined numerical parameters. Parameter definitions
!       for non-integer constants are used throughout the code to facilitate
!       changes to the default floating point precision if required:
!
!       ZERO      = 0.0
!       ONE       = 1.0
!       TWO       = 2.0
!       THREE     = 3.0
!       POINT_5   = 0.5
!       POINT_75  = 0.75
!
!       PI                 = 3.14159265
!       DEGREES_TO_RADIANS = PI / 180.0
!
!       TOA_PRESSURE = 0.005
!
!
!       Pseudo-parameters
!       -----------------
!
!       These values are not parameters in the Fortran sense in that they are
!       defined at run-time based on user inputs but once defined, they are
!       (or should be) invariant.
!
!       MAX_N_CHANNELS:              INTEGER defining the maximum number of
!                                    instrument channels. This defines the
!                                    valid channels for the valid satellites
!                                    that the user has selected.
!                                    Upon RTM initialisation and destruction
!                                    the value is set to -1.
!                                    The value of MAX_N_CHANNELS can only be
!                                    accessed through its own methods:
!                                      - set_max_n_channels()
!                                      - reset_max_n_channels()
!                                      - get_max_n_channels()
!       
!
! MODULES:
!       None.
!
! CONTAINS:
!       None.
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       None.
!
! RESTRICTIONS:
!       None.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS@NOAA/NCEP 31-Jul-2000
!                       pvandelst@ncep.noaa.gov
!
!
!  Copyright (C) 2000 Paul van Delst
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!M-
!------------------------------------------------------------------------------

MODULE parameters

  ! ---------------------
  ! Module use statements
  ! ---------------------

  USE type_kinds, ONLY : fp_kind


  ! ------------------
  ! Default visibility
  ! ------------------

  PRIVATE


  ! ---------------------------
  ! Current number of absorbers
  ! ---------------------------

  INTEGER, PUBLIC, PARAMETER :: MAX_N_ABSORBERS = 3


  ! --------------------
  ! Number of predictors.
  ! --------------------

  INTEGER, PUBLIC, PARAMETER :: MAX_N_PREDICTORS_USED = 5

  INTEGER, PUBLIC, PARAMETER :: MAX_N_STANDARD_PREDICTORS   = 9
  INTEGER, PUBLIC, PARAMETER :: MAX_N_INTEGRATED_PREDICTORS = 6

  INTEGER, PUBLIC, PARAMETER :: MAX_N_PREDICTORS = MAX_N_STANDARD_PREDICTORS + &
                                                   ( MAX_N_ABSORBERS * MAX_N_INTEGRATED_PREDICTORS )



  ! --------------------------------------
  ! Number of absorber layers in algorithm
  ! --------------------------------------

  INTEGER, PUBLIC, PARAMETER :: MAX_N_ABSORBER_LAYERS = 300


  ! ----------------------------------------------------------
  ! Number of channels (for ALL satellites - really the number
  ! of satellites x number of channels USED per satellite)
  !
  ! This is also the number of lines in the satellite 
  ! information file.
  !
  ! Eventually the MAX_N_PROFILES and MAX_N_LAYERS values
  ! will be dynamic, i.e they will be defined by user inputs.
  ! For now, however, they're hardwired.
  ! --------------------------------------------------------

  ! -- Accessed via SET_<name>, RESET_<name>, and GET_<name> routines
  INTEGER, PRIVATE, PARAMETER :: RESET_VALUE = -1
  INTEGER, PRIVATE, SAVE      :: MAX_N_CHANNELS = RESET_VALUE

  INTEGER, PUBLIC, PARAMETER :: MAX_N_PROFILES = 128
  INTEGER, PUBLIC, PARAMETER :: MAX_N_LAYERS   = 100
!!!  INTEGER, PRIVATE :: MAX_N_PROFILES = RESET_VALUE
!!!  INTEGER, PRIVATE :: MAX_N_LAYERS   = RESET_VALUE


  ! -----
  ! Flags
  ! -----

  ! -- Direction flags for transmittance calculation
  INTEGER, PUBLIC, PARAMETER :: DOWN = 0
  INTEGER, PUBLIC, PARAMETER :: UP   = 1


  ! --------------------
  ! Numerical parameters
  ! --------------------

  ! -- Numbers
  REAL( fp_kind ), PUBLIC, PARAMETER :: ZERO      = 0.0_fp_kind
  REAL( fp_kind ), PUBLIC, PARAMETER :: ONE       = 1.0_fp_kind
  REAL( fp_kind ), PUBLIC, PARAMETER :: TWO       = 2.0_fp_kind
  REAL( fp_kind ), PUBLIC, PARAMETER :: THREE     = 3.0_fp_kind
  REAL( fp_kind ), PUBLIC, PARAMETER :: POINT_5   = 0.5_fp_kind
  REAL( fp_kind ), PUBLIC, PARAMETER :: POINT_75  = 0.75_fp_kind

  ! -- Precision/tolerance
  REAL( fp_kind ), PUBLIC, PARAMETER :: TOLERANCE = EPSILON( ONE )

  ! -- Constant to allow degrees->radians conversion
  REAL( fp_kind ), PUBLIC, PARAMETER :: PI = 3.14159265358979323_fp_kind
  REAL( fp_kind ), PUBLIC, PARAMETER :: DEGREES_TO_RADIANS = PI / 180.0_fp_kind

  ! -- Top-Of-Atmosphere pressure in hPa
  REAL( fp_kind ), PUBLIC, PARAMETER :: TOA_PRESSURE = 0.005_fp_kind

  ! -- Reciprocal gravity (scaled by 100 for use with pressure in hPa)
  REAL( fp_kind ), PUBLIC, PARAMETER :: RECIPROCAL_GRAVITY = ONE / 980.665_fp_kind

  ! -- Diffusivity angle secant = ACOS( 3/5 ) in degrees (~53.13)
  REAL( fp_kind ), PUBLIC, PARAMETER :: SECANT_DIFFUSIVITY_ANGLE = 5.0_fp_kind / THREE

  ! -- Maximum solar zenith angle secant definition. Should be determined
  ! -- by the maximum angle secant used in generating the transmittance
  ! -- model coefficients, i.e. a secant of 2.25 => 63.6deg. Users have
  ! -- requested the value be 85deg => secant of ~11.47.
  REAL( fp_kind ), PUBLIC, PARAMETER :: MAX_SOLAR_ANGLE = 85.0_fp_kind
  REAL( fp_kind ), PUBLIC, PARAMETER :: MAX_SECANT_SOLAR_ANGLE = 11.473711738554476_fp_kind

!!!! The following is preferred but not allowed in Fortran 90          !!!!
!!!! Use of non-integer result intrinsics allowed in Fortran 95 though !!!!
!  REAL( fp_kind ), PRIVATE, PARAMETER :: MAX_SOLAR_ANGLE = 85.0_fp_kind
!  REAL( fp_kind ), PUBLIC,  PARAMETER :: MAX_SECANT_SOLAR_ANGLE = ONE / COS( DEGREES_TO_RADIANS * MAX_SOLAR_ANGLE )
                                                                 

  ! ---------------------
  ! Subprogram visibility
  ! ---------------------

  PUBLIC :: set_max_n_channels, reset_max_n_channels, get_max_n_channels


CONTAINS


  ! -------------------------------------------
  ! Subroutines to SET and GET the value of the 
  ! "pseudo-parameter" MAX_N_CHANNELS
  ! -------------------------------------------

  ! -- Set the value
  SUBROUTINE set_max_n_channels( value )
    INTEGER, INTENT( IN ) :: value
    MAX_N_CHANNELS = value
  END SUBROUTINE set_max_n_channels

  ! -- REset the value
  SUBROUTINE reset_max_n_channels()
    MAX_N_CHANNELS = RESET_VALUE
  END SUBROUTINE reset_max_n_channels

  ! -- Get the value and test if it's been set
  SUBROUTINE get_max_n_channels( value, is_set )
    INTEGER, INTENT( OUT )           :: value
    LOGICAL, INTENT( OUT ), OPTIONAL :: is_set
    value = MAX_N_CHANNELS
    IF ( PRESENT( is_set ) ) THEN
      is_set = .FALSE.
      IF ( value /= RESET_VALUE ) is_set = .TRUE.
    END IF
  END SUBROUTINE get_max_n_channels

END MODULE parameters


!-------------------------------------------------------------------------------
!                          -- MODIFICATION HISTORY --
!-------------------------------------------------------------------------------
!
! $Id$
!
! $Date$
!
! $Revision$
!
! $State$
!
! $Log$
! Revision 1.8  2001/08/31 21:09:45  paulv
! - Added the secant of the maximum solar angle as a parameter. The secant
!   value was calculated and expressed as an explicit number since Fortran 90
!   does not allow intrinsics that have other than integer results in a
!   parameter initialisation expression.
! - Changed  MAX_SOLAR_ZENITH_ANGLE name to MAX_SOLAR_ANGLE.
!
! Revision 1.7  2001/08/16 16:44:14  paulv
! - Updated documentation.
! - Changed MAX_N_CHANNELS attributes from PUBLIC to PRIVATE, SAVE. The value
!   of MAX_N_CHANNELS is now accessed via its public methods:
!     set_max_n_channels
!     reset_max_n_channels()
!    get_max_n_channels()
! - Added RESET_VALUE parameter.
! - Removed POINT_333 parameter.
!
! Revision 1.6  2001/07/12 16:46:12  paulv
! - Added PRIVATE statement to prevent definitions in module TYPE_KINDS
!   being available outside the scope of this module.
!
! Revision 1.5  2001/05/29 17:42:55  paulv
! - Now use TYPE_KINDS module parameter FP_KIND to set the floating point
!   data type. All REAL declarations are now typed with FP_KIND.
! - Added direction flags for transmittance calculation.
!
! Revision 1.4  2000/11/09 20:32:11  paulv
! - Removed MAX_N_CHANNELS as a parameter. It is now a "pseudo" parameter
!   in that it is determined by the number of channels for which coefficients
!   are defined.
! - Added some more numerical parameters.
!
! Revision 1.3  2000/08/31 19:36:33  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 1.2  2000/08/24 15:39:45  paulv
! - Changed the parameter name that references how many predictors of the
!   total set to use from MAX_N_PREDICTORS_TO_USE to MAX_N_PREDICTORS_USED.
!   I felt this would clarify (for me at least) that while the maximum
!   number of predictors is set, the number that is actually used can be
!   less than that.
! - Current maximum number of layers is 100. This is a temporary limit for
!   testing purposes.
! - The parameter RECIPROCAL_GRAVITY was removed from this module and placed
!   in the ABSORBER_PROFILE module where it is used.
!
! Revision 1.1  2000/08/08 16:57:21  paulv
! Initial checkin
!
!
!
!
