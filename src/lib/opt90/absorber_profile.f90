!------------------------------------------------------------------------------
!M+
! NAME:
!       absorber_profile
!
! PURPOSE:
!       Module containing routines to compute and assemble the integrated
!       absorber profiles.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE absorber_profile
!
! OUTPUTS:
!       None.
!
! MODULES:
!       parameters:  Module containing parameter definitions for the
!                    RT model.
!
! CONTAINS:
!       compute_absorber_amount:     PUBLIC subroutine to compute the integrated
!                                    absorber profiles. Currently the absorbers
!                                    are:
!                                      - Water vapor
!                                      - Dry/fixed gases
!                                      - Ozone
!
!       compute_absorber_amount_TL:  PUBLIC subroutine to compute the tangent-
!                                    linear form of the integrated absorber 
!                                    profiles.
!
!       compute_absorber_amount_AD:  PUBLIC subroutine to compute the adjoint of
!                                    the integrated absorber profiles.
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
! COMMENTS:
!       All of the array documentation lists the dimensions by a single letter.
!       Throughout the RTM code these are:
!         I: Array dimension is of I predictors (Istd and Iint are variants).
!         J: Array dimension is of J absorbing species.
!         K: Array dimension is of K atmospheric layers.
!         L: Array dimension is of L spectral channels.
!         M: Array dimension is of M profiles.
!       Not all of these dimensions will appear in every module.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS@NOAA/NCEP 01-Aug-2000
!                       pvandelst@ncep.noaa.gov
!
!       Adapted from code written by: Thomas J.Kleespies
!                                     NOAA/NESDIS/ORA
!                                     tkleespies@nesdis.noaa.gov
!
!  Copyright (C) 2000 Thomas Kleespies, Paul van Delst
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

MODULE absorber_profile


  ! ---------------------
  ! Module use statements
  ! ---------------------

  USE type_kinds, ONLY : fp_kind
  USE parameters
  USE transmittance_coefficients, ONLY : alpha, &
                                         absorber_space_levels


  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------

  IMPLICIT NONE


  ! ------------------
  ! Default visibility
  ! ------------------

  PRIVATE


  ! ----------------------------------
  ! Explicit visibility of subprograms
  ! ----------------------------------

  PUBLIC :: compute_absorber_amount
  PUBLIC :: compute_absorber_amount_TL
  PUBLIC :: compute_absorber_amount_AD
  PUBLIC :: find_absorber_layer_index


CONTAINS



!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_absorber_amount
!
! PURPOSE:
!       PUBLIC subroutine to compute the integrated profiles for all the
!       absorbers. Currently the number of absorbers are:
!         - Water vapor
!         - Dry/fixed gases (pressure == absorber amount)
!         - Ozone
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_absorber_amount( pressure,    &  ! Input,  K
!                                     water_vapor, &  ! Input,  K
!                                     ozone,       &  ! Input,  K
!                                     absorber     )  ! Output, 0:K x J
!
! INPUT ARGUMENTS:
!       pressure:     Profile LEVEL pressure array.
!                     UNITS:      hPa
!                     TYPE:       REAL( fp_kind )
!                     DIMENSION:  K
!                     ATTRIBUTES: INTENT( IN )
!
!       water_vapor:  Profile LAYER water vapor mixing ratio array.
!                     UNITS:      g/kg
!                     TYPE:       REAL( fp_kind )
!                     DIMENSION:  K
!                     ATTRIBUTES: INTENT( IN )
!
!       ozone:        Profile LAYER ozone mixing ratio array.
!                     UNITS:      ppmv
!                     TYPE:       REAL( fp_kind )
!                     DIMENSION:  K
!                     ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       absorber:     Profile LEVEL integrated absorber amount array.
!                     UNITS:      Varies with absorber
!                     TYPE:       REAL( fp_kind )
!                     DIMENSION:  0:K x J
!                     ATTRIBUTES: INTENT( OUT )
!
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! CALLS:
!       None.
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       None known.
!
! RESTRICTIONS:
!       None.
!
! PROCEDURE:
!       The function calculates and accumulates the integrated path length
!       through k atmospheric layers,
!
!                     __ k
!                 1  \
!         u(k) = ---  >  q(i).dp(i)
!                 g  /__
!                       i=1
!
!       with the units of u dependent on those of the input mixing ratio, q.
!
!       The exception is the dry gas absorber amount which is represented
!       by the pressure (which is already an integrated quantity).
!
!       Also,
!
!         u(0) = 0.0          for water vapor and ozone
!              = TOA_PRESSURE for dry/fixed gases
!
!       The routine loop over layers, k, with each absorber amount calculated
!       independently (like looping over absorber, j). This is opposite
!       to what is recommended (since the arrays are dimensioned [k,j]) but 
!       later use of the arrays require the [k,j] ordering where the absorber
!       loop over j is the OUTSIDE loop (see COMPUTE_TRANSMITTANCE()).
!S-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_absorber_amount( pressure,    &  ! Input,  K
                                      water_vapor, &  ! Input,  K
                                      ozone,       &  ! Input,  K
                                      absorber     )  ! Output, 0:K x J


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: pressure     ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: water_vapor  ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: ozone        ! Input,  K

    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( OUT ) :: absorber     ! Output, 0:K x J


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: k
    REAL( fp_kind ) :: dp



    !#--------------------------------------------------------------------------#
    !#              -- Initialise 0'th LEVEL absorber amounts --                #
    !#                                                                          #
    !# This is done so that layer differences and averages can be calculated    #
    !# simply in the predictor and transmittance routines.                      #
    !#--------------------------------------------------------------------------#

    absorber( 0, : ) = ZERO



    !#--------------------------------------------------------------------------#
    !#             -- Assemble the NADIR LEVEL absorber profiles --             #
    !#--------------------------------------------------------------------------#

    ! ---------------------------------------------
    ! TOA layer. Pressure is dimensioned as K hence
    ! the top input layer is treated separately
    ! ---------------------------------------------

    dp = pressure( 1 )

    ! -- Integrated absorber amount
    absorber( 1, 1 ) = RECIPROCAL_GRAVITY * dp * water_vapor( 1 )
    absorber( 1, 2 ) = dp
    absorber( 1, 3 ) = RECIPROCAL_GRAVITY * dp * ozone( 1 )


    ! --------------------------------
    ! Loop over layers, TOA - 1 -> SFC
    ! --------------------------------

    k_layer_loop: DO k = 2, SIZE( pressure )

      ! -- Layer pressure difference
      dp = pressure( k ) - pressure( k-1 )

      ! -- Integrated absorber amounts
      absorber( k, 1 ) = absorber( k-1, 1 ) + ( RECIPROCAL_GRAVITY * dp * water_vapor( k ) )
      absorber( k, 2 ) = pressure( k )
      absorber( k, 3 ) = absorber( k-1, 3 ) + ( RECIPROCAL_GRAVITY * dp * ozone( k ) )

    END DO k_layer_loop

  END SUBROUTINE compute_absorber_amount





!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_absorber_amount_TL
!
! PURPOSE:
!       PUBLIC subroutine to compute the tangent linear form of the integrated
!       profiles for all the absorbers. Currently the number of absorbers
!       are:
!         - Water vapor
!         - Dry/fixed gases (pressure == absorber amount)
!         - Ozone
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_absorber_amount_TL( &
!                                        ! -- Forward input
!                                        pressure,       &  ! Input, K
!                                        water_vapor,    &  ! Input, K
!                                        ozone,          &  ! Input, K
!
!                                        ! -- Tangent-linear input
!                                        pressure_TL,    &  ! Input, K
!                                        water_vapor_TL, &  ! Input, K
!                                        ozone_TL,       &  ! Input, K
!
!                                        ! -- Tangent-linear output
!                                        absorber_TL     )  ! Output, 0:K x J
!
! INPUT ARGUMENTS:
!       pressure:        Profile LEVEL pressure array.
!                        UNITS:      hPa
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN )
!
!       water_vapor:     Profile LAYER water vapor mixing ratio array.
!                        UNITS:      g/kg
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN )
!
!       ozone:           Profile LAYER ozone mixing ratio array.
!                        UNITS:      ppmv
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN )
!
!       pressure_TL:     Profile LEVEL tangent-linear pressure array,
!                        i.e. the pressure perturbation.
!                        UNITS:      hPa
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K, number of levels - 1
!                        ATTRIBUTES: INTENT( IN )
!
!       water_vapor_TL:  Profile LAYER tangent-linear water vapor mixing
!                        ratio array, i.e. the water vapor mixing
!                        ratio perturbation.
!                        UNITS:      g/kg
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN )
!
!       ozone_TL:        Profile LAYER tangent-linear ozone mixing ratio
!                        array i.e. the ozone mixing ratio perturbation.
!                        UNITS:      ppmv
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       absorber_TL:     Profile LEVEL tangent-linear average integrated 
!                        absorber amount array.
!                        UNITS:      Varies with absorber.
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  0:K x J
!                        ATTRIBUTES: INTENT( OUT )
!
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! CALLS:
!       None
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       None known.
!
! RESTRICTIONS:
!       None.
!
! PROCEDURE:
!       The function calculates and accumulates the tangent-linear of the
!       integrated path length through k atmospheric layers,
!                        __ k
!                    1  \
!         u_TL(k) = ---  >  [ q_TL(i).dp(i) + q(i).dp_TL(i) ]
!                    g  /__
!                          i=1
!
!       with the units of u dependent on those of the input mixing ratio, q.
!
!       The exception is the dry gas absorber amount which is represented
!       by the tangent linear of the pressure (i.e. a perturbation).
!
!       Also,
!
!         u(0)    = 0.0          for water vapor and ozone
!                 = TOA_PRESSURE for dry/fixed gases
!
!         u_TL(0) = 0.0 for all absorbers
!
!       The routine loop over layers, k, with each absorber amount calculated
!       independently (like looping over absorber, j). This is opposite
!       to what is recommended (since the arrays are dimensioned [k,j]) but 
!       later use of the arrays require the [k,j] ordering where the absorber
!       loop over j is the OUTSIDE loop (see COMPUTE_TRANSMITTANCE_TL()).
!S-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_absorber_amount_TL( &
                                         ! -- Forward input
                                         pressure,       &  ! Input,  K
                                         water_vapor,    &  ! Input,  K
                                         ozone,          &  ! Input,  K

                                         ! -- Tangent-linear input
                                         pressure_TL,    &  ! Input,  K
                                         water_vapor_TL, &  ! Input,  K
                                         ozone_TL,       &  ! Input,  K

                                         ! -- Tangent-linear output
                                         absorber_TL     )  ! Output, 0:K x J


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward input
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: pressure        ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: water_vapor     ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: ozone           ! Input,  K

    ! -- Tangent-linear input
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: pressure_TL     ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: water_vapor_TL  ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: ozone_TL        ! Input,  K

    ! -- Tangent-linear output
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( OUT ) :: absorber_TL     ! Output, 0:K x J



    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: k

    REAL( fp_kind ) :: dp
    REAL( fp_kind ) :: dp_TL



    !#--------------------------------------------------------------------------#
    !#            -- Initialise 0'th LEVEL TL absorber amounts --               #
    !#                                                                          #
    !# This is done so that layer differences and averages can be calculated    #
    !# simply in the predictor and transmittance routines.                      #
    !#--------------------------------------------------------------------------#

    absorber_TL( 0, : ) = ZERO



    !#--------------------------------------------------------------------------#
    !#        -- Assemble the NADIR tangent-linear absorber profiles --         #
    !#--------------------------------------------------------------------------#

    ! -----------------------------------------------
    ! TOA layer. Profile inputs are dimensioned as K
    ! hence the top input layer is treated separately
    ! -----------------------------------------------

    absorber_TL( 1, 1 ) = RECIPROCAL_GRAVITY * (( water_vapor_TL( 1 ) * pressure( 1 )    ) + &
                                                ( water_vapor( 1 )    * pressure_TL( 1 ) ))
    absorber_TL( 1, 2 ) = pressure_TL( 1 )
    absorber_TL( 1, 3 ) = RECIPROCAL_GRAVITY * (( ozone_TL( 1 ) * pressure( 1 )    ) + &
                                                ( ozone( 1 )    * pressure_TL( 1 ) ))


    ! --------------------------------
    ! Loop over layers, TOA - 1 -> SFC
    ! --------------------------------

    k_layer_loop: DO k = 2, SIZE( pressure )

      ! -- Layer pressure differences
      dp    = pressure( k )    - pressure( k-1 )
      dp_TL = pressure_TL( k ) - pressure_TL( k-1 )

      ! -- Integrated TL absorber amounts
      absorber_TL( k, 1 ) = absorber_TL( k-1, 1 ) + &
                            ( RECIPROCAL_GRAVITY * (( water_vapor_TL( k ) * dp    ) + &
                                                    ( water_vapor( k )    * dp_TL )) )
      absorber_TL( k, 2 ) = pressure_TL( k )

      absorber_TL( k, 3 ) = absorber_TL( k-1, 3 ) + &
                            ( RECIPROCAL_GRAVITY * (( ozone_TL( k ) * dp    ) + &
                                                    ( ozone( k )    * dp_TL )) )

    END DO k_layer_loop

  END SUBROUTINE compute_absorber_amount_TL





!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_absorber_amount_AD
!
! PURPOSE:
!       PUBLIC subroutine to compute the adjoint of the integrated
!       profiles for all the absorbers. Currently the number of absorbers
!       are:
!         - Water vapor
!         - Dry/fixed gases (pressure == absorber amount)
!         - Ozone
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_absorber_amount_AD( &
!                                        ! -- Forward input
!                                        pressure,       &  ! Input, K
!                                        water_vapor,    &  ! Input, K
!                                        ozone,          &  ! Input, K
!
!                                        ! -- Adjoint input
!                                        absorber_AD,    &  ! In/Output, 0:K x J
!
!                                        ! -- Adjoint output
!                                        pressure_AD,    &  ! In/Output, K
!                                        water_vapor_AD, &  ! In/Output, K
!                                        ozone_AD        )  ! In/Output, K
!
! INPUT ARGUMENTS:
!       pressure:        Profile LEVEL pressure array.
!                        UNITS:      hPa
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN )
!
!       water_vapor:     Profile LAYER water vapor mixing ratio array.
!                        UNITS:      g/kg
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN )
!
!       ozone:           Profile LAYER ozone mixing ratio array.
!                        UNITS:      ppmv
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN )
!
!       absorber_AD:     Profile LEVEL adjoint average integrated 
!                        absorber amount array.
!                        ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                        UNITS:      Varies with absorber
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  0:K x J
!                        ATTRIBUTES: INTENT( OUT )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       pressure_AD:     Profile LAYER adjoint pressure array.
!                        UNITS:      hPa
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN OUT )
!
!       water_vapor_AD:  Profile LAYER adjoint water vapor array.
!                        UNITS:      g/kg
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN OUT )
!
!       ozone_AD:        Profile LAYER adjoint ozone water vapor array.
!                        UNITS:      ppmv
!                        TYPE:       REAL( fp_kind )
!                        DIMENSION:  K
!                        ATTRIBUTES: INTENT( IN OUT )
!
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! CALLS:
!       None
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       The input argument ABSORBER_AD is set to zero on output.
!
! RESTRICTIONS:
!       None.
!
! PROCEDURE:
!       The function calculates the adjoints of the integrated path length
!       through k atmospheric layers. Given the tangent-linear expression,
!
!                        __ k
!                    1  \
!         u_TL(k) = ---  >  [ q_TL(i).dp(i) + q(i).dp_TL(i) ]
!                    g  /__
!                          i=1
!
!       for the profile inputs (water vapor and ozone), the adjoints are
!       given by,
!
!                              1
!         dp_AD     = dp_AD + --- q(k).u_AD(k)
!                              g
!
!                                1
!         q_AD(k)   = q_AD(k) + --- dp.u_AD(k)
!                                g
!
!
!         u_AD(k-1) = u_AD(k-1) + u_AD(k)
!
!         u_AD(k)   = 0.0
!
!       looping over the layer index, k, from the surface (max.) to the
!       top of the atmosphere. For the dry absorber component, which uses
!       pressure as the absorber amount proxy, the tangent-linear expression
!       is simply,
!
!         u_TL(k) = p_TL(k)
!
!       The adjoint is thus,
!
!         p_AD(k) = p_AD(k) + u_AD(k)
!
!         u_AD(k) = 0.0
!
!       The actual execution of the above expressions in the code are done to
!       prevent recalculation of the the same quantities in individual loops.
!
!       Typically, the adjoint quantities that correspond to INPUTS in the
!       forward model (pressure_AD, water_vapor_AD, and ozone_AD) are set to
!       0.0 before entering this routine, and the adjoint quantities 
!       corresponding to OUTPUTS in the forward model (absorber_AD) are set to
!       1.0 before entering this routine. This will return the gradient vector
!       of the output variable with respect to the input variable (partial
!       derivative.) E.g. if the forward problem is defined as,
!
!
!         u = f(q,p)
!
!       then the TL form would be,
!
!                 df          df
!         u_TL = ---- q_TL + ---- p_TL   (NOTE: the df/dX is partial)
!                 dq          dp
!
!       If u_AD = 1.0 and q_AD,p_AD = 0.0 on input, then on output,
!
!                 df            df
!         q_AD = ---- , p_AD = ---- , and u_AD = 0.0
!                 dq            dp 
! 
!S-
!--------------------------------------------------------------------------------

  SUBROUTINE  compute_absorber_amount_AD( &
                                          ! -- Forward input
                                          pressure,       &  ! Input, K
                                          water_vapor,    &  ! Input, K
                                          ozone,          &  ! Input, K

                                          ! -- Adjoint input
                                          absorber_AD,    &  ! In/Output, 0:K x J

                                          ! -- Adjoint output
                                          pressure_AD,    &  ! In/Output, K
                                          water_vapor_AD, &  ! In/Output, K
                                          ozone_AD        )  ! In/Output, K



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward input
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )     :: pressure        ! Input, K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )     :: water_vapor     ! Input, K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )     :: ozone           ! Input, K

    ! -- Adjoint input
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN OUT ) :: absorber_AD     ! In/Output, 0:K x J

    ! -- Adjoint output
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN OUT ) :: pressure_AD     ! In/Output, K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN OUT ) :: water_vapor_AD  ! In/Output, K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN OUT ) :: ozone_AD        ! In/Output, K


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: k

    REAL( fp_kind ) :: dp
    REAL( fp_kind ) :: dp_AD



    !#--------------------------------------------------------------------------#
    !#           -- Assemble the NADIR adjoint absorber profiles --             #
    !#--------------------------------------------------------------------------#

    ! --------------------------------
    ! Loop over layers, SFC -> TOA - 1
    ! --------------------------------

    k_layer_loop: DO k = SIZE( pressure ), 2, -1

      ! -- Layer pressure differences
      dp = pressure( k ) - pressure( k-1 )

      ! -- Ozone absorber adjoint
      ozone_AD( k )       = ozone_AD( k ) + &
                            ( RECIPROCAL_GRAVITY * dp * absorber_AD( k, 3 ) )

      ! -- Pressure absorber adjoint
      pressure_AD( k )    = pressure_AD( k ) + absorber_AD( k, 2 )

      ! -- Water vapor absorber adjoint
      water_vapor_AD( k ) = water_vapor_AD( k ) + &
                            ( RECIPROCAL_GRAVITY * dp * absorber_AD( k, 1 ) )

      ! -- Layer pressure difference adjoint
      dp_AD = RECIPROCAL_GRAVITY * ( ( water_vapor( k ) * absorber_AD( k, 1 ) ) + &
                                     ( ozone( k )       * absorber_AD( k, 3 ) )   )
      pressure_AD( k )   = pressure_AD( k )   + dp_AD
      pressure_AD( k-1 ) = pressure_AD( k-1 ) - dp_AD
      ! NOTE: No dp_AD = ZERO here as it is reinitialised every layer

      ! -- Previous layer absorber amounts
      absorber_AD( k-1, 3 ) = absorber_AD( k-1, 3 ) + absorber_AD( k, 3 )
      absorber_AD( k,   3 ) = ZERO

      absorber_AD( k,   2 ) = ZERO

      absorber_AD( k-1, 1 ) = absorber_AD( k-1, 1 ) + absorber_AD( k, 1 )
      absorber_AD( k,   1 ) = ZERO

    END DO k_layer_loop


    ! -----------------------------------------------
    ! TOA layer. Profile inputs are dimensioned as K
    ! hence the top input layer is treated separately
    ! -----------------------------------------------

    ! -- Ozone
    ozone_AD( 1 )       = ozone_AD( 1 ) + &
                          ( RECIPROCAL_GRAVITY * pressure( 1 ) * absorber_AD( 1, 3 ) )

    ! -- Pressure
    pressure_AD( 1 )    = pressure_AD( 1 ) + &
                          ( RECIPROCAL_GRAVITY * ( ( ozone( 1 )       * absorber_AD( 1, 3 ) ) + &
                                                   ( water_vapor( 1 ) * absorber_AD( 1, 1 ) )   ) ) + &
                          absorber_AD( 1, 2 )
    ! -- Water vapor
    water_vapor_AD( 1 ) = water_vapor_AD( 1 ) + &
                          ( RECIPROCAL_GRAVITY * pressure( 1 ) * absorber_AD( 1, 1 ) )

    ! -- Zero absorber amount adjoints
    absorber_AD( 0:1, : ) = ZERO

  END SUBROUTINE compute_absorber_amount_AD






  SUBROUTINE find_absorber_layer_index( absorber,   &  ! Input,  0:K x J
                                        layer_index )  ! Output, K x J

    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN )  :: absorber
    INTEGER,         DIMENSION(  :, : ), INTENT( OUT ) :: layer_index

    INTEGER :: j,  n_absorbers
    INTEGER :: k,  n_layers

    REAL( fp_kind ) :: ave_absorber


    INTRINSIC ABS, &
              SIZE



    !#--------------------------------------------------------------------------#
    !#                  -- Determine the array dimension --                     #
    !#--------------------------------------------------------------------------#

    n_layers    = SIZE( absorber, DIM = 1 ) - 1
    n_absorbers = SIZE( absorber, DIM = 2 )



    !#--------------------------------------------------------------------------#
    !#                       -- Loop over absorbers --                          #
    !#--------------------------------------------------------------------------#

    j_absorber_loop: DO j = 1, n_absorbers


      ! ------------------------------
      ! Loop over the user layer space
      ! ------------------------------

      k_user_space_loop: DO k = 1, n_layers


        ! ---------------------------------------
        ! Calculate layer average absorber amount
        ! ---------------------------------------

        ave_absorber = POINT_5 * ( absorber( k, j ) + absorber( k-1, j ) )


        ! -------------------------------------------------------
        ! Calculate the absorber space bracket layer, ka, for:
        !
        !   abs_spc(ka-1) < absorber < abs_spc(ka)
        !
        ! This is done using the exponential parameter, alpha,
        ! used in determining the absorber space layers to 
        ! begin with:
        !
        !                   EXP( ka.alpha ) - 1
        !   A(ka) = A(1) . ---------------------
        !                    EXP( alpha ) - 1
        !
        ! so given an absorber amount, A(k), the required ka
        ! value is found using:
        !
        !          1        (  A(k)                        )
        !   k = ------- . LN( ------( EXP(alpha)-1 )  +  1 )
        !        alpha      (  A(1)                        )
        !
        ! The use of the CEILING intrinsic function provides the
        ! integer value equal to or greater than the argument.
        ! Then use of MAX( MIN( ka, MAX_N_ABSORBERS_LAYERS ), 1 )
        ! is to ensure that the calculated value never exceeds the 
        ! maximum number of absorber space layers - which could
        ! happen if an input absorber profile contained amounts
        ! which exceeded the maximum absorber space amount - and
        ! is never less than 1 - which could happen if the input
        ! absorber profile has adjacent layers with zero amounts.
        ! -------------------------------------------------------

        layer_index( k, j ) = MAX( MIN( INT( CEILING( ( ONE / alpha(j) ) * &
                                                 LOG( ( ave_absorber / absorber_space_levels(1,j) ) * &
                                                      ( EXP( alpha(j) ) - ONE ) + ONE ) ) ), &
                                   MAX_N_ABSORBER_LAYERS ), 1 )


      END DO k_user_space_loop

    END DO j_absorber_loop

  END SUBROUTINE find_absorber_layer_index

END MODULE absorber_profile


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
! $Name$
!
! $State$
!
! $Log$
! Revision 1.11  2001/09/25 15:51:29  paulv
! - Changed the calculation of the bracketing absorber space layer in
!   sbroutine FIND_ABSORBER_LAYER_INDEX from
!     MIN( ka, MAX_N_ABSORBERS_LAYERS )
!   to
!     MAX( MIN( ka, MAX_N_ABSORBERS_LAYERS ), 1 )
!   so as to avoid the result ever being zero - which could happen before if
!   adjacent layers of the input absorber profile were zero.
!
! Revision 1.10  2001/08/31 20:41:18  paulv
! - Altered method of searching for bracketing absorber space layers in
!   FIND_ABSORBER_LAYER_INDEX. Previosuly a trickle down search was performed.
!   Now the actual corresponding layer is calculated using the exponential
!   factor used in generating the absorber space.
!
! Revision 1.9  2001/08/16 16:30:38  paulv
! - Updated documentation.
!
! Revision 1.8  2001/08/01 16:36:34  paulv
! - Removed use of module ABSORBER_SPACE and replaced it with
!     USE transmittance_coefficients, ONLY : absorber_space_levels
!   to reflect changes in code. The absorber space levels are no longer
!   calculated during model initialisation, but are precalculated and stored
!   in the transmittance coefficient data file.
!
! Revision 1.7  2001/06/05 21:18:10  paulv
! - Changed adjoint routine slightly to make adjoint calcs a bit clearer
!   when looking at the tangent-linear code.
! - Corrected bug in TOA layer pressure_AD calculation.
!
! Revision 1.6  2001/05/29 17:32:51  paulv
! - Some cosmetic changes
! - Removed subtraction of the TOA_PRESSURE parameter from the DRY absorber
!   calculation. This was causing the upper level channels to produce
!   spurious numbers in the forward calculation.
! - Added the  FIND_ABSORBER_LAYER_INDEX routine. Removed it from the FORWARD_MODEL
!   module. It seemed more appropriate in this one.
! - Using pressure array data directly in first layer calcs rather than
!   DP variable.
!
! Revision 1.5  2001/03/26 18:45:59  paulv
! - Now use TYPE_KINDS module parameter FP_KIND to set the floating point
!   data type.
! - Module parameter RECIPROCAL_GRAVITY moved to PARAMETERS module.
! - ONLY clause used in USE PARAMETERS statement. Only parameters available
!   in ABSORBER_PROFILE module are ZERO, TOA_PRESSURE, and RECIPROCAL_GRAVITY.
! - Output ABSORBER argument is now dimensioned as 0:K. This eliminates the
!   need for using an ABSORBER_KM1 variable in other routines that use the
!   ABSORBER array variable where the layer loop always goes from 1 -> n_layers.
! - Removed output arguments of AVE_ABSORBER and DELTA_ABSORBER due to the
!   ABSORBER dimension change to 0:K. Calculating the average and layer
!   difference absorber amounts in other routines can now be done simply
!   by averaging or differencing ABOSRBER(K) and ABSORBER(K-1) even for
!   layer #1.
! - Integration of absorber amount for the TOA layer is done OUTSIDE of the
!   layer loop. This avoids the need for a PRESSURE_KM1 variable since
!   pressure is dimensioned 1:K.
! - Layer loop, thus, goes from 2 -> n_layers.
! - Changed order of argument list in COMPUTE_PREDICTORS_TL and its
!   associated routines. All forward arguments are listed followed by
!   the tangent-linear arguments rather than interspersing them as before.
! - Added adjoint routine COMPUTE_ABSORBER_AMOUNT_AD.
!
! Revision 1.4  2001/03/26 18:30:54  paulv
! - Removed integrate_absorber_profile and integrate_absorber_profile_tl
!   functions. Integration is now done in-line in the main routines.
!
! Revision 1.3  2000/11/09 20:29:40  paulv
! - Added tangent linear form of the absorber profile routines.
!
! Revision 1.2  2000/08/31 19:36:31  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 1.1  2000/08/24 13:11:27  paulv
! Initial checkin.
!
!
!
!
