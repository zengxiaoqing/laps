!------------------------------------------------------------------------------
!M+
! NAME:
!       transmittance
!
! PURPOSE:
!       RT model transmittance module
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE transmittance
!
! OUTPUTS:
!       None.
!
! MODULES:
!       type_kinds:                 Module containing data type kind definitions.
!
!       parameters:                 Module containing parameter definitions for the
!                                   RT model.
!
!       absorber_space:             Module containing the RT model absorber space
!                                   definitions.
!
!       transmitance_coefficients:  Module containing the RT model absorber space
!                                   definitions and transmittance coefficients.
!
!       predictors:                 Module containing the predictor calculation
!                                   routines and data.
!
! CONTAINS:
!       compute_transmittance:      PUBLIC subroutine to calculate the layer
!                                   transmittances of an input atmospheric profile.
!
!       compute_transmittance_TL:   PUBLIC subroutine to calculate the tangent-linear
!                                   layer transmittances of an input atmospheric profile.
!
!       compute_transmittance_AD:   PUBLIC subroutine to calculate the layer transmittance
!                                   adjoints of an input atmospheric profile.
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
!       Written by:     Paul van Delst, CIMSS/SSEC 11-Jul-2000
!                       paul.vandelst@ssec.wisc.edu
!
!       Adapted from code written by: Thomas J.Kleespies
!                                     NOAA/NESDIS/ORA
!                                     thomas.j.kleespies@noaa.gov
!
!                                     and
!
!                                     John Derber
!                                     NOAA/NCEP/EMC
!                                     john.derber@noaa.gov
!
!
!  Copyright (C) 2000 Thomas Kleespies, John Derber, Paul van Delst
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

MODULE transmittance

  ! ------------
  ! Module usage
  ! ------------

  USE type_kinds, ONLY : fp_kind
  USE parameters
  USE transmittance_coefficients


  ! -----------------------
  ! Disable implicit typing
  ! -----------------------

  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------

  PRIVATE
  PUBLIC :: compute_transmittance
  PUBLIC :: compute_transmittance_TL
  PUBLIC :: compute_transmittance_AD


CONTAINS


!------------------------------------------------------------------------------
!S+
! NAME:
!       compute_transmittance
!
! PURPOSE:
!       PUBLIC subroutine to calculate the layer transmittances given an
!       input atmospheric profile.
!
! CALLING SEQUENCE:
!       CALL compute_transmittance( absorber,      &   ! Input, 0:K x J
!                                   predictor,     &   ! Input, I x K
!                                   layer_index,   &   ! Input, K x J
!                                   channel_index, &   ! Input, scalar
!                                   direction,     &   ! Input, scalar
!
!                                   tau            )   ! Output, K
!
! INPUT ARGUMENTS:
!       absorber:         Profile LEVEL integrated absorber amount array.
!                         UNITS:      Varies with absorber.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  0:K x J
!                         ATTRIBUTES: INTENT( IN )
!
!       predictor:        Profile LAYER predictors array.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  I x K
!                         ATTRIBUTES: INTENT( IN )
!
!       layer_index:      Index array array associating the input absorber
!                         layer amount to the bracketing absorber space levels.
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  K x J
!                         ATTRIBUTES: INTENT( IN )
!
!       channel_index:    Channel index id. This is a unique index associated
!                         with a (supported) sensor channel.
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN )
!
!       direction:        Direction identifier.
!                         If = 0, calculate layer->surface transmittances (i.e. down)
!                            = 1, calculate layer->space   transmittances (i.e. up)
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!        None.
!
! OUTPUT ARGUMENTS:
!        tau:             Layer to boundary transmittances for the input atmosphere
!                         and channel.
!                         UNITS:      None
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( OUT )
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
!       McMillin, L.M., L.J. Crone, M.D. Goldberg, and T.J. Kleespies,
!         "Atmospheric transmittance of an absorbing gas. 4. OPTRAN: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", Applied Optics, 1995, v34, pp6269-6274.
!
!       The layer absorption coefficient is calculated using,
!
!                                   __ I
!                                  \   _
!         abs_coeff(k) = b(0,k,l) + >  b(i,k,l).X(k,l)
!                                  /__
!                                     i=1
!             _
!       where b == interpolated coefficients
!             X == predictors
!             i == predictor index
!             k == layer index
!             l == channel index
!
!       The input coefficients are linearly interpolated in absorber space
!       from the precalculated absorber levels to that defined by the
!       input profile.
!
!       The absorber layer optical depth is then calculated using,
!
!         optical_depth(k) = abs_coeff(k) * dA(k)
!
!       where dA == layer absorber difference.
!
!       The layer transmittance is then calculated using,
!
!         tau(k) = EXP( -optical_depth(k) )
!S-
!------------------------------------------------------------------------------

  SUBROUTINE compute_transmittance( absorber,      &   ! Input, 0:K x J
                                    predictor,     &   ! Input, I x K
                                    layer_index,   &   ! Input, K x J
                                    channel_index, &   ! Input, scalar
                                    direction,     &   ! Input, scalar, 0==down, 1==up

                                    tau            )   ! Output, K



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN )  :: absorber         ! Input, 0:K x J
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )  :: predictor        ! Input, I x K
    INTEGER,         DIMENSION( :, : ),  INTENT( IN )  :: layer_index      ! Input, K x J
    INTEGER,                             INTENT( IN )  :: channel_index    ! Input, scalar
    INTEGER,                             INTENT( IN )  :: direction        ! Input, scalar, 0==down, 1==up

    REAL( fp_kind ), DIMENSION( : ),     INTENT( OUT ) :: tau              ! Output, K


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_TRANSMITTANCE'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: l
    INTEGER :: k, k1, k2, dk, n_layers
    INTEGER :: j, n_absorbers
    INTEGER :: i, ip, n_predictors

    REAL( fp_kind ) :: ave_absorber, d_absorber
    REAL( fp_kind ) :: gradient, b, b1, b2
    REAL( fp_kind ) :: absorption_coefficient
    REAL( fp_kind ) :: total_od, od_tolerance

    REAL( fp_kind ), DIMENSION( SIZE( tau ) ) :: optical_depth


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC EXP, &
              PRESENT, &
              SIZE



    !#--------------------------------------------------------------------------#
    !#                   -- Determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    ! -- Number of atmospheric layers and absorbers. The "-1"
    ! -- for the layer assign is because absorber is based on
    ! -- LEVELS.
    n_layers    = SIZE( absorber, DIM = 1 ) - 1
    n_absorbers = SIZE( absorber, DIM = 2 )

    ! -- Number of predictors *TO USE*. Currently this is
    ! -- pretty much hardwired to 5 but it could change
    ! -- in the future. This value is explicitly
    ! -- checked as the TAU_COEFFICIENTS array contained
    ! -- in the TRANSMITTANCE_COEFFICIENTS module is
    ! -- allocated dynamically when initialising the 
    ! -- radiative transfer model.
 
    ! -- The "-1" is used because the first element of the
    ! -- predictor index array is used to determine if there
    ! -- is *any* absorption for a particular absorber species.
    n_predictors = SIZE( tau_coefficients, DIM = 1 ) - 1



    !#--------------------------------------------------------------------------#
    !#                -- Calculate the layer optical depths --                  #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------
    ! Assign the channel index to a short name
    ! ----------------------------------------

    l = channel_index


    ! ---------------------------
    ! Initilise the optical depth
    ! ---------------------------

    optical_depth( : ) = ZERO


    ! -----------------------------------------------------
    ! Loop over each absorber for optical depth calculation
    ! -----------------------------------------------------

    j_absorber_loop: DO j = 1, n_absorbers


      ! -----------------------------------------
      ! Check if there is any absorption for this
      ! absorber/channel combination.
      !
      ! This check is the reason why all channels
      ! cannot be processed at once and why the
      ! layer loop is within the absorber loop.
      ! -----------------------------------------

      IF ( predictor_index( 0, l, j ) == 0 ) CYCLE j_absorber_loop



      !#------------------------------------------------------------------------#
      !#                    -- Begin loop over layers --                        #
      !#------------------------------------------------------------------------#

      k_layer_od_loop: DO k = 1, n_layers


        ! -----------------------------------
        ! Calculate the current layer average
        ! absorber amount and difference
        ! -----------------------------------

        ave_absorber = POINT_5 * ( absorber( k, j ) + absorber( k-1, j ) )
        d_absorber   = absorber( k, j ) - absorber( k-1, j )


        ! -----------------------------------------------------------
        ! To linearly interpolate the tau_coeffs to the actual user
        ! space absorber amount, need the gradient across the layer
        ! -----------------------------------------------------------

        k2 = layer_index( k, j )
        k1 = k2 - 1

        gradient = ( absorber_space_levels( k2, j ) - ave_absorber                   ) / &
        !          -------------------------------------------------------------------
                   ( absorber_space_levels( k2, j ) - absorber_space_levels( k1, j ) )


        ! -------------------------------
        ! Calculate absorption coeficient
        ! -------------------------------

        ! -- Offset term
        b1 = tau_coefficients( 0, k1, l, j )
        b2 = tau_coefficients( 0, k2, l, j )
        absorption_coefficient = b2 + ( gradient * ( b1 - b2 ) )


        i_predictor_loop: DO i = 1, n_predictors

          ! -- Current predictor interpolated coefficient
          b1 = tau_coefficients( i, k1, l, j )
          b2 = tau_coefficients( i, k2, l, j )
          b  = b2 + ( gradient * ( b1 - b2 ) )

          !  - Sum current predictor's contribution
          ip = predictor_index( i, l, j )
          absorption_coefficient = absorption_coefficient + ( b * predictor( ip, k ) ) 

        END DO i_predictor_loop


        ! -- **** IT WOULD BE NICE TO NOT NEED THIS AT ALL! ****
        absorption_coefficient = MAX( absorption_coefficient, ZERO )


        ! -----------------------
        ! Calculate optical_depth
        ! -----------------------

        optical_depth( k ) = optical_depth( k ) + &
                             ( absorption_coefficient * d_absorber )


      END DO k_layer_od_loop

    END DO j_absorber_loop



    !#--------------------------------------------------------------------------#
    !#           -- Calculate the layer->boundary transmittances --             #
    !#                                                                          #
    !# This step involves another loop over layers. One *could* reverse the     #
    !# order of the j absorber loop and k layer od loop and calculate the tau   #
    !# values outside the absorber loop. However, this would involve an IF test #
    !# for every layer even if PREDICTOR_INDEX( 0, l ) == 0. I figured it is    #
    !# better to skip the layer loop altogether if there is no absorption for   #
    !# a particular absorber, despite the need for an extra layer loop here to  #
    !# calculate the transmittance.                                             #
    !#--------------------------------------------------------------------------#

    ! ----------------------
    ! Determine loop indices
    ! ----------------------

    IF ( direction == UP ) THEN
      ! -- Transmittance going up, e.g. UPWELLING. Layer 1 => TOA
      k1 = 1
      k2 = n_layers
      dk = 1
    ELSE
      ! -- Transmittance going down, e.g. DOWNWELLING FLUX, SOLAR. Layer N_LAYERS => SFC
      k1 =  n_layers
      k2 =  1
      dk = -1
    END IF


    ! ------------------------------------------
    ! Initialise the total optical depth and set
    ! its tolerance value. The tolerance value
    ! should prevent floating point underflows.
    ! ------------------------------------------

    total_od     = ZERO
    od_tolerance = ABS( LOG( TOLERANCE ) )


    ! ----------------------------------------------
    ! Loop over layers for transmittance calculation
    ! ----------------------------------------------

    k_layer_tau_loop: DO k = k1, k2, dk

      ! -- Update total optical depth
      total_od = total_od + optical_depth( k )

      ! -- If optical depth is < than tolerance, calculate transmittance
      IF ( total_od < od_tolerance ) THEN

        tau( k ) = EXP( -total_od )

      ELSE

        ! -- The following *could* be slow if dk is negative
        ! -- If so, then tau( k ) = ZERO would work.
        tau( k:k2:dk ) = ZERO
        EXIT k_layer_tau_loop

      END IF

    END DO k_layer_tau_loop

  END SUBROUTINE compute_transmittance




!------------------------------------------------------------------------------
!S+
! NAME:
!       compute_transmittance_TL
!
! PURPOSE:
!       PUBLIC subroutine to calculate the tangent-linear layer transmittances
!       of an input atmospheric profile.
!
! CALLING SEQUENCE:
!        CALL compute_transmittance_TL( &
!                                       ! -- Forward input
!                                       absorber,      &   ! Input, 0:K x J
!                                       predictor,     &   ! Input, I x K
!                                       tau,           &   ! Input, K
!
!                                       ! -- Tangent-liner input
!                                       absorber_TL,   &   ! Input, 0:K x J
!                                       predictor_TL,  &   ! Input, I x K
!
!                                       ! -- Other input
!                                       layer_index,   &   ! Input, K x J
!                                       channel_index, &   ! Input, scalar
!                                       direction,     &   ! Input, scalar, 0==down, 1==up
!
!                                       ! -- Tangent-liner output
!                                       tau_TL         )   ! Output, K
!
! INPUT ARGUMENTS:
!       absorber:         Profile LEVEL integrated absorber amount array.
!                         UNITS:      Varies with absorber.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  0:K x J
!                         ATTRIBUTES: INTENT( IN )
!
!       predictor:        Profile LAYER predictors array.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  I x K
!                         ATTRIBUTES: INTENT( IN )
!
!       absorber_TL:      Profile LEVEL tangent-linear integrated absorber
!                         amount array.
!                         UNITS:      Varies with absorber.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  0:K x J
!                         ATTRIBUTES: INTENT( IN )
!
!       predictor_TL:     Profile LAYER tangent-linear predictors array.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  I x K
!                         ATTRIBUTES: INTENT( IN )
!
!       layer_index:      Index array array associating the input absorber
!                         layer amount to the bracketing absorber space levels.
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  K x J
!                         ATTRIBUTES: INTENT( IN )
!
!       channel_index:    Channel index id. This is a unique index associated
!                         with a (supported) sensor channel.
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN )
!
!       direction:        Direction identifier.
!                         If = 0, calculate layer->surface transmittances (i.e. down)
!                            = 1, calculate layer->space   transmittances (i.e. up)
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!        None.
!
! OUTPUT ARGUMENTS:
!        tau_TL:          Layer to boundary tangent-linear transmittances for the
!                         input atmosphere and channel.
!                         UNITS:      None
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( OUT )
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
!       McMillin, L.M., L.J. Crone, M.D. Goldberg, and T.J. Kleespies,
!         "Atmospheric transmittance of an absorbing gas. 4. OPTRAN: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", Applied Optics, 1995, v34, pp6269-6274.
!
!       The layer absorption coefficient is calculated using,
!
!                                   __ I
!                        _         \   _
!         abs_coeff(k) = b(0,k,l) + >  b(i,k,l).X(k,l)
!                                  /__
!                                     i=1
!
!       and the tangent-linear absorption coefficient from,
!
!                                         __ I
!                           _            \   _                    _
!         abs_coeff_TL(k) = b_TL(0,k,l) + >  b_TL(i,k,l).X(k,l) + b(i,k,l).X_TL(k,l)
!                                        /__
!                                           i=1
!             _  _
!       where b, b_TL == interpolated coefficients
!             X, X_TL == predictors
!             i       == predictor index
!             k       == layer index
!             l       == channel index
!
!       The input coefficients are linearly interpolated in absorber space
!       from the precalculated absorber levels to that defined by the
!       input profile.
!
!       The absorber layer optical depth is then calculated using,
!
!         optical_depth(k) = abs_coeff(k).dA(k)
!
!       and the tangent-linear term is determined from,
!
!         optical_depth_TL(k) = abs_coeff_TL(k).dA(k) + abs_coeff(k).dA_TL(k)
!
!       where dA, dA_TL == layer absorber difference.
!
!       The tangent-linear transmittance is then calculated using,
!
!                               __ K
!                              \
!         tau_TL(k) = -tau(k) . >  optical_depth_TL(i)
!                              /__
!                                 i=k
!
!S-
!------------------------------------------------------------------------------

  SUBROUTINE compute_transmittance_TL( &
                                       ! -- Forward input
                                       absorber,      &   ! Input, 0:K x J
                                       predictor,     &   ! Input, I x K
                                       tau,           &   ! Input, K

                                       ! -- Tangent-liner input
                                       absorber_TL,   &   ! Input, 0:K x J
                                       predictor_TL,  &   ! Input, I x K

                                       ! -- Other input
                                       layer_index,   &   ! Input, K x J
                                       channel_index, &   ! Input, scalar
                                       direction,     &   ! Input, scalar, 0==down, 1==up

                                       ! -- Tangent-liner output
                                       tau_TL         )   ! Output, K



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward input
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN )  :: absorber         ! 0:K x J
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )  :: predictor        ! I x K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: tau              ! K

    ! -- Tangent_linear input
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN )  :: absorber_TL      ! 0:K x J
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )  :: predictor_TL     ! I x K

    ! -- Other input
    INTEGER,         DIMENSION( :, : ),  INTENT( IN )  :: layer_index      ! K x J
    INTEGER,                             INTENT( IN )  :: channel_index    ! scalar
    INTEGER,                             INTENT( IN )  :: direction        ! scalar, 0==down, 1==up

    ! -- Tangent_linear output
    REAL( fp_kind ), DIMENSION( : ),     INTENT( OUT ) :: tau_TL           ! K


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_TRANSMITTANCE_TL'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: l
    INTEGER :: k, k1, k2, dk, n_layers
    INTEGER :: j, n_absorbers
    INTEGER :: i, ip, n_predictors

    REAL( fp_kind ) :: ave_absorber,    d_absorber
    REAL( fp_kind ) :: ave_absorber_TL, d_absorber_TL

    REAL( fp_kind ) :: b1, b2
    REAL( fp_kind ) :: gradient,    b
    REAL( fp_kind ) :: gradient_TL, b_TL

    REAL( fp_kind ) :: absorption_coefficient
    REAL( fp_kind ) :: absorption_coefficient_TL

    REAL( fp_kind ) :: total_od_TL

    REAL( fp_kind ), DIMENSION( SIZE( tau ) ) :: optical_depth_TL


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC EXP, &
              PRESENT, &
              SIZE



    !#--------------------------------------------------------------------------#
    !#                   -- Determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    ! -- Number of atmospheric layers and absorbers. The "-1"
    ! -- for the layer assign is because absorber is based on
    ! -- LEVELS.
    n_layers    = SIZE( absorber, DIM = 1 ) - 1
    n_absorbers = SIZE( absorber, DIM = 2 )

    ! -- Number of predictors *TO USE*. Currently this is
    ! -- pretty much hardwired to 5 but it could change
    ! -- in the future. This value is explicitly
    ! -- checked as the TAU_COEFFICIENTS array contained
    ! -- in the TRANSMITTANCE_COEFFICIENTS module is
    ! -- allocated dynamically when initialising the 
    ! -- radiative transfer model.
 
    ! -- The "-1" is used because the first element of the
    ! -- predictor index array is used to determine if there
    ! -- is *any* absorption for a particular absorber species.
    n_predictors = SIZE( tau_coefficients, DIM = 1 ) - 1



    !#--------------------------------------------------------------------------#
    !#                -- Calculate the layer optical depths --                  #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------
    ! Assign the channel index to a short name
    ! ----------------------------------------

    l = channel_index


    ! ------------------------------------------
    ! Initilise the tangent-linear optical depth
    ! ------------------------------------------

    optical_depth_TL( : ) = ZERO


    ! -----------------------------------------------------
    ! Loop over each absorber for optical depth calculation
    ! -----------------------------------------------------

    j_absorber_loop: DO j = 1, n_absorbers


      ! -----------------------------------------
      ! Check if there is any absorption for this
      ! absorber/channel combination.
      !
      ! This check is the reason why all channels
      ! cannot be processed at once and why the
      ! layer loop is within the absorber loop.
      ! -----------------------------------------

      IF ( predictor_index( 0, l, j ) == 0 ) CYCLE j_absorber_loop



      !#------------------------------------------------------------------------#
      !#                    -- Begin loop over layers --                        #
      !#------------------------------------------------------------------------#

      k_layer_od_loop: DO k = 1, n_layers


        ! -----------------------------------
        ! Calculate the current layer average
        ! absorber amounts and differences
        ! -----------------------------------

        ave_absorber    = POINT_5 * ( absorber(    k, j ) + absorber(    k-1, j ) )
        ave_absorber_TL = POINT_5 * ( absorber_TL( k, j ) + absorber_TL( k-1, j ) )

        d_absorber    = absorber(    k, j ) - absorber(    k-1, j )
        d_absorber_TL = absorber_TL( k, j ) - absorber_TL( k-1, j )


        ! -----------------------------------------------------------
        ! To linearly interpolate the tau_coeffs to the actual user
        ! space absorber amount, need the gradient across the layer
        ! -----------------------------------------------------------

        k2 = layer_index( k, j )
        k1 = k2 - 1

        gradient = ( absorber_space_levels( k2, j ) - ave_absorber                   ) / &
        !          -------------------------------------------------------------------
                   ( absorber_space_levels( k2, j ) - absorber_space_levels( k1, j ) )

        gradient_TL =                        ( -ave_absorber_TL )                         / &
        !             -------------------------------------------------------------------
                      ( absorber_space_levels( k2, j ) - absorber_space_levels( k1, j ) )


        ! --------------------------------
        ! Calculate absorption coeficients
        ! --------------------------------

        ! -- Offset term
        b1 = tau_coefficients( 0, k1, l, j )
        b2 = tau_coefficients( 0, k2, l, j )

        absorption_coefficient    = b2 + ( gradient    * ( b1 - b2 ) )
        absorption_coefficient_TL =        gradient_TL * ( b1 - b2 )


        i_predictor_loop: DO i = 1, n_predictors

          ! -- Current predictor interpolated coefficient
          b1 = tau_coefficients( i, k1, l, j )
          b2 = tau_coefficients( i, k2, l, j )

          b     = b2 + ( gradient    * ( b1 - b2 ) )
          b_TL  =        gradient_TL * ( b1 - b2 )

          !  - Sum current predictor's contribution
          ip = predictor_index( i, l, j )

          absorption_coefficient    = absorption_coefficient + ( b * predictor( ip, k ) ) 
          absorption_coefficient_TL = absorption_coefficient_TL        + &
                                      ( b_TL * predictor(    ip, k ) ) + &
                                      ( b    * predictor_TL( ip, k ) )

        END DO i_predictor_loop


        ! -- **** IT WOULD BE NICE TO NOT NEED THIS AT ALL! ****
        absorption_coefficient = MAX( absorption_coefficient, ZERO )


        ! --------------------------------------
        ! Calculate tangent-linear optical depth
        ! --------------------------------------

        optical_depth_TL( k ) = optical_depth_TL( k ) + &
                                ( absorption_coefficient_TL * d_absorber    ) + &
                                ( absorption_coefficient    * d_absorber_TL )

      END DO k_layer_od_loop

    END DO j_absorber_loop



    !#--------------------------------------------------------------------------#
    !#    -- Calculate the layer->boundary tangent-linear transmittances --     #
    !#                                                                          #
    !# This step involves another loop over layers. One *could* reverse the     #
    !# order of the j absorber loop and k layer od loop and calculate the tau   #
    !# values outside the absorber loop. However, this would involve an IF test #
    !# for every layer even if PREDICTOR_INDEX( 0, l ) == 0. I figured it is    #
    !# better to skip the layer loop altogether if there is no absorption for   #
    !# a particular absorber, despite the need for an extra layer loop here to  #
    !# calculate the transmittance.                                             #
    !#--------------------------------------------------------------------------#

    ! ----------------------
    ! Determine loop indices
    ! ----------------------

    IF ( direction == UP ) THEN
      ! -- Transmittance going up, e.g. UPWELLING. Layer 1 => TOA
      k1 = 1
      k2 = n_layers
      dk = 1
    ELSE
      ! -- Transmittance going down, e.g. DOWNWELLING FLUX, SOLAR. Layer N_LAYERS => SFC
      k1 =  n_layers
      k2 =  1
      dk = -1
    END IF


    ! ----------------------------------------------
    ! Loop over layers for transmittance calculation
    ! ----------------------------------------------

    total_od_TL = ZERO

    k_layer_tau_loop: DO k = k1, k2, dk

      total_od_TL = total_od_TL + optical_depth_TL( k )
      tau_TL( k ) = -tau( k ) * total_od_TL

    END DO k_layer_tau_loop

  END SUBROUTINE compute_transmittance_TL





!------------------------------------------------------------------------------
!S+
! NAME:
!       compute_transmittance_AD
!
! PURPOSE:
!       PUBLIC subroutine to calculate the adjoint of the layer transmittances
!       of an input atmospheric profile.
!
! CALLING SEQUENCE:
!       CALL compute_transmittance_AD( &
!                                      ! -- Forward input
!                                      absorber,      &   ! Input, 0:K x J
!                                      predictor,     &   ! Input, I x K
!                                      tau,           &   ! Input, K
!
!                                      ! -- Adjoint input
!                                      tau_AD,        &   ! In/Output, K
!
!                                      ! -- Other input
!                                      layer_index,   &   ! Input, K x J
!                                      channel_index, &   ! Input, scalar
!                                      direction,     &   ! Input, scalar, 0==down, 1==up
!
!                                      ! -- Adjoint output
!                                      absorber_AD,   &   ! In/Output, 0:K x J
!                                      predictor_AD   )   ! In/Output, I x K
!
! INPUT ARGUMENTS:
!       absorber:         Profile LEVEL integrated absorber amount array.
!                         UNITS:      Varies with absorber.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  0:K x J
!                         ATTRIBUTES: INTENT( IN )
!
!       predictor:        Profile LAYER predictors array.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  I x K
!                         ATTRIBUTES: INTENT( IN )
!
!       tau:              Layer to boundary transmittances for the input
!                         atmosphere and channel.
!                         UNITS:      None
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( OUT )
!
!       tau_AD:           Adjoint of the layer to boundary transmittances
!                         for the input atmosphere and channel.
!                         **THIS ARGUMENT IS SET TO ZERO ON OUTPUT.**
!                         UNITS:      None
!                         TYPE:       Real
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN OUT )
!
!       layer_index:      Index array array associating the input absorber
!                         layer amount to the bracketing absorber space levels.
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  K x J
!                         ATTRIBUTES: INTENT( IN )
!
!       channel_index:    Channel index id. This is a unique index associated
!                         with a (supported) sensor channel.
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN )
!
!       direction:        Direction identifier.
!                         If = 0, calculate layer->surface transmittances (i.e. down)
!                            = 1, calculate layer->space   transmittances (i.e. up)
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       absorber_AD:      Adjoint of the profile LEVEL integrated absorber
!                         amount array.
!                         UNITS:      Varies with absorber.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  0:K x J
!                         ATTRIBUTES: INTENT( IN OUT )
!
!       predictor_AD:     Adjoint of the profile LAYER predictors.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  I x K
!                         ATTRIBUTES: INTENT( IN OUT )
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
!       The input argument TAU_AD is set to zero upon output.
!
! RESTRICTIONS:
!       None.
!
! PROCEDURE:
!       McMillin, L.M., L.J. Crone, M.D. Goldberg, and T.J. Kleespies,
!         "Atmospheric transmittance of an absorbing gas. 4. OPTRAN: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", Applied Optics, 1995, v34, pp6269-6274.
!
!S-
!------------------------------------------------------------------------------

  SUBROUTINE compute_transmittance_AD( &
                                       ! -- Forward input
                                       absorber,      &   ! Input, 0:K x J
                                       predictor,     &   ! Input, I x K
                                       tau,           &   ! Input, K

                                       ! -- Adjoint input
                                       tau_AD,        &   ! In/Output, K

                                       ! -- Other input
                                       layer_index,   &   ! Input, K x J
                                       channel_index, &   ! Input, scalar
                                       direction,     &   ! Input, scalar, 0==down, 1==up

                                       ! -- Adjoint output
                                       absorber_AD,   &   ! In/Output, 0:K x J
                                       predictor_AD   )   ! In/Output, I x K



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward input
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN )     :: absorber       ! 0:K x J
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )     :: predictor      ! I x K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )     :: tau            ! K

    ! -- Adjoint input
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN OUT ) :: tau_AD         ! K

    ! -- Other input
    INTEGER,         DIMENSION( :, : ),  INTENT( IN )     :: layer_index    ! K x J
    INTEGER,                             INTENT( IN )     :: channel_index  ! scalar
    INTEGER,                             INTENT( IN )     :: direction      ! scalar, 0==down, 1==up

    ! -- Adjoint output
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN OUT ) :: absorber_AD    ! 0:K x J
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN OUT ) :: predictor_AD   ! I x K


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_TRANSMITTANCE_AD'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: l
    INTEGER :: k, k1, k2, dk, n_layers
    INTEGER :: j, n_absorbers
    INTEGER :: i, ip, n_predictors

    REAL( fp_kind ) :: ave_absorber,    d_absorber
    REAL( fp_kind ) :: ave_absorber_AD, d_absorber_AD

    REAL( fp_kind ) :: d_absorber_space

    REAL( fp_kind ) :: b1o, b2o
    REAL( fp_kind ) :: b1,  b2
    REAL( fp_kind ) :: gradient,    b
    REAL( fp_kind ) :: gradient_AD, b_AD

    REAL( fp_kind ) :: absorption_coefficient
    REAL( fp_kind ) :: absorption_coefficient_AD

    REAL( fp_kind ) :: total_od_AD

    REAL( fp_kind ), DIMENSION( SIZE( tau ) ) :: optical_depth_AD


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC SIZE



    !#--------------------------------------------------------------------------#
    !#                   -- Determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    ! -- Number of atmospheric layers and absorbers. The "-1"
    ! -- for the layer assign is because absorber is based on
    ! -- LEVELS.
    n_layers    = SIZE( absorber, DIM = 1 ) - 1
    n_absorbers = SIZE( absorber, DIM = 2 )

    ! -- Number of predictors *TO USE*. Currently this is
    ! -- pretty much hardwired to 5 but it could change
    ! -- in the future. This value is explicitly
    ! -- checked as the TAU_COEFFICIENTS array contained
    ! -- in the TRANSMITTANCE_COEFFICIENTS module is
    ! -- allocated dynamically when initialising the 
    ! -- radiative transfer model.
 
    ! -- The "-1" is used because the first element of the
    ! -- predictor index array is used to determine if there
    ! -- is *any* absorption for a particular absorber species.
    n_predictors = SIZE( tau_coefficients, DIM = 1 ) - 1



    !#--------------------------------------------------------------------------#
    !#                       -- SOME INITIALISATIONS --                         #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------
    ! Assign the channel index to a short name
    ! ----------------------------------------

    l = channel_index


    ! -------------------------------------
    ! Initilise the local adjoint variables
    ! -------------------------------------

    optical_depth_AD( : ) = ZERO

    total_od_AD = ZERO
    gradient_AD = ZERO



    !#--------------------------------------------------------------------------#
    !#       -- CALCULATE THE LAYER->BOUNDARY TRANSMITTANCE ADJOINT --          #
    !#                                                                          #
    !# This step involves another loop over layers. One *could* reverse the     #
    !# order of the j absorber loop and k layer od loop and calculate the tau   #
    !# values outside the absorber loop. However, this would involve an IF test #
    !# for every layer even if PREDICTOR_INDEX( 0, l ) == 0. I figured it is    #
    !# better to skip the layer loop altogether if there is no absorption for   #
    !# a particular absorber, despite the need for an extra layer loop here to  #
    !# calculate the transmittance.                                             #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------
    ! Determine loop indices. Same direction
    ! as FORWARD and TANGENT-LINEAR models
    ! --------------------------------------

    IF ( direction == UP ) THEN
      ! -- Transmittance going up, e.g. UPWELLING. Layer 1 => TOA
      k1 = 1
      k2 = n_layers
      dk = 1
    ELSE
      ! -- Transmittance going down, e.g. DOWNWELLING FLUX, SOLAR. Layer N_LAYERS => SFC
      k1 =  n_layers
      k2 =  1
      dk = -1
    END IF


    ! -----------------------------------------------
    ! Loop over layers for transmittance calculation.
    ! Note the loop index reversal.
    ! -----------------------------------------------

    k_layer_tau_loop: DO k = k2, k1, -dk

      total_od_AD = total_od_AD - ( tau( k ) * tau_AD( k ) )
      tau_AD( k ) = ZERO

      optical_depth_AD( k ) = optical_depth_AD( k ) + total_od_AD
      ! Note: No total_od_AD = ZERO here because
      !       total_od_TL = total_od_TL + (....)

    END DO k_layer_tau_loop


    !#--------------------------------------------------------------------------#
    !#                        -- LOOP OVER ABSORBERS --                         #
    !#--------------------------------------------------------------------------#

    j_absorber_loop: DO j = 1, n_absorbers


      ! -----------------------------------------
      ! Check if there is any absorption for this
      ! absorber/channel combination.
      !
      ! This check is the reason why all channels
      ! cannot be processed at once and why the
      ! layer loop is within the absorber loop.
      ! -----------------------------------------

      IF ( predictor_index( 0, l, j ) == 0 ) CYCLE j_absorber_loop



      !#------------------------------------------------------------------------#
      !#                        -- LOOP OVER LAYERS --                          #
      !#------------------------------------------------------------------------#

      k_layer_od_loop: DO k = n_layers, 1, -1


        ! -----------------------------------
        ! Calculate the current layer average
        ! absorber amounts and differences
        ! -----------------------------------

        ave_absorber = POINT_5 * ( absorber( k, j ) + absorber( k-1, j ) )
        d_absorber   = absorber( k, j ) - absorber( k-1, j )


        ! -----------------------------------
        ! Assign absorber space layer indices
        ! -----------------------------------

        k2 = layer_index( k, j )
        k1 = k2 - 1



        !#----------------------------------------------------------------------#
        !#           -- HERE REPEAT THE FORWARD CALCULATION OF THE   --         #
        !#           -- ABSORPTION COEFFICIENT FOR THE CURRENT LAYER --         #
        !#----------------------------------------------------------------------#

        ! ---------------------------------------------------------
        ! To linearly interpolate the tau_coeffs to the actual user
        ! space absorber amount, need the gradient across the layer
        ! ---------------------------------------------------------

        d_absorber_space = absorber_space_levels( k2, j ) - absorber_space_levels( k1, j )

        gradient = ( absorber_space_levels( k2, j ) - ave_absorber ) / &
        !          -------------------------------------------------
                                   d_absorber_space


        ! --------------------------------
        ! Calculate absorption coeficients
        ! --------------------------------

        ! -- Offset term (save these values for later)
        b1o = tau_coefficients( 0, k1, l, j )
        b2o = tau_coefficients( 0, k2, l, j )

        absorption_coefficient = b2o + ( gradient * ( b1o - b2o ) )


        i_predictor_loop_FWD: DO i = 1, n_predictors

          ! -- Current predictor interpolated coefficient
          b1 = tau_coefficients( i, k1, l, j )
          b2 = tau_coefficients( i, k2, l, j )

          b = b2 + ( gradient * ( b1 - b2 ) )

          !  - Sum current predictor's contribution
          ip = predictor_index( i, l, j )
          absorption_coefficient = absorption_coefficient + ( b * predictor( ip, k ) ) 

        END DO i_predictor_loop_FWD


        ! -- **** IT WOULD BE NICE TO NOT NEED THIS AT ALL! ****
        absorption_coefficient = MAX( absorption_coefficient, ZERO )



        !#----------------------------------------------------------------------#
        !#                  -- BEGIN ADJOINT CALCULATIONS --                    #
        !#----------------------------------------------------------------------#

        ! ---------------------------------------------------------------
        ! Adjoints of the optical depth.
        !
        ! These quantities are local to the k_layer_od_loop
        ! and are equal to zero at this point so a straight
        ! initialisation is used, i.e. there is no
        !   d_adbsorber_AD            = d_absorber_AD + (...)
        !   absorption_coefficient_AD = absorption_coefficient_AD + (...)
        ! This also eliminates the need to zero out the two
        ! quanitities later in the loop once they no longer
        ! have an impact on the gradient vector result.
        !
        ! Also not that there is no
        !   optical_depth_AD( k ) = ZERO
        ! because
        !   optical_depth_TL( k ) = optical_depth_TL( k ) + (....)
        ! ---------------------------------------------------------------

        d_absorber_AD = absorption_coefficient * optical_depth_AD( k )   ! .... (1)
        absorption_coefficient_AD = d_absorber * optical_depth_AD( k )


        ! ---------------------------------
        ! Adjoint of absorption coefficient
        ! ---------------------------------

        ! -- Loop over predictors
        ! -- Note that b_AD is local to this loop only
        ! -- hence rather than
        ! --   b_AD = b_AD + (....)
        ! -- and, at end of loop,
        ! --   b_AD = ZERO
        ! -- it is simply initialised each iteration
        i_predictor_loop_AD: DO i = n_predictors, 1, -1

          ! -- Current predictor interpolated coefficient
          b1 = tau_coefficients( i, k1, l, j )
          b2 = tau_coefficients( i, k2, l, j )

          b = b2 + ( gradient * ( b1 - b2 ) )

          ! -- Adjoints of the absorption coefficient
          ip = predictor_index( i, l, j )

          predictor_AD( ip, k ) = predictor_AD( ip, k ) + ( b * absorption_coefficient_AD )
          b_AD = predictor( ip, k ) * absorption_coefficient_AD
          ! NOTE: No absorption_coefficient_AD = ZERO here because
          !       absorption_coefficient_TL = absorption_coefficient_TL + (....)

          ! -- Coefficient adjoint
          gradient_AD = gradient_AD + ( ( b1 - b2 ) * b_AD )

        END DO i_predictor_loop_AD


        ! -- Offset term.
        ! -- Note that absorption_coefficient_AD is not set to zero
        ! -- after this statement due to its initialisation each layer.
        gradient_AD = gradient_AD + ( ( b1o - b2o ) * absorption_coefficient_AD )


        ! -------------------------------------------------
        ! Adjoint of coefficient layer gradient.
        !
        ! The ave_absorber_AD quantity is also local to the
        ! k_layer_od_loop so a straight initialisation can
        ! be used rather than,
        !   ave_absorber_AD = ave_absorber_AD - (....)
        ! -------------------------------------------------

        ave_absorber_AD = -( gradient_AD / d_absorber_space )   ! .... (2)
        gradient_AD = ZERO


        ! ---------------------------------------------------
        ! Adjoints of the current layer average
        ! absorber amount and difference.
        !
        ! Neither d_absorber_AD nor ave_absorber_AD need
        ! to be set to zero after this as they are explicitly
        ! reassigned each layer iteration at (1) and (2) above
        ! respectively.
        ! ---------------------------------------------------

        absorber_AD( k-1, j ) = absorber_AD( k-1, j ) - d_absorber_AD + ( POINT_5 * ave_absorber_AD )
        absorber_AD( k,   j ) = absorber_AD( k,   j ) + d_absorber_AD + ( POINT_5 * ave_absorber_AD )

      END DO k_layer_od_loop

    END DO j_absorber_loop

  END SUBROUTINE compute_transmittance_AD

END MODULE transmittance


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
! Revision 1.8  2001/08/16 17:19:11  paulv
! - Updated documentation
!
! Revision 1.7  2001/08/01 17:04:00  paulv
! - The absorber space levels are no longer calculated during model
!   initialisation, but are precalculated and stored in the transmittance
!   coefficient data file. This means that,
!     USE absorber_space, ONLY : absorber_space_levels
!   was deleted as the absorber space level array is now available from
!   the TRANSMITTANCE_COEFFICIENTS module.
!
! Revision 1.6  2001/07/12 18:38:28  paulv
! - Use of ABSORBER_SPACE module now includes an ONLY clause so that the
!   absorber_space_levels is all that is available.
! - Direction specification changed from
!     IF ( direction == 0 ) THEN
!       ...do DOWNWELLING STUFF...
!     ELSE
!       ...do UPWELLING STUFF...
!     END IF
!   to
!     IF ( direction == UP ) THEN
!       ...do UPWELLING STUFF...
!     ELSE
!       ...do DOWNWELLING STUFF...
!     END IF
!   since the upwelling case is required for every call, but the downwelling
!   may not be. Also, the parameter UP is now used in the IF rather than
!   an actual number (0 in this case).
! - Changed
!      od_tolerance = ABS( ALOG( TOLERANCE ) )
!   to
!      od_tolerance = ABS( LOG( TOLERANCE ) )
! - Corrected bug in the forward calculation of the absorption coefficient
!   in TRANSMITTANCE_AD. The offset coefficients are defined as
!     b1o = tau_coefficients( 0, k1, l, j )
!     b2o = tau_coefficients( 0, k2, l, j )
!   and the offset term was initialised as
!     absorption_coefficient = b2 + ( gradient * ( b1o - b2o ) )
!   instead of
!     absorption_coefficient = b2o + ( gradient * ( b1o - b2o ) )
!   where in the former, B2 was specified rather than B2O
!
! Revision 1.5  2001/05/29 18:00:08  paulv
! - Added adjoint form of the transmittance calculation.
! - Removed the FIND_ABSORBER_SPACE_LAYER  routine. Now resides in the
!   ABSORBER_PROFILE module. The absorber space bracket layer indices are
!   now passed as arguments from the calling routine.
! - The predictor indices and transmittance coefficients are no longer passed
!   as arguments but read from the transmittance_coefficients module.
!
! Revision 1.4  2000/11/14 18:42:32  paulv
! - Merged branch incorporating tangent-linear code into main truck.
!   Optical depth debug code still present.
!
! Revision 1.3.1.2  2000/11/14 18:34:56  paulv
! - Finished adding tangent-linear code. Optical depth debug code still
!   present - output sent to unit numbers 51 and 61.
!
! Revision 1.3.1.1  2000/11/09 20:49:35  paulv
! - Adding tangent linear forms of the optical depth and transmittance
!   computation. IN PROGRESS AND INCOMPLETE.
! - Removed code that finds the absorber space bracket layers into its own
!   subroutine. Both the forward and tangent linear routines use the same
!   search method.
!
! Revision 1.3  2000/08/31 19:36:33  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 1.2  2000/08/24 16:55:33  paulv
! - Added optional NO_STANDARD input argument to the COMPUTE_TRANSMITTANCE
!   subprogram to prevent the (angle independent) standard predictors from being
!   recalculated when only the path angle has changed in the calling procedure.
! - The  profile data integration has been removed from the COMPUTE_TRANSMITTANCE
!   subprogram and is now performed outside of this module in the ABSORBER_PROFILE
!   module. This has a number of consequences:
!   o The VIEW_ANGLE input argument was removed and replaced with the path-angle
!     scaled ABSORBER_AMOUNTS argument.
!   o The interface pressure and ozone profile data are no longer required
!     and have been removed from the input argument list.
! - All profile integration and predictor calculation code has been removed
!   from the COMPUTE_TRANSMITTANCE subprogram.
! - Updated module and subprogram documentation.
!
! Revision 1.1  2000/08/22 15:57:26  paulv
! Initial checkin.
!
!
!
!
