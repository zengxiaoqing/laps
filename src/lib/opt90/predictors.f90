!------------------------------------------------------------------------------
!M+
! NAME:
!       predictors
!
! PURPOSE:
!       RT model predictor module
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE predictors
!
! OUTPUTS:
!       None.
!
! MODULES:
!       type_kinds:     Module containing data type kind definitions.
!
!       parameters:     Module containing parameter definitions for the
!                       RT model.
!
!       error_handler:  Module to define error codes and handle error
!                       conditions
!
! CONTAINS:
!       compute_predictors:             PUBLIC subroutine to calculate all
!                                       the predictors. This routine calls the
!                                       private functions that follow.
!
!       compute_std_predictors:         PRIVATE function to compute the standard,
!                                       absorber independent predictor set.
!
!       compute_int_predictors:         PRIVATE function to compute the integrated,
!                                       absorber dependent predictor set.
!
!       compute_predictors_TL:          PUBLIC subroutine to calculate all
!                                       tangent-linear predictors. This routine
!                                       calls the private functions that follow.
!
!       compute_std_predictors_TL:      PRIVATE function to compute the standard,
!                                       absorber independent tangent-linear 
!                                       predictor set.
!
!       compute_int_predictors_TL:      PRIVATE function to compute the integrated,
!                                       absorber dependent tangent-linear predictor
!                                       set.
!
!       compute_predictors_AD:          PUBLIC subroutine to calculate the
!                                       predictor adjoints. This routine calls
!                                       the functions that follow.
!
!       compute_std_predictors_AD:      PRIVATE function to compute the adjoint of
!                                       the standard, absorber independent predictor
!                                       set.
!
!       compute_int_predictors_AD:      PRIVATE function to compute the adjoint of
!                                       integrated, absorber dependent predictor set.
!
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
!       Written by:     Paul van Delst, CIMSS@NOAA/NCEP 11-Jul-2000
!                       pvandelst@ncep.noaa.gov
!
!       Adapted from code written by: Thomas J.Kleespies
!                                     NOAA/NESDIS/ORA
!                                     tkleespies@nesdis.noaa.gov
!
!                                     and modified by:
!
!                                     John Derber
!                                     NOAA/NCEP/EMC
!                                     jderber@ncep.noaa.gov
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

MODULE predictors


  ! ---------------------
  ! Module use statements
  ! ---------------------

  USE type_kinds, ONLY : fp_kind
  USE parameters
  USE error_handler


  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------

  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------

  PRIVATE
  PUBLIC :: compute_predictors
  PUBLIC :: compute_predictors_TL
  PUBLIC :: compute_predictors_AD


CONTAINS



!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_predictors
!
! PURPOSE:
!       PUBLIC routine to calculate the forward transmittance model predictors.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_predictors ( pressure,    &  ! Input,  K
!                                 temperature, &  ! Input,  K
!                                 water_vapor, &  ! Input,  K
!                                 absorber,    &  ! Input,  0:K x J
!
!                                 predictor,   &  ! Output, I x K
!
!                                 no_standard  )  ! Optional input
!
! INPUT ARGUMENTS:
!       pressure:         Profile LAYER average pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       temperature:      Profile LAYER average temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       water_vapor:      Profile LAYER average water vapor mixing ratio array.
!                         UNITS:      g/kg
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       absorber:         Profile LEVEL integrated absorber amount array.
!                         UNITS:      Varies with absorber.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  0:K x J
!                         ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       no_standard:      If present, the standard predictors are not calculated.
!                         This prevents recalculation of the standard predictors
!                         is only the view angle has changed - which only affects
!                         the integrated predictors.
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       predictor:        Profile LAYER predictors array.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  I x K
!                         ATTRIBUTES: INTENT( OUT ), TARGET
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! CALLS:
!       compute_std_predictors:   PRIVATE function to compute the standard
!                                 (absorber independent) predictor set.
!
!       compute_int_predictors:   PRIVATE function to compute the absorber
!                                 integrated predictor set.
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
!       The predictors used in the NCEP transmittance model are organised in
!       the following manner:
!
!       --------------------------------------------
!       | 1 | 2 | 3 | ... | 9 | 10 | 11 | ... | 27 |
!       --------------------------------------------
!
!       \                    / \                   /
!        \                  /   \                 /
!         ------------------     -----------------
!                  |                      |
!                  v                      v
!
!              Standard               Integrated
!             Predictors            Predictors for
!                                   each absorber
!
!       Pointers are used to reference the module scope predictor data array
!       before the calls are made to the standard and integrated predictor
!       calculation functions. This eliminates array copying.
!S-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_predictors ( pressure,    &  ! Input,  K
                                  temperature, &  ! Input,  K
                                  water_vapor, &  ! Input,  K
                                  absorber,    &  ! Input,  0:K x J

                                  predictor,   &  ! Output, I x K

                                  no_standard  )  ! Optional input



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )           :: pressure      ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )           :: temperature   ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )           :: water_vapor   ! K
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN )           :: absorber      ! 0:K x J

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ), TARGET  :: predictor     ! I x K

    INTEGER,                             INTENT( IN ), OPTIONAL :: no_standard


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_PREDICTORS'


    ! ---------------
    ! Local variables
    ! ---------------

    CHARACTER( 80 ) :: message

    INTEGER :: i1, i2, j

    REAL( fp_kind ), POINTER, DIMENSION( :, : ) :: std_predictors   ! Istd x K
    REAL( fp_kind ), POINTER, DIMENSION( :, : ) :: int_predictors   ! Iint x K



    !#--------------------------------------------------------------------------#
    !#            -- Calculate the standard predictors if needed --             #
    !#--------------------------------------------------------------------------#

    IF ( .NOT. PRESENT( no_standard ) ) THEN

      ! -- Alias the input predictor array
      std_predictors => predictor( 1:MAX_N_STANDARD_PREDICTORS, : )

      CALL compute_std_predictors( pressure,      &
                                   temperature,   &
                                   water_vapor,   &
                                   std_predictors )

    END IF


                                                   
    !#--------------------------------------------------------------------------#
    !#                -- Calculate the integrated predictors --                 #
    !#--------------------------------------------------------------------------#

    j_absorber_loop: DO j = 1, SIZE( absorber, DIM = 2 )

      ! -- Determine indices of current absorber predictors
      i1 = MAX_N_STANDARD_PREDICTORS + ( ( j - 1 ) * MAX_N_INTEGRATED_PREDICTORS ) + 1
      i2 = i1 + MAX_N_INTEGRATED_PREDICTORS - 1

      ! -- Alias the input predictor array
      int_predictors => predictor( i1:i2, : )

      ! -- Compute the predictors for the current absorber
      CALL compute_int_predictors( pressure,          &
                                   temperature,       &
                                   absorber( 0:, j ), &
                                   int_predictors     )

    END DO j_absorber_loop

  END SUBROUTINE compute_predictors



!--------------------------------------------------------------------------------
!P+
! NAME:
!       compute_std_predictors
!
! PURPOSE:
!       PRIVATE function to compute the standard, absorber independent
!       predictor set for the forward transmittance model.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_std_predictors( pressure,    &
!                                    temperature, &
!                                    water_vapor, &
!                                    predictor    )
!
! INPUT ARGUMENTS:
!       pressure:         Profile LAYER average pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       temperature:      Profile LAYER average temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       water_vapor:      Profile LAYER average water vapor mixing ratio array.
!                         UNITS:      g/kg
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None
!
! OUTPUT ARGUMENTS:
!       predictor:        Array containing the calculated standard predictor set.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  Istd x K
!                         ATTRIBUTES: INTENT( OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None
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
!       McMillin, L.M., L.J. Crone, M.D. Goldberg, and T.J. Kleespies,
!         "Atmospheric transmittance of an absorbing gas. 4. OPTRAN: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", Applied Optics, 1995, v34, pp6269-6274.
!
!       The standard predictors are the following:
!
!         1) Temperature, T
!         2) Pressure, P
!         3) T^2
!         4) P^2
!         5) T.P
!         6) T^2.P
!         7) T.P^2
!         8) T^2.P^2
!         9) Water vapor mixing ratio, Q
!P-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_std_predictors( p,        &  ! Input,  K
                                     t,        &  ! Input,  K
                                     w,        &  ! Input,  K
                                     predictor )  ! Output, Istd x K



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: p          ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: t          ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: w          ! K

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: predictor  ! Istd x K


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_STD_PREDICTORS'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: k

    REAL( fp_kind ) :: p2
    REAL( fp_kind ) :: t2



    !#--------------------------------------------------------------------------#
    !#                -- Calculate the standard predictor set --                #
    !#--------------------------------------------------------------------------#

    k_layer_loop: DO k = 1, SIZE( p )


      ! -- Precalculate the squared terms
      p2 = p( k ) * p( k )
      t2 = t( k ) * t( k )

      ! -- Calculate and assign the absorber independent predictors
      predictor( 1, k ) = t( k )
      predictor( 2, k ) = p( k )
      predictor( 3, k ) = t2
      predictor( 4, k ) = p2
      predictor( 5, k ) = t( k ) * p( k )
      predictor( 6, k ) = t2     * p( k )
      predictor( 7, k ) = t( k ) * p2
      predictor( 8, k ) = t2     * p2
      predictor( 9, k ) = w( k )

    END DO k_layer_loop

  END SUBROUTINE compute_std_predictors






!--------------------------------------------------------------------------------
!P+
! NAME:
!       compute_int_predictors
!
! PURPOSE:
!       PRIVATE function to compute the integrated, absorber dependent
!       predictor set for the forward transmittance model.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_int_predictors( pressure,    &
!                                    temperature, &
!                                    absorber,    &
!                                    predictor    )
!
! INPUT ARGUMENTS:
!       pressure:         Profile LAYER average pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       temperature:      Profile LAYER average temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       absorber:         Profile LEVEL integrated absorber amount array.
!                         UNITS:      Varies with absorber
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  0:K
!                         ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None
!
! OUTPUT ARGUMENTS:
!       predictor:        Array containing the calculated integrated predictor set
!                         for the passed absorber.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  Iint x K
!                         ATTRIBUTES: INTENT( OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None
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
!       McMillin, L.M., L.J. Crone, M.D. Goldberg, and T.J. Kleespies,
!         "Atmospheric transmittance of an absorbing gas. 4. OPTRAN: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", Applied Optics, 1995, v34, pp6269-6274.
!
!       The integrated predictors consist of six that are repeated for every
!       absorber:
!
!         1) T*
!         2) P*
!         3) T**
!         4) P**
!         5) T***
!         6) P***
!
!       The reference above details the predictor value for input LEVEL values. In
!       the equations below, the "X" quantities can be replaced with temperature, T, 
!       of pressure, P, to determine the required predictor:
!
!                           __ k
!                    1     \ 
!         X*(k) = --------  >  ( X(i) + X(i-1) ) ( A(i) - A(i-1) )
!                  c.A(k)  /__
!                             i=1
!
!                            __ k
!                     1     \
!        X**(k) = ---------  >  ( X(i)A(i) + X(i-1)A(i-1) ) ( A(i) - A(i-1) )
!                  c.A^2(k) /__
!                              i=1
!
!                            __ k
!                     1     \ 
!       X***(k) = ---------  >  ( X(i)A^2(i) + X(i-1)A^2(i-1) ) ( A(i) - A(i-1) )
!                  c.A^3(k) /__
!                              i=1
!
!       To accomodate input LAYER values (i.e. layer average values) from the
!       NCEP GDAS, the predictor formulations were modified to:
!
!                   __ k                        __k-1
!                  \                           \ 
!                   >  X(i)( A(i)-A(i-1) )      >  X(i)( A(i)-A(i-1) )
!                  /__                         /__ 
!                     i=1                         i=1
!         X*(k) = ------------------------- + -------------------------
!                            2.A(k)                    2.A(k-1)
!
!                   __ k
!                  \ 
!                   >  X(i)( A(i) + A(i-1) ) ( A(i) - A(i-1) )
!                  /__
!                     i=1
!        X**(k) = --------------------------------------------- +
!                                       2.A^2(k) 
!                   
!                   __k-1
!                  \ 
!                   >  X(i)( A(i) + A(i-1) ) ( A(i) - A(i-1) )
!                  /__
!                     i=1
!                 ---------------------------------------------
!                                       2.A^2(k-1)
!
!                      __ k
!                     \ 
!                 3 .  >  X(i)( A^2(i) + A^2(i-1) ) ( A(i) - A(i-1) )
!                     /__
!                        i=1
!       X***(k) = ---------------------------------------------------- +
!                                       4.A^3(k)                    
!
!                      __k-1
!                     \ 
!                 3 .  >  X(i)( A^2(i) + A^2(i-1) ) ( A(i) - A(i-1) )
!                     /__
!                        i=1
!                 ----------------------------------------------------
!                                       4.A^3(k-1)
!
!       Thus the transmittance model coefficients calculated using the LEVEL
!       predictor formulation are used with the predictors constructed above 
!       with LAYER values.
!P-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_int_predictors( pressure,    &  ! Input,  K
                                     temperature, &  ! Input,  K
                                     absorber,    &  ! Input,  0:K
                                     predictor    )  ! Output, Iint x K



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: pressure     ! K
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: temperature  ! K
    REAL( fp_kind ), DIMENSION( 0: ),   INTENT( IN )  :: absorber     ! 0:K

    REAL( fp_kind ), DIMENSION( :, : ), INTENT( OUT ) :: predictor    ! Iint x K


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_INT_PREDICTORS'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: i, k
    INTEGER :: n_layers
    INTEGER :: n_predictors

    REAL( fp_kind ) :: d_absorber
    REAL( fp_kind ) :: factor_1
    REAL( fp_kind ) :: factor_2
    REAL( fp_kind ) :: inverse_1
    REAL( fp_kind ) :: inverse_2
    REAL( fp_kind ) :: inverse_3
    REAL( fp_kind ) :: absorber_3


    ! -- Square of the absorber amount. 0:K
    REAL( fp_kind ), DIMENSION( 0:SIZE( pressure ) ) :: absorber_2

    ! -- Intermediate summation array. Iint
    REAL( fp_kind ), DIMENSION( SIZE( predictor, DIM=1 ) ) :: s

    ! -- LEVEL predictor, Iint x 0:K
    REAL( fp_kind ), DIMENSION( SIZE( predictor, DIM=1 ), 0:SIZE( pressure ) ) :: x



    !#--------------------------------------------------------------------------#
    !#                  -- Determine the number predictors --                   #
    !#--------------------------------------------------------------------------#

    n_predictors = SIZE( predictor, DIM = 1 )



    !#--------------------------------------------------------------------------#
    !#                         -- Initialise values --                          #
    !#--------------------------------------------------------------------------#

    absorber_2( 0 ) = ZERO

    s( : )    = ZERO
    x( :, 0 ) = ZERO



    !#--------------------------------------------------------------------------#
    !#               -- Calculate the integrated predictor set --               #
    !#--------------------------------------------------------------------------#

    k_layer_loop: DO k = 1, SIZE( pressure )


      ! -----------------------------------------
      ! Calculate absorber multiplicative factors
      ! -----------------------------------------

      absorber_2( k ) = absorber( k ) * absorber( k )

      d_absorber = absorber( k ) - absorber( k-1 )                      ! For * terms
      factor_1   = ( absorber( k )   + absorber( k-1 )   ) * d_absorber ! For ** terms
      factor_2   = ( absorber_2( k ) + absorber_2( k-1 ) ) * d_absorber ! For *** terms


      ! -------------------------------
      ! Calculate the intermediate sums
      ! -------------------------------

      s( 1 ) = s( 1 ) + ( temperature( k ) * d_absorber )  ! T*
      s( 2 ) = s( 2 ) + ( pressure( k )    * d_absorber )  ! P*

      s( 3 ) = s( 3 ) + ( temperature( k ) * factor_1 )    ! T**
      s( 4 ) = s( 4 ) + ( pressure( k )    * factor_1 )    ! P**

      s( 5 ) = s( 5 ) + ( temperature( k ) * factor_2 )    ! T***
      s( 6 ) = s( 6 ) + ( pressure( k )    * factor_2 )    ! P***


      ! -------------------------------------------------------
      ! Calculate the normalising factors for the integrated
      ! predictors. Note that the checks below, the IF tests to
      ! determine if the absorber products are represenatble
      ! are to minimise the number of calcs. I.e if inverse_1
      ! is toast because absorber(k) is too small there's no
      ! need to check any further.
      ! -------------------------------------------------------

      ! -- Is inverse_1 representable?
      inverse_1_check: IF ( absorber( k ) > TOLERANCE ) THEN

        inverse_1 = ONE / absorber( k )

        ! -- Is inverse_2 representable
        inverse_2_check: IF ( absorber_2( k ) > TOLERANCE ) THEN

          inverse_2  = inverse_1 * inverse_1
          absorber_3 = absorber( k ) * absorber_2( k )
         
          ! -- Is inverse_3 representable
          inverse_3_check: IF ( absorber_3 > TOLERANCE ) THEN

            inverse_3 = inverse_2 * inverse_1

          ELSE

            inverse_3 = ZERO

          END IF inverse_3_check

        ELSE

          inverse_2 = ZERO
          inverse_3 = ZERO

        END IF inverse_2_check

      ELSE

        inverse_1 = ZERO
        inverse_2 = ZERO
        inverse_3 = ZERO

      END IF inverse_1_check


      ! ---------------------------------------------
      ! Scale and normalise the integrated predictors
      ! ---------------------------------------------

      x( 1, k ) = POINT_5  * s( 1 ) * inverse_1  ! T*
      x( 2, k ) = POINT_5  * s( 2 ) * inverse_1  ! P*

      x( 3, k ) = POINT_5  * s( 3 ) * inverse_2  ! T**
      x( 4, k ) = POINT_5  * s( 4 ) * inverse_2  ! P**

      x( 5, k ) = POINT_75 * s( 5 ) * inverse_3  ! T***
      x( 6, k ) = POINT_75 * s( 6 ) * inverse_3  ! P***


      ! ----------------------------
      ! Sum predictors across layers
      ! ----------------------------

      DO i = 1, n_predictors

        predictor( i, k ) = x( i, k ) + x( i, k-1 )

      END DO

    END DO k_layer_loop

  END SUBROUTINE compute_int_predictors






!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_predictors_TL
!
! PURPOSE:
!       PUBLIC routine to calculate the tangent-linear transmittance model
!       predictors.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_predictors_TL ( pressure,       &  ! Input,  K
!                                    temperature,    &  ! Input,  K
!                                    water_vapor,    &  ! Input,  K
!                                    absorber,       &  ! Input,  0:K x J
!
!                                    pressure_TL,    &  ! Input,  K
!                                    temperature_TL, &  ! Input,  K
!                                    water_vapor_TL, &  ! Input,  K
!                                    absorber_TL,    &  ! Input,  0:K x J
!
!                                    predictor_TL,   &  ! Output, I x K
!
!                                    no_standard     )  ! Optional input
!
! INPUT ARGUMENTS:
!       pressure:            Profile LAYER average pressure array.
!                            UNITS:      hPa
!                            TYPE:       REAL( fp_kind )
!                            DIMENSION:  K
!                            ATTRIBUTES: INTENT( IN )
!
!       temperature:         Profile LAYER average temperature array.
!                            UNITS:      Kelvin
!                            TYPE:       REAL( fp_kind )
!                            DIMENSION:  K
!                            ATTRIBUTES: INTENT( IN )
!
!       water_vapor:         Profile LAYER average water vapor mixing ratio array.
!                            UNITS:      g/kg
!                            TYPE:       REAL( fp_kind )
!                            DIMENSION:  K
!                            ATTRIBUTES: INTENT( IN )
!
!       absorber:            Profile LEVEL integrated absorber amount array.
!                            UNITS:      Varies with absorber
!                            TYPE:       REAL( fp_kind )
!                            DIMENSION:  K
!                            ATTRIBUTES: INTENT( IN )
!
!       pressure_TL:         Profile tangent-linear pressure array, i.e. the
!                            pressure perturbation.
!                            UNITS:      hPa
!                            TYPE:       REAL( fp_kind )
!                            DIMENSION:  K
!                            ATTRIBUTES: INTENT( IN )
!
!       temperature_TL:      Profile tangent-linear LAYER average temperature array
!                            i.e. the temperature perturbation.
!                            UNITS:      Kelvin
!                            TYPE:       REAL( fp_kind )
!                            DIMENSION:  K
!                            ATTRIBUTES: INTENT( IN )
!
!       water_vapor_TL:      Profile tangent-linear water vapor mixing
!                            ratio array, i.e. the water vapor mixing
!                            ratio perturbation.
!                            UNITS:      g/kg
!                            TYPE:       REAL( fp_kind )
!                            DIMENSION:  K
!                            ATTRIBUTES: INTENT( IN )
!
!       absorber_TL:         Profile LEVEL integrated tangent-linear 
!                            absorber amount array.
!                            UNITS:      Varies with absorber
!                            TYPE:       REAL( fp_kind )
!                            DIMENSION:  K
!                            ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       no_standard:         If present, the standard predictors are not calculated.
!                            This prevents recalculation of the standard predictors
!                            is only the view angle has changed - which only affects
!                            the integrated predictors.
!                            UNITS:      None
!                            TYPE:       Integer
!                            DIMENSION:  Scalar
!                            ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       predictor_TL:        Profile LAYER tangent-linear predictors array.
!                            UNITS:      Varies with predictor type.
!                            TYPE:       REAL( fp_kind )
!                            DIMENSION:  Iint x K
!                            ATTRIBUTES: INTENT( OUT ), TARGET
!
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! CALLS:
!       compute_std_predictors_TL:   PRIVATE function to compute the tangent-
!                                    linear form of the standard (absorber
!                                    independent) predictor set.
!
!       compute_int_predictors_TL:   PRIVATE function to compute the tangent-
!                                    linear form of the absorber integrated
!                                    predictor set.
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
!       The predictors used in the NCEP transmittance model are organised in
!       the following manner:
!
!       --------------------------------------------
!       | 1 | 2 | 3 | ... | 9 | 10 | 11 | ... | 27 |
!       --------------------------------------------
!
!       \                    / \                   /
!        \                  /   \                 /
!         ------------------     -----------------
!                  |                      |
!                  v                      v
!
!              Standard               Integrated
!             Predictors            Predictors for
!                                   each absorber
!
!       Pointers are used to reference the module scope predictor data array
!       before the calls are made to the standard and integrated predictor
!       calculation functions. This eliminates array copying.
!S-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_predictors_TL ( pressure,       &  ! Input,  K
                                     temperature,    &  ! Input,  K
                                     water_vapor,    &  ! Input,  K
                                     absorber,       &  ! Input,  0:K x J

                                     pressure_TL,    &  ! Input,  K
                                     temperature_TL, &  ! Input,  K
                                     water_vapor_TL, &  ! Input,  K
                                     absorber_TL,    &  ! Input,  0:K x J

                                     predictor_TL,   &  ! Output, I x K

                                     no_standard     )  ! Optional input



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )           :: pressure        ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )           :: temperature     ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )           :: water_vapor     ! K
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN )           :: absorber        ! 0:K x J

    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )           :: pressure_TL     ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )           :: temperature_TL  ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )           :: water_vapor_TL  ! K
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN )           :: absorber_TL     ! 0:K x J

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ), TARGET  :: predictor_TL    ! I x K

    INTEGER,                             INTENT( IN ), OPTIONAL :: no_standard


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_PREDICTORS_TL'


    ! ---------------
    ! Local variables
    ! ---------------

    CHARACTER( 80 ) :: message

    INTEGER :: i1, i2, j

    REAL( fp_kind ), POINTER, DIMENSION( :, : ) :: std_predictors_TL   ! Istd x K
    REAL( fp_kind ), POINTER, DIMENSION( :, : ) :: int_predictors_TL   ! Iint x K




    !#--------------------------------------------------------------------------#
    !#     -- Calculate the tangent-linear standard predictors if needed --     #
    !#--------------------------------------------------------------------------#

    IF ( .NOT. PRESENT( no_standard ) ) THEN

      ! -- Alias the input predictor array
      std_predictors_TL => predictor_TL( 1:MAX_N_STANDARD_PREDICTORS, : )

      CALL compute_std_predictors_TL( pressure,         &
                                      temperature,      &
                                      water_vapor,      &

                                      pressure_TL,      &
                                      temperature_TL,   &
                                      water_vapor_TL,   &

                                      std_predictors_TL )

    END IF


                                                   
    !#--------------------------------------------------------------------------#
    !#                -- Calculate the integrated predictors --                 #
    !#--------------------------------------------------------------------------#

    j_absorber_loop: DO j = 1, SIZE( absorber_TL, DIM = 2 )

      ! -- Determine indices of current absorber predictors
      i1 = MAX_N_STANDARD_PREDICTORS + ( ( j - 1 ) * MAX_N_INTEGRATED_PREDICTORS ) + 1
      i2 = i1 + MAX_N_INTEGRATED_PREDICTORS - 1

      ! -- Alias the input predictor array
      int_predictors_TL => predictor_TL( i1:i2, : )

      ! -- Calculate tangent-linear predictors for current absorber
      CALL compute_int_predictors_TL( pressure,             &
                                      temperature,          &
                                      absorber( 0:, j ),    &

                                      pressure_TL,          &
                                      temperature_TL,       &
                                      absorber_TL( 0:, j ), &

                                      int_predictors_TL     )

    END DO j_absorber_loop

  END SUBROUTINE compute_predictors_TL





!--------------------------------------------------------------------------------
!P+
! NAME:
!       compute_std_predictors_TL
!
! PURPOSE:
!       PRIVATE function to compute the standard, absorber independent
!       predictor set for the tangent-linear transmittance model.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_std_predictors_TL( pressure,       &  ! Input,  K
!                                       temperature,    &  ! Input,  K
!                                       water_vapor,    &  ! Input,  K
!
!                                       pressure_TL,    &  ! Input,  K
!                                       temperature_TL, &  ! Input,  K
!                                       water_vapor_TL, &  ! Input,  K
!
!                                       predictors_TL   )  ! Output, Istd x K
!
! INPUT ARGUMENTS:
!       pressure:         Profile LAYER average pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       temperature:      Profile LAYER average temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       water_vapor:      Profile LAYER average water vapor mixing ratio array.
!                         UNITS:      g/kg
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       pressure_TL:      Profile LAYER tangent-linear pressure array, i.e. the
!                         pressure perturbation.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       temperature_TL:   Profile LAYER tangent-linear LAYER average temperature array
!                         i.e. the temperature perturbation.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       water_vapor_TL:   Profile LAYER tangent-linear water vapor mixing
!                         ratio array, i.e. the water vapor mixing
!                         ratio perturbation.
!                         UNITS:      g/kg
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None
!
! OUTPUT ARGUMENTS:
!       predictor_TL:     Array containing the calculated tangent-linear
!                         standard predictor set.
!                         UNITS:      Varies depending on predictor index.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  Istd x K
!                         ATTRIBUTES: INTENT( OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None
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
!       McMillin, L.M., L.J. Crone, M.D. Goldberg, and T.J. Kleespies,
!         "Atmospheric transmittance of an absorbing gas. 4. OPTRAN: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", Applied Optics, 1995, v34, pp6269-6274.
!
!       The standard predictors are the following:
!
!         1) Temperature, T
!         2) Pressure, P
!         3) T^2
!         4) P^2
!         5) T.P
!         6) T^2.P
!         7) T.P^2
!         8) T^2.P^2
!         9) Water vapor mixing ratio, Q
!
!       Thus the tangent-linear form of these are
!
!         1) dT
!         2) dP
!         3) 2T.dT
!         4) 2P.dP
!         5) P.dT + T.dP
!         6) 2TP.dT + T^2.dP
!         7) 2TP.dP + P^2.dT
!         8) 2T(P^2).dT + 2(T^2)P.dP
!         9) dQ
!
!P-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_std_predictors_TL( p,           &  ! Input,  K
                                        t,           &  ! Input,  K
                                        w,           &  ! Input,  K

                                        p_TL,        &  ! Input,  K
                                        t_TL,        &  ! Input,  K
                                        w_TL,        &  ! Input,  K

                                        predictor_TL  )  ! Output, Istd x K



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: p              ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: t             ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: w             ! Input,  K

    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: p_TL          ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: t_TL          ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: w_TL          ! Input,  K

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: predictor_TL  ! Output, Istd x K


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_STD_PREDICTORS_TL'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: k

    REAL( fp_kind ) :: p2, p2_TL
    REAL( fp_kind ) :: t2, t2_TL



    !#--------------------------------------------------------------------------#
    !#        -- Calculate the tangent-linear standard predictor set --         #
    !#--------------------------------------------------------------------------#

    k_layer_loop: DO k = 1, SIZE( p )


      ! -- Precalculate the squared terms
      p2 = p( k ) * p( k )
      t2 = t( k ) * t( k )

      ! -- Tangent-linear of squared terms
      p2_TL = TWO * p( k ) * p_TL( k )
      t2_TL = TWO * t( k ) * t_TL( k )
      
      ! -- Calculate and assign the absorber independent predictors
      predictor_TL( 1, k ) = t_TL( k )
      predictor_TL( 2, k ) = p_TL( k )
      predictor_TL( 3, k ) = t2_TL
      predictor_TL( 4, k ) = p2_TL
      predictor_TL( 5, k ) = ( t( k ) * p_TL( k ) ) + ( p( k ) * t_TL( k ) )
      predictor_TL( 6, k ) = ( p( k ) * t2_TL     ) + ( t2     * p_TL( k ) )
      predictor_TL( 7, k ) = ( t( k ) * p2_TL     ) + ( p2     * t_TL( k ) )
      predictor_TL( 8, k ) = ( t2     * p2_TL     ) + ( p2     * t2_TL     )
      predictor_TL( 9, k ) = w_TL( k )

    END DO k_layer_loop

  END SUBROUTINE compute_std_predictors_TL





!--------------------------------------------------------------------------------
!P+
! NAME:
!       compute_int_predictors_TL
!
! PURPOSE:
!       PRIVATE function to compute the integrated, absorber dependent
!       predictor set for the tangent-linear transmittance model.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_int_predictors_TL( pressure,       &  ! Input, K
!                                       temperature,    &  ! Input, K
!                                       absorber,       &  ! Input, 0:K
!
!                                       pressure_TL,    &  ! Input, K
!                                       temperature_TL, &  ! Input, K
!                                       absorber_TL,    &  ! Input, 0:K
!
!                                       predictor_TL    )  ! Input, Iint x K
!
! INPUT ARGUMENTS:
!       pressure:         Profile LAYER pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       temperature:      Profile LAYER temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       absorber_TL:      Profile LEVEL integrated absorber amount array.
!                         UNITS:      Varies with absorber
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       pressure_TL:      Profile LAYER tangent-linear pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       temperature_TL:   Profile LAYER tangent-linear temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       absorber_TL:      Profile LEVEL tangent-linear integrated absorber
!                         amount array.
!                         UNITS:      Varies with absorber
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None
!
! OUTPUT ARGUMENTS:
!       predictor_TL:     Array containing the calculated tangent-linear
!                         integrated predictor set for every absorber.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  Iint x K
!                         ATTRIBUTES: INTENT( OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None
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
!       McMillin, L.M., L.J. Crone, M.D. Goldberg, and T.J. Kleespies,
!         "Atmospheric transmittance of an absorbing gas. 4. OPTRAN: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", Applied Optics, 1995, v34, pp6269-6274.
!
!       The integrated predictors consist of six that are repeated for every
!       absorber:
!
!         1) T*
!         2) P*
!         3) T**
!         4) P**
!         5) T***
!         6) P***
!
!       The reference above details the predictor value for input LEVEL values. In
!       the equations below, the "X" quantities can be replaced with temperature, T, 
!       of pressure, P, to determine the required predictor, and A is the LEVEL
!       absorber amount:
!
!                           __ k
!                    1     \ 
!         X*(k) = --------  >  ( X(i) + X(i-1) ) ( A(i) - A(i-1) )
!                  c.A(k)  /__
!                             i=1
!
!                            __ k
!                     1     \
!        X**(k) = ---------  >  ( X(i)A(i) + X(i-1)A(i-1) ) ( A(i) - A(i-1) )
!                  c.A^2(k) /__
!                              i=1
!
!                            __ k
!                     1     \ 
!       X***(k) = ---------  >  ( X(i)A^2(i) + X(i-1)A^2(i-1) ) ( A(i) - A(i-1) )
!                  c.A^3(k) /__
!                              i=1
!
!       To accomodate input LAYER values (i.e. layer average values) from the
!       NCEP GDAS, the predictor formulations were modified to:
!
!                   __ k                        __k-1
!                  \                           \ 
!                   >  X(i)( A(i)-A(i-1) )      >  X(i)( A(i)-A(i-1) )
!                  /__                         /__ 
!                     i=1                         i=1
!         X*(k) = ------------------------- + -------------------------
!                            2.A(k)                    2.A(k-1)
!
!                   __ k
!                  \ 
!                   >  X(i)( A(i) + A(i-1) ) ( A(i) - A(i-1) )
!                  /__
!                     i=1
!        X**(k) = --------------------------------------------- +
!                                       2.A^2(k) 
!                   
!                   __k-1
!                  \ 
!                   >  X(i)( A(i) + A(i-1) ) ( A(i) - A(i-1) )
!                  /__
!                     i=1
!                 ---------------------------------------------
!                                       2.A^2(k-1)
!
!                      __ k
!                     \ 
!                 3 .  >  X(i)( A^2(i) + A^2(i-1) ) ( A(i) - A(i-1) )
!                     /__
!                        i=1
!       X***(k) = ---------------------------------------------------- +
!                                       4.A^3(k)                    
!
!                      __k-1
!                     \ 
!                 3 .  >  X(i)( A^2(i) + A^2(i-1) ) ( A(i) - A(i-1) )
!                     /__
!                        i=1
!                 ----------------------------------------------------
!                                       4.A^3(k-1)
!
!
!       The tangent-linear form of these LAYER predictor formulations are:
!
!                    __ k                                             
!                   \                                                 
!                    >  dX(i)( A(i)-A(i-1) ) + X(i)( dA(i)-dA(i-1) )  
!                   /__                                               
!                      i=1                                            
!         dX*(k) = -------------------------------------------------- -
!                                       2.A(k)                        
!
!                           __ k
!                          \
!                   dA(k) . >  X(i)( A(i)-A(i-1) )
!                          /__
!                             i=1
!                  -------------------------------- +
!                               2.A^2(k)
!                 
!                    __k-1                                              
!                   \                                                   
!                    >  dX(i)( A(i)-A(i-1) ) + X(i)( dA(i)-dA(i-1) )    
!                   /__                                                 
!                      i=1                                              
!                  -------------------------------------------------- - 
!                                      2.A(k-1)                         
!
!                             __k-1
!                            \
!                   dA(k-1) . >  X(i)( A(i)-A(i-1) )
!                            /__
!                               i=1
!                  ----------------------------------
!                              2.A^2(k-1)
!
!                 
!                 
!                    __ k                                                                                                       
!                   \                                                                                                           
!                    >  [dX(i)( A(i)+A(i-1) ) + X(i)( dA(i)+dA(i-1) )]( A(i)-A(i-1) ) + X(i)( A(i)+A(i-1) )( dA(i)-dA(i-1) )    
!                   /__                                                                                                         
!                      i=1                                                                                                      
!        dX**(k) = ---------------------------------------------------------------------------------------------------------- - 
!                                                                  2.A^2(k)                                                     
!
!                           __ k
!                          \
!                   dA(k) . >  X(i)( A(i)+A(i-1) )( A(i)-A(i-1) )
!                          /__
!                             i=1
!                  ----------------------------------------------- +
!                                      A^3(k)
!
!                    __k-1                                                                                                      
!                   \                                                                                                           
!                    >  [dX(i)( A(i)+A(i-1) ) + X(i)( dA(i)+dA(i-1) )]( A(i)-A(i-1) ) + X(i)( A(i)+A(i-1) )( dA(i)-dA(i-1) )    
!                   /__                                                                                                         
!                      i=1                                                                                                      
!                  ---------------------------------------------------------------------------------------------------------- - 
!                                                                 2.A^2(k-1)                                                    
!
!                             __k-1
!                            \
!                   dA(k-1) . >  X(i)( A(i)+A(i-1) )( A(i)-A(i-1) )
!                            /__
!                               i=1
!                  -------------------------------------------------
!                                      A^3(k-1)
!
!
!                 
!                      __ k                                                                                                       
!                     \            2    2                 2     2                                2    2                           
!                  3 . >  [dX(i)( A(i)+A(i-1) ) + X(i)( dA(i)+dA(i-1) )]( A(i)-A(i-1) ) + X(i)( A(i)+A(i-1) )( dA(i)-dA(i-1) )    
!                     /__                                                                                                         
!                        i=1                                                                                                      
!       dX***(k) = ------------------------------------------------------------------------------------------------------------ - 
!                                                                  4.A^3(k)                                                       
!
!                            __ k
!                           \          2    2
!                  9.dA(k) . >  X(i)( A(i)+A(i-1) )( A(i)-A(i-1) )
!                           /__
!                              i=1
!                  ------------------------------------------------ +
!                                      4.A^4(k)
!
!                      __k-1                                                                                                      
!                     \            2    2                 2     2                                2    2                           
!                  3 . >  [dX(i)( A(i)+A(i-1) ) + X(i)( dA(i)+dA(i-1) )]( A(i)-A(i-1) ) + X(i)( A(i)+A(i-1) )( dA(i)-dA(i-1) )    
!                     /__                                                                                                         
!                        i=1                                                                                                      
!                  ------------------------------------------------------------------------------------------------------------ - 
!                                                                 4.A^3(k-1)                                                      
!
!                            __k-1
!                           \          2    2
!                  9.dA(k) . >  X(i)( A(i)+A(i-1) )( A(i)-A(i-1) )
!                           /__
!                              i=1
!                  ------------------------------------------------
!                                      4.A^4(k-1)
!
!
!       Thus the transmittance model coefficients calculated using the LEVEL
!       predictor formulation are used with the tangent-linear predictors
!       constructed above with LAYER values.
!
!P-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_int_predictors_TL( pressure,       &  ! Input,  K
                                        temperature,    &  ! Input,  K
                                        absorber,       &  ! Input,  0:K

                                        pressure_TL,    &  ! Input,  K
                                        temperature_TL, &  ! Input,  K
                                        absorber_TL,    &  ! Input,  0:K

                                        predictor_TL    )  ! Output, Iint x K
 


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: pressure        ! K
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: temperature     ! K
    REAL( fp_kind ), DIMENSION( 0: ),   INTENT( IN )  :: absorber        ! 0:K

    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: pressure_TL     ! K
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: temperature_TL  ! K
    REAL( fp_kind ), DIMENSION( 0: ),   INTENT( IN )  :: absorber_TL     ! 0:K

    REAL( fp_kind ), DIMENSION( :, : ), INTENT( OUT ) :: predictor_TL    ! Iint x K


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_INT_PREDICTORS_TL'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: i, k
    INTEGER :: n_predictors

    REAL( fp_kind ) :: d_absorber
    REAL( fp_kind ) :: d_absorber_TL
    REAL( fp_kind ) :: factor_1
    REAL( fp_kind ) :: factor_1_TL
    REAL( fp_kind ) :: factor_2
    REAL( fp_kind ) :: factor_2_TL
    REAL( fp_kind ) :: inverse_1
    REAL( fp_kind ) :: inverse_2
    REAL( fp_kind ) :: inverse_3
    REAL( fp_kind ) :: inverse_4
    REAL( fp_kind ) :: absorber_3
    REAL( fp_kind ) :: absorber_4
    REAL( fp_kind ) :: inverse_1_TL
    REAL( fp_kind ) :: inverse_2_TL
    REAL( fp_kind ) :: inverse_3_TL


    ! -- Square of the absorber amount. 0:K
    REAL( fp_kind ), DIMENSION( 0:SIZE( pressure ) ) :: absorber_2

    ! -- Intermediate summation arrays. Iint
    REAL( fp_kind ), DIMENSION( SIZE( predictor_TL, DIM=1 ) ) :: s
    REAL( fp_kind ), DIMENSION( SIZE( predictor_TL, DIM=1 ) ) :: s_TL

    ! -- LEVEL predictor, Iint x 0:K
    REAL( fp_kind ), DIMENSION( SIZE( predictor_TL, DIM=1 ), 0:SIZE( pressure ) ) :: x_TL



    !#--------------------------------------------------------------------------#
    !#          -- Determine the number of layers and predictors --             #
    !#--------------------------------------------------------------------------#

    n_predictors = SIZE( predictor_TL, DIM = 1 )



    !#--------------------------------------------------------------------------#
    !#                         -- Initialise values --                          #
    !#--------------------------------------------------------------------------#

    absorber_2( 0 ) = ZERO

    s( : )       = ZERO
    s_TL( : )    = ZERO
    x_TL( :, 0 ) = ZERO



    !#--------------------------------------------------------------------------#
    !#               -- Calculate the integrated predictor set --               #
    !#--------------------------------------------------------------------------#

    k_layer_loop: DO k = 1, SIZE( pressure )


      ! --------------------------------
      ! Calculate multiplicative factors
      ! --------------------------------

      absorber_2( k ) = absorber( k ) * absorber( k )

      ! -- For the * terms
      d_absorber    = absorber( k )    - absorber( k-1 )
      d_absorber_TL = absorber_TL( k ) - absorber_TL( k-1 )

      ! -- For the ** terms
      factor_1    = ( absorber( k ) + absorber( k-1 ) ) * d_absorber
      factor_1_TL = ( ( absorber( k )    + absorber( k-1 )    ) * d_absorber_TL ) + &
                    ( ( absorber_TL( k ) + absorber_TL( k-1 ) ) * d_absorber    )

      ! -- For the *** terms       
      factor_2    = ( absorber_2( k ) + absorber_2( k-1 ) ) * d_absorber
      factor_2_TL = ( ( absorber_2( k ) + absorber_2( k-1 ) ) * d_absorber_TL ) + &
                    ( TWO * ( ( absorber( k )   * absorber_TL( k )   ) + &
                              ( absorber( k-1 ) * absorber_TL( k-1 ) ) ) * d_absorber )


      ! -------------------------------
      ! Calculate the intermediate sums
      ! -------------------------------

      ! -- T*
      s( 1 )    = s( 1 )    + ( temperature( k )    * d_absorber    )     ! Forward predictor
      s_TL( 1 ) = s_TL( 1 ) + ( temperature_TL( k ) * d_absorber    ) + &
                              ( temperature( k )    * d_absorber_TL )

      ! -- P*
      s( 2 )    = s( 2 )    + ( pressure( k )       * d_absorber    )     ! Forward predictor
      s_TL( 2 ) = s_TL( 2 ) + ( pressure_TL( k )    * d_absorber    ) + &
                              ( pressure( k )       * d_absorber_TL )

      ! -- T**
      s( 3 )    = s( 3 )    + ( temperature( k )    * factor_1    )       ! Forward predictor
      s_TL( 3 ) = s_TL( 3 ) + ( temperature_TL( k ) * factor_1    ) + &
                              ( temperature( k )    * factor_1_TL )

      ! -- P**
      s( 4 )    = s( 4 )    + ( pressure( k )       * factor_1    )       ! Forward predictor
      s_TL( 4 ) = s_TL( 4 ) + ( pressure_TL( k )    * factor_1    ) + &
                              ( pressure( k )       * factor_1_TL )

      ! -- T***
      s( 5 )    = s( 5 )    + ( temperature( k )    * factor_2    )       ! Forward predictor
      s_TL( 5 ) = s_TL( 5 ) + ( temperature_TL( k ) * factor_2    ) + &
                              ( temperature( k )    * factor_2_TL )

      ! -- P***
      s( 6 )    = s( 6 )    + ( pressure( k )       * factor_2    )       ! Forward predictor
      s_TL( 6 ) = s_TL( 6 ) + ( pressure_TL( k )    * factor_2    ) + &
                              ( pressure( k )       * factor_2_TL )


      ! ------------------------------------------------------
      ! Calculate the normalising factors for the integrated
      ! tangent-linear predictors. Note that the checks below,
      ! the IF tests to determine if the absorber products are
      ! represenatble are to minimise the number of calcs. I.e
      ! if inverse_1 is toast because absorber(k) is too small
      ! there's no need to check any further.
      ! ------------------------------------------------------

      ! -- Is inverse_1 representable?
      inverse_1_check: IF ( absorber( k ) > TOLERANCE ) THEN

        inverse_1 = ONE / absorber( k )

        ! -- Is inverse_2 representable
        inverse_2_check: IF ( absorber_2( k ) > TOLERANCE ) THEN

          inverse_2    =  inverse_1 * inverse_1
          inverse_1_TL = -inverse_2 * absorber_TL( k )
          absorber_3   =  absorber( k ) * absorber_2( k )
         
          ! -- Is inverse_3 representable
          inverse_3_check: IF ( absorber_3 > TOLERANCE ) THEN

            inverse_3    =  inverse_2 * inverse_1
            inverse_2_TL = -inverse_3 * absorber_TL( k ) * TWO
            absorber_4   =  absorber( k ) * absorber_3

            ! -- Is inverse_4 represenatble?
            inverse_4_check: IF ( absorber_4 > TOLERANCE ) THEN

              inverse_4    =  inverse_3 * inverse_1
              inverse_3_TL = -inverse_4 * absorber_TL( k ) * THREE

            ELSE

              inverse_3_TL = ZERO

            END IF inverse_4_check

          ELSE

            inverse_3 = ZERO

            inverse_2_TL = ZERO
            inverse_3_TL = ZERO

          END IF inverse_3_check

        ELSE

          inverse_2 = ZERO
          inverse_3 = ZERO

          inverse_1_TL = ZERO
          inverse_2_TL = ZERO
          inverse_3_TL = ZERO

        END IF inverse_2_check

      ELSE

        inverse_1 = ZERO
        inverse_2 = ZERO
        inverse_3 = ZERO

        inverse_1_TL = ZERO
        inverse_2_TL = ZERO
        inverse_3_TL = ZERO

      END IF inverse_1_check


      ! ------------------------------------------------------------
      ! Scale and normalise the tangent-linear integrated predictors
      ! ------------------------------------------------------------

      ! -- T*
      x_TL( 1, k ) = POINT_5  * ( ( s_TL( 1 ) * inverse_1    ) + &
                                  ( s( 1 )    * inverse_1_TL ) )

      ! -- P*
      x_TL( 2, k ) = POINT_5  * ( ( s_TL( 2 ) * inverse_1    ) + &
                                  ( s( 2 )    * inverse_1_TL ) )

      ! -- T**
      x_TL( 3, k ) = POINT_5  * ( ( s_TL( 3 ) * inverse_2    ) + &
                                  ( s( 3 )    * inverse_2_TL ) )

      ! -- P**
      x_TL( 4, k ) = POINT_5  * ( ( s_TL( 4 ) * inverse_2    ) + &
                                  ( s( 4 )    * inverse_2_TL ) )

      ! -- T***
      x_TL( 5, k ) = POINT_75 * ( ( s_TL( 5 ) * inverse_3    ) + &
                                  ( s( 5 )    * inverse_3_TL ) )

      ! -- P***
      x_TL( 6, k ) = POINT_75 * ( ( s_TL( 6 ) * inverse_3    ) + &
                                  ( s( 6 )    * inverse_3_TL ) )


      ! ----------------------------
      ! Sum predictors across layers
      ! ----------------------------

      DO i = 1, n_predictors

        predictor_TL( i, k ) = x_TL( i, k ) + x_TL( i, k - 1 )

      END DO

    END DO k_layer_loop

  END SUBROUTINE compute_int_predictors_TL




!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_predictors_AD
!
! PURPOSE:
!       PUBLIC routine to calculate the adjoint transmittance model
!       predictors.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_predictors_AD ( &
!                                    ! -- Forward input
!                                    pressure,       &  ! Input,  K
!                                    temperature,    &  ! Input,  K
!                                    water_vapor,    &  ! Input,  K
!                                    absorber,       &  ! Input,  0:K x J
!
!                                    ! -- Adjoint input
!                                    predictor_AD,   &  ! In/Output, I x K
!
!                                    ! -- Adjoint output
!                                    pressure_AD,    &  ! In/Output,  K
!                                    temperature_AD, &  ! In/Output,  K
!                                    water_vapor_AD, &  ! In/Output,  K
!                                    absorber_AD,    &  ! In/Output,  0:K x J
!
!                                    ! -- Optional input
!                                    no_standard     )  ! Optional input
!
! INPUT ARGUMENTS:
!       pressure:         Profile LAYER average pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       temperature:      Profile LAYER average temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       water_vapor:      Profile LAYER average water vapor mixing ratio array.
!                         UNITS:      g/kg
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       absorber:         Profile LEVEL integrated absorber amount array.
!                         UNITS:      Varies with absorber
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       predictor_AD:     Profile LAYER predictor adjoint array.
!                         ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  I x K
!                         ATTRIBUTES: INTENT( IN OUT ), TARGET
!
!
! OPTIONAL INPUT ARGUMENTS:
!       no_standard:      If present, the standard predictors are not calculated.
!                         This prevents recalculation of the standard predictors
!                         is only the view angle has changed - which only affects
!                         the integrated predictors.
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       pressure_AD:      Profile LAYER adjoint pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN OUT )
!
!       temperature_AD:   Profile LAYER adjoint temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN OUT )
!
!       water_vapor_AD:   Profile LAYER adjoint water vapor array.
!                         UNITS:      g/kg
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN OUT )
!
!       absorber_AD:      Profile LEVEL adjoint integrated absorber
!                         amount array.
!                         UNITS:      Varies with absorber
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  0:K
!                         ATTRIBUTES: INTENT( IN OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! CALLS:
!       compute_std_predictors_TL:   PRIVATE function to compute the tangent-
!                                    linear form of the standard (absorber
!                                    independent) predictor set.
!
!       compute_int_predictors_TL:   PRIVATE function to compute the tangent-
!                                    linear form of the absorber integrated
!                                    predictor set.
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
!       The predictors used in the NCEP transmittance model are organised in
!       the following manner:
!
!       --------------------------------------------
!       | 1 | 2 | 3 | ... | 9 | 10 | 11 | ... | 27 |
!       --------------------------------------------
!
!       \                    / \                   /
!        \                  /   \                 /
!         ------------------     -----------------
!                  |                      |
!                  v                      v
!
!              Standard               Integrated
!             Predictors            Predictors for
!                                   each absorber
!
!       Pointers are used to reference the module scope predictor data array
!       before the calls are made to the standard and integrated predictor
!       calculation functions. This eliminates array copying.
!S-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_predictors_AD ( &
                                     ! -- Forward input
                                     pressure,       &  ! Input,  K
                                     temperature,    &  ! Input,  K
                                     water_vapor,    &  ! Input,  K
                                     absorber,       &  ! Input,  0:K x J

                                     ! -- Adjoint input
                                     predictor_AD,   &  ! In/Output, I x K

                                     ! -- Adjoint output
                                     pressure_AD,    &  ! In/Output,  K
                                     temperature_AD, &  ! In/Output,  K
                                     water_vapor_AD, &  ! In/Output,  K
                                     absorber_AD,    &  ! In/Output,  0:K x J

                                     ! -- Optional input
                                     no_standard     )  ! Optional input



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward input
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )             :: pressure        ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )             :: temperature     ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )             :: water_vapor     ! K
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN )             :: absorber        ! 0:K x J

    ! -- Adjoint input
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN OUT ), TARGET :: predictor_AD    ! I x K

    ! -- Adjoint output
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN OUT )         :: pressure_AD     ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN OUT )         :: temperature_AD  ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN OUT )         :: water_vapor_AD  ! K
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN OUT )         :: absorber_AD     ! 0:K x J

    ! -- Optional input
    INTEGER,                             INTENT( IN ), OPTIONAL   :: no_standard


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_PREDICTORS_AD'


    ! ---------------
    ! Local variables
    ! ---------------

    CHARACTER( 80 ) :: message

    INTEGER :: i1, i2, j

    REAL( fp_kind ), POINTER, DIMENSION( :, : ) :: std_predictors_AD   ! Istd x K
    REAL( fp_kind ), POINTER, DIMENSION( :, : ) :: int_predictors_AD   ! Iint x K



    !#--------------------------------------------------------------------------#
    !#         -- Calculate the adjoint of the integrated predictors --         #
    !#--------------------------------------------------------------------------#

    j_absorber_loop: DO j = 1, SIZE( absorber, DIM = 2 )

      ! -- Determine indices of current absorber predictors
      i1 = MAX_N_STANDARD_PREDICTORS + ( ( j - 1 ) * MAX_N_INTEGRATED_PREDICTORS ) + 1
      i2 = i1 + MAX_N_INTEGRATED_PREDICTORS - 1

      ! -- Alias the input predictor array
      int_predictors_AD => predictor_AD( i1:i2, : )

      ! -- Compute the predictor adjoints for the current absorber
      CALL compute_int_predictors_AD( &
                                      ! -- Forward input
                                      pressure,            &  ! Input,  K
                                      temperature,         &  ! Input,  K
                                      absorber( 0:, j ),   &  ! Input,  0:K

                                      ! -- Adjoint input
                                      int_predictors_AD,   &  ! In/Output, Iint x K
                                                              
                                      ! -- Adjoint output
                                      pressure_AD,         &  ! In/Output,  K
                                      temperature_AD,      &  ! In/Output,  K
                                      absorber_AD( 0:, j ) )  ! In/Output,  0:K
                                                              
    END DO j_absorber_loop                                    




    !#--------------------------------------------------------------------------#
    !#     -- Calculate the adjoint of the standard predictors if needed --     #
    !#--------------------------------------------------------------------------#

    IF ( .NOT. PRESENT( no_standard ) ) THEN

      ! -- Alias the input predictor array
      std_predictors_AD => predictor_AD( 1:MAX_N_STANDARD_PREDICTORS, : )

      ! -- Compute the predictor adjoints
      CALL compute_std_predictors_AD( &
                                      ! -- Forward input
                                      pressure,          &  ! Input,  K
                                      temperature,       &  ! Input,  K
                                      water_vapor,       &  ! Input,  K

                                      ! -- Adjoint input
                                      std_predictors_AD, &  ! In/Output, Istd x K

                                      ! -- Adjoint output
                                      pressure_AD,       &  ! In/Output,  K
                                      temperature_AD,    &  ! In/Output,  K
                                      water_vapor_AD     )  ! In/Output,  K
    ELSE

      predictor_AD( 1:MAX_N_STANDARD_PREDICTORS, : ) = ZERO

    END IF

  END SUBROUTINE compute_predictors_AD





!--------------------------------------------------------------------------------
!P+
! NAME:
!       compute_std_predictors_AD
!
! PURPOSE:
!       PRIVATE function to compute the standard, absorber independent
!       predictor set for the adjoint transmittance model.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_std_predictors_AD( &
!                                       ! -- Forward input
!                                       pressure,       &  ! Input,  K
!                                       temperature,    &  ! Input,  K
!                                       water_vapor,    &  ! Input,  K
!
!                                       ! -- Adjoint input
!                                       predictor_AD,   &  ! In/Output, Istd x K
!
!                                       ! -- Adjoint output
!                                       pressure_AD,    &  ! In/Output,  K
!                                       temperature_AD, &  ! In/Output,  K
!                                       water_vapor_AD  )  ! In/Output,  K
!
! INPUT ARGUMENTS:
!       pressure:         Profile LAYER average pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       temperature:      Profile LAYER average temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       water_vapor:      Profile LAYER average water vapor mixing ratio array.
!                         UNITS:      g/kg
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       predictor_AD:     Adjoint of the LAYER predictor arrays.
!                         ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  Iint x K
!                         ATTRIBUTES: INTENT( IN OUT )
!
! OPTIONAL INPUT ARGUMENTS:
!       None
!
! OUTPUT ARGUMENTS:
!       pressure_AD:      Profile LAYER adjoint pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN OUT )
!
!       temperature_AD:   Profile LAYER adjoint temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN OUT )
!
!       water_vapor_AD:   Profile LAYER adjoint water vapor array.
!                         UNITS:      g/kg
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None
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
!       McMillin, L.M., L.J. Crone, M.D. Goldberg, and T.J. Kleespies,
!         "Atmospheric transmittance of an absorbing gas. 4. OPTRAN: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", Applied Optics, 1995, v34, pp6269-6274.
!
!       The standard predictors are the following:
!
!         pred(1) = T, Temperature
!         pred(2) = P, Pressure
!         pred(3) = T^2
!         pred(4) = P^2
!         pred(5) = T.P
!         pred(6) = T^2.P
!         pred(7) = T.P^2
!         pred(8) = T^2.P^2
!         pred(9) = W, Water vapor mixing ratio
!
!       The tangent-linear form of these are
!
!         dpred(1) = dT
!         dpred(2) = dP
!         dpred(3) = 2T.dT
!         dpred(4) = 2P.dP
!         dpred(5) = P.dT + T.dP
!         dpred(6) = 2TP.dT + T^2.dP
!         dpred(7) = 2TP.dP + P^2.dT
!         dpred(8) = 2T(P^2).dT + 2(T^2)P.dP
!         dpred(9) = dW
!
!       Thus, the adjoint of the predictors are
!
!         d#T = 2T( P^2.d#pred(8) + P.d#pred(6) + d#pred(3) ) + 
!                   P^2.d#pred(7) + P.d#pred(5) + d#pred(1)
!
!         d#P = 2P( T^2.d#pred(8) + T.d#pred(7) + d#pred(4) ) + 
!                   T^2.d#pred(6) + T.d#pred(5) + d#pred(2)
!
!         d#W = d#pred(9)
!
!       where the "#"ed terms represent the adjoints.
!
!P-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_std_predictors_AD( &
                                        ! -- Forward input
                                        p,            &  ! Input,  K
                                        t,            &  ! Input,  K
                                        w,            &  ! Input,  K

                                        ! -- Adjoint input
                                        predictor_AD, &  ! In/Output, Istd x K

                                        ! -- Adjoint output
                                        p_AD,         &  ! In/Output,  K
                                        t_AD,         &  ! In/Output,  K
                                        w_AD          )  ! In/Output,  K




    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward input
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )     :: p             ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )     :: t             ! Input,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )     :: w             ! Input,  K

    ! -- Adjoint input
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN OUT ) :: predictor_AD  ! In/Output, Istd x K

    ! -- Adjoint output
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN OUT ) :: p_AD          ! In/Output,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN OUT ) :: t_AD          ! In/Output,  K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN OUT ) :: w_AD          ! In/Output,  K



    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_STD_PREDICTORS_AD'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: i, n_predictors
    INTEGER :: k

    REAL( fp_kind ) :: p2, p2_AD
    REAL( fp_kind ) :: t2, t2_AD



    !#--------------------------------------------------------------------------#
    !#              -- Determine the number of predictors --                    #
    !#--------------------------------------------------------------------------#

    n_predictors = SIZE( predictor_AD, DIM = 1 )



    !#--------------------------------------------------------------------------#
    !#        -- Calculate the adjoints of the standard predictor set --        #
    !#                                                                          #
    !# Don't have to loop backwards here as this is a parallel loop.            #
    !#                                                                          #
    !# Pressure and temperature squared adjoint terms are not zeroed out every  #
    !# loop iteration as they are local to each iteration and can be simply     #
    !# re-assigned.                                                             #
    !#--------------------------------------------------------------------------#

    k_layer_loop: DO k = 1, SIZE( p )


      ! -- Precalculate the squared terms
      p2 = p( k ) * p( k )
      t2 = t( k ) * t( k )

      ! -- Pressure squared adjoint
      p2_AD = ( t2   * predictor_AD( 8, k ) ) + &   ! Predictor #8, T^2.P^2
              ( t(k) * predictor_AD( 7, k ) ) + &   ! Predictor #7, T.P^2
                       predictor_AD( 4, k )         ! Predictor #4, P^2

      ! -- Temperature squared adjoint
      t2_AD = ( p2   * predictor_AD( 8, k ) ) + &   ! Predictor #8, T^2.P^2
              ( p(k) * predictor_AD( 6, k ) ) + &   ! Predictor #6, T^2.P
                       predictor_AD( 3, k )         ! Predictor #3, T^2

      ! -- Water vapor adjoint
      w_AD( k ) = w_AD( k ) + predictor_AD( 9, k )  ! Predictor #9, W

      ! -- Temperature adjoint
      t_AD( k ) = t_AD( k ) + &
                  ( p2   * predictor_AD( 7, k ) ) + &   ! Predictor #7, T.P^2
                  ( p(k) * predictor_AD( 5, k ) ) + &   ! Predictor #5, T.P
                           predictor_AD( 1, k )   + &   ! Predictor #1, T
                  ( TWO * t(k) * t2_AD )                ! T^2 term

      ! -- Pressure adjoint
      p_AD( k ) = p_AD( k ) + &
                  ( t2   * predictor_AD( 6, k ) ) + &   ! Predictor #6, T^2.P
                  ( t(k) * predictor_AD( 5, k ) ) + &   ! Predictor #5, T.P
                           predictor_AD( 2, k )   + &   ! Predictor #2, P
                  ( TWO * p(k) * p2_AD )                ! P^2 term

      ! -- Zero output adjoint
      predictor_AD( :, k ) = ZERO

    END DO k_layer_loop

  END SUBROUTINE compute_std_predictors_AD





!--------------------------------------------------------------------------------
!P+
! NAME:
!       compute_int_predictors_AD
!
! PURPOSE:
!       PRIVATE function to compute the integrated, absorber dependent predictor
!       set for the adjoint transmittance model.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_int_predictors_AD( &
!                                       ! -- Forward input
!                                       pressure,       &  ! Input,  K
!                                       temperature,    &  ! Input,  K
!                                       absorber,       &  ! Input,  0:K
!
!                                       ! -- Adjoint input
!                                       predictor_AD,   &  ! In/Output, Iint x K
!
!                                       ! -- Adjoint output
!                                       pressure_AD,    &  ! In/Output,  K
!                                       temperature_AD, &  ! In/Output,  K
!                                       absorber_AD     )  ! In/Output,  0:K
!
! INPUT ARGUMENTS:
!       pressure:         Profile LAYER pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       temperature:      Profile LAYER temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN )
!
!       absorber :        Profile LEVEL integrated absorber amount array.
!                         UNITS:      Varies with absorber.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  0:K
!                         ATTRIBUTES: INTENT( IN )
!
!       predictor_AD:     Adjoint of the LAYER predictor arrays.
!                         ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                         UNITS:      Varies with predictor type.
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  Iint x K
!                         ATTRIBUTES: INTENT( IN OUT )
!
! OPTIONAL INPUT ARGUMENTS:
!       None
!
! OUTPUT ARGUMENTS:
!       pressure_AD:      Profile LAYER adjoint pressure array.
!                         UNITS:      hPa
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN OUT )
!
!       temperature_AD:   Profile LAYER adjoint temperature array.
!                         UNITS:      Kelvin
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  K
!                         ATTRIBUTES: INTENT( IN OUT )
!
!       absorber_AD:      Profile LEVEL adjoint integrated absorber
!                         amount array.
!                         UNITS:      Varies with absorber
!                         TYPE:       REAL( fp_kind )
!                         DIMENSION:  0:K
!                         ATTRIBUTES: INTENT( IN OUT )
!
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None
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
!       The input argument PREDICTOR_AD is set to zero on output.
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
!       The integrated predictors consist of six that are repeated for every
!       absorber:
!
!         1) T*
!         2) P*
!         3) T**
!         4) P**
!         5) T***
!         6) P***
!
!P-
!--------------------------------------------------------------------------------

  SUBROUTINE compute_int_predictors_AD( &
                                        ! -- Forward input
                                        pressure,       &  ! Input,  K
                                        temperature,    &  ! Input,  K
                                        absorber,       &  ! Input,  0:K

                                        ! -- Adjoint input
                                        predictor_AD,   &  ! In/Output, Iint x K

                                        ! -- Adjoint output
                                        pressure_AD,    &  ! In/Output,  K
                                        temperature_AD, &  ! In/Output,  K
                                        absorber_AD     )  ! In/Output,  0:K

 


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward input
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )     :: pressure        ! K
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )     :: temperature     ! K
    REAL( fp_kind ), DIMENSION( 0: ),   INTENT( IN )     :: absorber        ! 0:K

    ! -- Adjoint input
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( IN OUT ) :: predictor_AD    ! Iint x K

    ! -- Adjoint output
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN OUT ) :: pressure_AD     ! K
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN OUT ) :: temperature_AD  ! K
    REAL( fp_kind ), DIMENSION( 0: ),   INTENT( IN OUT ) :: absorber_AD     ! 0:K



    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_INT_PREDICTORS_AD'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Square of the absorber amount. 0:K
    REAL( fp_kind ), DIMENSION( 0:SIZE( pressure ) ) :: absorber_2

    ! -- Multiplicative factors, K
    REAL( fp_kind ), DIMENSION( SIZE( pressure ) ) :: d_absorber
    REAL( fp_kind ), DIMENSION( SIZE( pressure ) ) :: factor_1
    REAL( fp_kind ), DIMENSION( SIZE( pressure ) ) :: factor_2

    ! -- Intermediate summation arrays, Iint x 0:K and Iint
    REAL( fp_kind ), DIMENSION( SIZE( predictor_AD, DIM=1 ), 0:SIZE( pressure ) ) :: s
    REAL( fp_kind ), DIMENSION( SIZE( predictor_AD, DIM=1 ) ) :: s_AD

    ! -- LEVEL predictor, Iint x 0:K
    REAL( fp_kind ), DIMENSION( SIZE( predictor_AD, DIM=1 ), 0:SIZE( pressure ) ) :: x_AD

    ! -- Scalars
    INTEGER :: i, n_predictors
    INTEGER :: k, n_layers

    REAL( fp_kind ) :: d_absorber_AD
    REAL( fp_kind ) :: factor_1_AD
    REAL( fp_kind ) :: factor_2_AD

    REAL( fp_kind ) :: inverse_1
    REAL( fp_kind ) :: inverse_2
    REAL( fp_kind ) :: inverse_3
    REAL( fp_kind ) :: inverse_4
    REAL( fp_kind ) :: absorber_3
    REAL( fp_kind ) :: absorber_4
    REAL( fp_kind ) :: inverse_1_AD
    REAL( fp_kind ) :: inverse_2_AD
    REAL( fp_kind ) :: inverse_3_AD



    !#--------------------------------------------------------------------------#
    !#                       -- ASSIGN THE DIMENSIONS --                        #
    !#--------------------------------------------------------------------------#

    n_predictors = SIZE( predictor_AD, DIM=1 )
    n_layers     = SIZE( pressure )

!    IF ( n_predictors /= MAX_N_INTEGRATED_PREDICTORS ) THEN....



    !#--------------------------------------------------------------------------#
    !#          -- RECALCULATE THE INTERMEDIATE FORWARD MODEL SUMS --           #
    !#--------------------------------------------------------------------------#

    ! -----------------
    ! Initialise arrays
    ! -----------------

    absorber_2( 0 ) = ZERO
    s( :, 0: )      = ZERO


    ! ----------------
    ! Loop over layers
    ! ----------------

    k_layer_loop_forward: DO k = 1, n_layers


      ! -----------------------------------------
      ! Calculate absorber multiplicative factors
      ! and save for adjoint calculation.
      ! -----------------------------------------

      absorber_2( k ) = absorber( k ) * absorber( k )

      d_absorber( k ) = absorber( k ) - absorber( k-1 )                           ! For * terms
      factor_1( k )   = ( absorber( k )   + absorber( k-1 )   ) * d_absorber( k ) ! For ** terms
      factor_2( k )   = ( absorber_2( k ) + absorber_2( k-1 ) ) * d_absorber( k ) ! For *** terms


      ! ----------------------------------------
      ! Calculate and save the intermediate sums
      ! ----------------------------------------

      s( 1, k ) = s( 1, k-1 ) + ( temperature( k ) * d_absorber( k ) )  ! T*
      s( 2, k ) = s( 2, k-1 ) + ( pressure( k )    * d_absorber( k ) )  ! P*

      s( 3, k ) = s( 3, k-1 ) + ( temperature( k ) * factor_1( k ) )    ! T**
      s( 4, k ) = s( 4, k-1 ) + ( pressure( k )    * factor_1( k ) )    ! P**

      s( 5, k ) = s( 5, k-1 ) + ( temperature( k ) * factor_2( k ) )    ! T***
      s( 6, k ) = s( 6, k-1 ) + ( pressure( k )    * factor_2( k ) )    ! P***

    END DO k_layer_loop_forward



    !#--------------------------------------------------------------------------#
    !#                -- INITIALISE LOCAL ADJOINT VARIABLES --                  #
    !#--------------------------------------------------------------------------#

    x_AD( :, 0: ) = ZERO
    s_AD( : )     = ZERO

    absorber_AD( 0 ) = ZERO



    !#--------------------------------------------------------------------------#
    !#            -- CALCULATE THE INTEGRATED PREDICTOR ADJOINTS --             #
    !#--------------------------------------------------------------------------#


    ! -----------------------------------
    ! Here loop order does matter as this
    ! is a sequential loop
    ! -----------------------------------

    k_layer_loop_adjoint: DO k = n_layers, 1, -1


      ! -------------------------------------------------------
      ! Calculate the normalising factors for the integrated
      ! predictors. Note that the checks below, the IF tests to
      ! determine if the absorber products are represenatble
      ! are to minimise the number of calcs. I.e if inverse_1
      ! is toast because absorber(k) is too small there's no
      ! need to check any further.
      ! -------------------------------------------------------

      ! -- Is inverse_1 representable?
      inverse_1_check: IF ( absorber( k ) > TOLERANCE ) THEN

        inverse_1 = ONE / absorber( k )

        ! -- Is inverse_2 representable
        inverse_2_check: IF ( absorber_2( k ) > TOLERANCE ) THEN

          inverse_2  =  inverse_1 * inverse_1
          absorber_3 =  absorber( k ) * absorber_2( k )
         
          ! -- Is inverse_3 representable
          inverse_3_check: IF ( absorber_3 > TOLERANCE ) THEN

            inverse_3  =  inverse_2 * inverse_1
            absorber_4 =  absorber( k ) * absorber_3

            ! -- Is inverse_4 represenatble?
            inverse_4_check: IF ( absorber_4 > TOLERANCE ) THEN

              inverse_4 =  inverse_3 * inverse_1

            ELSE

              inverse_4 = ZERO

            END IF inverse_4_check

          ELSE

            inverse_3 = ZERO
            inverse_4 = ZERO

          END IF inverse_3_check

        ELSE

          inverse_2 = ZERO
          inverse_3 = ZERO
          inverse_4 = ZERO

        END IF inverse_2_check

      ELSE

        inverse_1 = ZERO
        inverse_2 = ZERO
        inverse_3 = ZERO
        inverse_4 = ZERO

      END IF inverse_1_check


      ! --------------------------------------------
      ! Adjoint of predictor summation across layers
      ! --------------------------------------------

      DO i = 1, n_predictors

        x_AD( i, k )   = x_AD( i, k )   + predictor_AD( i, k )
        x_AD( i, k-1 ) = x_AD( i, k-1 ) + predictor_AD( i, k )
        predictor_AD( i, k ) = ZERO

      END DO


      ! ----------------------------------------------------------------
      ! Adjoint of the scaled and normalised LEVEL integrated predictors
      !
      ! Note that the adjoint variables inverse_X_AD are local to this
      ! loop iteration so they are simply assigned when they are first
      ! used.
      ! ----------------------------------------------------------------

      ! -- P*** and T***, predictor indices #6 and 5
      s_AD( 6 )    = s_AD( 6 ) + ( POINT_75 * inverse_3 * x_AD( 6, k ) )
      s_AD( 5 )    = s_AD( 5 ) + ( POINT_75 * inverse_3 * x_AD( 5, k ) )
      inverse_3_AD = POINT_75 * ( ( s( 6, k ) * x_AD( 6, k ) ) + &
                                  ( s( 5, k ) * x_AD( 5, k ) ) )

      ! -- P** and T**, predictor indices #4 and 3
      s_AD( 4 )    = s_AD( 4 ) + ( POINT_5 * inverse_2 * x_AD( 4, k ) )
      s_AD( 3 )    = s_AD( 3 ) + ( POINT_5 * inverse_2 * x_AD( 3, k ) )
      inverse_2_AD = POINT_5 * ( ( s( 4, k ) * x_AD( 4, k ) ) + &
                                 ( s( 3, k ) * x_AD( 3, k ) ) )

      ! -- P* and T*, predictor indices #2 and 1
      ! -- Simply assign a value for inverse_1_AD
      s_AD( 2 )    = s_AD( 2 ) + ( POINT_5 * inverse_1 * x_AD( 2, k ) )
      s_AD( 1 )    = s_AD( 1 ) + ( POINT_5 * inverse_1 * x_AD( 1, k ) )
      inverse_1_AD = POINT_5 * ( ( s( 2, k ) * x_AD( 2, k ) ) + &
                                 ( s( 1, k ) * x_AD( 1, k ) ) )

      ! -- Zero out the LEVEL predictor adjoint array for current k.
      x_AD( :, k ) = ZERO

      ! -- Adjoint of inverse terms. Note that the inverse_X_AD
      ! -- terms are *not* zeroed out as they are re-assigned values
      ! -- each loop iteration above.
      absorber_AD( k ) = absorber_AD( k ) - (         inverse_2 * inverse_1_AD ) - &
                                            ( TWO *   inverse_3 * inverse_2_AD ) - &
                                            ( THREE * inverse_4 * inverse_3_AD )


      ! ----------------------------------------------------------
      ! Pressure and temperature adjoints of the intermediate sums
      ! ----------------------------------------------------------

      ! -- Pressure
      pressure_AD( k ) = pressure_AD( k ) + ( factor_2( k )   * s_AD( 6 ) ) + &  ! P***
                                            ( factor_1( k )   * s_AD( 4 ) ) + &  ! P**
                                            ( d_absorber( k ) * s_AD( 2 ) )      ! P*


      ! -- Temperature
      temperature_AD( k ) = temperature_AD( k ) + ( factor_2( k )   * s_AD( 5 ) ) + &  ! T***
                                                  ( factor_1( k )   * s_AD( 3 ) ) + &  ! T**
                                                  ( d_absorber( k ) * s_AD( 1 ) )      ! T*


      ! -------------------------------------
      ! Adjoint of the multiplicative factors
      ! -------------------------------------

      ! --------------------------------------------------
      ! Note that the adjoint variables factor_X_AD and
      ! d_absorber_AD are local to this loop iteration
      ! so they are simply assigned when they are first
      ! used.
      !
      ! Note there are no
      !   s_AD() = 0
      ! because all the tangent-linear forms are
      !   s_TL() = s_TL() + (...)
      ! summing from the previous layer.
      !
      ! Note that the factor_X_AD and d_absorber_AD
      ! terms are *not* zeroed out as they are re-assigned
      ! values each loop iteration.
      ! --------------------------------------------------

      ! -- Multiplicative factors
      factor_2_AD = ( pressure( k )    * s_AD( 6 ) ) + &
                    ( temperature( k ) * s_AD( 5 ) )

      factor_1_AD = ( pressure( k )    * s_AD( 4 ) ) + &
                    ( temperature( k ) * s_AD( 3 ) )

      d_absorber_AD = ( pressure( k )    * s_AD( 2 ) ) + &
                      ( temperature( k ) * s_AD( 1 ) )

      ! -- Adjoint of factor_2
      absorber_AD( k-1 ) = absorber_AD( k-1 ) + ( TWO * absorber( k-1 ) * d_absorber( k ) * factor_2_AD )
      absorber_AD(  k  ) = absorber_AD(  k  ) + ( TWO * absorber(  k  ) * d_absorber( k ) * factor_2_AD )
      d_absorber_AD      = d_absorber_AD      + ( ( absorber_2( k ) + absorber_2( k-1 ) ) * factor_2_AD )

      ! -- Adjoint of factor_1
      absorber_AD( k-1 ) = absorber_AD( k-1 ) + ( d_absorber( k ) * factor_1_AD )
      absorber_AD(  k  ) = absorber_AD(  k  ) + ( d_absorber( k ) * factor_1_AD )
      d_absorber_AD      = d_absorber_AD      + ( ( absorber( k ) + absorber( k-1 ) ) * factor_1_AD )

      ! -- Adjoint of d_absorber
      absorber_AD( k-1 ) = absorber_AD( k-1 ) - d_absorber_AD
      absorber_AD(  k  ) = absorber_AD(  k  ) + d_absorber_AD

    END DO k_layer_loop_adjoint

  END SUBROUTINE compute_int_predictors_AD

END MODULE predictors


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
! Revision 2.10  2001/08/16 16:45:05  paulv
! - Updated documentation
!
! Revision 2.9  2001/07/12 18:22:33  paulv
! - Added more robust code for calculating the powers of the inverse absorber
!   amount. The squares, cubes and fourth powers of absorber amounts were
!   causing floating point underflows which, when used in a denominator were
!   greatly inflating any precision errors. So, the solution I adopted was,
!   prior to each inverse calculation, to check the value of the absorber
!   quantity (e.g. absorber**2, absorber**3, or absorber**4). If they are less
!   than a tolerance value (defined using the EPSILON intrinsic) then their
!   inverse is set to zero - as are any higher power inverses. This prevents
!   too small values from being used.
!   This is mostly a problem for water vapor near at the top of the defined
!   atmosphere, although some low ozone values can be found (rarely).
! - Changed all loop initialisation to vector expressions, e.g. from
!     DO i = 1, n_predictors
!       s( i )       = ZERO
!       s_TL( i )    = ZERO
!       x_TL( i, 0 ) = ZERO
!     END DO
!   to
!     s( : )       = ZERO
!     s_TL( : )    = ZERO
!     x_TL( :, 0 ) = ZERO
! - Corrected bug in definition of the intermediate summation array in
!   COMPUTE_INT_PREDICTORS_AD from
!     REAL( fp_kind ), DIMENSION( SIZE( predictor_AD, DIM=1 ), SIZE( pressure ) ) :: s
!   to
!     REAL( fp_kind ), DIMENSION( SIZE( predictor_AD, DIM=1 ), 0:SIZE( pressure ) ) :: s
! - Added initialisation of absorber_2( 0 )
!     absorber_2( 0 ) = ZERO
!   in COMPUTE_INT_PREDICTORS_AD.
! - Corrected bug in intermediate summation loop. The sums were being calculated
!   like:
!     DO k = 1, n_layers
!       s( i, k ) = s( i, k ) + .......
!   rather than:
!     DO k = 1, n_layers
!       s( i, k ) = s( i, k-1 ) + .......
!   hence the need for the redefintion of s(K) to s(0:K).
!
! Revision 2.8  2001/05/29 17:51:34  paulv
! - Corrected some documentation errors.
!
! Revision 2.7  2001/05/04 14:39:43  paulv
! - Removed shared predictor arrays from module.
! - Now use TYPE_KINDS module parameter FP_KIND to set the floating point
!   data type.
! - Added adjoint form of routines to module. Use of NO_STANDARD optional
!   keyword in COMPUTE_PREDICTORS_AD has not yet been looked into.
! - Changed names of PRIVATE integrated predictor routines from
!   COMPUTE_STARD_PREDICTORS to COMPUTE_INT_PREDICTORS. The difference between
!   "stard" and "std" is small enough to cause mild confusion when scanning
!   the code.
! - Changed references to parameter maximums to input array sizes, e.g.
!     j_absorber_loop: DO j = 1, MAX_N_ABSORBERS
!   becomes
!     j_absorber_loop: DO j = 1, SIZE( absorber, DIM = 2 )
! - Shortened variable names in standard predictor routines, e.g. p instead
!   of pressure. These calcs are simple enough that short names are clear.
! - ABSORBER arrays are now dimensioned as 0:K. This eliminates the
!   need for using an ABSORBER_KM1 variable in computing the absorber layer
!   difference, D_ABSORBER, and average, AVE_ABSORBER where the layer loop
!   always goes from 1 -> n_layers.
! - Simplified COMPUTE_INT_PREDICTORS_TL. This is a combination of the
!   change in the absorber array dimensioning and just thinking about it
!   for a while.
! - Updated header documentation. Adjoint routines not yet fully documented.
!
! Revision 2.6  2001/04/03 20:02:56  paulv
! - Commented out shared predictor data arrays. Predictor arrays are now
!   passed arguments.
! - Removed reference to profile number. Calls to routines are now a simgle
!   profile passed per call.
! - Removed planned allocation of predictor arrays.
! - Correted bug in 1st and 2nd order predictor calculation.
!
! Revision 2.5  2001/01/24 20:14:21  paulv
! - Latest test versions.
!
! Revision 2.4  2000/11/09 20:36:07  paulv
! - Added tangent linear form of routines in this module.
! - Changed names of PRIVATE routines used to compute predictors. The
!   appearance of "standard" and "integrated" in routine names have been
!   replaced with "std" and "stard" respectively.
! - Initial setup in place for dynamic allocation of predictor arrays. These
!   arrays are still statically allocated with parameterised dimensions.
!
! Revision 2.3  2000/08/31 19:36:33  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 2.2  2000/08/24 16:43:29  paulv
! - Added optional NO_STANDARD argument to COMPUTE_PREDICTORS subprogram.
!   Using this argument prevents the standard predictors (which are angle
!   independent) from being recalculated when only the path angle has changed
!   in the calling procedure.
! - Updated module and subprogram documentation.
!
! Revision 2.1  2000/08/21 21:03:16  paulv
! - Standard and integrated predictor sets calculated in separate
!   functions.
! - Predictor values saved in linear store. Simplifies application of
!   predictors when calculating absorption coefficients.
! - Wrapper function "compute_predictors" added to simplify predictor
!   calculation.
!
! Revision 1.1  2000/08/08 16:34:17  paulv
! Initial checkin
!
!
!
