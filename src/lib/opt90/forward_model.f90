!------------------------------------------------------------------------------
!M+
! NAME:
!       forward_model
!
! PURPOSE:
!       Module containing the NCEP RT forward model function.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE forward_model
!
! OUTPUTS:
!       None.
!
! MODULES:
!       type_kinds:            Module to define kind types for variable declaration.
!
!       error_handler:         Module to define error codes and handle error conditions
!
!       parameters:            Module containing parameter definitions for the
!                              RT model.
!
!       spectral_coefficients: Module containing the RT model spectral coefficients.
!
!       absorber_profile:      Module containing routines for generating the absorber
!                              profiles.
!
!       predictors:            Module containing routines for generating the predictor
!                              profiles.
!
!       transmittance:         Module containing transmittance calculation routines.
!
!       radiance:              Module containing radiance calculation routines.
!
! CONTAINS:
!       compute_rtm:           PUBLIC function that calculates the forward model 
!                              top-of-atmosphere (TOA) radiances and brightness 
!                              temperatures for an input atmospheric profile set and
!                              user specified satellites/channels.
!
!                              This function is simply a wrapper around the FORWARD
!                              model so that the user doesn't have to declare the
!                              absorber/predictor/etc arrays in the calling routine.
!
!       forward_rtm:           PUBLIC function that calculates top-of-atmosphere (TOA)
!                              radiances and brightness temperatures for user specified
!                              profiles and satellites/channels.
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
!       Written by:     Paul van Delst, CIMSS@NOAA/NCEP 15-July-2000
!                       pvandelst@ncep.noaa.gov
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

MODULE forward_model


  ! ------------
  ! Module usage
  ! ------------

  USE type_kinds,            ONLY : fp_kind
  USE error_handler
  USE parameters
  USE spectral_coefficients, ONLY : is_solar_channel, &
                                    is_microwave_channel
  USE absorber_profile,      ONLY : compute_absorber_amount, &
                                    find_absorber_layer_index
  USE predictors,            ONLY : compute_predictors
  USE transmittance,         ONLY : compute_transmittance
  USE radiance,              ONLY : compute_radiance


  ! -----------------------
  ! Disable implicit typing
  ! -----------------------

  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------

  PRIVATE
  PUBLIC :: compute_rtm
  PUBLIC :: forward_rtm


  ! --------------------
  ! Function overloading
  ! --------------------

  INTERFACE compute_rtm
    MODULE PROCEDURE compute_rtm_rank1
    MODULE PROCEDURE compute_rtm_rank2
  END INTERFACE ! compute_rtm

  INTERFACE forward_rtm
    MODULE PROCEDURE forward_rtm_rank1
    MODULE PROCEDURE forward_rtm_rank2
  END INTERFACE ! forward_rtm

CONTAINS


!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_rtm
!
! PURPOSE:
!       PUBLIC function that calculates the forward model top-of-atmosphere (TOA)
!       radiances and brightness temperatures for an input atmospheric profile
!       set and user specified satellites/channels.
!
!       This function is simply a wrapper around the FORWARD model so that the
!       user doesn't have to declare the absorber/predictor/etc arrays in the
!       calling routine.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       result = compute_rtm( &
!                             ! -- Forward inputs
!                             level_p, layer_p, layer_t, layer_w, layer_o, &  ! Input, K x M
!
!                             surface_temperature,                         &  ! Input, M
!                             surface_emissivity,                          &  ! Input, L*M
!                             surface_reflectivity,                        &  ! Input, L*M
!
!                             ! -- Other inputs
!                             secant_view_angle,                           &  ! Input, M
!                             secant_solar_angle,                          &  ! Input, M
!                             n_channels_per_profile,                      &  ! Input, M
!                             channel_index,                               &  ! Input, L*M
!
!                             ! -- Forward output
!                             tau,                                         &  ! Input, K x L*M
!                             flux_tau,                                    &  ! Input, K x L*M
!                             solar_tau,                                   &  ! Input, K x L*M
!
!                             upwelling_radiance,                          &  ! Input, L*M
!                             brightness_temperature,                      &  ! Input, L*M
!
!                             ! Optional inputs
!                             n_input_profiles = n_input_profiles,         &
!                             message_log      = message_log               )
!
! INPUT ARGUMENTS:
!
!       level_p:                   Profile set layer interface pressure array. The TOA
!                                  pressure is not included. TOA pressure is parameterised
!                                  in the PARAMETERS module.
!                                  UNITS:      hPa
!                                  TYPE:       Real
!                                  DIMENSION:  K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_p:                   Profile set layer average pressure array.
!                                  UNITS:      hPa
!                                  TYPE:       Real
!                                  DIMENSION:  K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_t:                   Profile set layer average temperature array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_w:      .            Profile set layer average water vapor mixing ratio array
!                                  UNITS:      g/kg
!                                  TYPE:       Real
!                                  DIMENSION:  K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_o:                   Profile set layer average ozone mixing ratio array.
!                                  UNITS:      ppmv
!                                  TYPE:       Real
!                                  DIMENSION:  K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_temperature:       Profile set surface temperature array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  M; M > or = 1 (i.e. scalar)
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_emissivity:        Profile set surface emissivity array
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_reflectivity:      Profile set surface reflectivity array
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!       secant_view_angle:         Secant of the satellite view angle measured
!                                  from nadir for each profile in the set.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  M; M > or = 1 (i.e. scalar)
!                                  ATTRIBUTES: INTENT( IN )
!
!       secant_solar_angle:        Secant of the solar zenith angle for each
!                                  profile in the set.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  M; M > or = 1 (i.e. scalar)
!                                  ATTRIBUTES: INTENT( IN )
!
!       n_channels_per_profile:    The number of channels for each profile in the
!                                  set for which radiances are required.
!                                  UNITS:      None
!                                  TYPE:       Integer
!                                  DIMENSION:  M; M > or = 1 (i.e. scalar)
!                                  ATTRIBUTES: INTENT( IN )
!
!       channel_index:             Channel index id array. Each element is a unique
!                                  index to a (supported) sensor channel.
!                                  UNITS:      None
!                                  TYPE:       Integer
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!
! OPTIONAL INPUT ARGUMENTS:
!
!       n_input_profiles:          The number of profiles in the passed arrays to process.
!                                  If not specified, the default value is the SECOND dimension
!                                  of the pressure array determined using the SIZE intrinsic,
!                                  N_PROFILES. If N_INPUT_PROFILES is specified and is < 1 or
!                                  greater than N_PROFILES, the default value is set to N_PROFILES.
!                                  This argument is ignored if the input profile arrays are
!                                  vectors, i.e. a single profile.
!                                  UNITS:      None
!                                  TYPE:       Integer
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN ), OPTIONAL
!
!       message_log:               Character string specifying a filename in which any
!                                  messages will be logged. If not specified, or if an
!                                  error occurs opening the log file, the default action
!                                  is to output messages to the screen.
!                                  UNITS:      None
!                                  TYPE:       Character
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! OUTPUT ARGUMENTS:
!
!       tau:                       Layer->TOA transmittance for the satellite
!                                  view angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( OUT )
!
!       flux_tau:                  Layer->SFC transmittance for the default
!                                  diffusivity angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( OUT )
!
!       solar_tau:                 Layer->SFC transmittance for the solar
!                                  zenith angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( OUT )
!
!       upwelling_radiance:        TOA radiances for each channel/profile.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( OUT )
!
!       brightness_temperature:    TOA brightness temperatures corresponding
!                                  to the TOA radiances.
!                                  N.B.: Set to ZERO upon output.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result = SUCCESS => Calculation was successful
!              = FAILURE => Error occurred
!
! CALLS:
!      display_message:            Subroutine to output messages
!                                  SOURCE: error_handler module
!
!      get_max_n_channels:         Routine to retrieve the value of the
!                                  MAX_N_CHANNELS "pseudo-parameter".
!                                  SOURCE: parameters module
!
!      forward_rtm:                Function to construct the forward model and calculate
!                                  the transmittance profiles and TOA radiance/temperatures.
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       All input adjoint quantities are set to ZERO on output.
!
! RESTRICTIONS:
!       Note the restrictions on the input array dimensions:
!         K == n_layers > 1
!         L == n_channels > 1
!         M == n_profiles > or = 1
!
! PROCEDURE:
!       See individual module function documentation.
!S-
!--------------------------------------------------------------------------------

  FUNCTION compute_rtm_rank2( &
             ! -- Forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o, &  ! Input, K x M

             surface_temperature,                         &  ! Input, M
             surface_emissivity,                          &  ! Input, L*M
             surface_reflectivity,                        &  ! Input, L*M

             ! -- Other inputs
             secant_view_angle,                           &  ! Input, M
             secant_solar_angle,                          &  ! Input, M
             n_channels_per_profile,                      &  ! Input, M
             channel_index,                               &  ! Input, L*M

             ! -- Forward output
             tau,                                         &  ! Output, K x L*M
             flux_tau,                                    &  ! Output, K x L*M
             solar_tau,                                   &  ! Output, K x L*M

             upwelling_radiance,                          &  ! Output, L*M
             brightness_temperature,                      &  ! Output, L*M

             ! -- Optional inputs
             n_input_profiles,                            &
             message_log )                                &

           RESULT ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- TYPE DECLARATIONS --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward inputs
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( IN )  :: level_p                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( IN )  :: layer_p                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( IN )  :: layer_t                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( IN )  :: layer_w                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( IN )  :: layer_o                 ! K x M

    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: surface_temperature     ! M
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: surface_emissivity      ! L*M
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: surface_reflectivity    ! L*M

    ! -- Other inputs
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: secant_view_angle       ! M
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: secant_solar_angle      ! M
    INTEGER,         DIMENSION( : ),    INTENT( IN )  :: n_channels_per_profile  ! M
    INTEGER,         DIMENSION( : ),    INTENT( IN )  :: channel_index           ! L*M

    ! -- Forward outputs
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( OUT ) :: tau                     ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( OUT ) :: flux_tau                ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( OUT ) :: solar_tau               ! K x L*M

    REAL( fp_kind ), DIMENSION( : ),    INTENT( OUT ) :: upwelling_radiance      ! L*M
    REAL( fp_kind ), DIMENSION( : ),    INTENT( OUT ) :: brightness_temperature  ! L*M

    ! -- Optional input
    INTEGER,        OPTIONAL,           INTENT( IN )  :: n_input_profiles        ! Scalar
    CHARACTER( * ), OPTIONAL,           INTENT( IN )  :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_RTM_RANK2'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Array for integrated absorber amounts, 0:K x J x M
    REAL( fp_kind ), DIMENSION( 0:SIZE( layer_p, DIM = 1 ), &
                                  MAX_N_ABSORBERS,          &
                                  SIZE( layer_p, DIM = 2 )  ) :: absorber

    ! -- Arrays for absorber space indexing, K x J x M
    INTEGER,         DIMENSION( SIZE( layer_p, DIM = 1 ), &
                                MAX_N_ABSORBERS,          &
                                SIZE( layer_p, DIM = 2 )  ) :: tau_layer_index,      &
                                                               flux_tau_layer_index, &
                                                               solar_tau_layer_index

    ! -- Arrays for predictors, Imax x K x M
    REAL( fp_kind ), DIMENSION( MAX_N_PREDICTORS,         &
                                SIZE( layer_p, DIM = 1 ), &
                                SIZE( layer_p, DIM = 2 )  ) :: tau_predictor,      &
                                                               flux_tau_predictor, &
                                                               solar_tau_predictor

    ! -- Array for layer Planck radiance term, K x L*M
    REAL( fp_kind ), DIMENSION( SIZE( layer_p, DIM = 1 ),  &
                                SIZE( upwelling_radiance ) ) :: layer_radiance
      

    ! -- Array for downwelling radiance (flux + solar), L*M
    REAL( fp_kind ), DIMENSION( SIZE( upwelling_radiance ) ) :: downwelling_radiance



    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC ADJUSTL, &
              MAXVAL,  &
              SIZE,    &
              TRIM



    !#--------------------------------------------------------------------------#
    !#         -- Compute the forward radiances and temperatures --             #
    !#--------------------------------------------------------------------------#

    error_status = forward_rtm_rank2( &
                     ! -- Forward inputs
                     level_p, layer_p, layer_t, layer_w, layer_o,  &  ! Input,  K x M

                     surface_temperature,                          &  ! Input,  M
                     surface_emissivity,                           &  ! Input,  L*M
                     surface_reflectivity,                         &  ! Input,  L*M

                     ! -- Other inputs
                     secant_view_angle,                            &  ! Input,  M
                     secant_solar_angle,                           &  ! Input,  M
                     n_channels_per_profile,                       &  ! Input,  M
                     channel_index,                                &  ! Input,  L*M

                     ! -- Outputs
                     absorber,                                     &  ! Output, 0:K x J x M

                     tau_layer_index,                              &  ! Output, K x J x M
                     flux_tau_layer_index,                         &  ! Output, K x J x M
                     solar_tau_layer_index,                        &  ! Output, K x J x M

                     tau_predictor,                                &  ! Output, Imax x K x M
                     flux_tau_predictor,                           &  ! Output, Imax x K x M
                     solar_tau_predictor,                          &  ! Output, Imax x K x M

                     tau,                                          &  ! Output, K x L*M
                     flux_tau,                                     &  ! Output, K x L*M
                     solar_tau,                                    &  ! Output, K x L*M

                     layer_radiance,                               &  ! Output, K x L*M
                     downwelling_radiance,                         &  ! Output, L*M
                     upwelling_radiance,                           &  ! Output, L*M

                     brightness_temperature,                       &  ! Output, L*M

                     ! -- Optional inputs
                     n_input_profiles = n_input_profiles,          &
                     message_log      = message_log                )


    ! -------------------------------
    ! Check for successful completion
    ! -------------------------------

    IF ( error_status /= SUCCESS ) THEN

      CALL display_message( ROUTINE_NAME, &
                            'Error occured in FORWARD_RTM_RANK2', &
                            error_status, &
                            message_log = message_log )
      RETURN

    END IF



    !#--------------------------------------------------------------------------#
    !#                              -- DONE --                                  #
    !#--------------------------------------------------------------------------#

    error_status = SUCCESS


  END FUNCTION compute_rtm_rank2



  FUNCTION compute_rtm_rank1( &
             ! -- Forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o, &  ! Input, K

             surface_temperature,                         &  ! Input, Scalar
             surface_emissivity,                          &  ! Input, L
             surface_reflectivity,                        &  ! Input, L

             ! -- Other inputs
             secant_view_angle,                           &  ! Input, Scalar
             secant_solar_angle,                          &  ! Input, Scalar
             n_channels,                                  &  ! Input, Scalar
             channel_index,                               &  ! Input, L

             ! -- Forward output
             tau,                                         &  ! Output, K x L
             flux_tau,                                    &  ! Output, K x L
             solar_tau,                                   &  ! Output, K x L

             upwelling_radiance,                          &  ! Output, L
             brightness_temperature,                      &  ! Output, L

             ! -- Optional inputs
             n_input_profiles,                            &
             message_log )                                &

           RESULT ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- TYPE DECLARATIONS --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward inputs
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: level_p                 ! K
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: layer_p                 ! K
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: layer_t                 ! K
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: layer_w                 ! K
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: layer_o                 ! K

    REAL( fp_kind ),                    INTENT( IN )  :: surface_temperature     ! Scalar
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: surface_emissivity      ! L
    REAL( fp_kind ), DIMENSION( : ),    INTENT( IN )  :: surface_reflectivity    ! L

    ! -- Other inputs
    REAL( fp_kind ),                    INTENT( IN )  :: secant_view_angle       ! Scalar
    REAL( fp_kind ),                    INTENT( IN )  :: secant_solar_angle      ! Scalar
    INTEGER,                            INTENT( IN )  :: n_channels              ! Scalar
    INTEGER,         DIMENSION( : ),    INTENT( IN )  :: channel_index           ! L

    ! -- Forward outputs
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( OUT ) :: tau                     ! K x L
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( OUT ) :: flux_tau                ! K x L
    REAL( fp_kind ), DIMENSION( :, : ), INTENT( OUT ) :: solar_tau               ! K x L

    REAL( fp_kind ), DIMENSION( : ),    INTENT( OUT ) :: upwelling_radiance      ! L
    REAL( fp_kind ), DIMENSION( : ),    INTENT( OUT ) :: brightness_temperature  ! L

    ! -- Optional input. Note that N_INPUT_PROFILES is not used in this
    ! -- function. It is included here so if a user specifies it by mistake
    ! -- for rank-1 profile input the code won't (hopefully) fall in a heap.
    INTEGER,        OPTIONAL,           INTENT( IN )  :: n_input_profiles        ! Scalar
    CHARACTER( * ), OPTIONAL,           INTENT( IN )  :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_RTM_RANK1'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Array for integrated absorber amounts, 0:K x J
    REAL( fp_kind ), DIMENSION( 0:SIZE( layer_p, DIM = 1 ), &
                                  MAX_N_ABSORBERS           ) :: absorber

    ! -- Arrays for absorber space indexing, K x J
    INTEGER,         DIMENSION( SIZE( layer_p, DIM = 1 ), &
                                MAX_N_ABSORBERS           ) :: tau_layer_index,      &
                                                               flux_tau_layer_index, &
                                                               solar_tau_layer_index

    ! -- Arrays for predictors, Imax x K
    REAL( fp_kind ), DIMENSION( MAX_N_PREDICTORS,         &
                                SIZE( layer_p, DIM = 1 )  ) :: tau_predictor,      &
                                                               flux_tau_predictor, &
                                                               solar_tau_predictor

    ! -- Array for layer Planck radiance term, K x L
    REAL( fp_kind ), DIMENSION( SIZE( layer_p, DIM = 1 ),  &
                                SIZE( upwelling_radiance ) ) :: layer_radiance
      

    ! -- Array for downwelling radiance (flux + solar), L
    REAL( fp_kind ), DIMENSION( SIZE( upwelling_radiance ) ) :: downwelling_radiance



    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC ADJUSTL, &
              MAXVAL,  &
              SIZE,    &
              TRIM



    !#--------------------------------------------------------------------------#
    !#         -- COMPUTE THE FORWARD RADIANCES AND TEMPERATURES --             #
    !#--------------------------------------------------------------------------#

    error_status = forward_rtm_rank1( &
                     ! -- Forward inputs
                     level_p, layer_p, layer_t, layer_w, layer_o,  &  ! Input,  K

                     surface_temperature,                          &  ! Input,  Scalar
                     surface_emissivity,                           &  ! Input,  L
                     surface_reflectivity,                         &  ! Input,  L

                     ! -- Other inputs
                     secant_view_angle,                            &  ! Input,  Scalar
                     secant_solar_angle,                           &  ! Input,  Scalar
                     n_channels,                                   &  ! Input,  Scalar
                     channel_index,                                &  ! Input,  L

                     ! -- Outputs
                     absorber,                                     &  ! Output, 0:K x J

                     tau_layer_index,                              &  ! Output, K x J
                     flux_tau_layer_index,                         &  ! Output, K x J
                     solar_tau_layer_index,                        &  ! Output, K x J

                     tau_predictor,                                &  ! Output, Imax x K
                     flux_tau_predictor,                           &  ! Output, Imax x K
                     solar_tau_predictor,                          &  ! Output, Imax x K

                     tau,                                          &  ! Output, K x L
                     flux_tau,                                     &  ! Output, K x L
                     solar_tau,                                    &  ! Output, K x L

                     layer_radiance,                               &  ! Output, K x L
                     downwelling_radiance,                         &  ! Output, L
                     upwelling_radiance,                           &  ! Output, L

                     brightness_temperature,                       &  ! Output, L

                     message_log = message_log )


    ! -------------------------------
    ! Check for successful completion
    ! -------------------------------

    IF ( error_status /= SUCCESS ) THEN

      CALL display_message( ROUTINE_NAME, &
                            'Error occured in FORWARD_RTM_RANK1', &
                            error_status, &
                            message_log = message_log )
      RETURN

    END IF



    !#--------------------------------------------------------------------------#
    !#                              -- DONE --                                  #
    !#--------------------------------------------------------------------------#

    error_status = SUCCESS


  END FUNCTION compute_rtm_rank1







!--------------------------------------------------------------------------------
!S+
! NAME:
!       forward_rtm
!
! PURPOSE:
!       PUBLIC function that calculates top-of-atmosphere (TOA) radiances
!       and brightness temperatures for an input atmospheric profile or
!       profile set and user specified satellites/channels.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!
!       result = forward_rtm( &
!                             ! -- Inputs
!                             level_p, layer_p, layer_t, layer_w, layer_o, &  ! Input, K x M
!
!                             surface_temperature,                         &  ! Input,  M
!                             surface_emissivity,                          &  ! Input,  L*M
!                             surface_reflectivity,                        &  ! Input,  L*M
!
!                             secant_view_angle,                           &  ! Input,  M
!                             secant_solar_angle,                          &  ! Input,  M
!                             n_channels_per_profile,                      &  ! Input,  M
!                             channel_index,                               &  ! Input,  L*M
!
!                             ! -- Outputs
!                             absorber,                                    &  ! Output, 0:K x J x M
!
!                             tau_layer_index,                             &  ! Output, K x J x M
!                             flux_tau_layer_index,                        &  ! Output, K x J x M
!                             solar_tau_layer_index,                       &  ! Output, K x J x M
!
!                             tau_predictor,                               &  ! Output, Imax x K x M
!                             flux_tau_predictor,                          &  ! Output, Imax x K x M
!                             solar_tau_predictor,                         &  ! Output, Imax x K x M
!
!                             tau,                                         &  ! Output, K x L*M
!                             flux_tau,                                    &  ! Output, K x L*M
!                             solar_tau,                                   &  ! Output, K x L*M
!
!                             layer_radiance,                              &  ! Output, K x L*M
!                             downwelling_radiance,                        &  ! Output, L*M
!                             upwelling_radiance,                          &  ! Output, L*M
!
!                             brightness_temperature,                      &  ! Output, L*M
!
!                             ! -- Optional inputs
!                             n_input_profiles = n_input_profiles,         &
!                             message_log      = message_log               )
!                             
! INPUT ARGUMENTS:
!
!       level_p:                 Profile set layer interface pressure array. The TOA
!                                pressure is not included. TOA pressure is parameterised
!                                in the PARAMETERS module.
!                                UNITS:      hPa
!                                TYPE:       Real
!                                DIMENSION:  K x M; K > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( IN )
!
!       layer_p:                 Profile set layer average pressure array.
!                                UNITS:      hPa
!                                TYPE:       Real
!                                DIMENSION:  K x M; K > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( IN )
!
!       layer_t:                 Profile set layer average temperature array.
!                                UNITS:      Kelvin
!                                TYPE:       Real
!                                DIMENSION:  K x M; K > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( IN )
!
!       layer_w:      .          Profile set layer average water vapor mixing ratio array
!                                UNITS:      g/kg
!                                TYPE:       Real
!                                DIMENSION:  K x M; K > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( IN )
!
!       layer_o:                 Profile set layer average ozone mixing ratio array.
!                                UNITS:      ppmv
!                                TYPE:       Real
!                                DIMENSION:  K x M; K > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( IN )
!
!       surface_temperature:     Profile set surface temperature array.
!                                UNITS:      Kelvin
!                                TYPE:       Real
!                                DIMENSION:  M; M > or = 1 (i.e. scalar)
!                                ATTRIBUTES: INTENT( IN )
!
!       surface_emissivity:      Profile set surface emissivity array
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  L*M; L > 1, and M > or = 1
!                                            NB: This is a 1-D array.
!                                ATTRIBUTES: INTENT( IN )
!
!       surface_reflectivity:    Profile set surface reflectivity array
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  L*M; L > 1, and M > or = 1
!                                            NB: This is a 1-D array.
!                                ATTRIBUTES: INTENT( IN )
!
!       secant_view_angle:       Secant of the satellite view angle measured
!                                from nadir for each profile in the set.
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  M; M > or = 1 (i.e. scalar)
!                                ATTRIBUTES: INTENT( IN )
!
!       secant_solar_angle:      Secant of the solar zenith angle for each
!                                profile in the set.
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  M; M > or = 1 (i.e. scalar)
!                                ATTRIBUTES: INTENT( IN )
!
!       n_channels_per_profile:  The number of channels for each profile in the
!                                set for which radiances are required.
!                                UNITS:      None
!                                TYPE:       Integer
!                                DIMENSION:  M; M > or = 1 (i.e. scalar)
!                                ATTRIBUTES: INTENT( IN )
!
!       channel_index:           Channel index id array. Each element is a unique
!                                index to a (supported) sensor channel.
!                                UNITS:      None
!                                TYPE:       Integer
!                                DIMENSION:  L*M; L > 1, and M > or = 1
!                                            NB: This is a 1-D array.
!                                ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!
!       n_input_profiles:        The number of profiles in the passed arrays to process.
!                                If not specified, the default value is the SECOND dimension
!                                of the pressure array determined using the SIZE intrinsic,
!                                N_PROFILES. If N_INPUT_PROFILES is specified and is < 1 or
!                                greater than N_PROFILES, the default value is set to N_PROFILES.
!                                This argument is ignored if the input profile arrays are
!                                vectors, i.e. a single profile.
!                                UNITS:      None
!                                TYPE:       Integer
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( IN ), OPTIONAL
!
!       message_log:             Character string specifying a filename in which any
!                                messages will be logged. If not specified, or if an
!                                error occurs opening the log file, the default action
!                                is to output messages to the screen.
!                                UNITS:      None
!                                TYPE:       Character
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! OUTPUT ARGUMENTS:
!
!       absorber:                Array of absorber amount for nadir view.
!                                UNITS:      Absorber dependent.
!                                TYPE:       Real
!                                DIMENSION:  0:K x J x M
!                                ATTRIBUTES: INTENT( OUT )
!
!       tau_layer_index:         Array of absorber space layer indices of the input
!                                absorber amounts at the satellite view angle.
!                                UNITS:      None.
!                                TYPE:       Integer
!                                DIMENSION:  K x J x M
!                                ATTRIBUTES: INTENT( OUT )
!
!       flux_tau_layer_index:    Array of absorber space layer indices of the input
!                                absorber amounts at the default diffusivity angle.
!                                UNITS:      None.
!                                TYPE:       Integer
!                                DIMENSION:  K x J x M
!                                ATTRIBUTES: INTENT( OUT )
!
!       solar_tau_layer_index:   Array of absorber space layer indices of the input
!                                absorber amounts at the solar zenith angle.
!                                UNITS:      None.
!                                TYPE:       Integer
!                                DIMENSION:  K x J x M
!                                ATTRIBUTES: INTENT( OUT )
!
!       tau_predictor:           Predictor profiles for the layer->TOA transmittance.
!                                UNITS:      None.
!                                TYPE:       Real
!                                DIMENSION:  I x K x M; K > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( OUT )
!
!       flux_tau_predictor:      Predictor profiles for the thermal flux transmittance.
!                                UNITS:      None.
!                                TYPE:       Real
!                                DIMENSION:  I x K x M; K > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( OUT )
!
!       solar_tau_predictor:     Predictor profiles for the solar transmittance.
!                                UNITS:      None.
!                                TYPE:       Real
!                                DIMENSION:  I x K x M; K > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( OUT )
!
!       tau:                     Layer->TOA transmittance for the satellite
!                                view angle.
!                                UNITS:      None.
!                                TYPE:       Real
!                                DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( OUT )
!
!       flux_tau:                Layer->SFC transmittance for the default
!                                diffusivity angle.
!                                UNITS:      None.
!                                TYPE:       Real
!                                DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( OUT )
!
!       solar_tau:               Layer->SFC transmittance for the solar
!                                zenith angle.
!                                UNITS:      None.
!                                TYPE:       Real
!                                DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( OUT )
!
!       layer_radiance:          Layer Planck radiances at every layer for
!                                each channel/profile.
!                                UNITS:      mW/(m^2.sr.cm^-1)
!                                TYPE:       Real
!                                DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                ATTRIBUTES: INTENT( OUT )
!
!       downwelling_radiance:    TOA->SFC radiances for each channel/profile due
!                                to thermal flux and solar components.
!                                UNITS:      mW/(m^2.sr.cm^-1)
!                                TYPE:       Real
!                                DIMENSION:  L*M; L > 1, and M > or = 1
!                                            NB: This is a 1-D array.
!                                ATTRIBUTES: INTENT( OUT )
!
!       upwelling_radiance:      TOA radiances for each channel/profile.
!                                UNITS:      mW/(m^2.sr.cm^-1)
!                                TYPE:       Real
!                                DIMENSION:  L*M; L > 1, and M > or = 1
!                                            NB: This is a 1-D array.
!                                ATTRIBUTES: INTENT( OUT )
!
!       brightness_temperature:  Temperatures corresponding to the TOA radiances
!                                for each channel/profile.
!                                UNITS:      Kelvin
!                                TYPE:       Real
!                                DIMENSION:  L*M; L > 1, and M > or = 1
!                                            NB: This is a 1-D array.
!                                ATTRIBUTES: INTENT( OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result = SUCCESS => Calculation was successful
!              = FAILURE => Error occurred checking input arguments
!
! CALLS:
!      display_message:            Subroutine to output messages
!                                  SOURCE: error_handler module
!
!      get_max_n_channels:         Routine to retrieve the value of the
!                                  MAX_N_CHANNELS "pseudo-parameter".
!                                  SOURCE: parameters module
!
!      compute_absorber_amount:    Subroutine to integrate the absorber profiles
!                                  SOURCE: absorber_profile module
!
!      find_absorber_layer_index:  Subroutine to find the absorber space levels
!                                  that bracket the calculated absorber amounts
!                                  SOURCE: absorber_profile module
!                                  
!      compute_predictors:         Subroutine to compute the transmittance predictor
!                                  profiles.
!                                  SOURCE: predictor module
!
!      compute_transmittance:      Subroutine to compute the transmittance profiles.
!                                  SOURCE: transmittance module
!
!      compute_radiance:           Subroutine to compute the TOA radiances and
!                                  brightness temperatures.
!                                  SOURCE: radiance module
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
!       Note the restrictions on the input array dimensions:
!         K == n_layers > 1
!         L == n_channels > 1
!         M == n_profiles > or = 1
!
! PROCEDURE:
!       See individual module function documentation.
!S-
!--------------------------------------------------------------------------------

  FUNCTION forward_rtm_rank2( &

             ! -- Inputs
             level_p, layer_p, layer_t, layer_w, layer_o, &  ! Input, K x M

             surface_temperature,                         &  ! Input,  M
             surface_emissivity,                          &  ! Input,  L*M
             surface_reflectivity,                        &  ! Input,  L*M

             secant_view_angle,                           &  ! Input,  M
             secant_solar_angle,                          &  ! Input,  M
             n_channels_per_profile,                      &  ! Input,  M
             channel_index,                               &  ! Input,  L*M

             ! -- Outputs
             absorber,                                    &  ! Output, 0:K x J x M

             tau_layer_index,                             &  ! Output, K x J x M
             flux_tau_layer_index,                        &  ! Output, K x J x M
             solar_tau_layer_index,                       &  ! Output, K x J x M

             tau_predictor,                               &  ! Output, Imax x K x M
             flux_tau_predictor,                          &  ! Output, Imax x K x M
             solar_tau_predictor,                         &  ! Output, Imax x K x M

             tau,                                         &  ! Output, K x L*M
             flux_tau,                                    &  ! Output, K x L*M
             solar_tau,                                   &  ! Output, K x L*M

             layer_radiance,                              &  ! Output, K x L*M
             downwelling_radiance,                        &  ! Output, L*M
             upwelling_radiance,                          &  ! Output, L*M

             brightness_temperature,                      &  ! Output, L*M

             ! -- Optional inputs
             n_input_profiles,                            &
             message_log )                                &

           RESULT ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- TYPE DECLARATIONS --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Inputs
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: level_p                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_p                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_t                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_w                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_o                 ! K x M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: surface_temperature     ! M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: surface_emissivity      ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: surface_reflectivity    ! L*M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: secant_view_angle       ! M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: secant_solar_angle      ! M
    INTEGER,         DIMENSION( : ),        INTENT( IN )  :: n_channels_per_profile  ! M
    INTEGER,         DIMENSION( : ),        INTENT( IN )  :: channel_index           ! L*M

    ! -- Outputs
    REAL( fp_kind ), DIMENSION( 0:, :, : ), INTENT( OUT ) :: absorber                ! 0:K x J x M

    INTEGER,         DIMENSION(  :, :, : ), INTENT( OUT ) :: tau_layer_index         ! K x J x M
    INTEGER,         DIMENSION(  :, :, : ), INTENT( OUT ) :: flux_tau_layer_index    ! K x J x M
    INTEGER,         DIMENSION(  :, :, : ), INTENT( OUT ) :: solar_tau_layer_index   ! K x J x M

    REAL( fp_kind ), DIMENSION( :, :, : ),  INTENT( OUT ) :: tau_predictor           ! Imax x K x M
    REAL( fp_kind ), DIMENSION( :, :, : ),  INTENT( OUT ) :: flux_tau_predictor      ! Imax x K x M
    REAL( fp_kind ), DIMENSION( :, :, : ),  INTENT( OUT ) :: solar_tau_predictor     ! Imax x K x M

    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT ) :: tau                     ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT ) :: flux_tau                ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT ) :: solar_tau               ! K x L*M

    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT ) :: layer_radiance          ! K x L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( OUT ) :: downwelling_radiance    ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( OUT ) :: upwelling_radiance      ! L*M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( OUT ) :: brightness_temperature  ! L*M

    ! -- Optional input
    INTEGER,        OPTIONAL,               INTENT( IN )  :: n_input_profiles        ! Scalar
    CHARACTER( * ), OPTIONAL,               INTENT( IN )  :: message_log
    

    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'FORWARD_RTM_RANK2'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Scalars
    CHARACTER( 100 ) :: message
    CHARACTER( 10 )  :: value_in, value_allowed

    INTEGER :: m, n_profiles          ! Profile loop variables
    INTEGER :: l, l1, l2, n_channels  ! Channel loop/index variables

    INTEGER :: valid_solar

    ! -- Maximum channels pseudo parameter
    INTEGER :: MAX_N_CHANNELS
    LOGICAL :: is_set


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC ADJUSTL, &
              ANY,     &
              COS,     &
              MAXVAL,  &
              PRESENT, &
              SIZE,    &
              TRIM



    !#--------------------------------------------------------------------------#
    !#          -- DETERMINE ARRAY DIMENSIONS AND CHECK THE VALUES --           #
    !#--------------------------------------------------------------------------#

    ! ----------------------------
    ! Check the number of profiles
    ! ----------------------------

    ! -- Number of atmospheric profiles. Default size is the SECOND
    ! -- dimension of the LAYER_P argument.
    n_profiles = SIZE( layer_p, DIM = 2 )
    IF ( PRESENT( n_input_profiles ) ) THEN
      IF ( n_input_profiles > 0 .AND. n_input_profiles <= n_profiles ) THEN
        n_profiles = n_input_profiles
      ELSE
        WRITE( message, '( "Invalid N_INPUT_PROFILES value: ", i5, &
                          &". Using pressure array dimension value of ", i5, "." )' ) &
                        n_input_profiles, n_profiles
        CALL display_message( ROUTINE_NAME,    &
                              TRIM( message ), &
                              WARNING,         &
                              message_log = message_log )
      END IF
    END IF

    ! -- Check that the number of profiles is not greater than
    ! -- MAX_N_PROFILES. This is simply a limit to restrict the
    ! -- size of the input arrays so they're not TOO big.
    IF ( n_profiles > MAX_N_PROFILES ) THEN

      error_status = FAILURE
      WRITE( value_in,      '( i5 )' ) n_profiles
      WRITE( value_allowed, '( i5 )' ) MAX_N_PROFILES
      CALL display_message( ROUTINE_NAME, &
                            'Number of passed profiles ('// &
                            TRIM( ADJUSTL( value_in ) )// &
                            ') > maximum number of profiles allowed ('// &
                            TRIM( ADJUSTL( value_allowed ) )//').', &
                            error_status, &
                            message_log = message_log )
      RETURN

    END IF


    ! ----------------------------
    ! Check the number of channels
    ! ----------------------------

    ! -- Check for a negative number of channels
    IF ( ANY( n_channels_per_profile < 0 ) ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Elements of N_CHANNELS_PER_PROFILE are negative.', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF

    ! -- Check that the maximum number of channels for any profile
    ! -- is not greater than than the number of channels with which
    ! -- the model was initialised
    n_channels = MAXVAL( n_channels_per_profile )

    CALL get_max_n_channels( MAX_N_CHANNELS, is_set )

    IF ( .NOT. is_set ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'MAX_N_CHANNELS value not set. Check that RTM is initialised.', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF

    IF ( n_channels > MAX_N_CHANNELS ) THEN
      error_status = FAILURE
      WRITE( value_in,      '( i5 )' ) n_channels
      WRITE( value_allowed, '( i5 )' ) MAX_N_CHANNELS
      CALL display_message( ROUTINE_NAME, &
                            'Number of requested channels ('// &
                            TRIM( ADJUSTL( value_in ) )// &
                            ') > number of initialisation channels ('// &
                            TRIM( ADJUSTL( value_allowed ) )//').', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF




    !#--------------------------------------------------------------------------#
    !#                           -- PROFILE LOOP --                             #
    !#--------------------------------------------------------------------------#

    ! -------------------------------
    ! Initialise channel index finder
    ! -------------------------------

    l1 = 1


    ! ------------------
    ! Begin profile loop
    ! ------------------

    m_profile_loop: DO m = 1, n_profiles


      ! -----------------------------------------------
      ! Check for the "no channel" case. If no channels
      ! required for this profile, go to the next one.
      ! -----------------------------------------------

      IF ( n_channels_per_profile( m ) == 0 ) CYCLE m_profile_loop


      ! -------------------------------------------
      ! Determine the end channel index index range
      ! -------------------------------------------

      l2 = l1 + n_channels_per_profile( m ) - 1


      ! ------------------------
      ! Call the rank-1 function
      ! ------------------------

      error_status = forward_rtm_rank1( &

                       ! -- Forward inputs
                       level_p( :, m ),                   &  ! Input,  K
                       layer_p( :, m ),                   &  ! Input,  K
                       layer_t( :, m ),                   &  ! Input,  K
                       layer_w( :, m ),                   &  ! Input,  K
                       layer_o( :, m ),                   &  ! Input,  K

                       surface_temperature( m ),          &  ! Input,  Scalar
                       surface_emissivity( l1:l2 ),       &  ! Input,  L
                       surface_reflectivity( l1:l2 ),     &  ! Input,  L

                       ! -- Other inputs
                       secant_view_angle( m ),            &  ! Input,  Scalar
                       secant_solar_angle( m ),           &  ! Input,  Scalar
                       n_channels_per_profile( m ),       &  ! Input,  Scalar
                       channel_index( l1:l2 ),            &  ! Input,  L

                       ! -- Outputs
                       absorber( 0:, :, m ),              &  ! Output, 0:K x J

                       tau_layer_index( :, :, m ),        &  ! Output, K x J
                       flux_tau_layer_index( :, :, m ),   &  ! Output, K x J
                       solar_tau_layer_index( :, :, m ),  &  ! Output, K x J

                       tau_predictor( :, :, m ),          &  ! Output, Imax x K
                       flux_tau_predictor( :, :, m ),     &  ! Output, Imax x K
                       solar_tau_predictor( :, :, m ),    &  ! Output, Imax x K

                       tau( :, l1:l2 ),                   &  ! Output, K x L
                       flux_tau( :, l1:l2 ),              &  ! Output, K x L
                       solar_tau( :, l1:l2 ),             &  ! Output, K x L

                       layer_radiance( :, l1:l2 ),        &  ! Output, K x L
                       downwelling_radiance( l1:l2 ),     &  ! Output, L
                       upwelling_radiance( l1:l2 ),       &  ! Output, L

                       brightness_temperature( l1:l2 ),   &  ! Output, L

                       message_log = message_log )


      ! -------------------------------
      ! Check for successful completion
      ! -------------------------------

      IF ( error_status /= SUCCESS ) THEN

        CALL display_message( ROUTINE_NAME, &
                              'Error occured in FORWARD_RTM_RANK1', &
                              error_status, &
                              message_log = message_log )
        RETURN

      END IF


      ! ------------------------------
      ! Update the channel begin index
      ! ------------------------------

      l1 = l2 + 1

    END DO m_profile_loop



    !#--------------------------------------------------------------------------#
    !#                              -- DONE --                                  #
    !#--------------------------------------------------------------------------#

    error_status = SUCCESS

  END FUNCTION forward_rtm_rank2


  FUNCTION forward_rtm_rank1( &

             ! -- Inputs
             level_p, layer_p, layer_t, layer_w, layer_o, &  ! Input, K

             surface_temperature,                         &  ! Input,  Scalar
             surface_emissivity,                          &  ! Input,  L
             surface_reflectivity,                        &  ! Input,  L

             secant_view_angle,                           &  ! Input,  Scalar
             secant_solar_angle,                          &  ! Input,  Scalar
             n_channels,                                  &  ! Input,  Scalar
             channel_index,                               &  ! Input,  L

             ! -- Outputs
             absorber,                                    &  ! Output, 0:K x J

             tau_layer_index,                             &  ! Output, K x J
             flux_tau_layer_index,                        &  ! Output, K x J
             solar_tau_layer_index,                       &  ! Output, K x J

             tau_predictor,                               &  ! Output, Imax x K
             flux_tau_predictor,                          &  ! Output, Imax x K
             solar_tau_predictor,                         &  ! Output, Imax x K

             tau,                                         &  ! Output, K x L
             flux_tau,                                    &  ! Output, K x L
             solar_tau,                                   &  ! Output, K x L

             layer_radiance,                              &  ! Output, K x L
             downwelling_radiance,                        &  ! Output, L
             upwelling_radiance,                          &  ! Output, L

             brightness_temperature,                      &  ! Output, L

             ! -- Optional inputs
             n_input_profiles,                            &
             message_log )                                &

           RESULT ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- TYPE DECLARATIONS --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Inputs
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: level_p                 ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_p                 ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_t                 ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_w                 ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_o                 ! K

    REAL( fp_kind ),                     INTENT( IN )  :: surface_temperature     ! Scalar
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: surface_emissivity      ! L
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: surface_reflectivity    ! L

    REAL( fp_kind ),                     INTENT( IN )  :: secant_view_angle       ! Scalar
    REAL( fp_kind ),                     INTENT( IN )  :: secant_solar_angle      ! Scalar
    INTEGER,                             INTENT( IN )  :: n_channels              ! Scalar
    INTEGER,         DIMENSION( : ),     INTENT( IN )  :: channel_index           ! L

    ! -- Outputs
    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( OUT ) :: absorber                ! 0:K x J

    INTEGER,         DIMENSION(  :, : ), INTENT( OUT ) :: tau_layer_index         ! K x J
    INTEGER,         DIMENSION(  :, : ), INTENT( OUT ) :: flux_tau_layer_index    ! K x J
    INTEGER,         DIMENSION(  :, : ), INTENT( OUT ) :: solar_tau_layer_index   ! K x J

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: tau_predictor           ! Imax x K
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: flux_tau_predictor      ! Imax x K
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: solar_tau_predictor     ! Imax x K

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: tau                     ! K x L
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: flux_tau                ! K x L
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: solar_tau               ! K x L

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: layer_radiance          ! K x L
    REAL( fp_kind ), DIMENSION( : ),     INTENT( OUT ) :: downwelling_radiance    ! L
    REAL( fp_kind ), DIMENSION( : ),     INTENT( OUT ) :: upwelling_radiance      ! L

    REAL( fp_kind ), DIMENSION( : ),     INTENT( OUT ) :: brightness_temperature  ! L

    ! -- Optional input. Note that N_INPUT_PROFILES is not used in this
    ! -- function. It is specified so if a user specifies it by mistake
    ! -- for rank-1 profile input the code won't (hopefully) fall in a heap.
    INTEGER,        OPTIONAL,            INTENT( IN )  :: n_input_profiles        ! Scalar
    CHARACTER( * ), OPTIONAL,            INTENT( IN )  :: message_log
    

    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'FORWARD_RTM_RANK1'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Scalars
    CHARACTER( 100 ) :: message
    CHARACTER( 10 )  :: value_in, value_allowed

    INTEGER :: n_layers   ! Layer dimension
    INTEGER :: l          ! Channel loop/index variables

    INTEGER :: valid_solar

    ! -- Maximum channels pseudo parameter
    INTEGER :: MAX_N_CHANNELS
    LOGICAL :: is_set

    ! -- Arrays for integrated absorber amounts.
    REAL( fp_kind ), DIMENSION( 0:SIZE( absorber, DIM = 1 )-1, &
                                  SIZE( absorber, DIM = 2 )    ) :: tau_absorber,      &
                                                                    flux_tau_absorber, &
                                                                    solar_tau_absorber

    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC ADJUSTL, &
              ANY,     &
              COS,     &
              MAXVAL,  &
              SIZE,    &
              TRIM



    !#--------------------------------------------------------------------------#
    !#           -- DETERMINE ARRAY DIMENSIONS AND CHECK INPUT --               #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------
    ! Check the number of channels - if zero
    ! then simply RETURN.
    ! --------------------------------------

    IF ( n_channels == 0 ) RETURN


    ! ------------------
    ! Get the dimensions
    ! ------------------

    n_layers = SIZE( layer_p )


    ! ------------------------------------------------------------
    ! Check that the number of layers is not greater than
    ! MAX_N_ABSORBER_LAYERS. If there are more input layers
    ! there may be non-unique assignation of absorber->atmospheric
    ! layers
    ! ------------------------------------------------------------

    IF ( n_layers > MAX_N_ABSORBER_LAYERS ) THEN

      error_status = FAILURE
      WRITE( value_in,      '( i5 )' ) n_layers
      WRITE( value_allowed, '( i5 )' ) MAX_N_ABSORBER_LAYERS
      CALL display_message( ROUTINE_NAME, &
                            'Number of passed layers ('// &
                            TRIM( ADJUSTL( value_in ) )// &
                            ') > maximum number of absorber layers ('// &
                            TRIM( ADJUSTL( value_allowed ) )//').', &
                            error_status, &
                            message_log = message_log )
      RETURN

    END IF


    ! -----------------------------------
    ! Perform a simple check on the input
    ! data for negative values
    ! -----------------------------------

    ! -- Profile data
    IF ( ANY( level_p < ZERO ) .OR. &
         ANY( layer_p < ZERO ) .OR. &
         ANY( layer_t < ZERO ) .OR. &
         ANY( layer_w < ZERO ) .OR. &
         ANY( layer_o < ZERO )      ) THEN

       write (6,*) 'level_p', level_p
       write (6,*) 'layer_p', layer_p
       write (6,*) 'layer_t', layer_t
       write (6,*) 'layer_w', layer_w
       write (6,*) 'layer_o', layer_o
       
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Negative values found in input profile data.', &
                            error_status, &
                            message_log = message_log )
      
    END IF

    ! -- Surface properties
    IF (      surface_temperature  < ZERO   .OR. &
         ANY( surface_emissivity   < ZERO ) .OR. &
         ANY( surface_reflectivity < ZERO )      ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Negative values found in surface properties (Tsfc,esfc,rsfc).', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF

    ! -- Number of channels
    IF ( n_channels < 0 ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Negative number of channels passed..', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#           -- CALCULATE THE PROFILE GENERIC ABSORBER AMOUNTS --           #
    !#--------------------------------------------------------------------------#

    CALL compute_absorber_amount( level_p( : ),     &  ! Input,  K
                                  layer_w( : ),     &  ! Input,  K
                                  layer_o( : ),     &  ! Input,  K

                                  absorber( 0:, : ) )  ! Output, 0:K x J



    !#--------------------------------------------------------------------------#
    !#      -- CALCULATE THE PREDICTORS FOR THE UPWELLING TRANSMITTANCE --      #
    !#--------------------------------------------------------------------------#

    ! ------------------------------------------------------
    ! Modify absorber quantities by the angle secant
    ! Could put a loop here but here's hoping the compiler
    ! recognises this as a group of loops over layer.
    ! ------------------------------------------------------

    tau_absorber( 0:, : ) = secant_view_angle * absorber( 0:, : )

    ! -- Subtract the top pressure from the dry absorber profile
    tau_absorber( 1:, 2 ) = tau_absorber( 1:, 2 ) - TOA_PRESSURE


    ! -----------------------------------------------------
    ! Calculate the predictors for the satellite view angle
    ! -----------------------------------------------------

    CALL compute_predictors( layer_p( : ),          &  ! Input,  K
                             layer_t( : ),          &  ! Input,  K
                             layer_w( : ),          &  ! Input,  K
                             tau_absorber( 0:, : ), &  ! Input,  0:K x J

                             tau_predictor( :, : )  )  ! Output, I x K


    ! ------------------------------------------------
    ! Determine the absorber space levels that bracket
    ! the "average absorber" amounts at the view angle
    ! ------------------------------------------------

    CALL find_absorber_layer_index( tau_absorber( 0:, : ),  &  ! Input, 0:K x J
                                    tau_layer_index( :, : ) )  ! Output, K x J



    !#--------------------------------------------------------------------------#
    !#        -- CALCULATE THE PREDICTORS FOR THE FLUX TRANSMITTANCE --         #
    !#--------------------------------------------------------------------------#

    ! ---------------------------------------------
    ! Have any INFRARED channels been specified for
    ! the current profile? (Microwave channels are
    ! flagged as == 1, so IR == 0).
    !
    ! For microwave channels the downwelling flux
    ! transmission is assumed == upwelling view
    ! angle transmission.
    ! ---------------------------------------------

    IF ( ANY( is_microwave_channel( channel_index( 1:n_channels ) ) == 0 ) ) THEN


      ! ---------------------------------
      ! Modify the nadir absorber amounts
      ! ---------------------------------

      flux_tau_absorber( 0:, : ) = SECANT_DIFFUSIVITY_ANGLE * absorber( 0:, : )
      
      ! -- Subtract the top pressure from the dry absorber profile
      flux_tau_absorber( 1:, 2 ) = flux_tau_absorber( 1:, 2 ) - TOA_PRESSURE


      ! --------------------------------
      ! Calculate the predictors for the
      ! diffusivity angle
      ! --------------------------------

      ! -- Calculate the integrated predictors only
      CALL compute_predictors( layer_p( : ),               &  ! Input,  K
                               layer_t( : ),               &  ! Input,  K
                               layer_w( : ),               &  ! Input,  K
                               flux_tau_absorber( 0:, : ), &  ! Input,  0:K x J

                               flux_tau_predictor( :, : ), &  ! Output, I x K

                               no_standard = 1             )  ! Optional input

      ! -- Copy the angle independent (standard) predictors
      flux_tau_predictor( 1:MAX_N_STANDARD_PREDICTORS, : ) = &
           tau_predictor( 1:MAX_N_STANDARD_PREDICTORS, : ) 


      ! ------------------------------------------------
      ! Determine the absorber space levels that bracket
      ! the absorber amounts at the diffusivity angle
      ! ------------------------------------------------

      CALL find_absorber_layer_index( flux_tau_absorber( 0:, : ),  &  ! Input,  0:K x J
                                      flux_tau_layer_index( :, : ) )  ! Output, K x J

    END IF
     


    !#--------------------------------------------------------------------------#
    !#       -- CALCULATE THE PREDICTORS FOR THE SOLAR TRANSMITTANCE --         #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------------------
    ! Have *any* SOLAR sensitive channels been specified
    ! for the current profile (Flagged as == 1)?
    !
    ! AND
    !
    ! Is the specified solar zenith angle valid?
    ! --------------------------------------------------

    IF ( ( ANY( is_solar_channel( channel_index( 1:n_channels ) ) == 1 ) ) .AND. &
         secant_solar_angle < MAX_SECANT_SOLAR_ANGLE ) THEN


      ! --------------------------------
      ! Modify the nadir absorber amount
      ! --------------------------------

      solar_tau_absorber( 0:, : ) = secant_solar_angle * absorber( 0:, : )

      ! -- Subtract the top pressure from the dry absorber profile
      solar_tau_absorber( 1:, 2 ) = solar_tau_absorber( 1:, 2 ) - TOA_PRESSURE


      ! --------------------------------
      ! Calculate the predictors for the
      ! solar zenith angle
      ! --------------------------------

      ! -- Calculate the integrated predictors only
      CALL compute_predictors( layer_p( : ),                &  ! Input,  K
                               layer_t( : ),                &  ! Input,  K
                               layer_w( : ),                &  ! Input,  K
                               solar_tau_absorber( 0:, : ), &  ! Input,  0:K x J

                               solar_tau_predictor( :, : ), &  ! Output, I x K

                               no_standard = 1              )  ! Optional input

      ! -- Copy the angle independent predictors
      solar_tau_predictor( 1:MAX_N_STANDARD_PREDICTORS, : ) = &
            tau_predictor( 1:MAX_N_STANDARD_PREDICTORS, : ) 


      ! ------------------------------------------------
      ! Determine the absorber space levels that bracket
      ! the absorber amounts at the solar zenith angle
      ! ------------------------------------------------

      CALL find_absorber_layer_index( solar_tau_absorber( 0:, : ),  &  ! Input, 0:K x J
                                      solar_tau_layer_index( :, : ) )  ! Output, K x J

    END IF
     


    !#--------------------------------------------------------------------------#
    !#                           -- CHANNEL LOOP --                             #
    !#--------------------------------------------------------------------------#

    l_channel_loop: DO l = 1, n_channels


      ! --------------------------------------------------
      ! Calculate the current channel layer transmittances
      ! for the satellite view angle
      ! --------------------------------------------------

      CALL compute_transmittance( tau_absorber( 0:, : ),   &   ! Input, 0:K x J
                                  tau_predictor( :, : ),   &   ! Input, I x K
                                  tau_layer_index( :, : ), &   ! Input, K x J
                                  channel_index( l ),      &   ! Input, scalar
                                  UP,                      &   ! Input, scalar

                                  tau( :, l )              )   ! Output, K


      ! -----------------------------------------------------
      ! If the current channel is an INFRARED channel,
      ! then calculate the downwelling flux transmittance.
      !
      ! If the current channel is a MICROWAVE channel,
      ! then use the predictors for the upwelling
      ! transmittance calculations.
      !
      ! Two things:
      ! - Currently, the predictors and coefficients are
      !   the same for the up- and downwelling cases,
      !   hence the simple "upending" of the transmittances
      !   for the microwave (assumed specular) case.
      ! - The "upending", or assumption that the downwelling
      !   transmittances can be derived directly from the
      !   upwelling transmittances, is only valid for
      !   monochromatic transmittances. For broadband sensors
      !   this is an approximation - and depending on the
      !   spectral structure of the transmittances within
      !   a channel's response - not usually a good one.
      ! -----------------------------------------------------

      IF ( is_microwave_channel( channel_index( l ) ) == 0 ) THEN

        ! -- IR channel
        CALL compute_transmittance( flux_tau_absorber( 0:, : ),   &   ! Input, 0:K x J
                                    flux_tau_predictor( :, : ),   &   ! Input, I x K
                                    flux_tau_layer_index( :, : ), &   ! Input, K x J
                                    channel_index( l ),           &   ! Input, scalar
                                    DOWN,                         &   ! Input, scalar

                                    flux_tau( :, l )              )   ! Output, K

      ELSE


        ! -- uW channel

!  This can be considered a "hook" for future versions where
!  downwelling will not be derived from the upwelling.
!
!        CALL compute_transmittance( tau_absorber( 0:, : ),   &   ! Input, 0:K x J
!                                    tau_predictor( :, : ),   &   ! Input, I x K
!                                    tau_layer_index( :, : ), &   ! Input, K x J
!                                    channel_index( l ),      &   ! Input, scalar
!                                    DOWN,                    &   ! Input, scalar
!
!                                    flux_tau( :, l )         )   ! Output, K

        ! -- This gives the identical result as a separate call
        ! -- but without the extra computational burden.
        flux_tau( :, l ) = ZERO
        IF ( tau( n_layers, l ) > TOLERANCE ) THEN
          flux_tau( 1, l ) = tau( n_layers, l )
          flux_tau( 2:n_layers, l ) =     flux_tau( 1, l ) / &
          !                           ------------------------
                                       tau( 1:n_layers-1, l )
        END IF

      END IF



      ! ----------------------------------------------------
      ! If the current channel is a SOLAR SENSITIVE channel,
      !   AND
      ! the solar angle is valid, then calculate the
      ! transmittance for direct solar.
      ! ----------------------------------------------------

      IF ( is_solar_channel( channel_index( l ) ) == 1 .AND. &
           secant_solar_angle < MAX_SECANT_SOLAR_ANGLE ) THEN

        valid_solar = 1

        CALL compute_transmittance( solar_tau_absorber( 0:, : ),   &   ! Input, 0:K x J
                                    solar_tau_predictor( :, : ),   &   ! Input, I x K
                                    solar_tau_layer_index( :, : ), &   ! Input, K x J
                                    channel_index( l ),            &   ! Input, scalar
                                    DOWN,                          &   ! Input, scalar

                                    solar_tau( :, l )              )   ! Output, K

      ELSE

        valid_solar = 0

        solar_tau( :, l ) = ZERO

      END IF


      ! ---------------------------------------
      ! Calculate the profile/channel radiances
      ! ---------------------------------------

      CALL compute_radiance( layer_t( : ),               &  ! Input, K

                             surface_temperature,        &  ! Input, scalar
                             surface_emissivity( l ),    &  ! Input, scalar
                             surface_reflectivity( l ),  &  ! Input, scalar

                             tau(      :, l ),           &  ! Input, K
                             flux_tau( :, l ),           &  ! Input, K
                             solar_tau( n_layers, l ),   &  ! Input, scalar

                             secant_solar_angle,         &  ! Input, scalar
                             valid_solar,                &  ! Input, scalar
                             channel_index( l ),         &  ! Input, scalar

                             layer_radiance( :, l ),     &  ! Output, K
                             downwelling_radiance( l ),  &  ! Output, scalar
                             upwelling_radiance( l ),    &  ! Output, scalar

                             brightness_temperature( l ) )  ! Output, scalar

    END DO l_channel_loop



    !#--------------------------------------------------------------------------#
    !#                              -- DONE --                                  #
    !#--------------------------------------------------------------------------#

    error_status = SUCCESS

  END FUNCTION forward_rtm_rank1

END MODULE forward_model


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
! Revision 2.12  2001/11/07 15:03:15  paulv
! - Added check for negative number of channels to FORWARD_RTM_RANK2() function.
! - Added profile loop CYCLE statement to FORWARD_RTM_RANK2() function.
! - Added check for negative number of channels to FORWARD_RTM_RANK1() function.
! - Changed the logical IF test for the number of input profiles in the
!   FORWARD_RTM_RANK2() function from
!     IF ( n_input_profiles >= 1 .AND. n_input_profiles <= n_profiles ) THEN
!   to
!     IF ( n_input_profiles > 0 .AND. n_input_profiles <= n_profiles ) THEN
!   The use of both the ">=" and "<=" realtional operators with the .AND. I
!   found confusing.
!
! Revision 2.11  2001/09/28 22:44:24  paulv
! - Overloaded the COMPUTE_RTM and FORWARD_RTM functions to accept both a
!   single or group of profiles. Contained PRIVATE functions are now
!     o COMPUTE_RTM_RANK1 and FORWARD_RTM_RANK1 for single profile input
!     o COMPUTE_RTM_RANK2 and FORWARD_RTM_RANK2 for multiple profile input
! - Put N_INPUT_PROFILES optional argument back in the COMPUTE_RTM and FORWARD_RTM
!   argument lists. This allows more control over how many profiles are to be
!   processed rather than simply relying on the dimension of the input arrays.
!   Now, as before,
!     n_profiles = SIZE( layer_p, DIM = 2 )
!   but also,
!     IF ( PRESENT( n_input_profiles ) ) THEN
!       IF ( n_input_profiles >= 1 .AND. n_input_profiles <= n_profiles ) THEN
!         n_profiles = n_input_profiles
!     ....
! - Changed SURFACE_TEMPERATURE argument check from
!     IF ( ANY( surface_temperature < ZERO ) )
!   to
!     IF (      surface_temperature < ZERO )
!   as the check is now done in the rank-1 FORWARD_RTM function. This eliminates
!   the need to fully populate the input arrays with data when only certain
!   "chunks" may be processed. Previously, the use of ANY() could generate
!   and error if the full surface_temperature array was not initialised.
! - Added "Name" to RCS keyword list.
!
! Revision 2.10  2001/09/04 21:29:11  paulv
! - Updated documentation.
!
! Revision 2.9  2001/08/31 21:22:23  paulv
! - Removed input data checks from COMPUTE_RTM. The same checks are performed
!   in the main routine, FORWARD_RTM, so there was no need to replicate them.
! - Added check for negative profile and surface data in FORWARD_RTM.
! - Maximum solar angle secant is no longer calculated in FORWARD_RTM but
!   is declared as a parameter in the PARAMETERS module.
!
! Revision 2.8  2001/08/16 16:36:29  paulv
! - Updated documentation.
! - The channel dimension is now obtained by:
!     n_channels = MAXVAL( n_channels_per_profile )
!   rather than
!     n_channels = SIZE( channel_index ) / n_profiles
!   since the latter assumes the same number of channels will be processed
!   for each profile - which may not be the case. The new method determines
!   the largest number of channels to be processed for any particular
!   profile.
! - The comparison of n_channels and MAX_N_CHANNELS is now done via the
!   MAX_N_CHANNELS methods in the PARAMETERS module. And extra check to
!   see if MAX_N_CHANNELS has been set was added.
!
! Revision 2.7  2001/08/01 16:47:59  paulv
! - Removed USE of ABSORBER_SPACE module.
! - ONLY clauses added to other USE statements so only those module members
!   required in this module are visible.
! - Added COMPUTE_RTM function. This is a wrapper for the FORWARD_RTM function.
! - Updated input argument checking. Now consistent with other model
!   components.
!
! Revision 2.6  2001/07/12 18:41:37  paulv
! - Commented out informational message output at start of function.
!
! Revision 2.5  2001/05/29 18:21:02  paulv
! - Now use TYPE_KINDS module parameter FP_KIND to set the floating point
!   data type.
! - Added ABSORBER_SPACE and PREDICTORS module use to provide access to
!   the absorber space definitions and the predictor calculation routines.
! - All angle parameters placed in the PARAMETERS module.
! - Argument list increased to provide data for other component RTM calls
!   (e.g. TL or AD).
! - For full listing of differences, do a cvs diff between revisions 2.4
!   and 2.5 (this one).
!
! Revision 2.4  2001/01/24 20:14:21  paulv
! - Latest test versions.
!
! Revision 2.3  2000/11/09 20:58:28  paulv
! - Made radiance function call consistent with changes to that module.
!   (specifically regarding the solar, flux, and reflectivity terms).
!
! Revision 2.2  2000/08/31 19:36:32  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 2.1  2000/08/24 18:08:28  paulv
! - Many changes from initial version. Too many to list here.
! - Updated module and subprogram documentation.
!
!
!
!
!
