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
!                               ! -- Forward inputs
!                               level_p, layer_p, layer_t, layer_w, layer_o, &  ! Input, K x M
!
!                               surface_temperature,                         &  ! Input, M
!                               surface_emissivity,                          &  ! Input, L*M
!                               surface_reflectivity,                        &  ! Input, L*M
!
!                               ! -- Other inputs
!                               secant_view_angle,                           &  ! Input, M
!                               secant_solar_angle,                          &  ! Input, M
!                               n_channels_per_profile,                      &  ! Input, M
!                               channel_index,                               &  ! Input, L*M
!
!                               ! -- Forward output
!                               tau,                                         &  ! Input, K x L*M
!                               flux_tau,                                    &  ! Input, K x L*M
!                               solar_tau,                                   &  ! Input, K x L*M
!
!                               upwelling_radiance,                          &  ! Input, L*M
!                               brightness_temperature,                      &  ! Input, L*M
!
!                               ! Optional inputs
!                               message_log = message_log )
!
! INPUT ARGUMENTS:
!
!       level_p:                   Profile set layer interface pressure array. The TOA
!                                  pressure is not included. TOA pressure is parameterised
!                                  in the PARAMETERS module.
!                                  UNITS:      hPa
!                                  TYPE:       Real
!                                  DIMENSION:  K x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_p:                   Profile set layer average pressure array.
!                                  UNITS:      hPa
!                                  TYPE:       Real
!                                  DIMENSION:  K x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_t:                   Profile set layer average temperature array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  K x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_w:      .            Profile set layer average water vapor mixing ratio array
!                                  UNITS:      g/kg
!                                  TYPE:       Real
!                                  DIMENSION:  K x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_o:                   Profile set layer average ozone mixing ratio array.
!                                  UNITS:      ppmv
!                                  TYPE:       Real
!                                  DIMENSION:  K x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_temperature:       Profile set surface temperature array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  M
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_emissivity:        Profile set surface emissivity array
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_reflectivity:      Profile set surface reflectivity array
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!       secant_view_angle:         Secant of the satellite view angle measured
!                                  from nadir for each profile in the set.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  M
!                                  ATTRIBUTES: INTENT( IN )
!
!       secant_solar_angle:        Secant of the solar zenith angle for each
!                                  profile in the set.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  M
!                                  ATTRIBUTES: INTENT( IN )
!
!       n_channels_per_profile:    The number of channels for each profile in the
!                                  set for which radiances are required.
!                                  UNITS:      None
!                                  TYPE:       Integer
!                                  DIMENSION:  M
!                                  ATTRIBUTES: INTENT( IN )
!
!       channel_index:             Channel index id array. Each element is a unique
!                                  index to a (supported) sensor channel.
!                                  UNITS:      None
!                                  TYPE:       Integer
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!
! OPTIONAL INPUT ARGUMENTS:
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
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( OUT )
!
!       flux_tau:                  Layer->SFC transmittance for the default
!                                  diffusivity angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( OUT )
!
!       solar_tau:                 Layer->SFC transmittance for the solar
!                                  zenith angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( OUT )
!
!       upwelling_radiance:        TOA radiances for each channel/profile.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( OUT )
!
!       brightness_temperature:    TOA brightness temperatures corresponding
!                                  to the TOA radiances.
!                                  N.B.: Set to ZERO upon output.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
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
!       The code has is not overloaded for scalar input so the input
!       arguments must be dimensioned accordingly, even if only one
!       profile or channel is being passed.
!
! PROCEDURE:
!       See individual module function documentation.
!S-
!--------------------------------------------------------------------------------

  FUNCTION compute_rtm( &
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
             tau,                                         &  ! Input, K x L*M
             flux_tau,                                    &  ! Input, K x L*M
             solar_tau,                                   &  ! Input, K x L*M

             upwelling_radiance,                          &  ! Input, L*M
             brightness_temperature,                      &  ! Input, L*M

             ! Optional inputs
             message_log )                                &

           RESULT ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward inputs
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )     :: level_p                    ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )     :: layer_p                    ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )     :: layer_t                    ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )     :: layer_w                    ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )     :: layer_o                    ! K x M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )     :: surface_temperature        ! M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )     :: surface_emissivity         ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )     :: surface_reflectivity       ! L*M

    ! -- Other inputs
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )     :: secant_view_angle          ! M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )     :: secant_solar_angle         ! M
    INTEGER,         DIMENSION( : ),        INTENT( IN )     :: n_channels_per_profile     ! M
    INTEGER,         DIMENSION( : ),        INTENT( IN )     :: channel_index              ! L*M

    ! -- Forward outputs
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT )    :: tau                        ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT )    :: flux_tau                   ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT )    :: solar_tau                  ! K x L*M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( OUT )    :: upwelling_radiance         ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( OUT )    :: brightness_temperature     ! L*M

    ! -- Optional input
    CHARACTER( * ), OPTIONAL,               INTENT( IN )     :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_RTM'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Scalars
    CHARACTER( 100 ) :: message
    CHARACTER( 5 )   :: value_in, value_allowed

    INTEGER :: m, n_profiles          ! Profile loop variables
    INTEGER :: l, l1, l2, n_channels  ! Channel loop/index variables
    INTEGER :: k, n_layers            ! Layer loop variables

    ! -- Maximum channels pseudo parameter
    INTEGER :: MAX_N_CHANNELS
    LOGICAL :: is_set

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

    error_status = forward_rtm( &
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

                                message_log = message_log )


    ! -------------------------------
    ! Check for successful completion
    ! -------------------------------

    IF ( error_status /= SUCCESS ) THEN

      CALL display_message( ROUTINE_NAME, &
                            'Error occured in FORWARD_RTM', &
                            error_status, &
                            message_log = message_log )
      RETURN

    END IF



    !#--------------------------------------------------------------------------#
    !#                              -- Done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = SUCCESS


  END FUNCTION compute_rtm





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
!                             message_log = message_log )
!                             
! INPUT ARGUMENTS:
!
!       level_p:                 Profile set layer interface pressure array. The TOA
!                                pressure is not included. TOA pressure is parameterised
!                                in the PARAMETERS module.
!                                UNITS:      hPa
!                                TYPE:       Real
!                                DIMENSION:  K x M
!                                ATTRIBUTES: INTENT( IN )
!
!       layer_p:                 Profile set layer average pressure array.
!                                UNITS:      hPa
!                                TYPE:       Real
!                                DIMENSION:  K x M
!                                ATTRIBUTES: INTENT( IN )
!
!       layer_t:                 Profile set layer average temperature array.
!                                UNITS:      Kelvin
!                                TYPE:       Real
!                                DIMENSION:  K x M
!                                ATTRIBUTES: INTENT( IN )
!
!       layer_w:      .          Profile set layer average water vapor mixing ratio array
!                                UNITS:      g/kg
!                                TYPE:       Real
!                                DIMENSION:  K x M
!                                ATTRIBUTES: INTENT( IN )
!
!       layer_o:                 Profile set layer average ozone mixing ratio array.
!                                UNITS:      ppmv
!                                TYPE:       Real
!                                DIMENSION:  K x M
!                                ATTRIBUTES: INTENT( IN )
!
!       surface_temperature:     Profile set surface temperature array.
!                                UNITS:      Kelvin
!                                TYPE:       Real
!                                DIMENSION:  M
!                                ATTRIBUTES: INTENT( IN )
!
!       surface_emissivity:      Profile set surface emissivity array
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  L*M
!                                            NB: This is a 1-D array.
!                                ATTRIBUTES: INTENT( IN )
!
!       surface_reflectivity:    Profile set surface reflectivity array
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  L*M
!                                            NB: This is a 1-D array.
!                                ATTRIBUTES: INTENT( IN )
!
!       secant_view_angle:       Secant of the satellite view angle measured
!                                from nadir for each profile in the set.
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  M
!                                ATTRIBUTES: INTENT( IN )
!
!       secant_solar_angle:      Secant of the solar zenith angle for each
!                                profile in the set.
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  M
!                                ATTRIBUTES: INTENT( IN )
!
!       n_channels_per_profile:  The number of channels for each profile in the
!                                set for which radiances are required.
!                                UNITS:      None
!                                TYPE:       Integer
!                                DIMENSION:  M
!                                ATTRIBUTES: INTENT( IN )
!
!       channel_index:           Channel index id array. Each element is a unique
!                                index to a (supported) sensor channel.
!                                UNITS:      None
!                                TYPE:       Integer
!                                DIMENSION:  L*M
!                                            NB: This is a 1-D array.
!                                ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!
!       message_log:           Character string specifying a filename in which any
!                              messages will be logged. If not specified, or if an
!                              error occurs opening the log file, the default action
!                              is to output messages to the screen.
!                              UNITS:      None
!                              TYPE:       Character
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! OUTPUT ARGUMENTS:
!
!       absorber:               Array of absorber amount for nadir view.
!                               UNITS:      Absorber dependent.
!                               TYPE:       Real
!                               DIMENSION:  0:K x J x M
!                               ATTRIBUTES: INTENT( OUT )
!
!       tau_layer_index:        Array of absorber space layer indices of the input
!                               absorber amounts at the satellite view angle.
!                               UNITS:      None.
!                               TYPE:       Integer
!                               DIMENSION:  K x J x M
!                               ATTRIBUTES: INTENT( OUT )
!
!       flux_tau_layer_index:   Array of absorber space layer indices of the input
!                               absorber amounts at the default diffusivity angle.
!                               UNITS:      None.
!                               TYPE:       Integer
!                               DIMENSION:  K x J x M
!                               ATTRIBUTES: INTENT( OUT )
!
!       solar_tau_layer_index:  Array of absorber space layer indices of the input
!                               absorber amounts at the solar zenith angle.
!                               UNITS:      None.
!                               TYPE:       Integer
!                               DIMENSION:  K x J x M
!                               ATTRIBUTES: INTENT( OUT )
!
!       tau_predictor:          Predictor profiles for the layer->TOA transmittance.
!                               UNITS:      None.
!                               TYPE:       Real
!                               DIMENSION:  I x K x M
!                               ATTRIBUTES: INTENT( OUT )
!
!       flux_tau_predictor:     Predictor profiles for the thermal flux transmittance.
!                               UNITS:      None.
!                               TYPE:       Real
!                               DIMENSION:  I x K x M
!                               ATTRIBUTES: INTENT( OUT )
!
!       solar_tau_predictor:    Predictor profiles for the solar transmittance.
!                               UNITS:      None.
!                               TYPE:       Real
!                               DIMENSION:  I x K x M
!                               ATTRIBUTES: INTENT( OUT )
!
!       tau:                    Layer->TOA transmittance for the satellite
!                               view angle.
!                               UNITS:      None.
!                               TYPE:       Real
!                               DIMENSION:  K x L*M
!                               ATTRIBUTES: INTENT( OUT )
!
!       flux_tau:               Layer->SFC transmittance for the default
!                               diffusivity angle.
!                               UNITS:      None.
!                               TYPE:       Real
!                               DIMENSION:  K x L*M
!                               ATTRIBUTES: INTENT( OUT )
!
!       solar_tau:              Layer->SFC transmittance for the solar
!                               zenith angle.
!                               UNITS:      None.
!                               TYPE:       Real
!                               DIMENSION:  K x L*M
!                               ATTRIBUTES: INTENT( OUT )
!
!       layer_radiance:         Layer Planck radiances at every layer for
!                               each channel/profile.
!                               UNITS:      mW/(m^2.sr.cm^-1)
!                               TYPE:       Real
!                               DIMENSION:  K x L*M
!                               ATTRIBUTES: INTENT( OUT )
!
!       downwelling_radiance:   TOA->SFC radiances for each channel/profile due
!                               to thermal flux and solar components.
!                               UNITS:      mW/(m^2.sr.cm^-1)
!                               TYPE:       Real
!                               DIMENSION:  L*M
!                                           NB: This is a 1-D array.
!                               ATTRIBUTES: INTENT( OUT )
!
!       upwelling_radiance:     TOA radiances for each channel/profile.
!                               UNITS:      mW/(m^2.sr.cm^-1)
!                               TYPE:       Real
!                               DIMENSION:  L*M
!                                           NB: This is a 1-D array.
!                               ATTRIBUTES: INTENT( OUT )
!
!       brightness_temperature: Temperatures corresponding to the TOA radiances
!                               for each channel/profile.
!                               UNITS:      Kelvin
!                               TYPE:       Real
!                               DIMENSION:  L*M
!                                           NB: This is a 1-D array.
!                               ATTRIBUTES: INTENT( OUT )
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
!       The code has is not overloaded for scalar input so the input
!       arguments must be dimensioned accordingly, even if only one
!       profile or channel is being passed.
!
! PROCEDURE:
!       See individual module function documentation.
!S-
!--------------------------------------------------------------------------------

  FUNCTION forward_rtm( &

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
             message_log )                                &

           RESULT ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
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
    CHARACTER( * ), OPTIONAL,               INTENT( IN )  :: message_log
    

    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'FORWARD_RTM'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Scalars
    CHARACTER( 100 ) :: message
    CHARACTER( 10 )  :: value_in, value_allowed

    INTEGER :: m, n_profiles          ! Profile loop variables
    INTEGER :: l, l1, l2, n_channels  ! Channel loop/index variables
    INTEGER :: k, n_layers            ! Layer loop variables

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
    !#                   -- Determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    ! ------------------
    ! Get the dimensions
    ! ------------------

    n_layers   = SIZE( layer_p, DIM = 1 )
    n_profiles = SIZE( layer_p, DIM = 2 )
    n_channels = MAXVAL( n_channels_per_profile ) ! Result is scalar for rank-1 array


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


    ! ------------------------------------------------------
    ! Check that the number of profiles is not greater than
    ! MAX_N_PROFILES. This is simply a limit to restrict the
    ! size of the input arrays so they're not TOO big.
    ! ------------------------------------------------------

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


    ! -----------------------------------------------------
    ! Check that the number of channels is not greater than
    ! than the number of channels with which the model was
    ! initialised
    ! -----------------------------------------------------

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
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Negative values found in input profile data.', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF

    ! -- Surface properties
    IF ( ANY( surface_temperature  < ZERO ) .OR. &
         ANY( surface_emissivity   < ZERO ) .OR. &
         ANY( surface_reflectivity < ZERO )      ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Negative values found in surface properties (Tsfc,esfc,rsfc).', &
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


      ! -------------------------------------------
      ! Determine the end channel index index range
      ! -------------------------------------------

      l2 = l1 + n_channels_per_profile( m ) - 1



      !#------------------------------------------------------------------------#
      !#          -- Calculate the profile generic absorber amounts --          #
      !#------------------------------------------------------------------------#

      CALL compute_absorber_amount( level_p( :, m ),     &  ! Input,  K
                                    layer_w( :, m ),     &  ! Input,  K
                                    layer_o( :, m ),     &  ! Input,  K

                                    absorber( 0:, :, m ) )  ! Output, 0:K x J



      !#------------------------------------------------------------------------#
      !#     -- Calculate the predictors for the upwelling transmittance --     #
      !#------------------------------------------------------------------------#

      ! ------------------------------------------------------
      ! Modify absorber quantities by the angle secant
      ! Could put a loop here but here's hoping the compiler
      ! recognises this as a group of loops over layer.
      ! ------------------------------------------------------

      tau_absorber( 0:, : ) = secant_view_angle( m ) * absorber( 0:, :, m )

      ! -- Subtract the top pressure from the dry absorber profile
      tau_absorber( 1:, 2 ) = tau_absorber( 1:, 2 ) - TOA_PRESSURE


      ! -----------------------------------------------------
      ! Calculate the predictors for the satellite view angle
      ! -----------------------------------------------------

      CALL compute_predictors( layer_p( :, m ),         &  ! Input,  K
                               layer_t( :, m ),         &  ! Input,  K
                               layer_w( :, m ),         &  ! Input,  K
                               tau_absorber( 0:, : ),   &  ! Input,  0:K x J

                               tau_predictor( :, :, m ) )  ! Output, I x K


      ! ------------------------------------------------
      ! Determine the absorber space levels that bracket
      ! the "average absorber" amounts at the view angle
      ! ------------------------------------------------

      CALL find_absorber_layer_index( tau_absorber( 0:, : ),     &  ! Input, 0:K x J
                                      tau_layer_index( :, :, m ) )  ! Output, K x J



      !#------------------------------------------------------------------------#
      !#       -- Calculate the predictors for the flux transmittance --        #
      !#------------------------------------------------------------------------#

      ! ---------------------------------------------
      ! Have any INFRARED channels been specified for
      ! the current profile? (Microwave channels are
      ! flagged as == 1, so IR == 0).
      !
      ! For microwave channels the downwelling flux
      ! transmission is assumed == upwelling view
      ! angle transmission.
      ! ---------------------------------------------

      IF ( ANY( is_microwave_channel( channel_index( l1:l2 ) ) == 0 ) ) THEN


        ! ---------------------------------
        ! Modify the nadir absorber amounts
        ! ---------------------------------

        flux_tau_absorber( 0:, : ) = SECANT_DIFFUSIVITY_ANGLE * absorber( 0:, :, m )
        
        ! -- Subtract the top pressure from the dry absorber profile
        flux_tau_absorber( 1:, 2 ) = flux_tau_absorber( 1:, 2 ) - TOA_PRESSURE


        ! --------------------------------
        ! Calculate the predictors for the
        ! diffusivity angle
        ! --------------------------------

        ! -- Calculate the integrated predictors only
        CALL compute_predictors( layer_p( :, m ),               &  ! Input,  K
                                 layer_t( :, m ),               &  ! Input,  K
                                 layer_w( :, m ),               &  ! Input,  K
                                 flux_tau_absorber( 0:, : ),    &  ! Input,  0:K x J

                                 flux_tau_predictor( :, :, m ), &  ! Output, I x K

                                 no_standard = 1                )  ! Optional input

        ! -- Copy the angle independent (standard) predictors
        flux_tau_predictor( 1:MAX_N_STANDARD_PREDICTORS, :, m ) = &
             tau_predictor( 1:MAX_N_STANDARD_PREDICTORS, :, m ) 


        ! ------------------------------------------------
        ! Determine the absorber space levels that bracket
        ! the absorber amounts at the diffusivity angle
        ! ------------------------------------------------

        CALL find_absorber_layer_index( flux_tau_absorber( 0:, : ),     &  ! Input,  0:K x J
                                        flux_tau_layer_index( :, :, m ) )  ! Output, K x J

      END IF
       


      !#------------------------------------------------------------------------#
      !#       -- Calculate the predictors for the solar transmittance --       #
      !#------------------------------------------------------------------------#

      ! --------------------------------------------------
      ! Have *any* SOLAR sensitive channels been specified
      ! for the current profile (Flagged as == 1)?
      !
      ! AND
      !
      ! Is the specified solar zenith angle valid (<85)?
      ! --------------------------------------------------

      IF ( ( ANY( is_solar_channel( channel_index( l1:l2 ) ) == 1 ) ) .AND. &
           secant_solar_angle( m ) < MAX_SECANT_SOLAR_ANGLE ) THEN


        ! --------------------------------
        ! Modify the nadir absorber amount
        ! --------------------------------

        solar_tau_absorber( 0:, : ) = secant_solar_angle( m ) * absorber( 0:, :, m )

        ! -- Subtract the top pressure from the dry absorber profile
        solar_tau_absorber( 1:, 2 ) = solar_tau_absorber( 1:, 2 ) - TOA_PRESSURE


        ! --------------------------------
        ! Calculate the predictors for the
        ! solar zenith angle
        ! --------------------------------

        ! -- Calculate the integrated predictors only
        CALL compute_predictors( layer_p( :, m ),                &  ! Input,  K
                                 layer_t( :, m ),                &  ! Input,  K
                                 layer_w( :, m ),                &  ! Input,  K
                                 solar_tau_absorber( 0:, : ),    &  ! Input,  0:K x J

                                 solar_tau_predictor( :, :, m ), &  ! Output, I x K

                                 no_standard = 1                 )  ! Optional input

        ! -- Copy the angle independent predictors
        solar_tau_predictor( 1:MAX_N_STANDARD_PREDICTORS, :, m ) = &
              tau_predictor( 1:MAX_N_STANDARD_PREDICTORS, :, m ) 


        ! ------------------------------------------------
        ! Determine the absorber space levels that bracket
        ! the absorber amounts at the solar zenith angle
        ! ------------------------------------------------

        CALL find_absorber_layer_index( solar_tau_absorber( 0:, : ),     &  ! Input, 0:K x J
                                        solar_tau_layer_index( :, :, m ) )  ! Output, K x J

      END IF
       


      !#------------------------------------------------------------------------#
      !#                           -- CHANNEL LOOP --                           #
      !#------------------------------------------------------------------------#

      l_channel_loop: DO l = l1, l2


        ! --------------------------------------------------
        ! Calculate the current channel layer transmittances
        ! for the satellite view angle
        ! --------------------------------------------------

        CALL compute_transmittance( tau_absorber( 0:, : ),      &   ! Input, 0:K x J
                                    tau_predictor( :, :, m ),   &   ! Input, I x K
                                    tau_layer_index( :, :, m ), &   ! Input, K x J
                                    channel_index( l ),         &   ! Input, scalar
                                    UP,                         &   ! Input, scalar

                                    tau( :, l )                 )   ! Output, K


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
          CALL compute_transmittance( flux_tau_absorber( 0:, : ),      &   ! Input, 0:K x J
                                      flux_tau_predictor( :, :, m ),   &   ! Input, I x K
                                      flux_tau_layer_index( :, :, m ), &   ! Input, K x J
                                      channel_index( l ),              &   ! Input, scalar
                                      DOWN,                            &   ! Input, scalar

                                      flux_tau( :, l )                 )   ! Output, K

        ELSE


          ! -- uW channel

!  This can be considered a "hook" for future versions where
!  downwelling will not be derived from the upwelling.
!
!          CALL compute_transmittance( tau_absorber( 0:, : ),      &   ! Input, 0:K x J
!                                      tau_predictor( :, :, m ),   &   ! Input, I x K
!                                      tau_layer_index( :, :, m ), &   ! Input, K x J
!                                      channel_index( l ),         &   ! Input, scalar
!                                      DOWN,                       &   ! Input, scalar
!
!                                      flux_tau( :, l )            )   ! Output, K

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
             secant_solar_angle( m ) < MAX_SECANT_SOLAR_ANGLE ) THEN

          valid_solar = 1

          CALL compute_transmittance( solar_tau_absorber( 0:, : ),      &   ! Input, 0:K x J
                                      solar_tau_predictor( :, :, m ),   &   ! Input, I x K
                                      solar_tau_layer_index( :, :, m ), &   ! Input, K x J
                                      channel_index( l ),               &   ! Input, scalar
                                      DOWN,                             &   ! Input, scalar

                                      solar_tau( :, l )                 )   ! Output, K

        ELSE

          valid_solar = 0

          solar_tau( :, l ) = ZERO

        END IF


        ! ---------------------------------------
        ! Calculate the profile/channel radiances
        ! ---------------------------------------

        CALL compute_radiance( layer_t( :, m ),            &  ! Input, K

                               surface_temperature( m ),   &  ! Input, scalar
                               surface_emissivity( l ),    &  ! Input, scalar
                               surface_reflectivity( l ),  &  ! Input, scalar

                               tau(      :, l ),           &  ! Input, K
                               flux_tau( :, l ),           &  ! Input, K
                               solar_tau( n_layers, l ),   &  ! Input, scalar

                               secant_solar_angle( m ),    &  ! Input, scalar
                               valid_solar,                &  ! Input, scalar
                               channel_index( l ),         &  ! Input, scalar

                               layer_radiance( :, l ),     &  ! Output, K
                               downwelling_radiance( l ),  &  ! Output, scalar
                               upwelling_radiance( l ),    &  ! Output, scalar

                               brightness_temperature( l ) )  ! Output, scalar

      END DO l_channel_loop


      ! ------------------------------
      ! Update the channel begin index
      ! ------------------------------

      l1 = l2 + 1

    END DO m_profile_loop



    !#--------------------------------------------------------------------------#
    !#                              -- Done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = SUCCESS

  END FUNCTION forward_rtm

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
! $State$
!
! $Log$
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
