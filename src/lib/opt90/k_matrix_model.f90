!------------------------------------------------------------------------------
!M+
! NAME:
!       k_matrix_model
!
! PURPOSE:
!       Module containing the NCEP RT K-matrix model functions
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE k_matrix_model
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
!       forward_model:         Module containing the forward model function.
!
! CONTAINS:
!       compute_rtm_K:         PUBLIC function that calculates the K-matrix of the 
!                              top-of-atmosphere (TOA) radiances and brightness 
!                              temperatures for an input atmospheric profile set and
!                              user specified satellites/channels.
!
!                              This function is simply a wrapper around both the FORWARD
!                              model and the K-MATRIX model so that the user doesn't have
!                              to declare the absorber/predictor/etc. arrays in the calling
!                              routine.
!
!       k_matrix_rtm:          PUBLIC function that calculates the K-matrix of the
!                              top-of-atmosphere (TOA) radiances and brightness
!                              temperatures for user specified profiles and
!                              satellite/channel transmittance profiles, radiances
!                              and brightness temperatures.
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
!       Written by:     Paul van Delst, CIMSS@NOAA/NCEP 18-Jul-2001
!                       pvandelst@ncep.noaa.gov
!
!  Copyright (C) 2001 Paul van Delst
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

MODULE k_matrix_model


  ! ------------
  ! Module usage
  ! ------------

  USE type_kinds,            ONLY : fp_kind
  USE error_handler
  USE parameters
  USE spectral_coefficients, ONLY : is_solar_channel, &
                                    is_microwave_channel
  USE absorber_profile,      ONLY : compute_absorber_amount_AD
  USE predictors,            ONLY : compute_predictors_AD
  USE transmittance,         ONLY : compute_transmittance_AD
  USE radiance,              ONLY : compute_radiance_AD
  USE forward_model,         ONLY : forward_rtm


  ! -----------------------
  ! Disable implicit typing
  ! -----------------------

  IMPLICIT NONE


  ! --------------------
  ! Default visibilities
  ! --------------------

  PRIVATE
  PUBLIC :: compute_rtm_K
  PUBLIC :: k_matrix_rtm


CONTAINS





!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_rtm_K
!
! PURPOSE:
!       PUBLIC function that calculates the K-matrix of the top-of-atmosphere (TOA)
!       radiances and brightness temperatures for an input atmospheric profile
!       set and user specified satellites/channels.
!
!       This function is simply a wrapper around both the FORWARD model and the
!       K-MATRIX model so that the user doesn't have to declare the absorber/
!       predictor/etc arrays in the calling routine.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       result = compute_rtm_K( &
!                               ! -- Forward inputs
!                               level_p, layer_p, layer_t, layer_w, layer_o,           &  ! Input, K x M
!
!                               surface_temperature,                                   &  ! Input, M
!                               surface_emissivity,                                    &  ! Input, L*M
!                               surface_reflectivity,                                  &  ! Input, L*M
!
!                               ! -- K-matrix inputs
!                               tau_K,                                                 &  ! In/Output, K x L*M
!                               flux_tau_K,                                            &  ! In/Output, K x L*M
!                               solar_tau_K,                                           &  ! In/Output, K x L*M
!
!                               upwelling_radiance_K,                                  &  ! In/Output, L*M
!                               brightness_temperature_K,                              &  ! In/Output, L*M
!
!                               ! -- Other inputs
!                               secant_view_angle,                                     &  ! Input, M
!                               secant_solar_angle,                                    &  ! Input, M
!                               n_channels_per_profile,                                &  ! Input, M
!                               channel_index,                                         &  ! Input, L*M
!
!                               ! -- Forward output
!                               tau,                                                   &  ! Input, K x L*M
!                               flux_tau,                                              &  ! Input, K x L*M
!                               solar_tau,                                             &  ! Input, K x L*M
!
!                               upwelling_radiance,                                    &  ! Input, L*M
!                               brightness_temperature,                                &  ! Input, L*M
!
!                               ! -- K-matrix outputs
!                               level_p_K, layer_p_K, layer_t_K, layer_w_K, layer_o_K, &  ! In/Output, K x L*M
!
!                               surface_temperature_K,                                 &  ! In/Output, L*M
!                               surface_emissivity_K,                                  &  ! In/Output, L*M
!                               surface_reflectivity_K,                                &  ! In/Output, L*M
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
!       tau_K:                     Layer->TOA adjoint transmittance for the satellite
!                                  view angle.
!                                  ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       flux_tau_K:                Layer->SFC adjoint transmittance for the default
!                                  diffusivity angle.
!                                  ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       solar_tau_K:               Layer->SFC adjoint transmittance for the solar
!                                  zenith angle.
!                                  ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       upwelling_radiance_K:      TOA adjoint radiances for each channel/profile.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       brightness_temperature_K:  Adjoint temperatures corresponding to the
!                                  TOA adjoint radiances.
!                                  ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
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
!       level_p_K:                 Profile set layer interface pressure k-matrix
!                                  adjoint array.
!                                  UNITS:      hPa
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                              NB: This is a 2-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       layer_p_K:                 Profile set layer average pressure k-matrix
!                                  adjoint array.
!                                  UNITS:      hPa
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                              NB: This is a 2-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       layer_t_K:                 Profile set layer average temperature k-matrix
!                                  adjoint array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                              NB: This is a 2-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       layer_w_K:      .          Profile set layer average water vapor mixing ratio
!                                  k-matrix adjoint array.
!                                  UNITS:      g/kg
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                              NB: This is a 2-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       layer_o_K:                 Profile set layer average ozone mixing ratio
!                                  k-matrix adjoint array.
!                                  UNITS:      ppmv
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                              NB: This is a 2-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       surface_temperature_K:     Profile set surface temperature k-matrix adjoint
!                                  array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       surface_emissivity_K:      Profile set surface emissivity k-matrix adjoint
!                                  array.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       surface_reflectivity_K:    Profile set surface reflectivity k-matrix adjoint
!                                  array.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
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
!                                  SOURCE: forward_model module
!
!      k_matrix_rtm:               Function that calculates the K-matrix of the (TOA)
!                                  radiances/temperatures.
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       All input adjoint arguments are set to ZERO on output.
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

  FUNCTION compute_rtm_K( &
             ! -- Forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o,           &  ! Input, K x M

             surface_temperature,                                   &  ! Input, M
             surface_emissivity,                                    &  ! Input, L*M
             surface_reflectivity,                                  &  ! Input, L*M

             ! -- K-matrix inputs
             tau_K,                                                 &  ! In/Output, K x L*M
             flux_tau_K,                                            &  ! In/Output, K x L*M
             solar_tau_K,                                           &  ! In/Output, K x L*M

             upwelling_radiance_K,                                  &  ! In/Output, L*M
             brightness_temperature_K,                              &  ! In/Output, L*M

             ! -- Other inputs
             secant_view_angle,                                     &  ! Input, M
             secant_solar_angle,                                    &  ! Input, M
             n_channels_per_profile,                                &  ! Input, M
             channel_index,                                         &  ! Input, L*M

             ! -- Forward output
             tau,                                                   &  ! Input, K x L*M
             flux_tau,                                              &  ! Input, K x L*M
             solar_tau,                                             &  ! Input, K x L*M

             upwelling_radiance,                                    &  ! Input, L*M
             brightness_temperature,                                &  ! Input, L*M

             ! -- K-matrix outputs
             level_p_K, layer_p_K, layer_t_K, layer_w_K, layer_o_K, &  ! In/Output, K x L*M

             surface_temperature_K,                                 &  ! In/Output, L*M
             surface_emissivity_K,                                  &  ! In/Output, L*M
             surface_reflectivity_K,                                &  ! In/Output, L*M

             ! Optional inputs
             message_log )                                          &

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

    ! -- K-matrix inputs
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: tau_K                      ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: flux_tau_K                 ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: solar_tau_K                ! K x L*M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: upwelling_radiance_K       ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: brightness_temperature_K   ! L*M

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

    ! -- K-matrix outputs
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: level_p_K                  ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: layer_p_K                  ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: layer_t_K                  ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: layer_w_K                  ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: layer_o_K                  ! K x L*M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: surface_temperature_K      ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: surface_emissivity_K       ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: surface_reflectivity_K     ! L*M

    ! -- Optional input
    CHARACTER( * ), OPTIONAL,               INTENT( IN )     :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_RTM_K'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Scalars
    CHARACTER( 100 ) :: message
    CHARACTER( 5 )   :: value_in, value_allowed

    INTEGER :: m, n_profiles          ! Profile loop variables
    INTEGER :: l, l1, l2, n_channels  ! Channel loop/index variables
    INTEGER :: k, n_layers            ! Layer loop variables

    INTEGER :: error_status_FWD
    INTEGER :: error_status_K

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

    ! -- Array for forward and K-matrix layer Planck radiance term, K x L*M
    REAL( fp_kind ), DIMENSION( SIZE( layer_p, DIM = 1 ),  &
                                SIZE( upwelling_radiance ) ) :: layer_radiance,  &
                                                                layer_radiance_K
      

    ! -- Array for forward and K-matrix downwelling radiance (flux + solar), L*M
    REAL( fp_kind ), DIMENSION( SIZE( upwelling_radiance ) ) :: downwelling_radiance,  &
                                                                downwelling_radiance_K 



    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC ADJUSTL, &
              MAXVAL,  &
              SIZE,    &
              TRIM



    !#--------------------------------------------------------------------------#
    !#                   -- Compute the forward radiances --                    #
    !#--------------------------------------------------------------------------#

    error_status_FWD = forward_rtm( &
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

    IF ( error_status_FWD /= SUCCESS ) THEN

      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Error occured in FORWARD_RTM', &
                            error_status, &
                            message_log = message_log )
      RETURN

    END IF



    !#--------------------------------------------------------------------------#
    !#                   -- Compute the K-matrix profiles --                    #
    !#--------------------------------------------------------------------------#

    ! -----------------------------------------------
    ! Initialise all local adjoint/K-matrix variables
    ! -----------------------------------------------

    layer_radiance_K( :, : )    = ZERO
    downwelling_radiance_K( : ) = ZERO


    ! -----------------------
    ! Call the k-matrix model
    ! -----------------------

    error_status_K = k_matrix_rtm( &
                                   ! -- Forward inputs
                                   level_p, layer_p, layer_t, layer_w, layer_o,           &  ! Input,  K x M

                                   surface_temperature,                                   &  ! Input, M
                                   surface_emissivity,                                    &  ! Input, L*M
                                   surface_reflectivity,                                  &  ! Input, L*M

                                   absorber,                                              &  ! Input, 0:K x J x M

                                   tau_layer_index,                                       &  ! Input, K x J x M
                                   flux_tau_layer_index,                                  &  ! Input, K x J x M
                                   solar_tau_layer_index,                                 &  ! Input, K x J x M

                                   tau_predictor,                                         &  ! Input, Imax x K x M
                                   flux_tau_predictor,                                    &  ! Input, Imax x K x M
                                   solar_tau_predictor,                                   &  ! Input, Imax x K x M

                                   tau,                                                   &  ! Input, K x L*M
                                   flux_tau,                                              &  ! Input, K x L*M
                                   solar_tau,                                             &  ! Input, K x L*M

                                   layer_radiance,                                        &  ! Input, K x L*M
                                   downwelling_radiance,                                  &  ! Input, L*M
                                   upwelling_radiance,                                    &  ! Input, L*M

                                   ! -- K-matrix inputs
                                   tau_K,                                                 &  ! In/Output, K x L*M
                                   flux_tau_K,                                            &  ! In/Output, K x L*M
                                   solar_tau_K,                                           &  ! In/Output, K x L*M

                                   layer_radiance_K,                                      &  ! In/Output, K x L*M
                                   downwelling_radiance_K,                                &  ! In/Output, L*M
                                   upwelling_radiance_K,                                  &  ! In/Output, L*M

                                   brightness_temperature_K,                              &  ! In/Output, L*M

                                   ! -- Other inputs
                                   secant_view_angle,                                     &  ! Input, M
                                   secant_solar_angle,                                    &  ! Input, M
                                   n_channels_per_profile,                                &  ! Input, M
                                   channel_index,                                         &  ! Input, L*M

                                   ! -- K-matrix outputs
                                   level_p_K, layer_p_K, layer_t_K, layer_w_K, layer_o_K, &  ! In/Output,  K x L*M

                                   surface_temperature_K,                                 &  ! In/Output, L*M
                                   surface_emissivity_K,                                  &  ! In/Output, L*M
                                   surface_reflectivity_K,                                &  ! In/Output, L*M

                                   message_log = message_log )


    ! -------------------------------
    ! Check for successful completion
    ! -------------------------------

    IF ( error_status_K /= SUCCESS ) THEN

      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Error occured in K_MATRIX_RTM', &
                            error_status, &
                            message_log = message_log )
      RETURN

    END IF


    !#--------------------------------------------------------------------------#
    !#                              -- Done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = SUCCESS


  END FUNCTION compute_rtm_K





!--------------------------------------------------------------------------------
!S+
! NAME:
!       k_matrix_rtm
!
! PURPOSE:
!       PUBLIC function that calculates the K-matrix of the top-of-atmosphere (TOA)
!       radiances and brightness temperatures for an input atmospheric profile
!       set and user specified satellites/channels.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       result = k_matrix_rtm( &
!                              ! -- Forward inputs
!                              level_p, layer_p, layer_t, layer_w, layer_o,           &  ! Input, K x M
!
!                              surface_temperature,                                   &  ! Input, M
!                              surface_emissivity,                                    &  ! Input, L*M
!                              surface_reflectivity,                                  &  ! Input, L*M
!
!                              absorber,                                              &  ! Input, 0:K x J x M
!
!                              tau_layer_index,                                       &  ! Input, K x J x M
!                              flux_tau_layer_index,                                  &  ! Input, K x J x M
!                              solar_tau_layer_index,                                 &  ! Input, K x J x M
!
!                              tau_predictor,                                         &  ! Input, Imax x K x M
!                              flux_tau_predictor,                                    &  ! Input, Imax x K x M
!                              solar_tau_predictor,                                   &  ! Input, Imax x K x M
!
!                              tau,                                                   &  ! Input, K x L*M
!                              flux_tau,                                              &  ! Input, K x L*M
!                              solar_tau,                                             &  ! Input, K x L*M
!
!                              layer_radiance,                                        &  ! Input, K x L*M
!                              downwelling_radiance,                                  &  ! Input, L*M
!                              upwelling_radiance,                                    &  ! Input, L*M
!
!                              ! -- K-matrix inputs
!                              tau_K,                                                 &  ! In/Output, K x L*M
!                              flux_tau_K,                                            &  ! In/Output, K x L*M
!                              solar_tau_K,                                           &  ! In/Output, K x L*M
!
!                              layer_radiance_K,                                      &  ! In/Output, K x L*M
!                              downwelling_radiance_K,                                &  ! In/Output, L*M
!                              upwelling_radiance_K,                                  &  ! In/Output, L*M
!
!                              brightness_temperature_K,                              &  ! In/Output, L*M
!
!                              ! -- Other inputs
!                              secant_view_angle,                                     &  ! Input, M
!                              secant_solar_angle,                                    &  ! Input, M
!                              n_channels_per_profile,                                &  ! Input, M
!                              channel_index,                                         &  ! Input, L*M
!
!                              ! -- K-matrix outputs
!                              level_p_K, layer_p_K, layer_t_K, layer_w_K, layer_o_K, &  ! In/Output, K x L*M
!
!                              surface_temperature_K,                                 &  ! In/Output, L*M
!                              surface_emissivity_K,                                  &  ! In/Output, L*M
!                              surface_reflectivity_K,                                &  ! In/Output, L*M
!
!                              ! Optional inputs
!                              message_log = message_log )
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
!       absorber:                  Array of absorber amount for nadir view.
!                                  UNITS:      Absorber dependent.
!                                  TYPE:       Real
!                                  DIMENSION:  0:K x J x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       tau_layer_index:           Array of absorber space layer indices of the input
!                                  absorber amounts at the satellite view angle.
!                                  UNITS:      None.
!                                  TYPE:       Integer
!                                  DIMENSION:  K x J x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       flux_tau_layer_index:      Array of absorber space layer indices of the input
!                                  absorber amounts at the default diffusivity angle.
!                                  UNITS:      None.
!                                  TYPE:       Integer
!                                  DIMENSION:  K x J x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       solar_tau_layer_index:     Array of absorber space layer indices of the input
!                                  absorber amounts at the solar zenith angle.
!                                  UNITS:      None.
!                                  TYPE:       Integer
!                                  DIMENSION:  K x J x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       tau_predictor:             Predictor profiles for the layer->TOA transmittance.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  I x K x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       flux_tau_predictor:        Predictor profiles for the thermal flux transmittance.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  I x K x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       solar_tau_predictor:       Predictor profiles for the solar transmittance.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  I x K x M
!                                  ATTRIBUTES: INTENT( IN )
!
!       tau:                       Layer->TOA transmittance for the satellite
!                                  view angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN )
!
!       flux_tau:                  Layer->SFC transmittance for the default
!                                  diffusivity angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN )
!
!       solar_tau:                 Layer->SFC transmittance for the solar
!                                  zenith angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_radiance:            Layer Planck radiances at every layer for
!                                  each channel/profile.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN )
!
!       downwelling_radiance:      TOA->SFC radiances for each channel/profile due
!                                  to thermal flux and solar components.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!       upwelling_radiance:        TOA radiances for each channel/profile.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!       tau_K:                     Layer->TOA adjoint transmittance for the satellite
!                                  view angle.
!                                  N.B.: Set to ZERO upon output.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       flux_tau_K:                Layer->SFC adjoint transmittance for the default
!                                  diffusivity angle.
!                                  N.B.: Set to ZERO upon output.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       solar_tau_K:               Layer->SFC adjoint transmittance for the solar
!                                  zenith angle.
!                                  N.B.: Set to ZERO upon output.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       layer_radiance_K:          Layer Planck adjoint radiances at every layer for
!                                  each channel/profile.
!                                  N.B.: Set to ZERO upon output.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       downwelling_radiance_K:    TOA->SFC adjoint radiances for each channel/profile due
!                                  to thermal flux and solar components.
!                                  N.B.: Set to ZERO upon output.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       upwelling_radiance_K:      TOA adjoint radiances for each channel/profile.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  N.B.: Set to ZERO upon output.
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       brightness_temperature_K:  Adjoint temperatures corresponding to the
!                                  TOA adjoint radiances.
!                                  N.B.: Set to ZERO upon output.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
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
!       level_p_K:                 Profile set layer interface pressure k-matrix
!                                  adjoint array.
!                                  UNITS:      hPa
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                              NB: This is a 2-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       layer_p_K:                 Profile set layer average pressure k-matrix
!                                  adjoint array.
!                                  UNITS:      hPa
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                              NB: This is a 2-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       layer_t_K:                 Profile set layer average temperature k-matrix
!                                  adjoint array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                              NB: This is a 2-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       layer_w_K:      .          Profile set layer average water vapor mixing ratio
!                                  k-matrix adjoint array.
!                                  UNITS:      g/kg
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                              NB: This is a 2-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       layer_o_K:                 Profile set layer average ozone mixing ratio
!                                  k-matrix adjoint array.
!                                  UNITS:      ppmv
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M
!                                              NB: This is a 2-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       surface_temperature_K:     Profile set surface temperature k-matrix adjoint
!                                  array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       surface_emissivity_K:      Profile set surface emissivity k-matrix adjoint
!                                  array.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
!
!       surface_reflectivity_K:    Profile set surface reflectivity k-matrix adjoint
!                                  array.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  L*M
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN OUT )
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
!      compute_absorber_amount_AD: Subroutine to calculate the adjoint
!                                  absorber profiles
!                                  SOURCE: absorber_profile module
!
!      compute_predictors_AD:      Subroutine to compute the adjoint 
!                                  transmittance predictor profiles.
!                                  SOURCE: predictor module
!
!      compute_transmittance_AD:   Subroutine to compute the adjoint
!                                  transmittance profiles.
!                                  SOURCE: transmittance module
!
!      compute_radiance_AD:        Subroutine to compute the TOA adjoint 
!                                  radiances and brightness temperatures.
!                                  SOURCE: radiance module
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

  FUNCTION k_matrix_rtm( &

             ! -- Forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o,           &  ! Input, K x M

             surface_temperature,                                   &  ! Input, M
             surface_emissivity,                                    &  ! Input, L*M
             surface_reflectivity,                                  &  ! Input, L*M

             absorber,                                              &  ! Input, 0:K x J x M

             tau_layer_index,                                       &  ! Input, K x J x M
             flux_tau_layer_index,                                  &  ! Input, K x J x M
             solar_tau_layer_index,                                 &  ! Input, K x J x M

             tau_predictor,                                         &  ! Input, Imax x K x M
             flux_tau_predictor,                                    &  ! Input, Imax x K x M
             solar_tau_predictor,                                   &  ! Input, Imax x K x M

             tau,                                                   &  ! Input, K x L*M
             flux_tau,                                              &  ! Input, K x L*M
             solar_tau,                                             &  ! Input, K x L*M

             layer_radiance,                                        &  ! Input, K x L*M
             downwelling_radiance,                                  &  ! Input, L*M
             upwelling_radiance,                                    &  ! Input, L*M

             ! -- K-matrix inputs
             tau_K,                                                 &  ! In/Output, K x L*M
             flux_tau_K,                                            &  ! In/Output, K x L*M
             solar_tau_K,                                           &  ! In/Output, K x L*M

             layer_radiance_K,                                      &  ! In/Output, K x L*M
             downwelling_radiance_K,                                &  ! In/Output, L*M
             upwelling_radiance_K,                                  &  ! In/Output, L*M

             brightness_temperature_K,                              &  ! In/Output, L*M

             ! -- Other inputs
             secant_view_angle,                                     &  ! Input, M
             secant_solar_angle,                                    &  ! Input, M
             n_channels_per_profile,                                &  ! Input, M
             channel_index,                                         &  ! Input, L*M

             ! -- K-matrix outputs
             level_p_K, layer_p_K, layer_t_K, layer_w_K, layer_o_K, &  ! In/Output, K x L*M

             surface_temperature_K,                                 &  ! In/Output, L*M
             surface_emissivity_K,                                  &  ! In/Output, L*M
             surface_reflectivity_K,                                &  ! In/Output, L*M

             ! Optional inputs
             message_log )                                          &

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

    REAL( fp_kind ), DIMENSION( 0:, :, : ), INTENT( IN )     :: absorber                   ! 0:K x J x M

    INTEGER,         DIMENSION(  :, :, : ), INTENT( IN )     :: tau_layer_index            ! K x J x M
    INTEGER,         DIMENSION(  :, :, : ), INTENT( IN )     :: flux_tau_layer_index       ! K x J x M
    INTEGER,         DIMENSION(  :, :, : ), INTENT( IN )     :: solar_tau_layer_index      ! K x J x M

    REAL( fp_kind ), DIMENSION( :, :, : ),  INTENT( IN )     :: tau_predictor              ! Imax x K x M
    REAL( fp_kind ), DIMENSION( :, :, : ),  INTENT( IN )     :: flux_tau_predictor         ! Imax x K x M
    REAL( fp_kind ), DIMENSION( :, :, : ),  INTENT( IN )     :: solar_tau_predictor        ! Imax x K x M

    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )     :: tau                        ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )     :: flux_tau                   ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )     :: solar_tau                  ! K x L*M

    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )     :: layer_radiance             ! K x L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )     :: downwelling_radiance       ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )     :: upwelling_radiance         ! L*M

    ! -- K-matrix inputs
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: tau_K                      ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: flux_tau_K                 ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: solar_tau_K                ! K x L*M

    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: layer_radiance_K           ! K x L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: downwelling_radiance_K     ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: upwelling_radiance_K       ! L*M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: brightness_temperature_K   ! L*M

    ! -- Other inputs
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )     :: secant_view_angle          ! M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )     :: secant_solar_angle         ! M
    INTEGER,         DIMENSION( : ),        INTENT( IN )     :: n_channels_per_profile     ! M
    INTEGER,         DIMENSION( : ),        INTENT( IN )     :: channel_index              ! L*M

    ! -- K-matrix outputs
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: level_p_K                  ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: layer_p_K                  ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: layer_t_K                  ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: layer_w_K                  ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN OUT ) :: layer_o_K                  ! K x L*M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: surface_temperature_K      ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: surface_emissivity_K       ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN OUT ) :: surface_reflectivity_K     ! L*M

    ! -- Optional input
    CHARACTER( * ), OPTIONAL,               INTENT( IN )     :: message_log
    
    
    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'K_MATRIX_RTM'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Scalars
    CHARACTER( 100 ) :: message
    CHARACTER( 5 )   :: value_in, value_allowed

    INTEGER :: m, n_profiles          ! Profile loop variables
    INTEGER :: l, l1, l2, n_channels  ! Channel loop/index variables
    INTEGER :: k, n_layers            ! Layer loop variables

    INTEGER :: valid_solar

    ! -- Maximum channels pseudo parameter
    INTEGER :: MAX_N_CHANNELS
    LOGICAL :: is_set

    ! -- Arrays for integrated absorber amounts, 0:K x J
    REAL( fp_kind ), DIMENSION( 0:SIZE( absorber, DIM = 1 )-1, &
                                  SIZE( absorber, DIM = 2 )    ) :: tau_absorber,      &
                                                                    flux_tau_absorber, &
                                                                    solar_tau_absorber

    ! -- Arrays for K-matrix adjoints of integrated absorber amounts, 0:K x J
    REAL( fp_kind ), DIMENSION( 0:SIZE( absorber, DIM = 1 )-1, &
                                  SIZE( absorber, DIM = 2 )    ) :: absorber_K, &
                                                                    tau_absorber_K,      &
                                                                    flux_tau_absorber_K, &
                                                                    solar_tau_absorber_K

    ! -- Arrays for K-matrix adjoint predictors, Imax x K
    REAL( fp_kind ), DIMENSION( SIZE( tau_predictor, DIM = 1 ), &
                                SIZE( tau_predictor, DIM = 2 )  ) :: tau_predictor_K,      &
                                                                     flux_tau_predictor_K, &
                                                                     solar_tau_predictor_K

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
    ! profile data for negative values
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


      ! ----------------------------------------------------
      ! Modify absorber quantities by the angle secant
      ! Could put a loop here but here's hoping the compiler
      ! recognises this as a group of loops over layer.
      ! ----------------------------------------------------

      ! -- Upwelling transmittance
      tau_absorber( 0:, : ) = secant_view_angle( m ) * absorber( 0:, :, m )
      tau_absorber( 1:, 2 ) = tau_absorber( 1:, 2 ) - TOA_PRESSURE

      ! -- Flux transmittance
      IF ( ANY( is_microwave_channel( channel_index( l1:l2 ) ) == 0 ) ) THEN
        flux_tau_absorber( 0:, : ) = SECANT_DIFFUSIVITY_ANGLE * absorber( 0:, :, m )
        flux_tau_absorber( 0:, : ) = secant_view_angle( m ) * absorber( 0:, :, m )
        flux_tau_absorber( 1:, 2 ) = flux_tau_absorber( 1:, 2 ) - TOA_PRESSURE
      END IF

      ! -- Solar transmittance
      IF ( ( ANY( is_solar_channel( channel_index( l1:l2 ) ) == 1 ) ) .AND. &
           secant_solar_angle( m ) < MAX_SECANT_SOLAR_ANGLE ) THEN
        solar_tau_absorber( 0:, : ) = secant_solar_angle( m ) * absorber( 0:, :, m )
        solar_tau_absorber( 1:, 2 ) = solar_tau_absorber( 1:, 2 ) - TOA_PRESSURE
      END IF


      !#------------------------------------------------------------------------#
      !#                           -- CHANNEL LOOP --                           #
      !#------------------------------------------------------------------------#

      l_channel_loop: DO l = l1, l2


        ! ---------------------------------------------------------
        ! Initialise local, channel dependent, adjoint variables
        !
        ! This initialisation may slow things down a bit if 
        ! a) the compiler isn't too clever about optimising and
        !    generates code that loops over each array separately
        ! b) there are a large number of channels.
        !
        ! The solution to (a) above is put in explicit loops, which
        ! may be necessary.
        !
        ! The alternative to (b) above is to define the arrays
        ! with a channel dimension which could be problematical
        ! memorywise for lots of channels.
        !
        ! These will eventually be moved outside both the channel
        ! and profile loop as once they are initialised, they are
        ! set to zero in the relevant adjoint routine on output.
        !----------------------------------------------------------

        ! -- Absorber arrays, 0:K x J
        absorber_K( 0:, : )           = ZERO
        tau_absorber_K( 0:, : )       = ZERO
        flux_tau_absorber_K( 0:, : )  = ZERO
        solar_tau_absorber_K( 0:, : ) = ZERO


        ! -- Predictor arrays, Imax x K
        tau_predictor_K( :, : )       = ZERO
        flux_tau_predictor_K( :, : )  = ZERO
        solar_tau_predictor_K( :, : ) = ZERO


        ! ----------------------------------------------------
        ! Set the "this is a channel influenced by solar" flag
        ! ----------------------------------------------------

        valid_solar = 0

        IF ( is_solar_channel( channel_index( l ) ) == 1 .AND. &
             secant_solar_angle( m ) < MAX_SECANT_SOLAR_ANGLE        ) valid_solar = 1



        ! ------------------------------------------------
        ! Calculate the adjoint of the current channel TOA
        ! radiance or brightness temperature
        ! ------------------------------------------------

        CALL compute_radiance_AD( &
                                  ! -- Forward input
                                  layer_t( :, m ),               &  ! Input, K

                                  surface_temperature( m ),      &  ! Input, scalar
                                  surface_emissivity( l ),       &  ! Input, scalar
                                  surface_reflectivity( l ),     &  ! Input, scalar

                                  tau(      :, l ),              &  ! Input, K
                                  flux_tau( :, l ),              &  ! Input, K
                                  solar_tau( n_layers, l ),      &  ! Input, scalar

                                  layer_radiance( :, l ),        &  ! Input, K
                                  downwelling_radiance( l ),     &  ! Input, scalar
                                  upwelling_radiance( l ),       &  ! Input, scalar

                                  ! -- K-matrix input
                                  layer_radiance_K( :, l ),      &  ! In/Output, K
                                  downwelling_radiance_K( l ),   &  ! In/Output, scalar
                                  upwelling_radiance_K( l ),     &  ! In/Output, scalar

                                  brightness_temperature_K( l ), &  ! In/Output, scalar

                                  ! -- Other input
                                  secant_solar_angle( m ),       &  ! Input, scalar
                                  valid_solar,                   &  ! Input, scalar
                                  channel_index( l ),            &  ! Input, scalar

                                  ! -- K-matrix output
                                  layer_t_K( :, l ),             &  ! In/Output, K

                                  surface_temperature_K( l ),    &  ! In/Output, scalar
                                  surface_emissivity_K( l ),     &  ! In/Output, scalar
                                  surface_reflectivity_K( l ),   &  ! In/Output, scalar

                                  tau_K( :, l ),                 &  ! In/Output, K
                                  flux_tau_K( :, l ),            &  ! In/Output, K
                                  solar_tau_K( n_layers, l )     )  ! In/Output, scalar



        !#----------------------------------------------------------------------#
        !# --    Calculate the K-matrix adjoint result for the solar term    -- #
        !#----------------------------------------------------------------------#

        ! ----------------------------------------------------
        ! If the current channel is a SOLAR SENSITIVE channel,
        !   AND
        ! the solar angle is valid, then calculate the adjoint
        ! of the transmittance for direct solar.
        ! ----------------------------------------------------

        solar_if_block: IF ( valid_solar == 1 ) THEN


          ! -----------------------------------------
          ! Adjoint of the direct solar transmittance
          ! -----------------------------------------

          CALL compute_transmittance_AD( &
                                         ! -- Forward input
                                         solar_tau_absorber( 0:, : ),      &   ! Input, 0:K x J
                                         solar_tau_predictor( :, :, m ),   &   ! Input, I x K
                                         solar_tau( :, l ),                &   ! Input, K

                                         ! -- K-matrix input
                                         solar_tau_K( :, l ),              &   ! In/Output, K

                                         ! -- Other input
                                         solar_tau_layer_index( :, :, m ), &   ! Input, K x J
                                         channel_index( l ),               &   ! Input, scalar
                                         DOWN,                             &   ! Input, scalar

                                         ! -- K-matrix output
                                         solar_tau_absorber_K( 0:, : ),    &   ! In/Output, 0:K x J
                                         solar_tau_predictor_K( :, : )     )   ! In/Output, I x K


          ! -----------------------------------------------
          ! K-matrix adjoint of the standard predictor copy
          ! -----------------------------------------------

          tau_predictor_K( 1:MAX_N_STANDARD_PREDICTORS, : ) = &
                  tau_predictor_K( 1:MAX_N_STANDARD_PREDICTORS, : ) + &
            solar_tau_predictor_K( 1:MAX_N_STANDARD_PREDICTORS, : )

          solar_tau_predictor_K( 1:MAX_N_STANDARD_PREDICTORS, : ) = ZERO


          ! ----------------------------------------------
          ! K-matrix adjoint of the integrateed predictors
          ! ----------------------------------------------

          CALL compute_predictors_AD( &
                                      ! -- Forward input
                                      layer_p( :, m ),               &  ! Input,  K
                                      layer_t( :, m ),               &  ! Input,  K
                                      layer_w( :, m ),               &  ! Input,  K
                                      solar_tau_absorber( 0:, : ),   &  ! Input,  0:K x J

                                      ! -- K-matrix input
                                      solar_tau_predictor_K( :, : ), &  ! In/Output, I x K

                                      ! -- K-matrix output
                                      layer_p_K( :, l ),             &  ! In/Output,  K
                                      layer_t_K( :, l ),             &  ! In/Output,  K
                                      layer_w_K( :, l ),             &  ! In/Output,  K
                                      solar_tau_absorber_K( 0:, : ), &  ! In/Output,  0:K x J

                                      no_standard = 1                )  ! Optional input


          ! ----------------------------------------------------------
          ! K-matrix adjoint of the nadir absorber amount modification
          ! ----------------------------------------------------------

          absorber_K( 0:, : ) = absorber_K( 0:, : ) + &
                                  ( secant_solar_angle( m ) * solar_tau_absorber_K( 0:, : ) )

        END IF solar_if_block



        !#----------------------------------------------------------------------#
        !# --    Calculate the K-matrix adjoint result for the flux term     -- #
        !#----------------------------------------------------------------------#

        ! --------------------------------------------------
        ! If the current channel is an INFRARED channel,
        ! then calculate the adjoint of the downwelling flux
        ! transmittance using the flux absorber amounts.
        !
        ! If the current channel is a MICROWAVE channel,
        ! then calculate the adjoint of the flux transmittance
        ! using the upwelling absorber amounts.
        ! --------------------------------------------------

        flux_if_block: IF ( is_microwave_channel( channel_index( l ) ) == 0 ) THEN


          ! ------------------------------------
          ! Adjoint of the IR flux transmittance
          ! ------------------------------------

          CALL compute_transmittance_AD( &
                                         ! -- Forward input
                                         flux_tau_absorber( 0:, : ),      &   ! Input, 0:K x J
                                         flux_tau_predictor( :, :, m ),   &   ! Input, I x K
                                         flux_tau( :, l ),                &   ! Input, K

                                         ! -- K-matrix input
                                         flux_tau_K( :, l ),              &   ! In/Output, K

                                         ! -- Other input
                                         flux_tau_layer_index( :, :, m ), &   ! Input, K x J
                                         channel_index( l ),              &   ! Input, scalar
                                         DOWN,                            &   ! Input, scalar

                                         ! -- K-matrix output
                                         flux_tau_absorber_K( 0:, : ),    &   ! In/Output, 0:K x J
                                         flux_tau_predictor_K( :, : )     )   ! In/Output, I x K


          ! -----------------------------------------------
          ! K-matrix adjoint of the standard predictor copy
          ! -----------------------------------------------

          tau_predictor_K( 1:MAX_N_STANDARD_PREDICTORS, : ) = &
                 tau_predictor_K( 1:MAX_N_STANDARD_PREDICTORS, : ) + &
            flux_tau_predictor_K( 1:MAX_N_STANDARD_PREDICTORS, : )

          flux_tau_predictor_K( 1:MAX_N_STANDARD_PREDICTORS, : ) = ZERO


          ! ----------------------------------------------
          ! K-matrix adjoint of the integrateed predictors
          ! ----------------------------------------------

          CALL compute_predictors_AD( &
                                      ! -- Forward input
                                      layer_p( :, m ),              &  ! Input,  K
                                      layer_t( :, m ),              &  ! Input,  K
                                      layer_w( :, m ),              &  ! Input,  K
                                      flux_tau_absorber( 0:, : ),   &  ! Input,  0:K x J

                                      ! -- K-matrix input
                                      flux_tau_predictor_K( :, : ), &  ! In/Output, I x K

                                      ! -- K-matrix output
                                      layer_p_K( :, l ),            &  ! In/Output,  K
                                      layer_t_K( :, l ),            &  ! In/Output,  K
                                      layer_w_K( :, l ),            &  ! In/Output,  K
                                      flux_tau_absorber_K( 0:, : ), &  ! In/Output,  0:K x J

                                      no_standard = 1               )  ! Optional input


          ! ----------------------------------------------------------
          ! K-matrix adjoint of the nadir absorber amount modification
          ! ----------------------------------------------------------

          absorber_K( 0:, : ) = absorber_K( 0:, : ) + &
                                   ( SECANT_DIFFUSIVITY_ANGLE * flux_tau_absorber_K( 0:, : ) )


        ELSE  ! We have a microwave channel....

          ! ----------------------------------------------
          ! If total transmittance /= 0, calculate adjoint
          ! flux transmittance for the microwave
          ! ----------------------------------------------
        
          IF ( tau( n_layers, l ) > TOLERANCE ) THEN
            DO k = n_layers, 2, -1
              flux_tau_K( 1, l ) = flux_tau_K( 1, l ) + ( flux_tau_K( k, l ) / &
              !                                           ------------------
                                                             tau( k-1, l )   )

              tau_K( k-1, l ) = tau_K( k-1, l ) - ( flux_tau_K( k, l ) * flux_tau( 1, l ) / &
              !                                     -------------------------------------
                                                                tau( k-1, l )**2          )
              flux_tau_K( k, l ) = ZERO
            END DO
            tau_K( n_layers, l ) = tau_K( n_layers, l ) + flux_tau_K( 1, l )
            flux_tau_K( 1, l ) = ZERO
          ELSE
            flux_tau_K( :, l ) = ZERO
          END IF


        END IF flux_if_block



        !#----------------------------------------------------------------------#
        !# --   Calculate the K-matrix adjoint result for the transmittance  -- #
        !#----------------------------------------------------------------------#

        ! -----------------------------------------------------
        ! Calculate the adjoint of the upwelling transmittances
        ! for the satellite view angle
        ! -----------------------------------------------------

        CALL compute_transmittance_AD( &
                                       ! -- Forward input
                                       tau_absorber( 0:, : ),      &   ! Input, 0:K x J
                                       tau_predictor( :, :, m ),   &   ! Input, I x K
                                       tau( :, l ),                &   ! Input, K

                                       ! -- K-matrix input
                                       tau_K( :, l ),              &   ! In/Output, K

                                       ! -- Other input
                                       tau_layer_index( :, :, m ), &   ! Input, K x J
                                       channel_index( l ),         &   ! Input, scalar
                                       UP,                         &   ! Input, scalar

                                       ! -- K-matrix output
                                       tau_absorber_K( 0:, : ),    &   ! In/Output, 0:K x J
                                       tau_predictor_K( :, : )     )   ! In/Output, I x K



        ! --------------------------------------
        ! K-matrix adjoint of all the predictors
        ! --------------------------------------

        CALL compute_predictors_AD( &
                                    ! -- Forward input
                                    layer_p( :, m ),         &  ! Input,  K
                                    layer_t( :, m ),         &  ! Input,  K
                                    layer_w( :, m ),         &  ! Input,  K
                                    tau_absorber( 0:, : ),   &  ! Input,  0:K x J

                                    ! -- K-matrix input
                                    tau_predictor_K( :, : ), &  ! In/Output, I x K

                                    ! -- K-matrix output
                                    layer_p_K( :, l ),       &  ! In/Output,  K
                                    layer_t_K( :, l ),       &  ! In/Output,  K
                                    layer_w_K( :, l ),       &  ! In/Output,  K
                                    tau_absorber_K( 0:, : )  )  ! In/Output,  0:K x J


        ! ----------------------------------------------------------
        ! K-matrix adjoint of the nadir absorber amount modification
        ! ----------------------------------------------------------

        absorber_K( 0:, : ) = absorber_K( 0:, : ) + &
                                 ( secant_view_angle( m ) * tau_absorber_K( 0:, : ) )



        !#----------------------------------------------------------------------#
        !#            -- Calculate the absorber K-matrix adjoints --            #
        !#----------------------------------------------------------------------#

        CALL compute_absorber_amount_AD( &
                                         ! -- Forward input
                                         level_p( :, m ),     &  ! Input,  K
                                         layer_w( :, m ),     &  ! Input,  K
                                         layer_o( :, m ),     &  ! Input,  K

                                         ! -- K-matrix input
                                         absorber_K( 0:, : ), &  ! In/Output, 0:K x J

                                         ! -- K-matrix output
                                         level_p_K( :, l ),   &  ! In/Ouput,  K
                                         layer_w_K( :, l ),   &  ! In/Ouput,  K
                                         layer_o_K( :, l )    )  ! In/Ouput,  K


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

  END FUNCTION k_matrix_rtm

END MODULE k_matrix_model


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
! Revision 1.4  2001/09/04 21:29:11  paulv
! - Updated documentation.
!
! Revision 1.3  2001/08/31 21:28:52  paulv
! - Removed input data checks from COMPUTE_RTM_K. The same checks are performed
!   in the main routine, K_MATRIX_RTM, so there was no need to replicate them.
! - Added check for negative profile and surface data in K_MATRIX_RTM.
! - Maximum solar angle secant is no longer calculated in K_MATRIX_RTM but
!   is declared as a parameter in the PARAMETERS module.
!
! Revision 1.2  2001/08/16 16:38:57  paulv
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
! Revision 1.1  2001/08/01 16:34:02  paulv
! Initial checkin.
!
!
!
!
