!------------------------------------------------------------------------------
!M+
! NAME:
!       tangent_linear_model
!
! PURPOSE:
!       Module containing the NCEP RT tangent linear model function.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE tangent_linear_model
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
!       tangent_linear_rtm:    PUBLIC function that calculates top-of-atmosphere (TOA)
!                              tangent-linear radiances and brightness temperatures for
!                              user specified profiles and satellites/channels. 
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

MODULE tangent_linear_model


  ! ------------
  ! Module usage
  ! ------------

  USE type_kinds, ONLY : fp_kind
  USE error_handler
  USE parameters
  USE spectral_coefficients
  USE absorber_profile
  USE predictors
  USE transmittance
  USE radiance


  ! -----------------------
  ! Disable implicit typing
  ! -----------------------

  IMPLICIT NONE


  ! --------------------
  ! Default visibilities
  ! --------------------

  PRIVATE
  PUBLIC :: tangent_linear_rtm


  ! --------------------
  ! Function overloading
  ! --------------------

  INTERFACE tangent_linear_rtm
    MODULE PROCEDURE tangent_linear_rtm_rank1
    MODULE PROCEDURE tangent_linear_rtm_rank2
  END INTERFACE ! tangent_linear_rtm


CONTAINS


!--------------------------------------------------------------------------------
!S+
! NAME:
!       tangent_linear_rtm
!
! PURPOSE:
!       PUBLIC function that calculates top-of-atmosphere (TOA) tangent-linear
!       radiances and brightness temperatures for an input atmospheric profile
!       set and user specified satellites/channels.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!
!       result = tangent_linear_rtm( &
!                                    ! Forward inputs
!                                    level_p, layer_p, layer_t, layer_w, layer_o,                &  ! Input, K x M
!
!                                    surface_temperature,                                        &  ! Input, M
!                                    surface_emissivity,                                         &  ! Input, L*M
!                                    surface_reflectivity,                                       &  ! Input, L*M
!
!                                    absorber,                                                   &  ! Input, 0:K x J x M
!
!                                    tau_layer_index,                                            &  ! Input, K x J x M
!                                    flux_tau_layer_index,                                       &  ! Input, K x J x M
!                                    solar_tau_layer_index,                                      &  ! Input, K x J x M
!
!                                    tau_predictor,                                              &  ! Input, Imax x K x M
!                                    flux_tau_predictor,                                         &  ! Input, Imax x K x M
!                                    solar_tau_predictor,                                        &  ! Input, Imax x K x M
!
!                                    tau,                                                        &  ! Input, K x L*M
!                                    flux_tau,                                                   &  ! Input, K x L*M
!                                    solar_tau,                                                  &  ! Input, K x L*M
!
!                                    layer_radiance,                                             &  ! Input, K x L*M
!                                    downwelling_radiance,                                       &  ! Input, L*M
!                                    upwelling_radiance,                                         &  ! Input, L*M
!
!                                    ! -- Tangent-linear inputs
!                                    level_p_TL, layer_p_TL, layer_t_TL, layer_w_TL, layer_o_TL, &  ! Input, K x M
!
!                                    surface_temperature_TL,                                     &  ! Input, M
!                                    surface_emissivity_TL,                                      &  ! Input, L*M
!                                    surface_reflectivity_TL,                                    &  ! Input, L*M
!
!                                    ! -- Other inputs
!                                    secant_view_angle,                                          &  ! Input, M
!                                    secant_solar_angle,                                         &  ! Input, M
!                                    n_channels_per_profile,                                     &  ! Input, M
!                                    channel_index,                                              &  ! Input, L*M
!
!                                    ! -- Tangent-linear outputs
!                                    tau_TL,                                                     &  ! Output, K x L*M
!                                    flux_tau_TL,                                                &  ! Output, K x L*M
!                                    solar_tau_TL,                                               &  ! Output, K x L*M
!
!                                    layer_radiance_TL,                                          &  ! Output, K x L*M
!                                    downwelling_radiance_TL,                                    &  ! Output, L*M
!                                    upwelling_radiance_TL,                                      &  ! Output, L*M
!
!                                    brightness_temperature_TL,                                  &  ! Output, L*M
!
!                                    ! Optional inputs
!                                    message_log = message_log )
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
!                                  DIMENSION:  I x K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       flux_tau_predictor:        Predictor profiles for the thermal flux transmittance.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  I x K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       solar_tau_predictor:       Predictor profiles for the solar transmittance.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  I x K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       tau:                       Layer->TOA transmittance for the satellite
!                                  view angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       flux_tau:                  Layer->SFC transmittance for the default
!                                  diffusivity angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       solar_tau:                 Layer->SFC transmittance for the solar
!                                  zenith angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_radiance:            Layer Planck radiances at every layer for
!                                  each channel/profile.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       downwelling_radiance:      TOA->SFC radiances for each channel/profile due
!                                  to thermal flux and solar components.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!       upwelling_radiance:        TOA radiances for each channel/profile.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!       level_p_TL:                Profile set layer interface pressure tangent-linear
!                                  array.
!                                  UNITS:      hPa
!                                  TYPE:       Real
!                                  DIMENSION:  K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_p_TL:                Profile set layer average pressure tangent-linear
!                                  array.
!                                  UNITS:      hPa
!                                  TYPE:       Real
!                                  DIMENSION:  K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_t_TL:                Profile set layer average temperature tangent-linear
!                                  array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_w_TL:      .         Profile set layer average water vapor mixing ratio
!                                  tangent-linear array.
!                                  UNITS:      g/kg
!                                  TYPE:       Real
!                                  DIMENSION:  K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_o_TL:                Profile set layer average ozone mixing ratio
!                                  tangent-linear array.
!                                  UNITS:      ppmv
!                                  TYPE:       Real
!                                  DIMENSION:  K x M; K > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_temperature_TL:    Profile set surface temperature tangent-linear array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  M; M > or = 1 (i.e. scalar)
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_emissivity_TL:     Profile set surface emissivity tangent-linear array
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_reflectivity_TL:   Profile set surface reflectivity tangent-linear array
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
!       tau_TL:                    Layer->TOA tangent-linear transmittance for the satellite
!                                  view angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( OUT )
!
!       flux_tau_TL:               Layer->SFC tangent-linear transmittance for the default
!                                  diffusivity angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( OUT )
!
!       solar_tau_TL:              Layer->SFC tangent-linear transmittance for the solar
!                                  zenith angle.
!                                  UNITS:      None.
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( OUT )
!
!       layer_radiance_TL:         Layer Planck tangent-linear radiances at every layer for
!                                  each channel/profile.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  K x L*M; K > 1, L > 1, and M > or = 1
!                                  ATTRIBUTES: INTENT( OUT )
!
!       downwelling_radiance_TL:   TOA->SFC tangent-linear radiances for each channel/profile due
!                                  to thermal flux and solar components.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( OUT )
!
!       upwelling_radiance_TL:     TOA tangent-linear radiances for each channel/profile.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( OUT )
!
!       brightness_temperature_TL: Tangent-linear temperatures corresponding to the
!                                  TOA tangent-linear radiances.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  L*M; L > 1, and M > or = 1
!                                              NB: This is a 1-D array.
!                                  ATTRIBUTES: INTENT( OUT )
!
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
!      compute_absorber_amount_TL: Subroutine to integrate the tangent-linear
!                                  absorber profiles
!                                  SOURCE: absorber_profile module
!
!      compute_predictors_TL:      Subroutine to compute the tangent-linear 
!                                  transmittance predictor profiles.
!                                  SOURCE: predictor module
!
!      compute_transmittance_TL:   Subroutine to compute the tangent-linear
!                                  transmittance profiles.
!                                  SOURCE: transmittance module
!
!      compute_radiance_TL:        Subroutine to compute the TOA tangent-linear 
!                                  radiances and brightness temperatures.
!                                  SOURCE: radiance module
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
!       Note the restrictions on the input array dimensions:
!         K == n_layers > 1
!         L == n_channels > 1
!         M == n_profiles > or = 1
!
! PROCEDURE:
!       See individual module function documentation.
!S-
!--------------------------------------------------------------------------------

  FUNCTION tangent_linear_rtm_rank2( &

             ! -- Forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o,                &  ! Input, K x M

             surface_temperature,                                        &  ! Input, M
             surface_emissivity,                                         &  ! Input, L*M
             surface_reflectivity,                                       &  ! Input, L*M

             absorber,                                                   &  ! Input, 0:K x J x M

             tau_layer_index,                                            &  ! Input, K x J x M
             flux_tau_layer_index,                                       &  ! Input, K x J x M
             solar_tau_layer_index,                                      &  ! Input, K x J x M

             tau_predictor,                                              &  ! Input, Imax x K x M
             flux_tau_predictor,                                         &  ! Input, Imax x K x M
             solar_tau_predictor,                                        &  ! Input, Imax x K x M

             tau,                                                        &  ! Input, K x L*M
             flux_tau,                                                   &  ! Input, K x L*M
             solar_tau,                                                  &  ! Input, K x L*M

             layer_radiance,                                             &  ! Input, K x L*M
             downwelling_radiance,                                       &  ! Input, L*M
             upwelling_radiance,                                         &  ! Input, L*M

             ! -- Tangent-linear inputs
             level_p_TL, layer_p_TL, layer_t_TL, layer_w_TL, layer_o_TL, &  ! Input, K x M

             surface_temperature_TL,                                     &  ! Input, M
             surface_emissivity_TL,                                      &  ! Input, L*M
             surface_reflectivity_TL,                                    &  ! Input, L*M

             ! -- Other inputs
             secant_view_angle,                                          &  ! Input, M
             secant_solar_angle,                                         &  ! Input, M
             n_channels_per_profile,                                     &  ! Input, M
             channel_index,                                              &  ! Input, L*M

             ! -- Tangent-linear outputs
             tau_TL,                                                     &  ! Output, K x L*M
             flux_tau_TL,                                                &  ! Output, K x L*M
             solar_tau_TL,                                               &  ! Output, K x L*M

             layer_radiance_TL,                                          &  ! Output, K x L*M
             downwelling_radiance_TL,                                    &  ! Output, L*M
             upwelling_radiance_TL,                                      &  ! Output, L*M

             brightness_temperature_TL,                                  &  ! Output, L*M

             ! -- Optional inputs
             n_input_profiles,                                           &
             message_log )                                               &

           RESULT ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward inputs
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: level_p                    ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_p                    ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_t                    ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_w                    ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_o                    ! K x M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: surface_temperature        ! M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: surface_emissivity         ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: surface_reflectivity       ! L*M

    REAL( fp_kind ), DIMENSION( 0:, :, : ), INTENT( IN )  :: absorber                   ! 0:K x J x M

    INTEGER,         DIMENSION(  :, :, : ), INTENT( IN )  :: tau_layer_index            ! K x J x M
    INTEGER,         DIMENSION(  :, :, : ), INTENT( IN )  :: flux_tau_layer_index       ! K x J x M
    INTEGER,         DIMENSION(  :, :, : ), INTENT( IN )  :: solar_tau_layer_index      ! K x J x M

    REAL( fp_kind ), DIMENSION( :, :, : ),  INTENT( IN )  :: tau_predictor              ! Imax x K x M
    REAL( fp_kind ), DIMENSION( :, :, : ),  INTENT( IN )  :: flux_tau_predictor         ! Imax x K x M
    REAL( fp_kind ), DIMENSION( :, :, : ),  INTENT( IN )  :: solar_tau_predictor        ! Imax x K x M

    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: tau                        ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: flux_tau                   ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: solar_tau                  ! K x L*M

    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_radiance             ! K x L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: downwelling_radiance       ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: upwelling_radiance         ! L*M

    ! -- Tangent-linear inputs
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: level_p_TL                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_p_TL                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_t_TL                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_w_TL                 ! K x M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( IN )  :: layer_o_TL                 ! K x M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: surface_temperature_TL     ! M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: surface_emissivity_TL      ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: surface_reflectivity_TL    ! L*M

    ! -- Other inputs
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: secant_view_angle          ! M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( IN )  :: secant_solar_angle         ! M
    INTEGER,         DIMENSION( : ),        INTENT( IN )  :: n_channels_per_profile     ! M
    INTEGER,         DIMENSION( : ),        INTENT( IN )  :: channel_index              ! L*M

    ! -- Tangent-linear outputs
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT ) :: tau_TL                     ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT ) :: flux_tau_TL                ! K x L*M
    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT ) :: solar_tau_TL               ! K x L*M

    REAL( fp_kind ), DIMENSION( :, : ),     INTENT( OUT ) :: layer_radiance_TL          ! K x L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( OUT ) :: downwelling_radiance_TL    ! L*M
    REAL( fp_kind ), DIMENSION( : ),        INTENT( OUT ) :: upwelling_radiance_TL      ! L*M

    REAL( fp_kind ), DIMENSION( : ),        INTENT( OUT ) :: brightness_temperature_TL  ! L*M

    ! -- Optional input
    INTEGER,        OPTIONAL,               INTENT( IN )  :: n_input_profiles          ! Scalar
    CHARACTER( * ), OPTIONAL,               INTENT( IN )  :: message_log
    
    
    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'TANGENT_LINEAR_RTM_RANK2'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Scalars
    CHARACTER( 100 ) :: message
    CHARACTER( 5 )   :: value_in, value_allowed

    INTEGER :: m, n_profiles          ! Profile loop/index variables
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

      error_status = tangent_linear_rtm_rank1( &
                       ! -- Forward inputs
                       level_p( :, m ), layer_p( :, m ), layer_t( :, m ), layer_w( :, m ), layer_o( :, m ), &  ! Input, K

                       surface_temperature( m ),           &  ! Input,  Scalar
                       surface_emissivity( l1:l2 ),        &  ! Input,  L
                       surface_reflectivity( l1:l2 ),      &  ! Input,  L

                       absorber( 0:, :, m ),               &  ! Input, 0:K x J

                       tau_layer_index( :, :, m ),         &  ! Input, K x J
                       flux_tau_layer_index( :, :, m ),    &  ! Input, K x J
                       solar_tau_layer_index( :, :, m ),   &  ! Input, K x J

                       tau_predictor( :, :, m ),           &  ! Input, Imax x K
                       flux_tau_predictor( :, :, m ),      &  ! Input, Imax x K
                       solar_tau_predictor( :, :, m ),     &  ! Input, Imax x K

                       tau( :, l1:l2 ),                    &  ! Input, K x L
                       flux_tau( :, l1:l2 ),               &  ! Input, K x L
                       solar_tau( :, l1:l2 ),              &  ! Input, K x L

                       layer_radiance( :, l1:l2 ),         &  ! Input, K x L
                       downwelling_radiance( l1:l2 ),      &  ! Input, L
                       upwelling_radiance( l1:l2 ),        &  ! Input, L

                       ! -- Tangent-linear inputs
                       level_p_TL( :, m ), layer_p_TL( :, m ), layer_t_TL( :, m ), layer_w_TL( :, m ), layer_o_TL( :, m ), &  ! Input, K

                       surface_temperature_TL,             &  ! Input, Scalar
                       surface_emissivity_TL( l1:l2 ),     &  ! Input, L
                       surface_reflectivity_TL( l1:l2 ),   &  ! Input, L

                       ! -- Other inputs
                       secant_view_angle( m ),             &  ! Input,  Scalar
                       secant_solar_angle( m ),            &  ! Input,  Scalar
                       n_channels_per_profile( m ),        &  ! Input,  Scalar
                       channel_index( l1:l2 ),             &  ! Input,  L

                       ! -- Tangent-linear outputs
                       tau_TL( :, l1:l2 ),                 &  ! Output, K x L*M
                       flux_tau_TL( :, l1:l2 ),            &  ! Output, K x L*M
                       solar_tau_TL( :, l1:l2 ),           &  ! Output, K x L*M

                       layer_radiance_TL( :, l1:l2 ),      &  ! Output, K x L*M
                       downwelling_radiance_TL( l1:l2 ),   &  ! Output, L*M
                       upwelling_radiance_TL( l1:l2 ),     &  ! Output, L*M

                       brightness_temperature_TL( l1:l2 ), &  ! Output, L*M

                       ! -- Optional inputs
                       message_log = message_log )


      ! -------------------------------
      ! Check for successful completion
      ! -------------------------------

      IF ( error_status /= SUCCESS ) THEN

        CALL display_message( ROUTINE_NAME, &
                              'Error occured in TANGENT_LINEAR_RTM_RANK1', &
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

  END FUNCTION tangent_linear_rtm_rank2



  FUNCTION tangent_linear_rtm_rank1( &

             ! Forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o,                &  ! Input, K

             surface_temperature,                                        &  ! Input, Scalar
             surface_emissivity,                                         &  ! Input, L
             surface_reflectivity,                                       &  ! Input, L

             absorber,                                                   &  ! Input, 0:K x J

             tau_layer_index,                                            &  ! Input, K x J
             flux_tau_layer_index,                                       &  ! Input, K x J
             solar_tau_layer_index,                                      &  ! Input, K x J

             tau_predictor,                                              &  ! Input, Imax x K
             flux_tau_predictor,                                         &  ! Input, Imax x K
             solar_tau_predictor,                                        &  ! Input, Imax x K

             tau,                                                        &  ! Input, K x L
             flux_tau,                                                   &  ! Input, K x L
             solar_tau,                                                  &  ! Input, K x L

             layer_radiance,                                             &  ! Input, K x L
             downwelling_radiance,                                       &  ! Input, L
             upwelling_radiance,                                         &  ! Input, L

             ! -- Tangent-linear inputs
             level_p_TL, layer_p_TL, layer_t_TL, layer_w_TL, layer_o_TL, &  ! Input, K

             surface_temperature_TL,                                     &  ! Input, Scalar
             surface_emissivity_TL,                                      &  ! Input, L
             surface_reflectivity_TL,                                    &  ! Input, L

             ! -- Other inputs
             secant_view_angle,                                          &  ! Input, Scalar
             secant_solar_angle,                                         &  ! Input, Scalar
             n_channels,                                                 &  ! Input, Scalar
             channel_index,                                              &  ! Input, L

             ! -- Tangent-linear outputs
             tau_TL,                                                     &  ! Output, K x L
             flux_tau_TL,                                                &  ! Output, K x L
             solar_tau_TL,                                               &  ! Output, K x L

             layer_radiance_TL,                                          &  ! Output, K x L
             downwelling_radiance_TL,                                    &  ! Output, L
             upwelling_radiance_TL,                                      &  ! Output, L

             brightness_temperature_TL,                                  &  ! Output, L

             ! -- Optional inputs
             n_input_profiles,                                           &
             message_log )                                               &

           RESULT ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward inputs
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: level_p                    ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_p                    ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_t                    ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_w                    ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_o                    ! K

    REAL( fp_kind ),                     INTENT( IN )  :: surface_temperature        ! Scalar
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: surface_emissivity         ! L
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: surface_reflectivity       ! L

    REAL( fp_kind ), DIMENSION( 0:, : ), INTENT( IN )  :: absorber                   ! 0:K x J

    INTEGER,         DIMENSION(  :, : ), INTENT( IN )  :: tau_layer_index            ! K x J
    INTEGER,         DIMENSION(  :, : ), INTENT( IN )  :: flux_tau_layer_index       ! K x J
    INTEGER,         DIMENSION(  :, : ), INTENT( IN )  :: solar_tau_layer_index      ! K x J

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )  :: tau_predictor              ! Imax x K
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )  :: flux_tau_predictor         ! Imax x K
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )  :: solar_tau_predictor        ! Imax x K

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )  :: tau                        ! K x L
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )  :: flux_tau                   ! K x L
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )  :: solar_tau                  ! K x L

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( IN )  :: layer_radiance             ! K x L
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: downwelling_radiance       ! L
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: upwelling_radiance         ! L

    ! -- Tangent-linear inputs
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: level_p_TL                 ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_p_TL                 ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_t_TL                 ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_w_TL                 ! K
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: layer_o_TL                 ! K

    REAL( fp_kind ),                     INTENT( IN )  :: surface_temperature_TL     ! Scalar
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: surface_emissivity_TL      ! L
    REAL( fp_kind ), DIMENSION( : ),     INTENT( IN )  :: surface_reflectivity_TL    ! L

    ! -- Other inputs
    REAL( fp_kind ),                     INTENT( IN )  :: secant_view_angle          ! Scalar
    REAL( fp_kind ),                     INTENT( IN )  :: secant_solar_angle         ! Scalar
    INTEGER,                             INTENT( IN )  :: n_channels                 ! Scalar
    INTEGER,         DIMENSION( : ),     INTENT( IN )  :: channel_index              ! L

    ! -- Tangent-linear outputs
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: tau_TL                     ! K x L
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: flux_tau_TL                ! K x L
    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: solar_tau_TL               ! K x L

    REAL( fp_kind ), DIMENSION( :, : ),  INTENT( OUT ) :: layer_radiance_TL          ! K x L
    REAL( fp_kind ), DIMENSION( : ),     INTENT( OUT ) :: downwelling_radiance_TL    ! L
    REAL( fp_kind ), DIMENSION( : ),     INTENT( OUT ) :: upwelling_radiance_TL      ! L

    REAL( fp_kind ), DIMENSION( : ),     INTENT( OUT ) :: brightness_temperature_TL  ! L

    ! -- Optional input. Note that N_INPUT_PROFILES is not used in this
    ! -- function. It is specified so if a user specifies it by mistake
    ! -- for rank-1 profile input the code won't (hopefully) fall in a heap.
    INTEGER,        OPTIONAL,            INTENT( IN )  :: n_input_profiles           ! Scalar
    CHARACTER( * ), OPTIONAL,            INTENT( IN )  :: message_log
    
    
    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'TANGENT_LINEAR_RTM_RANK1'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Scalars
    CHARACTER( 100 ) :: message
    CHARACTER( 5 )   :: value_in, value_allowed

    INTEGER :: k, n_layers  ! Layer loop/index variables
    INTEGER :: l            ! Channel loop/index variable

    INTEGER :: valid_solar

    ! -- Maximum channels pseudo parameter
    INTEGER :: MAX_N_CHANNELS
    LOGICAL :: is_set

    ! -- Arrays for integrated absorber amounts, 0:K x J
    REAL( fp_kind ), DIMENSION( 0:SIZE( absorber, DIM = 1 )-1, &
                                  SIZE( absorber, DIM = 2 )    ) :: tau_absorber,      &
                                                                    flux_tau_absorber, &
                                                                    solar_tau_absorber

    REAL( fp_kind ), DIMENSION( 0:SIZE( absorber, DIM = 1 )-1, &
                                  SIZE( absorber, DIM = 2 )    ) :: absorber_TL, &
                                                                    tau_absorber_TL,      &
                                                                    flux_tau_absorber_TL, &
                                                                    solar_tau_absorber_TL

    ! -- Arrays for tangent-linear predictors, Imax x K
    REAL( fp_kind ), DIMENSION( SIZE( tau_predictor, DIM = 1 ), &
                                SIZE( tau_predictor, DIM = 2 )  ) :: tau_predictor_TL,      &
                                                                     flux_tau_predictor_TL, &
                                                                     solar_tau_predictor_TL

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
    !#             -- DETERMINE ARRAY DIMENSIONS AND CHECK INPUT --             #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------
    ! Check the number of channels - if zero
    ! then simply RETURN.
    ! --------------------------------------

    IF ( n_channels == 0 ) RETURN


    ! ------------------
    ! Get the dimensions
    ! ------------------

    n_layers   = SIZE( layer_p )


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
    !#   -- CALCULATE THE PROFILE GENERIC TANGENT-LINEAR ABSORBER AMOUNTS --    #
    !#--------------------------------------------------------------------------#

    CALL compute_absorber_amount_TL( &
                                     ! -- Forward input
                                     level_p( : ),        &  ! Input,  K
                                     layer_w( : ),        &  ! Input,  K
                                     layer_o( : ),        &  ! Input,  K

                                     ! -- Tangent-linear input
                                     level_p_TL( : ),     &  ! Input,  K
                                     layer_w_TL( : ),     &  ! Input,  K
                                     layer_o_TL( : ),     &  ! Input,  K

                                     ! -- Tangent-linear output
                                     absorber_TL( 0:, : ) )  ! Output, 0:K x J



    !#--------------------------------------------------------------------------#
    !#      -- CALCULATE THE PREDICTORS FOR THE UPWELLING TRANSMITTANCE --      #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------------------
    ! Modify absorber quantities by the angle secant
    ! Could put a loop here but here's hoping the compiler
    ! recognises this as a group of loops over layer.
    ! ----------------------------------------------------

    tau_absorber(    0:, : ) = secant_view_angle * absorber(    0:, : )
    tau_absorber_TL( 0:, : ) = secant_view_angle * absorber_TL( 0:, : )

    ! -- Subtract the top pressure from the dry absorber profile
    tau_absorber( 1:, 2 ) = tau_absorber( 1:, 2 ) - TOA_PRESSURE


    ! ---------------------------------------
    ! Calculate the tangent-linear predictors
    ! for the satellite view angle
    ! ---------------------------------------

    CALL compute_predictors_TL( &
                                ! -- Forward input
                                layer_p( : ),              &  ! Input,  K
                                layer_t( : ),              &  ! Input,  K
                                layer_w( : ),              &  ! Input,  K
                                tau_absorber( 0:, : ),     &  ! Input,  0:K x J

                                ! -- Tangent-linear input
                                layer_p_TL( : ),           &  ! Input,  K
                                layer_t_TL( : ),           &  ! Input,  K
                                layer_w_TL( : ),           &  ! Input,  K
                                tau_absorber_TL( 0:, : ),  &  ! Input,  0:K x J

                                ! -- Tangent-linear output
                                tau_predictor_TL( :, : )   )  ! Output, I x K
 


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

      flux_tau_absorber(    0:, : ) = SECANT_DIFFUSIVITY_ANGLE * absorber(    0:, : )
      flux_tau_absorber_TL( 0:, : ) = SECANT_DIFFUSIVITY_ANGLE * absorber_TL( 0:, : )
      
      ! -- Subtract the top pressure from the dry absorber profile
      flux_tau_absorber( 1:, 2 ) = flux_tau_absorber( 1:, 2 ) - TOA_PRESSURE


      ! ---------------------------------------
      ! Calculate the tangent-linear predictors
      ! for the diffusivity angle
      ! ---------------------------------------

      ! -- Calculate the integrated predictors only
      CALL compute_predictors_TL( &
                                  ! -- Forward input
                                  layer_p( : ),                  &  ! Input,  K
                                  layer_t( : ),                  &  ! Input,  K
                                  layer_w( : ),                  &  ! Input,  K
                                  flux_tau_absorber( 0:, : ),    &  ! Input,  0:K x J

                                  ! -- Tangent-linear input
                                  layer_p_TL( : ),               &  ! Input,  K
                                  layer_t_TL( : ),               &  ! Input,  K
                                  layer_w_TL( : ),               &  ! Input,  K
                                  flux_tau_absorber_TL( 0:, : ), &  ! Input,  0:K x J

                                  ! -- Tangent-linear output
                                  flux_tau_predictor_TL( :, : ), &  ! Output, I x K

                                  no_standard = 1                )  ! Optional input

      ! -- Copy the angle independent (standard) predictors
      flux_tau_predictor_TL( 1:MAX_N_STANDARD_PREDICTORS, : ) = &
           tau_predictor_TL( 1:MAX_N_STANDARD_PREDICTORS, : ) 

    END IF



    !#--------------------------------------------------------------------------#
    !#        -- CALCULATE THE PREDICTORS FOR THE SOLAR TRANSMITTANCE --        #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------------------
    ! Have *any* SOLAR sensitive channels been specified
    ! for the current profile (Flagged as == 1)?
    !
    ! AND
    !
    ! Is the specified solar zenith angle valid (<85)?
    ! --------------------------------------------------

    IF ( ( ANY( is_solar_channel( channel_index( 1:n_channels ) ) == 1 ) ) .AND. &
         secant_solar_angle < MAX_SECANT_SOLAR_ANGLE ) THEN


      ! --------------------------------
      ! Modify the nadir absorber amount
      ! --------------------------------

      solar_tau_absorber(    0:, : ) = secant_solar_angle * absorber(    0:, : )
      solar_tau_absorber_TL( 0:, : ) = secant_solar_angle * absorber_TL( 0:, : )

      ! -- Subtract the top pressure from the dry absorber profile
      solar_tau_absorber( 1:, 2 ) = solar_tau_absorber( 1:, 2 ) - TOA_PRESSURE


      ! ---------------------------------------
      ! Calculate the tangent-linear predictors
      ! for the solar zenith angle
      ! ---------------------------------------

      ! -- Calculate the integrated predictors only
      CALL compute_predictors_TL( &
                                  ! -- Forward input
                                  layer_p( : ),                   &  ! Input,  K
                                  layer_t( : ),                   &  ! Input,  K
                                  layer_w( : ),                   &  ! Input,  K
                                  solar_tau_absorber( 0:, : ),    &  ! Input,  0:K x J

                                  ! -- Tangent-linear input
                                  layer_p_TL( : ),                &  ! Input,  K
                                  layer_t_TL( : ),                &  ! Input,  K
                                  layer_w_TL( : ),                &  ! Input,  K
                                  solar_tau_absorber_TL( 0:, : ), &  ! Input,  0:K x J

                                  ! -- Tangent-linear output
                                  solar_tau_predictor_TL( :, : ), &  ! Output, I x K

                                  no_standard = 1                 )  ! Optional input

      ! -- Copy the angle independent predictors
      solar_tau_predictor_TL( 1:MAX_N_STANDARD_PREDICTORS, : ) = &
            tau_predictor_TL( 1:MAX_N_STANDARD_PREDICTORS, : ) 

    END IF
     

    !#--------------------------------------------------------------------------#
    !#                            -- CHANNEL LOOP --                            #
    !#--------------------------------------------------------------------------#

    l_channel_loop: DO l = 1, n_channels


      ! --------------------------------------------
      ! Calculate the current channel tangent-linear
      ! transmittances for the satellite view angle
      ! --------------------------------------------

      CALL compute_transmittance_TL( &
                                     ! -- Forward input
                                     tau_absorber( 0:, : ),    &   ! Input, 0:K x J
                                     tau_predictor( :, : ),    &   ! Input, I x K
                                     tau( :, l ),              &   ! Input, K

                                     ! -- Tangent-linear input
                                     tau_absorber_TL( 0:, : ), &   ! Input, 0:K x J
                                     tau_predictor_TL( :, : ), &   ! Input, I x K

                                     ! -- Other input
                                     tau_layer_index( :, : ),  &   ! Input, K x J
                                     channel_index( l ),       &   ! Input, scalar
                                     UP,                       &   ! Input, scalar

                                     ! -- Tangent-linear output
                                     tau_TL( :, l )            )   ! Output, K



      ! --------------------------------------------------
      ! If the current channel is an INFRARED channel,
      ! then calculate the downwelling flux transmittance.
      !
      ! If the current channel is a MICROWAVE channel,
      ! then use the predictors for the upwelling
      ! transmittance calculations.
      !
      ! This situation is a toss-up. One could return the
      ! layer optical depth rather than the transmittance
      ! and use that to derive the microwave flux
      ! transmittance but then there would be two sets
      ! of quantities to be carried through - the layer
      ! to boundary transmittances and the individual
      ! layer optical depths. If the speed penalty of
      ! calculating the microwave flux transmittances in
      ! this manner is too high, and memory is available,
      ! then this is an option.
      ! --------------------------------------------------

      IF ( is_microwave_channel( channel_index( l ) ) == 0 ) THEN

        ! -- IR channel
        CALL compute_transmittance_TL( &
                                       ! -- Forward input
                                       flux_tau_absorber( 0:, : ),    &   ! Input, 0:K x J
                                       flux_tau_predictor( :, : ),    &   ! Input, I x K
                                       flux_tau( :, l ),              &   ! Input, K

                                       ! -- Tangent-linear input
                                       flux_tau_absorber_TL( 0:, : ), &   ! Input, 0:K x J
                                       flux_tau_predictor_TL( :, : ), &   ! Input, I x K

                                       ! -- Other input
                                       flux_tau_layer_index( :, : ),  &   ! Input, K x J
                                       channel_index( l ),            &   ! Input, scalar
                                       DOWN,                          &   ! Input, scalar

                                       ! -- Tangent-linear output
                                       flux_tau_TL( :, l )            )   ! Output, K

      ELSE

        ! -- uW channel
!        CALL compute_transmittance_TL( &
!                                       ! -- Forward input
!                                       tau_absorber( 0:, : ),    &   ! Input, 0:K x J
!                                       tau_predictor( :, : ),    &   ! Input, I x K
!                                       flux_tau( :, l ),         &   ! Input, K
!
!                                       ! -- Tangent-linear input
!                                       tau_absorber_TL( 0:, : ), &   ! Input, 0:K x J
!                                       tau_predictor_TL( :, : ), &   ! Input, I x K
!
!                                       ! -- Other input
!                                       tau_layer_index( :, : ),  &   ! Input, K x J
!                                       channel_index( l ),       &   ! Input, scalar
!                                       DOWN,                     &   ! Input, scalar
!
!                                       ! -- Tangent-linear output
!                                       flux_tau_TL( :, l )       )   ! Output, K

        ! -- If total transmittance /= 0, calculate 
        ! -- tangent-linear flux transmittance
        flux_tau_TL( :, l ) = ZERO
        IF ( tau( n_layers, l ) > TOLERANCE ) THEN
          flux_tau_TL( 1, l ) = tau_TL( n_layers, l )
          DO k = 2, n_layers

            flux_tau_TL( k, l ) = ( flux_tau_TL( 1, l ) / &
            !                       -------------------   
                                      tau( k-1, l )     ) - &

                                  ( flux_tau( 1, l ) * tau_TL( k-1, l ) / &
            !                       -----------------------------------
                                              tau( k-1, l )**2          )
          END DO
        END IF

      END IF



      ! ----------------------------------------------------
      ! If the current channel is a SOLAR SENSITIVE channel,
      !   AND
      ! the solar angle is valid, then calculate the
      ! transmittance for direct solar.
      ! ----------------------------------------------------

      IF ( is_solar_channel( channel_index( l ) ) == 1 .AND. &
           secant_solar_angle < MAX_SECANT_SOLAR_ANGLE       ) THEN

        valid_solar = 1

        CALL compute_transmittance_TL( &
                                       ! -- Forward input
                                       solar_tau_absorber( 0:, : ),    &   ! Input, 0:K x J
                                       solar_tau_predictor( :, : ),    &   ! Input, I x K
                                       solar_tau( :, l ),              &   ! Input, K

                                       ! -- Tangent-linear input
                                       solar_tau_absorber_TL( 0:, : ), &   ! Input, 0:K x J
                                       solar_tau_predictor_TL( :, : ), &   ! Input, I x K

                                       ! -- Other input
                                       solar_tau_layer_index( :, : ),  &   ! Input, K x J
                                       channel_index( l ),             &   ! Input, scalar
                                       DOWN,                           &   ! Input, scalar

                                       ! -- Tangent-linear output
                                       solar_tau_TL( :, l )            )   ! Output, K

      ELSE

        valid_solar = 0
        solar_tau_TL( :, l ) = ZERO

      END IF


      ! ------------------------------------------------------
      ! Calculate the profile/channel tangent-linear radiances
      ! ------------------------------------------------------

      CALL compute_radiance_TL( &
                                ! -- Forward input
                                layer_t( : ),                  &  ! Input, K

                                surface_temperature,           &  ! Input, scalar
                                surface_emissivity( l ),       &  ! Input, scalar
                                surface_reflectivity( l ),     &  ! Input, scalar

                                tau(      :, l ),              &  ! Input, K
                                flux_tau( :, l ),              &  ! Input, K
                                solar_tau( n_layers, l ),      &  ! Input, scalar

                                layer_radiance( :, l ),        &  ! Input, K
                                downwelling_radiance( l ),     &  ! Input, scalar
                                upwelling_radiance( l ),       &  ! Input, scalar

                                ! -- Tangent-linear input
                                layer_t_TL( : ),               &  ! Input, K

                                surface_temperature_TL,        &  ! Input, scalar
                                surface_emissivity_TL( l ),    &  ! Input, scalar
                                surface_reflectivity_TL( l ),  &  ! Input, scalar

                                tau_TL(      :, l ),           &  ! Input, K
                                flux_tau_TL( :, l ),           &  ! Input, K
                                solar_tau_TL( n_layers, l ),   &  ! Input, scalar

                                ! -- Other input
                                secant_solar_angle,            &  ! Input, scalar
                                valid_solar,                   &  ! Input, scalar
                                channel_index( l ),            &  ! Input, scalar

                                ! -- Tangent-linear output
                                layer_radiance_TL( :, l ),     &  ! Output, K
                                downwelling_radiance_TL( l ),  &  ! Output, scalar
                                upwelling_radiance_TL( l ),    &  ! Output, scalar

                                brightness_temperature_TL( l ) )  ! Output, scalar


    END DO l_channel_loop



    !#--------------------------------------------------------------------------#
    !#                              -- DONE --                                  #
    !#--------------------------------------------------------------------------#

    error_status = SUCCESS

  END FUNCTION tangent_linear_rtm_rank1

END MODULE tangent_linear_model


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
! Revision 1.9  2001/11/07 15:08:07  paulv
! - Changed the logical IF test for the number of input profiles in the
!   [<TANGENT_LINEAR><ADJOINT><K_MATRIX>]_RTM_RANK2() functions from
!     IF ( n_input_profiles >= 1 .AND. n_input_profiles <= n_profiles ) THEN
!   to
!     IF ( n_input_profiles > 0 .AND. n_input_profiles <= n_profiles ) THEN
!   The use of both the ">=" and "<=" relational operators with the .AND. I
!   found confusing.
!
! Revision 1.8  2001/11/07 14:57:46  paulv
! - Added check for negative number of channels to TANGENT_LINEAR_RTM_RANK2() function.
! - Added profile loop CYCLE statement to TANGENT_LINEAR_RTM_RANK2() function.
! - Added check for negative number of channels to TANGENT_LINEAR_RTM_RANK1() function.
!
! Revision 1.7  2001/10/01 20:26:20  paulv
! - Overloaded the TANGENT_LINEAR_RTM function to accept both a
!   single or group of profiles. Contained PRIVATE functions are now
!     o TANGENT_LINEAR_RTM_RANK1 for single profile input
!     o TANGENT_LINEAR_RTM_RANK2 for multiple profile input
! - Put N_INPUT_PROFILES optional argument back in TANGENT_LINEAR_RTM
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
!   as the check is now done in the rank-1 TANGENT_LINEAR_RTM function. This eliminates
!   the need to fully populate the input arrays with data when only certain
!   "chunks" may be processed. Previously, the use of ANY() could generate
!   and error if the full surface_temperature array was not initialised.
! - Added "Name" to RCS keyword list.
!
! Revision 1.6  2001/09/04 21:29:11  paulv
! - Updated documentation.
!
! Revision 1.5  2001/08/31 21:33:04  paulv
! - Added check for negative profile and surface data in TANGENT_LINEAR_RTM.
! - Maximum solar angle secant is no longer calculated in TANGENT_LINEAR_RTM but
!   is declared as a parameter in the PARAMETERS module.
!
! Revision 1.4  2001/08/16 16:37:52  paulv
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
! Revision 1.3  2001/08/01 17:00:16  paulv
! - Altered function declaration to avoid more than the standard 39
!   continuation lines allowed in Fortran 90.
! - Updated input argument checking. Now consistent with other model
!   components.
!
! Revision 1.2  2001/07/12 18:46:06  paulv
! - Commented out informational message output at start of function.
! - Corrected bug in the calculation of the thermal flux transmittance for
!   the infrared channels. The call to COMPUTE_TRANSMITTANCE_TL used the
!   absorber bracketing index array TAU_LAYER_INDEX (for the upwelling
!   transmittance) rather than FLUX_TAU_LAYER_INDEX.
!
! Revision 1.1  2001/05/29 16:36:02  paulv
! Initial checkin
!
!
!
