!------------------------------------------------------------------------------
!M+
! NAME:
!       radiance
!
! PURPOSE:
!       RT model radiance module
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE radiance
!
! OUTPUTS:
!       None.
!
! MODULES:
!       type_kinds:              Module containing data type kind definitions.
!
!       parameters:              Module containing parameter definitions for
!                                the RT model.
!
!       spectral_coefficients:   Module containing the RT model spectral
!                                coefficients.
!
!       sensor_planck_routines:  Module containing all the forward, tangent-
!                                linear, and adjoint Planck radiance and
!                                temperature subroutines. 
!
! CONTAINS:
!       compute_radiance:        PUBLIC subroutine to calculate the channel TOA
!                                radiance and brightness temperature.
!
!       compute_radiance_TL:     PUBLIC subroutine to calculate the tangent-
!                                linear TOA radiance and brightness temperature.
!
!       compute_radiance_AD:     PUBLIC subroutine to calculate the adjoint of
!                                the TOA radiance and brightness temperature.
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

MODULE radiance


  ! ---------------------
  ! Module use statements
  ! ---------------------

  USE type_kinds, ONLY : fp_kind
  USE parameters
  USE spectral_coefficients
  USE sensor_planck_routines


  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------

  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------

  PRIVATE
  PUBLIC  :: compute_radiance
  PUBLIC  :: compute_radiance_TL
  PUBLIC  :: compute_radiance_AD


CONTAINS


!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_radiance
!
! PURPOSE:
!       PUBLIC subroutine to calculate the TOA radiance and brightness temperature.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_radiance( temperature,           &  ! Input, K
!
!                              surface_temperature,   &  ! Input, scalar
!                              surface_emissivity,    &  ! Input, scalar
!                              surface_reflectivity,  &  ! Input, scalar
!
!                              tau,                   &  ! Input, K
!                              flux_tau,              &  ! Input, K
!                              solar_tau,             &  ! Input, scalar
!
!                              secant_solar_angle,    &  ! Input, scalar
!                              valid_solar,           &  ! Input, scalar
!                              channel_index,         &  ! Input, scalar
!
!                              layer_radiance,        &  ! Output, K
!                              downwelling_radiance,  &  ! Output, scalar
!                              upwelling_radiance,    &  ! Output, scalar
!
!                              brightness_temperature )  ! Output, scalar
!
!
! INPUT ARGUMENTS:
!       temperature:             Profile LAYER average temperature array.
!                                UNITS:      Kelvin
!                                TYPE:       Float
!                                DIMENSION:  K
!                                ATTRIBUTES: INTENT( IN )
!
!       surface_temperature:     Surface boundary temperature.
!                                UNITS:      Kelvin
!                                TYPE:       Real
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( IN )
!
!       surface_emissivity:      Surface boundary emissivity
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( IN )
!
!       surface_reflectivity:    Surface boundary reflectivity
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( IN )
!
!       tau:                     Layer-to-space transmittance profile for
!                                a particular satellite view angle.
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  K
!                                ATTRIBUTES: INTENT( IN )
!
!       flux_tau:                Layer-to-surface transmittance profile for
!                                either the diffuse approximation angle (IR)
!                                or the satellite view angle (MW). The latter
!                                assumes specular reflectivity.
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  K
!                                ATTRIBUTES: INTENT( IN )
!
!       solar_tau:               Total space-to-surface transmittance at the
!                                solar zenith angle.
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( IN )
!
!       secant_solar_angle:      Secant of the solar zenith angle corresponding
!                                to that used in calculating the total solar
!                                transmittance.
!                                UNITS:      None
!                                TYPE:       Real
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( IN )
!
!       valid_solar:             Flag indicating if the solar component should
!                                be included.
!                                If = 0, no solar (if sensor channel frequency
!                                        is less than a preset cutoff or if solar
!                                        zenith angle is greater than its preset
!                                         cutoff.)
!                                   = 1, include solar
!                                UNITS:      None.
!                                TYPE:       Integer
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( IN )
!
!       channel_index:           Channel index id. This is a unique index
!                                to a (supported) sensor channel.
!                                UNITS:      None
!                                TYPE:       Integer
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       layer_radiance:          Channel Planck radiance for every input layer.
!                                UNITS:      mW/(m^2.sr.cm^-1)
!                                TYPE:       Real
!                                DIMENSION:  K.
!                                ATTRIBUTES: INTENT( OUT )
!
!       downwelling_radiance:    Channel radiance at surface due to downwelling
!                                flux and solar components.
!                                UNITS:      mW/(m^2.sr.cm^-1)
!                                TYPE:       Real
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( OUT )
!
!       upwelling_radiance:      Channel TOA radiance simulating the satellite
!                                sensor measurement. This is composed of the
!                                reflected downwelling propagated through the
!                                atmosphere as well as the upwelling only component.
!                                UNITS:      mW/(m^2.sr.cm^-1)
!                                TYPE:       Real
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( OUT )
!
!       brightness_temperature:  Channel TOA brightness temperature.
!                                UNITS:      Kelvin
!                                TYPE:       Real
!                                DIMENSION:  Scalar
!                                ATTRIBUTES: INTENT( OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! CALLS:
!       sensor_planck_radiance:    Function to compute the Planck radiance
!                                  for a specified channel given the temperature.
!                                  SOURCE: SENSOR_PLANCK_ROUTINES module
!
!       sensor_planck_temperature: Function to compute the Planck temperature
!                                  for a specified channel given the radiance.
!                                  SOURCE: SENSOR_PLANCK_ROUTINES module
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
!       The downwelling radiance is first initialised to the space emissisio
!       boundary term using precalculated cosmic background radiances,
!
!         R_down = CBR * flux_tau(1)
!
!       where the emissivity of space is implicitly assumed to be 1.0 and
!       flux_tau(1) is the space-to-ground transmittance.
!
!       The contributions of all the layers EXCEPT THE SURFACE LAYER to the
!       downwelling flux is accumulated,
!
!                            __K-1
!                           \
!         R_down = R_down +  >  B(k) * dflux_tau(k)
!                           /__
!                              k=1
!
!       The surface layer contribution is then added explicitly,
!
!         R_down = R_down + ( B(K) * ( 1 - flux_tau(K) ) )
!
!       to avoid exceeding the arrays bounds of flux_tau 
!       (i.e. flux_tau(K+1) == 1.0 ) or requiring flux_tau to be
!       dimensioned 0:K.
!
!       The solar term is then added if required,
!
!         R_down = R_down + ( solar_irradiance * solar_tau * COS(solar_theta) )
!
!       The downwelling radiance is then reflected off the surface, added
!       to the surface emission term, propagated upwards through the atmosphere
!       and used to initialise the upwelling radiance term,
!
!         R_up = ( ( e_sfc * B_sfc ) + ( r_sfc * R_down ) ) * tau(K)
!
!       The contributions of all the layers EXCEPT THE TOP LAYER to the
!       upwelling radiance is accumulated,
!
!                        __ 2
!                       \
!         R_up = R_up +  >  B(k) * dtau(k)
!                       /__
!                          k=K
!
!       The top layer contribution is then added explicitly,
!
!         R_up = R_up + ( B(1) * ( 1 - tau(1) ) )
!
!       to avoid exceeding the arrays bounds of tau (i.e. tau(0) == 1.0 )
!       or requiring tau to be dimensioned 0:K.
!
!       The final upwelling radiance is then converted to a brightness
!       temperature.
!
!S-      
!--------------------------------------------------------------------------------


  SUBROUTINE compute_radiance( temperature,           &  ! Input, K      

                               surface_temperature,   &  ! Input, scalar 
                               surface_emissivity,    &  ! Input, scalar 
                               surface_reflectivity,  &  ! Input, scalar 

                               tau,                   &  ! Input, K      
                               flux_tau,              &  ! Input, K      
                               solar_tau,             &  ! Input, scalar 

                               secant_solar_angle,    &  ! Input, scalar 
                               valid_solar,           &  ! Input, scalar 
                               channel_index,         &  ! Input, scalar 

                               layer_radiance,        &  ! Output, K     
                               downwelling_radiance,  &  ! Output, scalar
                               upwelling_radiance,    &  ! Output, scalar

                               brightness_temperature )  ! Output, scalar



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: temperature

    REAL( fp_kind ),                 INTENT( IN )  :: surface_temperature
    REAL( fp_kind ),                 INTENT( IN )  :: surface_emissivity
    REAL( fp_kind ),                 INTENT( IN )  :: surface_reflectivity

    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: tau
    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: flux_tau
    REAL( fp_kind ),                 INTENT( IN )  :: solar_tau

    REAL( fp_kind ),                 INTENT( IN )  :: secant_solar_angle
    INTEGER,                         INTENT( IN )  :: valid_solar
    INTEGER,                         INTENT( IN )  :: channel_index

    REAL( fp_kind ), DIMENSION( : ), INTENT( OUT ) :: layer_radiance
    REAL( fp_kind ),                 INTENT( OUT ) :: downwelling_radiance
    REAL( fp_kind ),                 INTENT( OUT ) :: upwelling_radiance

    REAL( fp_kind ),                 INTENT( OUT ) :: brightness_temperature

 

    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_RADIANCE'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: k, n_layers
    INTEGER :: l

    REAL( fp_kind ) :: surface_B


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC SIZE



    !#--------------------------------------------------------------------------#
    !#                   -- Determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    n_layers = SIZE( temperature )



    !#--------------------------------------------------------------------------#
    !#              -- Calculate the downwelling thermal flux --                #
    !#--------------------------------------------------------------------------#

    ! -- Assign the channel index to a short name
    l = channel_index


    ! --------------------------------------------
    ! Initialise the downwelling radiance to the
    ! space emission boundary term reaching the
    ! surface. Thhe cosmic background radiance is
    ! zero for infrared channels and precalculated
    ! for microwave channels. The emissivity of
    ! space is assumed to be 1.0.
    !
    ! Cosmic background data from the
    ! SPECTRAL_COEFFICIENTS module
    ! --------------------------------------------

    downwelling_radiance = cosmic_background_radiance( l ) * flux_tau( 1 )


    ! --------------------------------
    ! Loop over layers from TOA->SFC-1
    ! --------------------------------

    k_down_layer_loop: DO k = 1, n_layers - 1

      ! -- Calculate the Planck layer radiance
      CALL sensor_planck_radiance( l,                  &  ! Input
                                   temperature( k ),   &  ! Input
                                   layer_radiance( k ) )  ! Output

      ! -- Accumulate absorption and emission for current layer.
      ! -- LTE assumed.
      downwelling_radiance = downwelling_radiance + &
                             ( layer_radiance( k ) * ( flux_tau( k+1 ) - flux_tau( k ) ) )

    END DO k_down_layer_loop


    ! ----------------------------------
    ! Flux bottom layer (closest to SFC)
    ! ----------------------------------

    ! -- Lowest layer Planck radiance
    CALL sensor_planck_radiance( l,                         &  ! Input
                                 temperature( n_layers ),   &  ! Input
                                 layer_radiance( n_layers ) )  ! Output

    ! -- Contribution of lowest layer. Note that at the
    ! -- surface, a transmittance of 1.0 is used.
    downwelling_radiance = downwelling_radiance + &
                           ( layer_radiance( n_layers ) * ( ONE - flux_tau( n_layers ) ) )



    !#--------------------------------------------------------------------------#
    !#                -- Calculate the downwelling solar terms --               #
    !#--------------------------------------------------------------------------#

    solar_term: IF ( valid_solar == 1 ) THEN

      downwelling_radiance = downwelling_radiance + &

                             ( solar_irradiance( l ) * solar_tau / &
      !                        ---------------------------------
                                      secant_solar_angle         )

    END IF solar_term



    !#--------------------------------------------------------------------------#
    !#   -- Reflect the downwelling radiance and add the surface emission --    #
    !#   -- and use it to initialise the upwelling radiance               --    #
    !#--------------------------------------------------------------------------#

    ! -- Calculate the surface term
    CALL sensor_planck_radiance( l,                   &  ! Input
                                 surface_temperature, &  ! Input
                                 surface_B            )  ! Output

    ! -- Initialise upwelling radiance
    upwelling_radiance = ( ( surface_emissivity   * surface_B            ) + &
                           ( surface_reflectivity * downwelling_radiance )   ) * tau( n_layers )



    !#--------------------------------------------------------------------------#
    !#                  -- Calculate the upwelling radiance --                  #
    !#--------------------------------------------------------------------------#

    ! --------------------------------
    ! Loop over layers from SFC->TOA-1
    ! --------------------------------

    k_up_layer_loop: DO k = n_layers, 2, -1

      upwelling_radiance = upwelling_radiance + &
                           ( layer_radiance( k ) * ( tau( k-1 ) - tau( k ) ) )

    END DO k_up_layer_loop


    ! --------------------------
    ! Top layer (closest to TOA)
    ! --------------------------

    upwelling_radiance = upwelling_radiance + &
                         ( layer_radiance( 1 ) * ( ONE - tau( 1 ) ) )



    !#--------------------------------------------------------------------------#
    !#           -- Convert the radiances to brightness temperatures --         #
    !#--------------------------------------------------------------------------#

    CALL sensor_planck_temperature( l,                     &  ! Input
                                    upwelling_radiance,    &  ! Input
                                    brightness_temperature )  ! Output

  END SUBROUTINE compute_radiance



!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_radiance_TL
!
! PURPOSE:
!       PUBLIC subroutine to calculate the tangent-linear TOA radiance and
!       brightness temperature.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_radiance_TL( &
!
!                                 ! -- Forward inputs
!                                 temperature,              &  ! Input, K
!
!                                 surface_temperature,      &  ! Input, scalar
!                                 surface_emissivity,       &  ! Input, scalar
!                                 surface_reflectivity,     &  ! Input, scalar
!
!                                 tau,                      &  ! Input, K
!                                 flux_tau,                 &  ! Input, K
!                                 solar_tau,                &  ! Input, scalar
!
!                                 layer_radiance,           &  ! Input, K
!                                 downwelling_radiance,     &  ! Input, scalar
!                                 upwelling_radiance,       &  ! Input, scalar
!
!                                 ! -- Tangent-linear inputs
!                                 temperature_TL,           &  ! Input, K
!
!                                 surface_temperature_TL,   &  ! Input, scalar
!                                 surface_emissivity_TL,    &  ! Input, scalar
!                                 surface_reflectivity_TL,  &  ! Input, scalar
!
!                                 tau_TL,                   &  ! Input, K
!                                 flux_tau_TL,              &  ! Input, K
!                                 solar_tau_TL,             &  ! Input, scalar
!
!                                 ! -- Other inputs
!                                 secant_solar_angle,       &  ! Input, scalar
!                                 valid_solar,              &  ! Input, scalar
!                                 channel_index,            &  ! Input, scalar
!
!                                 ! -- Tangent-linear outputs
!                                 layer_radiance_TL,        &  ! Output, K
!                                 downwelling_radiance_TL,  &  ! Output, scalar
!                                 upwelling_radiance_TL,    &  ! Output, scalar
!
!                                 brightness_temperature_TL )  ! Output, scalar
!
!
! INPUT ARGUMENTS:
!       temperature:               Profile LAYER average temperature array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Float
!                                  DIMENSION:  K
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_temperature:       Surface boundary temperature.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_emissivity:        Surface boundary emissivity
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_reflectivity:      Surface boundary reflectivity
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       tau:                       Layer-to-space transmittance profile for
!                                  a particular satellite view angle.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  K.
!                                  ATTRIBUTES: INTENT( IN )
!
!       flux_tau:                  Layer-to-surface transmittance profile for
!                                  either the diffuse approximation angle (IR)
!                                  or the satellite view angle (MW). The latter
!                                  assumes specular reflectivity.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  K.
!                                  ATTRIBUTES: INTENT( IN )
!
!       solar_tau:                 Total space-to-surface transmittance at the
!                                  solar zenith angle.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_radiance:            Channel Planck radiance for every layer.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  K.
!                                  ATTRIBUTES: INTENT( IN )
!
!       downwelling_radiance:      Channel radiance at surface due to downwelling
!                                  flux and solar components.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       upwelling_radiance:        Channel TOA radiance simulating the satellite
!                                  sensor measurement. This is composed of the
!                                  reflected downwelling propagated through the
!                                  atmosphere as well as the upwelling only component.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       temperature_TL:            Tangent-linear temperature profile.
!                                  UNITS:      Kelvin
!                                  TYPE:       Float
!                                  DIMENSION:  K.
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_temperature_TL:    Tangent-linear surface boundary temperature.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_emissivity_TL:     Tangent-linear surface boundary emissivity
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_reflectivity_TL:   Tangent-linear surface boundary reflectivity
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       tau_TL:                    Tangent-linear layer-to-space transmittance
!                                  profile.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  K.
!                                  ATTRIBUTES: INTENT( IN )
!
!       flux_tau_TL:               Tangent-linear layer-to-surface flux transmittance
!                                  profile.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  K.
!                                  ATTRIBUTES: INTENT( IN )
!
!       solar_tau_TL:              Tangent-linear total space-to-surface solar
!                                  transmittance.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       secant_solar_angle:        Secant of the solar zenith angle corresponding
!                                  to that used in calculating the total solar
!                                  transmittance.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       valid_solar:               Flag indicating if the solar component should
!                                  be included.
!                                  If = 0, no solar (if sensor channel frequency
!                                          is less than a preset cutoff or if solar
!                                          zenith angle is greater than its preset
!                                           cutoff.)
!                                     = 1, include solar
!                                  UNITS:      None.
!                                  TYPE:       Integer
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       channel_index:             Channel index id. This is a unique index
!                                  to a (supported) sensor channel.
!                                  UNITS:      None
!                                  TYPE:       Integer
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       layer_radiance_TL:         Tangent-linear channel Planck radiance for
!                                  every layer.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  K.
!                                  ATTRIBUTES: INTENT( OUT )
!
!       downwelling_radiance_TL:   Tangent-linear channel radiance at surface
!                                  due to downwelling flux and solar components.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( OUT )
!
!       upwelling_radiance_TL:     Tangent-linear channel TOA radiance.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( OUT )
!
!       brightness_temperature_TL: Tangent-linear channel TOA brightness
!                                  temperature.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! CALLS:
!       sensor_planck_radiance_TL:    Function to compute the tangent-linear
!                                     Planck radiance.
!                                     SOURCE: SENSOR_PLANCK_ROUTINES module
!
!       sensor_planck_temperature_TL: Function to compute the tangent-linear
!                                     Planck temperature.
!                                     SOURCE: SENSOR_PLANCK_ROUTINES module
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
!       The downwelling radiance is first initialised to the space emission
!       boundary term using precalculated cosmic background radiances,
!
!         R_down_TL = CBR * flux_tau_TL(1)
!
!       and the emissivity of space is implicitly assumed to be 1.0 and
!       flux_tau_TL(1) is the tangent-linear form of the space-to-ground
!       transmittance.
!
!       The contributions of all the layers EXCEPT THE SURFACE LAYER to the
!       downwelling flux is accumulated,
!
!                                  __K-1
!                                 \
!         R_down_TL = R_down_TL +  >  ( B(k) * dflux_tau_TL(k) ) + ( B_TL(k) * dflux_tau(k) )
!                                 /__
!                                    k=1
!
!       The surface layer contribution is then added explicitly,
!
!         R_down_TL = R_down_TL + ( B(K) * ( -flux_tau_TL(K) ) ) + ( B_TL(K) * ( 1 - flux_tau(K) ) )
!
!       to avoid exceeding the arrays bounds of flux_tau and flux_tau_TL
!       (i.e. flux_tau(K+1) == 1.0 ) or requiring them to be
!       dimensioned 0:K.
!
!       The solar term is then added if required,
!
!         R_down_TL = R_down_TL + ( solar_irradiance * solar_tau_TL * COS(solar_theta) )
!
!       The downwelling radiance is then reflected off the surface, added
!       to the surface emission term, propagated upwards through the atmosphere
!       and used to initialise the upwelling radiance term,
!
!         R_up_TL = ( e_sfc    * B_sfc_TL  * tau(K) ) + &
!                   ( e_sfc_TL * B_sfc     * tau(K) ) + &
!                   ( r_sfc    * R_down_TL * tau(K) ) + &
!                   ( r_sfc_TL * R_down    * tau(K) ) + &
!                   ( ( ( e_sfc * B_sfc ) + ( r_sfc * R_down ) ) * tau_TL(K) )
!
!       The contributions of all the layers EXCEPT THE TOP LAYER to the
!       upwelling radiance is accumulated,
!
!                              __ 2
!                             \
!         R_up_TL = R_up_TL +  >  ( B(k) * dtau_TL(k) ) + ( B_TL(k) * dtau(k) )
!                             /__
!                                k=K
!
!       The top layer contribution is then added explicitly,
!
!         R_up_TL = R_up_TL + ( B(1) * ( -tau_TL(1) ) ) + ( B_TL(1) * ( 1 - tau(1) ) )
!
!       to avoid exceeding the arrays bounds of tau (i.e. tau(0) == 1.0 )
!       or tau_TL or requiring them to be dimensioned 0:K.
!
!       The final tangent-linear upwelling radiance is then converted to a
!       tangent-linear  brightness temperature.
!
!S-      
!--------------------------------------------------------------------------------


  SUBROUTINE compute_radiance_TL( &
                                  ! -- Forward inputs
                                  temperature,              &  ! Input, K

                                  surface_temperature,      &  ! Input, scalar
                                  surface_emissivity,       &  ! Input, scalar
                                  surface_reflectivity,     &  ! Input, scalar

                                  tau,                      &  ! Input, K
                                  flux_tau,                 &  ! Input, K
                                  solar_tau,                &  ! Input, scalar

                                  layer_radiance,           &  ! Input, K
                                  downwelling_radiance,     &  ! Input, scalar
                                  upwelling_radiance,       &  ! Input, scalar

                                  ! -- Tangent-linear inputs
                                  temperature_TL,           &  ! Input, K

                                  surface_temperature_TL,   &  ! Input, scalar
                                  surface_emissivity_TL,    &  ! Input, scalar
                                  surface_reflectivity_TL,  &  ! Input, scalar

                                  tau_TL,                   &  ! Input, K
                                  flux_tau_TL,              &  ! Input, K
                                  solar_tau_TL,             &  ! Input, scalar

                                  ! -- Other inputs
                                  secant_solar_angle,       &  ! Input, scalar
                                  valid_solar,              &  ! Input, scalar
                                  channel_index,            &  ! Input, scalar

                                  ! -- Tangent-linear outputs
                                  layer_radiance_TL,        &  ! Output, K
                                  downwelling_radiance_TL,  &  ! Output, scalar
                                  upwelling_radiance_TL,    &  ! Output, scalar

                                  brightness_temperature_TL )  ! Output, scalar



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward input
    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: temperature

    REAL( fp_kind ),                 INTENT( IN )  :: surface_temperature
    REAL( fp_kind ),                 INTENT( IN )  :: surface_emissivity
    REAL( fp_kind ),                 INTENT( IN )  :: surface_reflectivity

    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: tau
    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: flux_tau
    REAL( fp_kind ),                 INTENT( IN )  :: solar_tau

    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: layer_radiance
    REAL( fp_kind ),                 INTENT( IN )  :: downwelling_radiance
    REAL( fp_kind ),                 INTENT( IN )  :: upwelling_radiance

    ! -- Tangent-linear input
    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: temperature_TL

    REAL( fp_kind ),                 INTENT( IN )  :: surface_temperature_TL
    REAL( fp_kind ),                 INTENT( IN )  :: surface_emissivity_TL
    REAL( fp_kind ),                 INTENT( IN )  :: surface_reflectivity_TL

    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: tau_TL
    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )  :: flux_tau_TL
    REAL( fp_kind ),                 INTENT( IN )  :: solar_tau_TL

    ! -- Other input
    REAL( fp_kind ),                 INTENT( IN )  :: secant_solar_angle
    INTEGER,                         INTENT( IN )  :: valid_solar
    INTEGER,                         INTENT( IN )  :: channel_index

    ! -- Tangent-linear output
    REAL( fp_kind ), DIMENSION( : ), INTENT( OUT ) :: layer_radiance_TL
    REAL( fp_kind ),                 INTENT( OUT ) :: downwelling_radiance_TL
    REAL( fp_kind ),                 INTENT( OUT ) :: upwelling_radiance_TL

    REAL( fp_kind ),                 INTENT( OUT ) :: brightness_temperature_TL

 

    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_RADIANCE_TL'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: k, n_layers
    INTEGER :: l

    REAL( fp_kind ) :: surface_B
    REAL( fp_kind ) :: surface_B_TL


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC SIZE



    !#--------------------------------------------------------------------------#
    !#                   -- Determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    n_layers = SIZE( temperature )



    !#--------------------------------------------------------------------------#
    !#       -- Calculate the tangent-linear downwelling thermal flux --        #
    !#--------------------------------------------------------------------------#

    ! -- Assign the channel index to a short name
    l = channel_index


    ! ---------------------------------------------
    ! Initialise the tangent-linear downwelling
    ! radiance to the space emission boundary term
    ! reaching the surface. The cosmic background
    ! radiance is zero for infrared channels and
    ! precalculated for microwave channels. The
    ! emissivity of space is assumed to be 1.0.
    !
    ! Cosmic background data from the
    ! SPECTRAL_COEFFICIENTS module
    ! ---------------------------------------------

    downwelling_radiance_TL = cosmic_background_radiance( l ) * flux_tau_TL( 1 )


    ! --------------------------------
    ! Loop over layers from TOA->SFC-1
    ! --------------------------------

    k_down_layer_loop: DO k = 1, n_layers - 1

      ! -- Calculate the tangent-linear layer Planck radiance
      CALL sensor_planck_radiance_TL( l,                     &  ! Input
                                      temperature( k ),      &  ! Input
                                      temperature_TL( k ),   &  ! Input
                                      layer_radiance_TL( k ) )  ! Output

      ! -- Accumulate tangent-linear absorption and emission for current layer.
      ! -- LTE assumed.
      downwelling_radiance_TL = downwelling_radiance_TL + &
                                ( layer_radiance(k)    * ( flux_tau_TL(k+1) - flux_tau_TL(k) ) ) + &
                                ( layer_radiance_TL(k) * ( flux_tau(k+1)    - flux_tau(k)    ) )

    END DO k_down_layer_loop


    ! ----------------------------------
    ! Flux bottom layer (closest to SFC)
    ! ----------------------------------

    ! -- Lowest layer tangent-linear Planck radiance
    CALL sensor_planck_radiance_TL( l,                            &  ! Input
                                    temperature( n_layers ),      &  ! Input
                                    temperature_TL( n_layers ),   &  ! Input
                                    layer_radiance_TL( n_layers ) )  ! Output

    ! -- Contribution of lowest layer. Note that at the
    ! -- surface, a transmittance of 1.0 is used.
    downwelling_radiance_TL = downwelling_radiance_TL + &
                              ( layer_radiance(n_layers)    * (      -flux_tau_TL(n_layers) ) ) + &
                              ( layer_radiance_TL(n_layers) * ( ONE - flux_tau(n_layers)    ) )



    !#--------------------------------------------------------------------------#
    !#               -- Calculate the tangent-linear solar term --              #
    !#--------------------------------------------------------------------------#

    solar_term: IF ( valid_solar == 1 ) THEN

      downwelling_radiance_TL = downwelling_radiance_TL + &

                                ( solar_irradiance(l) * solar_tau_TL / &
      !                           ------------------------------------
                                           secant_solar_angle          )

    END IF solar_term



    !#--------------------------------------------------------------------------#
    !#      -- Reflect the tangent-linear downwelling radiance, add the  --     #
    !#      -- surface emission and use it to initialise the tangent-    --     #
    !#      -- linear upwelling radiance                                 --     #
    !#--------------------------------------------------------------------------#

    ! -- Calculate the surface terms
    CALL sensor_planck_radiance( l,                   &  ! Input
                                 surface_temperature, &  ! Input
                                 surface_B            )  ! Output

    CALL sensor_planck_radiance_TL( l,                      &  ! Input
                                    surface_temperature,    &  ! Input
                                    surface_temperature_TL, &  ! Input
                                    surface_B_TL            )  ! Output

    ! -- Initialise the tangent-linear upwelling radiance
    upwelling_radiance_TL = ( tau(n_layers) * surface_emissivity   * surface_B_TL            ) + &
                            ( tau(n_layers) * surface_reflectivity * downwelling_radiance_TL ) + &

                            ( ((surface_emissivity   * surface_B           ) + &
                               (surface_reflectivity * downwelling_radiance)) * tau_TL(n_layers) ) + &

                            ( tau(n_layers) * surface_B            * surface_emissivity_TL   ) + &
                            ( tau(n_layers) * downwelling_radiance * surface_reflectivity_TL )



    !#--------------------------------------------------------------------------#
    !#           -- Calculate the tangent-linear upwelling radiance --          #
    !#--------------------------------------------------------------------------#

    ! --------------------------------
    ! Loop over layers from SFC->TOA-1
    ! --------------------------------

    k_up_layer_loop: DO k = n_layers, 2, -1

      upwelling_radiance_TL = upwelling_radiance_TL + &
                              ( layer_radiance(k)    * ( tau_TL(k-1) - tau_TL(k) ) ) + &
                              ( layer_radiance_TL(k) * ( tau(k-1)    - tau(k)    ) )

    END DO k_up_layer_loop


    ! --------------------------
    ! Top layer (closest to TOA)
    ! --------------------------

    upwelling_radiance_TL = upwelling_radiance_TL + &
                            ( layer_radiance(1)    * (      -tau_TL(1) ) ) + &
                            ( layer_radiance_TL(1) * ( ONE - tau(1)    ) )



    !#--------------------------------------------------------------------------#
    !#    -- Convert the tangent-linear radiance to brightness temperature --   #
    !#--------------------------------------------------------------------------#

    CALL sensor_planck_temperature_TL( l,                        &  ! Input
                                       upwelling_radiance,       &  ! Input
                                       upwelling_radiance_TL,    &  ! Input
                                       brightness_temperature_TL )  ! Output

  END SUBROUTINE compute_radiance_TL






!--------------------------------------------------------------------------------
!S+
! NAME:
!       compute_radiance_AD
!
! PURPOSE:
!       PUBLIC subroutine to calculate the adjoint of the TOA radiance and
!       brightness temperature.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL compute_radiance_AD( &
!                                 ! -- Forward inputs
!                                 temperature,               &  ! Input, K
!
!                                 surface_temperature,       &  ! Input, scalar
!                                 surface_emissivity,        &  ! Input, scalar
!                                 surface_reflectivity,      &  ! Input, scalar
!
!                                 tau,                       &  ! Input, K
!                                 flux_tau,                  &  ! Input, K
!                                 solar_tau,                 &  ! Input, scalar
!
!                                 layer_radiance,            &  ! Input, K
!                                 downwelling_radiance,      &  ! Input, scalar
!                                 upwelling_radiance,        &  ! Input, scalar
!
!                                 ! -- Adjoint inputs
!                                 layer_radiance_AD,         &  ! In/Output, K
!                                 downwelling_radiance_AD,   &  ! In/Output, scalar
!                                 upwelling_radiance_AD,     &  ! In/Output, scalar
!
!                                 brightness_temperature_AD, &  ! In/Output, scalar
!
!                                 ! -- Other inputs
!                                 secant_solar_angle,        &  ! Input, scalar
!                                 valid_solar,               &  ! Input, scalar
!                                 channel_index,             &  ! Input, scalar
!
!                                 ! -- Adjoint outputs
!                                 temperature_AD,            &  ! In/Output, K
!
!                                 surface_temperature_AD,    &  ! In/Output, scalar
!                                 surface_emissivity_AD,     &  ! In/Output, scalar
!                                 surface_reflectivity_AD,   &  ! In/Output, scalar
!
!                                 tau_AD,                    &  ! In/Output, K
!                                 flux_tau_AD,               &  ! In/Output, K
!                                 solar_tau_AD               )  ! In/Output, scalar
!
!
! INPUT ARGUMENTS:
!       temperature:               Profile LAYER average temperature array.
!                                  UNITS:      Kelvin
!                                  TYPE:       Float
!                                  DIMENSION:  K
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_temperature:       Surface boundary temperature.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_emissivity:        Surface boundary emissivity
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_reflectivity:      Surface boundary reflectivity
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       tau:                       Layer-to-space transmittance profile for
!                                  a particular satellite view angle.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  K
!                                  ATTRIBUTES: INTENT( IN )
!
!       flux_tau:                  Layer-to-surface transmittance profile for
!                                  either the diffuse approximation angle (IR)
!                                  or the satellite view angle (MW). The latter
!                                  assumes specular reflectivity.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  K
!                                  ATTRIBUTES: INTENT( IN )
!
!       solar_tau:                 Total space-to-surface transmittance at the
!                                  solar zenith angle.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_radiance:            Channel Planck radiance for every layer.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  K
!                                  ATTRIBUTES: INTENT( IN )
!
!       downwelling_radiance:      Channel radiance at surface due to downwelling
!                                  flux and solar components.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       upwelling_radiance:        Channel TOA radiance simulating the satellite
!                                  sensor measurement. This is composed of the
!                                  reflected downwelling propagated through the
!                                  atmosphere as well as the upwelling only component.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       layer_radiance_AD:         Adjoint of the channel Planck radiance for every
!                                  layer.
!                                  ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  K
!                                  ATTRIBUTES: INTENT( IN )
!
!       downwelling_radiance_AD:   Adjoint of the channel radiance at surface due
!                                  to downwelling flux and solar components.
!                                  ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       upwelling_radiance_AD:     Adjoint of the channel TOA radiance simulating
!                                  the satellite sensor measurement.
!                                  ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                                  UNITS:      mW/(m^2.sr.cm^-1)
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       brightness_temperature_AD: Adjoint of the channel TOA brightness temperature.
!                                  ** THIS ARGUMENT IS SET TO ZERO ON OUTPUT **.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       secant_solar_angle:        Secant of the solar zenith angle corresponding
!                                  to that used in calculating the total solar
!                                  transmittance.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       valid_solar:               Flag indicating if the solar component should
!                                  be included.
!                                  If = 0, no solar (if sensor channel frequency
!                                          is less than a preset cutoff or if solar
!                                          zenith angle is greater than its preset
!                                           cutoff.)
!                                     = 1, include solar
!                                  UNITS:      None.
!                                  TYPE:       Integer
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       channel_index:             Channel index id. This is a unique index
!                                  to a (supported) sensor channel.
!                                  UNITS:      None
!                                  TYPE:       Integer
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       temperature_AD:            Adjoint of the profile LAYER average temperature.
!                                  UNITS:      Kelvin
!                                  TYPE:       Float
!                                  DIMENSION:  K
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_temperature_AD:    Adjoint of the surface boundary temperature.
!                                  UNITS:      Kelvin
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_emissivity_AD:     Adjoint of the surface boundary emissivity
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       surface_reflectivity_AD:   Adjoint of the surface boundary reflectivity
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
!       tau_AD:                    Adjoint of the layer-to-space transmittance
!                                  profile for a particular satellite view angle.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  K
!                                  ATTRIBUTES: INTENT( IN )
!
!       flux_tau_AD:               Adjoint of the layer-to-surface transmittance
!                                  profile for either the diffuse approximation
!                                  angle (IR) or the satellite view angle (MW).
!                                  The latter assumes specular reflectivity.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  K
!                                  ATTRIBUTES: INTENT( IN )
!
!       solar_tau_AD:              Adjoint of the total space-to-surface 
!                                  transmittance at the solar zenith angle.
!                                  UNITS:      None
!                                  TYPE:       Real
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT( IN )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! CALLS:
!       sensor_planck_radiance_AD:    Function to compute the adjoint
!                                     Planck radiance.
!                                     SOURCE: SENSOR_PLANCK_ROUTINES module
!
!       sensor_planck_temperature_AD: Function to compute the adjoint
!                                     Planck temperature.
!                                     SOURCE: SENSOR_PLANCK_ROUTINES module
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       All input adjoint arguments are set to zero on output.
!
! RESTRICTIONS:
!       None.
!
! PROCEDURE:
!       This adjoint code is derived directly from the transpose of the
!       tangent-linear model. Starting with the final statement of the
!       tangent-linear model,
!
!         R_up_TL = R_up_TL + ( B(1) * ( -tau_TL(1) ) ) + ( B_TL(1) * ( 1 - tau(1) ) )
!
!       the adjoint of the components is given by,
!
!         tau_AD(1) = -B(1) * R_up_AD
!         B_AD(1)   = ( 1 - tau(1) ) * R_up_AD
!
!       R_up_AD is not set to zero at this point as it is being summed. The 
!       surface to space tangent-linear summation is given by,
!
!                              __ 2
!                             \
!         R_up_TL = R_up_TL +  >  ( B(k) * dtau_TL(k) ) + ( B_TL(k) * dtau(k) )
!                             /__
!                                k=K
!
!       The adjoint of the components are,
!
!         B_AD(k)     = dtau(k) * R_up_AD
!         tau_AD(k)   = tau_AD(k)   - ( B(k) * R_up_AD )
!         tau_AD(k-1) = tau_AD(k-1) + ( B(k) * R_up_AD )
!         T_AD(k)     = planck_AD(T(k),B_AD(k))
!         B_AD(k)     = 0.0
!
!       Next comes the tangent-linea surface term,
!
!         R_up_TL = ( e_sfc    * B_sfc_TL  * tau(K) ) + &
!                   ( e_sfc_TL * B_sfc     * tau(K) ) + &
!                   ( r_sfc    * R_down_TL * tau(K) ) + &
!                   ( r_sfc_TL * R_down    * tau(K) ) + &
!                   ( ( ( e_sfc * B_sfc ) + ( r_sfc * R_down ) ) * tau_TL(K) )
!
!       with the adjoints,
!
!         r_AD      = tau(K) * R_down * R_up_AD
!         e_AD      = tau(K) * B(sfc) * R_up_AD
!         tau_AD(K) = tau_AD(K) + ( e.B(sfc) + r.R_down ) * R_up_AD
!         R_down_AD = R_down_AD + ( tau(K) * r * R_up_AD )
!         B_AD(sfc) =               tau(K) * e * R_up_AD
!         R_up_AD   = 0.0
!         T_AD(sfc) = planck_AD(T(sfc),B_AD(sfc))
!
!       The tangent-linear solar term, if used, is
!
!         R_down_TL = R_down_TL + ( solar_irradiance * solar_tau_TL * COS(solar_theta) )
!
!       with its adjoint being,
!
!         solar_tau_AD = solar_tau_AD + ( solar_irradiance * COS(solar_theta) * R_down_AD )
!
!       The tangent-linear surface layer contribution to the downwelling flux
!       is calculated separately,
!
!         R_down_TL = R_down_TL + ( B(K) * ( -flux_tau_TL(K) ) ) + ( B_TL(K) * ( 1 - flux_tau(K) ) )
!         
!       and its component adjoints are,
!
!         B_AD(K)        = B_AD(K) + ( 1 - flux_tau(K) ) * R_down_AD
!         flux_tau_AD(K) = flux_tau_AD(K) - ( B(K) * R_down_AD )
!         T_AD(K)        = planck_AD(T(K),B_AD(K))
!         B_AD(K)        = 0.0
!
!       As with the upwelling tangent-linear summation, the downwelling is
!       given by,
!
!                                  __K-1
!                                 \
!         R_down_TL = R_down_TL +  >  ( B(k) * dflux_tau_TL(k) ) + ( B_TL(k) * dflux_tau(k) )
!                                 /__
!                                    k=1
!
!       The adjoint of the components are,
!
!         B_AD(k)          = B_AD(k ) + ( dflux_tau(k) * R_down_AD )
!         flux_tau_AD(k)   = flux_tau_AD(k)   - ( B(k) * R_down_AD )
!         flux_tau_AD(k-1) = flux_tau_AD(k-1) + ( B(k) * R_down_AD )
!         T_AD(k)          = planck_AD(T(k),B_AD(k))
!         B_AD(k)          = 0.0
!
!       The final step is to determine the adjoints of the cosmic background
!       tangent-linear term,
!
!         R_down_TL = CBR * flux_tau_TL(1)
!
!       which is
!
!         flux_tau_AD(1) = flux_tau_AD(1) + ( CBR * R_down_AD )
!         R_down_AD      = 0.0
!
!S-      
!--------------------------------------------------------------------------------


  SUBROUTINE compute_radiance_AD( &
                                  ! -- Forward inputs
                                  temperature,               &  ! Input, K

                                  surface_temperature,       &  ! Input, scalar
                                  surface_emissivity,        &  ! Input, scalar
                                  surface_reflectivity,      &  ! Input, scalar

                                  tau,                       &  ! Input, K
                                  flux_tau,                  &  ! Input, K
                                  solar_tau,                 &  ! Input, scalar

                                  layer_radiance,            &  ! Input, K
                                  downwelling_radiance,      &  ! Input, scalar
                                  upwelling_radiance,        &  ! Input, scalar

                                  ! -- Adjoint inputs
                                  layer_radiance_AD,         &  ! In/Output, K
                                  downwelling_radiance_AD,   &  ! In/Output, scalar
                                  upwelling_radiance_AD,     &  ! In/Output, scalar

                                  brightness_temperature_AD, &  ! In/Output, scalar

                                  ! -- Other inputs
                                  secant_solar_angle,        &  ! Input, scalar
                                  valid_solar,               &  ! Input, scalar
                                  channel_index,             &  ! Input, scalar

                                  ! -- Adjoint outputs
                                  temperature_AD,            &  ! In/Output, K

                                  surface_temperature_AD,    &  ! In/Output, scalar
                                  surface_emissivity_AD,     &  ! In/Output, scalar
                                  surface_reflectivity_AD,   &  ! In/Output, scalar

                                  tau_AD,                    &  ! In/Output, K
                                  flux_tau_AD,               &  ! In/Output, K
                                  solar_tau_AD               )  ! In/Output, scalar




    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Forward input
    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )     :: temperature

    REAL( fp_kind ),                 INTENT( IN )     :: surface_temperature
    REAL( fp_kind ),                 INTENT( IN )     :: surface_emissivity
    REAL( fp_kind ),                 INTENT( IN )     :: surface_reflectivity

    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )     :: tau
    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )     :: flux_tau
    REAL( fp_kind ),                 INTENT( IN )     :: solar_tau

    REAL( fp_kind ), DIMENSION( : ), INTENT( IN )     :: layer_radiance
    REAL( fp_kind ),                 INTENT( IN )     :: downwelling_radiance
    REAL( fp_kind ),                 INTENT( IN )     :: upwelling_radiance

    ! -- Adjoint inputs
    REAL( fp_kind ), DIMENSION( : ), INTENT( IN OUT ) :: layer_radiance_AD
    REAL( fp_kind ),                 INTENT( IN OUT ) :: downwelling_radiance_AD
    REAL( fp_kind ),                 INTENT( IN OUT ) :: upwelling_radiance_AD

    REAL( fp_kind ),                 INTENT( IN OUT ) :: brightness_temperature_AD

    ! -- Other input
    REAL( fp_kind ),                 INTENT( IN )     :: secant_solar_angle
    INTEGER,                         INTENT( IN )     :: valid_solar
    INTEGER,                         INTENT( IN )     :: channel_index

    ! -- Adjoint outputs
    REAL( fp_kind ), DIMENSION( : ), INTENT( IN OUT ) :: temperature_AD

    REAL( fp_kind ),                 INTENT( IN OUT ) :: surface_temperature_AD
    REAL( fp_kind ),                 INTENT( IN OUT ) :: surface_emissivity_AD
    REAL( fp_kind ),                 INTENT( IN OUT ) :: surface_reflectivity_AD

    REAL( fp_kind ), DIMENSION( : ), INTENT( IN OUT ) :: tau_AD
    REAL( fp_kind ), DIMENSION( : ), INTENT( IN OUT ) :: flux_tau_AD
    REAL( fp_kind ),                 INTENT( IN OUT ) :: solar_tau_AD


 

    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'COMPUTE_RADIANCE_AD'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: k, n_layers
    INTEGER :: l

    REAL( fp_kind ) :: surface_B
    REAL( fp_kind ) :: surface_B_AD

    REAL( fp_kind ) :: B_ur_AD
    REAL( fp_kind ) :: B_dr_AD


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC SIZE



    !#--------------------------------------------------------------------------#
    !#                   -- Determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    n_layers = SIZE( temperature )



    !#--------------------------------------------------------------------------#
    !#         -- Calculate the adjoint of the brightness temperature --        #
    !#--------------------------------------------------------------------------#

    ! -- Assign the channel index to a short name
    l = channel_index

    ! -- upwelling radiance adjoint
    CALL sensor_planck_temperature_AD( l,                         &  ! Input
                                       upwelling_radiance,        &  ! Input
                                       brightness_temperature_AD, &  ! Input
                                       upwelling_radiance_AD      )  ! In/Output
    brightness_temperature_AD = ZERO



    !#--------------------------------------------------------------------------#
    !#      -- Calculate the adjoints of the upwelling radiance term --         #
    !#--------------------------------------------------------------------------#

    ! --------------------------
    ! Top layer (closest to TOA)
    ! --------------------------

    ! -- Adjoint of top layer radiance
    layer_radiance_AD( 1 ) = layer_radiance_AD( 1 ) + ( ( ONE - tau( 1 ) ) * upwelling_radiance_AD )

    ! -- Adjoint of top layer dtau
    tau_AD( 1 ) = tau_AD( 1 ) + ( -layer_radiance( 1 ) * upwelling_radiance_AD )
    ! NOTE: No upwelling_radiance_AD = 0 here since
    !       upwelling_radiance_tl = upwelling_radiance_tl + (...)

    ! -- Adjoint of top layer temperature
    CALL sensor_planck_radiance_AD( l,                    &  ! Input
                                    temperature(1),       &  ! Input
                                    layer_radiance_AD(1), &  ! Input
                                    temperature_AD(1)     )  ! In/Output
    layer_radiance_AD( 1 ) = ZERO


    ! --------------------------------
    ! Loop over layers from TOA-1->SFC
    ! --------------------------------

    k_down_layer_loop: DO k = 2, n_layers

      ! -- Adjoint of layer radiance
      layer_radiance_AD( k ) = layer_radiance_AD( k ) + &
                               ( ( tau( k-1 ) - tau( k ) ) * upwelling_radiance_AD )

      ! -- Adjoint of dtau
      B_ur_AD = layer_radiance( k ) * upwelling_radiance_AD
      tau_AD( k )   = tau_AD(k)     - B_ur_AD
      tau_AD( k-1 ) = tau_AD( k-1 ) + B_ur_AD
      ! NOTE: No upwelling_radiance_AD = 0 here since
      !       upwelling_radiance_tl = upwelling_radiance_tl + (...)

      ! -- Adjoint of layer temperature
      CALL sensor_planck_radiance_AD( l,                      &  ! Input
                                      temperature( k ),       &  ! Input
                                      layer_radiance_AD( k ), &  ! Input
                                      temperature_AD( k )     )  ! In/Output
      layer_radiance_AD( k ) = ZERO
 
    END DO k_down_layer_loop



    !#--------------------------------------------------------------------------#
    !#     -- Calculate the adjoints of the surface tangent-linear terms --     #
    !#--------------------------------------------------------------------------#

    ! -- Recalculate surface Planck radiance
    CALL sensor_planck_radiance( l,                   &  ! Input
                                 surface_temperature, &  ! Input
                                 surface_B            )  ! Output


    ! -- Surface property adjoints
    surface_reflectivity_AD = surface_reflectivity_AD + &
                              ( tau( n_layers ) * downwelling_radiance * upwelling_radiance_AD )
 
    surface_emissivity_AD   = surface_emissivity_AD + &
                              ( tau( n_layers ) * surface_B * upwelling_radiance_AD )
 
    ! -- Total transmittance adjoint
    tau_AD( n_layers ) = tau_AD( n_layers ) + &
                         ( ( ( surface_emissivity   * surface_B            ) + &
                             ( surface_reflectivity * downwelling_radiance ) ) * upwelling_radiance_AD )
 
    ! -- Downweling radiance adjoint
    downwelling_radiance_AD = downwelling_radiance_AD + &
                              ( tau( n_layers ) * surface_reflectivity * upwelling_radiance_AD )
 
    ! -- Surface emission adjoint. This quantity
    ! -- is initialised each call.
    surface_B_AD = ( tau( n_layers ) * surface_emissivity * upwelling_radiance_AD )

    ! -- Set upwelling adjoint to zero.
    ! --  No more impact on gradient vector
    upwelling_radiance_AD = ZERO

    ! -- Surface temperature adjoint
    CALL sensor_planck_radiance_AD( l,                     & ! Input
                                    surface_temperature,   & ! Input
                                    surface_B_AD,          & ! Input
                                    surface_temperature_AD ) ! In/Output
    ! -- No need to zero surface_B_AD as it is initialised each call.
  


    !#--------------------------------------------------------------------------#
    !#      -- Calculate the adjoints of the solar transmittance term --        #
    !#--------------------------------------------------------------------------#

    solar_term: IF ( valid_solar == 1 ) THEN

      solar_tau_AD = solar_tau_AD + &
                     ( ( solar_irradiance( l ) / secant_solar_angle ) * downwelling_radiance_AD )

      ! NOTE: No downwelling_radiance_AD = 0 here since
      !       downwelling_radiance_TL = downwelling_radiance_TL + (TL solar term)

    END IF solar_term




    !#--------------------------------------------------------------------------#
    !#        -- Calculate the adjoints of the downwelling flux term --         #
    !#--------------------------------------------------------------------------#

    ! -----------------------------
    ! Bottom layer (closest to SFC)
    ! -----------------------------

    ! -- Adjoint of bottom layer radiance
    layer_radiance_AD( n_layers ) = layer_radiance_AD( n_layers ) + &
                                    ( ( ONE - flux_tau( n_layers ) ) * downwelling_radiance_AD )

    ! -- Adjoint of flux transmittance
    flux_tau_AD( n_layers ) = flux_tau_AD( n_layers ) - &
                              ( layer_radiance( n_layers ) * downwelling_radiance_AD )
    ! NOTE: No downwelling_radiance_AD = 0 here since
    !       downwelling_radiance_tl = downwelling_radiance_tl + (TL flux term)
 
    ! -- Adjoint of layer temperature
    CALL sensor_planck_radiance_AD( l,                             &  ! Input
                                    temperature( n_layers ),       &  ! Input
                                    layer_radiance_AD( n_layers ), &  ! Input
                                    temperature_AD( n_layers )     )  ! In/Output
    layer_radiance_AD( n_layers ) = ZERO


    ! --------------------------------
    ! Loop over layers from SFC-1->TOA
    ! --------------------------------

    k_up_layer_loop: DO k = n_layers - 1, 1, -1

      ! -- Adjoint of layer radiance
      layer_radiance_AD( k ) = layer_radiance_AD( k ) + &
                               ( ( flux_tau( k+1 ) - flux_tau( k ) ) * downwelling_radiance_AD )

      ! -- Adjoint of dtau
      B_dr_AD = layer_radiance( k ) * downwelling_radiance_AD
      flux_tau_AD( k )   = flux_tau_AD( k )   - B_dr_AD
      flux_tau_AD( k+1 ) = flux_tau_AD( k+1 ) + B_dr_AD
      ! NOTE: No downwelling_radiance_AD = 0 here since
      !       downwelling_radiance_tl = downwelling_radiance_tl + (...)

      ! -- Adjoint of layer temperature
      CALL sensor_planck_radiance_AD( l,                 &  ! Input
                                 temperature( k ),       &  ! Input
                                 layer_radiance_AD( k ), &  ! Input
                                 temperature_AD( k )     )  ! In/Output
      layer_radiance_AD( k ) = ZERO

    END DO k_up_layer_loop


    ! --------------------------------------------
    ! Background term. Note that the emissivity of
    ! space is implicitly assumed to 1.0.
    ! --------------------------------------------

    flux_tau_AD( 1 ) = flux_tau_AD( 1 ) + &
                       ( cosmic_background_radiance( l ) * downwelling_radiance_AD )

    downwelling_radiance_AD = ZERO

  END SUBROUTINE compute_radiance_AD

END MODULE radiance


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
! Revision 2.4  2001/08/16 17:11:57  paulv
! - Updated documentation
!
! Revision 2.3  2001/08/08 20:02:12  paulv
! - Removed sensor Planck function routines and placed them in their own
!   module, SENSOR_PLANCK_ROUTINES. Some routines were required for other
!   uses so their PRIVATE subprogram status wasn't amenable to code-sharing.
! - Updated header documentation.
!
! Revision 2.2  2001/08/01 16:58:32  paulv
! - Corrected bug in COMPUTE_RADIANCE() function. The initialisation
!   statement of the downwelling radiance was,
!      downwelling_radiance = cosmic_background_radiance( l )
!   and has been changed to,
!      downwelling_radiance = cosmic_background_radiance( l ) * flux_tau( 1 )
!   i.e. the transmission of the space emission term to the surface. This
!   was a holdover from earlier versions of the functions when the transmittances
!   were calculated and passed as *layer* rather than layer-to-surface
!   transmittances.
! - Removed initialisation and zeroing of the adjoint of the surface emissioni
!   term SURFACE_B_AD. This is used in only one place so there is no need to
!   do,
!     surface_B_AD = ZERO
!     surface_B_AD = surface_B_AD + &
!                    ( tau( n_layers ) * surface_emissivity * upwelling_radiance_AD )
!     ....use surface_B_AD...
!     surface_B_AD = ZERO
!   when,
!     surface_B_AD = ( tau( n_layers ) * surface_emissivity * upwelling_radiance_AD )
!   will do.
! - Updated documentation.
!
! Revision 2.1  2001/05/29 18:05:29  paulv
! - All tangent-linear and adjoint routines included.
! - No more optional arguments of downwelling flux and surface reflectivity -
!   they are expected.
! - Altered the method of calculating the layer contributions. Changed code
!   from:
!     layer_radiance(k) = (1-tau)*B(T) + tau*layer_radiance(k-1)
!   to:
!     layer_radiance(k) = B(T) * dtau
!
! Revision 1.5  2001/01/24 20:14:21  paulv
! - Latest test versions.
!
! Revision 1.4  2000/11/09 20:46:07  paulv
! - Added solar term.
! - Downwelling flux transmittance term is now an optional argument.
!   If not specified, the layer radiances are calculated during the
!   upwelling radiance integration.
! - Surface reflectivity is an optional argument. If not specified the
!   surface emissivity is used to generate an isotropic reflectivity.
!
! Revision 1.3  2000/08/31 19:36:33  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 1.2  2000/08/24 15:48:34  paulv
! - Replaced "regular" reflectivity for reflected downwelling thermal with
!   the isotropic reflectivity in the COMPUTE_RADIANCE subprogram.
! - Updated module and subprogram documentation.
!
! Revision 1.1  2000/08/21 20:59:34  paulv
! Initial checkin.
!
!
!
!
