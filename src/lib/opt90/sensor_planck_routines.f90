!------------------------------------------------------------------------------
!M+
! NAME:
!       sensor_planck_routines
!
! PURPOSE:
!       Module containing the sensor Planck function routines.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE sensor_planck_routines
!
! OUTPUTS:
!       None.
!
! MODULES:
!       type_kinds:                    Module containing data type kind definitions.
!
!       parameters:                    Module containing parameter definitions for
!                                      the RT model.
!
!       spectral_coefficients:         Module containing the RT model spectral
!                                      coefficients
!
! CONTAINS:
!       sensor_planck_radiance:        PUBLIC subroutine to calculate the instrument
!                                      channel radiance.
!
!       sensor_planck_radiance_TL:     PUBLIC subroutine to calculate the tangent-linear
!                                      instrument channel radiance.
!
!       sensor_planck_radiance_AD:     PUBLIC subroutine to calculate the adjoint
!                                      instrument channel radiance.
!
!       sensor_planck_temperature:     PUBLIC subroutine to calculate the instrument
!                                      channel brightness temperature.
!
!       sensor_planck_temperature_TL:  PUBLIC subroutine to calculate the tangent-linear
!                                      instrument channel brightness temperature.
!
!       sensor_planck_temperature_AD:  PUBLIC subroutine to calculate the adjoint
!                                      instrument channel brightness temperature.
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! RESTRICTIONS:
!       These functions are called frequently so no input checking is
!       performed.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 08-Aug-2001
!                       paul.vandelst@ssec.wisc.edu
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

MODULE sensor_planck_routines


  ! ---------------------
  ! Module use statements
  ! ---------------------

  USE type_kinds, ONLY : fp_kind
  USE parameters
  USE spectral_coefficients


  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------

  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------

  PRIVATE

  PUBLIC  :: sensor_planck_radiance
  PUBLIC  :: sensor_planck_radiance_TL
  PUBLIC  :: sensor_planck_radiance_AD

  PUBLIC  :: sensor_planck_temperature
  PUBLIC  :: sensor_planck_temperature_TL
  PUBLIC  :: sensor_planck_temperature_AD


CONTAINS


!--------------------------------------------------------------------------------
!P+
! NAME:
!       sensor_planck_radiance
!
! PURPOSE:
!       Subroutine to calculate the instrument channel radiance.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL sensor_planck_radiance( channel,     &  ! Input
!                                    temperature, &  ! Input
!                                    radiance     )  ! Output
!
! INPUT ARGUMENTS:
!       channel:     Channel index id. This is a unique index
!                    to a (supported) sensor channel.
!                    UNITS:      None
!                    TYPE:       Integer
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: INTENT( IN )
!
!       temperature: Temperature for which the Planck radiance is
!                    to be calculated.
!                    UNITS:      Kelvin, K
!                    TYPE:       Real
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       radiance:    Channel Planck radiance.
!                    UNITS:      mW/(m^2.sr.cm^-1)
!                    TYPE:       Real
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: INTENT( OUT )
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
!       None.
!
! RESTRICTIONS:
!       Spectral coefficients are obtained from the SPECTRAL_COEFFICIENTS module
!       so only radiances for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! PROCEDURE:
!       First a polychromatic correction is applied to give an effective
!       temperature,
!
!         T_eff = bc1 + ( bc2 * T )
!
!       The sensor radiance is then calculated using the effective temperature:
!
!                       pc1
!         R = ------------------------
!              EXP( pc2 / T_eff ) - 1
!
!       The bc1, bc2, pc1, and pc2 values are obtained from the 
!       SPECTRAL_COEFFICIENTS module which is filled during the initialisation
!       phase.
!P-
!--------------------------------------------------------------------------------

  SUBROUTINE sensor_planck_radiance( channel,     &  ! Input
                                     temperature, &  ! Input
                                     radiance     )  ! Output

    ! -- Arguments
    INTEGER,         INTENT( IN )  :: channel
    REAL( fp_kind ), INTENT( IN )  :: temperature
    REAL( fp_kind ), INTENT( OUT ) :: radiance

    ! -- Local
    REAL( fp_kind ) :: effective_temperature

    INTRINSIC EXP


    ! -------------------------------------
    ! Apply the polychromaticity correction
    ! to obtain an effective temperature
    ! -------------------------------------

    effective_temperature = band_c1( channel ) + ( band_c2( channel ) * temperature )


    ! -----------------------------
    ! Calculate the Planck radiance
    ! -----------------------------

    radiance =                  planck_c1( channel )  / &
    !          -------------------------------------------------------------
               ( EXP( planck_c2( channel ) / effective_temperature ) - ONE )

  END SUBROUTINE sensor_planck_radiance




!--------------------------------------------------------------------------------
!P+
! NAME:
!       sensor_planck_radiance_TL
!
! PURPOSE:
!       Subroutine to calculate the tangent-linear instrument channel radiance.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL sensor_planck_radiance_TL( channel,        &  ! Input
!                                       temperature,    &  ! Input
!                                       temperature_TL, &  ! Input
!                                       radiance_TL     )  ! Output
!
! INPUT ARGUMENTS:
!       channel:        Channel index id. This is a unique index
!                       to a (supported) sensor channel.
!                       UNITS:      None
!                       TYPE:       Integer
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
!       temperature:    Temperature for which the tangent-linear Planck radiance
!                       is to be calculated.
!                       UNITS:      Kelvin, K
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
!       temperature_TL: Tangent-linear temperature for which the tangent-linear
!                       Planck radiance is required.
!                       UNITS:      Kelvin, K
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       radiance_TL:    Tangent-linear Planck radiance.
!                       UNITS:      mW/(m^2.sr.cm^-1)
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( OUT )
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
!       None.
!
! RESTRICTIONS:
!       Spectral coefficients are obtained from the SPECTRAL_COEFFICIENTS module
!       so only radiances for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! PROCEDURE:
!       First a polychromatic correction is applied to give an effective
!       temperature,
!
!         T_eff = bc1 + ( bc2 . T )
!
!       The sensor tangent-linear radiance is then calculated by first computing
!       the exponent term,
!
!          exponent = EXP( pc2 / T_eff )
!
!       and then the actual operator,
!
!                 pc1 . pc2 . bc1 . exponent
!         F = ---------------------------------
!              ( T_eff . ( exponent - 1 ) )^2
!
!       which is the derivate of the Planck equation wrt temperature. The
!       tangent-linear radiance is then determined by,
!
!         dR = F . dT
!
!       where dT is the input tangent-linear temeprature.
!
!       The bc1, bc2, pc1, and pc2 values are obtained from the 
!       SPECTRAL_COEFFICIENTS module which is filled during the initialisation
!       phase.
!P-
!--------------------------------------------------------------------------------

  SUBROUTINE sensor_planck_radiance_TL( channel,        &  ! Input
                                        temperature,    &  ! Input
                                        temperature_TL, &  ! Input
                                        radiance_TL     )  ! Output

    ! -- Arguments
    INTEGER,         INTENT( IN )  :: channel
    REAL( fp_kind ), INTENT( IN )  :: temperature
    REAL( fp_kind ), INTENT( IN )  :: temperature_TL
    REAL( fp_kind ), INTENT( OUT ) :: radiance_TL

    ! -- Local
    REAL( fp_kind ) :: effective_temperature
    REAL( fp_kind ) :: exponent
    REAL( fp_kind ) :: F

    INTRINSIC EXP


    ! -------------------------------------
    ! Apply the polychromaticity correction
    ! -------------------------------------

    effective_temperature = band_c1( channel ) + ( band_c2( channel ) * temperature )


    ! --------------------------------------
    ! Calculate the Planck function operator
    ! --------------------------------------

    ! -- The exponent term
    exponent = EXP( planck_c2( channel ) / effective_temperature )

    ! -- The operator, call it F
    F =  planck_c1( channel ) * planck_c2( channel ) * exponent * band_c2( channel ) / &
    !   -----------------------------------------------------------------------------
                      ( effective_temperature * ( exponent - ONE ) )**2


    ! -------------------------------------
    ! Calculate the tangent-linear radiance
    ! -------------------------------------

    radiance_TL = F * temperature_TL

  END SUBROUTINE sensor_planck_radiance_TL




!--------------------------------------------------------------------------------
!P+
! NAME:
!       sensor_planck_radiance_AD
!
! PURPOSE:
!       Subroutine to calculate the adjoint instrument channel radiance.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL sensor_planck_radiance_AD( channel,       &  ! Input
!                                       temperature,   &  ! Input
!                                       radiance_AD,   &  ! Input
!                                       temperature_AD )  ! In/Output
!
! INPUT ARGUMENTS:
!       channel:        Channel index id. This is a unique index
!                       to a (supported) sensor channel.
!                       UNITS:      None
!                       TYPE:       Integer
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
!       temperature:    Temperature for which the tangent-linear Planck radiance
!                       is to be calculated.
!                       UNITS:      Kelvin
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
!       radiance_AD:    Adjoint Planck radiance.
!                       UNITS:      mW/(m2.sr.cm-1)
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       temperature_AD: Adjoint Planck temperature
!                       UNITS:      Kelvin
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN OUT )
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
!       The input adjoint radiance argument, radiance_AD, is NOT set to zero
!       before returning to the calling routine.
!
! RESTRICTIONS:
!       Spectral coefficients are obtained from the SPECTRAL_COEFFICIENTS module
!       so only radiances for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! PROCEDURE:
!       First a polychromatic correction is applied to give an effective
!       temperature,
!
!         T_eff = bc1 + ( bc2 . T )
!
!       The sensor tangent-linear radiance is then calculated by first computing
!       the exponent term,
!
!          exponent = EXP( pc2 / T_eff )
!
!       and then the actual operator,
!
!                 pc1 . pc2 . bc1 . exponent
!         F = ---------------------------------
!              ( T_eff . ( exponent - 1 ) )^2
!
!       which is the derivate of the Planck equation wrt temperature. The
!       adjoint temperature is then determined from,
!
!         T_AD = T_AD + ( F . R_AD )
!
!       where T_AD and R_AD on the LHS are the input adjoint temperature and
!       radiance respectively.
!
!       The bc1, bc2, pc1, and pc2 values are obtained from the 
!       SPECTRAL_COEFFICIENTS module which is filled during the initialisation
!       phase.
!P-
!--------------------------------------------------------------------------------

  SUBROUTINE sensor_planck_radiance_AD( channel,       &  ! Input
                                        temperature,   &  ! Input
                                        radiance_AD,   &  ! Input
                                        temperature_AD )  ! In/Output

    ! -- Arguments
    INTEGER,         INTENT( IN )     :: channel
    REAL( fp_kind ), INTENT( IN )     :: temperature
    REAL( fp_kind ), INTENT( IN )     :: radiance_AD
    REAL( fp_kind ), INTENT( IN OUT ) :: temperature_AD

    ! -- Local
    REAL( fp_kind ) :: effective_temperature
    REAL( fp_kind ) :: exponent
    REAL( fp_kind ) :: F

    INTRINSIC EXP


    ! -------------------------------------
    ! Apply the polychromaticity correction
    ! -------------------------------------

    effective_temperature = band_c1( channel ) + ( band_c2( channel ) * temperature )


    ! --------------------------------------
    ! Calculate the Planck function operator
    ! --------------------------------------

    ! -- The exponent term
    exponent = EXP( planck_c2( channel ) / effective_temperature )

    ! -- The operator, call it F
    F =  planck_c1( channel ) * planck_c2( channel ) * exponent * band_c2( channel ) / &
    !   -----------------------------------------------------------------------------
                      ( effective_temperature * ( exponent - ONE ) )**2


    ! ---------------------------------
    ! Calculate the adjoint temperature
    ! ---------------------------------

    temperature_AD = temperature_AD + ( F * radiance_AD )

  END SUBROUTINE sensor_planck_radiance_AD




!--------------------------------------------------------------------------------
!P+
! NAME:
!       sensor_planck_temperature
!
! PURPOSE:
!       Subroutine to calculate the instrument channel brightness temperature.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL sensor_planck_temperature( channel,    &  ! Input
!                                       radiance,   &  ! Input
!                                       temperature )  ! Output
!
! INPUT ARGUMENTS:
!       channel:     Channel index id. This is a unique index
!                    to a (supported) sensor channel.
!                    UNITS:      None
!                    TYPE:       Integer
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: INTENT( IN )
!
!       radiance:    Radiance for which the Planck temperature is desired.
!                    UNITS:      mW/(m^2.sr.cm^-1)
!                    TYPE:       Real
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       temperature: Planck temperature.
!                    UNITS:      Kelvin, K
!                    TYPE:       Real
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: INTENT( IN )
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
!       None.
!
! RESTRICTIONS:
!       Spectral coefficients are obtained from the SPECTRAL_COEFFICIENTS module
!       so only temperatures for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! PROCEDURE:
!       First the effective temperature is calculated from the inverse Planck function,
!
!                        pc2
!         T_eff = ------------------
!                  LOG( pc1/R + 1 )
!
!       The polychromatic correction is then removed to provide the brightness
!       temperature,
!
!              T_eff - bc1
!         T = -------------
!                  bc2
!
!       The bc1, bc2, pc1, and pc2 values are obtained from the 
!       SPECTRAL_COEFFICIENTS module which is filled during the initialisation
!       phase.
!P-
!--------------------------------------------------------------------------------

  SUBROUTINE sensor_planck_temperature( channel,    &  ! Input
                                        radiance,   &  ! Input
                                        temperature )  ! Output

    ! -- Arguments
    INTEGER,         INTENT( IN )  :: channel
    REAL( fp_kind ), INTENT( IN )  :: radiance
    REAL( fp_kind ), INTENT( OUT ) :: temperature

    ! -- Local
    REAL( fp_kind ) :: effective_temperature

    INTRINSIC LOG


    ! -----------------------------------
    ! Calculate the effective temperature
    ! -----------------------------------

    effective_temperature =              planck_c2( channel )  / &
    !                       ------------------------------------------------
                            LOG( ( planck_c1( channel ) / radiance ) + ONE )

    ! -------------------------------------
    ! Apply the polychromatic correction to 
    ! obtain the true temperature
    ! -------------------------------------

    temperature = ( effective_temperature - band_c1( channel ) ) / &
    !             ----------------------------------------------
                                band_c2( channel )

  END SUBROUTINE sensor_planck_temperature




!--------------------------------------------------------------------------------
!P+
! NAME:
!       sensor_planck_temperature_TL
!
! PURPOSE:
!       Subroutine to calculate the tangent-linear instrument channel
!       brightness temperature.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL sensor_planck_temperature_TL( channel,       &  ! Input
!                                          radiance,      &  ! Input
!                                          radiance_TL,   &  ! Input
!                                          temperature_TL )  ! Output
!
! INPUT ARGUMENTS:
!       channel:        Channel index id. This is a unique index
!                       to a (supported) sensor channel.
!                       UNITS:      None
!                       TYPE:       Integer
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
!       radiance:       Radiance at which the tangent-linear Planck temperature
!                       is desired.
!                       UNITS:      mW/(m^2.sr.cm^-1)
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
!       radiance_TL:    Tangent-linear radiance for which the tangent-linear
!                       Planck temperature is desired.
!                       UNITS:      mW/(m^2.sr.cm^-1)
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       temperature_TL: Tangent-linear Planck temperature.
!                       UNITS:      Kelvin, K
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
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
!       None.
!
! RESTRICTIONS:
!       Spectral coefficients are obtained from the SPECTRAL_COEFFICIENTS module
!       so only temperatures for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! PROCEDURE:
!       First the logarithm argument is calculated,
!
!         a = pc1/R + 1
!
!       The inverse Planck function operator is then calculated,
!
!                      pc1 . pc2
!         F = ------------------------------
!              bc2 . a . ( R . LOG( a ) )^2
!
!       and the tangent-linear temperature is then given by,
!
!         dT = F . dR
!
!       The bc1, bc2, pc1, and pc2 values are obtained from the 
!       SPECTRAL_COEFFICIENTS module which is filled during the initialisation
!       phase.
!P-
!--------------------------------------------------------------------------------


  SUBROUTINE sensor_planck_temperature_TL( channel,       &  ! Input
                                           radiance,      &  ! Input
                                           radiance_TL,   &  ! Input
                                           temperature_TL )  ! Output

    ! -- Arguments
    INTEGER,         INTENT( IN )  :: channel
    REAL( fp_kind ), INTENT( IN )  :: radiance
    REAL( fp_kind ), INTENT( IN )  :: radiance_TL
    REAL( fp_kind ), INTENT( OUT ) :: temperature_TL

    ! -- Local
    REAL( fp_kind ) :: argument
    REAL( fp_kind ) :: F

    INTRINSIC LOG


    ! --------------------------------------
    ! Calculate the Planck function operator
    ! --------------------------------------

    ! -- The logarithm argument
    argument = ( planck_c1( channel ) / radiance ) + ONE

    ! -- The operator, call it F
    F =             planck_c1( channel ) * planck_c2( channel ) / &
    !   -----------------------------------------------------------------------
         ( band_c2( channel ) * argument * ( radiance * LOG( argument ) )**2 )


    ! ----------------------------------------
    ! Calculate the tangent-linear temperature
    ! ----------------------------------------

    temperature_TL = F * radiance_TL

  END SUBROUTINE sensor_planck_temperature_TL



!--------------------------------------------------------------------------------
!P+
! NAME:
!       sensor_planck_temperature_AD
!
! PURPOSE:
!       Subroutine to calculate the adjoint instrument channel
!       brightness temperature.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       CALL sensor_planck_temperature_AD( channel,        &  ! Input
!                                          radiance,       &  ! Input
!                                          temperature_AD, &  ! Input
!                                          radiance_AD     )  ! In/Output
!
! INPUT ARGUMENTS:
!       channel:        Channel index id. This is a unique index
!                       to a (supported) sensor channel.
!                       UNITS:      None
!                       TYPE:       Integer
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
!       radiance:       Radiance at which the adjoint radiance is desired.
!                       UNITS:      mW/(m^2.sr.cm^-1)
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
!       temperature_AD: Adjoint Planck temperature.
!                       UNITS:      Kelvin, K
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       radiance_AD:    Adjoint radiance.
!                       UNITS:      K.m^2.sr.cm^-1/mW
!                       TYPE:       Real
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN OUT )
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
!       The input adjoint temperature argument, temperature_AD, is NOT set to zero
!       before returning to the calling routine.
!
! RESTRICTIONS:
!       Spectral coefficients are obtained from the SPECTRAL_COEFFICIENTS module
!       so only temperatures for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! PROCEDURE:
!       First the logarithm argument is calculated,
!
!         a = pc1/R + 1
!
!       The inverse Planck function operator is then calculated,
!
!                      pc1 . pc2
!         F = ------------------------------
!              bc2 . a . ( R . LOG( a ) )^2
!
!       which is the derivate of the Planck temperature wrt radiance. The
!       adjoint radiance is then determined from,
!
!         R_AD = R_AD + ( F . T_AD )
!
!       where R_AD and T_AD on the LHS are the input adjoint radiance and
!       temperature respectively.
!
!       The bc1, bc2, pc1, and pc2 values are obtained from the 
!       SPECTRAL_COEFFICIENTS module which is filled during the initialisation
!       phase.
!P-
!--------------------------------------------------------------------------------


  SUBROUTINE sensor_planck_temperature_AD( channel,        &  ! Input
                                           radiance,       &  ! Input
                                           temperature_AD, &  ! Input
                                           radiance_AD     )  ! In/Output

    ! -- Arguments
    INTEGER,         INTENT( IN )     :: channel
    REAL( fp_kind ), INTENT( IN )     :: radiance
    REAL( fp_kind ), INTENT( IN )     :: temperature_AD
    REAL( fp_kind ), INTENT( IN OUT ) :: radiance_AD

    ! -- Local
    REAL( fp_kind ) :: argument
    REAL( fp_kind ) :: F

    INTRINSIC LOG


    ! --------------------------------------
    ! Calculate the Planck function operator
    ! --------------------------------------

    ! -- The logarithm argument
    argument = ( planck_c1( channel ) / radiance ) + ONE

    ! -- The operator, call it F
    F =             planck_c1( channel ) * planck_c2( channel ) / &
    !   -----------------------------------------------------------------------
         ( band_c2( channel ) * argument * ( radiance * LOG( argument ) )**2 )


    ! ------------------------------
    ! Calculate the adjoint radiance
    ! ------------------------------

    radiance_AD = radiance_AD + ( F * temperature_AD )

  END SUBROUTINE sensor_planck_temperature_AD

END MODULE sensor_planck_routines


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
! Revision 1.1  2001/08/08 20:04:03  paulv
! Initial checkin.
! - Routines were extracted from the RADIANCE module and placed in their
!   own module to facilitate code-sharing.
!
!
!
