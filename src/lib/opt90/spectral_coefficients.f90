!------------------------------------------------------------------------------
!M+
! NAME:
!       spectral_coefficients
!
! PURPOSE:
!       Module to hold the RT model spectral coefficients and their access
!       routines. 
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE spectral_coefficients
!
! OUTPUTS:
!       frequency:                     Channel frequency
!                                      UNITS:      GHz
!                                      TYPE:       Double
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       wavenumber:                    Channel frequency
!                                      UNITS:      cm^-1
!                                      TYPE:       Double
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       planck_c1:                     First Planck function coefficient
!                                      UNITS:      mW/(m^2.sr.cm^-2)
!                                      TYPE:       Double
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       planck_c2:                     Second Planck function coefficient
!                                      UNITS:      K/cm^1
!                                      TYPE:       Double
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       band_c1:                       Polychromatic band correction offset.
!                                      UNITS:      K
!                                      TYPE:       Double
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       band_c2:                       Polychromatic band correction slope.
!                                      UNITS:      K/K
!                                      TYPE:       Double
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       is_microwave_channel:          Flag indicating if the particular satellite
!                                      channel is a microwave channel.
!                                      If = 0, IR channel
!                                         = 1, uW channel
!                                      UNITS:      None
!                                      TYPE:       Long
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       cosmic_background_temperature: Cosmic background temperature in Kelvin to
!                                      use in the radiative transfer. These values
!                                      are frequency dependent for microwave channels
!                                      to account for instrument calibration using the
!                                      Rayleigh-Jeans approximation, but the radiative
!                                      transfer using the full Planck function
!                                      expression. IR channels use a fixed value of
!                                      2.736K.
!                                      UNITS:      K
!                                      TYPE:       Double
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       cosmic_background_radiance:    Cosmic background radiance used to initialise
!                                      the radiative transfer. These values
!                                      are frequency dependent for microwave channels
!                                      to account for instrument calibration using the
!                                      Rayleigh-Jeans approximation, but the radiative
!                                      transfer using the full Planck function
!                                      expression. IR channel values are all set to 0.0.
!                                      UNITS:      mW/(m^2.sr.cm^-1)
!                                      TYPE:       Double
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       solar_irradiance:              Exterrestrial solar irradiance source values
!                                      (based on Kurucz).
!                                      UNITS:      mW/(m^2.cm^-1)
!                                      TYPE:       Double
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       blackbody_irradiance:          Exterrestrial solar irradiance source values
!                                      for a blackbody source at the equivalent solar
!                                      temperature of 5783K.
!                                      UNITS:      mW/(m^2.cm^-1)
!                                      TYPE:       Double
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
!       is_solar_channel:              Array of integer flag values indicating if the
!                                      particular satellite channel possesses valid
!                                      solar irradiance values.
!                                      If = 0, No
!                                         = 1, Yes
!                                      UNITS:      None
!                                      TYPE:       Long
!                                      DIMENSION:  L
!                                      ATTRIBUTES: PUBLIC, SAVE
!
! MODULES:
!       type_kinds:           Module containing data type kind definitions.
!
!       error_handler:        Module to define error codes and handle error
!                             conditions
!
!       parameters:           Module containing parameter definitions for the
!                             RT model.
!
!       coefficient_utility:  Module containing coefficient file utility subprograms.
!
! CONTAINS:
!       read_spectral_coefficients:     PUBLIC function to read the spectral coefficients
!                                       for the satellite/channels used by the RT model
!                                       and fill the PUBLIC data arrays.
!
!       destroy_spectral_coefficients:  PUBLIC function to release the memory that
!                                       was allocated for the spectral coefficient data.
!
!       write_spectral_coefficients:    PUBLIC function to write the new format spectral
!                                       coefficient data file.
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
!       Written by:     Paul van Delst, CIMSS@NOAA/NCEP 12-Jun-2000
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

MODULE spectral_coefficients


  ! ---------------------
  ! Module use statements
  ! ---------------------

  USE type_kinds
  USE file_utility
  USE error_handler
  USE parameters
  USE coefficient_utility


  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------

  IMPLICIT NONE


  ! ------------------
  ! Default visibilities
  ! ------------------

  PRIVATE
  PUBLIC :: read_spectral_coefficients
  PUBLIC :: destroy_spectral_coefficients
  PUBLIC :: write_spectral_coefficients



  ! -----------------
  ! Module parameters
  ! -----------------

  ! -- Spectral coefficient file version numbers
  INTEGER( Long ), PARAMETER, PRIVATE :: MIN_VALID_RELEASE = 1_Long
  INTEGER( Long ), PARAMETER, PRIVATE :: MAX_VALID_RELEASE = 1_Long

  INTEGER( Long ), PARAMETER, PRIVATE :: MIN_VALID_VERSION = 1_Long
  INTEGER( Long ), PARAMETER, PRIVATE :: MAX_VALID_VERSION = 1_Long

  ! -- Spectral coefficient file descriptors:
  ! -- Number of items in the spectral coefficient file
  INTEGER( Long ), PARAMETER, PRIVATE :: N_SPECTRAL_ITEMS = 12_Long

  ! -- Data types of the spectral coefficient data
  !    5 = Double (i.e. 8-byte float)
  !    4 = Single (i.e. 4-byte float)
  !    3 = Long   (i.e. 4-byte integer)
  INTEGER( Long ), PARAMETER, PRIVATE, &
                   DIMENSION( N_SPECTRAL_ITEMS ) :: SPECTRAL_DATA_TYPE = &
                                                    (/ 5_Long, 5_Long, 5_Long, 5_Long, &
                                                       5_Long, 5_Long, 3_Long, 5_Long, &
                                                       5_Long, 5_Long, 5_Long, 3_Long /)

  ! -- Names of the data items (for error processing)
  CHARACTER( * ), PARAMETER, PRIVATE, &
                  DIMENSION( N_SPECTRAL_ITEMS ) :: SPECTRAL_DATA_NAME = &
                                                     (/ 'FREQUENCY                    ', &
                                                        'WAVENUMBER                   ', &
                                                        'PLANCK_C1                    ', &
                                                        'PLANCK_C2                    ', &
                                                        'BAND_C1                      ', &
                                                        'BAND_C2                      ', &
                                                        'IS_MICROWAVE_CHANNEL         ', &
                                                        'COSMIC_BACKGROUND_TEMPERATURE', &
                                                        'COSMIC_BACKGROUND_RADIANCE   ', &
                                                        'SOLAR_IRRADIANCE             ', &
                                                        'BLACKBODY_IRRADIANCE         ', &
                                                        'IS_SOLAR_CHANNEL             ' /)

  ! ----------------
  ! Module variables
  ! ----------------

  CHARACTER( 128 ) :: message


  ! ---------------------------------------------------
  ! Definitions of shared data.
  !
  ! Note that the SAVE attribute is specified to ensure
  ! that the data is retained even when this module is
  ! not being directly accessed.
  ! ---------------------------------------------------

  ! -- Channel frequency in GHz (MW only, 0.0 for IR)
  REAL( Double ),  SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: frequency

  ! -- Channel frequency in cm^-1
  REAL( Double ),  SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: wavenumber

  ! -- First channel Planck function coefficient in mW/(m^2.sr.cm^-2)
  REAL( Double ),  SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: planck_c1

  ! -- Second channel Planck function coefficient in K/cm^1
  REAL( Double ),  SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: planck_c2

  ! -- Band correction (for polychromaticity) offset in K
  REAL( Double ),  SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: band_c1

  ! -- Band correction (for polychromaticity) slope in K/K
  REAL( Double ),  SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: band_c2

  ! -- Flag indicating whether channel is microwave or infrared.
  INTEGER( Long ), SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: is_microwave_channel

  ! -- Effective cosmic background temperature in Kelvin.
  REAL( Double ),  SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: cosmic_background_temperature

  ! -- Cosmic background radiance in mW/(m^2.sr.cm^-1).
  REAL( Double ),  SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: cosmic_background_radiance

  ! -- Kurucz solar irradiance source function in mW/(m^2.cm^-1)
  REAL( Double ),  SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: solar_irradiance

  ! -- Equivalent blackbody solar irradiance source function in mW/(m^2.cm^-1)
  REAL( Double ),  SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: blackbody_irradiance

  ! -- Flag indicating whether channel is sensitive to solar contribution.
  INTEGER( Long ), SAVE, PUBLIC, ALLOCATABLE, DIMENSION( : ) :: is_solar_channel



CONTAINS


!------------------------------------------------------------------------------
!S+
! NAME:
!       read_spectral_coefficients
!
! PURPOSE:
!       PUBLIC function to read the spectral coefficients for the satellite/
!       channels used by the RT model and fill the PUBLIC data arrays..
!
! CALLING SEQUENCE:
!       result = read_spectral_coefficients( coefficient_file, &
!                                            message_log  = message_log )
!
! INPUT ARGUMENTS:
!       coefficient_file: Name of the file containing the spectral coefficient data.
!                         UNITS:      None
!                         TYPE:       Character
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN )
!
!
! OPTIONAL INPUT ARGUMENTS:
!       message_log:      Character string specifying a filename in which any
!                         messages will be logged. If not specified, or if an
!                         error occurs opening the log file, the default action
!                         is to output messages to the screen.
!                         UNITS:      None
!                         TYPE:       Character
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       None.
!
! OPTIONAL OUTUPT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result = SUCCESS => coefficient read was successful
!              = FAILURE => error occurred opening or accessing coefficient
!                           data file
!
! CALLS:
!      display_message:         Subroutine to output messages
!                               SOURCE: error_handler module
!
!      open_coefficient_file:   Function to open the transmittance coefficient
!                               data file.
!                               SOURCE: coefficient_utility module
!
!      get_max_n_channels:      Routine to retrieve the value of the
!                               MAX_N_CHANNELS "pseudo-parameter".
!                               SOURCE: parameters module
!
!      set_max_n_channels:      Routine to set the value of the
!                               MAX_N_CHANNELS "pseudo-parameter".
!                               SOURCE: parameters module
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
!       All spectral coefficients are read in and stored by CHANNEL.
!S-
!------------------------------------------------------------------------------

  FUNCTION read_spectral_coefficients( coefficient_file, &
                                       message_log ) &
                                     RESULT ( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    CHARACTER( * ), INTENT( IN )           :: coefficient_file
    CHARACTER( * ), INTENT( IN ), OPTIONAL :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'READ_SPECTRAL_COEFFICIENTS'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- These are dimensioned ( Long ) as they are read in
    INTEGER( Long ) :: file_release
    INTEGER( Long ) :: file_version

    INTEGER( Long ) :: n_channels
    INTEGER( Long ) :: n_items
    INTEGER( Long ), DIMENSION( N_SPECTRAL_ITEMS ) :: data_type

    INTEGER :: l
    INTEGER :: io_status
    INTEGER :: file_id
    INTEGER :: allocate_status

    ! -- Maximum channels pseudo parameter
    INTEGER :: MAX_N_CHANNELS
    LOGICAL :: is_set


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC TRIM                          


    !#--------------------------------------------------------------------------#
    !#              -- Open the spectral coefficient data file --               #
    !#--------------------------------------------------------------------------#

    error_status = open_coefficient_file( coefficient_file, &
                                          file_id, &
                                          message_log = message_log )

    IF ( error_status /= SUCCESS ) THEN
      CALL display_message( ROUTINE_NAME, &
                            'Error opening '//TRIM( coefficient_file ), &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#                      -- Read the "file header" --                        #
    !#--------------------------------------------------------------------------#

    ! ------------------------------------
    ! Read the release/version information
    ! ------------------------------------

    READ( file_id, IOSTAT = io_status ) file_release, file_version

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading file release/version information. IOSTAT = ", i5 )' ) &
                      io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF

    ! -- Check that this file is not an OLD release/version
    IF ( ( file_release + file_version ) < ( MIN_VALID_RELEASE + MIN_VALID_VERSION ) ) THEN
      error_status = FAILURE
      WRITE( message, '( "Need to update the coefficient file. ", &
                        &"File version is ", i1, ".", i2.2, &
                        &". Oldest valid release is ",i1,".",i2.2,"." )' ) &
                      file_release,  file_version, &
                      MIN_VALID_RELEASE, MIN_VALID_VERSION
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF


    ! -- Check that this file is not a TOO NEW release/version
    ! -- i.e. update the software!
    IF ( ( file_release + file_version ) > ( MAX_VALID_RELEASE + MAX_VALID_VERSION ) ) THEN
      error_status = FAILURE
      WRITE( message, '( "Need to update the coefficient read software. ", &
                        &"File version is ", i1, ".", i2.2, &
                        &". Newest valid release is ",i1,".",i2.2,"." )' ) &
                      file_release,  file_version, &
                      MAX_VALID_RELEASE, MAX_VALID_VERSION
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF


    ! ---------------------------
    ! Read the number of channels
    ! ---------------------------

    READ( file_id, IOSTAT = io_status ) n_channels

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading total number of channels. IOSTAT = ", i5 )' ) &
                      io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF


    ! -----------------------------        
    ! Read the number of data items
    ! -----------------------------

    READ( file_id, IOSTAT = io_status ) n_items

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading total number of data items. IOSTAT = ", i5 )' ) &
                      io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF

    IF ( n_items /= N_SPECTRAL_ITEMS ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Number of data items in '//&
                            TRIM( coefficient_file )//&
                            ' inconsistent with definition.', &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF


    ! -------------------
    ! Read the data types
    ! -------------------

    READ( file_id, IOSTAT = io_status ) data_type
    
    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading data types. IOSTAT = ", i5 )' ) &
                      io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF


    DO l = 1, n_items

      IF ( data_type( l ) /= SPECTRAL_DATA_TYPE( l ) ) THEN
        error_status = FAILURE
        WRITE( message, '( "Invalid type for ", a, " data item." )' ) &
                        TRIM( SPECTRAL_DATA_NAME( l ) )
        CALL display_message( ROUTINE_NAME, &
                              TRIM( message ), &
                              error_status, &
                              message_log = message_log )
        CLOSE( file_id )
        RETURN
      END IF

    END DO


    ! --------------------------------------------
    ! Set the global definition for MAX_N_CHANNELS
    ! --------------------------------------------

    CALL get_max_n_channels( MAX_N_CHANNELS, is_set )

    IF ( is_set  ) THEN
      IF ( MAX_N_CHANNELS /= n_channels ) THEN
        error_status = WARNING
        WRITE( message, '( "MAX_N_CHANNELS set to different value, ", i4, ", ", &
                          &"than defined in coefficient file, ", i4, ". Overwriting." )' ) &
                        MAX_N_CHANNELS, n_channels
        CALL display_message( ROUTINE_NAME, &
                              TRIM( message ), &
                              error_status, &
                              message_log = message_log )
        CALL set_max_n_channels( n_channels )
      END IF
    ELSE
      CALL set_max_n_channels( n_channels )
    END IF



    !#--------------------------------------------------------------------------#
    !#          -- Allocate arrays for spectral coefficient data --             #
    !#--------------------------------------------------------------------------#

    ! -- Check if arrays are already allocated
    IF ( ALLOCATED( frequency                     ) .OR. &
         ALLOCATED( wavenumber                    ) .OR. &
         ALLOCATED( planck_c1                     ) .OR. &
         ALLOCATED( planck_c2                     ) .OR. &
         ALLOCATED( band_c1                       ) .OR. &
         ALLOCATED( band_c2                       ) .OR. &
         ALLOCATED( is_microwave_channel          ) .OR. &
         ALLOCATED( cosmic_background_temperature ) .OR. &
         ALLOCATED( cosmic_background_radiance    ) .OR. &
         ALLOCATED( solar_irradiance              ) .OR. &
         ALLOCATED( blackbody_irradiance          ) .OR. &
         ALLOCATED( is_solar_channel              )      ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Spectral coefficient data arrays already allocated.', &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    ENDIF

    ! -- If not, allocate them
    ALLOCATE( frequency( n_channels ),                     &
              wavenumber( n_channels ),                    &
              planck_c1( n_channels ),                     &
              planck_c2( n_channels ),                     &
              band_c1( n_channels ),                       &
              band_c2( n_channels ),                       &
              is_microwave_channel( n_channels ),          &
              cosmic_background_temperature( n_channels ), &
              cosmic_background_radiance( n_channels ),    &
              solar_irradiance( n_channels ),              &
              blackbody_irradiance( n_channels ),          &
              is_solar_channel( n_channels ),              &
              STAT = allocate_status                       )

    IF ( allocate_status /= 0 ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Unable to allocate arrays for spectral coefficient data.', &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF
   


    !#--------------------------------------------------------------------------#
    !#                        -- Loop over channels --                          #
    !#--------------------------------------------------------------------------#

    l_channel_loop: DO l = 1, n_channels


      ! -------------------
      ! Read channel record
      ! -------------------

      READ( file_id, IOSTAT = io_status )   &
        frequency( l ),                     &  ! Frequency
        wavenumber( l ),                    &  ! Wavenumber
        planck_c1( l ),                     &  ! First Planck constant
        planck_c2( l ),                     &  ! Second Planck constant
        band_c1( l ),                       &  ! Band correction offset
        band_c2( l ),                       &  ! Band correction slope
        is_microwave_channel( l ),          &  ! Microwave channel flag
        cosmic_background_temperature( l ), &  ! Effective cosmic background temperature
        cosmic_background_radiance( l ),    &  ! Cosmic background radiance
        solar_irradiance( l ),              &  ! Solar irradiance source
        blackbody_irradiance( l ),          &  ! Blackbody solar irradiance source
        is_solar_channel( l )                  ! Solar use flag


      IF ( io_status /= 0 ) THEN
        error_status = FAILURE
        WRITE( message, '( "Error reading channel ", i4, " coefficients. IOSTAT = ", i5 )' ) &
                        l, io_status
        CALL display_message( ROUTINE_NAME, &
                              TRIM( message ), &
                              error_status, &
                              message_log = message_log )
        RETURN
      END IF

    END DO l_channel_loop



    !#--------------------------------------------------------------------------#
    !#                          -- Close the file --                            #
    !#--------------------------------------------------------------------------#

    CLOSE( file_id, IOSTAT = io_status )

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error closing ", a, ". IOSTAT = ", i5 )' ) &
                      TRIM( coefficient_file ), io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#                      -- Successful completion --                         #
    !#--------------------------------------------------------------------------#

    ! -- Output an info message
!    WRITE( message, '( "FILE VERSION: ", i1, ".", i2.2, 2x, &
!                      &"N_CHANNELS=",i4 )' ) &
!                    file_release, file_version, &
!                    n_channels
!    CALL display_message( ROUTINE_NAME, &
!                          TRIM( message ), &
!                          INFORMATION, &
!                          message_log = message_log )
 
    error_status = SUCCESS

  END FUNCTION read_spectral_coefficients 




!------------------------------------------------------------------------------
!S+
! NAME:
!       destroy_spectral_coefficients
!
! PURPOSE:
!       PUBLIC function to deallocate the spectral coefficients PUBLIC
!       data arrays..
!
! CALLING SEQUENCE:
!       result = destroy_spectral_coefficients( message_log = message_log )
!
! INPUT ARGUMENTS:
!       None.
!
! OPTIONAL INPUT ARGUMENTS:
!       message_log:      Character string specifying a filename in which any
!                         messages will be logged. If not specified, or if an
!                         error occurs opening the log file, the default action
!                         is to output messages to the screen.
!                         UNITS:      None
!                         TYPE:       Character
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       None.
!
! OPTIONAL OUTUPT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result = SUCCESS => deallocation was successful
!              = FAILURE => error occurred deallocating coefficient data arrays
!
! CALLS:
!      display_message:         Subroutine to output messages
!                               SOURCE: error_handler module
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
!S-
!------------------------------------------------------------------------------

  FUNCTION destroy_spectral_coefficients( message_log ) &
                                        RESULT ( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    CHARACTER( * ), INTENT( IN ), OPTIONAL :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'DESTROY_SPECTRAL_COEFFICIENTS'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: allocate_status



    !#--------------------------------------------------------------------------#
    !#           -- Deallocate arrays for spectral coefficient data --          #
    !#--------------------------------------------------------------------------#

    DEALLOCATE( frequency,                     &
                wavenumber,                    &
                planck_c1,                     &
                planck_c2,                     &
                band_c1,                       &
                band_c2,                       &
                is_microwave_channel,          &
                cosmic_background_temperature, &
                cosmic_background_radiance,    &
                solar_irradiance,              &
                blackbody_irradiance,          &
                is_solar_channel,              &
                STAT = allocate_status         )

    IF ( allocate_status /= 0 ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Error occurred deallocating spectral coefficient data arrays.', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#           -- Reset the MAX_N_CHANNELS "pseudo-parameter" --              #
    !#--------------------------------------------------------------------------#

    CALL reset_max_n_channels()


    !#--------------------------------------------------------------------------#
    !#                      -- Successful completion --                         #
    !#--------------------------------------------------------------------------#
 
    error_status = SUCCESS

  END FUNCTION destroy_spectral_coefficients 





!------------------------------------------------------------------------------
!
! NAME:
!       write_spectral_coefficients
!
! PURPOSE:
!       Function to write the spectral coefficient data to file.
!
! CATEGORY:
!       NCEP RTM
!
! LANGUAGE:
!       Fortran-90
!
! CALLING SEQUENCE:
!       result = write_spectral_coefficients( coefficient_file,              &  ! Input
!                                             frequency,                     &  ! Input, L
!                                             wavenumber,                    &  ! Input, L
!                                             planck_c1,                     &  ! Input, L
!                                             planck_c2,                     &  ! Input, L
!                                             band_c1,                       &  ! Input, L
!                                             band_c2,                       &  ! Input, L
!                                             is_microwave_channel,          &  ! Input, L
!                                             cosmic_background_temperature, &  ! Input, L
!                                             cosmic_background_radiance     &  ! Input, L
!                                             solar_irradiance,              &  ! Input, L
!                                             blackbody_irradiance,          &  ! Input, L
!                                             is_solar_channel,              &  ! Input, L
!                                             message_log                    )  ! Optional input
!
! INPUT ARGUMENTS:
!       coefficient_file:              Name of the file to which the spectral
!                                      data is to be written.
!                                      UNITS:      None
!                                      TYPE:       Character
!                                      DIMENSION:  Scalar
!                                      ATTRIBUTES: INTENT( IN )
!
!       frequency:                     Channel frequency
!                                      UNITS:      GHz
!                                      TYPE:       Double
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       wavenumber:                    Channel frequency
!                                      UNITS:      cm^-1
!                                      TYPE:       Double
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       planck_c1:                     First Planck function coefficient
!                                      UNITS:      mW/(m^2.sr.cm^-2)
!                                      TYPE:       Double
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       planck_c2:                     Second Planck function coefficient
!                                      UNITS:      K/cm^1
!                                      TYPE:       Double
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       band_c1:                       Polychromatic band correction offset.
!                                      UNITS:      K
!                                      TYPE:       Double
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       band_c2:                       Polychromatic band correction slope.
!                                      UNITS:      K/K
!                                      TYPE:       Double
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       is_microwave_channel:          Flag indicating if the particular satellite
!                                      channel is a microwave channel.
!                                      If = 0, IR channel
!                                         = 1, uW channel
!                                      UNITS:      None
!                                      TYPE:       Long integer
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       cosmic_background_temperature: Cosmic background temperature in Kelvin to
!                                      use in the radiative transfer. These values
!                                      may be frequency dependent for microwave channels
!                                      to account for instrument calibration using the
!                                      Rayleigh-Jeans approximation, but the radiative
!                                      transfer using the full Planck function
!                                      expression. IR channels use a fixed value of
!                                      2.736K.
!                                      UNITS:      K
!                                      TYPE:       Double
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       cosmic_background_radiance:    Cosmic background radiance used to initialise
!                                      the radiative transfer. These values
!                                      may be frequency dependent for microwave channels
!                                      to account for instrument calibration using the
!                                      Rayleigh-Jeans approximation, but the radiative
!                                      transfer using the full Planck function
!                                      expression. IR channel values are all set to 0.0.
!                                      UNITS:      mW/(m^2.sr.cm^-1)
!                                      TYPE:       Double
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       solar_irradiance:              Exterrestrial solar irradiance source values
!                                      (based on Kurucz).
!                                      UNITS:      mW/(m^2.cm^-1)
!                                      TYPE:       Double
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       blackbody_irradiance:          Exterrestrial solar irradiance source values
!                                      for a blackbody source at the equivalent solar
!                                      temperature of 5783K.
!                                      UNITS:      mW/(m^2.cm^-1)
!                                      TYPE:       Double
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
!       is_solar_channel:              Array of integer flag values indicating if the
!                                      particular satellite channel possesses valid
!                                      solar irradiance values.
!                                      If = 0, No
!                                         = 1, Yes
!                                      UNITS:      None
!                                      TYPE:       Long
!                                      DIMENSION:  L, n_channels
!                                      ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       message_log:                   Character string specifying a filename in which
!                                      any messages will be logged. If not specified,
!                                      or if an error occurs opening the log file,
!                                      the default action is to output messages to
!                                      the screen.
!                                      UNITS:      None
!                                      TYPE:       Character
!                                      DIMENSION:  Scalar
!                                      ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       None.
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result = SUCCESS => data file write was successful
!              = FAILURE => error occurred opening or writing to the
!                           data file
!
! CALLS:
!       display_message:        Subroutine to output messages
!                               SOURCE: error_handler module
!
! CONTAINS:
!       None.
!
! INCLUDE FILES:
!       None.
!
! EXTERNALS:
!       None.
!
! COMMON BLOCKS:
!       None.
!
! FILES ACCESSED:
!       None.
!
! SIDE EFFECTS:
!       If a coefficient file with the specified name exists, it is
!       overwritten.
!
! RESTRICTIONS:
!       None.
!
!------------------------------------------------------------------------------

  FUNCTION write_spectral_coefficients( coefficient_file,              &  ! Input
                                        frequency,                     &  ! Input
                                        wavenumber,                    &  ! Input
                                        planck_c1,                     &  ! Input
                                        planck_c2,                     &  ! Input
                                        band_c1,                       &  ! Input
                                        band_c2,                       &  ! Input
                                        is_microwave_channel,          &  ! Input
                                        cosmic_background_temperature, &  ! Input
                                        cosmic_background_radiance,    &  ! Input
                                        solar_irradiance,              &  ! Input
                                        blackbody_irradiance,          &  ! Input
                                        is_solar_channel,              &  ! Input

                                        release,                       &  ! Optional input
                                        version,                       &  ! Optional input
                                        message_log )                  &  ! Optional input
                                      RESULT ( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Input
    CHARACTER( * ),                  INTENT( IN ) :: coefficient_file              ! Input
    REAL( Double ),  DIMENSION( : ), INTENT( IN ) :: frequency                     ! Input, L
    REAL( Double ),  DIMENSION( : ), INTENT( IN ) :: wavenumber                    ! Input, L
    REAL( Double ),  DIMENSION( : ), INTENT( IN ) :: planck_c1                     ! Input, L
    REAL( Double ),  DIMENSION( : ), INTENT( IN ) :: planck_c2                     ! Input, L
    REAL( Double ),  DIMENSION( : ), INTENT( IN ) :: band_c1                       ! Input, L
    REAL( Double ),  DIMENSION( : ), INTENT( IN ) :: band_c2                       ! Input, L
    INTEGER( Long ), DIMENSION( : ), INTENT( IN ) :: is_microwave_channel          ! Input, L
    REAL( Double ),  DIMENSION( : ), INTENT( IN ) :: cosmic_background_temperature ! Input, L
    REAL( Double ),  DIMENSION( : ), INTENT( IN ) :: cosmic_background_radiance    ! Input, L
    REAL( Double ),  DIMENSION( : ), INTENT( IN ) :: solar_irradiance              ! Input, L
    REAL( Double ),  DIMENSION( : ), INTENT( IN ) :: blackbody_irradiance          ! Input, L
    INTEGER( Long ), DIMENSION( : ), INTENT( IN ) :: is_solar_channel              ! Input, L


    INTEGER,        OPTIONAL,        INTENT( IN ) :: release
    INTEGER,        OPTIONAL,        INTENT( IN ) :: version

    CHARACTER( * ), OPTIONAL,        INTENT( IN ) :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ),  PARAMETER :: ROUTINE_NAME = 'WRITE_SPECTRAL_COEFFICIENTS'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- Typed ( Long ) for output
    INTEGER( Long ) :: n_channels

    INTEGER( Long ) :: file_release
    INTEGER( Long ) :: file_version

    INTEGER :: file_id
    INTEGER :: io_status

    INTEGER :: i, j, k, l


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC SIZE, &
              TRIM



    !#--------------------------------------------------------------------------#
    !#                           -- Check input --                              #
    !#--------------------------------------------------------------------------#

    n_channels = SIZE( frequency )

    IF ( SIZE( wavenumber                    ) /= n_channels .OR. &
         SIZE( planck_c1                     ) /= n_channels .OR. &
         SIZE( planck_c2                     ) /= n_channels .OR. &
         SIZE( band_c1                       ) /= n_channels .OR. &
         SIZE( band_c1                       ) /= n_channels .OR. &
         SIZE( is_microwave_channel          ) /= n_channels .OR. &
         SIZE( cosmic_background_temperature ) /= n_channels .OR. &
         SIZE( cosmic_background_radiance    ) /= n_channels .OR. &
         SIZE( solar_irradiance              ) /= n_channels .OR. &
         SIZE( blackbody_irradiance          ) /= n_channels .OR. &
         SIZE( is_solar_channel              ) /= n_channels      ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Inconsistent vector channel dimensions.', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF


    ! ----------------------------------
    ! Assign release and version numbers
    ! ----------------------------------

    IF ( PRESENT( release ) ) THEN
      file_release = INT( release, Long )
    ELSE
      file_release = MAX_VALID_RELEASE
    END IF

    IF ( PRESENT( version ) ) THEN
      file_version = INT( version, Long )
    ELSE
      file_version = MAX_VALID_VERSION
    END IF



    !#--------------------------------------------------------------------------#
    !#        -- Open the spectral coefficient data file for output --          #
    !#--------------------------------------------------------------------------#

    error_status = open_coefficient_file( coefficient_file,         &
                                          file_id,                  &
                                          for_output  = 1,          &
                                          message_log = message_log )

    IF ( error_status /= SUCCESS ) THEN
      CALL display_message( ROUTINE_NAME, &
                            'Error opening '//TRIM( coefficient_file )//' for output.', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#                     -- Write the "file header" --                        #
    !#--------------------------------------------------------------------------#

    ! -------------------------------------
    ! Write the release/version information
    ! -------------------------------------

    WRITE( file_id, IOSTAT = io_status ) file_release, file_version 

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error writing file release/version information to ", a, &
                        &". IOSTAT = ", i5 )' ) &
                      TRIM( coefficient_file ), io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF


    ! --------------------
    ! Write the dimensions
    ! --------------------

    WRITE( file_id, IOSTAT = io_status ) n_channels

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error writing dimension values to ", a, ". IOSTAT = ", i5 )' ) &
                      TRIM( coefficient_file ), io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF


    ! ----------------------------------------------
    ! Write the number of data items per channel and
    ! the type of the data items.
    ! ----------------------------------------------

    WRITE( file_id, IOSTAT = io_status ) N_SPECTRAL_ITEMS

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error writing data item count to ", a, ". IOSTAT = ", i5 )' ) &
                      TRIM( coefficient_file ), io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF


    WRITE( file_id, IOSTAT = io_status ) SPECTRAL_DATA_TYPE

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error writing data item types to ", a, ". IOSTAT = ", i5 )' ) &
                      TRIM( coefficient_file ), io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF


    !#--------------------------------------------------------------------------#
    !#                     -- Write the data by channel --                      #
    !#--------------------------------------------------------------------------#

    DO l = 1, n_channels

      WRITE( file_id, IOSTAT = io_status ) frequency( l ),                     &                    
                                           wavenumber( l ),                    &
                                           planck_c1( l ),                     &
                                           planck_c2( l ),                     &
                                           band_c1( l ),                       &
                                           band_c2( l ),                       &
                                           is_microwave_channel( l ),          &
                                           cosmic_background_temperature( l ), &
                                           cosmic_background_radiance( l ),    &
                                           solar_irradiance( l ),              &
                                           blackbody_irradiance( l ),          &
                                           is_solar_channel( l )

      IF ( io_status /= 0 ) THEN
        error_status = FAILURE
        WRITE( message, '( "Error writing channel ", i4, 1x, &
                          &"spectral data. IOSTAT = ", i5 )' ) &
                        l, io_status
        CALL display_message( ROUTINE_NAME, &
                              TRIM( message ), &
                              error_status, &
                              message_log = message_log )
        CLOSE( file_id )
        RETURN
      END IF

    END DO



    !#--------------------------------------------------------------------------#
    !#                          -- Close the file --                            #
    !#--------------------------------------------------------------------------#

    CLOSE( file_id, IOSTAT = io_status )

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error closing ", a, ". IOSTAT = ", i5 )' ) &
                      TRIM( coefficient_file ), io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#                      -- Successful completion --                         #
    !#--------------------------------------------------------------------------#

    ! -- Output an info message
    ! -- Output an info message
    WRITE( message, '( "FILE VERSION: ", i1, ".", i2.2, 2x, &
                      &"N_CHANNELS=",i4 )' ) &
                    file_release, file_version, &
                    n_channels
    CALL display_message( ROUTINE_NAME, &
                          TRIM( message ), &
                          INFORMATION, &
                          message_log = message_log )

    error_status = SUCCESS

  END FUNCTION write_spectral_coefficients

END MODULE spectral_coefficients


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
! Revision 1.1  2002/11/15 15:21:33  birk
! Added to cvs mainly to see how this compiles on other platforms, it currently
! seems to compile on the IBM
!
! Revision 1.11  2001/08/31 21:11:41  paulv
! - Added MIN and MAX release/version parameters to allow for valid use of
!   data files within a specified range.
!
! Revision 1.10  2001/08/16 17:16:30  paulv
! - Updated documentation
! - The comparison of n_channels and MAX_N_CHANNELS is now done via the
!   MAX_N_CHANNELS methods in the PARAMETERS module.
!
! Revision 1.9  2001/08/09 20:45:33  paulv
! - Added the WRITE_SPECTRAL_COEFFICIENTS function.
! - Moved all the spectral data type and name definitions from the
!   COEFFICIENT_UTILITY module to this one. Altered USE statement of the
!   COEFFICIENT_UTILITY module to reflect this change.
! - Added VALID_RELEASE and VALID_VERSION parameters for data file version
!   checking.
! - Added data file release and version number read/write to the requisite
!   read/write function.
!
! Revision 1.8  2001/07/12 17:46:31  paulv
! - Removed definitions of the number, type, and name of the spectral items
!   and moved them into the COEFFICIENT_UTILITY module. They are now available
!   via:
!     USE coefficient_utility, ONLY: open_coefficient_file, &
!                                    N_SPECTRAL_ITEMS,      &
!                                    SPECTRAL_DATA_TYPE,    &
!                                    SPECTRAL_DATA_NAME
!   This was done to allow the COEFFICIENT_UTILITY module to be used for
!   reformatting. However, this may change - now definitions for the contents
!   of the spectral coefficient data file are distributed in two different
!   modules. I don't like that.
! - Compressed the coefficient READ statement.
!
! Revision 1.7  2001/05/29 17:49:32  paulv
! - Made ALL real valued spectral coefficients double precision. Previously
!   the cosmic background and solar terms were single precision.
! - Added precalculated cosmic background radiance term.
! - Some cosmetic changes.
! - Changed USE_SOLAR array name to IS_SOLAR_CHANNEL.
! - Added DESTROY_SPECTRAL_COEFFICIENTS function. This deallocates the
!   data arrays used to store the spectral coefficient data.
!
! Revision 1.6  2001/01/09 21:23:27  paulv
! - Added IS_MICROWAVE_CHANNEL and COSMIC_BACKGROUND_TEMPERATURE arrays to
!   read.
! - Updated module documentation.
!
! Revision 1.5  2000/11/09 20:39:49  paulv
! - Coefficient arrays are now ALLOCATABLE.
! - Input file format has changed to contain data dimension and type
!   information for file data checks and array allocation.
!
! Revision 1.4  2000/08/31 19:36:33  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 1.3  2000/08/24 16:06:28  paulv
! - Removed references to the record length parameter. No longer needed as
!   file access is SEQUENTIAL rather than DIRECT.
! - Replaced error check after OPEN_COEFFICIENT_FILE call with a simple
!   test for error_status /= SUCCESS. The open function no longer returns
!   any error status other than SUCCESS or FAILURE (used to return WARNING
!   in some circumstances.)
! - "REC =" keyword removed from file READ statement.
! - Channel loop construct name changed from "channel_loop" to "l_channel_loop"
!   to indicate the loop counter variable is "l". This is not a big deal for
!   this situation but has proven useful in other modules with a high degree
!   of nested loops.
! - Updated module and subprogram documentation.
!
! Revision 1.2  2000/08/08 17:04:02  paulv
! Module modified to:
! - Read the spectral coefficients correctly! and
! - To use the PARAMETERS module rather than the CONSTANTS module.
!
! Revision 1.1  2000/07/12 16:08:10  paulv
! Initial checked in version
!
!
!
