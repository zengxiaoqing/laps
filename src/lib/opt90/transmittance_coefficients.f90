!------------------------------------------------------------------------------
!M+
! NAME:
!       transmittance_coefficients
!
! PURPOSE:
!       Module to hold RT model transmittance coefficients 
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE transmittance_coefficients
!
! OUTPUTS:
!       absorber_space_levels:     Array containing teh absorber space levels
!                                  used in generating the transmittance
!                                  coefficients.
!                                  UNITS:      Depends on absorber type
!                                  TYPE:       Double
!                                  DIMENSION:  0:Ka x J
!                                  ATTRIBUTES: SAVE, PUBLIC
!
!       predictor_index:           Array containing the predictor indices used to
!                                  identify predictors in the transmittance model.
!                                  UNITS:      None.
!                                  TYPE:       Integer
!                                  DIMENSION:  0:I x L x J
!                                  ATTRIBUTES: SAVE, PUBLIC
!
!       tau_coefficients:          Array containing the transmittance model
!                                  coefficients.
!                                  UNITS:      Absorber, predictor dependent
!                                  TYPE:       Double
!                                  DIMENSION:  0:I x 0:Ka x L x J
!                                  ATTRIBUTES: SAVE, PUBLIC
!
! MODULES:
!       type_kinds:                Module containing data type kind definitions.
!
!       error_handler:             Module to define error codes and handle error
!                                  conditions
!
!       parameters:                Module containing parameter definitions for the
!                                  RT model.
!
!       coefficient_utility:       Module containing coefficient file utility
!                                  subprograms.
!
! CONTAINS:
!       read_tau_coefficients:     PUBLIC function to read the transmittance model
!                                  coefficients and predictor indices and fill the
!                                  PUBLIC data arrays.
!
!       destroy_tau_coefficients:  PUBLIC function to release the memory that
!                                  was allocated for the transmittance coefficient
!                                  data.
!
!       write_tau_coefficients:    PUBLIC function to write the new format
!                                  transmittance coefficient data file.
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

MODULE transmittance_coefficients


  ! ---------------------
  ! Module use statements
  ! ---------------------

  USE type_kinds
  USE error_handler
  USE parameters
  USE coefficient_utility


  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------

  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------

  PRIVATE
  PUBLIC :: read_tau_coefficients
  PUBLIC :: destroy_tau_coefficients
  PUBLIC :: write_tau_coefficients



  ! -----------------
  ! Module parameters
  ! -----------------

  ! -- Transmittance coefficient file version numbers
  INTEGER( Long ), PARAMETER, PRIVATE :: MIN_VALID_RELEASE = 1_Long
  INTEGER( Long ), PARAMETER, PRIVATE :: MAX_VALID_RELEASE = 1_Long

  INTEGER( Long ), PARAMETER, PRIVATE :: MIN_VALID_VERSION = 4_Long
  INTEGER( Long ), PARAMETER, PRIVATE :: MAX_VALID_VERSION = 4_Long

  ! -- Transmittance coefficient file descriptors:
  ! -- Number of items in the transmittance coefficient file
  INTEGER( Long ), PARAMETER, PRIVATE :: N_TRANSMITTANCE_ITEMS = 4_Long

  ! -- Data types of the transmittance coefficient data
  !    5 = Double (i.e. 8-byte float)
  !    4 = Single (i.e. 4-byte float)
  !    3 = Long   (i.e. 4-byte integer)
  INTEGER( Long ), PARAMETER, PRIVATE, &
                   DIMENSION( N_TRANSMITTANCE_ITEMS ) :: TRANSMITTANCE_DATA_TYPE = &
                                                         (/ 5_Long, 5_Long, 3_Long, 5_Long /)

  ! -- Names of the data items (for error processing)
  CHARACTER( * ), PARAMETER, PRIVATE, &
                  DIMENSION( N_TRANSMITTANCE_ITEMS ) :: TRANSMITTANCE_DATA_NAME = &
                                                        (/ 'ALPHA                    ', &
                                                           'ABSORBER_SPACE_LEVELS    ', &
                                                           'PREDICTOR_INDEX          ', &
                                                           'TRANSMITTANCE_COEFFICIENT' /)

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

  ! -- Alpha value
  !
  ! Dimension is: number of absorbers
  REAL( Double ), SAVE, &
                  PUBLIC, &
                  ALLOCATABLE, &
                  DIMENSION( : ) :: alpha


  ! -- Absorber space levels
  !
  ! Eventual dimensions are:
  !   #1: 0->maximum number of absorber layers  (fixed by transmittance algorithm)
  !   #2: 1->number of absorbers                (fixed by transmittance algorithm)
  !
  ! Although both dimensions are fixed by the transmittance
  ! algorithm, there is the possibility of two different sets
  ! of transmittance coefficients existing: when different
  ! profile sets are used to generate the coefficients. By
  ! including the data with the transmittance coefficients
  ! rather than calculating them every time the model is
  ! initialised using different coeficient sets is a simple
  ! matter of reading a different file.
  REAL( Double ), SAVE, &
                  PUBLIC, &
                  ALLOCATABLE, &
                  DIMENSION( :, : ) :: absorber_space_levels



  ! -- Absorber predictor index arrays
  !
  ! Remember the 0'th index is an indicator of absorption for
  ! a particular channel. If, for any channel, l, 
  ! predictor_index( 0,l ) = 0 then there is no absorption due
  ! to the absorber in channel l.
  !
  ! Eventual dimensions are:
  !   #1: 0->maximum number of predictors used - (fixed by transmittance algorithm)
  !   #2: 1->number of channels requested      - (user selectable)
  !   #3: 1->number of absorbers               - (fixed by transmittance algorithm)
  INTEGER( Long ), SAVE, &
                   PUBLIC, &
                   ALLOCATABLE, &
                   DIMENSION( :, :, : ) :: predictor_index



  ! -- Coefficients for calculating transmittances
  !
  ! The 0'th predictor coefficient is the offset term used in:
  !
  !                                  n
  ! absorption_coefficient = b(0) + sum{ b(i)*X(i) }
  !                                 i=1
  !
  ! where n = n_predictors_used
  !
  ! Eventual dimensions are:
  !   #1: 0->maximum number of predictors used - (fixed by transmittance algorithm)
  !   #2: 0->number of absorber layers         - (fixed by transmittance algorithm)
  !   #3: 1->number of channels requested      - (user selectable)
  !   #4: 1->number of absorbers               - (fixed by transmittance algorithm)
  REAL( Double ), SAVE, &
                  PUBLIC, &
                  ALLOCATABLE, &
                  DIMENSION( :, :, :, : ) :: tau_coefficients



CONTAINS


!------------------------------------------------------------------------------
!S+
! NAME:
!       read_tau_coefficients
!
! PURPOSE:
!       PUBLIC function to read the transmittance model coefficients and
!       predictor indices and fill the PUBLIC data arrays.
!
! CALLING SEQUENCE:
!       result = read_tau_coefficients( coefficient_file, &
!                                       message_log = message_log )
!
! INPUT ARGUMENTS:
!       coefficient_file: Name of the file containing the coefficient data.
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
!S-
!------------------------------------------------------------------------------

  FUNCTION read_tau_coefficients( coefficient_file, &
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

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'READ_TAU_COEFFICIENTS'


    ! ---------------
    ! Local variables
    ! ---------------

    ! -- These are dimensioned ( Long ) as they are read in.
    INTEGER( Long ) :: file_release
    INTEGER( Long ) :: file_version

    INTEGER( Long ) :: n_predictors_to_use
    INTEGER( Long ) :: n_absorber_layers
    INTEGER( Long ) :: n_channels
    INTEGER( Long ) :: n_absorbers
    INTEGER( Long ) :: n_items
    INTEGER, DIMENSION( N_TRANSMITTANCE_ITEMS ) :: data_type

    INTEGER :: i, j, k, l
    INTEGER :: io_status
    INTEGER :: file_id
    INTEGER :: allocate_status

    LOGICAL :: min_relver
    LOGICAL :: max_relver

    ! -- Maximum channels pseudo parameter
    INTEGER :: MAX_N_CHANNELS
    LOGICAL :: is_set


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC TRIM



    !#--------------------------------------------------------------------------#
    !#           -- Open the transmittance coefficient data file --             #
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


    ! -------------------
    ! Read the dimensions
    ! -------------------

    READ( file_id, IOSTAT = io_status ) n_predictors_to_use, &
                                        n_absorber_layers, &
                                        n_channels, &
                                        n_absorbers

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading data dimensions. IOSTAT = ", i5 )' ) &
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

    IF ( n_items /= N_TRANSMITTANCE_ITEMS ) THEN
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

      IF ( data_type( l ) /= TRANSMITTANCE_DATA_TYPE( l ) ) THEN
        error_status = FAILURE
        WRITE( message, '( "Invalid type for ", a, " data item." )' ) &
                        TRIM( TRANSMITTANCE_DATA_NAME( l ) )
        CALL display_message( ROUTINE_NAME, &
                              TRIM( message ), &
                              error_status, &
                              message_log = message_log )
        CLOSE( file_id )
        RETURN
      END IF

    END DO


    !#--------------------------------------------------------------------------#
    !#                     -- Check the dimension values --                     #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------------------------------
    ! Check the dimensions that are transmittance algorithm dependent:
    !   n_predictors_to_use, n_absorber_layers, n_absorbers
    ! ----------------------------------------------------------------

    ! -- Number of predictors to use
    IF ( n_predictors_to_use /= MAX_N_PREDICTORS_USED ) THEN
      error_status = FAILURE
      WRITE( message, '( "Predictor index dimension inconsistent. Value is ", i3, &
                        &" but should be ", i1 )' ) &
                      n_predictors_to_use, MAX_N_PREDICTORS_USED
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF

    ! -- Number of absorber layers
    IF ( n_absorber_layers /= MAX_N_ABSORBER_LAYERS ) THEN
      error_status = FAILURE
      WRITE( message, '( "Absorber layer dimension inconsistent. Value is ", i3, &
                        &" but should be ", i3 )' ) &
                      n_absorber_layers, MAX_N_ABSORBER_LAYERS
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF

    ! -- Number of absorbers
    IF ( n_absorbers /= MAX_N_ABSORBERS ) THEN
      error_status = FAILURE
      WRITE( message, '( "Absorber dimension inconsistent. Value is ", i3, &
                        &" but should be ", i1 )' ) &
                      n_absorbers, MAX_N_ABSORBERS
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF


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
    !#          -- Allocate arrays for transmittance coefficient data --        #
    !#--------------------------------------------------------------------------#

    ! -- Check if arrays are already allocated
    IF ( ALLOCATED( alpha                 ) .OR. &
         ALLOCATED( absorber_space_levels ) .OR. &
         ALLOCATED( predictor_index       ) .OR. &
         ALLOCATED( tau_coefficients      )      ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Transmittance coefficient data arrays already allocated.', &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF

    ! -- If not, allocate them
    ALLOCATE( alpha( n_absorbers ), &

              absorber_space_levels( 0:n_absorber_layers, &
                                     n_absorbers          ), &

              predictor_index( 0:n_predictors_to_use, &
                               n_channels,            &
                               n_absorbers            ), &

              tau_coefficients( 0:n_predictors_to_use, &
                                0:n_absorber_layers,   &
                                n_channels,            &
                                n_absorbers            ), &
              STAT = allocate_status )

    IF ( allocate_status /= 0 ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Unable to allocate arrays for transmittance coefficient data.', &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF
   


    !#--------------------------------------------------------------------------#
    !#                 -- Read the absorber space data --                       #
    !#--------------------------------------------------------------------------#

    ! -- Alpha values
    READ( file_id, IOSTAT = io_status ) alpha

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading alpha values. IOSTAT = ", i5 )' ) &
                      io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF

    ! -- Absorber space levels
    READ( file_id, IOSTAT = io_status ) absorber_space_levels

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading absorber space levels. IOSTAT = ", i5 )' ) &
                      io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#                        -- Loop over channels --                          #
    !#--------------------------------------------------------------------------#

    l_channel_loop: DO l = 1, n_channels


      ! ----------------------
      ! Read predictor indices
      ! ----------------------

      READ( file_id, IOSTAT = io_status ) predictor_index( :, l, : )

      IF ( io_status /= 0 ) THEN
        error_status = FAILURE
        WRITE( message, '( "Error reading channel ", i4, " predictor indices. IOSTAT = ", i5 )' ) &
                        l, io_status
        CALL display_message( ROUTINE_NAME, &
                              TRIM( message ), &
                              error_status, &
                              message_log = message_log )
        RETURN
      END IF


      ! -------------------------------
      ! Read transmittance coefficients
      ! -------------------------------

      READ( file_id, IOSTAT = io_status ) tau_coefficients( :, :, l, : )

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
                      coefficient_file, io_status
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
!                      &"N_PREDICTORS_TO_USE=",i1,2x,&
!                      &"N_ABSORBER_LAYERS=",i3,2x,&
!                      &"N_CHANNELS=",i4,2x,&
!                      &"N_ABSORBERS=",i1 )' ) &
!                    file_release, file_version, &
!                    n_predictors_to_use, &
!                    n_absorber_layers, &
!                    n_channels, &
!                    n_absorbers
!    CALL display_message( ROUTINE_NAME, &
!                          TRIM( message ), &
!                          INFORMATION, &
!                          message_log = message_log )

    error_status = SUCCESS

  END FUNCTION read_tau_coefficients 



!------------------------------------------------------------------------------
!S+
! NAME:
!       destroy_tau_coefficients
!
! PURPOSE:
!       PUBLIC function to deallocate the alpha, absorber space level, 
!       predictor index and transmittance model coefficient PUBLIC data arrays.
!
! CALLING SEQUENCE:
!       result = destroy_tau_coefficients( message_log = message_log )
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
!              = FAILURE => error occurred deallocating coefficient arrays.
!
! CALLS:
!      display_message:         Subroutine to output messages
!                               SOURCE: error_handler module
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

  FUNCTION destroy_tau_coefficients( message_log ) &
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

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'DESTROY_TAU_COEFFICIENTS'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: allocate_status



    !#--------------------------------------------------------------------------#
    !#         -- Deallocate arrays for transmittance coefficient data --       #
    !#--------------------------------------------------------------------------#

    DEALLOCATE( alpha,                 &
                absorber_space_levels, &
                predictor_index,       &
                tau_coefficients,      &
                STAT = allocate_status )

    IF ( allocate_status /= 0 ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Error occurred deallocating transmittance coefficient data arrays.', &
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

  END FUNCTION destroy_tau_coefficients 





!------------------------------------------------------------------------------
!
! NAME:
!       write_tau_coefficients
!
! PURPOSE:
!       Function to write the transmittance coefficient data to file.
!
! CATEGORY:
!       NCEP RTM
!
! LANGUAGE:
!       Fortran-90
!
! CALLING SEQUENCE:
!       result = write_tau_coefficients( coefficient_file,           &  ! Input
!                                        alpha,                      &  ! Input
!                                        absorber_space_levels,      &  ! Input
!                                        predictor_index,            &  ! Input, 0:Iuse        x L x J
!                                        transmittance_coefficients, &  ! Input, 0:Iuse x 0:Ka x L x J
!                                        message_log                 )  ! Optional input
!
! INPUT ARGUMENTS:
!       coefficient_file:            Name of the file to which the transmittance
!                                    coefficient data is to be written.
!                                    UNITS:      None
!                                    TYPE:       Character
!                                    DIMENSION:  Scalar
!                                    ATTRIBUTES: INTENT( IN )
!
!       alpha:                       Vector containing the exponential parameter
!                                    for each absorber used to construct the
!                                    absorber space.
!                                    ** NOTE: This variable is in a different scope
!                                             from the module data variable of the
!                                             same name.
!                                    UNITS:      None
!                                    TYPE:       Double
!                                    DIMENSION:  J
!                                    ATTRIBUTES: INTENT( IN )
!
!       absorber_space_levels:       Array containing the absorber space levels
!                                    used in the transmittance model.
!                                    ** NOTE: This variable is in a different scope
!                                             from the module data variable of the
!                                             same name.
!                                    UNITS:      Varies with absorber
!                                    TYPE:       Double
!                                    DIMENSION:  0:Ka x J
!                                    ATTRIBUTES: INTENT( IN )
!
!       predictor_index:             Array containing the predictor indices
!                                    used to identify which predictors to use
!                                    in the transmittance model.
!                                    ** NOTE: This variable is in a different scope
!                                             from the module data variable of the
!                                             same name.
!
!                                    Remember the 0'th index is an indicator
!                                    of absorption for a particular channel.
!                                    If, for any channel, l,
!                                      predictor_index( 0,l ) = 0
!                                    then there is no absorption due to that
!                                    absorber in channel l. If the 0'th index
!                                    is non-zero, then it is the number of
!                                    valid predictor indices to follow (may
!                                    not always need all the predictors).
!                                    UNITS:      None
!                                    TYPE:       Long
!                                    DIMENSION:  0:Iuse x L x J
!                                    ATTRIBUTES: INTENT( IN )
!
!       transmittance_coefficients:  Array containing the transmittance
!                                    model coefficients in absorber, and
!                                    predictor dependent units.
!                                    ** NOTE: This variable is in a different scope
!                                             from the module data variable of the
!                                             same name.
!
!                                    The 0'th predictor coefficient is the
!                                    offset term used in:
!                                                         __ n
!                                                        \
!                                      abs_coeff = b(0) + >  b(i)*X(i)
!                                                        /__
!                                                           i=1
!
!                                    where n = n_predictors_used (or the 0'th
!                                    term of predictor_index if it is non-zero).
!                                    UNITS:      None
!                                    TYPE:       Double
!                                    DIMENSION:  0:Iuse x 0:Ka x L x J
!                                    ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       message_log:                 Character string specifying a filename in which
!                                    any messages will be logged. If not specified,
!                                    or if an error occurs opening the log file,
!                                    the default action is to output messages to
!                                    the screen.
!                                    UNITS:      None
!                                    TYPE:       Character
!                                    DIMENSION:  Scalar
!                                    ATTRIBUTES: INTENT( IN ), OPTIONAL
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

  FUNCTION write_tau_coefficients( coefficient_file,           &  ! Input
                                   alpha,                      &  ! Input
                                   absorber_space_levels,      &  ! Input
                                   predictor_index,            &  ! Input
                                   transmittance_coefficients, &  ! Input

                                   release,                    &  ! Optional input
                                   version,                    &  ! Optional input
                                   message_log )               &  ! Optional input
                                 RESULT ( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    ! -- Input
    CHARACTER( * ),                             INTENT( IN ) :: coefficient_file            ! Input
    REAL( Double ),  DIMENSION(            : ), INTENT( IN ) :: alpha                       ! Input, J
    REAL( Double ),  DIMENSION(     0:,    : ), INTENT( IN ) :: absorber_space_levels       ! Input, 0:Ka x J
    INTEGER( Long ), DIMENSION( 0:,     :, : ), INTENT( IN ) :: predictor_index             ! Input, 0:Iuse x L x J
    REAL( Double ),  DIMENSION( 0:, 0:, :, : ), INTENT( IN ) :: transmittance_coefficients  ! Input, 0:Iuse x 0:Ka x L x J
    !                            ^   ^  ^  ^
    !                            I   K  L  J

    INTEGER,        OPTIONAL,                   INTENT( IN ) :: release
    INTEGER,        OPTIONAL,                   INTENT( IN ) :: version

    CHARACTER( * ), OPTIONAL,                   INTENT( IN ) :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ),  PARAMETER :: ROUTINE_NAME = 'WRITE_TAU_COEFFICIENTS'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: n_alpha_absorbers 

    INTEGER :: n_as_layers
    INTEGER :: n_as_absorbers 

    INTEGER :: n_pi_predictors
    INTEGER :: n_pi_channels  
    INTEGER :: n_pi_absorbers 

    ! -- These are typed ( Long ) as they will be output
    INTEGER( Long ) :: n_tc_predictors
    INTEGER( Long ) :: n_tc_layers    
    INTEGER( Long ) :: n_tc_channels  
    INTEGER( Long ) :: n_tc_absorbers 

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

    ! -------------------------------------------
    ! Check the predictor index and transmittance
    ! coefficient input arrays
    ! -------------------------------------------

    n_alpha_absorbers = SIZE( alpha )

    n_as_layers       = SIZE( absorber_space_levels, DIM = 1 ) - 1
    n_as_absorbers    = SIZE( absorber_space_levels, DIM = 2 )

    n_pi_predictors   = SIZE( predictor_index, DIM = 1 ) - 1
    n_pi_channels     = SIZE( predictor_index, DIM = 2 )
    n_pi_absorbers    = SIZE( predictor_index, DIM = 3 )

    n_tc_predictors   = SIZE( transmittance_coefficients, DIM = 1 ) - 1
    n_tc_layers       = SIZE( transmittance_coefficients, DIM = 2 ) - 1
    n_tc_channels     = SIZE( transmittance_coefficients, DIM = 3 )
    n_tc_absorbers    = SIZE( transmittance_coefficients, DIM = 4 )

    ! -- Check number of predictors
    IF ( n_pi_predictors /= n_tc_predictors ) THEN 
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Inconsistent PREDICTOR_INDEX and '//&
                            'TRANSMITTANCE_COEFFICIENTS predictor dimension.', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF

    ! -- Check number of layers
    IF ( n_as_layers /= n_tc_layers ) THEN 
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Inconsistent ABSORBER_SPACE_LEVELS and '//&
                            'TRANSMITTANCE_COEFFICIENTS level dimension.', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF

    ! -- Check number of channels
    IF ( n_pi_channels /= n_tc_channels ) THEN 
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Inconsistent PREDICTOR_INDEX and '//&
                            'TRANSMITTANCE_COEFFICIENTS channel dimension.', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF

    ! -- Check number of absorbers
    IF ( n_alpha_absorbers /= n_tc_absorbers .OR. &
         n_as_absorbers    /= n_tc_absorbers .OR. &
         n_pi_absorbers    /= n_tc_absorbers      ) THEN 
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            'Inconsistent ALPHA/ABSORBER_SPACE_LEVELS/PREDICTOR_INDEX '//&
                            'and TRANSMITTANCE_COEFFICIENTS absorber dimension.', &
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
    !#      -- Open the transmittance coefficient data file for output --       #
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

    WRITE( file_id, IOSTAT = io_status ) n_tc_predictors, &
                                         n_tc_layers,     &
                                         n_tc_channels,   &
                                         n_tc_absorbers

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

    WRITE( file_id, IOSTAT = io_status ) N_TRANSMITTANCE_ITEMS

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


    WRITE( file_id, IOSTAT = io_status ) TRANSMITTANCE_DATA_TYPE

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
    !#                 -- Write the absorber space data --                      #
    !#--------------------------------------------------------------------------#

    ! -- Alpha values
    WRITE( file_id, IOSTAT = io_status ) alpha

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error writing alpha values. IOSTAT = ", i5 )' ) &
                      io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      CLOSE( file_id )
      RETURN
    END IF

    ! -- Absorber space level data
    WRITE( file_id, IOSTAT = io_status ) absorber_space_levels

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error writing absorber space level data. IOSTAT = ", i5 )' ) &
                      io_status
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

    DO l = 1, n_tc_channels


      ! ---------------
      ! Predictor index
      ! ---------------

      WRITE( file_id, IOSTAT = io_status ) predictor_index( :, l, : )

      IF ( io_status /= 0 ) THEN
        error_status = FAILURE
        WRITE( message, '( "Error writing channel ", i4, 1x, &
                          &"predictor indices. IOSTAT = ", i5 )' ) &
                        l, io_status
        CALL display_message( ROUTINE_NAME, &
                              TRIM( message ), &
                              error_status, &
                              message_log = message_log )
        CLOSE( file_id )
        RETURN
      END IF


      ! --------------------------
      ! Transmittance coefficients
      ! --------------------------

      WRITE( file_id, IOSTAT = io_status ) transmittance_coefficients( :, :, l, : )

      IF ( io_status /= 0 ) THEN
        error_status = FAILURE
        WRITE( message, '( "Error writing channel ", i4, 1x, &
                          &"transmittance coefficients. IOSTAT = ", i5 )' ) &
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
!    WRITE( message, '( "FILE VERSION: ", i1, ".", i2.2, 2x, &
!                      &"N_PREDICTORS=",i1,2x,&
!                      &"N_ABSORBER_LAYERS=",i3,2x,&
!                      &"N_CHANNELS=",i4,2x,&
!                      &"N_ABSORBERS=",i1 )' ) &
!                    file_release, file_version, &
!                    n_tc_predictors, &
!                    n_tc_layers, &
!                    n_tc_channels, &
!                    n_tc_absorbers
!    CALL display_message( ROUTINE_NAME, &
!                          TRIM( message ), &
!                          INFORMATION, &
!                          message_log = message_log )

    error_status = SUCCESS

  END FUNCTION write_tau_coefficients

END MODULE transmittance_coefficients


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
! Revision 1.11  2001/08/31 21:14:36  paulv
! - Added MIN and MAX release/version parameters to allow for valid use of
!   data files within a specified range.
! - Added the absorber space exponential parameter ALPHA to the transittance
!   file data. READ_ and WRITE_ functions updated as well as the valid version
!   number.
!
! Revision 1.10  2001/08/16 17:17:54  paulv
! - Updated documentation
! - The comparison of n_channels and MAX_N_CHANNELS is now done via the
!   MAX_N_CHANNELS methods in the PARAMETERS module.
!
! Revision 1.9  2001/08/09 20:47:25  paulv
! - Added the WRITE_TRANSMITTANCE_COEFFICIENTS function.
! - Moved all the transmittance data type and name definitions from the
!   COEFFICIENT_UTILITY module to this one. Altered USE statement of the
!   COEFFICIENT_UTILITY module to reflect this change.
! - Added VALID_RELEASE and VALID_VERSION parameters for data file version
!   checking.
! - Added data file release and version number read/write to the requisite
!   read/write function.
!
! Revision 1.8  2001/08/01 17:07:02  paulv
! - The absorber space levels are no longer calculated during model
!   initialisation, but are precalculated and stored in the transmittance
!   coefficient data file. Thus a shared data array, ABSORBER_SPACE_LEVELS,
!   was added to this module as were array allocation/deallocation and
!   data read statements.
!
! Revision 1.7  2001/07/12 17:49:04  paulv
! - Removed definitions of the number, type, and name of the transmittance items
!   and moved them into the COEFFICIENT_UTILITY module. They are now available
!   via:
!     USE coefficient_utility, ONLY: open_coefficient_file, &
!                                    N_TRANSMITTANCE_ITEMS,      &
!                                    TRANSMITTANCE_DATA_TYPE,    &
!                                    TRANSMITTANCE_DATA_NAME
!   This was done to allow the COEFFICIENT_UTILITY module to be used for
!   reformatting. However, this may change - now definitions for the contents
!   of the transmittance coefficient data file are distributed in two different
!   modules. I don't like that.
!
! Revision 1.6  2001/05/29 17:50:31  paulv
! - Added DESTROY_TRANSMITTANCE_COEFFICIENTS function. This deallocates the
!   data arrays used to store the transmittance coefficient data.
!
! Revision 1.5  2000/11/09 20:42:11  paulv
! - Coefficient arrays are now ALLOCATABLE.
! - Input file format has changed to contain data dimension and type
!   information for file data checks and array allocation.
! - Coefficient data is now read directly into the correct channel position
!   straight from the file. Previosuly a dummy chunk was read and then placed
!   into the shared data array in a separate loop. The new method is faster
!   and a lot easier to follow.
!
! Revision 1.4  2000/08/31 19:36:34  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 1.3  2000/08/24 16:17:02  paulv
! - Dimension order of PREDICTOR_INDEX changed from:
!
!     0:MAX_N_PREDICTORS_TO_USE x MAX_N_ABSORBERS x MAX_N_CHANNELS
!   to:
!     0:MAX_N_PREDICTORS_TO_USE x MAX_N_CHANNELS  x MAX_N_ABSORBERS
!
!   and for TAU_COEFFICIENTS from
!
!     0:MAX_N_PREDICTORS_TO_USE x MAX_N_ABSORBERS x 0:MAX_N_ABSORBER_LAYERS x MAX_N_CHANNELS
!   to
!     0:MAX_N_PREDICTORS_TO_USE x 0:MAX_N_ABSORBER_LAYERS x MAX_N_CHANNELS x MAX_N_ABSORBERS
!
!   This allowed for more efficent access to the various data on an absorber
!   by absorber basis (since if there is no significant channel absorption by
!   any particular absorber go to the next absorber).
! - Removed references to the record length parameter. No longer needed as
!   file access is SEQUENTIAL rather than DIRECT.
! - Replaced error check after OPEN_COEFFICIENT_FILE call with a simple
!   test for error_status /= SUCCESS. The open function no longer returns
!   any error status other than SUCCESS or FAILURE (used to return WARNING
!   in some circumstances.)
! - "REC =" keyword removed from file READ statement.
! - TAU_COEFFICIENTS array filled by looping over the number of absorber layers.
!   This is only marginally faster than using array syntax but it makes the code
!   a bit easier to understand (IMO).
! - Channel loop construct name changed from "channel_loop" to "l_channel_loop"
!   to indicate the loop counter variable is "l". This is not a big deal for
!   this situation but has proven useful in other modules with a high degree
!   of nested loops.
! - Updated module and subprogram documentation.
!
! Revision 1.2  2000/08/08 17:03:38  paulv
! Module modified to:
! - Read the transmittance coefficients correctly! and
! - To use the PARAMETERS module rather than the CONSTANTS module.
!
! Revision 1.1  2000/07/12 16:08:10  paulv
! Initial checked in version
!
!
!
