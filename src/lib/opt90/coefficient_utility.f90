!------------------------------------------------------------------------------
!M+
! NAME:
!       coefficient_utility
!
! PURPOSE:
!       Module to hold utility data and routines for reading and writing RT
!       model coefficient data
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE coefficient_utility
!
! OUTPUTS:
!       None.
!
! MODULES:
!       error_handler:          Module to define error codes and handle error
!                               conditions
!
!       file_utility:   .       Module containing global file utility routines
!
! CONTAINS:
!       open_coefficient_file:  PUBLIC function to open the sequential access
!                               coefficient files.
!
!       check_coefficient_file: PRIVATE function to determine if the coefficient
!                               file is in the correct format, endian-wise.
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
!                       paul.vandelst@ssec.wisc.edu
!                         or
!                       pvandelst@ncep.noaa.gov
!
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

MODULE coefficient_utility


  ! ----------
  ! Module use
  ! ----------

  USE type_kinds
  USE file_utility
  USE error_handler


  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------

  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------

  PRIVATE
  PUBLIC :: open_coefficient_file


  ! -----------------
  ! Module parameters
  ! -----------------

  INTEGER( Long ), PARAMETER, PRIVATE :: MAGIC_NUMBER = 123456789_Long


  ! ----------------
  ! Module variables
  ! ----------------

  CHARACTER( 128 ) :: message


  ! -----------------
  ! Module intrinsics
  ! -----------------

  INTRINSIC PRESENT, &
            TRIM



CONTAINS


!------------------------------------------------------------------------------
!S+
! NAME:
!       open_coefficient_file
!
! PURPOSE:
!       PUBLIC function to open the sequential access coefficient files
!
! CALLING SEQUENCE:
!       result = open_coefficient_file( coefficient_file, &
!                                       file_id, &
!                                       message_log = message_log )
!
! INPUT ARGUMENTS:
!       coefficient_file: Name of the file containing the coefficient data.
!                         UNITS:      None
!                         TYPE:       Character
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( IN )
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
!       file_id:          File logical unit number.
!                         UNITS:      None
!                         TYPE:       Integer
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT( OUT )
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result = SUCCESS => file open was successful
!              = FAILURE => error occurred during file open OR
!                           error occurred during file check
!
! CALLS:
!      file_exists:             Function to determine if a named file exists.
!                               SOURCE: file_utility module
!
!      get_lun:                 Function to return a free logical unit number
!                               for file access.
!                               SOURCE: file_utility module
!
!      display_message:         Subroutine to output messages
!                               SOURCE: error_handler module
!
!      check_coefficient_file:  Function to check the file for the correct
!                               endian-ness.
!                               SOURCE: PRIVATE module subprogram
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
!       The file is inquired to determine if it exists. If so, it is opened
!       and the endian-ness is checked.
!S-
!------------------------------------------------------------------------------

  FUNCTION open_coefficient_file( coefficient_file, &
                                  file_id,          &
                                  for_output,       &
                                  no_check,         &
                                  message_log )     &
                                RESULT( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    CHARACTER( * ),           INTENT( IN )  :: coefficient_file
    INTEGER,                  INTENT( OUT ) :: file_id

    INTEGER,        OPTIONAL, INTENT( IN )  :: for_output
    INTEGER,        OPTIONAL, INTENT( IN )  :: no_check

    CHARACTER( * ), OPTIONAL, INTENT( IN )  :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'OPEN_COEFFICIENT_FILE'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: file_check
    INTEGER :: file_output
    INTEGER :: io_status

    CHARACTER( 7 ) :: file_status
    CHARACTER( 5 ) :: file_action



    !#--------------------------------------------------------------------------#
    !#                      -- Check optional arguments --                      #
    !#--------------------------------------------------------------------------#

    IF ( PRESENT( no_check ) ) THEN
      file_check = 0
    ELSE
      file_check = 1
    END IF

    IF ( PRESENT( for_output ) ) THEN
      file_output = 1
      file_check  = 0
    ELSE
      file_output = 0
    END IF


    !#--------------------------------------------------------------------------#
    !#                      -- Check data file existence --                     #
    !#--------------------------------------------------------------------------#

    IF ( file_output == 0 ) THEN

      ! -- If data file does not exist, return an error
      IF ( .NOT. file_exists( coefficient_file ) ) THEN
        error_status = FAILURE
        WRITE( message, '( "Coefficient file, ", a, " not found." )' ) &
                        TRIM( coefficient_file )
        CALL display_message( ROUTINE_NAMe, &
                              TRIM( message ), &
                              error_status, &
                              message_log = message_log )
        RETURN
      END IF

      ! -- Set OPEN keywords for reading
      file_status = 'OLD   '
      file_action = 'READ '

    ELSE

      ! -- If data file does exist, output a warning message
      IF ( file_exists( coefficient_file ) ) THEN
        WRITE( message, '( "Coefficient file, ", a, " will be overwritten." )' ) &
                        TRIM( coefficient_file )
        CALL display_message( ROUTINE_NAMe, &
                              TRIM( message ), &
                              WARNING, &
                              message_log = message_log )
      END IF

      ! -- Set OPEN keywords for writing
      file_status = 'REPLACE'
      file_action = 'WRITE'

    END IF



    !#--------------------------------------------------------------------------#
    !#                        -- Open the data file --                          #
    !#--------------------------------------------------------------------------#

    file_id = get_lun()
    OPEN( file_id, FILE   = coefficient_file, &
                   STATUS = TRIM( file_status ), &
                   ACTION = TRIM( file_action ), &
                   ACCESS = 'SEQUENTIAL', &
                   FORM   = 'UNFORMATTED', &
                   IOSTAT = io_status )

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error opening ", a, ". IOSTAT = ", i5 )' ) &
                      coefficient_file, io_status
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#        -- If file is opened for output, write the magic number --        #
    !#--------------------------------------------------------------------------#

    IF ( file_output == 1 ) THEN

      WRITE( file_id, IOSTAT = io_status ) MAGIC_NUMBER

      IF ( io_status /= 0 ) THEN
        error_status = FAILURE
        WRITE( message, '( "Error writing magic number to ", a, ". IOSTAT = ", i5 )' ) &
                        TRIM( coefficient_file ), io_status
        CALL display_message( ROUTINE_NAME, &
                              TRIM( message ), &
                              error_status, &
                              message_log = message_log )
        CLOSE( file_id )
        RETURN
      END IF

    END IF


    
    !#--------------------------------------------------------------------------#
    !#             -- Check the coefficient data file if required --            #
    !#--------------------------------------------------------------------------#

    IF ( file_check == 0 ) THEN
      error_status = SUCCESS
      RETURN
    END IF

    error_status = check_coefficient_file( file_id, &
                                           message_log = message_log )

    IF ( error_status /= SUCCESS ) THEN
      CLOSE( file_id )
      WRITE( message, '( "Error checking ", a, ". File closed." )' ) &
                      coefficient_file
      CALL display_message( ROUTINE_NAME, &
                            TRIM( message ), &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF

  END FUNCTION open_coefficient_file



!------------------------------------------------------------------------------
!P+
! NAME:
!       check_coefficient_file
!
! PURPOSE:
!       PRIVATE function to determine if the coefficient file is in the correct
!       format, endian-wise.
!
! CALLING SEQUENCE:
!       result = check_coefficient_file( file_id, &
!                                        message_log = message_log )
!
! INPUT ARGUMENTS:
!       file_id:          File logical unit number for the open file that is
!                         to be checked.
!                         UNITS:      None
!                         TYPE:       Integer
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
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result = SUCCESS => file check was successful
!              = FAILURE => error occurred reading a file record OR
!                           8- and/or 32-bit integers not supported.
!
! CALLS:
!       display_message:  Subroutine to output messages
!                         SOURCE: error_handler module
!
! MODULES:
!       type_kinds:  Module containing data type kind definitions. Only the
!                    Byte and Long kind types, and the definition of the number
!                    of bytes used for each type are used.
!
! CONTAINS:
!       swap_endian_long_integer:  Function to byte-swap a long integer to
!                                  determine if input data byte-swapping may
!                                  be required.
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
!       The file is inquired to determine if it exists. If so, it is opened
!       and the endian-ness is checked.
!P-
!------------------------------------------------------------------------------

  FUNCTION check_coefficient_file( file_id,  &
                                   message_log ) &
                                 RESULT( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    INTEGER,        INTENT( IN )           :: file_id
    CHARACTER( * ), INTENT( IN ), OPTIONAL :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ),  PARAMETER :: ROUTINE_NAME = 'CHECK_COEFFICIENT_FILE'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER         :: io_status
    INTEGER( Long ) :: magic_number_read, dummy_long



    !#--------------------------------------------------------------------------#
    !#             -- Check that the current compilation supports --            #
    !#             -- 1- and 4-byte integer types                 --            #
    !#--------------------------------------------------------------------------#

    IF ( BIT_SIZE( 1_Long ) /= 32 .OR. &
         BIT_SIZE( 1_Byte ) /=  8      ) THEN
      error_status = FAILURE
      CALL display_message( ROUTINE_NAME, &
                            '8- and/or 32-bit integers not supported. '//&
                            'Unable to determine endian-ness', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#                        -- Read the magic number --                       #
    !#--------------------------------------------------------------------------#

    READ( file_id, IOSTAT = io_status ) magic_number_read

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading file. IOSTAT = ", i5 )' ) &
                      io_status
      CALL display_message( ROUTINE_NAME, &
                            message, &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#                      -- Compare the magic numbers --                     #
    !#--------------------------------------------------------------------------#

    IF ( magic_number_read /= MAGIC_NUMBER ) THEN

      ! -- Set the error status
      error_status = FAILURE

      ! -- Byte swap the file data.
      dummy_long = swap_endian_long_integer( magic_number_read )

      ! -- Check the file data again
      IF ( dummy_long /= MAGIC_NUMBER ) THEN
        CALL display_message( ROUTINE_NAME, &
                              'Unrecognised file format. Invalid magic number.', &
                              error_status, &
                              message_log = message_log )
        RETURN
      END IF

      ! -- If we get here then the data does need to be byte-swapped.
      CALL display_message( ROUTINE_NAME, &
                            'Data file needs to be byte-swapped.', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF

    error_status = SUCCESS


  CONTAINS


    !#--------------------------------------------------------------------------#
    !#          -- Internal subprogram to byte-swap a long integer --           #
    !#--------------------------------------------------------------------------#

    FUNCTION swap_endian_long_integer ( input ) &
                                      RESULT ( output )


      ! -----------------
      ! Type declarations
      ! -----------------

      ! -- Argument and result
      INTEGER( Long ), INTENT( IN ) :: input
      INTEGER( Long )               :: output


      ! -- Local variables
      INTEGER( Byte ), DIMENSION( n_bytes_for_Long_kind ) :: binput, boutput
      INTEGER( Long )                                     :: linput, loutput
      INTEGER :: i


      ! -------------------------------------------
      ! Equivalence the byte array and long integer
      ! -------------------------------------------

      EQUIVALENCE ( binput,  linput  ), &
                  ( boutput, loutput )


      ! ----------------------------------------------------------------
      ! Loop over the number of bytes for swapping.
      !
      ! Doing it this way is little bit faster (by about a factor of 4)
      ! than using the MVBITS intrinsic ( on the systems tested; Linux
      ! and AIX):
      !
      !  CALL MVBITS( input, 0,  8, output, 24 )  ! Bits  0-7  --> 24-31
      !  CALL MVBITS( input, 8,  8, output, 16 )  ! Bits  8-15 --> 16-23
      !  CALL MVBITS( input, 16, 8, output,  8 )  ! Bits 16-23 -->  8-15
      !  CALL MVBITS( input, 24, 8, output,  0 )  ! Bits 24-31 -->  0-8
      !
      ! but ONLY if the byte swap loop is inline (rather than by calling
      ! a generic byte swap routine with the number of bytes to swap is
      ! passed as an argument.)
      ! ----------------------------------------------------------------

      ! -- Reassign the input argument. Can't
      ! -- equivalence dummy arguments.
      linput = input

      ! -- Loop over the bytes and swap
      DO i = 1, n_bytes_for_Long_kind
        boutput( i ) = binput( n_bytes_for_Long_kind - ( i - 1 ) )
      END DO

      ! -- Assign the output argument. Can't
      ! -- equivalence dummy arguments.
      output = loutput

    END FUNCTION swap_endian_long_integer

  END FUNCTION check_coefficient_file

END MODULE coefficient_utility


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
! Revision 1.8  2001/08/16 16:39:54  paulv
! - Updated documentation
!
! Revision 1.7  2001/08/09 20:38:15  paulv
! - Changed magic number visibility attribute to PRIVATE. Ahhh....
! - Moved all the spectral and transmittance coefficient data type and name
!   definitions into their respective modules. Another ahhh.....
! - Added optional FOR_OUTPUT argument to OPEN_COEFFICIENT_FILE function
!   so that the same function can be used to open coefficient files for
!   writing. It also means the magic number write can be done in this module
!   totally encapsulating that functionality in this module only. Double ahhh....
!
! Revision 1.6  2001/08/01 16:43:05  paulv
! - Updated the definitions of data items and types in the transmittance
!   coefficient data file to reflect changes in code. The absorber space
!   levels are no longer calculated during model initialisation, but are
!   precalculated and stored in the transmittance coefficient data file.
!
! Revision 1.5  2001/07/12 16:58:18  paulv
! - Added USE of TYPE_KINDS module at top of this module. Previously it was
!   USEd only in the CHECK_COEFFICIENT_FILE() function.
! - Data file magic number definition now defined at top of module rather
!   than in the CHECK_COEFFICIENT_FILE() function.
! - Definitions for the number, type, and names of items in the transmittance
!   and spectral coefficient files moved from the TRANSMITTANCE_COEFFICIENTS
!   and SPECTRAL COEFFICIENTS module to this one. This was done to allow this
!   module to be used in both reading and writing/reformatting the coefficient
!   data files.
! - Module-wide error message character string defined.
! - Added NO_CHECK optional argument to the OPEN_COEFFICIENT_FILE() function.
!   This was done to allow the function to be used to open the old format
!   coefficient files for reformatting by not performing a magic number check.
!
! Revision 1.4  2000/08/31 19:36:31  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 1.3  2000/08/24 15:22:10  paulv
! - File access changed from DIRECT to SEQUENTIAL. Record length argument
!   no longer required by OPEN_COEFFICIENT_FILE and CHECK_COEFFICIENT_FILE
!   subprograms.
! - INQUIRE statement in OPEN_COEFFICIENT_FILE that checks for existence
!   of the file replaced by function FILE_EXISTS in module FILE_UTILITY.
! - CHECK_COEFFICIENT_FILE used to return a WARNING status if either 8- or
!   32-bit integers were not supported. This condition now returns a
!   FAILURE status as the magic number would not be read so any subsequent
!   attempt to read data would either fail or return junk.
! - The name of the SWAP_ENDIAN_FOURBYTE_INTEGER subprogram was changed to
!   SWAP_ENDIAN_LONG_INTEGER to remove any indication of how many bytes are
!   expected for this data type *apart* from the definition of
!   N_BYTES_FOR_LONG_KIND in the TYPE_KINDS module.
! - Updated module and subprogram documentation.
!
! Revision 1.2  2000/08/08 17:05:45  paulv
! Cosmetic changes to highlight parameters in the source by making them
! uppercase.
!
! Revision 1.1  2000/07/12 16:08:10  paulv
! Initial checked in version
!
!
!

