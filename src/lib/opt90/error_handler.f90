!------------------------------------------------------------------------------
!M+
! NAME:
!       error_handler
!
! PURPOSE:
!       Module to define simple error codes and handle error conditions
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE error_handler
!
! OUTPUTS:
!       SUCCESS:     Code specifying successful completion.
!                    UNITS:      None.
!                    TYPE:       Integer
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: PARAMETER, PUBLIC
!
!       INFORMATION: Code specifying information output.
!                    UNITS:      None.
!                    TYPE:       Integer
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: PARAMETER, PUBLIC
!
!       WARNING:     Code specifying warning state. Execution can
!                    continue but results may be incorrect.
!                    UNITS:      None.
!                    TYPE:       Integer
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: PARAMETER, PUBLIC
!
!       FAILURE:     Code specifying severe error. Execution cannot
!                    continue.
!                    UNITS:      None.
!                    TYPE:       Integer
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: PARAMETER, PUBLIC
!
!       UNDEFINED:   Code specifying undefined completion status.
!                    UNITS:      None.
!                    TYPE:       Integer
!                    DIMENSION:  Scalar
!                    ATTRIBUTES: PARAMETER, PUBLIC
!
!
! MODULES:
!       file_utility: Module containing global file utility routines.
!                     Only the get_lun() function is used in this
!                     module.
!
! CONTAINS:
!       display_message:  PUBLIC subroutine to display error/status messages
!                         either to standard output (default) or to a log file.
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
! EXAMPLE:
!       USE error_handler
!       error_status = calculate_widget_size()
!       IF ( error_status /= SUCCESS ) THEN
!         CALL display_message( routine_name, &
!                               'Error calculating widget size', &
!                               error_status, &
!                               message_log = 'error_log.txt' )
!         RETURN
!       END IF
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

MODULE error_handler


  ! ---------------------
  ! Module use statements
  ! ---------------------
  ! FSL Modification (BLS, 26 Nov 02):
  !  Commented out the USE file_utility statement here, because
  !  it causes problems when using PGF90.  It appears to pgf90
  !  that get_lun is defined in two different modules.  This
  !  USE statement is instead placed into the subroutine below
  !  that actually calls "get_lun"
  !USE file_utility, ONLY: get_lun


  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------

  IMPLICIT NONE


  ! ------------------
  ! Default visibility
  ! ------------------

  PRIVATE


  ! ------------------------------------
  ! Definitions of public parameter data
  ! ------------------------------------
 
  ! -- Integer values that define the error state.
  ! -- Note: These values are totally arbitrary. 
  INTEGER, PARAMETER, PUBLIC :: SUCCESS     = 0
  INTEGER, PARAMETER, PUBLIC :: INFORMATION = SUCCESS + 1
  INTEGER, PARAMETER, PUBLIC :: WARNING     = INFORMATION + 1
  INTEGER, PARAMETER, PUBLIC :: FAILURE     = WARNING + 1
  INTEGER, PARAMETER, PUBLIC :: UNDEFINED   = FAILURE + 1


  ! -----------------------------------
  ! Definitions of local parameter data
  ! -----------------------------------

  ! -- Character descriptors of the error states
  INTEGER,         PARAMETER :: MAX_N_STATES = 5
  CHARACTER( 11 ), PARAMETER, DIMENSION( 0:MAX_N_STATES-1 ) :: &
    STATE_DESCRIPTOR = (/ 'SUCCESS    ', &
                          'INFORMATION', &
                          'WARNING    ', &
                          'FAILURE    ', &
                          'UNDEFINED  ' /)


  ! ----------------------------------
  ! Explicit visibility of subprograms
  ! ----------------------------------

  PUBLIC :: display_message


CONTAINS



!------------------------------------------------------------------------------
!S+
! NAME:
!       display_message
!
! PURPOSE:
!       RECURSIVE PUBLIC routine to display messages.
!
! CALLING SEQUENCE:
!       CALL display_message( routine_name, &
!                             message,      &
!                             error_state,  &
!                             message_log  = message_log )
!
! INPUT ARGUMENTS:
!       routine_name: Name of the routine in which the message originated.
!                     UNITS:      None
!                     TYPE:       Character
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT( IN )
!
!       message:      Message text
!                     UNITS:      None
!                     TYPE:       Character
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT( IN )
!
!       error_state:  Flag corresponding to one of the defined error states.
!                     If not, the error state is set to UNDEFINED.
!                     UNITS:      None
!                     TYPE:       Integer
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       message_log:  Character string specifying a filename in which any
!                     messages will be logged. If not specified, or if an
!                     error occurs opening the log file, the default action
!                     is to output messages to the screen.
!                     UNITS:      None
!                     TYPE:       Character
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! CALLS:
!      get_lun:   Function to return a free logical unit number for
!                 file access.
!                 SOURCE: file_utility module
!
!      Routine calls itself if the optional argument message_log is passed and
!      an error occurs opening the output log file.
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
!       Output message format is:
!
!         "routine name"("state description") : "message"
!
!       For example, if an error occurs in this routine the output is:
!
!         "DISPLAY_MESSAGE(FAILURE) : Error opening message log file"
!S-
!------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE display_message ( routine_name, &
                                         message,      &
                                         error_state,  &
                                         message_log   )


    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#
    ! FSL Modification (BLS, 26 Nov 02)
    !  Added the USE file_utility statement originally at top of
    !  module to this private subroutine.
    USE file_utility, ONLY: get_lun
    ! ---------
    ! Arguments
    ! ---------

    CHARACTER( * ), INTENT( IN )           :: routine_name
    CHARACTER( * ), INTENT( IN )           :: message
    INTEGER,        INTENT( IN )           :: error_state
    CHARACTER( * ), INTENT( IN ), OPTIONAL :: message_log


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: THIS_ROUTINE_NAME = 'DISPLAY_MESSAGE'


    ! ---------------
    ! Local variables
    ! ---------------

    INTEGER :: error_state_to_use
    INTEGER :: log_to_file
    INTEGER :: file_id
    INTEGER :: io_status

    CHARACTER( 28 ) :: fmt_string


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC PRESENT, &
              TRIM


    !#--------------------------------------------------------------------------#
    !#                   -- Check the input error state --                      #
    !#--------------------------------------------------------------------------#

    error_state_to_use = error_state
    IF ( error_state < 0 .OR. error_state > MAX_N_STATES ) THEN
      error_state_to_use = UNDEFINED
    END IF



    !#--------------------------------------------------------------------------#
    !#      -- Set the message log. If not specified, output to screen --       #
    !#--------------------------------------------------------------------------#

    IF ( PRESENT( message_log ) ) THEN

      log_to_file = 1
      file_id     = get_lun()

      OPEN( file_id, FILE     = message_log,  &
                     ACCESS   = 'SEQUENTIAL', &
                     FORM     = 'FORMATTED',  &
                     STATUS   = 'UNKNOWN',    &
                     POSITION = 'APPEND',     &
                     ACTION   = 'READWRITE',  & ! Just READ may cause probs on some
                                                ! systems using POSITION = 'APPEND'
                     IOSTAT   = io_status )

      IF ( io_status /= 0 ) THEN
        CALL display_message( THIS_ROUTINE_NAME, &
                              'Error opening message log file', &
                              FAILURE )
        log_to_file = 0
      END IF

    ELSE

      log_to_file = 0

    END IF


    !#--------------------------------------------------------------------------#
    !#                         -- Output the message --                         #
    !#--------------------------------------------------------------------------#

    fmt_string = '( 1x, a, "(", a, ") : ", a )'

    log_message: IF ( log_to_file == 0 ) THEN
      WRITE( *, FMT = fmt_string ) &
                TRIM( routine_name ), &
                TRIM( STATE_DESCRIPTOR( error_state_to_use ) ), &
                TRIM( message )
    ELSE
      WRITE( file_id, FMT = fmt_string ) &
                      TRIM( routine_name ), &
                      TRIM( STATE_DESCRIPTOR( error_state_to_use ) ), &
                      TRIM( message )
      CLOSE( file_id )
    END IF log_message

  END SUBROUTINE display_message

END MODULE error_handler


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
! Revision 1.1  2002/11/15 15:21:31  birk
! Added to cvs mainly to see how this compiles on other platforms, it currently
! seems to compile on the IBM
!
! Revision 1.3  2000/08/31 19:36:32  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 1.2  2000/08/24 15:27:18  paulv
! - The DISPLAY_MESSAGE subprogram was made RECURSIVE so it can call itself
!   if an error occurs opening the message log file defined by the optional
!   input argument MESSAGE_LOG.
! - The message log file is now closed after the message is written (as it
!   should have always been...oops).
! - Updated module and subprogram documentation.
!
! Revision 1.1  2000/07/12 16:08:10  paulv
! Initial checked in version
!
!
!
