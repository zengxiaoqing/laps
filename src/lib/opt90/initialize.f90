!------------------------------------------------------------------------------
!M+
! NAME:
!       initialize
!
! PURPOSE:
!       Module for RT model initialisation. The absorber space layer values
!       are calculated and the transmittance and spectral coefficient data
!       are read from their respective data files.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       USE initialize
!
! OUTPUTS:
!       None.
!
! MODULES:
!       error_handler:               Module to define error codes and handle
!                                    error conditions.
!
!       spectral_coefficients:       Module containing the spectral coefficient
!                                    data and read methods.
!
!       transmittance_coefficients:  Module containing the transmittance
!                                    coefficient data and read methods.
!
! CONTAINS:
!       initialize_rtm:   PUBLIC function to initialise the RT model.
!
!       destroy_rtm:      PUBLIC function to destroy the RT model space.
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       Absorber space, transmittance coefficient, and spectral coefficient
!       data arrays are filled.
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
!       Written by:     Paul van Delst, CIMSS@NOAA/NCEP 31-Jul-2000
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

MODULE initialize

  ! ----------
  ! Module use
  ! ----------

  USE error_handler
  USE transmittance_coefficients
  USE spectral_coefficients


  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------

  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------

  PRIVATE
  PUBLIC :: initialize_rtm
  PUBLIC :: destroy_rtm


CONTAINS


!------------------------------------------------------------------------------
!S+
! NAME:
!       initialize_rtm
!
! PURPOSE:
!       PUBLIC function to initialise all the data arrays required by the
!       RT model.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       result = initialize_rtm( tau_file      = tau_file,      &
!                                spectral_file = spectral_file, &
!                                path          = path,          &
!                                message_log   = message_log    )
!
! INPUT ARGUMENTS:
!       None.
!
! OPTIONAL INPUT ARGUMENTS:
!       tau_file:      Character string specifying a file name for the
!                      transmittance coefficient data file. If not
!                      specified, "transmittance_coefficients" is the
!                      default.
!                      UNITS:      None
!                      TYPE:       Character
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT( IN ), OPTIONAL
!
!       spectral_file: Character string specifying a file name for the
!                      spectral coefficient data file. If not
!                      specified, "spectral_coefficients" is the
!                      default.
!                      UNITS:      None
!                      TYPE:       Character
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT( IN ), OPTIONAL
!
!       path:          Character string specifying a file path for the
!                      input coefficient files. If not specified, the
!                      current directory is the default.
!                      UNITS:      None
!                      TYPE:       Character
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT( IN ), OPTIONAL
!
!       message_log:   Character string specifying a filename in which any
!                      messages will be logged. If not specified, or if an
!                      error occurs opening the log file, the default action
!                      is to output messages to the screen.
!                      UNITS:      None
!                      TYPE:       Character
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT( IN ), OPTIONAL
!
! OUTPUT ARGUMENTS:
!       None.
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result = SUCCESS => initialisation was successful
!              = FAILURE => error occurred during initialisation
!
! CALLS:
!       display_message:             Subroutine to output messages
!                                    SOURCE: error_handler module
!
!       read_tau_coefficients:       Function to read the transmittance model
!                                    coefficients and predictor indices, and
!                                    absorber space definitions.
!                                    SOURCE: transmittance_coefficients module
!
!       read_spectral_coefficients:  Function to read the spectral coefficients
!                                    for the various satellite/channels.
!                                    SOURCE: spectral_coefficients module
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       All public data arrays accessed by this module and its dependencies
!       are overwritten.
!
! RESTRICTIONS:
!       If specified, the length of the combined path and filename strings
!       cannot exceed 255 characters.
!
!S-
!------------------------------------------------------------------------------

  FUNCTION initialize_rtm( tau_file,      &
                           spectral_file, &
                           path,          &
                           message_log    ) &
                         RESULT ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    CHARACTER( * ), OPTIONAL, INTENT( IN ) :: tau_file
    CHARACTER( * ), OPTIONAL, INTENT( IN ) :: spectral_file
    CHARACTER( * ), OPTIONAL, INTENT( IN ) :: path
    CHARACTER( * ), OPTIONAL, INTENT( IN ) :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'INITIALIZE_RTM'


    ! ---------------
    ! Local variables
    ! ---------------

    CHARACTER( 255 ) :: transmittance_coefficient_file
    CHARACTER( 255 ) :: spectral_coefficient_file


    ! ----------
    ! Intrinsics
    ! ----------

    INTRINSIC ADJUSTL, &
              PRESENT, &
              TRIM



    !#--------------------------------------------------------------------------#
    !#           -- Read the RT model transmittance coefficients --             #
    !#--------------------------------------------------------------------------#

    ! -- Construct filename
    transmittance_coefficient_file = 'transmittance_coefficients'

    IF ( PRESENT( tau_file ) ) &
      transmittance_coefficient_file = tau_file(1:len_trim(tau_file))

    IF ( PRESENT( path ) ) &
      transmittance_coefficient_file = path(1:len_trim(path)) // &
                                       transmittance_coefficient_file


    ! -- Read the data file
    error_status = read_tau_coefficients( &
         transmittance_coefficient_file(1:len_trim(transmittance_coefficient_file)) , &
         message_log = message_log        &
         )

    IF ( error_status /= SUCCESS ) THEN
      CALL display_message( ROUTINE_NAME, &
                            'Error reading transmittance model coefficients', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#               -- Read the RT model spectral coefficients --              #
    !#--------------------------------------------------------------------------#

    ! -- Construct filename
    spectral_coefficient_file = 'spectral_coefficients'

    IF ( PRESENT( spectral_file ) ) &
      spectral_coefficient_file = spectral_file(1:len_trim(spectral_file))

    IF ( PRESENT( path ) ) &
      spectral_coefficient_file = path(1:len_trim(path)) // &
                                  spectral_coefficient_file


    ! -- Read the data file
    error_status = read_spectral_coefficients( &
          spectral_coefficient_file(1:len_trim(spectral_coefficient_file)), &
                                               message_log = message_log          )

    IF ( error_status /= SUCCESS ) THEN
      CALL display_message( ROUTINE_NAME, &
                            'Error reading spectral coefficients', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF


    !#--------------------------------------------------------------------------#
    !#                               -- Done --                                 #
    !#--------------------------------------------------------------------------#

    error_status = SUCCESS

  END FUNCTION initialize_rtm


!------------------------------------------------------------------------------
!S+
! NAME:
!       destroy_rtm
!
! PURPOSE:
!       PUBLIC function to destroy all the public allocated data arrays
!       creating during initialisation of the RT model.
!
! CATEGORY:
!       NCEP RTM
!
! CALLING SEQUENCE:
!       result = destroy_rtm( message_log = message_log )
!
! INPUT ARGUMENTS:
!       None.
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
! OUTPUT ARGUMENTS:
!       None.
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result = SUCCESS => deallocations were successful
!              = FAILURE => error occurred during deallocations
!
! CALLS:
!       display_message:                Subroutine to output messages
!                                       SOURCE: error_handler module
!
!       destroy_tau_coefficients:       Function to deallocate the transmittance
!                                       model PUBLIC arrays.
!                                       SOURCE: transmittance_coefficients module
!
!       destroy_spectral_coefficients:  Function to deallocate the spectral
!                                       coefficient PUBLIC arrays.
!                                       SOURCE: spectral_coefficients module
!
! EXTERNALS:
!       None
!
! COMMON BLOCKS:
!       None.
!
! SIDE EFFECTS:
!       All public data arrays are deallocated.
!
! RESTRICTIONS:
!       None.
!
!S-
!------------------------------------------------------------------------------

  FUNCTION destroy_rtm( message_log ) RESULT ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- Type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! Arguments
    ! ---------

    CHARACTER( * ), OPTIONAL, INTENT( IN ) :: message_log


    ! ---------------
    ! Function result
    ! ---------------

    INTEGER :: error_status


    ! ----------------
    ! Local parameters
    ! ----------------

    CHARACTER( * ), PARAMETER :: ROUTINE_NAME = 'DESTROY_RTM'



    !#--------------------------------------------------------------------------#
    !#      -- Deallocate the RT model transmittance coefficient arrays --      #
    !#--------------------------------------------------------------------------#

    error_status = destroy_tau_coefficients( message_log = message_log )

    IF ( error_status /= SUCCESS ) THEN
      CALL display_message( ROUTINE_NAME, &
                            'Error deallocating transmittance coefficients', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF



    !#--------------------------------------------------------------------------#
    !#         -- Deallocate the RT model spectral coefficient arrays --        #
    !#--------------------------------------------------------------------------#

    error_status = destroy_spectral_coefficients( message_log = message_log )

    IF ( error_status /= SUCCESS ) THEN
      CALL display_message( ROUTINE_NAME, &
                            'Error deallocating spectral coefficients', &
                            error_status, &
                            message_log = message_log )
      RETURN
    END IF


    !#--------------------------------------------------------------------------#
    !#                               -- Done --                                 #
    !#--------------------------------------------------------------------------#

    error_status = SUCCESS

  END FUNCTION destroy_rtm

END MODULE initialize


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
! Revision 1.6  2001/08/31 21:17:25  paulv
! - Add TAU_FILE, SPECTRAL_FILE, and PATH optional arguments to INITIALIZE_RTM
!   to allow users to specify alternate file names and locations.
!
! Revision 1.5  2001/08/16 16:40:56  paulv
! - Updated documentation
!
! Revision 1.4  2001/08/01 16:51:52  paulv
! - Removed USE of module ABSORBER_SPACE to reflect changes in code.
!   The absorber space levels are no longer calculated during model
!   initialisation, but are precalculated and stored in the transmittance
!   coefficient data file.
! - Removed COMPUTE_ABSORBER_SPACE() function call for same reason.
!
! Revision 1.3  2001/05/29 17:40:22  paulv
! - Changed name of initialisation routine from RTM_INITIALIZE to INITIALIZE_RTM.
! - Added DESTROY_RTM function to deallocate memory used in coefficient read.
!
! Revision 1.2  2000/08/31 19:36:32  paulv
! - Added documentation delimiters.
! - Updated documentation headers.
!
! Revision 1.1  2000/08/08 16:39:35  paulv
! Initial checkin
!
!
!
