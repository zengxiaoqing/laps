!------------------------------------------------------------------------------
!
! NAME:
!       rtm_test_FWD
!
! PURPOSE:
!       Program to test the FWD RTM routine
!
! CATEGORY:
!       NCEP RTM
!
! LANGUAGE:
!       Fortran-90
!
! MODULES:
!       type_kinds:            Module to hold specification kinds for variable
!                              declaration.
!
!       file_utility:          Module containing generic file utility routines
!
!       error_handler:         Module to define simple error codes and handle
!                              error conditions
!
!       parameters:            Module to hold RT model parameter constants
!
!       initialize:            Module for RT model initialisation.
!
!       forward_model:         Module containing the RT forward model function.
!
!       spectral_coefficients: Module to hold the RT model spectral coefficients
!                              and their access routines.
!                              NOTE: This routine is USEd only to access a single
!                                    data entity for test purposes. Does not need 
!                                    to be USEd to call the RT routines.
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
!       Input:  Binary, direct access profile data file.
!       Output: ASCII forward model output file.
!
! SIDE EFFECTS:
!       Output file(s) are overwritten.
!
! RESTRICTIONS:
!       None.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 20-Aug-2001
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
!
!------------------------------------------------------------------------------


PROGRAM rtm_test_FWD


  !#----------------------------------------------------------------------------#
  !#                             -- MODULE USAGE --                             #
  !#----------------------------------------------------------------------------#

  USE type_kinds
  USE file_utility
  USE error_handler
  USE parameters
  USE initialize
  USE forward_model, ONLY: compute_rtm

  ! -- For access to IS_MICROWAVE_CHANNEL only.
  USE spectral_coefficients



  !#----------------------------------------------------------------------------#
  !#                           -- TYPE DECLARATIONS --                          #
  !#----------------------------------------------------------------------------#

  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------

  IMPLICIT NONE


  ! ------------------
  ! Program parameters
  ! ------------------

  ! -- Name
  CHARACTER( * ),  PARAMETER :: PROGRAM_NAME = 'RTM_TEST_FWD'

  ! -- Angle definitions
  INTEGER,         PARAMETER :: N_ANGLES  = 11
  REAL( fp_kind ), PARAMETER :: MIN_ANGLE =  0.0_fp_kind
  REAL( fp_kind ), PARAMETER :: MAX_ANGLE = 60.0_fp_kind

  ! -- Profile definitions
  CHARACTER( * ),  PARAMETER :: PROFILE_FILE        = 'profiles_LAYER.bin'
  INTEGER,         PARAMETER :: N_PROFILES          = 20
  INTEGER,         PARAMETER :: N_LAYERS            = 100

  ! -- Other dimension parameters
  INTEGER,         PARAMETER :: N_ABSORBERS  = MAX_N_ABSORBERS
  INTEGER,         PARAMETER :: N_PREDICTORS = MAX_N_PREDICTORS

  ! -- Emissivity parameters
  REAL( fp_kind ), PARAMETER :: DEFAULT_MW_EMISSIVITY = 0.6_fp_kind
  REAL( fp_kind ), PARAMETER :: DEFAULT_IR_EMISSIVITY = 0.96_fp_kind

  ! -- Default solar angle secant (> 11.47 means no solar)
  REAL( fp_kind ), PARAMETER :: DEFAULT_SECANT_SOLAR_ANGLE = 12.0_fp_kind


  ! -----------------
  ! Program variables
  ! -----------------

  ! -- Error message
  CHARACTER( 128 ) :: message 

  ! -- Status variables
  INTEGER :: error_status
  INTEGER :: allocate_status
  INTEGER :: io_status

  ! -- Loop counters, dimensions and file lun
  INTEGER :: i, j, k, l, m, il
  INTEGER :: n_available_channels, l1, l2, begin_channel, end_channel, n_channels
  INTEGER :: record_length, record_number
  INTEGER :: lun
  INTEGER :: lun_FWD

  ! -- Secant angle arrays
  INTEGER         :: i_angle
  REAL( fp_kind ) :: d_angle
  REAL( fp_kind ), DIMENSION( N_ANGLES ) :: view_angle
  REAL( fp_kind ), DIMENSION( N_ANGLES ) :: secant_view_angle
  REAL( fp_kind ), DIMENSION( N_ANGLES ) :: secant_solar_angle

  ! -- Profile read array
  REAL( fp_kind ), DIMENSION( N_LAYERS ) :: level_pressure
  REAL( fp_kind ), DIMENSION( N_LAYERS ) :: layer_pressure
  REAL( fp_kind ), DIMENSION( N_LAYERS ) :: layer_temperature
  REAL( fp_kind ), DIMENSION( N_LAYERS ) :: layer_water_vapor
  REAL( fp_kind ), DIMENSION( N_LAYERS ) :: layer_ozone

  ! -- Profile data arrays
  REAL( fp_kind ), DIMENSION( N_LAYERS, N_ANGLES ) :: level_p
  REAL( fp_kind ), DIMENSION( N_LAYERS, N_ANGLES ) :: layer_p
  REAL( fp_kind ), DIMENSION( N_LAYERS, N_ANGLES ) :: layer_t
  REAL( fp_kind ), DIMENSION( N_LAYERS, N_ANGLES ) :: layer_w
  REAL( fp_kind ), DIMENSION( N_LAYERS, N_ANGLES ) :: layer_o

  ! -- Number of channels processed for each profile
  INTEGER, DIMENSION( N_ANGLES ) :: n_channels_per_profile

  ! -- Allocatable arrays
  REAL( fp_kind ), DIMENSION( : ),    ALLOCATABLE :: surface_emissivity
  REAL( fp_kind ), DIMENSION( : ),    ALLOCATABLE :: surface_reflectivity
  INTEGER,         DIMENSION( : ),    ALLOCATABLE :: channel_index

  REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: tau
  REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: flux_tau
  REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: solar_tau

  REAL( fp_kind ), DIMENSION( : ),    ALLOCATABLE :: upwelling_radiance
  REAL( fp_kind ), DIMENSION( : ),    ALLOCATABLE :: brightness_temperature

  ! -- Stuff for emissivities
  REAL( fp_kind ) :: angle_modifier
  REAL( fp_kind ), DIMENSION( N_ANGLES ) :: mw_emissivity
  REAL( fp_kind ), DIMENSION( N_ANGLES ) :: ir_emissivity


  ! ----------
  ! Intrinsics
  ! ----------

  INTRINSIC COS, &
            MIN, MAX, &
            REAL, &
            TRIM



  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#                         -- INITIALIZE THE RTM --                           #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################

  WRITE( *, '( /5x, "Initialising RTM...." )' )

  error_status = initialize_rtm( tau_file = 'TEST_transmittance_coefficients', &
                                 spectral_file = 'TEST_spectral_coefficients'  )

  IF ( error_status /= SUCCESS ) THEN
    CALL display_message( PROGRAM_NAME, &
                          'Error initialzing RTM', &
                          error_status )
    STOP
  END IF




  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#               -- SETUP DATA AND ARRAYS FOR CALCULATIONS --                 #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################

  !#----------------------------------------------------------------------------#
  !#                    -- ALLOCATE ARRAYS FOR RTM MODEL --                     #
  !#----------------------------------------------------------------------------#

  ! ----------------------------------
  ! Get the number of channels read in
  ! from the coefficient data files
  ! ----------------------------------

  CALL get_max_n_channels( n_channels )

  begin_channel = 1
  end_channel   = n_channels

!!!Replace the above code with the following to process channel subsets
!!!
!!!  CALL get_max_n_channels( n_available_channels )
!!!
!!!  ! -- Get user channel limits
!!!  WRITE( *, '( /5x, "Number of available channels: ", i4 )' ) n_available_channels
!!!  WRITE( *, '( /5x, "Enter begin,end channels to process: " )', &
!!!            ADVANCE = 'NO' )
!!!  READ( *, * ) l1, l2
!!!
!!!  ! -- Only keep valid values
!!!  l1 = MAX( 1, MIN( l1, n_available_channels ) )
!!!  l2 = MAX( 1, MIN( l2, n_available_channels ) )
!!!  begin_channel = MIN( l1, l2 )
!!!  end_channel   = MAX( l1, l2 )
!!!  n_channels    = end_channel - begin_channel + 1

  WRITE( *, '( /5x, "Processing channels ", i4, " to ", i4, "...", / )' ) &
            begin_channel, end_channel


  ! ---------------------------------------------------
  ! Allocate the Forward model channel dependent arrays
  ! ---------------------------------------------------

  ALLOCATE( surface_emissivity(   n_channels * N_ANGLES ),   &  ! Input,  L*M
            surface_reflectivity( n_channels * N_ANGLES ),   &  ! Input,  L*M
            channel_index(        n_channels * N_ANGLES ),   &  ! Input,  L*M

            tau(       N_LAYERS, n_channels * N_ANGLES ),    &  ! Output, K x L*M
            flux_tau(  N_LAYERS, n_channels * N_ANGLES ),    &  ! Output, K x L*M
            solar_tau( N_LAYERS, n_channels * N_ANGLES ),    &  ! Output, K x L*M

            upwelling_radiance(     n_channels * N_ANGLES ), &  ! Output, L*M
            brightness_temperature( n_channels * N_ANGLES ), &  ! Output, L*M

            STAT = allocate_status )

  IF ( allocate_status /= 0 ) THEN
    WRITE( message, '( "Error allocating forward model channel ", &
                      &"dependent arrays. STAT = ", i5 )' ) &
                    allocate_status
    CALL display_message( PROGRAM_NAME,    &
                          TRIM( message ), &
                          FAILURE          )
    STOP
  END IF



  !#----------------------------------------------------------------------------#
  !#                 -- FILL PROFILE INDEPENDENT INPUT ARRAYS --                #
  !#----------------------------------------------------------------------------#

  ! -- Compute delta angle value
  d_angle = ( MAX_ANGLE - MIN_ANGLE ) / REAL( N_ANGLES-1, fp_kind )

  ! -- Initialise angle x channel counter
  il = 0


  ! ----------------
  ! Loop over angles
  ! ----------------

  DO i_angle = 1, N_ANGLES

    ! -- Fill angle arrays
    view_angle( i_angle )         = MIN_ANGLE + ( REAL( i_angle-1, fp_kind ) * d_angle )
    secant_view_angle( i_angle )  = ONE / COS( view_angle( i_angle ) * DEGREES_TO_RADIANS )
    secant_solar_angle( i_angle ) = DEFAULT_SECANT_SOLAR_ANGLE

    ! -- Fill channel count array, i.e. do them all
    n_channels_per_profile( i_angle ) = n_channels


    ! ------------------
    ! Loop over channels
    ! ------------------

    DO l = 1, n_channels

      ! -- Increment angle x channel counter
      ! -- and set the channel index
      il = il + 1
      channel_index( il ) = begin_channel + l - 1

      ! -- Assign pretend surface emissivity and reflectivity
      ! -- For uW, r = specular; for IR, r = isotropic
      ! -- The angle modifier is just something to provide
      ! -- a little bit of angular variation in the surface
      ! -- emissivities and reflectivities.
      angle_modifier = (COS( view_angle( i_angle ) * DEGREES_TO_RADIANS ))**(0.1_fp_kind)
      IF ( is_microwave_channel( channel_index( il ) ) == 1 ) THEN
        surface_emissivity( il )   = DEFAULT_MW_EMISSIVITY * angle_modifier
        surface_reflectivity( il ) = ONE - surface_emissivity( il )
        mw_emissivity( i_angle ) = surface_emissivity( il )
      ELSE
        surface_emissivity( il )   = DEFAULT_IR_EMISSIVITY * angle_modifier
        surface_reflectivity( il ) = ( ONE - surface_emissivity( il ) ) / PI
        ir_emissivity( i_angle ) = surface_emissivity( il )
      END IF

    END DO

  END DO



  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#                   -- OPEN THE INPUT PROFILE DATA FILE --                   #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################

  ! ------------------------------
  ! Set the record length in bytes
  ! ------------------------------

  record_length = N_LAYERS * n_bytes_for_fp_kind


  ! ---------
  ! Open file
  ! ---------

  lun = get_lun()
  OPEN( lun, FILE   = PROFILE_FILE,  &
             STATUS = 'OLD',         &
             FORM   = 'UNFORMATTED', &
             ACCESS = 'DIRECT',      &
             RECL   = record_length, &
             IOSTAT = io_status      )

  IF ( io_status /= 0 ) THEN
    error_status = FAILURE
    CALL display_message( PROGRAM_NAME, &
                          'Error occurred opening file '//PROFILE_FILE, &
                          error_status )
    STOP
  END IF


  ! ---------------------------------
  ! Read the interface pressure array
  ! ---------------------------------

  record_number = 1
  READ( lun, REC    = record_number, &
             IOSTAT = io_status      ) level_pressure

  IF ( io_status /= 0 ) THEN
    error_status = FAILURE
    CALL display_message( PROGRAM_NAME, &
                          'Error occurred reading interface pressure from '//&
                            PROFILE_FILE, &
                          error_status )
    STOP
  END IF



  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#                     -- OPEN THE OUTPUT DATA FILES --                       #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################

  lun_FWD = get_lun()

  OPEN( lun_FWD, FILE   = 'rtm_test.FWD.output', &
                 STATUS = 'REPLACE',             &
                 FORM   = 'FORMATTED',           &
                 ACCESS = 'SEQUENTIAL',          &
                 IOSTAT = io_status              )

  IF ( io_status /= 0 ) THEN
    error_status = FAILURE
    CALL display_message( PROGRAM_NAME, &
                          'Error occurred opening file rtm_test.FWD.output', &
                          error_status )
    STOP
  END IF


  ! -- Write some header information
  WRITE( lun_FWD, FMT = 100 ) 'MW', mw_emissivity
  WRITE( lun_FWD, FMT = 100 ) 'IR', ir_emissivity
  
  ! -- Output the requisite dimensions
  WRITE( lun_FWD, FMT = 110 ) N_PROFILES, N_ANGLES, N_CHANNELS

  ! -- Write a title
  WRITE( lun_FWD, FMT = 120 )

  ! -- Output format statements
  100 FORMAT( '! ', a, ' surface emissivity (fn. of angle) = ', 20( 2x, f5.3, : ) )
  110 FORMAT( 3( 1x, i5 ) )
  120 FORMAT( '  PROF   ANG    CH  ANGLE  RADIANCE     TEMP.', /, &
             &'---------------------------------------------' )



  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#                       -- BEGIN LOOP OVER PROFILES --                       #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################


  m_profile_loop: DO m = 1, N_PROFILES

    WRITE( *, '( 5x, "Processing profile: ", i2 )' ) m



    !#--------------------------------------------------------------------------#
    !#                           -- READ A PROFILE --                           #
    !#--------------------------------------------------------------------------#

    ! -- Layer pressure
    record_number = record_number + 1
    READ( lun, REC    = record_number, &
               IOSTAT = io_status      ) layer_pressure

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading layer pressure for profile #", i2.2, &
                        &" from ", a, ". IOSTAT = ", i5 )' ) &
                      m, PROFILE_FILE, io_status
      CALL display_message( PROGRAM_NAME,    &
                            TRIM( message ), &
                            error_status     )
      STOP
    END IF


    ! -- Layer temperature
    record_number = record_number + 1
    READ( lun, REC    = record_number, &
               IOSTAT = io_status      ) layer_temperature

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading layer temperature for profile #", i2.2, &
                        &" from ", a, ". IOSTAT = ", i5 )' ) &
                      m, PROFILE_FILE, io_status
      CALL display_message( PROGRAM_NAME,    &
                            TRIM( message ), &
                            error_status     )
      STOP
    END IF


    ! -- Layer water vapor
    record_number = record_number + 1
    READ( lun, REC    = record_number, &
               IOSTAT = io_status      ) layer_water_vapor

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading layer water vapor for profile #", i2.2, &
                        &" from ", a, ". IOSTAT = ", i5 )' ) &
                      m, PROFILE_FILE, io_status
      CALL display_message( PROGRAM_NAME,    &
                            TRIM( message ), &
                            error_status     )
      STOP
    END IF


    ! -- Layer ozone
    record_number = record_number + 1
    READ( lun, REC    = record_number, &
               IOSTAT = io_status      ) layer_ozone

    IF ( io_status /= 0 ) THEN
      error_status = FAILURE
      WRITE( message, '( "Error reading layer ozone for profile #", i2.2, &
                        &" from ", a, ". IOSTAT = ", i5 )' ) &
                      m, PROFILE_FILE, io_status
      CALL display_message( PROGRAM_NAME,    &
                            TRIM( message ), &
                            error_status     )
      STOP
    END IF


    ! ---------------------------------------------------------
    ! Load RTM input profile arrays
    ! (Using same profile for all angles - not efficient since
    !  all quantities are recalculated for each angle, but, eh,
    !  it's a test.)
    ! ---------------------------------------------------------

    DO i_angle = 1, N_ANGLES

      level_p( :, i_angle ) = level_pressure
      layer_p( :, i_angle ) = layer_pressure
      layer_t( :, i_angle ) = layer_temperature
      layer_w( :, i_angle ) = layer_water_vapor
      layer_o( :, i_angle ) = layer_ozone

    END DO


    !#--------------------------------------------------------------------------#
    !#                            -- FORWARD MODEL --                           #
    !#--------------------------------------------------------------------------#

    error_status = compute_rtm( &
                                ! -- Forward inputs
                                level_p,                &  ! Input,  K x M
                                layer_p,                &  ! Input,  K x M
                                layer_t,                &  ! Input,  K x M
                                layer_w,                &  ! Input,  K x M
                                layer_o,                &  ! Input,  K x M

                                layer_t(N_LAYERS,:),    &  ! Input, M
                                surface_emissivity,     &  ! Input, L*M
                                surface_reflectivity,   &  ! Input, L*M

                                secant_view_angle,      &  ! Input, M
                                secant_solar_angle,     &  ! Input, M
                                n_channels_per_profile, &  ! Input, M
                                channel_index,          &  ! Input, L*M

                                ! -- Forward outputs
                                tau,                    &  ! Output, K x L*M
                                flux_tau,               &  ! Output, K x L*M
                                solar_tau,              &  ! Output, K x L*M

                                upwelling_radiance,     &  ! Output, L*M
                                brightness_temperature  )  ! Output, L*M

    IF ( error_status /= SUCCESS ) THEN
      CALL display_message( PROGRAM_NAME, &
                            'Error occured in COMPUTE_RTM', &
                            error_status )
      STOP
    END IF



    ! --------------------------------
    ! Output the forward model results
    ! --------------------------------

    ! -- Initialise the channel/profile counter
    il = 0

    ! -- Loop over the number of ANGLES
    DO i_angle = 1, N_ANGLES

      ! -- Loop over the channels
      DO l = 1, n_channels_per_profile( i_angle )

        ! -- Increment the channel counter and output the
        ! -- channel radiance, and brightness temperature
        il = il + 1
        WRITE( lun_FWD, '( 3( 1x, i5 ), 1x, f6.3, 2( 1x, f9.5 ) )' ) &
                        m, i_angle, l, &
                        view_angle( i_angle ), &
                        upwelling_radiance( il ), &
                        brightness_temperature( il )

      END DO

    END DO


  END DO m_profile_loop

  CLOSE( lun_FWD )
  CLOSE( lun )



  !#----------------------------------------------------------------------------#
  !#                           -- DESTROY THE RTM --                            #
  !#----------------------------------------------------------------------------#

  ! ---------------------------------------
  ! Deallocate the channel dependent arrays
  ! ---------------------------------------

  ! -- Forward model
  DEALLOCATE( surface_emissivity,     &  ! Input,  L*M
              surface_reflectivity,   &  ! Input,  L*M
              channel_index,          &  ! Input,  L*M

              tau,                    &  ! Output, K x L*M
              flux_tau,               &  ! Output, K x L*M
              solar_tau,              &  ! Output, K x L*M

              upwelling_radiance,     &  ! Output, L*M
              brightness_temperature, &  ! Output, L*M

              STAT = allocate_status  )

  IF ( allocate_status /= 0 ) THEN
    WRITE( message, '( "Error deallocating forward model channel ", &
                      &"dependent arrays. STAT = ", i5 )' ) &
                    allocate_status
    CALL display_message( PROGRAM_NAME,    &
                          TRIM( message ), &
                          WARNING          )
  END IF


  ! ---------------------------------
  ! Deallocate the coefficient arrays
  ! ---------------------------------

  WRITE( *, '( /5x, "Destroying RTM...." )' )

  error_status = destroy_rtm()

  IF ( error_status /= SUCCESS ) THEN
    CALL display_message( PROGRAM_NAME, &
                          'Error destroying RTM', &
                          error_status )
    STOP
  END IF

END PROGRAM rtm_test_FWD


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
! Revision 1.1  2001/09/13 22:08:13  paulv
! Initial checkin.
!
!
!
!
