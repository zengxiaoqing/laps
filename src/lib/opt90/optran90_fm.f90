!------------------------------------------------------------------------------
!
! NAME:
!       optran90_fm ! initial interface to f90 code to the LAPS system
!                     built on optran90 test code
!                     and Paul van Delst (acknowledged below)
!                     Forecast systems laboratory
!
!  Copywrite (C) 2003 Daniel Birkenheuer
!
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 20-Aug-2001
!                       paul.vandelst@ssec.wisc.edu
!                       Modified by Dan Birkenheuer
!                       birk@fsl.noaa.gov for interface to LAPS
!                       Used in LAPS per permission of Paul van Delst
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


  subroutine optran90_fm ( kk, &
          lev_p, lay_p, lay_t, lay_w, lay_o, &
          laps_sfc_t, &
          sfc_emis, &
          sfc_refl, &
          sec_za, &
          sec_solar, &
          Mchan, &
          tau90, &
          flx_tau, &
          sol_tau, &
          up_radiance, &
          bright_temp, & 
          sndr_coeff,&
          sndr_trans, &
          sndr_coeff_len, &
          sndr_trans_len &
          )


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

  integer kk
  integer mchan

  real tau90 (kk,mchan)
  real sec_solar
  real flx_tau (kk,mchan)
  real bright_temp (mchan)
  real up_radiance (mchan)
  real sol_tau (kk,mchan)
  real sec_za
  real lay_w (kk)
  real lay_t (kk)
  real lev_p (kk+1)
  real lay_p (kk)
  real lay_o (kk)
  real sfc_refl
  real sfc_emis
  real laps_sfc_t

! optran 90 coefficient data information

  character*256 sndr_coeff, sndr_trans
  integer sndr_coeff_len, sndr_trans_len



  ! ------------------
  ! Program parameters
  ! ------------------

  ! path for LAPS

  character*256 fname
  integer len
  integer channels_used


  ! -- Name
  CHARACTER( * ),  PARAMETER :: PROGRAM_NAME = 'RTM_TEST_FWD'

  ! -- Angle definitions
  INTEGER,         PARAMETER :: N_ANGLES  = 1
  REAL( fp_kind ), PARAMETER :: MIN_ANGLE =  0.0_fp_kind
  REAL( fp_kind ), PARAMETER :: MAX_ANGLE = 0.0_fp_kind

  ! -- Profile definitions
  CHARACTER( * ),  PARAMETER :: PROFILE_FILE        = 'profiles_LAYER.bin'
  INTEGER,         PARAMETER :: N_PROFILES          = 1
  INTEGER N_LAYERS    ! arbitrary level number to be kk


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
  REAL( fp_kind ), DIMENSION( N_ANGLES ) :: surface_temperature

  ! -- Profile read array
  REAL( fp_kind ), DIMENSION( : ), ALLOCATABLE  :: level_pressure
  REAL( fp_kind ), DIMENSION( : ), ALLOCATABLE  :: layer_pressure
  REAL( fp_kind ), DIMENSION( : ), ALLOCATABLE  :: layer_temperature
  REAL( fp_kind ), DIMENSION( : ), ALLOCATABLE  :: layer_water_vapor
  REAL( fp_kind ), DIMENSION( : ), ALLOCATABLE  :: layer_ozone



  ! -- Profile data arrays
  REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: level_p
  REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: layer_p
  REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: layer_t
  REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: layer_w
  REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: layer_o
 
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
  integer ::  first_time = 1


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

  N_LAYERS = kk  ! profile dependent length


  if (first_time == 1) then ! first time it is called
     first_time = 0

     call get_directory ('static',fname,len)
     fname = fname(1:len)//'optranlib/'
     len = len_trim(fname)
     
     error_status = initialize_rtm( tau_file = sndr_trans(1:sndr_trans_len),&
                   path = fname(1:len),  &
               spectral_file = sndr_coeff(1:sndr_coeff_len)  )

     IF ( error_status /= SUCCESS ) THEN
        CALL display_message( PROGRAM_NAME, &
             'Error initialzing RTM', &
             error_status )
        STOP
     END IF

  endif ! first_time called


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

! ALLOCATE ARRAYS FOR VARIABLE PROFILE CONTENTS
!

  ALLOCATE (      level_p                  ( N_LAYERS,N_ANGLES), &
                  layer_p                  ( N_LAYERS,N_ANGLES), &    
                  layer_t                  ( N_LAYERS,N_ANGLES), & 
                  layer_w                  ( N_LAYERS,N_ANGLES), & 
                  layer_o                  ( N_LAYERS,N_ANGLES), &
                  level_pressure           ( N_LAYERS), &
                  layer_pressure           ( N_LAYERS), & 
                  layer_temperature        ( N_LAYERS), & 
                  layer_water_vapor        ( N_LAYERS), & 
                  layer_ozone              ( N_LAYERS)   )

  CALL get_max_n_channels( n_channels )

  begin_channel = 1
  end_channel   = n_channels


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





    !#--------------------------------------------------------------------------#
    !#                           -- Assign variable from LAPS input             #
    !#--------------------------------------------------------------------------#

      do k = 1, kk
         level_p(k, 1 ) = lev_p(k+1)
      enddo
      layer_p(1:kk, 1 ) = lay_p(1:kk)
      layer_t(1:kk, 1 ) = lay_t(1:kk)
      layer_w(1:kk, 1 ) = lay_w(1:kk)
      layer_o(1:kk, 1 ) = lay_o(1:kk)
      surface_temperature = laps_sfc_t
      surface_emissivity = sfc_emis    !array assignment to constant
      surface_reflectivity = sfc_refl  !array assignment to constant
      secant_view_angle = sec_za
      secant_solar_angle = sec_solar
      n_channels_per_profile  = n_channels  ! don't understand

      do k = 1,n_channels
         channel_index(k) = k  ! don't understand
      enddo

      ! note the (in,out) cautions on upwelling radiance and brightness_temperature
      ! this is preset here.

      upwelling_radiance = 0.
      brightness_temperature = 1.0  ! recommended preset


   

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

!                                layer_t(N_Layers,:),    &  ! Input, M
                                surface_temperature,    &  ! Input, M
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

! here the returned arrays (upwelling_radiance and brightness_temperature)
! have dimensions of the insturment channels.  The LAPS array however is a 
! fixed one of 18.  this takes and only uses the part of the LAPS array needed

    
    channels_used = size (upwelling_radiance)

    do i = 1,channels_used
       up_radiance(i) = upwelling_radiance(i)
       bright_temp(i) = brightness_temperature(i)
    enddo




  !#----------------------------------------------------------------------------#
  !#                           -- DESTROY THE RTM --                            #
  !#----------------------------------------------------------------------------#

  ! ---------------------------------------
  ! Deallocate the channel dependent arrays
  ! ---------------------------------------

  ! -- Forward model

  DEALLOCATE (  layer_p,  layer_t  ,layer_w  ,layer_o ,  layer_pressure,  layer_temperature , layer_water_vapor ,  layer_ozone,  &
                level_p, level_pressure  )

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

 ! this utility is now handled by the new module below.  called from 
 ! variational.f after all calls to the forward model are complete.

END SUBROUTINE  optran90_fm

SUBROUTINE optran_deallocate (istatus)

!  USE type_kinds
!  USE file_utility
!  USE error_handler
!  USE parameters

  USE initialize
  integer :: istatus
  
  istatus = destroy_rtm()
  
!  IF ( error_status /= SUCCESS ) THEN
!     CALL display_message( PROGRAM_NAME, &
!          'Error destroying RTM', &
!          error_status )
!     STOP
!  END IF
  

END SUBROUTINE optran_deallocate


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
! Revision 1.3  2003/07/10 17:09:24  birk
! ready for laps optran90
!
! Revision 1.2  2002/11/18 20:01:39  birk
! changes made to avoid compilation errors on jet, statement order specifics.
!
! Revision 1.1  2002/11/15 15:21:32  birk
! Added to cvs mainly to see how this compiles on other platforms, it currently
! seems to compile on the IBM
!
! Revision 1.1  2001/09/13 22:08:13  paulv
! Initial checkin.
!
!
!
!
