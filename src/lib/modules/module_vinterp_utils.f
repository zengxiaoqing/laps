!dis   
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis    
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis   
!dis 

MODULE vinterp_utils

  ! Module that contains arrays and routines for linear and logarithmic
  ! interpolation 

  
  IMPLICIT NONE
  
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE compute_lin_weights(p_mid_mb,p_lower_mb,p_upper_mb, &
                                  weight_bot, weight_top)

    ! Computes the weighting coefficient for the top bounding level
    ! of two pressure levels to interpolate to a level in-between using
    ! linear interpolation.
   
    IMPLICIT NONE
    REAL, INTENT(IN)             :: p_mid_mb  ! Desired pressure level
    REAL, INTENT(IN)             :: p_lower_mb ! Lower bounding pressure
    REAL, INTENT(IN)             :: p_upper_mb ! Upper bounding pressure
    REAL, INTENT(OUT)            :: weight_bot ! Weight given to bottom level
    REAL, INTENT(OUT)            :: weight_top ! Weight given to top level

    weight_bot = (p_mid_mb - p_upper_mb) / (p_lower_mb - p_upper_mb) 
    weight_top = 1.0 - weight_bot
    
  END SUBROUTINE compute_lin_weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE compute_log_weights(p_mid_mb,p_lower_mb,p_upper_mb, &
                                 weight_bot, weight_top)
    
    ! Computes weighting coefficient for upper pressure level that bounds
    ! a desired pressure level for logarithmic interpolation.
    
    ! NOTE: Pressures must be in mb!!!
    
    IMPLICIT NONE
    REAL, INTENT(IN)             :: p_mid_mb  ! Desired pressure level
    REAL, INTENT(IN)             :: p_lower_mb ! Lower bounding pressure
    REAL, INTENT(IN)             :: p_upper_mb ! Upper bounding pressure
    REAL, INTENT(OUT)            :: weight_bot ! Weight given to bottom level
    REAL, INTENT(OUT)            :: weight_top ! Weight given to top level
    
    weight_bot = ( ALOG(p_mid_mb) - ALOG(p_upper_mb) ) / &
                 ( ALOG(p_lower_mb) - ALOG(p_upper_mb) ) 
    weight_top = 1.0 - weight_bot

  END SUBROUTINE compute_log_weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE vinterp_3d(var_in, trap_bot_ind, trap_top_ind, &
                        weight_top, below_ground_value, var_out, &
                        nx, ny, nzin, nzout)

    ! Interpolates var_in to var_out using array of trapping indices
    ! and top weight values.  Allows the user to set a default
    ! 2d array to use for "below ground" values

    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: nx
    INTEGER, INTENT(IN)            :: ny
    INTEGER, INTENT(IN)            :: nzin
    INTEGER, INTENT(IN)            :: nzout
    REAL, INTENT(IN)               :: var_in (nx,ny,nzin)
    INTEGER, INTENT(IN)            :: trap_bot_ind(nx,ny,nzout)
    INTEGER, INTENT(IN)            :: trap_top_ind(nx,ny,nzout)
    REAL, INTENT(IN)               :: weight_top(nx,ny,nzout)
    REAL, INTENT(IN)               :: below_ground_value(nx,ny)
    REAL, INTENT(OUT)              :: var_out(nx,ny,nzout)
    INTEGER                        :: i,j,k

    DO j = 1,ny 
      DO i = 1,nx
        DO k= 1,nzout

          ! Is level below ground?  If so, zero it out
          IF (trap_bot_ind(i,j,k).LT. 1) THEN
            var_out(i,j,k) = below_ground_value(i,j)
          
          ! Is it above model top? If so, replicate top value
          ELSE IF (trap_bot_ind(i,j,k).EQ.nzin) THEN
            var_out(i,j,k) = var_in(i,j,nzin)

          ELSE
            var_out(i,j,k) = weight_top(i,j,k) * &
                             var_in(i,j,trap_top_ind(i,j,k)) + &
                             (1.-weight_top(i,j,k)) * &
                             var_in(i,j,trap_bot_ind(i,j,k))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
 
    RETURN
  END SUBROUTINE vinterp_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE vinterp_utils  
