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

SUBROUTINE model_pblhgt(theta,thsfc,psig,zsig,topo,nx,ny,nz,pblhgt)

  !  Subroutine to estimate height AGL in meters of PBL from native
  !  coordinate model data.  Adapted from the LAPS routine for 
  !  terrain-following model coordinates.

  IMPLICIT NONE

  INTEGER, INTENT(IN)   :: nx,ny,nz
  REAL, INTENT(IN)      :: theta(nx,ny,nz)
  REAL, INTENT(IN)      :: thsfc(nx,ny)
  REAL, INTENT(IN)      :: psig(nx,ny,nz)
  REAL, INTENT(IN)      :: zsig(nx,ny,nz)
  REAL, INTENT(IN)      :: topo(nx,ny)
  REAL, INTENT(OUT)     :: pblhgt(nx,ny)

  INTEGER  :: i,j,k,ktop
  REAL     :: thresh_k, topwgt, botwgt
  LOGICAL  :: found_pbl_top
  
  PRINT *, 'Generating PBL height using theta and surface temp.'
  loop_j:  DO j = 1 , ny
    loop_i:  DO i = 1 , nx

      ! Compute threshold value that theta needs to exceed
      ! to be above PBL.  We use surface theta plus an
      ! additional 3K for slop.

      thresh_k = thsfc(i,j) + 3.0  

      ! Now begin at the bottom and work our way up until
      ! we find the first level with a theta exceeding the
      ! threshold

      found_pbl_top = .false.
      loop_k: DO k = 1 , nz
 
        IF (theta(i,j,k) .GE. thresh_k) THEN
          ktop = k
          found_pbl_top = .true.
          EXIT loop_k
        ENDIF
      ENDDO loop_k

      ! If we did not find a good PBL, set PBL to first level
      ! and print out some diagnostics
      IF (.NOT. found_pbl_top) THEN
        PRINT *, 'PBL height not found at i/j = ',i,j
        PRINT *, 'Surface Theta = ', thsfc(i,j)
        PRINT *, 'Theta in the column:'
        PRINT *, 'Pressure Height  Theta'
        PRINT *, '-------- ------- --------'
        diag_loop: DO k = 1,nz
          PRINT '(F8.0,F8.0,F8.2)',psig(i,j,k),zsig(i,j,k),theta(i,j,k)
        ENDDO diag_loop
        ktop = 1
        pblhgt(i,j) = zsig(i,j,1) - topo(i,j)

      ELSE 
  
        ! We found the top k-level bounding the PBL so interpolate
        ! to the actual level
 
        IF (ktop .EQ. 1) THEN
          pblhgt(i,j) = zsig(i,j,1) - topo(i,j)
        ELSE
          ! Interpolate to get height at thresh_k
          botwgt = ( (theta(i,j,ktop)-thresh_k) / &
                     (theta(i,j,ktop)-theta(i,j,ktop-1)) )
          topwgt = 1.0 - botwgt
          pblhgt(i,j) = botwgt * zsig(i,j,ktop-1) + &
                        topwgt * zsig(i,j,ktop) - topo(i,j)
        ENDIF
      ENDIF
    ENDDO loop_i
  ENDDO loop_j
  RETURN
END SUBROUTINE model_pblhgt
