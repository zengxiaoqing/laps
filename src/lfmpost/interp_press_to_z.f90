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

  SUBROUTINE interp_press_to_z(press_levels, heights, new_level, press_out, &
                                  nx,ny,nz)

    ! Given a 1D array of pressure levels and a 3D array of heights (m) at
    ! those levels, this routine interpolates the pressure to the desired
    ! new_level.

    ! Pressures are in Pa, heights are in m!

    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: nx,ny,nz
    REAL, INTENT(IN)                  :: press_levels(nz)
    REAL, INTENT(IN)                  :: heights(nx,ny,nz)
    REAL, INTENT(IN)                  :: new_level
    REAL, INTENT(OUT)                 :: press_out (nx,ny)
    INTEGER                           :: i,j,k, ktop
    REAL                              :: ptop,pbot,ztop,zbot

    nsloop: DO j = 1, ny
      ewloop: DO i = 1, nx
        vertloop: DO k = 2, nz
          IF (heights(nx,ny,nz) .ge. new_level) THEN
             ktop = k
             EXIT vertloop
          ENDIF
        ENDDO vertloop
        ztop = heights(i,j,k)
        zbot = heights(i,j,k-1)
        ptop = press_levels(k)*0.01    ! Convert to mb!
        pbot = press_levels(k-1)*0.01
        press_out(i,j) = EXP ( &
                              ( new_level*ALOG(pbot/ptop) - &
                              ztop*ALOG(pbot) + &
                              zbot*ALOG(ptop)  ) / &
                              (zbot - ztop) ) * 100.  ! convert back to Pa
      ENDDO ewloop
    ENDDO nsloop
    RETURN
  END SUBROUTINE interp_press_to_z
 
    

          
