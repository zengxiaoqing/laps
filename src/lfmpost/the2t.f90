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

  SUBROUTINE the2t(thetae,p,tparcel)
 
    ! Calculates temperature at any pressure level along a saturation
    ! adiabat by iteratively solving eq 2.76 from Wallace and Hobbs (1977).

    ! Adpated from USAF routine, B. Shaw, NOAA/FSL

    USE constants
    IMPLICIT NONE

    INTEGER                          :: iter
    REAL, EXTERNAL                   :: mixsat
    REAL, INTENT(IN)                 :: p
    REAL                             :: tcheck
    REAL                             :: theta
    REAL, INTENT(IN)                 :: thetae
    REAL                             :: tovtheta
    REAL, INTENT(OUT)                :: tparcel
    LOGICAL                          :: converged

    converged = .false.
    tovtheta = (p / 1000.0) ** kappa
    tparcel = thetae / exp(lv * 0.012 / (cp *295.0)) *tovtheta

    iterloop: DO iter = 1 , 50 

      theta = thetae / exp(lv*mixsat(tparcel,p*100.)/(cp*tparcel))
      tcheck = theta * tovtheta

      IF (ABS(tparcel - tcheck) .lt. 0.05) THEN
        converged = .true.
        EXIT iterloop
      ENDIF
      tparcel = tparcel + (tcheck - tparcel) * 0.3
 
    ENDDO iterloop
    IF (.NOT.converged) THEN
      PRINT *, 'WARNING!! THETAE to TEMP calc did not converge!'
      print *, 'thetae and p:',thetae,p
    ENDIF

    RETURN
  END SUBROUTINE the2t
