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

  SUBROUTINE reflectivity(nx,ny,nz,rho,hgt,rainmr,icemr,snowmr,graupelmr, &
                          refl,max_refl,echo_tops)

    ! Subroutine to compute estimated radar reflectivity from
    ! the precipitation mixing ratios.  Will also return 
    ! column max reflectivity and echo tops.  The estimation
    ! is done using formulae from Kessler (1969) and 
    ! Rogers and Yau (1989).  
    !
    ! Adapted from USAF Weather Agency routine.  
    !
    !  Brent Shaw, NOAA Forecast System Lab, Dec 2000
    !
  
    IMPLICIT NONE

    INTEGER, INTENT(IN)             :: nx
    INTEGER, INTENT(IN)             :: ny
    INTEGER, INTENT(IN)             :: nz
    REAL, INTENT(IN)                :: rho(nx,ny,nz)
    REAL, INTENT(IN)                :: hgt(nx,ny,nz)
    REAL, INTENT(IN)                :: rainmr(nx,ny,nz)
    REAL, INTENT(IN)                :: icemr(nx,ny,nz)
    REAL, INTENT(IN)                :: snowmr(nx,ny,nz)
    REAL, INTENT(IN)                :: graupelmr(nx,ny,nz)
    REAL, INTENT(OUT)               :: refl(nx,ny,nz)
    REAL, INTENT(OUT)               :: max_refl(nx,ny)
    REAL, INTENT(OUT)               :: echo_tops(nx,ny)
   
    ! Locals
    REAL, PARAMETER                 :: svnfrth = 7.0/4.0
    REAL, PARAMETER                 :: max_top_thresh =  5.0
    INTEGER                         :: i,j,k
    REAL                            :: w
    max_refl = 0.0
    echo_tops = 1.0e37
    refl = 0.0

    DO j = 1,ny
      DO i = 1, nx
        DO k=1,nz

          ! Compute the basic reflectivity
          !refl(i,j,k) =17300.0 * &
          !            (rho(i,j,k) * 1000.0 * &
          !             MAX(0.0,rainmr(i,j,k)))**svnfrth

          ! Add the ice component
          !refl(i,j,k)=refl(i,j,k) + &
          !       38000.0*(rho(i,j,k) * 1000.0 * &
          !       MAX(0.0,icemr(i,j,k)+snowmr(i,j,k)+graupelmr(i,j,k)))**2.2

          ! Convert to dBZ
          !refl(i,j,k) = 10.*ALOG10(MAX(refl(i,j,k),1.0))

          ! Test RAMS reflectivity algorithm
          w = 264083.11 * (rainmr(i,j,k) + &
              0.2*(icemr(i,j,k) + snowmr(i,j,k)) + &
              2.0 * graupelmr(i,j,k))
          w = MAX(1.,w)
          refl(i,j,k) = 17.8 * alog10(w)
          ! Since we are going from the ground up, we can 
          ! check threshold and set echo top.

          IF (refl(i,j,k) .GE. max_top_thresh) echo_tops(i,j)=hgt(i,j,k) 
        
        ENDDO
        ! Compute the max value in the column
        max_refl(i,j) = MAXVAL(refl(i,j,:))
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE reflectivity
          

      
 
