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

  SUBROUTINE clouds ( nx, ny, nz, grid_spacing, cldliqmr, cldicemr, snowmr, &
                      height, topo, cldbase, cldtop, ceiling, cldamt)

    ! This routine is used to compute cloud ceiling (AGL), cloud base
    ! height (MSL), cloud top height (MSL), and coverage (fraction) given
    ! mixing ratios of the various microphysical species.  

    ! Adapted from the USAF Weather Agency MM5 post processer.
    !  Brent Shaw, NOAA Forecast Systems Lab, Dec 2000

    IMPLICIT NONE
    
    ! Arguments
    INTEGER, INTENT(IN)             :: nx      ! Array x-dimension
    INTEGER, INTENT(IN)             :: ny      ! Array y-dimension
    INTEGER, INTENT(IN)             :: nz      ! Array z-dimension
    REAL, INTENT(IN)                :: grid_spacing  ! In meters
    REAL, INTENT(IN)                :: cldliqmr (nx,ny,nz) ! Cloud Liq Mix Rat
    REAL, INTENT(IN)                :: cldicemr (nx,ny,nz) ! Cloud Ice Mix Rat
    REAL, INTENT(IN)                :: snowmr   (nx,ny,nz) ! Snow Mix Rat
    REAL, INTENT(IN)                :: height   (nx,ny,nz) ! Heights in meters
    REAL, INTENT(IN)                :: topo     (nx,ny)    ! Topography in m
    REAL, INTENT(OUT)               :: cldbase  (nx,ny)
    REAL, INTENT(OUT)               :: cldtop   (nx,ny)
    REAL, INTENT(OUT)               :: ceiling  (nx,ny)
    REAL, INTENT(OUT)               :: cldamt   (nx,ny)

    ! Locals
    INTEGER                         :: i,j,k
    REAL                            :: icethresh
    REAL                            :: liqthresh
    REAL                            :: snowthresh

    ! Intialize the arrays
     
    cldbase = 0.0
    cldtop = 0.0
    cldamt = 0.0
    ceiling = 1e+30  ! unlimited ceiling

    ! Set thresholds for cloud determination based on grid
    ! resolution.  Kind of hokey, but will get us by for now

    IF (grid_spacing .LE. 10000) THEN
      icethresh  = 0.000005
      snowthresh    = 0.000003
      liqthresh  = 0.000003
    ELSE
      icethresh  = 0.000005
      snowthresh    = 0.000025
      liqthresh  = 0.000025
    ENDIF

    ! Loop through using these thresholds to determine cloudiness

    DO j = 1 , ny
      DO i = 1, nx
         find_base: DO k=1,nz
           IF ( (cldliqmr(i,j,k) .GE. liqthresh) .OR. &
                (cldicemr(i,j,k) .GE. icethresh) .OR. &
                (snowmr(i,j,k)   .GE. snowthresh) ) THEN
              cldbase(i,j) = height(i,j,k)
              ceiling(i,j) = height(i,j,k) - topo(i,j)
              
              ! For now all we can do is use cldamt as a yes/no
              ! We should look at coming up with a fraction function

              cldamt(i,j) = 1.0
              EXIT find_base
            ENDIF
          ENDDO find_base
          find_top: DO k = nz,1,-1
            IF ( (cldliqmr(i,j,k) .GE. liqthresh) .OR. &
                 (cldicemr(i,j,k) .GE. icethresh) .OR. &
                 (snowmr(i,j,k)   .GE. snowthresh) ) THEN
               cldtop(i,j) = height(i,j,k)
               EXIT find_top
            ENDIF
          ENDDO find_top
       ENDDO
    ENDDO    
    
    RETURN
  END SUBROUTINE clouds
    

