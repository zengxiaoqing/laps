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

  SUBROUTINE integrated_liquid (nx,ny,nz,cond_mr,vapor_mr, rho,height,topo, &
                                intliqwater,totpcpwater)

    ! Computes integrated liquid water and total precip. water in a column.  
    !
    !  Adapted from USAF Weather Agency MM5 Post Processor
    !
    !  Brent Shaw, NOAA Forecast Systems Lab, Dec 2000

    IMPLICIT NONE
    ! Arguments
   
    INTEGER, INTENT(IN)                 :: nx
    INTEGER, INTENT(IN)                 :: ny
    INTEGER, INTENT(IN)                 :: nz
    REAL,    INTENT(IN)                 :: cond_mr(nx,ny,nz)
    REAL,    INTENT(IN)                 :: vapor_mr(nx,ny,nz)
    REAL,    INTENT(IN)                 :: rho(nx,ny,nz)
    REAL,    INTENT(IN)                 :: height(nx,ny,nz)
    REAL,    INTENT(IN)                 :: topo(nx,ny)
    REAL,    INTENT(OUT)                :: intliqwater(nx,ny)
    REAL,    INTENT(OUT)                :: totpcpwater(nx,ny)

    ! Locals
    REAL                                :: height_top,height_bot
    REAL                                :: dz
    INTEGER                             :: i,j,k

    DO j = 1, ny
      DO i = 1, nx
        intliqwater(i,j) = 0.0
        totpcpwater(i,j) = 0.0
        DO k = 1, nz
          ! Compute layer thickness
          IF (k.eq.1) THEN
            height_bot = topo(i,j)
            height_top = 0.5*(height(i,j,1)+height(i,j,2))
          ELSE IF (k.eq.nz) THEN
            height_bot = 0.5*(height(i,j,nz-1)+height(i,j,nz))
            height_top = 2*height(i,j,nz)-height_bot
          ELSE
            height_bot = 0.5*(height(i,j,k-1)+height(i,j,k))
            height_top = 0.5*(height(i,j,k)+height(i,j,k+1))
          ENDIF
          dz = height_top - height_bot
          intliqwater(i,j) = intliqwater(i,j)+cond_mr(i,j,k)*rho(i,j,k)*dz
          totpcpwater(i,j) = totpcpwater(i,j)+vapor_mr(i,j,k)*rho(i,j,k)*dz
        ENDDO
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE integrated_liquid
         
