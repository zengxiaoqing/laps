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

MODULE fire 

  ! Module to contain subroutines for computing various fire weather
  ! indices.

  IMPLICIT NONE

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ventilation(usig,vsig,zsig,pblhgt,topo,nx,ny,nz, &
                         upbl, vpbl, vent_ind)
  
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: nx
    INTEGER, INTENT(IN)    :: ny
    INTEGER, INTENT(IN)    :: nz
    REAL, INTENT(IN)       :: usig(nx,ny,nz)  ! U on sigma
    REAL, INTENT(IN)       :: vsig(nx,ny,nz)  ! V on sigma
    REAL, INTENT(IN)       :: zsig(nx,ny,nz)  ! Z on sigma
    REAL, INTENT(IN)       :: pblhgt(nx,ny)
    REAL, INTENT(IN)       :: topo(nx,ny)
    REAL, INTENT(OUT)      :: upbl(nx,ny)
    REAL, INTENT(OUT)      :: vpbl(nx,ny)
    REAL, INTENT(OUT)      :: vent_ind(nx,ny)

    INTEGER                :: i,j,k, nbl
    REAL                   :: usum, vsum, umean, vmean, spmean

    PRINT *, '---- Subroutine ventilation ----'
  
    DO j = 1 , ny
      DO i = 1 , nx
        
        IF (pblhgt(i,j) .GT. 0.) THEN

          ! Compute mean wind within the PBL
          nbl = 0
          usum = 0.
          vsum = 0.
          umean = 0.
          vmean = 0.
          sum_pbl: DO k = 1 , nz
            IF (zsig(i,j,k)-topo(i,j) .LE. pblhgt(i,j)) THEN
              nbl = nbl + 1
              usum = usum + usig(i,j,k)
              vsum = vsum + vsig(i,j,k)
            ELSE
              EXIT sum_pbl
            ENDIF
          ENDDO sum_pbl 
          IF (nbl .GT. 0) THEN
            umean = usum/FLOAT(nbl)
            vmean = vsum/FLOAT(nbl)
   
            ! Compute mean wind speed for this layer
            spmean = SQRT(umean**2 + vmean**2)
            
            ! Multiply mean speed by PBL depth to get index
            vent_ind(i,j) = pblhgt(i,j) * spmean
            upbl(i,j) = umean
            vpbl(i,j) = vmean
          ELSE
            ! PBL height is lower than the lowest model level...use
            ! lowest model wind
            spmean = SQRT(usig(i,j,1)**2 + vsig(i,j,1)**2)
            vent_ind(i,j) = pblhgt(i,j) * spmean
            upbl(i,j) = usig(i,j,1)
            vpbl(i,j) = vsig(i,j,1)
          ENDIF
        ELSE
          PRINT *, 'WARNING:  PBL Height <=0 in ventilation index'
          vent_ind(i,j) = 0.
          upbl(i,j) = 0.
          vpbl(i,j) = 0.
        ENDIF
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE ventilation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE haines_layer(p3d_mb, t3d_k, td3d_k, haines2d, nx, ny, nz, &
                    pmbbot, pmbtop)

    ! Computes haines index for layer bounded by top and bottom pressure 
    ! levels pmbtop and pmbbot (mb).

    IMPLICIT NONE

    INTEGER, INTENT(IN)      :: nx, ny, nz
    REAL, INTENT(IN)         :: p3d_mb(nx,ny,nz)  ! 3D pressure in mb
    REAL, INTENT(IN)         :: t3d_k (nx,ny,nz)  ! 3D Temp in K
    REAL, INTENT(IN)         :: td3d_k(nx,ny,nz)  ! 3D Dewpoint in K
    REAL, INTENT(OUT)        :: haines2d(nx,ny)   ! 2D haines index
    REAL, INTENT(IN)         :: pmbbot,pmbtop     ! Bounding mb levels

    ! Local Variables

    INTEGER :: i,j,k, kk, km1
    REAL :: tmkt, tmkb, tdkb, deltat, dpdep
    REAL :: factor1, factor2
    PRINT *, '---- Subroutine haines_layer ----'
    DO j = 1 , ny
      DO i = 1 , nx
      
        IF (p3d_mb(i,j,1) .lt. pmbbot) THEN
          haines2d(i,j) = 1e37  ! Cannot be computed
        ELSE
    
          DO k = 2, nz
        
            ! Account for flipped vertical coordinate from original 
            ! algorithm from Seattle
  
            kk = nz+1-k
            km1 = kk + 1 
            
            ! Find temperature at the top
            IF (p3d_mb(i,j,kk).GT.pmbtop.AND.p3d_mb(i,j,km1).LE.pmbtop) THEN 
              tmkt = t3d_k(i,j,km1) + (t3d_k(i,j,kk)-t3d_k(i,j,km1)) * &
                     (log(pmbtop)-log(p3d_mb(i,j,km1))) / &
                     (log(p3d_mb(i,j,kk))-log(p3d_mb(i,j,km1)))
            ENDIF

            ! Find Temp/dewpoint at the bottom of the layer
            
            IF (p3d_mb(i,j,kk).GT.pmbbot.AND.p3d_mb(i,j,km1).LE.pmbbot) THEN
              tmkb = t3d_k(i,j,km1) + (t3d_k(i,j,kk)-t3d_k(i,j,km1)) * &
                     (log(pmbbot) - log(p3d_mb(i,j,km1))) / &
                     (log(p3d_mb(i,j,kk))-log(p3d_mb(i,j,km1)))
              tdkb = td3d_k(i,j,km1) + (td3d_k(i,j,kk) - td3d_k(i,j,km1)) * &
                     (log(pmbbot) - log(p3d_mb(i,j,km1))) / &
                     (log(p3d_mb(i,j,kk))-log(p3d_mb(i,j,km1)))
            ENDIF

          ENDDO

          deltat = tmkb - tmkt
          dpdep  = tmkb - tdkb

          IF (NINT(pmbbot) .EQ. 700) THEN ! High Haines
            IF (deltat.LE. 17.5) THEN
              factor1 = 1.
            ELSEIF(deltat .GT. 17.5 .AND. deltat .LE. 21.5 ) THEN
              factor1 = 2.   ! deltat > 21.5
            ELSE
              factor1 = 3.
            ENDIF

            IF (dpdep.LE. 14.5) THEN
              factor2 = 1.
            ELSEIF(dpdep .GT. 14.5 .AND. dpdep .LE. 20.5) THEN
              factor2 = 2.
            ELSE    ! dpdep > 20.5
              factor2 = 3.
            ENDIF

          ELSEIF(NINT(pmbbot) .EQ. 850) THEN   ! Mid-level Haines
            IF (deltat .LE. 5.5) THEN
              factor1 = 1.
            ELSEIF(deltat .GT. 5.5 .AND. deltat .LE. 10.5) THEN
              factor1 = 2.
            ELSE   ! deltat > 10.5
              factor1 =3.
            ENDIF

            IF (dpdep .LE. 5.5 ) THEN
              factor2 = 1.
            ELSEIF( dpdep .GT. 5.5 .AND. dpdep .LE. 12.5) THEN
              factor2 = 2.
            ELSE   ! dpdep > 12.5
              factor2 = 3.
            ENDIF
   
          ELSE ! Cannot determine mid or high
  
            PRINT *, 'Haines_layer needs 850 or 700 as bottom layer'
            PRINT *, 'Bottom level (mb) specified: ', pmbbot
 
            STOP 'bad_haines_layer'

          ENDIF

          haines2d(i,j) = factor1 + factor2

        ENDIF 
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE haines_layer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE fosberg_fwi(t2k, rh2pct, p_sfc_mb, u10, v10, &
                         nx, ny, fwi)

    ! Computes the Fosberg Fire Weather Index
 
    IMPLICIT NONE
    INTEGER, INTENT(IN)              :: nx,ny
    REAL, INTENT(IN)                 :: t2k(nx,ny)   ! Sfc Temp (K)
    REAL, INTENT(IN)                 :: rh2pct(nx,ny)   ! Sfc RH (%)
    REAL, INTENT(IN)                 :: p_sfc_mb(nx,ny) ! Sfc Press (mb)
    REAL, INTENT(IN)                 :: u10(nx,ny)      ! 10 m U wind (m/s)
    REAL, INTENT(IN)                 :: v10(nx,ny)      ! 10 m V wind (m/s)
    REAL, INTENT(OUT)                :: fwi(nx,ny)      ! Fosberg Index

    ! Local Variables

    INTEGER :: i,j
    REAL :: uuu10, vvv10, m, n, t2f,rh2

    DO j = 1 , ny
      DO i = 1 , nx

        ! Convert Temperature from K to F
        t2f = 1.8 * (t2k(i,j)-273.15) + 32.0
        
        ! Convert u/v from m/s to mph
        uuu10 = u10(i,j) * 2.237
        vvv10 = v10(i,j) * 2.237
        rh2 = rh2pct(i,j)

        IF ( rh2 .LE. 10.5 ) THEN
           m = 0.03229 + (0.281073 * rh2) - (0.000578 * rh2 * t2f)
        ELSE IF( rh2 .GT. 10.5 .AND. rh2 .LE. 50.5 ) THEN
           m = 2.22749 + (0.160107 * rh2) - (0.014784 * t2f)
        ELSE IF( rh2 .GT. 50.5 .AND. rh2 .LE. 100 ) THEN
           m = 21.0606 + (0.005565 * rh2**2) - (0.00035 * rh2 * t2f) - &
              (0.483199 * rh2)
        ELSE         !   rh2 > 100
           m = 21.0606 + (0.005565 * 100**2) - (0.00035 * 100 * t2f) - &
               (0.483199 * 100)
           PRINT *, 'fwi calculation: rh2 > 100 ', j, i, rh2, t2f
        ENDIF

        n = 1.0 - 2.0*(m/30.) + 1.5*(m/30.)**2 - 0.5*(m/30.)**3
        fwi(i,j) = (n * SQRT(1.0 + uuu10**2 + vvv10**2)) / 0.3002

      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE fosberg_fwi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE fire
