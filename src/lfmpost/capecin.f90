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


  SUBROUTINE CAPECIN (PSIG, TSIG, THETAESIG, THETASIG, &
                      RHSIG, ZSIGFULL, &
                      TPRS, LI, & 
                      POSBUOYEN, NEGBUOYEN, K500, IMAX, JMAX, KSIG, KPRS)

!-----------------------------------------------------------------------
!--
!--   PURPOSE
!--   =======
!--   CALCULATE POSITIVE BUOYANT ENERGY (OR CONVECTIVE AVAILABLE
!--   ENERGY) AND NEGATIVE BUOYANT ENERGY (OR CONVECTIVE INHIBITION).
!--
!--   ALSO CALCULATE LIFTED INDEX.
!--
!--   REFERENCE
!--   =========
!--   DOSWELL AND RASMUSSEN (1994), WEA AND FCSTING, P 625.
!--
!--   UPDATES
!--   =======
!--   4 Jan 01  - Adapted from USAF Weather Agency routine
!--               B. Shaw, NOAA/FSL
!-----------------------------------------------------------------------

      USE constants
      IMPLICIT NONE

      INTEGER                     :: I
      INTEGER,   INTENT(IN)       :: IMAX
      INTEGER                     :: J
      INTEGER,   INTENT(IN)       :: JMAX
      INTEGER                     :: K
      INTEGER,   INTENT(IN)       :: K500
      INTEGER                     :: KMAX
      INTEGER,   INTENT(IN)       :: KSIG
      INTEGER,   INTENT(IN)       :: KPRS
      REAL,      INTENT(OUT)      :: LI        ( imax , jmax )
      REAL,      INTENT(OUT)      :: NEGBUOYEN ( imax , jmax )    
      REAL,      INTENT(OUT)      :: POSBUOYEN ( imax , jmax )    
      REAL                        :: PRSLCL
      REAL,      INTENT(IN)       :: PSIG      ( imax , jmax , ksig )
      REAL,      INTENT(IN)       :: RHSIG     ( imax , jmax , ksig )
      REAL,      INTENT(IN)       :: THETASIG  ( imax , jmax , ksig )
      REAL,      INTENT(IN)       :: THETAESIG ( imax , jmax , ksig )
      REAL                        :: THW
      REAL                        :: THWMAX
      REAL,      EXTERNAL         :: TLCL
      REAL                        :: TMPLCL
      REAL                        :: TPARCEL
      REAL,      INTENT(IN)       :: TPRS      ( imax , jmax , kprs )
      REAL,      INTENT(IN)       :: TSIG      ( imax , jmax , ksig )
      REAL,      EXTERNAL         :: WOBF
      REAL,      INTENT(IN)       :: ZSIGFULL  ( imax , jmax , ksig + 1 )
      REAL                        :: deltaz, dtheta, thetaparcel
      REAL, EXTERNAL              :: potential_temp
      REAL, PARAMETER             :: PREF = 1000.
      REAL                        :: psave
      LOGICAL                     :: compute_cin
      REAL, ALLOCATABLE           :: buoy(:)
      POSBUOYEN = 0.0
      NEGBUOYEN = 0.0
      ALLOCATE(buoy(ksig)) 
      DO J = 1, JMAX
        DO I = 1, IMAX

          THWMAX = -9999.0

          find_most_unstable: DO K = 1, KSIG

!           ------------------------------------------------------------
!           PICK THE MOST UNSTABLE PARCEL IN THE LOWEST 50 MB AS
!           INDICATED BY THE SIGMA LEVEL WITH THE HIGHEST WET BULB
!           POTENTIAL TEMPERATURE.  STORE INDEX IN KMAX.
!           ------------------------------------------------------------

            IF ( ((PSIG(I,J,1) - PSIG(I,J,K)) .LT. 50.0) ) THEN

              THW = THETAESIG(I,J,K) - WOBF(THETAESIG(I,J,K) - T0)

              IF (THW .GT. THWMAX) THEN
                KMAX = K
                THWMAX = THW
              END IF

            ELSE

              EXIT find_most_unstable

            END IF

          ENDDO find_most_unstable

!         --------------------------------------------------------------
!         CALCULATE LIFTED INDEX BY LIFTING THE MOST UNSTABLE
!         PARCEL.
!         --------------------------------------------------------------

          CALL THE2T (THETAESIG(I,J,KMAX), 500.0, TPARCEL)
          LI(I,J)  = TPRS(I,J,K500) - TPARCEL
    
!         --------------------------------------------------------------
!         CALCULATE THE TEMPERATURE AND PRESSURE OF THE LIFTING
!         CONDENSATION LEVEL.
!         --------------------------------------------------------------

          TMPLCL = TLCL ( TSIG(I,J,KMAX), RHSIG(I,J,KMAX) )
          PRSLCL = PSIG(I,J,KMAX) * (TMPLCL / TSIG(I,J,KMAX)) ** CPOR

!         --------------------------------------------------------------
!         CALCULATE THE BUOYANCY.
!         --------------------------------------------------------------
          posbuoyen(i,j) = 0.
          negbuoyen(i,j) = 0.
          DO K = KMAX, KSIG

!           ------------------------------------------------------------
!           ABOVE THE LCL, CALCULATE VIRTUAL TEMPERATURE OF THE
!           PARCEL AS IT MOVES ALONG A MOIST ADIABAT.  BELOW THE
!           LCL, LIFT PARCEL ALONG A DRY ADIABAT.
!           ------------------------------------------------------------
            
            IF (PSIG(I,J,K) .LE. PRSLCL) THEN
              
              CALL THE2T (THETAESIG(I,J,KMAX), PSIG(I,J,K), TPARCEL)
            ELSE
   
              TPARCEL   = THETASIG(I,J,KMAX) /(PREF / PSIG(I,J,K)) ** KAPPA

            END IF

            
            ! Compute the potential temperature of the parcel
            thetaparcel = potential_temp(tparcel, psig(i,j,k)*100.)
            dtheta = thetaparcel - thetasig(i,j,k)
            deltaz = zsigfull(i,j,k+1) - zsigfull(i,j,k)
            buoy(k) = deltaz * dtheta/thetasig(i,j,k)
          ENDDO
          
          ! Now loop through the column again, partitioning the buoyency
          ! into positive (CAPE) and negative (CIN) component.  We terminate
          ! the contribution to CIN when/if a layer of CAPE greater than
          ! 150 mb deep is found.
 
          compute_cin = .true.
          psave = -100.
          DO k = kmax, ksig
            IF (buoy(k) .GT. 0.) THEN
              IF (psave .LT. 0) THEN
                psave = psig(i,j,k)
              ELSE
                IF ( (psave-psig(i,j,k)) .GT. 150.) THEN
                  compute_cin = .false.
                ENDIF
              ENDIF
              posbuoyen(i,j) = posbuoyen(i,j) + buoy(k)
            ELSE IF (buoy(k).LT.0.) THEN
              psave = -100.
              IF (compute_cin) THEN
                negbuoyen(i,j) = negbuoyen(i,j) + buoy(k)
              ENDIF
            ENDIF 
          ENDDO
        ENDDO
      ENDDO
    
      posbuoyen = grav * posbuoyen
      negbuoyen = grav * negbuoyen
 
      ! Cap the negative buoyancy to a maximum value of 700 J/kg

      WHERE(negbuoyen .LT. -700) negbuoyen = -700. 
      DEALLOCATE(buoy)
      END SUBROUTINE CAPECIN
