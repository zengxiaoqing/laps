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
  SUBROUTINE HELICITY (USIG, VSIG, ZSIG, TERRAIN, IMAX, JMAX, KSIG, SRELHEL)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--
!--   NAME: HELICITY
!--
!--   PURPOSE
!--   =======
!--   CALCULATE STORM RELATIVE HELICITY.
!--
!--   REMARKS
!--   =======
!--   HELICITY IS EQUAL TO THE VERTICAL INTEGRAL OF:
!-- 
!--   (V - c) x dV/dz
!--
!--   WHERE V IS THE MODEL PREDICTED WIND AND C IS THE 
!--   ESTIMATED STORM MOTION VECTOR. 
!--
!--   UPDATES
!--   =======
!--   ??? 97   INITIAL VERSION......................................DNXM
!--   DEC 98   MODIFIED CALL TO WDIR FUNCTION.  AN ERROR IN WDIR WAS 
!--            CORRECTED............................................DNXM
!--   JAN 99   CHANGED VARIABLE NAMES USIGCRS AND VSIGCRS TO USIG AND
!--            AND VSIG AS WINDS ARE NOW AT THE CROSS POINT.........DNXM
!--   JAN 01   Adapted for use at FSL by B. Shaw, NOAA/FSL
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      USE constants
      IMPLICIT NONE

      REAL                         :: DZ
      INTEGER                      :: I
      INTEGER,     INTENT(IN)      :: IMAX
      INTEGER                      :: J
      INTEGER,     INTENT(IN)      :: JMAX
      INTEGER                      :: K
      INTEGER                      :: K3KM
      INTEGER                      :: K10KM
      INTEGER,     INTENT(IN)      :: KSIG
      REAL,        INTENT(OUT)     :: SRELHEL  ( imax , jmax )    
      REAL                         :: STORMDIR
      REAL                         :: STORMSPD
      REAL                         :: STORMU
      REAL                         :: STORMV
      REAL                         :: SUMDZ
      REAL                         :: SUMHEL
      REAL                         :: SUMU
      REAL                         :: SUMV
      REAL,        INTENT(IN)      :: TERRAIN  ( imax , jmax )
      REAL,        INTENT(IN)      :: USIG     ( imax , jmax , ksig )
      REAL,        INTENT(IN)      :: VSIG     ( imax , jmax , ksig )
      REAL,        EXTERNAL        :: WDIR
      REAL,        EXTERNAL        :: WSPD
      REAL,        INTENT(IN)      :: ZSIG     ( imax , jmax , ksig )

!-----------------------------------------------------------------------
!--   DETERMINE INDICES OF 3 AND 10 KM LAYERS (ABOVE GROUND LEVEL)
!-----------------------------------------------------------------------

      DO J = 1, JMAX
        DO I = 1, IMAX

          DO K = 1, KSIG

            IF (((ZSIG(I,J,K) - TERRAIN(I,J)) .GE. 3000.0)) THEN
              K3KM = K
              EXIT
            END IF
     
          ENDDO


          DO K = KSIG, K3KM, -1

            IF (((ZSIG(I,J,K) - TERRAIN(I,J)) .LE. 10000.0)) THEN
              K10KM = K
              EXIT 
            END IF

          ENDDO

!-----------------------------------------------------------------------
!--       ESTIMATE STORM MOTION VECTOR. SPEED IS 75 PERCENT OF THE MEAN
!--       WIND BETWEEN 3 AND 10 KM.  DIRECTION IS 30 DEGREES TO THE
!--       RIGHT OF THE MEAN WIND.
!-----------------------------------------------------------------------

          SUMU  = 0.0
          SUMV  = 0.0
          SUMDZ = 0.0

          DO K = K3KM, K10KM

            DZ    = ZSIG(I,J,K) - ZSIG(I,J,K-1)
            SUMDZ = SUMDZ + DZ
            SUMU  = SUMU + (0.5 * DZ * (USIG(I,J,K) + USIG(I,J,K-1)))
            SUMV  = SUMV + (0.5 * DZ * (VSIG(I,J,K) + VSIG(I,J,K-1)))

          ENDDO

          STORMU = SUMU / SUMDZ
          STORMV = SUMV / SUMDZ

          STORMSPD = WSPD (STORMU, STORMV) * 0.75

!-----------------------------------------------------------------------
!--       WHEN CALLING WDIR, SEND IN A CONE FACTOR OF ZERO
!--       SO THAT DIRECTION IS GRID RELATIVE.
!-----------------------------------------------------------------------

          STORMDIR = WDIR (STORMU, STORMV, 0.0, 0.0, 0.0) + 30.0
          IF (STORMDIR .GT. 360.0) STORMDIR = STORMDIR - 360.0

          STORMU = -STORMSPD * SIN(STORMDIR * DEG2RAD)
          STORMV = -STORMSPD * COS(STORMDIR * DEG2RAD)

          SUMHEL = 0.0
  
!-----------------------------------------------------------------------
!--       CALCULATE HELICITY.  INTEGRATE BETWEEN THE GROUND AND 3 KM,
!--       A DEPTH THAT IS FREQUENTLY USED IN THE LITERATURE.
!-----------------------------------------------------------------------

          DO K = 1, K3KM
     
            SUMHEL = SUMHEL +  &
                   ((USIG(I,J,K) - STORMU) *  &
                    (VSIG(I,J,K+1) - VSIG(I,J,K))) -  &
                   ((VSIG(I,J,K) - STORMV) *     &
                    (USIG(I,J,K+1) - USIG(I,J,K)))

      
          ENDDO
          
          SRELHEL(I,J) = -SUMHEL
  
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE HELICITY
