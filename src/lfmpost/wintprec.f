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
  SUBROUTINE WINTPREC (TSIG, ZSIG, ZPRS, PSFC, TSFC, &
                       TERRAIN, PRCPINC, IMAX, JMAX, &
                       KSIG, KPRS, K700, K850, &
                       K1000, PRCPCONV, PRECIPTYPE)         

!-----------------------------------------------------------------------  
!--   (WINTER) PRECIPITATION ALGORITHM
!--
!--   PURPOSE:  This product identifies areas of precipitation and the
!--   type expected, based on MM5 data.  The process is essentially two-
!--   fold.  First, areas of precipitation are identified when the MM5 
!--   precipitation array (PRCPINC) exceeds 0.01 inch. 
!--
!--   Second, thickness thresholds are used at the gridpoints for which 
!--   it has been determined that precipitation is occurring and the    
!--   surface pressure is greater than 850mb (i.e., non-mountainous    
!--   regions).  The thickness thresholds utilized are based on         
!--   meteorological research from the pertinent sources
!--   (e.g., MWR, W&F), and are as follows: 
!--
!--                               (THICK1)     (THICK2)   
!--                             1000mb-850mb  850mb-700mb  <--THICKNESS
!--
!--   PRECIPITATION TYPE:   RAIN   GT 1310      GT 1540 
!--
!--                FREEZING RAIN   GT 1310      GT 1540 [sig1 T < 0Co]
!--
!--                    ICE/MIXED   LE 1310      GT 1540 
!--
!--                         SNOW   LE 1310      LE 1540. 
!--
!--   Over mountainous terrain, precipitation type is limited to either 
!--   rain or snow.  This is consistent with climatic data presented in 
!--   "A Regional Climatology of Freezing Precipitation for the        
!--   Contiguous United States" (10th Conf. on Applied Climatology, 20-
!--   24 Oct 97).  Where a precipitation occurrence has been determined, 
!--   the temperatures in the lowest 1500 m are checked.  If all are 
!--   below freezing, SNOW is forecasted; otherwise RAIN is predicted.
!--
!--   MODIFICATION:  Added ability to predict regions where thunderstorm
!--   activity may occur.  Prior to exiting the main loop, a check is 
!--   made:  where rain is predicted AND the convective component of the
!--   precip exceeds 0.01", forecast for thunderstorms.  
!--
!--   UPDATES
!--   =======
!--   JAN 2001 1998  INITIAL VERSION, ADAPTED FROM USAF WEATHER AGENCY
!--       Brent Shaw, NOAA/Forecast Systems Lab
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------


      USE constants
      IMPLICIT NONE

      INTEGER                        :: I
      INTEGER,     INTENT(IN)        :: IMAX
      INTEGER                        :: J
      INTEGER,     INTENT(IN)        :: JMAX
      INTEGER                        :: K
      INTEGER,     INTENT(IN)        :: K700
      INTEGER,     INTENT(IN)        :: K850
      INTEGER,     INTENT(IN)        :: K1000
      INTEGER                        :: KPRS
      INTEGER                        :: KSIG
      INTEGER                        :: K1500
      REAL,        INTENT(IN)        :: PRCPCONV    ( imax , jmax )
      REAL,        INTENT(IN)        :: PRCPINC     ( imax , jmax )
      INTEGER,     INTENT(OUT)       :: PRECIPTYPE  ( imax , jmax ) 
      REAL,        INTENT(IN)        :: PSFC        ( imax , jmax )
      REAL,        INTENT(IN)        :: TERRAIN     ( imax , jmax )
      REAL,        INTENT(IN)        :: TSFC        ( imax , jmax )
      REAL                           :: TSFCF    
      REAL,        INTENT(IN)        :: TSIG        ( imax , jmax , ksig )
      REAL                           :: THICKHIGH
      REAL                           :: THICKLOW
      REAL                           :: TSIG1, TSIG2, TSIG3
      REAL,        INTENT(IN)        :: ZPRS        ( imax , jmax , kprs )
      REAL,        INTENT(IN)        :: ZSIG        ( imax , jmax , ksig )
      REAL,        EXTERNAL          :: fahren

      DO J = 1, JMAX
        DO I = 1, IMAX

!-----------------------------------------------------------------------
!--       THE THRESHOLD FOR CALCULATING PRECIP TYPE IS 0.0001 meter
!--       PER TIME PERIOD.
!-----------------------------------------------------------------------

          IF ( PRCPINC(I,J) .LE. 0.0001) THEN

            PRECIPTYPE(I,J) = 0

          ELSE

!-----------------------------------------------------------------------
!--         CHECK THE SURFACE PRESSURE TO DETERMINE WHETHER
!--         HIGH OR LOW-ELEVATION LOGIC IS USED. 850 MB IS
!--         THE CUTOFF.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!--         LOW ELEVATION GRID POINT.
!-----------------------------------------------------------------------

            IF (PSFC(I,J) .GT. 85000.0 ) THEN

!-----------------------------------------------------------------------
!--           CALCULATE THICKNESSES THAT WILL BE USED TO DETERMINE
!--           PRECIP TYPE.
!-----------------------------------------------------------------------

              THICKLOW   = ZPRS(I,J,K850) - ZPRS(I,J,K1000)

              THICKHIGH  = ZPRS(I,J,K700) - ZPRS(I,J,K850)

!-----------------------------------------------------------------------
!--           RAIN, OR IF SURFACE TEMPERATURE IS BELOW FREEZING,
!--           FREEZING RAIN.  
!-----------------------------------------------------------------------

              IF ((THICKLOW .GT. 1310.0) .AND. (THICKHIGH .GT. 1540.0)) THEN

                IF (TSFC(I,J) .GE. T0) THEN  
                  PRECIPTYPE(I,J) = 1  
                ELSE
                  PRECIPTYPE(I,J) = 3  
                ENDIF

!-----------------------------------------------------------------------
!--           ICE/MIXED.
!-----------------------------------------------------------------------

              ELSEIF ((THICKLOW .LE. 1310.0) .AND. &
                      (THICKHIGH .GT. 1540.0)) THEN

                PRECIPTYPE(I,J) = 4  

!-----------------------------------------------------------------------
!--           RAIN OR SNOW.
!-----------------------------------------------------------------------

              ELSEIF ((THICKLOW .LE. 1310.0) .AND. & 
                      (THICKHIGH .LE. 1540.0)) THEN

                TSFCF = fahren(tsfc(i,j)-t0) 

                IF (TSFCF .GE. 37.0) THEN
                  PRECIPTYPE(I,J) = 1    
                ELSE                                  
                  PRECIPTYPE(I,J) = 5   
                ENDIF

!-----------------------------------------------------------------------
!--           RAIN.
!-----------------------------------------------------------------------

              ELSE

                PRECIPTYPE(I,J) = 1

              ENDIF 

!-----------------------------------------------------------------------
!--         HIGH TERRAIN GRID POINT.
!-----------------------------------------------------------------------

            ELSE

!-----------------------------------------------------------------------
!--           FIND TOP OF 1500 M AGL LAYER.
!-----------------------------------------------------------------------
              
              DO K = 1, KSIG

                IF ( (ZSIG(I,J,K) - TERRAIN(I,J)) .GE. 1500.0 ) THEN 
 
                  K1500 = K
                  EXIT
 
                END IF
 
              END DO    
!-----------------------------------------------------------------------
!--           IF THE MODEL TOP IS EVER LOWERED, THE ABOVE CODE
!--           COULD FAIL ON THE TOP OF A HIGH MOUNTAIN.
!-----------------------------------------------------------------------

              K1500 = MAX (K1500, 1)

!-----------------------------------------------------------------------
!--           FIND TEMPERATURE AT THE BOTTOM, TOP AND MIDDLE OF
!--           1500 M AGL LAYER.  FOR MIDDLE LAYER, RECYCLE VARIABLE
!--           K1500.
!-----------------------------------------------------------------------

              TSIG1 = TSIG(I,J,1)

              TSIG3 = TSIG(I,J,K1500)

              K1500 = MAX (1, NINT(FLOAT(K1500)/2.0))

              TSIG2 = TSIG(I,J,K1500)

!----------------------------------------------------------------------
!--           SNOW.
!-----------------------------------------------------------------------

              IF ( (TSIG1 .LT. T0) .AND. &
                   (TSIG2 .LT. T0) .AND. &
                   (TSIG3 .LT. T0) ) THEN 

                PRECIPTYPE(I,J) = 5  

!-----------------------------------------------------------------------
!--           RAIN & CHECK FOR THUNDERSTORMS.
!-----------------------------------------------------------------------

              ELSE

                PRECIPTYPE(I,J) = 1    

              ENDIF 
 
            ENDIF
       
          END IF

          IF ( (PRECIPTYPE(I,J) .EQ. 1) .AND. & 
               (PRCPCONV(I,J) .GE. 0.001) ) THEN

             PRECIPTYPE(I,J) = 2         ! Thunderstorm

          END IF 

        END DO
      END DO
      RETURN
  END SUBROUTINE WINTPREC 

