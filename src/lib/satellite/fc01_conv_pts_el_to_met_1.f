C
C 5-11-99:  J. Smart  AFWA code modified for inclusion into LAPS satellite
C                     library to compute METEOSAT line/pixel values for
c                     given LAPS lat/lon.
C   4.4.1.2.2.11.7.1 FC01_CONV_PTS_EL_TO_MET_1 
CA
C   FC01_CONV_PTS_EL_TO_MET_1 converts latitude/longitude in radians to
C   pixel scanline number using algorithm described in annex "E"
C   to the meteosat-5 calibration report (April 1991). prepared by the
C   meteosat exploitation project at european space operations center.
CE
C---------------------  ARGUMENTS PASSED  ---------------------------
CCA
C
C      NAME                    I/O        DESCRIPTION
C
C      LATITUDE                 I   INPUT LATITUDE
C
C      LONGITUDE                I   INPUT LONGITUDE
C
C      X_COORD                  O   OUTPUT X COORDINATE
C
C      Y_COORD                  O   OUTPUT Y COORDINATE
C
C      CSTATUS                  O   ERROR CODE C
CCAE
C-------------------  SOFTWARE UNIT HISTORY  ------------------------
CH
C
C       06DEC95 : created, SCR95124, sfl;
C       26FEB97 : renamed RS to RS1 to avoid conflict with Ingest
C                  include file, sdhsup3, kla;
C
C      25MAR97 : SDHSU   Placed key in comments for
C                        SDD Sec 4 generation; RFB
CHE

C--------------------------------------------------------------------
C
C     Procedure FC01_CONV_PTS_EL_TO_MET_1 is
C
CP    begin
CE    end FC01_CONV_PTS_EL_TO_MET_1 ;
C
C--------------------------------------------------------------------

      SUBROUTINE FC01_CONV_PTS_EL_TO_MET_1(golonsbp,
     &                          stop_pix,fsci,
     &                          start_line,decimat,
     &                          ct,
     &                          LATITUDE,
     &                          LONGITUDE,
     &                          X_COORD,
     &                          Y_COORD,
     &                          CSTATUS)

      IMPLICIT NONE

c     INCLUDE 'data_type_code_tables.inc'
c     INCLUDE 'message_data_variables.inc'
c     INCLUDE 'is1_ingest_parameters_5.inc'

      CHARACTER ct*(*)

      REAL   LATITUDE,
     &       LONGITUDE,
     &       X_COORD,
     &       Y_COORD

      INTEGER   CSTATUS

      REAL   H
      PARAMETER (H = 35785.845)
      REAL   RE
      PARAMETER (RE = 6378.155)
      REAL   A
      PARAMETER (A = 1.0/297.0)
      REAL   RP
      PARAMETER (RP = RE/(1.0+A))
      REAL   PI
      PARAMETER (PI = 3.141592653)
      REAL   CDR
      PARAMETER (CDR = PI/180.0)
      REAL   CRD 
      PARAMETER (CRD = 180.0/PI)
      INTEGER   LPSI2
      PARAMETER (LPSI2 = 1)
      REAL   DELTAX
      PARAMETER (DELTAX = 18.0/2500.0)
      REAL   DELTAY 
      PARAMETER (DELTAY = 18.0/2500.0)

      REAL   XFI,
     &       XLA,
     &       ROM,
     &       Y,
     &       R1,
     &       R2,
     &       RS1,
     &       XR,
     &       YR,
     &       XT,
     &       YT,
     &       ZT,
     &       PX,
     &       PY,
     &       TETA
      REAL   RADJUST
      Real      golonsbp
      Integer   start_line
      Integer   stop_pix
      Integer   fsci
      Integer   decimat
      
c     INCLUDE 'data_areas_imd_hdr_offsets.inc'
c     INCLUDE 'data_areas_imd_hdr_data_names.inc'

c     COMMON /BLOCK_13/ IMAGE_LIST_CHAR

C-----------------------------------------------------------
C--   BEGIN
C-----------------------------------------------------------
C#
      CSTATUS = 0
      RADJUST=1.0
      IF(CT.EQ.'vis')RADJUST=2.0
      XFI = LATITUDE * CDR
      XLA = (LONGITUDE - GOLONSBP) * CDR
      ROM = (RE * RP)/SQRT(RP**2 * COS(XFI)**2 + RE**2 * SIN(XFI)**2)
      Y = SQRT(H**2 + ROM**2 - 2 * H * ROM * COS(XFI) * COS(XLA))
      R1 = Y**2 + ROM**2
      R2 = H**2

      IF (R1 .GT. R2) THEN
         X_COORD = 9999.
         Y_COORD = 9999.
         CSTATUS = 1
      ELSE
         RS1 = RE + H
         TETA = ATAN((RP / RE) * TAN(XFI))
         XT = RE * COS(TETA) * COS(XLA)
         YT = RE * COS(TETA) * SIN(XLA)
         ZT = RP * SIN(TETA)

         PX = ATAN(YT / (XT-RS1))

         PY = ATAN(ZT * (-1.0 / (XT-RS1))
     &                * COS(PX))

         PX = PX * CRD
         PY = PY * CRD

         XR = PX / (DELTAX * LPSI2)
         YR = PY / (DELTAY * LPSI2)

         IF(XR .GE. 0.0) XR = INT(PX / (DELTAX * LPSI2)) + 0.5
         IF(XR .LT. 0.0) XR = INT(PX / (DELTAX * LPSI2)) - 0.5
         IF(YR .GE. 0.0) YR = INT(PY / (DELTAY * LPSI2)) + 0.5
         IF(YR .LT. 0.0) YR = INT(PY / (DELTAY * LPSI2)) - 0.5

         X_COORD = (XR + 1250.0)*radjust
         Y_COORD = (YR + 1251.0)*radjust
c
c compute file relative x/y coordinate positions.
c
         X_COORD = X_COORD/DECIMAT
         X_COORD = STOP_PIX-X_COORD

         Y_COORD = Y_COORD - FSCI
         Y_COORD = Y_COORD/DECIMAT + 0.5
         Y_COORD = Y_COORD - START_LINE

C        Y_COORD = ((Y_COORD-(FSCI/radjust))/DECIMAT)+0.5-START_LINE

      END IF

99999 RETURN
      END
