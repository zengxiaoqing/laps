      SUBROUTINE FI7517 (IRET,IWORK,NPTS,ISTRTB,INRNGA,
     *                           MAXB,MINB,MXVALB,LWIDEB)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    FI7517      SCAN BACKWARD
C   PRGMMR: CAVANAUGH        ORG: W/NMC42    DATE: 94-01-21
C
C ABSTRACT: SCAN BACKWARDS UNTIL A VALUE EXCEEDS RANGE OF GROUP B
C           THIS MAY SHORTEN GROUP A
C
C PROGRAM HISTORY LOG:
C   94-01-21  CAVANAUGH
C   95-10-31  IREDELL     REMOVED SAVES AND PRINTS
C   98-06-17  IREDELL     REMOVED ALTERNATE RETURN
C
C USAGE:    CALL FI7517 (IRET,IWORK,NPTS,ISTRTB,INRNGA,
C    *                           MAXB,MINB,MXVALB,LWIDEB)
C   INPUT ARGUMENT LIST:
C     IWORK    -
C     ISTRTB   -
C     NPTS     -
C     INRNGA   -
C
C   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
C     IRET     -
C     JLAST    -
C     MAXB     -
C     MINB     -
C     LWIDTH   - NUMBER OF BITS TO CONTAIN MAX DIFF
C
C REMARKS: SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSING ENVIRONMENT.
C
C ATTRIBUTES:
C   LANGUAGE: IBM VS FORTRAN 77, CRAY CFT77 FORTRAN
C   MACHINE:  HDS, CRAY C916/256, Y-MP8/64, Y-MP EL92/256
C
C$$$
      INTEGER        IWORK(*),NPTS,ISTRTB,INRNGA
      INTEGER        MAXB,MINB,LWIDEB,MXVALB
      INTEGER        IBITS(31)
C
      DATA           IBITS/1,3,7,15,31,63,127,255,511,1023,2047,
     *               4095,8191,16383,32767,65535,131071,262143,
     *               524287,1048575,2097151,4194303,8388607,
     *               16777215,33554431,67108863,134217727,268435455,
     *               536870911,1073741823,2147483647/
C  ----------------------------------------------------------------
      IRET=0
C     PRINT *,'          FI7517'
      NPOS  = ISTRTB - 1
      ITST  = 0
      KSET  = INRNGA
C
 1000 CONTINUE
C     PRINT *,'TRY NPOS',NPOS,IWORK(NPOS),MAXB,MINB
      ITST  = ITST + 1
      IF (ITST.LE.KSET) THEN
          IF (IWORK(NPOS).GT.MAXB) THEN
              IF ((IWORK(NPOS)-MINB).GT.MXVALB) THEN
C                 PRINT *,'WENT OUT OF RANGE AT',NPOS
                  IRET=1
                  RETURN
              ELSE
                  MAXB    = IWORK(NPOS)
              END IF
          ELSE IF (IWORK(NPOS).LT.MINB) THEN
              IF ((MAXB-IWORK(NPOS)).GT.MXVALB) THEN
C                 PRINT *,'WENT OUT OF RANGE AT',NPOS
                  IRET=1
                  RETURN
              ELSE
                  MINB    = IWORK(NPOS)
              END IF
          END IF
          INRNGA  = INRNGA - 1
          NPOS  = NPOS - 1
          GO TO 1000
      END IF
C  ----------------------------------------------------------------
C
 9000 CONTINUE
      RETURN
      END
