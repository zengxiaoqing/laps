      SUBROUTINE FI7513 (IWORK,ISTART,NPTS,MAX,MIN,INRNGE)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    FI7513      SELECT BLOCK OF DATA FOR PACKING
C   PRGMMR: CAVANAUGH        ORG: W/NMC42    DATE: 94-01-21
C
C ABSTRACT: SELECT A BLOCK OF DATA FOR PACKING
C
C PROGRAM HISTORY LOG:
C   94-01-21  CAVANAUGH
C   95-10-31  IREDELL     REMOVED SAVES AND PRINTS
C
C USAGE:    CALL FI7513 (IWORK,ISTART,NPTS,MAX,MIN,INRNGE)
C   INPUT ARGUMENT LIST:
C     *        - RETURN ADDRESS IF ENCOUNTER SET OF SAME VALUES
C     IWORK    -
C     ISTART   -
C     NPTS     -
C
C   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
C     MAX      -
C     MIN      -
C     INRNGE   -
C
C REMARKS: SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSING ENVIRONMENT.
C
C ATTRIBUTES:
C   LANGUAGE: IBM VS FORTRAN 77, CRAY CFT77 FORTRAN
C   MACHINE:  HDS, CRAY C916/256, Y-MP8/64, Y-MP EL92/256
C
C$$$
      INTEGER        IWORK(*),NPTS,ISTART,INRNGE,INRNGA,INRNGB
      INTEGER        MAX,MIN,MXVAL,MAXB,MINB,MXVALB
      INTEGER        IBITS(31)
C
      DATA           IBITS/1,3,7,15,31,63,127,255,511,1023,2047,
     *               4095,8191,16383,32767,65535,131071,262143,
     *               524287,1048575,2097151,4194303,8388607,
     *               16777215,33554431,67108863,134217727,268435455,
     *               536870911,1073741823,2147483647/
C  ----------------------------------------------------------------
C                        IDENTIFY NEXT BLOCK OF DATA FOR PACKING AND
C                           RETURN TO CALLER
C  ********************************************************************
      ISTRTA  = ISTART
C
C                     GET BLOCK A
      CALL FI7516 (IWORK,NPTS,INRNGA,ISTRTA,
     *                                  MAX,MIN,MXVAL,LWIDE)
C  ********************************************************************
C
      ISTRTB  = ISTRTA + INRNGA
 2000 CONTINUE
C                         IF HAVE PROCESSED ALL DATA, RETURN
      IF (ISTRTB.GT.NPTS) THEN
C                         NO MORE DATA TO LOOK AT
          INRNGE  = INRNGA
          RETURN
      END IF
C                     GET BLOCK B
      CALL FI7502 (IWORK,ISTRTB,NPTS,ISAME)
      IF (ISAME.GE.15) THEN
C         PRINT *,'BLOCK B HAS ALL IDENTICAL VALUES'
C         PRINT *,'BLOCK A HAS INRNGE =',INRNGA
C                     BLOCK B CONTAINS ALL IDENTICAL VALUES
          INRNGE  = INRNGA
C                     EXIT WITH BLOCK A
          RETURN
      END IF
C                     GET BLOCK B
C
      ISTRTB  = ISTRTA + INRNGA
      CALL FI7516 (IWORK,NPTS,INRNGB,ISTRTB,
     *                                  MAXB,MINB,MXVALB,LWIDEB)
C     PRINT *,'BLOCK A',INRNGA,' BLOCK B',INRNGB
C  ********************************************************************
C                     PERFORM TREND ANALYSIS TO DETERMINE
C                     IF DATA COLLECTION CAN BE IMPROVED
C
      KTRND  = LWIDE - LWIDEB
C     PRINT *,'TREND',LWIDE,LWIDEB
      IF (KTRND.LE.0) THEN
C         PRINT *,'BLOCK A - SMALLER, SHOULD EXTEND INTO BLOCK B'
          MXVAL   = IBITS(LWIDE)
C
C                     IF BLOCK A REQUIRES THE SAME OR FEWER BITS
C                             LOOK AHEAD
C                        AND GATHER THOSE DATA POINTS THAT CAN
C                        BE RETAINED IN BLOCK A
C                        BECAUSE THIS BLOCK OF DATA
C                            USES FEWER BITS
C
          CALL FI7518 (IRET,IWORK,NPTS,ISTRTA,INRNGA,INRNGB,
     *                          MAX,MIN,LWIDE,MXVAL)
          IF(IRET.EQ.1) GO TO 8000
C         PRINT *,'18 INRNGA IS NOW ',INRNGA
          IF (INRNGB.LT.20) THEN
              RETURN
          ELSE
              GO TO 2000
          END IF
      ELSE
C         PRINT *,'BLOCK A - LARGER, B SHOULD EXTEND BACK INTO A'
          MXVALB  = IBITS(LWIDEB)
C
C                     IF BLOCK B REQUIRES FEWER BITS
C                             LOOK BACK
C                            SHORTEN BLOCK A BECAUSE NEXT BLOCK OF DATA
C                            USES FEWER BITS
C
          CALL FI7517 (IRET,IWORK,NPTS,ISTRTB,INRNGA,
     *                               MAXB,MINB,LWIDEB,MXVALB)
          IF(IRET.EQ.1) GO TO 8000
C         PRINT *,'17 INRNGA IS NOW ',INRNGA
      END IF
C
C                           PACK UP BLOCK A
C                           UPDATA POINTERS
 8000 CONTINUE
      INRNGE  = INRNGA
C                           GET NEXT BLOCK A
 9000 CONTINUE
      RETURN
      END
