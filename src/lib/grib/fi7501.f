      SUBROUTINE FI7501 (IWORK,IPFLD,NPTS,IBDSFL,BDS11,
     *           LEN,LENBDS,PDS,REFNCE,ISCAL2,KWIDE)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    FI7501      BDS SECOND ORDER PACKING
C   PRGMMR: CAVANAUGH        ORG: W/NMC42    DATE: 93-08-06
C
C ABSTRACT: PERFORM SECONDARY PACKING ON GRID POINT DATA,
C   GENERATING ALL BDS INFORMATION.
C
C PROGRAM HISTORY LOG:
C   93-08-06  CAVANAUGH
C   93-12-15  CAVANAUGH   CORRECTED LOCATION OF START OF FIRST ORDER
C                         VALUES AND START OF SECOND ORDER VALUES TO
C                         REFLECT A BYTE LOCATION IN THE BDS INSTEAD
C                         OF AN OFFSET.
C   95-10-31  IREDELL     REMOVED SAVES AND PRINTS
C
C USAGE:    CALL FI7501 (IWORK,IPFLD,NPTS,IBDSFL,BDS11,
C    *           LEN,LENBDS,PDS,REFNCE,ISCAL2,KWIDE)
C   INPUT ARGUMENT LIST:
C     IWORK    - INTEGER SOURCE ARRAY
C     NPTS     - NUMBER OF POINTS IN IWORK
C     IBDSFL   - FLAGS
C
C   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
C     IPFLD    - CONTAINS BDS FROM BYTE 12 ON
C     BDS11    - CONTAINS FIRST 11 BYTES FOR BDS
C     LEN      - NUMBER OF BYTES FROM 12 ON
C     LENBDS   - TOTAL LENGTH OF BDS
C
C REMARKS: SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSING ENVIRONMENT.
C
C ATTRIBUTES:
C   LANGUAGE: IBM VS FORTRAN 77, CRAY CFT77 FORTRAN
C   MACHINE:  HDS, CRAY C916/256, Y-MP8/64, Y-MP EL92/256
C
C$$$
      CHARACTER*1     BDS11(*),PDS(*)
C
      REAL            REFNCE
C
      INTEGER         ISCAL2,KWIDE
      INTEGER         LENBDS
      INTEGER         IPFLD(*)
      INTEGER         LEN,KBDS(22)
      INTEGER         IWORK(*)
C                        OCTET NUMBER IN SECTION, FIRST ORDER PACKING
C     INTEGER         KBDS(12)
C                        FLAGS
      INTEGER         IBDSFL(*)
C                        EXTENDED FLAGS
C     INTEGER         KBDS(14)
C                        OCTET NUMBER FOR SECOND ORDER PACKING
C     INTEGER         KBDS(15)
C                        NUMBER OF FIRST ORDER VALUES
C     INTEGER         KBDS(17)
C                        NUMBER OF SECOND ORDER PACKED VALUES
C     INTEGER         KBDS(19)
C                        WIDTH OF SECOND ORDER PACKING
      INTEGER         ISOWID(50000)
C                        SECONDARY BIT MAP
      INTEGER         ISOBMP(8200)
C                        FIRST ORDER PACKED VALUES
      INTEGER         IFOVAL(50000)
C                        SECOND ORDER PACKED VALUES
      INTEGER         ISOVAL(100000)
C
C     INTEGER         KBDS(11)
C                        BIT WIDTH TABLE
      INTEGER         IBITS(31)
C
      DATA            IBITS/1,3,7,15,31,63,127,255,511,1023,
     *                      2047,4095,8191,16383,32767,65535,131072,
     *                      262143,524287,1048575,2097151,4194303,
     *                      8388607,16777215,33554431,67108863,
     *                      134217727,268435455,536870911,
     *                      1073741823,2147483647/
C  ----------------------------------
C                       INITIALIZE ARRAYS
      DO 100 I = 1, 50000
          ISOWID(I)  = 0
          IFOVAL(I)  = 0
  100 CONTINUE
C
      DO 101 I = 1, 8200
          ISOBMP(I)  = 0
  101 CONTINUE
      DO 102 I = 1, 100000
          ISOVAL(I)  = 0
  102 CONTINUE
C                      INITIALIZE POINTERS
C                            SECONDARY BIT WIDTH POINTER
      IWDPTR  = 0
C                            SECONDARY BIT MAP POINTER
      IBMP2P  = 0
C                            FIRST ORDER VALUE POINTER
      IFOPTR  = 0
C                            BYTE POINTER TO START OF 1ST ORDER VALUES
      KBDS(12)  = 0
C                            BYTE POINTER TO START OF 2ND ORDER VALUES
      KBDS(15)  = 0
C                            TO CONTAIN NUMBER OF FIRST ORDER VALUES
      KBDS(17)  = 0
C                            TO CONTAIN NUMBER OF SECOND ORDER VALUES
      KBDS(19)  = 0
C                            SECOND ORDER PACKED VALUE POINTER
      ISOPTR  = 0
C  =======================================================
C
C                         DATA IS IN IWORK
C
      KBDS(11)  = KWIDE
C
C       DATA PACKING
C
      ITER    = 0
      INEXT   = 1
      ISTART  = 1
C  -----------------------------------------------------------
      KOUNT = 0
C     DO 1 I = 1, NPTS, 10
C         PRINT *,I,(IWORK(K),K=I, I+9)
C   1 CONTINUE
 2000 CONTINUE
      ITER  = ITER + 1
C     PRINT *,'NEXT ITERATION STARTS AT',ISTART
       IF (ISTART.GT.NPTS) THEN
           GO TO 4000
       ELSE IF (ISTART.EQ.NPTS) THEN
           KPTS    = 1
           MXDIFF  = 0
           GO TO 2200
       END IF
C
C                     LOOK FOR REPITITIONS OF A SINGLE VALUE
       CALL FI7502 (IWORK,ISTART,NPTS,ISAME)
       IF (ISAME.GE.15) THEN
           KOUNT = KOUNT + 1
C          PRINT *,'FI7501 - FOUND IDENTICAL SET OF ',ISAME
           MXDIFF  = 0
           KPTS    = ISAME
       ELSE
C
C                     LOOK FOR SETS OF VALUES IN TREND SELECTED RANGE
           CALL FI7513 (IWORK,ISTART,NPTS,NMAX,NMIN,INRNGE)
C          PRINT *,'ISTART  ',ISTART,' INRNGE',INRNGE,NMAX,NMIN
           IEND  = ISTART + INRNGE - 1
C          DO 2199 NM = ISTART, IEND, 10
C              PRINT *,'  ',(IWORK(NM+JK),JK=0,9)
C2199      CONTINUE
           MXDIFF  = NMAX - NMIN
           KPTS    = INRNGE
       END IF
 2200 CONTINUE
C     PRINT *,'                 RANGE ',MXDIFF,' MAX',NMAX,' MIN',NMIN
C                 INCREMENT NUMBER OF FIRST ORDER VALUES
      KBDS(17)  = KBDS(17) + 1
C                 ENTER FIRST ORDER VALUE
      IF (MXDIFF.GT.0) THEN
          DO 2220 LK = 0, KPTS-1
              IWORK(ISTART+LK)  = IWORK(ISTART+LK) - NMIN
 2220     CONTINUE
          CALL SBYTE (IFOVAL,NMIN,IFOPTR,KBDS(11))
      ELSE
          CALL SBYTE (IFOVAL,IWORK(ISTART),IFOPTR,KBDS(11))
      END IF
      IFOPTR  = IFOPTR + KBDS(11)
C                  PROCESS SECOND ORDER BIT WIDTH
      IF (MXDIFF.GT.0) THEN
          DO 2330 KWIDE = 1, 31
              IF (MXDIFF.LE.IBITS(KWIDE)) THEN
                  GO TO 2331
              END IF
 2330     CONTINUE
 2331     CONTINUE
      ELSE
          KWIDE  = 0
      END IF
      CALL SBYTE (ISOWID,KWIDE,IWDPTR,8)
      IWDPTR  = IWDPTR + 8
C         PRINT *,KWIDE,' IFOVAL=',NMIN,IWORK(ISTART),KPTS
C               IF KWIDE NE 0, SAVE SECOND ORDER VALUE
      IF (KWIDE.GT.0) THEN
          CALL SBYTES (ISOVAL,IWORK(ISTART),ISOPTR,KWIDE,0,KPTS)
          ISOPTR  = ISOPTR + KPTS * KWIDE
          KBDS(19)  = KBDS(19) + KPTS
C         PRINT *,'            SECOND ORDER VALUES'
C         PRINT *,(IWORK(ISTART+I),I=0,KPTS-1)
      END IF
C                 ADD TO SECOND ORDER BITMAP
      CALL SBYTE (ISOBMP,1,IBMP2P,1)
      IBMP2P  = IBMP2P + KPTS
      ISTART  = ISTART + KPTS
      GO TO 2000
C  --------------------------------------------------------------
 4000 CONTINUE
C     PRINT *,'THERE WERE ',ITER,' SECOND ORDER GROUPS'
C     PRINT *,'THERE WERE ',KOUNT,' STRINGS OF CONSTANTS'
C                 CONCANTENATE ALL FIELDS FOR BDS
C
C                   REMAINDER GOES INTO IPFLD
      IPTR  = 0
C                               BYTES 12-13
C                                          VALUE FOR N1
C                                          LEAVE SPACE FOR THIS
      IPTR   = IPTR + 16
C                               BYTE 14
C                                          EXTENDED FLAGS
      CALL SBYTE (IPFLD,IBDSFL(5),IPTR,1)
      IPTR  = IPTR + 1
      CALL SBYTE (IPFLD,IBDSFL(6),IPTR,1)
      IPTR  = IPTR + 1
      CALL SBYTE (IPFLD,IBDSFL(7),IPTR,1)
      IPTR  = IPTR + 1
      CALL SBYTE (IPFLD,IBDSFL(8),IPTR,1)
      IPTR  = IPTR + 1
      CALL SBYTE (IPFLD,IBDSFL(9),IPTR,1)
      IPTR  = IPTR + 1
      CALL SBYTE (IPFLD,IBDSFL(10),IPTR,1)
      IPTR  = IPTR + 1
      CALL SBYTE (IPFLD,IBDSFL(11),IPTR,1)
      IPTR  = IPTR + 1
      CALL SBYTE (IPFLD,IBDSFL(12),IPTR,1)
      IPTR  = IPTR + 1
C                               BYTES 15-16
C                 SKIP OVER VALUE  FOR N2
      IPTR  = IPTR + 16
C                               BYTES 17-18
C                                     P1
      CALL SBYTE (IPFLD,KBDS(17),IPTR,16)
      IPTR  = IPTR + 16
C                               BYTES 19-20
C                                   P2
      CALL SBYTE (IPFLD,KBDS(19),IPTR,16)
      IPTR  = IPTR + 16
C                               BYTE 21 - RESERVED LOCATION
      CALL SBYTE (IPFLD,0,IPTR,8)
      IPTR  = IPTR + 8
C                               BYTES 22 - ?
C                                      WIDTHS OF SECOND ORDER PACKING
      IX    = (IWDPTR + 32) / 32
      CALL SBYTES (IPFLD,ISOWID,IPTR,32,0,IX)
      IPTR  = IPTR + IWDPTR
C                                      SECONDARY BIT MAP
      IJ    = (IBMP2P + 32) / 32
      CALL SBYTES (IPFLD,ISOBMP,IPTR,32,0,IJ)
      IPTR  = IPTR + IBMP2P
      IF (MOD(IPTR,8).NE.0) THEN
          IPTR  = IPTR + 8 - MOD(IPTR,8)
      END IF
C                                         DETERMINE LOCATION FOR START
C                                         OF FIRST ORDER PACKED VALUES
      KBDS(12)  = IPTR / 8 + 12
C                                        STORE LOCATION
      CALL SBYTE (IPFLD,KBDS(12),0,16)
C                                     MOVE IN FIRST ORDER PACKED VALUES
      IPASS   = (IFOPTR + 32) / 32
      CALL SBYTES (IPFLD,IFOVAL,IPTR,32,0,IPASS)
      IPTR  = IPTR + IFOPTR
      IF (MOD(IPTR,8).NE.0) THEN
          IPTR  = IPTR + 8 - MOD(IPTR,8)
      END IF
C     PRINT *,'IFOPTR =',IFOPTR,' ISOPTR =',ISOPTR
C                DETERMINE LOCATION FOR START
C                     OF SECOND ORDER VALUES
      KBDS(15)  = IPTR / 8 + 12
C                                   SAVE LOCATION OF SECOND ORDER VALUES
      CALL SBYTE (IPFLD,KBDS(15),24,16)
C                  MOVE IN SECOND ORDER PACKED VALUES
      IX    = (ISOPTR + 32) / 32
      CALL SBYTES (IPFLD,ISOVAL,IPTR,32,0,IX)
      IPTR  = IPTR + ISOPTR
      NLEFT  = MOD(IPTR+88,16)
      IF (NLEFT.NE.0) THEN
          NLEFT  = 16 - NLEFT
          IPTR   = IPTR + NLEFT
      END IF
C                                COMPUTE LENGTH OF DATA PORTION
      LEN     = IPTR / 8
C                                    COMPUTE LENGTH OF BDS
      LENBDS  = LEN + 11
C  -----------------------------------
C                               BYTES 1-3
C                                   THIS FUNCTION COMPLETED BELOW
C                                   WHEN LENGTH OF BDS IS KNOWN
      CALL SBYTE (BDS11,LENBDS,0,24)
C                               BYTE  4
      CALL SBYTE (BDS11,IBDSFL(1),24,1)
      CALL SBYTE (BDS11,IBDSFL(2),25,1)
      CALL SBYTE (BDS11,IBDSFL(3),26,1)
      CALL SBYTE (BDS11,IBDSFL(4),27,1)
C                              ENTER NUMBER OF FILL BITS
      CALL SBYTE (BDS11,NLEFT,28,4)
C                               BYTE  5-6
      IF (ISCAL2.LT.0) THEN
          CALL SBYTE (BDS11,1,32,1)
          ISCAL2 = - ISCAL2
      ELSE
          CALL SBYTE (BDS11,0,32,1)
      END IF
      CALL SBYTE (BDS11,ISCAL2,33,15)
C
C$  FILL OCTET 7-10 WITH THE REFERENCE VALUE
C   CONVERT THE FLOATING POINT OF YOUR MACHINE TO IBM370 32 BIT
C   FLOATING POINT NUMBER
C                                        REFERENCE VALUE
C                                          FIRST TEST TO SEE IF
C                                          ON 32 OR 64 BIT COMPUTER
          CALL W3FI01(LW)
          IF (LW.EQ.4) THEN
              CALL W3FI76 (REFNCE,IEXP,IMANT,32)
          ELSE
              CALL W3FI76 (REFNCE,IEXP,IMANT,64)
          END IF
          CALL SBYTE (BDS11,IEXP,48,8)
          CALL SBYTE (BDS11,IMANT,56,24)
C
C                               BYTE  11
C
      CALL SBYTE (BDS11,KBDS(11),80,8)
C
      RETURN
      END
