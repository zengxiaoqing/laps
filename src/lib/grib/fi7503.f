      SUBROUTINE FI7503 (IWORK,IPFLD,NPTS,IBDSFL,BDS11,
     *           LEN,LENBDS,PDS,REFNCE,ISCAL2,KWIDE,IGDS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:    FI7501      ROW BY ROW, COL BY COL PACKING
C   PRGMMR: CAVANAUGH        ORG: W/NMC42    DATE: 94-05-20
C
C ABSTRACT: PERFORM ROW BY ROW OR COLUMN BY COLUMN PACKING
C   GENERATING ALL BDS INFORMATION.
C
C PROGRAM HISTORY LOG:
C   93-08-06  CAVANAUGH
C   95-10-31  IREDELL     REMOVED SAVES AND PRINTS
C
C USAGE:    CALL FI7503 (IWORK,IPFLD,NPTS,IBDSFL,BDS11,
C    *           LEN,LENBDS,PDS,REFNCE,ISCAL2,KWIDE,IGDS)
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
      INTEGER         IPFLD(*),IGDS(*)
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
C                         BUILD SECOND ORDER BIT MAP IN EITHER
C                         ROW BY ROW OR COL BY COL FORMAT
      IF (IAND(IGDS(13),32).NE.0) THEN
C                              COLUMN BY COLUMN
          KOUT  = IGDS(4)
          KIN   = IGDS(5)
C         PRINT *,'COLUMN BY COLUMN',KOUT,KIN
      ELSE
C                              ROW BY ROW
          KOUT  = IGDS(5)
          KIN   = IGDS(4)
C         PRINT *,'ROW BY ROW',KOUT,KIN
      END IF
      KBDS(17)  = KOUT
      KBDS(19)  = NPTS
C
C     DO 4100 J = 1, NPTS, 53
C         WRITE (6,4101) (IWORK(K),K=J,J+52)
 4101     FORMAT (1X,25I4)
C         PRINT *,' '
C4100 CONTINUE
C
C                             INITIALIZE BIT MAP POINTER
      IBMP2P = 0
C                             CONSTRUCT WORKING BIT MAP
      DO 2000 I = 1, KOUT
          DO 1000 J = 1, KIN
              IF (J.EQ.1) THEN
                  CALL SBYTE (ISOBMP,1,IBMP2P,1)
              ELSE
                  CALL SBYTE (ISOBMP,0,IBMP2P,1)
              END IF
              IBMP2P  = IBMP2P + 1
 1000     CONTINUE
 2000 CONTINUE
      LEN  = IBMP2P / 32 + 1
C     CALL BINARY(ISOBMP,LEN)
C
C                       PROCESS OUTER LOOP OF ROW BY ROW OR COL BY COL
C
      KPTR  = 1
      KBDS(11)  = KWIDE
      DO 6000 I = 1, KOUT
C                       IN CURRENT ROW OR COL
C                              FIND FIRST ORDER VALUE
          JPTR  = KPTR
          LOWEST  = IWORK(JPTR)
          DO 4000 J = 1, KIN
              IF (IWORK(JPTR).LT.LOWEST) THEN
                  LOWEST = IWORK(JPTR)
              END IF
              JPTR  = JPTR + 1
 4000     CONTINUE
C                            SAVE FIRST ORDER VALUE
          CALL SBYTE (IFOVAL,LOWEST,IFOPTR,KWIDE)
          IFOPTR  = IFOPTR + KWIDE
C         PRINT *,'FOVAL',I,LOWEST,KWIDE
C                            SUBTRACT FIRST ORDER VALUE FROM OTHER VALS
C                                         GETTING SECOND ORDER VALUES
          JPTR  = KPTR
          IBIG  = IWORK(JPTR) - LOWEST
          DO 4200 J = 1, KIN
              IWORK(JPTR)  = IWORK(JPTR) - LOWEST
              IF (IWORK(JPTR).GT.IBIG) THEN
                  IBIG  = IWORK(JPTR)
              END IF
              JPTR  = JPTR + 1
 4200     CONTINUE
C                            HOW MANY BITS TO CONTAIN LARGEST SECOND
C                                         ORDER VALUE IN SEGMENT
          CALL FI7505 (IBIG,NWIDE)
C                            SAVE BIT WIDTH
          CALL SBYTE (ISOWID,NWIDE,IWDPTR,8)
          IWDPTR  = IWDPTR + 8
C         PRINT *,I,'SOVAL',IBIG,' IN',NWIDE,' BITS'
C         WRITE (6,4101) (IWORK(K),K=KPTR,KPTR+52)
C                            SAVE SECOND ORDER VALUES OF THIS SEGMENT
          DO 5000 J = 0, KIN-1
              CALL SBYTE (ISOVAL,IWORK(KPTR+J),ISOPTR,NWIDE)
              ISOPTR  = ISOPTR + NWIDE
 5000     CONTINUE
          KPTR    = KPTR + KIN
 6000 CONTINUE
C  =======================================================
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
C     PRINT *,'ISOWID',IWDPTR,IX
C     CALL BINARY (ISOWID,IX)
C
C                     NO SECONDARY BIT MAP

C                                         DETERMINE LOCATION FOR START
C                                         OF FIRST ORDER PACKED VALUES
      KBDS(12)  = IPTR / 8 + 12
C                                        STORE LOCATION
      CALL SBYTE (IPFLD,KBDS(12),0,16)
C                                     MOVE IN FIRST ORDER PACKED VALUES
      IPASS   = (IFOPTR + 32) / 32
      CALL SBYTES (IPFLD,IFOVAL,IPTR,32,0,IPASS)
      IPTR  = IPTR + IFOPTR
C     PRINT *,'IFOVAL',IFOPTR,IPASS,KWIDE
C     CALL BINARY (IFOVAL,IPASS)
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
C     PRINT *,'ISOVAL',ISOPTR,IX
C     CALL BINARY (ISOVAL,IX)
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
      KLEN  = LENBDS / 4 + 1
C     PRINT *,'BDS11 LISTING',4,LENBDS
C     CALL BINARY (BDS11,4)
C     PRINT *,'IPFLD LISTING'
C     CALL BINARY (IPFLD,KLEN)
      RETURN
      END
