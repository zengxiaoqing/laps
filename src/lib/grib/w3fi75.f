      SUBROUTINE W3FI75 (IBITL,ITYPE,ITOSS,FLD,IFLD,IBMAP,IBDSFL,
     &  NPTS,BDS11,IPFLD,PFLD,LEN,LENBDS,IBERR,PDS,IGDS)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM:  W3FI75        GRIB PACK DATA AND FORM BDS OCTETS(1-11)
C   PRGMMR: FARLEY           ORG: NMC421      DATE:94-11-22
C
C ABSTRACT: THIS ROUTINE PACKS A GRIB FIELD AND FORMS OCTETS(1-11)
C   OF THE BINARY DATA SECTION (BDS).
C
C PROGRAM HISTORY LOG:
C   92-07-10  M. FARLEY   ORIGINAL AUTHOR
C   92-10-01  R.E.JONES   CORRECTION FOR FIELD OF CONSTANT DATA
C   92-10-16  R.E.JONES   GET RID OF ARRAYS FP AND INT
C   93-08-06  CAVANAUGH   ADDED ROUTINES FI7501, FI7502, FI7503
C                         TO ALLOW SECOND ORDER PACKING IN PDS.
C   93-07-21 STACKPOLE    ASSORTED REPAIRS TO GET 2ND DIFF PACK IN
C   93-10-28 CAVANAUGH    COMMENTED OUT NONOPERATIONAL PRINTS AND
C                         WRITE STATEMENTS
C   93-12-15  CAVANAUGH   CORRECTED LOCATION OF START OF FIRST ORDER
C                         VALUES AND START OF SECOND ORDER VALUES TO
C                         REFLECT A BYTE LOCATION IN THE BDS INSTEAD
C                         OF AN OFFSET IN SUBROUTINE FI7501.
C   94-01-27  CAVANAUGH   ADDED IGDS AS INPUT ARGUMENT TO THIS ROUTINE
C                         AND ADDED PDS AND IGDS ARRAYS TO THE CALL TO
C                         W3FI82 TO PROVIDE INFORMATION NEEDED FOR
C                         BOUSTROPHEDONIC PROCESSING.
C   94-05-25  CAVANAUGH   SUBROUTINE FI7503 HAS BEEN ADDED TO PROVIDE
C                         FOR ROW BY ROW OR COLUMN BY COLUMN SECOND
C                         ORDER PACKING.  THIS FEATURE CAN BE ACTIVATED
C                         BY SETTING IBDSFL(7) TO ZERO.
C   94-07-08  CAVANAUGH   COMMENTED OUT PRINT STATEMENTS USED FOR DEBUG
C   94-11-22  FARLEY      ENLARGED WORK ARRAYS TO HANDLE .5DEGREE GRIDS
C   95-06-01  R.E.JONES   CORRECTION FOR NUMBER OF UNUSED BITS AT END
C                         OF SECTION 4, IN BDS BYTE 4, BITS 5-8.
C   95-10-31  IREDELL     REMOVED SAVES AND PRINTS
C
C USAGE:    CALL W3FI75 (IBITL,ITYPE,ITOSS,FLD,IFLD,IBMAP,IBDSFL,
C    &              NPTS,BDS11,IPFLD,PFLD,LEN,LENBDS,IBERR,PDS,IGDS)
C   INPUT ARGUMENT LIST:
C     IBITL     - 0, COMPUTER COMPUTES PACKING LENGTH FROM POWER
C                    OF 2 THAT BEST FITS THE DATA.
C                 8, 12, ETC. COMPUTER RESCALES DATA TO FIT INTO
C                    SET NUMBER OF BITS.
C     ITYPE     - 0 = IF INPUT DATA IS FLOATING POINT (FLD)
C                 1 = IF INPUT DATA IS INTEGER (IFLD)
C     ITOSS     - 0 = NO BIT MAP IS INCLUDED (DON'T TOSS DATA)
C                 1 = TOSS NULL DATA ACCORDING TO IBMAP
C     FLD       - REAL ARRAY OF DATA TO BE PACKED IF ITYPE=0
C     IFLD      - INTEGER ARRAY TO BE PACKED IF ITYPE=1
C     IBMAP     - BIT MAP SUPPLIED FROM USER
C     IBDSFL    - INTEGER ARRAY CONTAINING TABLE 11 FLAG INFO
C                 BDS OCTET 4:
C                 (1) 0 = GRID POINT DATA
C                     1 = SPHERICAL HARMONIC COEFFICIENTS
C                 (2) 0 = SIMPLE PACKING
C                     1 = SECOND ORDER PACKING
C                 (3) 0 = ORIGINAL DATA WERE FLOATING POINT VALUES
C                     1 = ORIGINAL DATA WERE INTEGER VALUES
C                 (4) 0 = NO ADDITIONAL FLAGS AT OCTET 14
C                     1 = OCTET 14 CONTAINS FLAG BITS 5-12
C                 (5) 0 = RESERVED - ALWAYS SET TO 0
C                 (6) 0 = SINGLE DATUM AT EACH GRID POINT
C                     1 = MATRIX OF VALUES AT EACH GRID POINT
C                 (7) 0 = NO SECONDARY BIT MAPS
C                     1 = SECONDARY BIT MAPS PRESENT
C                 (8) 0 = SECOND ORDER VALUES HAVE CONSTANT WIDTH
C                     1 = SECOND ORDER VALUES HAVE DIFFERENT WIDTHS
C     NPTS      - NUMBER OF GRIDPOINTS IN ARRAY TO BE PACKED
C     IGDS      - ARRAY OF GDS INFORMATION
C
C   OUTPUT ARGUMENT LIST:
C     BDS11     - FIRST 11 OCTETS OF BDS
C     PFLD      - PACKED GRIB FIELD
C     LEN       - LENGTH OF PFLD
C     LENBDS    - LENGTH OF BDS
C     IBERR     - 1, ERROR CONVERTING IEEE F.P. NUMBER TO IBM370 F.P.
C
C REMARKS: SUBPROGRAM CAN BE CALLED FROM A MULTIPROCESSING ENVIRONMENT.
C
C ATTRIBUTES:
C   LANGUAGE: IBM VS FORTRAN 77, CRAY CFT77 FORTRAN
C   MACHINE:  HDS, CRAY C916/256, Y-MP8/64, Y-MP EL92/256
C
C$$$
C
      REAL            FLD(*)
C     REAL            FWORK(260000)
C
C     FWORK CAN USE DYNAMIC ALLOCATION OF MEMORY ON CRAY
C
      REAL            FWORK(NPTS)
      REAL            RMIN,REFNCE
C
      INTEGER         IPFLD(*)
      INTEGER         IBDSFL(*)
      INTEGER         IBMAP(*)
      INTEGER         IFLD(*),IGDS(*)
C     INTEGER         IWORK(260000)
C
C     IWORK CAN USE DYNAMIC ALLOCATION OF MEMORY ON CRAY
C
      INTEGER         IWORK(NPTS)
C
      LOGICAL         CONST
C
      CHARACTER * 1   BDS11(11),PDS(*)
      CHARACTER * 1   PFLD(*)
      CHARACTER * 1   CIEXP(8)
      CHARACTER * 1   CIMANT(8)
C
      EQUIVALENCE     (IEXP,CIEXP(1))
      EQUIVALENCE     (IMANT,CIMANT(1))
C
C            1.0   PACK THE FIELD.
C
C            1.1   TOSS DATA IF BITMAP BEING USED,
C                  MOVING 'DATA' TO WORK AREA...
C
      CONST = .FALSE.
      IBERR = 0
      IW    = 0
C
      IF (ITOSS .EQ. 1) THEN
        IF (ITYPE .EQ. 0) THEN
          DO 110 IT=1,NPTS
            IF (IBMAP(IT) .EQ. 1) THEN
              IW = IW + 1
              FWORK(IW) = FLD(IT)
            ENDIF
  110     CONTINUE
          NPTS = IW
        ELSE IF (ITYPE .EQ. 1) THEN
          DO 111 IT=1,NPTS
            IF (IBMAP(IT) .EQ. 1) THEN
              IW = IW + 1
              IWORK(IW) = IFLD(IT)
            ENDIF
  111     CONTINUE
          NPTS = IW
        ENDIF
C
C             ELSE, JUST MOVE DATA TO WORK ARRAY
C
      ELSE IF (ITOSS .EQ. 0) THEN
        IF (ITYPE .EQ. 0) THEN
          DO 112 IT=1,NPTS
            FWORK(IT) = FLD(IT)
  112     CONTINUE
        ELSE IF (ITYPE .EQ. 1) THEN
          DO 113 IT=1,NPTS
            IWORK(IT) = IFLD(IT)
  113     CONTINUE
        ENDIF
      ENDIF
C
C            1.2   CONVERT DATA IF NEEDED PRIOR TO PACKING.
C                  (INTEGER TO F.P. OR F.P. TO INTEGER)
C     ITYPE = 0...FLOATING POINT DATA
C       IBITL = 0...PACK IN LEAST # BITS...CONVERT TO INTEGER
C     ITYPE = 1...INTEGER DATA
C       IBITL > 0...PACK IN FIXED # BITS...CONVERT TO FLOATING POINT
C
      IF (ITYPE .EQ. 0 .AND. IBITL .EQ. 0) THEN
        DO 120 IF=1,NPTS
          IWORK(IF) = NINT(FWORK(IF))
  120   CONTINUE
      ELSE IF (ITYPE .EQ. 1 .AND. IBITL .NE. 0) THEN
        DO 123 IF=1,NPTS
          FWORK(IF) = FLOAT(IWORK(IF))
  123   CONTINUE
      ENDIF
C
C            1.3   PACK THE DATA.
C
      IF (IBDSFL(2).NE.0) THEN
C                                    SECOND ORDER PACKING
C
C            PRINT*,'  DOING SECOND ORDER PACKING...'
          IF (IBITL.EQ.0) THEN
C
C             PRINT*,'    AND VARIABLE BIT PACKING'
C
C                           WORKING WITH INTEGER VALUES
C                           SINCE DOING VARIABLE BIT PACKING
C
              MAX  = IWORK(1)
              MIN  = IWORK(1)
              DO 300 I = 2, NPTS
                  IF (IWORK(I).LT.MIN) THEN
                      MIN  = IWORK(I)
                  ELSE IF (IWORK(I).GT.MAX) THEN
                      MAX  = IWORK(I)
                  END IF
  300         CONTINUE
C                           EXTRACT MINIMA
              DO 400 I = 1, NPTS
C                 IF (IWORK(I).LT.0) THEN
C                     PRINT *,'MINIMA 400',I,IWORK(I),NPTS
C                 END IF
                  IWORK(I)  = IWORK(I) - MIN
  400         CONTINUE
              REFNCE  = MIN
              IDIFF   = MAX - MIN
C             PRINT *,'REFERENCE VALUE',REFNCE
C
C             WRITE (6,FMT='(''  MINIMA REMOVED      = '',/,
C    &              10(3X,10I10,/))') (IWORK(I),I=1,6)
C             WRITE (6,FMT='(''  END OF ARRAY  = '',/,
C    &              10(3X,10I10,/))') (IWORK(I),I=NPTS-5,NPTS)
C
C                      FIND BIT WIDTH OF IDIFF
C
              CALL FI7505 (IDIFF,KWIDE)
C             PRINT*,'  BIT WIDTH FOR ORIGINAL DATA', KWIDE
              ISCAL2 = 0
C
C             MULTIPLICATIVE SCALE FACTOR SET TO 1
C             IN ANTICIPATION OF POSSIBLE USE IN GLAHN 2DN DIFF
C
              SCAL2 = 1.
C
          ELSE
C
C             PRINT*,'   AND FIXED BIT PACKING, IBITL = ', IBITL
C                               FIXED BIT PACKING
C                               - LENGTH OF FIELD IN IBITL
C                               - MUST BE REAL DATA
C                            FLOATING POINT INPUT
C
              RMAX  = FWORK(1)
              RMIN  = FWORK(1)
              DO 100 I = 2, NPTS
                  IF (FWORK(I).LT.RMIN) THEN
                      RMIN  = FWORK(I)
                  ELSE IF (FWORK(I).GT.RMAX) THEN
                      RMAX  = FWORK(I)
                  END IF
  100         CONTINUE
              REFNCE  = RMIN
C             PRINT *,'100 REFERENCE',REFNCE
C                             EXTRACT MINIMA
              DO 200 I = 1, NPTS
                  FWORK(I)  = FWORK(I) - RMIN
  200         CONTINUE
C             PRINT *,'REFERENCE VALUE',REFNCE
C             WRITE (6,FMT='(''  MINIMA REMOVED      = '',/,
C    &              10(3X,10F8.2,/))') (FWORK(I),I=1,6)
C             WRITE (6,FMT='(''  END OF ARRAY  = '',/,
C    &              10(3X,10F8.2,/))') (FWORK(I),I=NPTS-5,NPTS)
C                                FIND LARGEST DELTA
              IDELT  = NINT(RMAX - RMIN)
C                                DO BINARY SCALING
C                                   FIND OUT WHAT BINARY SCALE FACTOR
C                                       PERMITS CONTAINMENT OF
C                                       LARGEST DELTA
              CALL FI7505 (IDELT,IWIDE)
C
C                                   BINARY SCALING
C
              ISCAL2  = IWIDE - IBITL
C             PRINT *,'SCALING NEEDED TO FIT =',ISCAL2
C             PRINT*,'  RANGE OF  = ',IDELT
C
C                                EXPAND DATA WITH BINARY SCALING
C                                CONVERT TO INTEGER
              SCAL2  = 2.0**ISCAL2
              SCAL2  = 1./ SCAL2
              DO 600 I = 1, NPTS
                  IWORK(I)  = NINT(FWORK(I) * SCAL2)
  600         CONTINUE
              KWIDE = IBITL
          END IF
C
C  *****************************************************************
C
C           FOLLOWING IS FOR GLAHN SECOND DIFFERENCING
C           NOT STANDARD GRIB
C
C            TEST FOR SECOND DIFFERENCE PACKING
C            BASED OF SIZE OF PDS - SIZE IN FIRST 3 BYTES
C
          CALL GBYTE (PDS,IPDSIZ,0,24)
          IF (IPDSIZ.EQ.50) THEN
C             PRINT*,'  DO SECOND DIFFERENCE PACKING '
C
C                   GLAHN PACKING TO 2ND DIFFS
C
C             WRITE (6,FMT='(''  CALL TO W3FI82 WITH = '',/,
C    &                  10(3X,10I6,/))') (IWORK(I),I=1,NPTS)
C
               CALL W3FI82 (IWORK,FVAL1,FDIFF1,NPTS,PDS,IGDS)
C
C             PRINT *,'GLAHN',FVAL1,FDIFF1
C             WRITE (6,FMT='(''  OUT FROM W3FI82 WITH = '',/,
C    &                  10(3X,10I6,/))') (IWORK(I),I=1,NPTS)
C
C             MUST NOW RE-REMOVE THE MINIMUM VALUE
C             OF THE SECOND DIFFERENCES TO ASSURE
C             ALL POSITIVE NUMBERS FOR SECOND ORDER GRIB PACKING
C
C             ORIGINAL REFERENCE VALUE ADDED TO FIRST POINT
C             VALUE FROM THE 2ND DIFF PACKER TO BE ADDED
C             BACK IN WHEN THE 2ND DIFF VALUES ARE
C             RECONSTRUCTED BACK TO THE BASIC VALUES
C
C             ALSO, THE REFERENCE VALUE IS
C             POWER-OF-TWO SCALED TO MATCH
C             FVAL1.  ALL OF THIS SCALING
C             WILL BE REMOVED AFTER THE
C             GLAHN SECOND DIFFERENCING IS UNDONE.
C             THE SCALING FACTOR NEEDED TO DO THAT
C             IS SAVED IN THE PDS AS A SIGNED POSITIVE
C             TWO BYTE INTEGER
C
C             THE SCALING FOR THE 2ND DIF PACKED
C             VALUES IS PROPERLY SET TO ZERO
C
              FVAL1 = FVAL1 + REFNCE*SCAL2
C                                          FIRST TEST TO SEE IF
C                                          ON 32 OR 64 BIT COMPUTER
              CALL W3FI01(LW)
              IF (LW.EQ.4) THEN
                  CALL W3FI76 (FVAL1,IEXP,IMANT,32)
              ELSE
                  CALL W3FI76 (FVAL1,IEXP,IMANT,64)
              END IF
              CALL SBYTE (PDS,IEXP,320,8)
              CALL SBYTE (PDS,IMANT,328,24)
C
              IF (LW.EQ.4) THEN
                  CALL W3FI76 (FDIFF1,IEXP,IMANT,32)
              ELSE
                  CALL W3FI76 (FDIFF1,IEXP,IMANT,64)
              END IF
              CALL SBYTE (PDS,IEXP,352,8)
              CALL SBYTE (PDS,IMANT,360,24)
C
C             TURN ISCAL2 INTO SIGNED POSITIVE INTEGER
C             AND STORE IN TWO BYTES
C
              IF(ISCAL2.GE.0)  THEN
                CALL SBYTE (PDS,ISCAL2,384,16)
              ELSE
                CALL SBYTE (PDS,1,384,1)
                ISCAL2 = - ISCAL2
                CALL SBYTE( PDS,ISCAL2,385,15)
              ENDIF
C
              MAX  = IWORK(1)
              MIN  = IWORK(1)
              DO 700 I = 2, NPTS
                  IF (IWORK(I).LT.MIN) THEN
                      MIN  = IWORK(I)
                  ELSE IF (IWORK(I).GT.MAX) THEN
                      MAX  = IWORK(I)
                  END IF
  700         CONTINUE
C                           EXTRACT MINIMA
              DO 710 I = 1, NPTS
                  IWORK(I)  = IWORK(I) - MIN
  710         CONTINUE
              REFNCE  = MIN
C             PRINT *,'710 REFERENCE',REFNCE
              ISCAL2 = 0
C
C             AND RESET VALUE OF KWIDE - THE BIT WIDTH
C             FOR THE RANGE OF THE VALUES
C
              IDIFF = MAX - MIN
              CALL FI7505 (IDIFF,KWIDE)
C
C             PRINT*,'BIT WIDTH (KWIDE) OF 2ND DIFFS', KWIDE
C
C  **************************** END OF GLAHN PACKING  ************
          ELSE IF (IBDSFL(2).EQ.1.AND.IBDSFL(7).EQ.0) THEN
C                        HAVE SECOND ORDER PACKING WITH NO SECOND ORDER
C                        BIT MAP. ERGO ROW BY ROW - COL BY COL
              CALL FI7503 (IWORK,IPFLD,NPTS,IBDSFL,BDS11,
     *              LEN,LENBDS,PDS,REFNCE,ISCAL2,KWIDE,IGDS)
              RETURN
          END IF
C         WRITE (6,FMT='(''  CALL TO FI7501 WITH = '',/,
C    &                  10(3X,10I6,/))') (IWORK(I),I=1,NPTS)
C         WRITE (6,FMT='(''  END OF ARRAY = '',/,
C    &                  10(3X,10I6,/))') (IWORK(I),I=NPTS-5,NPTS)
C         PRINT*,' REFNCE,ISCAL2, KWIDE AT CALL TO FI7501',
C    &             REFNCE, ISCAL2,KWIDE
C
C                         SECOND ORDER PACKING
C
          CALL FI7501 (IWORK,IPFLD,NPTS,IBDSFL,BDS11,
     *             LEN,LENBDS,PDS,REFNCE,ISCAL2,KWIDE)
C
C              BDS COMPLETELY ASSEMBLED IN FI7501 FOR SECOND ORDER
C              PACKING.
C
      ELSE
C                                      SIMPLE PACKING
C
C                PRINT*,'  SIMPLE FIRST ORDER PACKING...'
          IF (IBITL.EQ.0) THEN
C                PRINT*,' WITH VARIABLE BIT LENGTH'
C
C                  WITH VARIABLE BIT LENGTH, ADJUSTED
C                  TO ACCOMMODATE LARGEST VALUE
C                  BINARY SCALING ALWAYS = 0
C
              CALL W3FI58(IWORK,NPTS,IWORK,PFLD,NBITS,LEN,KMIN)
              RMIN   = KMIN
              REFNCE  = RMIN
              ISCALE = 0
C             PRINT*,'  BIT LENGTH CAME OUT AT ...',NBITS
C
C           SET CONST .TRUE. IF ALL VALUES ARE THE SAME
C
              IF (LEN.EQ.0.AND.NBITS.EQ.0) CONST = .TRUE.
C
          ELSE
C           PRINT*,' FIXED BIT LENGTH, IBITL = ', IBITL
C
C             FIXED BIT LENGTH PACKING (VARIABLE PRECISION)
C             VALUES SCALED BY POWER OF 2 (ISCALE) TO
C             FIT LARGEST VALUE INTO GIVEN BIT LENGTH (IBITL)
C
              CALL W3FI59(FWORK,NPTS,IBITL,IWORK,PFLD,ISCALE,LEN,RMIN)
              REFNCE = RMIN
C             PRINT *,' SCALING NEEDED TO FIT IS ...', ISCALE
              NBITS = IBITL
C
C           SET CONST .TRUE. IF ALL VALUES ARE THE SAME
C
              IF (LEN.EQ.0) THEN
                  CONST = .TRUE.
                  NBITS = 0
              END IF
          END IF
C
C$        COMPUTE LENGTH OF BDS IN OCTETS
C
          INUM  = NPTS * NBITS + 88
C         PRINT *,'NUMBER OF BITS BEFORE FILL ADDED',INUM
C
C                  NUMBER OF FILL BITS
          NFILL  = 0
          NLEFT  = MOD(INUM,16)
          IF (NLEFT.NE.0) THEN
              INUM  = INUM + 16 - NLEFT
              NFILL = 16 - NLEFT
          END IF
C         PRINT *,'NUMBER OF BITS AFTER FILL ADDED',INUM
C                  LENGTH OF BDS IN BYTES
          LENBDS = INUM / 8
C
C                2.0   FORM THE BINARY DATA SECTION (BDS).
C
C                 CONCANTENATE ALL FIELDS FOR BDS
C
C                               BYTES 1-3
          CALL SBYTE (BDS11,LENBDS,0,24)
C
C                               BYTE  4
C                                       FLAGS
          CALL SBYTE (BDS11,IBDSFL(1),24,1)
          CALL SBYTE (BDS11,IBDSFL(2),25,1)
          CALL SBYTE (BDS11,IBDSFL(3),26,1)
          CALL SBYTE (BDS11,IBDSFL(4),27,1)
C                                        NR OF FILL BITS
          CALL SBYTE (BDS11,NFILL,28,4)
C
C$      FILL OCTETS 5-6 WITH THE SCALE FACTOR.
C
C                               BYTE  5-6
          IF (ISCALE.LT.0) THEN
              CALL SBYTE (BDS11,1,32,1)
              ISCALE  = - ISCALE
              CALL SBYTE (BDS11,ISCALE,33,15)
          ELSE
              CALL SBYTE (BDS11,ISCALE,32,16)
          END IF
C
C$  FILL OCTET 7-10 WITH THE REFERENCE VALUE
C   CONVERT THE FLOATING POINT OF YOUR MACHINE TO IBM370 32 BIT
C   FLOATING POINT NUMBER
C
C                               BYTE  7-10
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
C
C$                        FILL OCTET 11 WITH THE NUMBER OF BITS.
C
C                               BYTE  11
          CALL SBYTE (BDS11,NBITS,80,8)
      END IF
C
      RETURN
      END
