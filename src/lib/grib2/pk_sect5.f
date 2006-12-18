      SUBROUTINE PK_SECT5(KFILDO,IPACK,ND5,IS5,NS5,MINA,XMINA,
     1                    MISSP,XMISSP,MISSS,XMISSS,L3264B,
     2                    LOCN,IPOS,LOCN5_20,IPOS5_20,LOCN5_32,
     3                    IPOS5_32,IER,ISEVERE,*)
C
C        MARCH    2000   GLAHN   TDL   FOR GRIB2
C        JANUARY  2001   GLAHN   COMMENTS; ELIMINATED TVALUE; ADDED
C                                CHECK ON SIZE OF IS5( )
C        JANUARY  2001   GLAHN/LAWRENCE RESTRUCTURED TO CHECK NS5
C                                BASED ON PACKING TYPE
C        FEBRUARY 2001   GLAHN/LAWRENCE CHANGED NS5.LT.10 TO LT.21
C                                FOR CASE (0); ELIMINATED PACKING
C                                REFERENCE VALUE AS INTEGER; COMMENTS
C        NOVEMBER 2001   GLAHN   CHANGED TEST ON NS5 FROM 10 TO 21 
C                                FOR SIMPLE PACKING.
C        DECEMBER 2001   GLAHN   MOVED TEST ON IS5(10) TO PREPR
C
C        PURPOSE
C            PACKS SECTION 5, THE DATA REPRESENTATION SECTION,
C            OF A GRIB2 MESSAGE.  THIS SUPPORTS ONLY SIMPLE,
C            COMPLEX, AND COMPLEX WITH SPATIAL DIFFERENCING.
C
C        DATA SET USE
C           KFILDO - UNIT NUMBER FOR OUTPUT (PRINT) FILE. (OUTPUT)
C
C        VARIABLES
C              KFILDO = UNIT NUMBER FOR OUTPUT (PRINT) FILE. (INPUT)
C            IPACK(J) = THE ARRAY THAT HOLDS THE ACTUAL PACKED MESSAGE
C                       (J=1,ND5). (INPUT/OUTPUT)
C                 ND5 = THE SIZE OF THE ARRAY IPACK( ). (INPUT)
C              IS5(J) = CONTAINS THE DATA REPRESENTATION DATA THAT
C                       WILL BE PACKED INTO IPACK( ) (J=1,NS5).
C                       (INPUT/OUTPUT)
C                 NS5 = SIZE OF IS5( ). (INPUT)
C                MINA = REFERENCE VALUE WHEN INTEGER.  (INPUT)
C               XMINA = REFERENCE VALUE WHEN FLOATING POINT.  (INPUT)
C               MISSP = PRIMARY MISSING VALUE, WHEN INTEGER.  (INPUT)
C              XMISSP = PRIMARY MISSING VALUE, WHEN FLOATING POINT.
C                       (INPUT)
C               MISSS = SECONDARY MISSING VALUE, WHEN INTEGER.  (INPUT)
C              XMISSS = SECONDARY MISSING VALUE, WHEN FLOATING POINT.
C                       (INPUT)
C              L3264B = THE INTEGER WORD LENGTH IN BITS OF THE MACHINE
C                       BEING USED. VALUES OF 32 AND 64 ARE
C                       ACCOMMODATED.  (INPUT)
C                LOCN = THE WORD POSITION TO PLACE THE NEXT VALUE.
C                       (INPUT/OUTPUT)
C                IPOS = THE BIT POSITION IN LOCN TO START PLACING
C                       THE NEXT VALUE. (INPUT/OUTPUT)
C           LOCN5_20 = LOCN FOR OCTET 20 IN SECTION 5.  (OUTPUT)
C           IPOS5_20 = IPOS FOR OCTET 20 IN SECTION 5.  (OUTPUT)
C           LOCN5_32 = LOCN FOR OCTET 32 IN SECTION 5.  (OUTPUT)
C           IPOS5_32 = IPOS FOR OCTET 32 IN SECTION 5.  (OUTPUT)
C                IER = RETURN STATUS CODE. (OUTPUT)
C                        0 = GOOD RETURN.
C                      1-4 = ERROR CODES GENERATED BY PKBG. SEE THE
C                            DOCUMENTATION IN THE PKBG ROUTINE.
C                      501 = IS5(5) DOES NOT INDICATE SECTION 5.
C                      502 = IS5( ) HAS NOT BEEN DIMENSIONED LARGE
C                            ENOUGH TO CONTINUE SECTION 5.
C             ISEVERE = THE SEVERITY LEVEL OF THE ERROR.  THE ONLY
C                       VALUE RETUNED IS:
C                       2 = A FATAL ERROR  (OUTPUT)
C                  * = ALTERNATE RETURN WHEN IER NE 0 FROM
C                      SUBROUTINES.
C
C             LOCAL VARIABLES
C              IVALUE = THE INTEGER EQUIVALENCE OF THE 'IEEE' REAL
C                       NUMBER CONTAINED IN RVALUE.
C               IZERO = CONTAINS THE VALUE '0'.  THIS IS USED IN THE
C                       CODE SIMPLY TO EMPHASIZE THAT A CERTAIN 
C                       GROUP OF OCTETS IN THE MESSAGE ARE BEING 
C                       ZEROED OUT.
C                   N = L3264B = THE INTEGER WORD LENGTH IN BITS OF
C                       THE MACHINE BEING USED. VALUES OF 32 AND
C                       64 ARE ACCOMMODATED.
C              RVALUE = CONTAINS A REAL NUMBER IN IEEE FORMAT.
C               ISIGN = SIGN OF VALUE BEING PACKED, 0 = POSITIVE,
C                       1 = NEGATIVE.  THE SIGN ALWAYS GOES IN THE
C                       LEFTMOST BIT OF THE AREA ASSIGNED TO THAT VALUE.
C
C        NON SYSTEM SUBROUTINES CALLED
C           FMKIEEE, PKBG
C
      DIMENSION IPACK(ND5),IS5(NS5)
C
      DATA IZERO/0/
C
      EQUIVALENCE(RVALUE,IVALUE)
C
      IER=0
      N=L3264B
C
C        ALL ERRORS GENERATED BY THIS ROUTINE ARE FATAL.
C 
      ISEVERE=2
C
C        CHECK TO MAKE SURE THAT IS5(5) SPECIFIES SECTION 5.
C
      IF(IS5(5).NE.5)THEN
         IER=501
      ELSE
C
C           CHECK MINIMUM SIZE OF IS5( ), WHICH DEPENDS ON
C           PACKING METHOD.
C
         SELECT CASE (IS5(10))
C
            CASE (0)
C                 THIS CASE IS FOR SIMPLE PACKING.
C
               IF(NS5.LT.21)THEN
                  IER=502
                  GO TO 900
               ENDIF
C
            CASE (2)
C                 THIS CASE IS FOR COMPLEX PACKING.
C
               IF(NS5.LT.47)THEN
                  IER=502
                  GO TO 900
               ENDIF
C
            CASE (3)
C                 THIS CASE IS FOR COMPLEX PACKING WITH SECOND
C                 ORDER DIFFERENCES.
C
               IF(NS5.LT.49)THEN
                  IER=502
                  GO TO 900
               ENDIF
C         
         END SELECT
C
C           PACK SECTION LENGTH.  IS5(1) IS INITIALIZED BY PREPR.
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IS5(1),32,N,IER,*900)
C
C           PACK THE SECTION NUMBER.  IS5(5) IS INPUT BY THE USER.
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IS5(5),8,N,IER,*900)
C
C           PACK THE NUMBER OF DATA POINTS.  IS5(6) IS INITIALIZED
C           BY PREPR.
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IS5(6),32,N,IER,*900)
C
C           PACK THE DATA REPRESENTATION TEMPLATE NUMBER.  IS5(10)
C           IS PROVIDED BY THE USER, BUT MAY BE OVERRIDDEN BY PREPR.
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IS5(10),16,N,IER,*900)
C
C           THE FOLLOWING SECTION IS FOR SIMPLE, COMPLEX, AND COMPLEX
C           WITH SECOND ORDER DIFFERENCE PACKING.
C
C           PACK THE REFERENCE VALUE.  IT MAY BE INTEGER OR
C           FLOATING POINT, BUT IS PACKED AS FLOATING POINT.
C           IS5(12) IS THE SCALED VALUE, AND MAY NOT BE OF ANY
C           USE TO THE USER, PARTICULARLY WHEN THE DATA WERE
C           FLOATING POINT.
C
         IF(IS5(21).EQ.0)THEN
C              THIS IS FLOATING POINT.
            IS5(12)=NINT(XMINA)
         ELSE
C              THIS IS INTEGER.
            XMINA=MINA
            IS5(12)=MINA
         ENDIF
C
         RVALUE=FMKIEEE(XMINA)
C           RVALUE AND IVALUE ARE EQUIVALENCED.
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IVALUE,32,N,IER,*900)
C
C            PACK THE BINARY AND DECIMAL SCALE FACTORS TAKING INTO
C            ACCOUNT THAT THESE VALUES MAY BE NEGATIVE.
C
         ISIGN=0
         IF(IS5(16).LT.0) ISIGN=1
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,ISIGN,1,N,IER,*900)
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,ABS(IS5(16)),15,N,
     1             IER,*900)
         ISIGN=0
         IF(IS5(18).LT.0) ISIGN=1
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,ISIGN,1,N,IER,*900)
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,ABS(IS5(18)),15,N,
     1             IER,*900)
C         
C           SAVE LOCATION FOR FIELD WIDTH FOR THE DATA.
C
         LOCN5_20=LOCN
         IPOS5_20=IPOS
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IZERO,8,N,IER,*900)
C       
C           PACK THE TYPE OF THE ORIGINAL FIELD VALUES.
C
         CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IS5(21),8,N,IER,*900)
C
C           THE FOLLOWING SECTION IS FOR ONLY COMPLEX, AND COMPLEX
C           WITH SECOND ORDER DIFFERENCE PACKING.
C
         IF(IS5(10).NE.0)THEN
C
C              THE PACKING METHOD IS COMPLEX OR SPATIAL DIFFERENCING.
C              PACK THE SPLITTING METHOD.
C
            CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IS5(22),8,N,IER,*900)
C
C              PACK THE PROVISION FOR MISSING VALUES.
C
            CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IS5(23),8,N,IER,*900)
C
C              PACK THE PRIMARY AND SECONDARY MISSING VALUES.  THEY
C              MAY BE PACKED EITHER INTEGER OR FLOATING POINT CONSISTENT
C              WITH THE ORIGINAL DATA (TEMPLATE 5.2, NOTE 7).
C
            IF(IS5(21).EQ.0)THEN
C
C                 THE MISSING VALUES ARE FLOATING POINT.
C
               IS5(24)=NINT(XMISSP)
               RVALUE=FMKIEEE(XMISSP)
               CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IVALUE,32,N,
     1                   IER,*900)
C
               IS5(28)=NINT(XMISSS)
               RVALUE=FMKIEEE(XMISSS)
               CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IVALUE,32,N,
     1                   IER,*900)
C
            ELSE
C
C                 THE MISSING VALUES ARE INTEGER.
C
               IS5(24)=MISSP
               ISIGN=0
               IF(MISSP.LT.0)ISIGN=1
               CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,ISIGN,1,N,IER,*900)
               CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,ABS(MISSP),31,N,
     1                   IER,*900)
               IS5(28)=MISSS
               ISIGN=0
               IF(MISSS.LT.0)ISIGN=1
               CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,ISIGN,1,N,IER,*900)
               CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,ABS(MISSS),31,N,
     1                   IER,*900)
            ENDIF
C         
C              SAVE THE LOCATIONS FOR THE REMAINDER OF SECTION 5.
C
            LOCN5_32=LOCN
            IPOS5_32=IPOS
C
C              SAVE SPACE FOR THE REMAINDER OF SECTION 5.  THESE WILL
C              BE FILLED BY PK_SECT7.
C                                  
            CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IZERO,32,N,IER,*900)
            CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IZERO,8,N,IER,*900)
            CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IZERO,8,N,IER,*900)
            CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IZERO,32,N,IER,*900)
            CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IZERO,8,N,IER,*900)
            CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IZERO,32,N,IER,*900)
            CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IZERO,8,N,IER,*900)
C
C              THE FOLLOWING SECTION IS FOR ONLY COMPLEX WITH SECOND
C              ORDER DIFFERENCE PACKING.
C
            IF(IS5(10).EQ.3)THEN
               CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IZERO,8,N,
     1                   IER,*900)
               CALL PKBG(KFILDO,IPACK,ND5,LOCN,IPOS,IZERO,8,N,
     1                   IER,*900)
C
            ENDIF
C
         ENDIF           
C
      ENDIF
C
C        ERROR RETURN SECTION
 900  IF(IER.NE.0)RETURN 1
C
      RETURN
      END
