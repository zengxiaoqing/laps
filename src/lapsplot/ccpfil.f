      subroutine ccpfil(ZREG,MREG,NREG)

C 
C Define error file, Fortran unit number, and workstation type,
C and workstation ID.
C 
      PARAMETER (IERRF=6, LUNIT=2, IWTYPE=1, IWKID=1)
      REAL XREG(MREG),YREG(NREG),ZREG(MREG,NREG)
      
      EXTERNAL COLOR2
C      
C Get data array
C
!     CALL GETDAT(XREG,YREG,ZREG,MREG,NREG)
C 
C Open GKS, open and activate a workstation.
C 
!     CALL GOPKS (IERRF, ISZDM)
!     CALL GOPWK (IWKID, LUNIT, IWTYPE)
!     CALL GACWK (IWKID)
C      
C Call Conpack color fill routine
C      
      CALL CCPFIL_SUB(ZREG,MREG,NREG,-15,COLOR2,IWKID)
C      
C Close frame
C      
!     CALL FRAME
C 
C Deactivate and close workstation, close GKS.
C 
!     CALL GDAWK (IWKID)
!     CALL GCLWK (IWKID)
!     CALL GCLKS

      call color
      
      return
      END

      
      SUBROUTINE CCPFIL_SUB(ZREG,MREG,NREG,NCL,COLOR2,IWKID)
      
      PARAMETER (LRWK=150000,LIWK=150000,LMAP=1000000,NWRK=150000
     1          ,NOGRPS=5)       
      REAL ZREG(MREG,NREG),RWRK(LRWK), XWRK(NWRK), YWRK(NWRK)
      INTEGER MREG,NREG,IWRK(LIWK)
      INTEGER MAP(LMAP),IAREA(NOGRPS),IGRP(NOGRPS)
      
      EXTERNAL FILL
      EXTERNAL COLOR2
C      
C Set up color table
C      
      CALL COLOR2(IWKID)
C      
C Initialize Areas
C      
      CALL ARINAM(MAP, LMAP)
C      
C Set number of contour levels and initialize Conpack
C      
!      CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',NCL)

      CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',+1)
      CALL CPSETR('CIS',0.1)
      CALL CPSETR('CMN',0.0)
      CALL CPSETR('CMX',1.0)

      CALL CPRECT(ZREG, MREG, MREG, NREG, RWRK, LRWK, IWRK, LIWK)
C      
C Add contours to area map
C      
      CALL CPCLAM(ZREG, RWRK, IWRK, MAP)
C      
C Set fill style to solid, and fill contours
C      
      CALL GSFAIS(1)
      CALL ARSCAM(MAP, XWRK, YWRK, NWRK, IAREA, IGRP, NOGRPS, FILL)
C      
C Draw Perimeter
C      
!     CALL CPBACK(ZREG, RWRK, IWRK)
C      
C Draw Labels
C      
!     CALL CPLBDR(ZREG,RWRK,IWRK)
C      
C Draw Contours
C      
!     CALL CPCLDR(ZREG,RWRK,IWRK)
      
      RETURN
      END
      
      SUBROUTINE FILL (XWRK,YWRK,NWRK,IAREA,IGRP,NGRPS)
C 
      DIMENSION XWRK(*),YWRK(*),IAREA(*),IGRP(*)
      
      DO 10, I=1,NGRPS
         IF (IGRP(I).EQ.3) IAREA3=IAREA(I)
 10   CONTINUE
      
      IF (IAREA3 .GT. 0) THEN
C      
C If the area is defined by 3 or more points, fill it
C      
         CALL GSFACI(IAREA3+1)
         CALL GFA(NWRK,XWRK,YWRK)
      ENDIF
C      
C Otherwise, do nothing
C      
      RETURN
      END


      
      SUBROUTINE COLOR2(IWKID)
C 
C BACKGROUND COLOR
C BLACK
C
      CALL GSCR(IWKID,0,0.,0.,0.)

      do i = 1,255
!         rintens = min(max(float(i-2) / 9.,0.),1.)
          rintens = min(max(float(i-2) / 9.,0.),1.)
          call GSCR(IWKID, i, rintens, rintens, rintens)
      enddo ! i

C 
      RETURN
C 
      END

