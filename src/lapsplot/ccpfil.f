      subroutine ccpfil(field_in,MREG,NREG,scale_l_in,scale_h_in
     1                 ,colortable)       

C 
C Define error file, Fortran unit number, and workstation type,
C and workstation ID.
C 
      PARAMETER (IERRF=6, LUNIT=2, IWTYPE=1, IWKID=1)
      REAL XREG(MREG),YREG(NREG),ZREG(MREG,NREG),field_in(MREG,NREG)
      character*(*)colortable
      
      write(6,*)' Subroutine ccpfil for solid fill plot...'

      if(scale_l_in .lt. scale_h_in)then
          ireverse = 0
          scale_l = scale_l_in
          scale_h = scale_h_in
      else
          ireverse = 1
          scale_l = scale_h_in
          scale_h = scale_l_in
      endif

!     Apply scaling to the array
      scale = scale_h - scale_l
      call addcon(field_in,-scale_l,ZREG,MREG,NREG)

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
      write(6,*)' Colortable is ',colortable,scale_l,scale_h,ireverse
      CALL CCPFIL_SUB(ZREG,MREG,NREG,-15,IWKID,scale,ireverse
     1                                              ,colortable)      
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

      jdot = 1
      
      return
      END

      
      SUBROUTINE CCPFIL_SUB(ZREG,MREG,NREG,NCL,IWKID,scale,ireverse
     1                                                    ,colortable)      
      
      PARAMETER (LRWK=150000,LIWK=150000,LMAP=1000000,NWRK=150000
     1          ,NOGRPS=5)       
      REAL ZREG(MREG,NREG),RWRK(LRWK), XWRK(NWRK), YWRK(NWRK)
      INTEGER MREG,NREG,IWRK(LIWK)
      INTEGER MAP(LMAP),IAREA(NOGRPS),IGRP(NOGRPS)
      character*(*) colortable
      
      EXTERNAL FILL

      ncols = 20
C      
C Set up color table
      write(6,*)' ccpfil_sub - scale = ',scale
C      
      CALL set_image_colortable(IWKID,ncols,ireverse,colortable)
C      
C Initialize Areas
C      
      CALL ARINAM(MAP, LMAP)
C      
C Set number of contour levels and initialize Conpack
C      
!      CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',NCL)

      CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',+1)
      CALL CPSETR('CIS', abs(scale) / float(ncols))
      CALL CPSETR('CMN',0.0 * abs(scale))
      CALL CPSETR('CMX',1.0 * abs(scale))

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


      
      SUBROUTINE set_image_colortable(IWKID,ncols,ireverse,colortable)    

      character*(*) colortable
C 
C BACKGROUND COLOR
C BLACK

      rcols = ncols - 1
C
      CALL GSCR(IWKID,0,0.,0.,0.)

      if(colortable .eq. 'linear')then

          do i = 1,255
!             rintens = min(max(float(i-2) / rcols,0.),1.)
              rintens = min(max(float(i-2) / rcols,0.),1.)
              if(ireverse .eq. 1)rintens = 1.0 - rintens
              call GSCR(IWKID, i, rintens, rintens, rintens)
          enddo ! i

      elseif(colortable .eq. 'hues')then
          call color_ramp(1,20,IWKID,1.0,0.7,0.7,1.0,0.7,0.7)          

      else
          write(6,*)' ERROR: Unknown color table ',colortable

      endif
C 
      RETURN
      END

      subroutine hsi_to_rgb(hue,sat,rintens,red,grn,blu)

!     Hue is 0:R, 1:B, 2:G, 3:R
      hue1 = hue + 3.0

      red = max(abs(hue1 - 3.0),0.0)
      grn = max(abs(hue  - 2.0),0.0)
      blu = max(abs(hue  - 1.0),0.0)

      red = (red*sat) + 1.0*(1.0-sat)
      grn = (grn*sat) + 1.0*(1.0-sat)
      blu = (blu*sat) + 1.0*(1.0-sat)

      red = red * rintens
      grn = grn * rintens
      blu = blu * rintens

      return
      end

      subroutine color_ramp(ncol1,ncol2,IWKID       ! I
     1                     ,hue1,sat1,rintens1      ! I
     1                     ,hue2,sat2,rintens2)     ! I

      do icol = ncol1,ncol2
          frac = float(icol-ncol1) / float(ncol2-ncol1)

          hue     = (1.0 - frac) * hue1     + frac * hue2  
          sat     = (1.0 - frac) * sat1     + frac * sat2
          rintens = (1.0 - frac) * rintens1 + frac * rintens2

          call hsi_to_rgb(hue,sat,rintens,red,grn,blu)
          
          call GSCR(IWKID,icol,red,grn,blu)
      enddo

      return
      end
