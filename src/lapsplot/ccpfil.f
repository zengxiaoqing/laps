cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
      subroutine ccpfil(field_in,MREG,NREG,scale_l_in,scale_h_in
     1                 ,colortable,n_image)       

C 
C Define error file, Fortran unit number, and workstation type,
C and workstation ID.
C 
      PARAMETER (IERRF=6, LUNIT=2, IWTYPE=1, IWKID=1)
      REAL XREG(MREG),YREG(NREG),ZREG(MREG,NREG),field_in(MREG,NREG)
      character*(*)colortable
      
      write(6,*)' Subroutine ccpfil for solid fill plot...'

      n_image = n_image + 1

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

      elseif(colortable .eq. 'hues' .or. colortable .eq. 'ref')then       
          ncols = 50
          call color_ramp(1,ncols/8,IWKID
     1                   ,0.5,0.15,0.6                 ! Pink
     1                   ,0.5,0.5,0.7)                ! Violet
          call color_ramp(ncols/8,59*ncols/120,IWKID
     1                   ,0.5,0.5,0.7                 ! Violet
     1                   ,1.5,1.0,0.7)                ! Aqua
          call color_ramp(59*ncols/120,73*ncols/120,IWKID
     1                   ,1.5,1.0,0.7                 ! Aqua
     1                   ,2.0,0.4,0.4)                ! Green
          call color_ramp(73*ncols/120,90*ncols/100,IWKID
     1                   ,2.0,0.4,0.4                 ! Green
     1                   ,3.0,0.9,0.7)                ! Red
          call color_ramp(90*ncols/100,ncols,IWKID
     1                   ,3.0,0.9,0.7                 ! Red
     1                   ,3.0,0.9,0.2)                ! Hot

          if(colortable .eq. 'ref')then
              do i = 1,3
                  call GSCR(IWKID, i, 0., 0., 0.)
              enddo 
          endif

      else
          write(6,*)' ERROR: Unknown color table ',colortable

      endif
C 
      RETURN
      END

      subroutine color_ramp(ncol1,ncol2,IWKID       ! I
     1                     ,hue1,sat1,rintens1      ! I
     1                     ,hue2,sat2,rintens2)     ! I

      write(6,*)' Subroutine color_ramp...'

      do icol = ncol1,ncol2
          frac = float(icol-ncol1) / float(ncol2-ncol1)

          hue     = (1.0 - frac) * hue1     + frac * hue2  
          sat     = (1.0 - frac) * sat1     + frac * sat2
          rintens = (1.0 - frac) * rintens1 + frac * rintens2

          call hsi_to_rgb(hue,sat,rintens,red,grn,blu)

          write(6,1)icol,hue,sat,rintens,red,grn,blu
 1        format(i5,6f8.3)
          
          call GSCR(IWKID,icol,red,grn,blu)
      enddo

      return
      end

      subroutine hsi_to_rgb(hue,sat,rintens,red,grn,blu)

!     Hue is 0:R, 1:B, 2:G, 3:R

      red1 = max(1.0 - abs(hue - 0.0),0.0)
      red2 = max(1.0 - abs(hue - 3.0),0.0)
      red = max(red1,red2)
      grn = max(1.0 - abs(hue  - 2.0),0.0)
      blu = max(1.0 - abs(hue  - 1.0),0.0)

!     Normalize to the max intensity
      colmax = max(red,grn,blu)
      if(colmax .gt. 0.)then
          red = red/colmax
          grn = grn/colmax
          blu = blu/colmax
      endif

      red = (red*sat) + 1.0*(1.0-sat)
      grn = (grn*sat) + 1.0*(1.0-sat)
      blu = (blu*sat) + 1.0*(1.0-sat)

      red = red * rintens
      grn = grn * rintens
      blu = blu * rintens

      return
      end
