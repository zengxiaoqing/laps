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
     1                 ,colortable,n_image,scale,c5_sect)       

C 
C Define error file, Fortran unit number, and workstation type,
C and workstation ID.
C 
      PARAMETER (IERRF=6, LUNIT=2, IWTYPE=1, IWKID=1)
      REAL XREG(MREG),YREG(NREG),ZREG(MREG,NREG),field_in(MREG,NREG)
      character*(*)colortable
      character*5 c5_sect

      logical log_scaling, l_discrete

      write(6,*)' Subroutine ccpfil for solid fill plot...'

      if(colortable(1:3) .eq. 'acc')then
          log_scaling = .true.
      else
          log_scaling = .false.
      endif

      if(colortable .eq. 'haines')then
          l_discrete = .true.
      else
          l_discrete = .false.
      endif

      n_image = n_image + 1

      call get_r_missing_data(r_missing_data,istatus)

!     if(n_image .gt. 1)then
!         write(6,*)' Image was already plotted - returning from ccpfil'
!         return
!     endif

      if(scale_l_in .lt. scale_h_in)then
          ireverse = 0
          scale_l = scale_l_in * scale
          scale_h = scale_h_in * scale
      else
          ireverse = 1
          scale_l = scale_h_in * scale
          scale_h = scale_l_in * scale
      endif

      write(6,*)' Colortable is ',colortable,scale_l,scale_h,ireverse

!     Apply scaling to the array
!     call addcon(field_in,-scale_l,ZREG,MREG,NREG)

      if(log_scaling)then ! e.g. for precip
          scale_l = .01 * scale

          do i = 1,MREG
          do j = 1,NREG
            ZREG(i,j) = alog10(max(field_in(i,j),scale_l))
          enddo ! j
          enddo ! i

          scale_l = alog10(scale_l)
          scale_h = alog10(scale_h)

          ZREG = ZREG - scale_l                          ! Array subtraction

          scale_loc = scale_h - scale_l

      else
          scale_loc = scale_h - scale_l

          ZREG = field_in - scale_l                      ! Array subtraction

!         Adjust field values
          do i = 1,MREG
          do j = 1,NREG
            if(field_in(i,j) .eq. r_missing_data)then

              if(c5_sect .eq. 'hsect')then

                if(colortable(1:3) .eq. 'lin')then
                  ZREG(i,j) = scale_loc * 0.50 ! e.g. CSC 

                elseif(colortable(1:3) .eq. 'cpe')then ! apply only to hsects
                  if(scale_h_in .eq. 7200. .or. scale_h_in .eq. 50.)then 
                      ZREG(i,j) = scale_loc * 0.00 ! ! cape/cin for hsect
                  endif

                elseif(l_discrete)then ! e.g. Haines
                  ZREG(i,j) = scale_loc * 100. ! 1.2 

                else

                  ZREG(i,j) = scale_loc * 0.96 ! e.g. CIN
                endif

              endif

            elseif(ireverse .eq. 1)then
              ZREG(i,j) = scale_loc - ZREG(i,j)

            endif

!           Prevent overshoot beyond colortable (except for CAPE/CIN/discrete)
            if(c5_sect .eq. 'hsect')then
              if(ZREG(i,j) .gt. scale_loc     .and. 
     1             colortable(1:3) .ne. 'cpe' .and.
     1             (.not. l_discrete)               )then
                ZREG(i,j) = scale_loc
              endif
            endif

            if(ZREG(i,j) .lt. 0.0      )then
              ZREG(i,j) = 0.0
            endif

          enddo ! j
          enddo ! i

      endif ! log_scaling


      ireverse_colorbar = ireverse

      ireverse = 0  ! Turn off later use of ireverse
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
      icol_offset = 40 ! Offset new colortable to preserve previous low end

      LMAP=MREG*NREG*256 ! 16000000
      LMAP = min(LMAP,32000000)
      CALL CCPFIL_SUB(ZREG,MREG,NREG,-15,IWKID,scale_loc,ireverse
     1                               ,r_missing_data
     1                               ,LMAP,log_scaling
     1                               ,colortable,ncols,icol_offset)      
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

c     Call local colorbar routine
      write(6,*)' Drawing colorbar: ',MREG,NREG
      call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
      call colorbar(MREG, NREG, ncols, ireverse_colorbar, log_scaling,
     1              scale_l, scale_h, colortable, scale,icol_offset,
     1              c5_sect)

      jdot = 1
      
      return
      END

      
      SUBROUTINE CCPFIL_SUB(ZREG,MREG,NREG,NCL,IWKID,scale,ireverse
     1                                ,r_missing_data
     1                                ,LMAP,log_scaling
     1                                ,colortable,ncols,icol_offset)      
      
      PARAMETER (LRWK=300000,LIWK=300000,NWRK=300000
     1          ,NOGRPS=5)       
      REAL ZREG(MREG,NREG),RWRK(LRWK), XWRK(NWRK), YWRK(NWRK)
      INTEGER MREG,NREG,IWRK(LIWK)
      INTEGER MAP(LMAP),IAREA(NOGRPS),IGRP(NOGRPS)
      character*(*) colortable
      logical log_scaling
      
      EXTERNAL FILL

      ncols = 20
C      
C Set up color table
      write(6,*)' ccpfil_sub - scale = ',scale
C      
      CALL set_image_colortable(IWKID,ncols,ireverse,colortable
     1                                     ,log_scaling,icol_offset)
C      
C Initialize Areas
C      
      CALL ARINAM(MAP, LMAP)

      col_offset = float(icol_offset) / float(ncols)
      write(6,*)' col_offset / scale = ',col_offset,scale

      do m = 1,MREG
      do n = 1,NREG
          ZREG(m,n) = ZREG(m,n) + (col_offset * scale)
      enddo ! n
      enddo ! m

C      
C Set number of contour levels and initialize Conpack
C      
!      CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',NCL)

      cis = abs(scale) / float(ncols)
      CALL CPSETI('CLS - CONTOUR LEVEL SELECTION FLAG',+1)
      CALL CPSETR('CIS', cis)
      CALL CPSETR('CMN',(0.0           ) * abs(scale) + 2.0*cis)
      CALL CPSETR('CMX',(1.0+col_offset) * abs(scale) + 2.0*cis)
      call cpsetr ('SPV',r_missing_data)

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
      
      SUBROUTINE set_image_colortable(IWKID,ncols,ireverse,colortable
     1                               ,log_scaling,icol_offset)    

      character*(*) colortable
      logical log_scaling
C 
C BACKGROUND COLOR
C BLACK

C
!     CALL GSCR(IWKID,0,0.,0.,0.)

      if(colortable(1:3) .eq. 'lin')then
          if(colortable .eq. 'linear_reduced')then
              ncols = 8
          else
              ncols = 20
          endif

          rcols = ncols - 1
          do i = 1,ncols+10 ! 255-icol_offset
              if(colortable .eq. 'linear_reduced')then
!                 rintens = min(max( (float(i) / rcols) - 0.4 ,0.),1.)       
                  rintens = min(max(float(i-1) / rcols,0.),1.)
              else
                  rintens = min(max(float(i-2) / rcols,0.),1.)
              endif

              if(ireverse .eq. 1)rintens = 1.0 - rintens
              call GSCR(IWKID, i+icol_offset, rintens, rintens, rintens)
          enddo ! i

      elseif(colortable .eq. 'hues' .or. colortable .eq. 'ref'
     1                              .or. colortable .eq. 'cpe')then       
          ncols1 = 50
          ncols = 60
          call color_ramp(1,ncols1/8,IWKID,icol_offset
     1                   ,0.5,0.15,0.6                ! Pink
     1                   ,0.5,0.5,0.7)                ! Violet
          call color_ramp(ncols1/8,59*ncols1/120,IWKID,icol_offset
     1                   ,0.5,0.5,0.7                 ! Violet
     1                   ,1.5,1.0,0.7)                ! Aqua
          call color_ramp(59*ncols1/120,73*ncols1/120,IWKID,icol_offset       
     1                   ,1.5,1.0,0.7                 ! Aqua
     1                   ,2.0,0.4,0.4)                ! Green
          call color_ramp(73*ncols1/120,70*ncols1/100,IWKID,icol_offset
     1                   ,2.0,0.4,0.4                 ! Green
!    1                   ,2.5,0.65,0.55)              ! Yellow Orig
     1                   ,2.45,0.95,0.60)             ! Yellow New
          call color_ramp(70*ncols1/100,82*ncols1/100,IWKID,icol_offset
     1                   ,2.45,0.95,0.60              ! Yellow 
     1                   ,2.72,0.9,0.7)               ! Orange
          call color_ramp(82*ncols1/100,90*ncols1/100,IWKID,icol_offset
     1                   ,2.72,0.9,0.7                ! Orange 
     1                   ,3.0,0.9,0.7)                ! Red
          call color_ramp(90*ncols1/100,ncols1,IWKID,icol_offset
     1                   ,3.0,0.9,0.7                 ! Red
     1                   ,3.0,0.9,0.2)                ! Dark Hot
          call color_ramp(100*ncols1/100,ncols,IWKID,icol_offset
     1                   ,3.0,0.9,0.2                 ! Dark Hot
     1                   ,3.0,0.15,0.6)               ! White Hot

          if(colortable .eq. 'ref')then
              do i = 1,3
                  call GSCR(IWKID, i+icol_offset, 0., 0., 0.)
              enddo 

              do i = ncols,ncols
                  call GSCR(IWKID, i+icol_offset, 0.3, 0.3, 0.3)
              enddo

          elseif(colortable .eq. 'cpe')then
              do i = 1,1
                  call GSCR(IWKID, i+icol_offset, 0., 0., 0.)
              enddo 

              ncols = 48

!             call GSCR(IWKID, ncols+icol_offset, 0.3, 0.3, 0.3)

          endif

      elseif(colortable .eq. 'vnt')then       
          ncols = 60
          call color_ramp(1,7,IWKID,icol_offset
     1                   ,0.5,0.15,0.6                ! Pink
     1                   ,0.5,0.5,0.7)                ! Violet
          call color_ramp(7,29,IWKID,icol_offset
     1                   ,0.5,0.5,0.7                 ! Violet
     1                   ,1.5,1.0,0.7)                ! Aqua
          call color_ramp(29,36,IWKID,icol_offset       
     1                   ,1.5,1.0,0.7                 ! Aqua
     1                   ,2.0,0.4,0.4)                ! Green
          call color_ramp(36,42,IWKID,icol_offset
     1                   ,2.0,0.4,0.4                 ! Green
!    1                   ,2.5,0.65,0.55)              ! Yellow Orig
     1                   ,2.45,0.95,0.60)             ! Yellow New
          call color_ramp(42,49,IWKID,icol_offset
     1                   ,2.45,0.95,0.60              ! Yellow 
     1                   ,2.72,0.9,0.7)               ! Orange
          call color_ramp(49,54,IWKID,icol_offset
     1                   ,2.72,0.9,0.7                ! Orange 
     1                   ,3.0,0.9,0.7)                ! Red
          call color_ramp(54,60,IWKID,icol_offset
     1                   ,3.0,0.9,0.7                 ! Red
     1                   ,3.0,0.9,0.2)                ! Dark Hot

      elseif(colortable .eq. 'tpw')then
          ncols = 60
          call color_ramp(1,8,IWKID,icol_offset
     1                   ,3.0,0.9,0.2                 ! Dark Hot
     1                   ,3.0,0.9,0.7)                ! Red
          call color_ramp(8,12,IWKID,icol_offset
     1                   ,3.0,0.9,0.7                 ! Red
     1                   ,2.75,0.9,0.7)               ! Orange
          call color_ramp(12,18,IWKID,icol_offset
     1                   ,2.75,0.9,0.7                ! Orange
     1                   ,2.5,0.95,0.60)              ! Yellow
          call color_ramp(18,29,IWKID,icol_offset
     1                   ,2.5,0.95,0.60               ! Yellow
     1                   ,2.0,0.4,0.4)                ! Green
          call color_ramp(29,36,IWKID,icol_offset       
     1                   ,2.0,0.4,0.4                 ! Green
     1                   ,1.5,1.0,0.7)                ! Aqua
          call color_ramp(36,54,IWKID,icol_offset
     1                   ,1.5,1.0,0.7                 ! Aqua
     1                   ,0.5,0.5,0.7)                ! Violet
          call color_ramp(54,60,IWKID,icol_offset
     1                   ,0.5,0.5,0.7                 ! Violet
     1                   ,0.5,0.15,0.6)               ! Pink

      elseif(colortable .eq. 'moist')then
          ncols = 60
          call color_ramp(1,12,IWKID,icol_offset
     1                   ,3.0,0.9,0.2                 ! Dark Hot
     1                   ,3.0,0.9,0.7)                ! Red
          call color_ramp(12,17,IWKID,icol_offset
     1                   ,3.0,0.9,0.7                 ! Red
     1                   ,2.75,0.9,0.7)               ! Orange
          call color_ramp(17,30,IWKID,icol_offset
     1                   ,2.75,0.9,0.7                ! Orange
     1                   ,2.5,0.95,0.60)              ! Yellow
          call color_ramp(30,47,IWKID,icol_offset
     1                   ,2.5,0.95,0.60               ! Yellow
     1                   ,2.0,0.4,0.4)                ! Green
          call color_ramp(47,60,IWKID,icol_offset       
     1                   ,2.0,0.4,0.4                 ! Green
     1                   ,1.5,1.0,0.7)                ! Aqua

      elseif(colortable .eq. 'spectral' .or. colortable .eq. 'acc')then       
          ncols = 40
          call color_ramp(1,35*ncols/100
     1                   ,IWKID,icol_offset
     1                   ,0.6,0.7,0.4                 ! Violet
     1                   ,1.5,1.0,0.7)                ! Aqua
          call color_ramp(35*ncols/100,55*ncols/100
     1                   ,IWKID,icol_offset       
     1                   ,1.5,1.0,0.7                 ! Aqua
     1                   ,2.0,0.4,0.4)                ! Green
          call color_ramp(55*ncols/100,74*ncols/100
     1                   ,IWKID,icol_offset
     1                   ,2.0,0.4,0.4                 ! Green
     1                   ,2.5,0.95,0.60)              ! Yellow
          call color_ramp(74*ncols/100,91*ncols/100
     1                   ,IWKID,icol_offset
     1                   ,2.5,0.95,0.60               ! Yellow
     1                   ,2.75,0.9,0.7)               ! Orange
          call color_ramp(91*ncols/100,ncols
     1                   ,IWKID,icol_offset
     1                   ,2.75,0.9,0.7                ! Orange
     1                   ,3.0,0.9,0.7)                ! Red

      elseif(colortable .eq. 'haines')then       
          ncols = 5 
          call color_ramp(1,1
     1                   ,IWKID,icol_offset
     1                   ,0.6,0.7,0.4                 ! Violet
     1                   ,0.6,0.7,0.4)                ! Violet
          call color_ramp(2,2
     1                   ,IWKID,icol_offset       
     1                   ,1.0,0.85,0.55               ! Blue
     1                   ,1.0,0.85,0.55)              ! Blue
          call color_ramp(3,3
     1                   ,IWKID,icol_offset
     1                   ,2.0,0.4,0.4                 ! Green
     1                   ,2.0,0.4,0.4)                ! Green
          call color_ramp(4,4
     1                   ,IWKID,icol_offset
     1                   ,2.5,0.95,0.60               ! Yellow
     1                   ,2.5,0.95,0.60)              ! Yellow
          call color_ramp(5,5
     1                   ,IWKID,icol_offset
     1                   ,3.0,0.9,0.7                 ! Red
     1                   ,3.0,0.9,0.7)                ! Red

!         Extra color at the top for r_missing_data
          call color_ramp(6,6
     1                   ,IWKID,icol_offset
     1                   ,1.0,0.0,1.0                 ! White
     1                   ,1.0,0.0,1.0)                ! White

      else
          write(6,*)' ERROR: Unknown color table ',colortable

      endif

      if(colortable .eq. 'acc')then ! Set colortable ends
          do i = 1,1
!         do i = 1,3
              call GSCR(IWKID, i+icol_offset, 0., 0., 0.)
          enddo 

          do i = ncols,ncols
              call GSCR(IWKID, i+icol_offset, 0.3, 0.3, 0.3)
          enddo
      endif
C 
      RETURN
      END

      subroutine color_ramp(ncol1,ncol2,IWKID,icol_offset       ! I
     1                     ,hue1,sat1,rintens1                  ! I
     1                     ,hue2,sat2,rintens2)                 ! I

      write(6,*)' Subroutine color_ramp.. ',ncol1,ncol2

      do icol = ncol1,ncol2
          if(ncol2 .ne. ncol1)then
              frac = float(icol-ncol1) / float(ncol2-ncol1)
          else
              frac = 0.0
          endif   

          hue     = (1.0 - frac) * hue1     + frac * hue2  
          sat     = (1.0 - frac) * sat1     + frac * sat2
          rintens = (1.0 - frac) * rintens1 + frac * rintens2

          call hsi_to_rgb(hue,sat,rintens,red,grn,blu)

          write(6,1)icol,hue,sat,rintens,red,grn,blu
 1        format(i5,6f8.3)
          
          call GSCR(IWKID,icol+icol_offset,red,grn,blu)
      enddo

      if(ncol2+2+icol_offset .le. 255)then
          call GSCR(IWKID,ncol2+1+icol_offset,red,grn,blu)
          call GSCR(IWKID,ncol2+2+icol_offset,red,grn,blu)
      endif

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


      subroutine colorbar(ni,nj,ncols,ireverse,log_scaling
     1                   ,scale_l,scale_h
     1                   ,colortable,scale,icol_offset,c5_sect)

      character*8 ch_low, ch_high, ch_mid, ch_frac
      character*(*)colortable
      character*5 c5_sect
      logical log_scaling,l_loop, l_discrete

      write(6,*)' colorbar: scale_l,scale_h,scale',scale_l,scale_h,scale

      range = abs(scale_h - scale_l) / scale

      if(scale_l .eq. -20. .or. scale_h .eq. 7200. 
     1                     .or. colortable .eq. 'haines' ! HAH, HAM
     1                     .or. colortable .eq. 'tpw'    ! TPW
     1                     .or. colortable .eq. 'vnt'    ! VNT
     1                     .or. range .eq. 100.    )then ! SFC T, Td, RH, CAPE
          l_loop = .true.
      else
          l_loop = .false.
      endif

      if(colortable(1:3) .eq. 'lin')l_loop = .false.

      if(colortable .eq. 'haines')then
          l_discrete = .true.
      else
          l_discrete = .false.
      endif

      call get_border(ni,nj,x_1,x_2,y_1,y_2)

      xlow =  0.40 ! 0.35
      xhigh = xlow + 0.50

      if(c5_sect .eq. 'xsect')then
          y_2 = 0.81
      endif

      ylow =  y_2 + .01
      yhigh = y_2 + .03

      ilow = 1
      ihigh = 999

      xrange = xhigh - xlow
      irange = ihigh - ilow

!     Put Colorbar
      do i = ilow,ihigh
          frac = float(i-ilow) / float(irange)
          x1   = xlow + frac*xrange 
          x2   = xlow + frac*xrange 

          if(ireverse .eq. 0)then
              rcol = 0.5 + float(ncols) * frac
          else
              rcol = 0.5 + float(ncols) * (1.0 - frac)
          endif

          icol = nint(rcol)

          call setusv_dum(2hIN,icol+icol_offset)

          y1 = ylow
          y2 = yhigh
          call line(x1,y1,x2,y2)
      enddo ! i

c     Restore original color table
!     call color

!     Write labels at middle and ends of colorbar
      call setusv_dum(2hIN,34) ! Gray

      call line(xlow,ylow,xhigh,ylow)
      call line(xlow,yhigh,xhigh,yhigh)
      call line(xlow,ylow,xlow,yhigh)
      call line(xhigh,ylow,xhigh,yhigh)

      call setusv_dum(2hIN,7)  ! Yellow

      rsize = .008
      iy = (y_2+.021) * 1024

      if(.not. l_discrete)then

!         Left Edge
          if(log_scaling)then
              rlow = 0.
          else
              rlow = scale_l/scale
          endif

          if(abs(rlow) .gt. 0.0 .and. 
     1       abs(rlow) .le. 0.5                  )then
              write(ch_low,3)rlow
              call right_justify(ch_low)
          else
              write(ch_low, 1)nint(rlow)
              call right_justify(ch_low)
          endif

!         ixl = 353 + nint(.05 * 1024.)
          ixl = nint((xlow - .005) * 1024.)
          CALL PCHIQU (  cpux(ixl),cpux(iy),ch_low,rsize ,0,+1.0)

!         Right Edge
          if(log_scaling)then
              rhigh = (10.**scale_h) / scale
          else
              rhigh = scale_h / scale
          endif

          if(abs(rhigh) .ge. 1.0)then
              if(abs(rhigh-float(nint(rhigh))) .lt. .05 .or. 
     1                                     colortable .eq. 'haines')then       
                  write(ch_high,1)nint(rhigh)
 1                format(i8)
              else
                  write(ch_high,2)rhigh
              endif
          else
              write(ch_high,3)rhigh
 3            format(f8.2)
          endif
          call left_justify(ch_high)

          ixh = ixl + 525 ! 878
          CALL PCHIQU (cpux(ixh),cpux(iy),ch_high,rsize,0,-1.0)

      endif ! l_discrete

      if(.not. l_loop .and. .not. l_discrete)then ! Plot Midpoint

          frac = 0.5
          x1   = xlow + frac*xrange 
          x2   = xlow + frac*xrange 
          call setusv_dum(2hIN,0)

          y1 = ylow
          y2 = yhigh
          call line(x1,y1,x2,y2)

!         Plot Number
          call setusv_dum(2hIN,7)  ! Yellow
          if(log_scaling)then
              rmid = (10.** ((scale_l+scale_h) / 2.0) ) / scale
          else
              rmid = ((scale_l+scale_h) / scale)/2.0
          endif

          if( (abs(rmid) .gt. 1.0 .or. abs(rlow) .gt. 1.0
     1                            .or. abs(rhigh) .gt. 1.0 )
     1                            .AND. 
     1                abs(rmid-float(nint(rmid))) .lt. .05   
     1                                                           )then       
              write(ch_mid,1)nint(rmid)
          elseif(abs(rhigh) .ge. 1.0)then
              write(ch_mid,2)rmid
2             format(f8.1)
          else
              write(ch_mid,3)rmid
          endif 
          call left_justify(ch_mid)
          call s_len(ch_mid,len_mid)

          ixm = (ixl+ixh)/2
          CALL PCHIQU (cpux(ixm),cpux(iy),ch_mid(1:len_mid),rsize,0,0.0)       

      endif

      if(colortable .eq. 'spectral' .or. colortable .eq. 'acc' 
     1                              .or. colortable(1:3) .eq. 'lin')then       

!         Other fractions
          do frac = 0.25,0.75,0.50

!             Plot Black Line
              x1   = xlow + frac*xrange 
              x2   = xlow + frac*xrange 
              call setusv_dum(2hIN,0)

              y1 = ylow
              y2 = yhigh
              call line(x1,y1,x2,y2)

!             Plot Number
              call setusv_dum(2hIN,7)  ! Yellow
              rarg = scale_l + (scale_h-scale_l) * frac
              if(log_scaling)then
                  rlabel = (10.**(rarg)) / scale
              else
                  rlabel = rarg / scale
              endif

              if(rlabel .lt. 0.999)then
                  write(ch_frac,3)rlabel
              elseif( (rlabel .lt. 2.0 .or. rlabel .ne. nint(rlabel))
     1                         .AND.   colortable .ne. 'haines'    )then       
                  write(ch_frac,2)rlabel
              else
                  write(ch_frac,1)nint(rlabel)
              endif

              call left_justify(ch_frac)
              call s_len(ch_frac,len_frac)

              ixm = ixl + (ixh-ixl)*frac
              CALL PCHIQU (cpux(ixm),cpux(iy),ch_frac(1:len_frac)
     1                    ,rsize,0 , 0.0)       
          enddo

      elseif(l_loop)then ! plot additional numbers

!         Set interval for writing numbers
          if(l_discrete)then
              colorbar_int = 0.5
          elseif(range .gt. 1000.)then
              colorbar_int = 1000.
          elseif(range .gt. 10.)then
              colorbar_int = 10.
          else ! range .le. 10
              colorbar_int = 1.0
          endif

!         Interval for writing lines
          colorbar_int = colorbar_int * scale / 1.0

          ixl = 409
          ixh = 924

          write(6,*)' Plotting colorbar',scale_l,colorbar_int,ixl,ixh          

          loop_count = 0

          do rarg = scale_l+colorbar_int,scale_h-.001,colorbar_int
              frac = (rarg - scale_l) / (scale_h - scale_l)

              loop_count = loop_count + 1

              if(.not. l_discrete)then ! Plot Black Line 
                  x1   = xlow + frac*xrange 
                  x2   = xlow + frac*xrange 
                  call setusv_dum(2hIN,0)

                  y1 = ylow
                  y2 = yhigh
                  call line(x1,y1,x2,y2)
              endif

              if(loop_count .eq. (loop_count/1) * 1 )then ! number every line
!                 Plot Number
                  call setusv_dum(2hIN,7)  ! Yellow
                  if(log_scaling)then
                      rlabel = (10.**(rarg)) / scale
                  else
                      rlabel = rarg / scale
                  endif

                  if(range .lt. 0.2)then
                      write(ch_frac,3)rlabel
                  elseif(range .lt. 2.0)then
                      write(ch_frac,2)rlabel
                  else
                      write(ch_frac,1)nint(rlabel)
                  endif

                  call left_justify(ch_frac)
                  call s_len(ch_frac,len_frac)

                  ixm = ixl + (ixh-ixl)*frac

                  if(l_discrete)then
                      if(rlabel .eq. float(nint(rlabel)))then
                          CALL PCHIQU (cpux(ixm),cpux(iy)
     1                        ,ch_frac(1:len_frac),rsize,0 , 0.0)       
                      endif
                  else
                      CALL PCHIQU (cpux(ixm),cpux(iy)
     1                    ,ch_frac(1:len_frac),rsize,0 , 0.0)       
                  endif

              endif ! loop_count
          enddo ! rarg
      endif

      return
      end 
