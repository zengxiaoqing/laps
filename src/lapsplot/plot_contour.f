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


      subroutine conrec_line(field,ni,nii,nj
     1                      ,clow_in,chigh_in,cint_in,plot_parms
     1                      ,i1,i2,i3,i4)

!     97-Aug-17     Ken Dritz     Commented out assignment to r_missing_data
!                                 and inserted call to get_r_missing_data

      include 'lapsplot.inc'

      data ihl/0/
      save ihl

      real field(ni,nj)

      character*48 c_blank

      c_blank = '   '

      NF = 51

      call gflas1(NF)
      call gflas2()

!     call color

!     r_missing_data = 1.e37
      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
         write (6,*) 'Error getting r_missing_data in conrec_line'
         stop
      endif

!     ihl = 2
      ihl = ihl + 1

      if(chigh_in .eq. clow_in)then ! special contouring (non-equal intervals)
          clow = clow_in
          chigh = chigh_in
          cint = cint_in
          lis = 1 ! labelling interval
!         ihl = ihl - 1
      else
!         ihl = 1
          clow = clow_in
          chigh = chigh_in
          cint = cint_in
          lis = 2
      endif

      write(6,*)' conrec_line: calling plot_contour'
     1                ,clow,chigh,cint,lis,ihl

      call plot_contour
     1  (NF,ni,nii,nj,field,c_blank,cint,cl,ch,c_blank,4,ihl
     1      ,clow,chigh,plot_parms,lis,r_missing_data)

      return
      end



      SUBROUTINE PLOT_CONTOUR(NF, MX, NX, NY, F, LAB1, FINT, FMIN, FMAX,
     + TIMELABEL,NC,IHL,clow,chigh,plot_parms,lis,r_missing_data)
C
C PURPOSE       Plot contours of F(NX,NY).
C
C ARGUMENTS     NF --- I  flash buffer to be used
C               MX --- I  maximum size of first dimension of F
C               NX --- I  actual size of first dimension of F
C               NY --- I  actual size of second dimension of F
C               F ---- I  two dimensional array to be contoured
C               LAB1 - I  user-supplied plot label
C               FINT - I  contour interval (if 0., .1*(FMAX-FMIN) is used)
C               FMIN -  O minimum grid value
C               FMAX -  O maximum grid value
C               TIMELABEL - I string
C               NC   - I  contour color index
C               IHL  - I  flag as whether to plot high/low labels
C
C


      include 'lapsplot.inc'

      PARAMETER (N2X=400)
      INTEGER NF, MX, NX, NY, IVAR
      REAL F1(N2X,N2X),HOLD(N2X,2)
      REAL F(MX,NY), FINT, FMIN, FMAX, machine_epsilon_p

      CHARACTER*48 LAB3
      CHARACTER*48 TIMELABEL
      CHARACTER*48 LAB1

      PARAMETER (LRWK=64000, LIWK=64000)

      REAL RWRK(LRWK)
      INTEGER IWRK(LIWK)

      integer hic,hlc

      REAL SCALE
      DATA SCALE /1.1/

      parameter (machine_epsilon_p = 1.19e-07)  ! from iftran.im - bj
      COMMON /SAVMAP/
     .  MX1, MX2, MY1, MY2,
     .  U1,  U2,  V1,  V2
      COMMON /FXFY1/ XA, YA, UA, VA, U1A, V1A, DUDI, DVDJ

        COMMON /NXNY/ NX1,NY1
        common /icol_index/ icol_current
        common /zoom/       zoom

c ... compare,spv def'n from iftran.im/b. jewett
      compare(a,b,c)=((int(sign(1.,c-abs(a-b)))+1)/2)
      spv(a) = compare(a,SPVAL_P,0.01)

        SPVAL_P=r_missing_data

        FMIN0 = FMIN
        FMAX0 = FMAX

        NX1 = NX
        NY1 = NY
        NX2X = NX*2-1
        NY2X = NY*2-1

!       MAXCHR = 'X'
!       MINCHR = 'N'

C
C     CALC MAX AND MIN OF F
C
      FMAX = -1.E8
      FMIN =  1.E8
      DO 55 J=1,NY
        DO 60 I=1,NX
!         IF (SPV(F(I,J)) .EQ. NO_P)
          IF (SPV(F(I,J)) .EQ. 0)
     .    THEN
!          if(f(i,j) .ne. r_missing_data)then
            IF (F(I,J).GT.FMAX) FMAX=F(I,J)
            IF (F(I,J).LT.FMIN) FMIN=F(I,J)
!          endif
          END IF
60      CONTINUE
55    CONTINUE
C
C     TEST FOR ZERO FINT
C
      IF (FINT .LE. MACHINE_EPSILON_P) THEN
        DF = .1*(FMAX-FMIN)
      ELSE
        DF = FINT
      END IF

      WRITE(6,5) FMIN,FMAX,DF
5     FORMAT(/,' MIN = ',G12.5,5X,'MAX = ',
     .  G12.5,5X,' PLOT INTERVAL = ',G12.5)
C
C    PLOT IF THERE IS ADEQUATE FIELD VARIATION
C

        write(LAB3,101)  FMAX,FMIN,DF
  101   FORMAT(2X,'MAX=',F10.3,2X,'MIN=',F10.3,2X,'INT=',F8.3,2X)

        call gflas3(NF)

!       CALL SET (0.,1.,0.,1.,0.,1.,0.,1.,1)

        IF (NC.EQ.1) CALL PLCHHQ (0.0,0.88,TIMELABEL,0.018,0.,-1.)
!       CALL PCSETI ('CC - CHARACTER COLOR', 7+NC)
!       CALL PLCHHQ (0.45,0.14-(NC-1)*0.02, LAB3   , .01,
!    1       0., -1.)
        CALL PLCHHQ (0.05,0.14-(NC-1)*0.02, LAB1   , .01,
     1       0., -1.)

!       FMIN = CEILING2(FMIN,DF,MACHINE_EPSILON_P)   ! Minimum contour
        SPVALU=SPVAL_P
!       CALL FXFY(NF)                               ! Initialize FX and FY fcns



        NULBLL = 1         ! number of unlabelled lines between labelled ones
        NHI = -1

C --- Do contouring
!       if (ihl.ge.1) then

        if (.true.) then
!         call CPSETR('HLS  - HIGH/LOW LABEL SIZE',.020)
          call CPSETR('HLS  - HIGH/LOW LABEL SIZE',.013)
          call CPSETC('ILT',' ')

          if(FMAX .eq. 100. .and. FMIN .le. 5.0)then ! RH case
              NSD = 3
              NLS = 1
          elseif(FMAX .ge. 1000. .and. FMAX .le. 1100.
     1     .and. FMIN .gt. 500.  .and. FMIN .lt. 1000.)then ! STNP case
              NSD = 3
              NLS = 1
          else ! general case
              NSD = 2
              NLS = 1
          endif

          call cpseti('NSD',NSD)
          call cpseti('NLS',NLS)

          write(6,*)' FMAX/NSD/NLS = ',FMAX,NSD,NLS

          CALL CPSETI ('HIC', icol_current)
          CALL CPSETI ('LOC', icol_current)
          CALL CPSETI ('HLC', icol_current)

!         if (ihl.ge.2) call cpsetc('HLT',' ')

        end if

!       call cpsetr ('LLS',.040)
        call cpseti ('CLS - contour level selection flag',20)
        call cpsetr ('CIS',DF)
        call cpseti ('LIS - label interval specifier',lis)
!       call cpseti ('LLP', 2)
!       call cpsetr ('CMN',FMIN)
!       call cpsetr ('CMX',FMAX)
        call cpsetr ('CMN',clow)
        call cpsetr ('CMX',chigh)
        call cpsetr ('SPV',SPVALU)
!       call cpsetr ('LLS - LINE LABEL SIZE',.025/sqrt(zoom))
        call cpsetr ('LLS - LINE LABEL SIZE',.025)
        call cpgetr ('LLS - LINE LABEL SIZE',clls)           
        call cpgetr ('HLS',hls)           
        call cpsetr ('CWM',1.00/sqrt(zoom))


        CALL CPRECT (F,NX,NX,NY,RWRK,LRWK,IWRK,LIWK)
        CALL CPPKCL (F,RWRK,IWRK)
       if (ihl.ge.1) then
        call cplbdr (f,rwrk,lwrk)
       end if


        CALL CPGETI ('NCL - NUMBER OF CONTOUR LEVELS', NCON)
        DO 111 I=1,NCON
          CALL CPSETI ('PAI - PARAMETER ARRAY INDEX', I)
          call cpgetr ('CLV - contour level values',cval)
          if (cval.lt.0.) then
            call cpseti ('CLD',61166)
          else
            call cpseti ('CLD',65535)
          end if
          call cpseti ('CLL',plot_parms%contour_line_width)       ! Line Width
!         CALL CPSETI ('CLC - CONTOUR LINE COLOR INDEX', 7+NC)
111     CONTINUE

        CALL CPCLDR (F,RWRK,IWRK)

        isize=7
        CALL PWRIT (MX1,MY2,' ',1,0,0,0)
        CALL SFLUSH                                  ! Flush frame buffer

        CALL cpgeti ('LIU', liu)
        CALL cpgeti ('LLP', llp)
        CALL cpgeti ('HIC', hic)
        CALL cpgeti ('HLC', hlc)
        CALL cpgeti ('LOC', loc)

        write(6,*)'NCON/IHL / HLS / HIC/HLC/LOC/LLP/LLS/LIS/LIU/ICOL = '
        write(6,*) NCON,IHL,hls,hic,hlc,loc,llp,clls,lis,liu
     1             ,icol_current


        WRITE(6,155)
155     FORMAT(' THE FIELD HAS BEEN PLOTTED')


      RETURN
      END

cccccccccccccccccccccccccccccc  ceiling2  cccccccccccccccccccccccccccc
c
c  ceiling2 - ceiling of x, over interval y, to accuracy z
c
        real function ceiling2(x,y,z)
        implicit none
        real x,y,z,trunc,a,b
        trunc(a,b) = (b*int(a/b))
c
        if (x.gt.0) then                ! ceiling(x,y,z) = trunc(x+y-z,y)
          ceiling2=trunc(x+y-z,y)
        else if (x.lt.0) then           ! ceiling(x,y,z) = -trunc(-x,y)
          ceiling2= -trunc(-x,y)
        else                            ! ceiling(x,y,z) = 0.
          ceiling2 = 0.
        endif
c
        return
        end


        SUBROUTINE FILL (XWRK,YWRK,NWRK,IAREA,IGRP,NGRPS)
C
        DIMENSION XWRK(*),YWRK(*),IAREA(*),IGRP(*)

        DO 10, I=1,NGRPS
          IF (IGRP(I).EQ.3) IAREA3=IAREA(I)
 10     CONTINUE

        IF (IAREA3 .GT. 0) THEN
C If the area is defined by 3 or more points, fill it
           CALL GSFACI(IAREA3+1)
           CALL GFA(NWRK,XWRK,YWRK)
        ENDIF
        return
        end

      SUBROUTINE COLOR(iwhite)
!     write(6,*)' White or Black background. [Enter 1/2]    ?   :'
!     accept 5, iwhite
!5     format (i)
!     iwhite = 2

C    BACKGROUND COLOR
C  The background is white here for better visibility on paper

      if (iwhite.eq.1) then
        write(6,*)' Color: White Background'
        CALL GSCR (1,0,1.,1.,1.)
        CALL GSCR (1,1,0.8,0.8,1.)
      else
        write(6,*)' Color: Black Background'
        CALL GSCR (1,0,0.,0.,0.)
        CALL GSCR (1,1,0.,0.,0.)
      end if

C
C     BACKGROUND COLOR
C     BLACK
c     CALL GSCR(1,0,0.,0.,0.)
C
C     FORGROUND COLORS
C White
      CALL GSCR(1,  1, 1.0, 1.0, 1.0)
        CALL GSCR (1,1,0.8,0.8,1.)
C white
      CALL GSCR(1,  2, 1.0, 1.0, 1.0)
C Red
      CALL GSCR(1,  3, 1.0, 0.0, 0.0)
C OrangeRed
      CALL GSCR(1,  4, 1.0, 0.30, 0.0)
C Orange
      CALL GSCR(1,  5, 1.0, 0.65, 0.0)
C Gold
      CALL GSCR(1,  6, 1.0, 0.85, 0.0)
C Yellow
      CALL GSCR(1,  7, 1.0, 1.0, 0.0)
C GreenYellow
      CALL GSCR(1,  8, 0.7, 1.0, 0.2)
C Chartreuse
      CALL GSCR(1,  9, 0.5, 1.0, 0.0)
C Celeste
      CALL GSCR(1, 10, 0.2, 1.0, 0.5)
C Green
!     CALL GSCR(1, 11, 0.2, 0.8, 0.2)
      CALL GSCR(1, 11, 0.1, 0.9, 0.1)
C Aqua
      CALL GSCR(1, 12, 0.0, 0.9, 1.0)
C DeepSkyBlue
      CALL GSCR(1, 13, 0.0, 0.75, 1.0)
!     CALL GSCR(1, 13, 0.0, 0.65, 0.9)
C RoyalBlue
      CALL GSCR(1, 14, 0.25, 0.45, 0.95)
C SlateBlue
      CALL GSCR(1, 15, 0.4, 0.35, 0.8)
C DarkViolet
      CALL GSCR(1, 16, 0.6, 0.0, 0.8)
C Lavender
      CALL GSCR(1, 17, 1.00, 0.65, 1.0)
C Black
      CALL GSCR (1,21,0.,0.,0.)
C Quasi Black
      CALL GSCR (1,22,0.05,0.05,0.05)
C Medium Gray
      CALL GSCR (1,23,.3,.3,.3)
! Unused...
      CALL GSCR (1,24,.0,.7,.0)
! Dark Gray
      CALL GSCR (1,25,.1,.1,.1)
! Pale Orange
      CALL GSCR (1,26,.6,.3,.1)
C Orchid
      CALL GSCR (1,27, 0.85, 0.45, 0.8)
! Unused...
      CALL GSCR (1,28,.0,.0,.95)
! Paler Yellow
      CALL GSCR (1,29,.50,.38,.0)
! Dark Pink
      CALL GSCR (1,30,.5,.30,.5)
      CALL GSCR (1,31,.2,.7,.7)

! Pale Yellow
      CALL GSCR (1,32,.6,.6,.0)
! Pale Red
      CALL GSCR (1,33,.5,.0,.0)
! Grey
      CALL GSCR (1,34,.5,.5,.5)
! DarkViolet
      CALL GSCR(1, 35, 0.6, 0.0, 0.8)

! Radar Vel Color Table (about 20 values)
      i_velcol_offset = 100
      do icol = i_velcol_offset,i_velcol_offset+10
          ivel = i_velcol_offset+10 - icol
          amp = float(ivel) / 10.
          amp_gry = 0.3
          amp_grn = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call GSCR(1,icol,amp_gry,amp_grn,amp_gry)
!         call GSCR(1,icol,0.,1.0,0.)
!         write(6,*)' Setting color ' ,icol
      enddo
      do icol = i_velcol_offset+11,i_velcol_offset+20
          ivel = icol - (i_velcol_offset+10)
          amp = float(ivel) / 10.
          amp_gry = 0.3
          amp_red = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call GSCR(1,icol,amp_red,amp_gry,amp_gry)
!         write(6,*)' Setting color ' ,icol
      enddo

! Radar Ref Color Table
      i_refcol_offset = 180
      do idbz = 0,20
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz) / 20.
          amp_gry = 0.3
          amp_red = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          amp_blu = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call GSCR(1,icol,amp_red,amp_gry,amp_blu)
      enddo
      do idbz = 21,30
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz-20) / 10.
          amp_gry = 0.3
          amp_grn = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call GSCR(1,icol,0.,amp_grn,0.)
      enddo
      do idbz = 31,40
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz-30) / 10.
          amp_gry = 0.3
          amp_blu = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call GSCR(1,icol,0.,0.,amp_blu)
      enddo
      do idbz = 41,50
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz-40) / 10.
          amp_gry = 0.3
          amp_red = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          amp_grn = amp_red * 0.5
          call GSCR(1,icol,amp_red,amp_grn,0.)
      enddo
      do idbz = 51,60
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz-50) / 10.
          amp_gry = 0.3
          amp_red = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          amp_grn = amp_red
          call GSCR(1,icol,amp_red,amp_grn,0.)
      enddo
      do idbz = 61,70
          icol = i_refcol_offset + idbz / 5
          amp = float(idbz-60) / 10.
          amp_gry = 0.3
          amp_red = min(amp,1.0) * (1.0 - amp_gry) + amp_gry
          call GSCR(1,icol,amp_red,0.,0.)
      enddo

!     White
      CALL GSCR(1, 251, 1.0, 1.0, 1.0)

C Done.
C
        RETURN
C
      END
