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
      subroutine draw_county_map(sw,ne,jproj,polat,polon,rrot,jdot
     1                          ,icol_sta,icol_cou,ni,nj,namelist_parms)       

      include 'lapsplot.inc'
c
      real sw(2),ne(2),pl1(2),pl2(2),pl3(2),pl4(2),
     +       polat,polon,rrot

c
!abdel      
      DOUBLE PRECISION GRSP,DPOLAT,DPOLON,DSW(2),DNE(2)
      INTEGER IRGL
      
      DATA GRSP,IRGL / 1.D0 , 0  /
      PARAMETER (NCRA=10000,NGPS=10,LRWK=2*NCRA)

      DIMENSION XCRA(NCRA),YCRA(NCRA)
      DIMENSION RWRK(LRWK)     
      EQUIVALENCE (RWRK(1),XCRA(1)),(RWRK(NCRA+1),YCRA(1))
! abdel     

      integer jproj,jjlts,jgrid,jus,jdot,ier
c
      COMMON/SUPMP9/DS,DI,DSRDI
      common /zoom/       zoom 
!     DI = 50.
!     polat=90.

!     rrot=0.
      pl1(1)=sw(1)
      pl2(1)=sw(2)
      pl3(1)=ne(1)
      pl4(1)=ne(2)
      jjlts=-2
      jgrid=0
 

!     call get_lapsplot_parms(namelist_parms,istatus)       

!     1 means use local version of supmap, 2 means use the NCARGlib version
!     3 means to try newer (ezmap) routines
      mode_supmap = namelist_parms%mode_supmap

      call get_grid_spacing(grid_spacing_m,istatus)
      if(istatus .ne. 1)then
          write(6,*)' no grid spacing, stop in draw_county_map'
          stop
      else
          write(6,*)
          write(6,*)' Subroutine draw_county_map...',mode_supmap,jproj
      endif

      domsize = (float(nj)-1.) * grid_spacing_m / zoom

!     Plot Counties
      if(jdot .eq. 1)then
          call gsln(3) ! Dotted
      else
          call gsln(1) ! Solid
      endif

      jgrid=namelist_parms%latlon_int        ! Draw lat/lon lines?

      if(domsize .le. 1500e3 .and. 
     1   namelist_parms%continent_line_width .gt. 0)then
          write(6,*)' Plotting Counties ',domsize,mode_supmap,jgrid
          call setusv_dum(2HIN,icol_cou)

          if(mode_supmap .eq. 1)then
              jus=-4
              call supmap_local(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                          jjlts,jgrid,jus,jdot,ier)
          elseif(mode_supmap .eq. 2)then
              iout = 0
              call supmap      (jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                          jjlts,jgrid*1000,iout,jdot,ier)
          elseif(mode_supmap .eq. 3)then
              write(6,*)' Calling MPLNDR, etc. for counties...'
!             call mapdrw()
              call mapint
!             call maplot
              CALL MPLNDR ('Earth..2',5)
              if(jgrid .gt. 0)then ! draw lat/lon lines
                  call mpsetr('GR',float(jgrid))
                  call mapgrd()
              endif
c abdel	      
          elseif(mode_supmap .eq. 4)then
              write(6,*)' Calling SUB submap=4 for europe...'
              DPOLAT=polat
              DPOLON=polon
              DSW(1)=sw(1)
              DSW(2)=sw(2)
              DNE(1)=ne(1)
              DNE(2)=ne(2)
              CALL MDPROJ ('ST',DPOLAT,DPOLON,0.D0)
              CALL MDPSET ('CO',DSW(1),DSW(2),DNE(1),DNE(2))
              CALL MAPSTD ('GR',GRSP)
              CALL MDRGOL (IRGL,RWRK,LRWK)    
          
          elseif(mode_supmap .eq. 5)then
              write(6,*)' Accessing RANGS database...'
              write(6,*)' Not yet supported - stop' 
              stop
          endif
          if(ier .ne. 0)write(6,*)' ier = ',ier

          call sflush

          jgrid=0                        ! Do not draw subsequent lat/lon lines

      else
          write(6,*)' Omitting counties ',domsize
     1             ,namelist_parms%continent_line_width

      endif

      call gsln(1)
      call setusv_dum('IN',namelist_parms%icol_state)

      call GSLWSC(namelist_parms%continent_line_width)

      if(mode_supmap .eq. 1)then
          write(6,*)' Plotting States From Counties ',mode_supmap,jgrid
          jus=-8
          call supmap_local(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid,jus,jdot,ier)
      elseif(mode_supmap .eq. 2)then
          write(6,*)' Plotting States From Counties ',mode_supmap,jgrid
          iout = 0
          call supmap      (jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid*1000,iout,jdot,ier)
      elseif(mode_supmap .eq. 3)then
          call mapint

          if(namelist_parms%state_line_width .gt. 0. .OR.
     1       namelist_parms%country_line_width .gt. 0.     )then
              write(6,*)' Calling MAPDRW, etc. for countries/states...'
     1                 ,namelist_parms%state_line_width
     1                 ,namelist_parms%icol_country
     1                 ,icol_sta
              call GSLWSC(namelist_parms%state_line_width)
              call setusv_dum('IN',namelist_parms%icol_state) 
              CALL MPLNDR ('Earth..2',4) ! states & countries
          else
              write(6,*)' Skip plot of countries/states'
              if(namelist_parms%continent_line_width .gt. 0.)then
                  write(6,*)' Calling MAPDRW just for continents...'
                  call GSLWSC(namelist_parms%continent_line_width)
                  call setusv_dum('IN',namelist_parms%icol_continent) 
                  CALL MPLNDR ('Earth..2',2) ! continents
              endif
          endif

          if( (namelist_parms%icol_country .ne. 
     1         namelist_parms%icol_state) 
     1                     .OR.
     1        (namelist_parms%country_line_width .eq. 0. 
     1                     .AND.
     1         namelist_parms%state_line_width .gt. 0.)       
     1                                                 )then
              write(6,*)
     1              ' Replotting countries in separate color/width: '
     1                 ,0
     1                 ,namelist_parms%state_line_width
              call GSLWSC(namelist_parms%state_line_width)
              if(namelist_parms%country_line_width .eq. 0.)then
                  call setusv_dum('IN',0) 
              else
                  call setusv_dum('IN',namelist_parms%icol_country) 
              endif
              CALL MPLNDR ('Earth..2',3) ! countries
          endif

          if(jgrid .gt. 0)then ! draw lat/lon lines
              call setusv_dum(2HIN,icol_cou)
              call mpsetr('GR',float(jgrid))
              call mapgrd()
          endif
c abdel	  
      elseif(mode_supmap .eq. 4)then
            write(6,*)' Calling SUB submap=4 for europe...'
            DPOLAT=polat
            DPOLON=polon
            DSW(1)=sw(1)
            DSW(2)=sw(2)
            DNE(1)=ne(1)
            DNE(2)=ne(2)
            CALL MDPROJ ('ST',DPOLAT,DPOLON,0.D0)
            CALL MDPSET ('CO',DSW(1),DSW(2),DNE(1),DNE(2))
            CALL MAPSTD ('GR',GRSP)
            CALL MDRGOL (IRGL,RWRK,LRWK)	  
	  
      else
          write(6,*)' Accessing RANGS database...'
          write(6,*)' Not yet supported - stop' 
          stop

      endif

      if(ier .ne. 0)write(6,*)' ier = ',ier

      call sflush

      call gsln(1)
      call setusv_dum(2HIN,icol_sta)

      jgrid=0                                ! Do not draw lat/lon lines
       
      if(mode_supmap .eq. 1)then
          write(6,*)' Plotting Continents ',mode_supmap,jgrid
          jus=-1
          call supmap_local(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid,jus,jdot,ier)
      elseif(mode_supmap .eq. 2)then
          write(6,*)' Plotting Continents ',mode_supmap,jgrid
          iout = 2
          call supmap      (jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid*1000,iout,jdot,ier)
      endif
      if(ier .ne. 0)write(6,*)' ier = ',ier

      call GSLWSC(1.0)

      call sflush

      write(6,*)

      return
      end

cabdel       
              SUBROUTINE MDRGDI (DINM)
C
C This is a user-replaceable routine that returns the name of the
C directory in which the RANGS/GSHHS data files have been placed.
C
       CHARACTER*(*) DINM

C Fitxer bo
C Return the name of the directory where the RANGS/GSHHS data reside.
C
       DINM='/usr/local/ncarg/lib/ncarg/database/RANGS_GSHHS'
C
C Done.
C
       RETURN

       END
