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
     1                          ,icol_sta,icol_cou,ni,nj)
c
c
      real*4 sw(2),ne(2),pl1(2),pl2(2),pl3(2),pl4(2),
     +       polat,polon,rrot

c
      integer*4 jproj,jjlts,jgrid,jus,jdot,ier
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

!     1 means use local version of supmap, 2 means use the NCARGlib version
      mode_supmap=1

      call get_grid_spacing(grid_spacing_m,istatus)
      if(istatus .ne. 1)then
          write(6,*)' no grid spacing, stop in draw_county_map'
          stop
      else
          write(6,*)' Subroutine draw_county_map...',jproj
      endif

      domsize = (max(ni,nj)-1.) * grid_spacing_m / zoom

      call get_lapsplot_parms(latlon_int,continent_line_width,istatus)       

!     Plot Counties
      if(jdot .eq. 1)then
          call gsln(3) ! Dotted
      else
          call gsln(1) ! Solid
      endif

      if(domsize .le. 2500e3)then
          write(6,*)' Plotting Counties'
          call setusv_dum(2HIN,icol_cou)
          jgrid=latlon_int                        ! Draw lat/lon lines?

          if(mode_supmap .eq. 1)then
              jus=-4
              call supmap_local(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                          jjlts,jgrid,jus,jdot,ier)
          else
              iout = 0
              call supmap      (jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                          jjlts,jgrid,iout,jdot,ier)
          endif
          if(ier .ne. 0)write(6,*)' ier = ',ier

          call sflush
      else
          write(6,*)' Large domain, omitting counties'

      endif

      write(6,*)' Plotting States From Counties'
      call gsln(1)
      call setusv_dum(2HIN,icol_sta)
      jgrid=0                                ! Do not draw lat/lon lines

      if(mode_supmap .eq. 1)then
          jus=-8
          call supmap_local(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid,jus,jdot,ier)
      else
          iout = 0
          call supmap      (jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid,iout,jdot,ier)
      endif
      if(ier .ne. 0)write(6,*)' ier = ',ier

      call sflush

      write(6,*)' Plotting Continents'
      call gsln(1)
      call setusv_dum(2HIN,icol_sta)

      jgrid=0                                ! Do not draw lat/lon lines
      call GSLWSC(continent_line_width)

      if(mode_supmap .eq. 1)then
          jus=-1
          call supmap_local(jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid,jus,jdot,ier)
      else
          iout = 2
          call supmap      (jproj,polat,polon,rrot,pl1,pl2,pl3,pl4,
     +                      jjlts,jgrid,iout,jdot,ier)
      endif
      if(ier .ne. 0)write(6,*)' ier = ',ier

      call GSLWSC(1.0)

      call sflush

      return
      end

 
       subroutine get_lapsplot_parms(latlon_int,continent_line_width
     1                              ,istatus)

       namelist /lapsplot_nl/ latlon_int,continent_line_width
 
       character*150 static_dir,filename

!      Set defaults
       latlon_int = 0
       continent_line_width = 1.0
 
       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/lapsplot.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,lapsplot_nl,err=901)
       close(1)

       print*,'success reading lapsplot_nl in ',filename
       write(*,lapsplot_nl)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading lapsplot_nl in ',filename
       write(*,lapsplot_nl)
       istatus = 0
       return

       end
