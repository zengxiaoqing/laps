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

 
       subroutine get_lapsplot_parms(namelist_parms,istatus)

       include 'lapsplot.inc'
 
       character*150 static_dir,filename
       character*3 c3_time_zone
       character*30 c_institution
       character*6 c_vnt_units
       character*7 c_tpw_units
       character*7 c_units_type
       character*7 c_pbl_depth_units
       character*10 c_ob_color
       character*20 btemp_colortable
       logical l_discrete,l_sphere,l_low_fill,l_high_fill       
       real time_zone

       namelist /lapsplot_nl/ latlon_int,continent_line_width
     1                       ,country_line_width
     1                       ,state_line_width
     1                       ,county_line_width
     1                       ,c3_time_zone,time_zone
     1                       ,c_institution,c_vnt_units,c_tpw_units
     1                       ,c_units_type,c_pbl_depth_units
     1                       ,chigh_sfcwind,chigh_3dwind,chigh_cape
     1                       ,chigh_tpw,power_tpw,scale_omega
     1                       ,l_discrete, l_sphere
     1                       ,l_low_fill, l_high_fill       
     1                       ,mode_supmap, iraster, icol_barbs
     1                       ,icol_continent,icol_country
     1                       ,icol_state,icol_county
     1                       ,dist_plot_ua, dist_plot_sfc
     1                       ,c_ob_color, i_background_color
     1                       ,btemp_colortable
     1                       ,i_pcp_sto_colorbar
     1                       ,i_sno_sto_colorbar
     1                       ,montage

!      Set defaults
       latlon_int = 0
       continent_line_width = 1.0
       country_line_width = 1.0
       state_line_width = 1.0
       county_line_width = 1.0
       c3_time_zone = 'UTC'
       time_zone = 0.0
       c_institution = 'NOAA/FSL LAPS'
       c_vnt_units = 'M**2/S'
       c_tpw_units = 'CM'
       mode_supmap = 3
       iraster = 0
       l_sphere = .false.
       icol_barbs = 0
       icol_continent = 7 ! yellow
       icol_country = 7   ! yellow
       icol_state = 7     ! yellow
       icol_county = 7    ! yellow
       chigh_cape = 7000.
       chigh_tpw = 7.
       power_tpw = 0.7
       scale_omega = 100. ! relative to Pa/S
       c_ob_color = 'default'
       i_background_color = 2
       btemp_colortable = 'linear'
       i_pcp_sto_colorbar = 3
       i_sno_sto_colorbar = 4 
       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/lapsplot.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,lapsplot_nl,err=901)
       close(1)

       print*,'success reading lapsplot_nl in ',filename
       write(*,lapsplot_nl)

!      Set namelist structure
       namelist_parms%latlon_int = latlon_int
       namelist_parms%continent_line_width = continent_line_width
       namelist_parms%country_line_width   = country_line_width
       namelist_parms%state_line_width     = state_line_width
       namelist_parms%county_line_width    = county_line_width
       namelist_parms%c3_time_zone = c3_time_zone
       namelist_parms%c_institution = c_institution
       namelist_parms%time_zone = time_zone
       namelist_parms%c_vnt_units = c_vnt_units
       namelist_parms%c_tpw_units = c_tpw_units
       namelist_parms%c_units_type = c_units_type
       namelist_parms%c_pbl_depth_units = c_pbl_depth_units
       namelist_parms%chigh_sfcwind = chigh_sfcwind
       namelist_parms%chigh_3dwind = chigh_3dwind
       namelist_parms%chigh_cape = chigh_cape
       namelist_parms%chigh_tpw = chigh_tpw
       namelist_parms%power_tpw = power_tpw
       namelist_parms%scale_omega = scale_omega
       namelist_parms%c_ob_color = c_ob_color
       namelist_parms%btemp_colortable = btemp_colortable
       namelist_parms%i_background_color = i_background_color
       namelist_parms%l_discrete = l_discrete
       namelist_parms%l_sphere = l_sphere
       namelist_parms%l_low_fill = l_low_fill
       namelist_parms%l_high_fill = l_high_fill
       namelist_parms%mode_supmap = mode_supmap
       namelist_parms%iraster = iraster
       namelist_parms%icol_barbs = icol_barbs
       namelist_parms%icol_continent = icol_continent
       namelist_parms%icol_country = icol_country
       namelist_parms%icol_state = icol_state
       namelist_parms%icol_county = icol_county
       namelist_parms%dist_plot_ua = dist_plot_ua
       namelist_parms%dist_plot_sfc = dist_plot_sfc
       namelist_parms%i_pcp_sto_colorbar = i_pcp_sto_colorbar
       namelist_parms%i_sno_sto_colorbar = i_sno_sto_colorbar

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
