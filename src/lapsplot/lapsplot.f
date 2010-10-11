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
        subroutine lapsplot(field,ni,nj,clow,chigh,cint,plot_parms
     1                     ,namelist_parms,lat,lon,c_metacode,jdot_in)       

        include 'lapsplot.inc'

        common /supmp1/ dummy,part
        common /supmp6/ umin,umax,vmin,vmax
        real dummy(8),part
        common /ERROR/ IFRAME, IERRR
        common /image/ n_image

        integer       idummy(6),ioffm,istatus

        integer ni,nj,i

        real field(ni,nj),lat(ni,nj),lon(ni,nj)

        real
     1  umin,umax,vmin,vmax

        integer N_CONTOURS
        parameter (N_CONTOURS = 26)
        real factor(N_CONTOURS)
        data factor/
     1  .0001,
     1  .0002,
     1  .0005,
     1  .001,
     1  .002,
     1  .005,
     1  .01,
     1  .02,
     1  .05,
     1  .1,
     1  .2,
     1  .5,
     1  1.,
     1  2.,
     1  5.,
     1  10.,
     1  20.,
     1  50.,
     1  100.,
     1  200.,
     1  500.,
     1  1000.,
     1  2000.,
     1  5000.,
     1  10000.,
     1  20000.
     1                  /


        real cint,cvalue,clow,chigh

        integer
     1  jlts,
     1  jnj,
     1  iusout,
     1  jdot,
     1  ier,iframe

        real
     1  tx,ty,
     1  polon,
     1  rot,
     1  plm1(2),plm2(2),plm3(2),plm4(2)

        character*2 c_metacode


        call get_r_missing_data(r_missing_data, istatus)
        if(istatus .ne. 1)then
            write(6,*)' Bad istatus in lapsplot.for'
            stop
        endif

        sizel = 2.0

        part = .90
        ioffm = 1

        if(c_metacode .eq. 'm ')then

            if(n_image .gt. 0)then
                write(6,*)' Skipping map plot - already done with image'       

            else ! set up supmap for plot
                call lapsplot_setup(ni,nj,lat,lon,jdot_in
     1                             ,namelist_parms,plot_parms)      

            endif

            return

        endif ! 'm'

        if(c_metacode .eq. 'c ')then

            call get_border(ni,nj,x_1,x_2,y_1,y_2)
            call set(x_1,x_2,y_1,y_2,0.05,0.95,0.05,0.95,1)

            part = 0.90
            ioffm = 1

            IOFFP = 1              ! we may need this
            SPVAL = r_missing_data ! we may need this

            write(6,*)' lapsplot: call "conrec", cint = ',cint

            if(cint .ge. 0.)then
                call conrec_line(field,ni,ni,nj,clow,chigh,cint
     1                          ,plot_parms,-1,0,-1848,0)

            else ! Special Contouring
                call conrec_line(field,ni,ni,nj,0.,1e8,1e8,plot_parms
     1                          ,-1,0,-1848,0)       
                do i = 1,N_CONTOURS
                    cvalue = factor(i)
                    if( cvalue .ge. abs(cint) )then
                        write(6,*)' Contouring at +/-',i,cvalue
                        call conrec_line(field,ni,ni,nj,cvalue
     1                                  ,cvalue,1e-6,plot_parms
     1                                  ,-1,0,-1848,0)
                        call conrec_line(field,ni,ni,nj,-cvalue
     1                                  ,-cvalue,1e-6,plot_parms
     1                                  ,-1,0,-1848,0)
                    else
                        write(6,*)' Skip contouring at +/-',i,cvalue
                    endif
                enddo
            endif

        endif

        return
        end


        subroutine lapsplot_setup(ni,nj,lat,lon,jdot_in
     1                           ,namelist_parms,plot_parms)       

!       implicit none

!       include 'lapsparms.cmn'

        include 'lapsplot.inc'

        common /supmp1/ dummy,part
        common /supmp6/ umin,umax,vmin,vmax
        real dummy(8),part
        common /ERROR/ IFRAME, IERRR
        common /image/ n_image

        integer       idummy(6),ioffm,istatus

        integer ni,nj,i

        real lat(ni,nj),lon(ni,nj)

        real
     1  umin,umax,vmin,vmax

        integer N_CONTOURS
        parameter (N_CONTOURS = 20)
        real factor(N_CONTOURS)
        data factor/
     1  .01,
     1  .02,
     1  .05,
     1  .1,
     1  .2,
     1  .5,
     1  1.,
     1  2.,
     1  5.,
     1  10.,
     1  20.,
     1  50.,
     1  100.,
     1  200.,
     1  500.,
     1  1000.,
     1  2000.,
     1  5000.,
     1  10000.,
     1  20000.
     1                  /


        real cint,cvalue,clow,chigh

        integer
     1  jlts,
     1  jnj,
     1  iusout,
     1  jdot,
     1  ier,iframe

        real
     1  tx,ty,
     1  polon,
     1  rot,
     1  plm1(2),plm2(2),plm3(2),plm4(2)

        character*2 c_metacode
        character*6 c6_maproj

        write(6,*)' lapsplot_setup: start'

        sizel = 2.0

        part = .90
        ioffm = 1

c       set up supmap for plot

        call get_c6_maproj(c6_maproj,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            stop 
        endif
        write(6,*)' c6_maproj = ',c6_maproj

        call get_standard_longitude(std_lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Bad istatus in lapsplot.for'
            stop
        endif
        write(6,*)' standard_lon = ',std_lon

        call get_standard_latitudes(std_lat1,std_lat2,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error calling laps routine'
            stop 
        endif
        write(6,*)' standard_lats = ',std_lat1,std_lat2

        rot=0
        jlts=-3
        jnj=0

!       write(6,*)' UMIN/UMAX/VMIN/VMAX',umin,umax,vmin,vmax

!       mpcol1 = 4 ! For Tektronix

        sizel = 2.0

        part = .90
        ioffm = 1

c       set up supmap for plot

        rot=0
        jlts=-3
        jnj=1000


        iusout = 3

        PLM1(1)=lat(1,nj)
        PLM1(2)=lon(1,nj)
        PLM2(1)=lat(ni,1)
        PLM2(2)=lon(ni,1)
        PLM3(1)=lat(1,1)
        PLM3(2)=lon(1,1)
        PLM4(1)=lat(ni,nj)
        PLM4(2)=lon(ni,nj)

        TX=0.5
        TY=0.9765
        IFRAME=IFRAME + 1

        write(6,*)' lapsplot_setup: IFRAME = ',IFRAME,' JDOT = '
     1                                   ,jdot_in,' ',c6_maproj

        call MAPINT

!       call GSCR(1,0,0.9,1.0,0.9)
!       call GSCR(1,1,0.,0.,0.)

        call MAPSTC('OU','PS')
        call MAPSTI('C6',5)
        call MAPSTI('DO',0)

        if(c6_maproj .eq. 'plrstr')then
            polat = std_lat2
            polon = std_lon
            rrot  = 0.
            call maproj('ST',polat,polon,rrot)
            jproj = 1

        elseif(c6_maproj .eq. 'lambrt')then
            jproj = 3
            polat = std_lat1
            polon = std_lon
            rrot  = std_lat2
            call maproj('LC',polat,polon,rrot)   

        elseif(c6_maproj .eq. 'merctr')then
            jproj = 9
            polat = 0.
            polon = std_lon
            rrot  = 0.
            call maproj('ME',polat,polon,rrot)

        elseif(c6_maproj .eq. 'latlon')then
            jproj = 9
            polat = 0.
            polon = std_lon
            rrot  = 0.
            call maproj('CE',polat,polon,rrot)

        else
            write(6,*)' lapsplot_setup: Error, maproj = ',c6_maproj

        endif

        call MAPINT

        map_mode = 2

!       Set up the colors, draw the county map
        if(n_image .gt. 0)then ! image present
            icol_sta = 7       ! yellow
            icol_cou = 17      ! lavender
            jdot = 0
        else                   ! no image present
            icol_sta = 7       ! yellow
            icol_cou = 34      ! grey
            jdot = 0
        endif

        call MAPSTI('C6',icol_sta)

        call MAPSET('PO',PLM1,PLM2,PLM3,PLM4)

        call MAPINT

        call draw_county_map(PLM3,PLM4,jproj,polat,polon,rrot,jdot
     1                      ,icol_sta,icol_cou,ni,nj,namelist_parms)

!       Set up colors, draw the state map?
 
!       if(map_mode .eq. 1)then
!           call setusv_dum(2HIN,32)
!       elseif(map_mode .eq. 2)then
!           call setusv_dum(2HIN,7)
!       endif

        call MAPINT
!       if(IFRAME .eq. 1)call MAPLOT

        IF(NERRO(IERR) .ne. 0)THEN
            call EPRIN
            call ERROF
        ENDIF

        call GSELNT(0)

        write(6,*)' UMIN/UMAX/VMIN/VMAX',umin,umax,vmin,vmax

        write(6,*)' lapsplot_setup: return'

        return
        end

