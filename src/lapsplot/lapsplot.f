cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
        subroutine lapsplot(field,ni,nj,clow,chigh,cint,lat,lon
     1                  ,c_metacode,c_file,c_domain,jdot_in)

!       implicit none

        common /supmp1/ dummy,part
        common /supmp6/ umin,umax,vmin,vmax
!       common /mapcol/ mpcol1,mpcol2,mpcol3,mpcol4
        real*4 dummy(8),part
!       common /CONRE1/IOFFP,SPVAL,EPSVAL,CNTMIN,CNTMAX,CNTINT,IOFFM
!       We need to set SIZEL
!       COMMON/LABS/IA(2),NC,NREP,NCRT,ILAB,NULBLL,SIZEL,SIZEM,SIZEP
        common /ERROR/ IFRAME, IERRR

        character*(*)   c_file
        character*(*)   c_domain

        integer*4       idummy(6),ioffm,istatus

        integer*4 ni,nj,i

        real*4 field(ni,nj),lat(ni,nj),lon(ni,nj)

        real*4
     1  umin,umax,vmin,vmax

        integer*4 N_CONTOURS
        parameter (N_CONTOURS = 20)
        real*4 factor(N_CONTOURS)
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


        real*4 cint,cbase,cvalue,clow,chigh

        integer*4
     1  jlts,
     1  jnj,
     1  iusout,
     1  jdot,
     1  ier,iframe

        real*4
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

c           set up supmap for plot

            call lapsplot_setup(ni,nj,lat,lon,jdot_in)

            return

        endif ! 'm'

        if(c_metacode .eq. 'c ')then

            call get_border(ni,nj,x_1,x_2,y_1,y_2)
            call set(x_1,x_2,y_1,y_2,0.05,0.95,0.05,0.95,1)

            part = 0.90
            ioffm = 1

            IOFFP = 1              ! we may need this
            SPVAL = r_missing_data ! we may need this

            write(6,*)' lapsplot: call "conrec"'

            if(cint .ge. 0.)then
                call conrec_line
     1  (field,ni,ni,nj,clow,chigh,cint,-1,0,-1848,0)

            else ! Special Contouring
                call conrec_line
     1    (field,ni,ni,nj,0.,1e8,1e8,-1,0,-1848,0)
                cbase = 1e-4

                do i = 1,N_CONTOURS
                    cvalue = factor(i)
                    if(cvalue .ge. abs(cint) .and.
     1                 cvalue .le. abs(chigh))then
                        call conrec_line
     1     (field,ni,ni,nj,cvalue,cvalue,1e-6,-1,0,-1848,0)
                        call conrec_line
     1     (field,ni,ni,nj,-cvalue,-cvalue,1e-6,-1,0,-1848,0)
                    endif
                enddo
            endif

        endif

        return
        end


        subroutine lapsplot_setup(ni,nj,lat,lon,jdot_in)

!       implicit none

!       include 'lapsparms.cmn'
        common /supmp1/ dummy,part
        common /supmp6/ umin,umax,vmin,vmax
        real*4 dummy(8),part
        common /ERROR/ IFRAME, IERRR

        integer*4       idummy(6),ioffm,istatus

        integer*4 ni,nj,i

        real*4 lat(ni,nj),lon(ni,nj)

        real*4
     1  umin,umax,vmin,vmax

        integer*4 N_CONTOURS
        parameter (N_CONTOURS = 20)
        real*4 factor(N_CONTOURS)
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


        real*4 cint,cbase,cvalue,clow,chigh

        integer*4
     1  jlts,
     1  jnj,
     1  iusout,
     1  jdot,
     1  ier,iframe

        real*4
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

!       call MAPINT

!       call GSCR(1,0,0.9,1.0,0.9)
!       call GSCR(1,1,0.,0.,0.)

        call MAPSTC('OU','US')

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

        else
            write(6,*)' lapsplot_setup: Error, maproj = ',c6_maproj

        endif


        map_mode = 2

!       Set up the colors, draw the county map
        if(map_mode .eq. 1)then
            icol_sta = 32
            icol_cou = 32
            jdot = jdot_in
        elseif(map_mode .eq. 2)then
            icol_sta = 7
            icol_cou = 34
            jdot = 0
        endif

        call draw_county_map(PLM3,PLM4,jproj,polat,polon,rrot,jdot
     1                      ,icol_sta,icol_cou,ni,nj)

!       Set up colors, draw the state map?
 
!       if(map_mode .eq. 1)then
!           call setusv_dum(2HIN,32)
!       elseif(map_mode .eq. 2)then
!           call setusv_dum(2HIN,7)
!       endif

        call MAPSET('PO',PLM1,PLM2,PLM3,PLM4)

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

