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

        subroutine xsect(c_display,i4time_ref,lun,l_atms
     1                  ,standard_longitude,NX_L,NY_L,NZ_L,NX_C,NZ_C
     1                  ,NX_P,r_missing_data,laps_cycle_time,maxstns)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, NZ_L as dummy arguments
!       97-Aug-14     Ken Dritz     Added NX_C, NZ_C as dummy arguments
!       97-Aug-14     Ken Dritz     Removed PARAMETER declarations for
!                                   NX_C, NZ_C (commented them out)
!       97-Aug-14     Ken Dritz     Added r_missing_data, laps_cycle_time
!                                   as dummy arguments
!       97-Aug-14     Ken Dritz     Added maxstns as dummy argument
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, NZ_L, NX_C, NZ_C
!                                   to interp_3d
!       97-Aug-14     Ken Dritz     Pass r_missing_data to interp_3d
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, NZ_L, NX_C, NZ_C
!                                   to interp_3dn
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, NZ_L, NX_C, NZ_C
!                                   to interp_3d_spread
!       97-Aug-14     Ken Dritz     Pass r_missing_data to interp_3d_spread
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, NX_C, r_missing_data
!                                   to interp_2d
!       97-Aug-14     Ken Dritz     Pass maxstns to label_other_stations
!       97-Aug-25     Steve Albers  Removed /lapsplot_cmn1/
!       97-Aug-25     Steve Albers  Removed /lapsplot_cmn2/
!       97-Aug-25     Steve Albers  Removed equivalence of pcp_type_3d,rh_3d.

        real*4 lat(NX_L,NY_L),lon(NX_L,NY_L),topo(NX_L,NY_L)
        real*4 rlaps_land_frac(NX_L,NY_L)
        real*4 dx(NX_L,NY_L)
        real*4 dy(NX_L,NY_L)

        integer*4 NX_C,NZ_C,NZ_B
!       parameter (NX_C = 61)   ! NX_L ! Number of horizontal points in X-Sect
!       parameter (NZ_C = NZ_L) ! Number of vertical levels in LAPS

        parameter (NZ_B = 5)    ! Bottom Level of ATMS X-Sect

        include 'laps_cloud.inc'

        real*4 clouds_3d(NX_L,NY_L,KCLOUD)

        common/lapsplot_omega/l_convert

        logical l_sta,l_convert,lapsplot_pregen,l_atms,l_pregen 
        logical l_radar_read, l_wind_read,l_arrival_gate
        logical iflag_mvd,iflag_icing_index,iflag_cloud_type
        logical iflag_snow_potential,iflag_bogus_w

        data lapsplot_pregen /.true./

        real*4 cld_pres(KCLOUD)

        character*2 c2_cloud_type,c2_cloud_types(0:10)
        data c2_cloud_types
     1  /'  ','St','Sc','Cu','Ns','Ac','As','Cs','Ci','Cc','Cb'/

        character*2 c2_precip_type,c2_precip_types(0:10)
        data c2_precip_types
     1  /'  ','Rn','Sn','Zr','Sl','Ha','  ','  ','  ','  ','  '/

        character*1 c1_precip_types(0:10)
        data c1_precip_types
     1       /' ','R','*','Z','I','H',' ',' ',' ',' ',' '/


        character*1 c_prodtype
        character*4 fcst_hhmm

        character*4 c4_log

        real*4 dum1_array(NX_L,1)
        real*4 dum2_array(NX_L,1)
        real*4 dum3_array(NX_L,1)
        real*4 dum4_array(NX_L,1)

        data mode_lwc/2/

!       character*255 c_filespec_wd/'*.lw3'/
!       character*255 c_filespec_wc/'*.lco'/
!       character*255 c_filespec_wb/'*.lba'/
        character*255 c_filespec_qg
        character*255 c_filespec

        data c_filespec_qg/'USER_DATA:*.lqo'/

        character*31 ext_wind

        character*3 var_2d
        character*150  directory
        character*31  ext
        character*10  units_2d
        character*125 comment_2d

      ! Used for "Potential" Precip Type
        logical l_mask_pcptype(NX_C,1)
        integer*4 ibase_array(NX_C,1)
        integer*4 itop_array(NX_C,1)

        real*4 u_3d(NX_L,NY_L,NZ_L)
        real*4 v_3d(NX_L,NY_L,NZ_L)
        real*4 field_3d(NX_L,NY_L,NZ_L)
        real*4 temp_3d(NX_L,NY_L,NZ_L)
        real*4 rh_3d(NX_L,NY_L,NZ_L)
        real*4 q_3d(NX_L,NY_L,NZ_L)
        real*4 slwc_3d(NX_L,NY_L,NZ_L)
        real*4 cice_3d(NX_L,NY_L,NZ_L)
        real*4 grid_ra_ref(NX_L,NY_L,NZ_L)
        real*4 grid_ra_vel(NX_L,NY_L,NZ_L)
!       real*4 grid_ra_rfill(NX_L,NY_L,NZ_L)

        real*4 pcp_type_3d(NX_L,NY_L,NZ_L)
!       equivalence(pcp_type_3d,rh_3d)

!       real*4 pcp_type_2d(NX_C,NZ_C)
!       equivalence(pcp_type_2d,rh_2d)

        real*4 field_2d(NX_C,NZ_C)

!       COMMON/LABS/IA(2),NC,NREP,NCRT,ILAB,NULBLL,SIZEL,SIZEM,SIZEP
        character c9_radarage*9

        real*4 lifted(NX_L,NY_L)

!       real*4 cloud_cvr_2d(NX_L,NY_L)
        real*4 cloud_ceil_2d(NX_L,NY_L)
        real*4 vis_2d(NX_L,NY_L)
        real*4 cloud_top_2d(NX_L,NY_L)

        real*4 clouds_vert(NX_C,KCLOUD)

!       real*4 cloud_cvr_1d(NX_C)
        real*4 cloud_ceil_1d(NX_C)
!       real*4 cloud_low_1d(NX_C)
        real*4 cloud_top_1d(NX_C)

        common /MCOLOR/mini,maxi

        real xcoord(NX_C),ycoord(NX_C)

!       COMMON /CONRE1/IOFFP,SPVAL,EPSVAL,CNTMIN,CNTMAX,CNTINT,IOFFM

        include 'icolors.inc'

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

        real*4 lat_1d(NX_C)
        real*4 lon_1d(NX_C)
        real*4 snow_1d(NX_C)

        real*4 pres_3d(NX_L,NY_L,NZ_L)
        real*4 pres_2d(NX_L,NY_L)
        real*4 pres_1d(NX_C)

        real*4 u_vert(NX_C,NZ_C)
        real*4 v_vert(NX_C,NZ_C)
        real*4 pres_vert(NX_C,NZ_C)
        real*4 temp_2d(NX_C,NZ_C)
        real*4 rh_2d(NX_C,NZ_C)
        real*4 heights_2d(NX_C,NZ_C)
        real*4 radar_2d(NX_C,NZ_C)
        real*4 slwc_2d(NX_C,NZ_C)
        real*4 cice_2d(NX_C,NZ_C)
        real*4 field_vert(NX_C,NZ_C)
        real*4 field_vert_buf(NX_C,NZ_C)
        real*4 field_vert_diff(NX_C,NZ_C)
        real*4 field_vert2(NX_C,NZ_C)
        real*4 field_vert3(NX_P,NX_P)
        real*4 w_2d(NX_C,NZ_C)
        character cldpcp_type_2d(NX_C,NZ_C)
        real*4 mvd_2d(NX_C,NZ_C)
        integer icing_index_2d(NX_C,NZ_C)
        real*4 terrain_vert(NX_C,NZ_C)
        real*4 terrain_vert1d(NX_C)
        real*4 lon_vert(NX_C)

        integer*4 iarg

        real*4 mspkt 
        data mspkt /.518/

        character*33 c33_label
        character*1 c_display
        character*1 c1_string
        character*2 c_metacode,c_wind
        character*3 c3_string,c3_sta,c3_type
        character*3 c3_ylow,c3_xlow
        character*5 c5_arrival_gate
        character*4 c4_string,c_field
        character*7 c7_string
        character*24 asc_tim_24
        character*9 a9time,c9_string
        character*20 c20_sta
        integer*4 ity

        integer*4 N_STATIONS
        parameter (N_STATIONS = 32)

        character*3 c3_sta_array(N_STATIONS)
        data c3_sta_array
     1  /'WIG','FTC','LOV','ELB','FLG','PTV','STP',
     1         'ELB','BJC','DEN','APA','COS','CYS','LAR',
     1         'LIC','AKO','GLD','LHX','BOU','KIO','GXY',
     1   'ERI','MHR','CP3','CHL','UND','OKC','MKC',
     1   'ICT','DSM','GRI','ASE'/

        real*4 sta_lat(N_STATIONS)
        data sta_lat
     1    /  40.29,  40.59,  40.59,  39.20,  39.36,  40.26,  39.75,
     1       39.23,  39.90,  39.75,  39.57,  38.82,  41.15,  41.32,
     1       39.18,  40.17,  39.37,  38.05,  40.01,  39.35,  40.42,
     1       40.10,  39.87,  39.95,  40.44,  40.10,  35.40,  39.12,
     1       37.65,  41.53,  40.97,  39.22/

        real*4 sta_lon(N_STATIONS)
        data sta_lon
     1    /-103.05,-105.14,-105.14,-104.50,-103.04,-104.87,-104.87,
     1     -104.63,-105.12,-104.87,-104.85,-104.72,-104.82,-105.68,
     1     -103.70,-103.22,-101.70,-103.52,-105.25,-104.42,-104.63,
     1     -105.03,-104.76,-105.19,-104.64,-104.34, -97.60, -94.60,
     1      -97.43, -93.65, -98.32,-106.87/

!       sizem = 1.0
        sizel = 2.0

        vxmin = .10
        vxmax = .90
        vymin = .20
        vymax = .80

        vymin2 = .50 - (.50-vymin) * 1.25
        vymax2 = .50 + (vymax-.50) * 1.25

!       vymin3 = .50 - (.50-vymin) * 1.20
!       vymax3 = .50 + (vymax-.50) * 1.20

        if(vymin .eq. .10)then
!           iyl_remap = 13
!           iyh_remap = 8
            iyl_remap = 10
            iyh_remap = 10
        elseif(vymin .eq. .20)then
            ixl_remap = nint(float(NX_P) * .0580)
            ixh_remap = nint(float(NX_P) * .0575)
!           iyl_remap = nint(float(NX_P) * .174)
!           iyh_remap = nint(float(NX_P) * .162)
            iyl_remap = nint(float(NX_P) * .168)
            iyh_remap = nint(float(NX_P) * .168)
        else
            write(6,*)' Error, invalid vymin ',vymin
        endif

        ioffm = 1 ! Don't plot label stuff in conrec

        lapsplot_pregen = .true.

c read in laps lat/lon and topo
        call get_laps_domain_95(NX_L,NY_L,lat,lon,topo,rlaps_land_frac
     1                      ,grid_spacing_dum,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            return
        endif

        i_graphics_overlay = 0
        i_label_overlay = 0
        i_map = 0
        i_initialize = 0

        l_wind_read = .false.
        l_radar_read = .false.

!       if(c_display .eq. 't')call setusv_dum(2hIN,16)

        rleft = 1
        right = NX_C

!       Decide whether the bottom of the X-Sect is at the bottom of the LAPS domain

!       if(l_atms)then
            topo_min = 1e10
            do j = 1,NY_L
            do i = 1,NX_L
                topo_min = min(topo_min,topo(i,j))
            enddo
            enddo

            ibottom_terrain = height_to_zcoord(topo_min,istatus)

            ibottom_terrain = max(ibottom_terrain,1)

!           ibottom_terrain = 1

            write(6,*)'    Lowest Displayed Level = ',ibottom_terrain

            bottom = ibottom_terrain
!       else
!           bottom = 1
!       endif

        ibottom = bottom

        top = NZ_C
        width = right - rleft
        r_height = top - bottom

!       This lets up plot outside the main box
!       call set(.00, 1.0, .00, 1.0, rleft - width/8., right + width/8.,
!       1                            bottom - r_height/8., top + r_height/8.,1)

        if(i_initialize .eq. 1)goto100

        i_initialize = 1

!       Define Segment for Cross Section on LAPS Grid
80      continue

        if(.true.)then
            write(6,102)
102         format(/
     1  '    Type of Xsect ',
     1  ' [we, sn, xxx (azimuth-true), arr (arrival gate)]   ? ',$)
        else
            write(6,103)
103         format(/
     1 '    Type of Xsect ',
     1  ' [we, sn, xxx (azimuth-true)]                       ? ',$)
        endif

        if(l_atms)write(6,*)' Reading X-sect type from lun = ',lun

        read(lun,1201)c3_type
1201    format(a3)

        if(l_atms)write(6,*)' Just read X-sect type  = ',c3_type,lun

        l_sta = .true.
        l_arrival_gate = .false.

        if(c3_type(1:2) .eq. 'we')then
            xlow = 1.
            xhigh = NX_L

            if(.true.)then
                write(6,111)NY_L,NY_L/2+1,NY_L
111             format(/'     N-S Location ',
     1        '[1 to ',i3,'; 1 = S Edge, '
     1        ,i3,' = Center, '
     1        ,i3,' = N Edge]   OR '//
     1        6x,' CLASS:       wig,ftc,lov,elb,flg'/
     1        6x,' RADIOMETERS: ptv,stp,elb,eri'/
     1        6x,' RADARS:      mhr,cp3,chl,und'/
     1        6x,' SAOs:        bjc,den,apa,cos,cys,lar,lic,ako,gld,lhx,
     1gxy'/
     1        6x,' VORs:        kiw'/
     1        ' ',6x,'MESONET:     ',
     1        'bou                                               ? ',$)

            else ! StormFest
                write(6,1110)NY_L,NY_L/2+1,NY_L
1110            format(/'     N-S Location ',
     1        '[1 to ',i3,'; 1 = S Edge, '
     1        ,i3,' = Center, '
     1        ,i3,' = N Edge]   OR '//
     1        '$',6x,'SAOs   :     ',
     1        'okc,mkc,ict,dsm,gri                               ? ')

            endif

            read(lun,120)c3_ylow
120         format(a3)

            call upcase(c3_ylow,c3_ylow)

            l_sta = .false.

            do i = 1,N_STATIONS
              if(c3_ylow .eq. c3_sta_array(i))then
                call latlon_to_rlapsgrid(sta_lat(i),sta_lon(i),lat,lon
     1                          ,NX_L,NY_L,xsta,ysta,istatus)

                if(xsta .lt. 1 .or. xsta .gt. NX_L .OR.
     1             ysta .lt. 1 .or. ysta .gt. NY_L)then
                    write(6,*)' Station is outside domain - try again...
     1'
                    goto80
                endif

                ylow = ysta
                l_sta = .true.
                i_sta = i

                pos_sta = 1. + (NX_C-1.) * (xsta-xlow)/(xhigh-xlow)

!               pos_sta = xsta

                c3_sta = c3_ylow
              endif

            enddo

            if(.not. l_sta)then
                read(c3_ylow,*,err=80)ylow
                write(6,*)
                write(6,*)'      J = ',nint(ylow)

                ylow = max(min(nint(ylow),NY_L),1)

!75              if(ylow .lt. 1 .or. ylow .gt. NY_L)then
!                    write(6,*)' Grid point is outside domain - try again
!     1...'
!                    goto80
!                endif

            else
                write(6,72)sta_lat(i_sta),sta_lon(i_sta)
     1                          ,nint(ylow),nint(pos_sta)
72              format(/7x,'lat/lon ',2f8.2,' J/I =',2i4)
            endif

            yhigh = ylow

        elseif(c3_type(1:2) .eq. 'sn')then
            ylow = 1.
            yhigh = NY_L

            if(.true.)then
                write(6,112)NX_L,NX_L/2+1,NX_L
112             format(/'     E-W Location ',
     1        '[1 to ',i3,'; 1 = W Edge, '
     1        ,i3,' = Center, '
     1        ,i3,' = E Edge]   OR '//
     1        6x,' CLASS:       wig,ftc,lov,elb,flg'/
     1        6x,' RADIOMETERS: ptv,stp,elb,eri'/
     1        6x,' RADARS:      mhr,cp3,chl,und'/
     1        6x,' SAOs:        bjc,den,apa,cos,cys,lar,lic,ako,gld,lhx,
     1gxy'/
     1        6x,' VORs:        kiw'/
     1        ' ',6x,'MESONET:     ',
     1        'bou                                               ? ',$)

            else ! StormFest
                write(6,1120)NX_L,NX_L/2+1,NX_L
1120            format(/'     E-W Location ',
     1        '[1 to ',i3,'; 1 = W Edge, '
     1        ,i3,' = Center, '
     1        ,i3,' = E Edge]   OR '//
     1        '$',6x,'SAOs   :     ',
     1        'okc,mkc,ict,dsm,gri                               ? ')

            endif

            read(lun,120)c3_xlow

            call upcase(c3_xlow,c3_xlow)

            l_sta = .false.

            do i = 1,N_STATIONS
              if(c3_xlow .eq. c3_sta_array(i))then
!               call latlon_grid(sta_lat(i),sta_lon(i),igrid,jgrid)
                call latlon_to_rlapsgrid(sta_lat(i),sta_lon(i),lat,lon
     1                          ,NX_L,NY_L,xsta,ysta,istatus)
                if(xsta .lt. 1 .or. xsta .gt. NX_L .OR.
     1           ysta .lt. 1 .or. ysta .gt. NY_L)then
                    write(6,*)' Station is outside domain - try again...
     1'
                    goto80
                endif

                xlow = xsta
                l_sta = .true.
                i_sta = i

                pos_sta = 1. + (NX_C-1.) * (ysta-ylow)/(yhigh-ylow)

!               pos_sta = ysta

                c3_sta = c3_xlow
              endif

            enddo

            if(.not. l_sta)then
                read(c3_xlow,*,err=80)xlow
                write(6,*)
                write(6,*)'      I = ',nint(xlow)

                xlow = max(min(nint(xlow),NX_L),1)

! 76             if(xlow .lt. 1 .or. xlow .gt. NX_L)then
!                    write(6,*)' Grid point is outside domain - try again
!     1...'
!                    goto80
!                endif

            else
                write(6,82)sta_lat(i_sta),sta_lon(i_sta)
     1                          ,nint(xlow),nint(pos_sta)
82              format(/7x,'lat/lon ',2f8.2,' I/J =',2i4)
            endif

            xhigh = xlow


        else ! Try to get an azimuth for the X-Sect
            read(c3_type,*,err=85)azi_xsect

            if(.true.)then
                write(6,113)
113             format(/'     Waypoint for X-Sect: '/
     1        6x,' CLASS:       wig,ftc,lov,elb,flg'/
     1        6x,' RADIOMETERS: ptv,stp,elb,eri'/
     1        6x,' RADARS:      mhr,cp3,chl,und'/
     1        6x,' SAOs:        bjc,den,apa,cos,cys,lar,lic,ako,gld,lhx,
     1gxy'/
     1        6x,' VORs:        kiw'/
     1        6x,' MESONET:     bou'/
     1        ' ',5x,'                     OR I,J location:',27x,'? ',$)       

            else ! Storm Fest
                write(6,1130)
1130             format(/'     Waypoint for X-Sect: '/
     1        6x,' SAOs:        okc,mkc,ict,dsm,gri'/
     1        '$',5x,'                     OR I,J location:',27x,'? ')

            endif

            read(lun,130)c20_sta
130         format(a20)
            c3_sta = c20_sta(1:3)
            call upcase(c3_sta,c3_sta)

            l_sta = .false.

            do i = 1,N_STATIONS
              if(c3_sta .eq. c3_sta_array(i))then
                call latlon_to_rlapsgrid(sta_lat(i),sta_lon(i),lat,lon
     1                          ,NX_L,NY_L,xsta,ysta,istatus)
                l_sta = .true.
!               i_sta = i
              endif
            enddo

            if(.not. l_sta)then ! Get I,J of Waypoint
                read(c20_sta,*,err=85)xsta,ysta
                write(6,86)xsta,ysta
86              format(/2x,'Waypt x,y = ',2f7.2)

            else
                write(6,87)c3_sta,xsta,ysta
87              format(/2x,a3,' x,y =   ',2f7.2)

            endif

!           Calculate endpoints of X-Sect from Waypoint and Azimuth
            call xsect_endpoints
     1  (xsta,ysta,azi_xsect,xlow,ylow,xhigh,yhigh,pos_sta,istatus,
     1   NX_L,NY_L,NX_C)

            goto90

85          write(6,*)' Try again...'
            goto80

90      endif ! Type of X-Sect


100    write(6,95)
95     format(
     1  /'  Field:  [WIND: di,sp,u,v,om,vo,va,vc (barbs)'
     1  /
     1  /'           TEMP: [t,pt,tb,pb] (T, Theta, T Blnc, Theta Blnc)'       
     1  /
     1  /'           HUMID: sh,rh,rl (Specific/Relative Humidity)'
     1  /
     1  /'           ts (Thetae Sat), tw (wetbulb)'
     1  /
     1  /'           cg/cf (3D Cloud Image),  tc (Cloud Type),  '
     1  ,'tp (Precip Type)'
!       1 /'           la (LWC - Adiabatic),         lj (LWC - Adjusted)'
!       1 /'                                         sj (SLWC - Adjusted)'
     1  /'           ls (cloud liquid)' 
                                          ! ,        ss (SLWC - Smith-Feddes)'
     1  /'           ci (cloud ice)'
     1  /
     1  /'           ic (icing index)    pc (precip conc)    mv (Mean Vo
     1l Diam)'
     1  /
     1  /'           cv (cloud cover contours)'
     1  /'           rf (reflectivity-graphic), ri (ref-image)'
     1  //'     Difference field: [dif-i] '
     1  /' ',49x,'q (quit/display)]   ? ',$)

        NULBLL = 3 ! for conrec (number of lines between labels)

        read(lun,1301)c_field
1301    format(a2)

        istatus = 1

        if(c_field .eq. 'q ')goto9999

        c4_log = 'x '//c_field
        if(lun .eq. 5)call logit(c4_log)

        write(6,*)' Generating Cross Section'

        i_image = 0

        call interp_2d(lat,lat_1d,xlow,xhigh,ylow,yhigh,
     1                 NX_L,NY_L,NX_C,r_missing_data)
        call interp_2d(lon,lon_1d,xlow,xhigh,ylow,yhigh,
     1                 NX_L,NY_L,NX_C,r_missing_data)

        call interp_2d(topo,terrain_vert1d,xlow,xhigh,ylow,yhigh,
     1                 NX_L,NY_L,NX_C,r_missing_data)

        if(c_field(1:3) .eq. 'dif')then
            write(6,*)' Plotting difference field of last two entries'       
            call diff(field_vert,field_vert_buf,field_vert_diff
     1               ,NX_C,NZ_C)       

            c33_label = 'difference field'

            scale = 1.

            if(c_field(4:4) .ne. 'i')then ! contour plot
                call contour_settings(field_vert_diff,NX_C,NZ_C
     1                               ,clow,chigh,cint,zoom,scale)

                call plot_cont(field_vert_diff,scale,clow,chigh,cint,
     1               asc9_tim_3dw,       
     1               c33_label,i_overlay,c_display,lat,lon,jdot,
     1               NX_C,NZ_C,r_missing_data,laps_cycle_time)

            else ! image plot
                call array_range(field_vert_diff,NX_C,NZ_C,rmin,rmax
     1                          ,r_missing_data)

                call ccpfil(field_vert_diff,NX_C,NZ_C,rmin,rmax,'hues'
     1                     ,n_image)    
                call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                call setusv_dum(2hIN,7)
                call write_label_lplot(NX_C,NZ_C,c33_label,asc9_tim_t
     1                                          ,i_overlay,'xsect')
                call lapsplot_setup(NX_C,NZ_C,lat,lon,jdot)

            endif

        else
            if(igrid .eq. 1)then
                write(6,*)' Copying field to field_buf'
                call move(field_vert,field_vert_buf,NX_C,NZ_C)       
            endif
        endif

        igrid = 1

        if(    c_field .eq. 'di' .or. c_field .eq. 'sp'
     1    .or. c_field .eq. 'u ' .or. c_field .eq. 'v '
     1    .or. c_field .eq. 'w ' .or. c_field .eq. 'dv'
     1    .or. c_field .eq. 'vo' .or. c_field .eq. 'va' 
     1    .or. c_field .eq. 'vc' .or. c_field .eq. 'om' )then

            call input_product_info(i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,3                       ! I
     1                             ,c_prodtype              ! O
     1                             ,ext                     ! O
     1                             ,directory               ! O
     1                             ,a9time                  ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                             ,istatus)                ! O

            i4time_3dw = i4_valid

            if(c_field .ne. 'w ' .and. c_field .ne. 'om')then
                if(c_prodtype .eq. 'N')then
                    c_wind = 'b'
                else
                    c_wind = 'k'
                endif

            elseif(c_field .eq. 'om' .and. c_prodtype .eq. 'A')then
                write(6,105)
105             format('  Omega field: Kinematic (lw3), '
     1                ,'Cloud Bogused'
     1                  ,' [k,c]  ',7x,'? ',$)

                read(lun,1301)c_wind
                if(c_wind .eq. 'k' .or. c_wind .eq. 'K')c_wind = 'k'
                if(c_wind .eq. 'c' .or. c_wind .eq. 'C')c_wind = 'c'

            endif

            if  (c_prodtype .eq. 'N')then
                ext_wind = 'balance'
                call get_directory(ext_wind,directory,len_dir)
                c_filespec = directory(1:len_dir)//'lw3/*.lw3'

            elseif(c_wind .eq. 'c')then
                ext_wind = 'lco'
                call get_directory(ext_wind,directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext_wind(1:3)

            elseif(c_prodtype .eq. 'A')then
                ext_wind = 'lw3'
                call get_directory(ext_wind,directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext_wind(1:3)

            else ! Background or Forecast
                ext_wind = ext
!               call get_directory(ext_wind,directory,len_dir)
                call s_len(directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext_wind(1:3)

            endif


            if(.not. l_wind_read)then
                write(6,*)
                write(6,*)' Looking for 3D wind data: ',ext_wind(1:10)
     1                   ,' ',ext,c_field,c_wind

                if(c_field .ne. 'w ' .and. c_field .ne. 'om')then ! Non-omega
                    write(6,*)' Reading U/V'
                    if(c_prodtype .eq. 'N')then
                        directory = directory(1:len_dir)//'lw3'
                        ext = 'lw3'

                        var_2d = 'U3'

                        call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_3dw
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,NX_L,NY_L,NZ_L,u_3d,istatus)       

                        var_2d = 'V3'

                        call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_3dw       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,NX_L,NY_L,NZ_L,v_3d,istatus)       

                    elseif(c_prodtype .eq. 'A')then
                        call get_file_time(c_filespec,i4time_ref
     1                                               ,i4time_3dw)
                        call get_uv_3d(i4time_3dw,NX_L,NY_L,NZ_L
     1                                  ,u_3d,v_3d,ext_wind,istatus)

                    else ! Background or Forecast
                        var_2d = 'U3'
                        call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,u_3d
     1                              ,istatus)

                        var_2d = 'V3'
                        call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,v_3d
     1                              ,istatus)

                    endif

                    call make_fnam_lp(i4time_3dw,a9time,istatus)
                    write(6,*)' a9time = ',a9time

                elseif(c_field .eq. 'w ' .or. c_field .eq. 'om')then ! Omega
                    write(6,*)' Reading Omega/W ',c_wind,c_prodtype
                    call get_file_time(c_filespec,i4time_ref,i4time_3dw)

                    if(c_wind .eq. 'c')then
                        call get_w_3d(i4time_3dw,NX_L,NY_L,NZ_L
     1                                  ,field_3d,ext_wind,istatus)

                    elseif(c_prodtype .eq. 'A')then
                        call get_w_3d(i4time_3dw,NX_L,NY_L,NZ_L
     1                                  ,field_3d,ext_wind,istatus)

                    elseif(c_prodtype .eq. 'N')then
                        directory = directory(1:len_dir)//'lw3'
                        ext = 'lw3'
                        var_2d = 'OM'

                        call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_3dw
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,NX_L,NY_L,NZ_L,field_3d,istatus)       

                    elseif(c_prodtype .eq. 'B' .or. 
     1                     c_prodtype .eq. 'F')then
                        var_2d = 'OM'
                        call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istatus)

                    endif

                    call make_fnam_lp(i4time_3dw,a9time,istatus)

                endif
            endif
!           l_wind_read = .true.

            if(istatus .ne. 1)then
                write(6,*)' Error reading in wind field'
                goto100
            endif

        endif

        if(c_field .eq. 'vc')then

!           Remap from 3d grid to Vert Xsect grid
            call interp_3d(u_3d,u_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            call interp_3d(v_3d,v_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            call interp_2d(lon,lon_vert,xlow,xhigh,ylow,yhigh,
     1                 NX_L,NY_L,NX_C,r_missing_data)

            i_contour = 2

            if       (c_prodtype .eq. 'N')then
                c33_label = 'LAPS Wind (Balanced)       knots '
            elseif   (c_prodtype .eq. 'A')then
                c33_label = 'LAPS Wind (Analyzed)       knots '
            elseif   (c_prodtype .eq. 'B')then
                c33_label = 'LAPS Wind (Background)     knots '
            elseif   (c_prodtype .eq. 'F')then
                c33_label = 'LAPS Wind (Forecast)       knots '
            else
                c33_label = 'LAPS Wind (??????????)     knots '
            endif


        elseif(c_field .eq. 'om' )then

            if(ext_wind .eq. 'lco')then ! Cloud Omega

!               Take out missing data values to ensure better interpolation
                do k = 1,NZ_L
                do j = 1,NY_L
                do i = 1,NX_L
                    if(field_3d(i,j,k) .eq. r_missing_data)then
                        field_3d(i,j,k) = -1e-30
                    endif
                enddo ! i
                enddo ! j
                enddo ! k

                call interp_3d(field_3d,field_vert,xlow,xhigh,ylow
     1                        ,yhigh,NX_L,NY_L,NZ_L,NX_C,NZ_C
     1                        ,r_missing_data)

                call get_pres_3d(i4time_ref,NX_L,NY_L,NZ_L,pres_3d
     1                                                    ,istatus)
                call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow
     1                            ,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

                do k = 1,NZ_C
                do i = 1,NX_C
                    if(field_vert(i,k) .ne. r_missing_data)then
                        field_vert(i,k) = field_vert(i,k)*10. ! Pa/S to ubar/S
                    else
                        field_vert(i,k) = -1e-30
                    endif
                enddo ! i
                enddo ! k

                cint = -1.

            else ! Not LCO field
                call interp_3d(field_3d,field_vert
     1                        ,xlow,xhigh,ylow,yhigh
     1                        ,NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)       


                call get_pres_3d(i4time_ref,NX_L,NY_L,NZ_L,pres_3d
     1                                                    ,istatus)
                call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow,yhigh,       
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)


                do k = NZ_C,1,-1
                do i = 1,NX_C
                    if(field_vert(i,k) .ne. r_missing_data)then
                        field_vert(i,k) = field_vert(i,k)*10. ! Pa/S to ubar/S
!                   else
!                       field_vert(i,k) = field_vert(i,min(k+1,NZ_C))
                    endif
                enddo ! i
                enddo ! k

                cint = -1.

            endif ! LCO field

            i_contour = 1

            if       (c_prodtype .eq. 'N')then
                c33_label = 'LAPS Omega (balanced)      ubar/s'
            else   if(c_wind .eq. 'c')then
                c33_label = 'LAPS Omega (cloud)         ubar/s'
            else   if(c_prodtype .eq. 'A')then
                c33_label = 'LAPS Omega (analyzed)      ubar/s'
            else   if(c_prodtype .eq. 'B')then
                c33_label = 'LAPS  Bkgnd   Omega  '//fcst_hhmm
     1                                             //'  ubar/s'
            else   if(c_prodtype .eq. 'F')then
                c33_label = 'LAPS  FUA     Omega  '//fcst_hhmm
     1                                             //'  ubar/s'
            else
                c33_label = '                                 '
            endif

        elseif(c_field .eq. 'w ' )then

            if(ext_wind .eq. 'lco')then ! Cloud Omega

!               Take out missing data values to insure better interpolation
                do k = 1,NZ_L
                do j = 1,NY_L
                do i = 1,NX_L
                    if(field_3d(i,j,k) .eq. r_missing_data)then
                        field_3d(i,j,k) = -1e-30
                    endif
                enddo ! i
                enddo ! j
                enddo ! k

                call interp_3d(field_3d,field_vert,xlow,xhigh,ylow
     1                        ,yhigh,NX_L,NY_L,NZ_L,NX_C,NZ_C
     1                        ,r_missing_data)

                call get_pres_3d(i4time_ref,NX_L,NY_L,NZ_L,pres_3d
     1                                                    ,istatus)
                call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow
     1                            ,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

                do k = 1,NZ_C
                do i = 1,NX_C
                    if(field_vert(i,k) .ne. r_missing_data)then
                        if(l_convert)then
                            field_vert(i,k) = ! Always should be .true.
     1           omega_to_w(field_vert(i,k),pres_vert(i,k)) * 100.
                        endif
                    else
                        field_vert(i,k) = -1e-30
                    endif
                enddo ! i
                enddo ! k

                cint = -1.

            else ! Not LCO field
                call interp_3d(field_3d,field_vert
     1                        ,xlow,xhigh,ylow,yhigh
     1                        ,NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)       


                call get_pres_3d(i4time_ref,NX_L,NY_L,NZ_L,pres_3d
     1                                                    ,istatus)
                call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow,yhigh,       
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)


                do k = NZ_C,1,-1
                do i = 1,NX_C
                    if(field_vert(i,k) .ne. r_missing_data)then
!                       l_convert is .false. if old 'W' data is read in (not OM)
                        if(l_convert)field_vert(i,k) =
     1             omega_to_w(field_vert(i,k),pres_vert(i,k)) * 100.       
                    else
                        field_vert(i,k) = field_vert(i,min(k+1,NZ_C))
                    endif
                enddo ! i
                enddo ! k

                cint = -1.

            endif ! LCO field

            i_contour = 1

            if       (c_prodtype .eq. 'N')then
                c33_label = 'LAPS W (bal)   Vert X-Sect (cm/s)'
            else   if(c_wind .eq. 'c')then
                c33_label = 'LAPS W (cloud) Vert X-Sect (cm/s)'
            else ! if(c_prodtype .eq. 'A')then
                c33_label = 'LAPS W (kinem) Vert X-Sect (cm/s)'
            endif


        elseif(c_field .eq. 'u ' )then
            call interp_3d(u_3d,u_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            call interp_3d(v_3d,v_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            call interp_2d(lon,lon_vert,xlow,xhigh,ylow,yhigh,
     1                 NX_L,NY_L,NX_C,r_missing_data)
            do k = NZ_C,1,-1
            do i = 1,NX_C
                if(u_vert(i,k) .ne. r_missing_data)then
                    field_vert(i,k) = u_vert(i,k)/mspkt
                else
                    field_vert(i,k) = field_vert(i,min(k+1,NZ_C))
                endif
            enddo ! i
            enddo ! k
            clow = -100.
            chigh = +1000.
            cint = 10.
            i_contour = 1

            if    (c_prodtype .eq. 'N')then
                c33_label = 'LAPS  U    (balanced)       (kt) '
            elseif(c_prodtype .eq. 'A')then   
                c33_label = 'LAPS  U    (analyzed)       (kt) '
            elseif(c_prodtype .eq. 'B')then   
                c33_label = 'LAPS  U    (background)     (kt) '
            elseif(c_prodtype .eq. 'F')then   
                c33_label = 'LAPS  U    (forecast)       (kt) '
            endif

        elseif(c_field .eq. 'v ' )then
            call interp_3d(u_3d,u_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            call interp_3d(v_3d,v_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            call interp_2d(lon,lon_vert,xlow,xhigh,ylow,yhigh,
     1                 NX_L,NY_L,NX_C,r_missing_data)
            do k = NZ_C,1,-1
            do i = 1,NX_C
                if(v_vert(i,k) .ne. r_missing_data)then
                    field_vert(i,k) = v_vert(i,k)/mspkt
                else
                    field_vert(i,k) = field_vert(i,min(k+1,NZ_C))
                endif
            enddo ! i
            enddo ! k
            clow = -100.
            chigh = +1000.
            cint = 10.
            i_contour = 1

            if    (c_prodtype .eq. 'N')then
                c33_label = 'LAPS  V    (balanced)       (kt) '
            elseif(c_prodtype .eq. 'A')then   
                c33_label = 'LAPS  V    (analyzed)       (kt) '
            elseif(c_prodtype .eq. 'B')then   
                c33_label = 'LAPS  V    (background)     (kt) '
            elseif(c_prodtype .eq. 'F')then   
                c33_label = 'LAPS  V    (forecast)       (kt) '
            endif

        elseif(c_field .eq. 'sp' )then
            call interp_3d(u_3d,u_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            call interp_3d(v_3d,v_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            call interp_2d(lon,lon_vert,xlow,xhigh,ylow,yhigh,
     1                 NX_L,NY_L,NX_C,r_missing_data)
            do k = NZ_C,1,-1
            do i = 1,NX_C
                if(v_vert(i,k) .ne. r_missing_data)then
                    call uv_to_disp(        u_vert(i,k),
     1                              v_vert(i,k),
     1                              di_dum,
     1                              speed_ms)
                    field_vert(i,k) = speed_ms/mspkt
                else
                    field_vert(i,k) = field_vert(i,min(k+1,NZ_C))
                endif
            enddo ! i
            enddo ! k
            clow = -100.
            chigh = +1000.
            cint = 10.
            i_contour = 1

            if       (c_prodtype .eq. 'N')then
                c33_label = 'LAPS Isotachs (Balanced)   knots '
            elseif   (c_prodtype .eq. 'A')then
                c33_label = 'LAPS Isotachs (Analyzed)   knots '
            elseif   (c_prodtype .eq. 'B')then
                c33_label = 'LAPS Isotachs (Background) knots '
            elseif   (c_prodtype .eq. 'F')then
                c33_label = 'LAPS Isotachs (Forecast)   knots '
            endif

        elseif(c_field .eq. 'di' )then
            call interp_3d(u_3d,u_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            call interp_3d(v_3d,v_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            call interp_2d(lon,lon_vert,xlow,xhigh,ylow,yhigh,
     1                 NX_L,NY_L,NX_C,r_missing_data)
            do k = NZ_C,1,-1
            do i = 1,NX_C
                if(u_vert(i,k) .ne. r_missing_data)then
                    call uv_to_disp(u_vert(i,k),
     1                              v_vert(i,k),
     1                              field_vert(i,k),
     1                              speed_dum)
                else
                    field_vert(i,k) = field_vert(i,min(k+1,NZ_C))
                endif
            enddo ! i
            enddo ! k
            clow = -100.
            chigh = +1000.
            cint = 10.
            i_contour = 1
            if       (c_prodtype .eq. 'N')then
                c33_label = 'LAPS Isogons (Balanced)    knots '
            elseif   (c_prodtype .eq. 'A')then
                c33_label = 'LAPS Isogons (Analysis)    knots '
            elseif   (c_prodtype .eq. 'B')then
                c33_label = 'LAPS Isogons (Background)  knots '
            elseif   (c_prodtype .eq. 'F')then
                c33_label = 'LAPS Isogons (Forecast)    knots '
            endif

        elseif(c_field .eq. 'dv' )then
            do k = 1,NZ_L
                call divergence(u_3d(1,1,k),v_3d(1,1,k),field_3d(1,1,k)       
     1                         ,lat,lon,NX_L,NY_L,.true.,r_missing_data)       
            enddo

            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            do k = NZ_C,1,-1
              do i = 1,NX_C
                if(field_vert(i,k) .eq. r_missing_data)then
                    field_vert(i,k) = field_vert(i,min(k+1,NZ_C))
                else
                    field_vert(i,k) = field_vert(i,k) * 1e5
                endif
              enddo ! i
            enddo ! k

            clow = -100.
            chigh = +1000.
            cint = 2.

            i_contour = 1

            if       (c_prodtype .eq. 'N')then
                c33_label = 'LAPS Divergence  (Bal)   [1e-5/s]'
            elseif   (c_prodtype .eq. 'A')then
                c33_label = 'LAPS Divergence  (Anal)  [1e-5/s]'
            elseif   (c_prodtype .eq. 'B')then
                c33_label = 'LAPS Divergence  (Bkgnd) [1e-5/s]'
            elseif   (c_prodtype .eq. 'F')then
                c33_label = 'LAPS Divergence  (Fcst)  [1e-5/s]'
            endif

        elseif(c_field .eq. 'vo' )then
            do k = 1,NZ_L
                call vorticity_abs(u_3d(1,1,k),v_3d(1,1,k)
     1                            ,field_3d(1,1,k),lat,lon
     1                            ,NX_L,NY_L,.true.,r_missing_data)       
            enddo

            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            do k = NZ_C,1,-1
              do i = 1,NX_C
                if(field_vert(i,k) .eq. r_missing_data)then
                    field_vert(i,k) = field_vert(i,min(k+1,NZ_C))
                else
                    field_vert(i,k) = field_vert(i,k) * 1e5
                endif
              enddo ! i
            enddo ! k

            clow = -100.
            chigh = +1000.
            cint = 2.

            i_contour = 1

            if       (c_prodtype .eq. 'N')then
                c33_label = 'LAPS Abs Vort  (Bal)     [1e-5/s]'
            elseif   (c_prodtype .eq. 'A')then
                c33_label = 'LAPS Abs Vort  (Anal)    [1e-5/s]'
            elseif   (c_prodtype .eq. 'B')then
                c33_label = 'LAPS Abs Vort  (Bkgnd)   [1e-5/s]'
            elseif   (c_prodtype .eq. 'F')then
                c33_label = 'LAPS Abs Vort  (Fcst)    [1e-5/s]'
            endif

        elseif(c_field .eq. 'va' )then
            call get_grid_spacing_array(lat,lon,NX_L,NY_L,dx,dy)
            do k = 1,NZ_L
                call vorticity_abs(u_3d(1,1,k),v_3d(1,1,k)
     1                            ,field_2d,lat,lon
     1                            ,NX_L,NY_L,.true.,r_missing_data)       

                call cpt_advection(field_2d,u_3d(1,1,k),v_3d(1,1,k)
     1                            ,dx,dy,NX_L,NY_L
     1                            ,field_3d(1,1,k))

            enddo

            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            do k = NZ_C,1,-1
              do i = 1,NX_C
                if(field_vert(i,k) .eq. r_missing_data)then
                    field_vert(i,k) = field_vert(i,min(k+1,NZ_C))
                else
                    field_vert(i,k) = field_vert(i,k) * 1e8
                endif
              enddo ! i
            enddo ! k

            clow = -100.
            chigh = +1000.
            cint = 5.

            i_contour = 1

            if       (c_prodtype .eq. 'N')then
                c33_label = 'LAPS Vort Adv  (Bal)   [1e-8/s^2]'
            elseif   (c_prodtype .eq. 'A')then
                c33_label = 'LAPS Vort Adv  (Anal)  [1e-8/s^2]'
            elseif   (c_prodtype .eq. 'B')then
                c33_label = 'LAPS Vort Adv  (Bkgnd) [1e-8/s^2]'
            elseif   (c_prodtype .eq. 'F')then
                c33_label = 'LAPS Vort Adv  (Fcst)  [1e-8/s^2]'
            endif

        elseif(c_field .eq. 'rf' .or. c_field .eq. 'rg'
     1                           .or. c_field .eq. 'rk')then
            if(c_field .ne. 'rg')then
                i4time_get = i4time_ref/laps_cycle_time 
     1                     * laps_cycle_time
                goto1300
            endif

1300        write(6,*)' Getting Radar data via get_laps_3dgrid'
!           var_2d = 'REF'
!           ext = 'lps'

!           call get_laps_3dgrid(i4time_ref,86400,i4time_radar,
!    1          NX_L,NY_L,NZ_L,ext,var_2d
!    1                  ,units_2d,comment_2d,grid_ra_ref,istatus)

            call input_product_info(i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,3                       ! I
     1                             ,c_prodtype              ! O
     1                             ,ext                     ! O
     1                             ,directory               ! O
     1                             ,a9time                  ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                             ,istatus)                ! O

            var_2d = 'REF'

            if(c_prodtype .eq. 'A')then
                write(6,*)' Getting Radar data via get_laps_3dgrid'
                ext = 'lps'

                call get_laps_3dgrid(i4time_get,86400,i4time_radar
     1              ,NX_L,NY_L,NZ_L,ext,var_2d
     1              ,units_2d,comment_2d,grid_ra_ref,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' Could not read lps via get_laps_3dgrid'       
                    goto100
                endif

                call make_fnam_lp(i4time_radar,a9time,istatus)
                c33_label = 'LAPS  Reflectivity  Vert X-Sect  '

            elseif(c_prodtype .eq. 'F')then
                call get_lapsdata_3d(i4_initial,i4_valid,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,grid_ra_ref
     1                              ,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' Could not read forecast ref'       
                    goto100
                endif
                c33_label = 'LAPS  FUA Reflectivity '//fcst_hhmm
     1                    //'   dbz'

            else
                goto100

            endif

            call interp_3d(grid_ra_ref,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)
            clow = 0.
            chigh = +100.
            cint = 10.
            i_contour = 1
!           c33_label = 'LAPS  Reflectivity  Vert X-Sect  '
!           call make_fnam_lp(i4time_radar,a9time,istatus)

        elseif(c_field .eq. 'ri' .or. c_field .eq. 'rj' 
     1    .or. c_field .eq. 'rs')then ! Reflectivity Image
            i_image = 1
            if(c_field .eq. 'ri')then
                i4time_get = i4time_ref/laps_cycle_time 
     1                                * laps_cycle_time
            else
                i4time_get = i4time_ref
            endif

            if(c_field .eq. 'ri')goto1500 ! Skip next part

            if(.not. l_radar_read)then
!               Obtain height field
                ext = 'lt1'
                var_2d = 'HT'
                call get_laps_3dgrid(
     1                   i4time_get,10000000,i4time_ht,
     1                   NX_L,NY_L,NZ_L,ext,var_2d
     1                  ,units_2d,comment_2d,heights_3d,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' Error locating height field'
                    go to 100
                endif

                call get_radar_ref(i4time_get,2000,i4time_radar,1
     1            ,1,NX_L,NY_L,NZ_L,lat,lon,topo,.true.,.true.
     1            ,heights_3d
     1            ,grid_ra_ref,n_ref
     1            ,rlat_radar,rlon_radar,rheight_radar,istat_2dref
     1            ,istat_3dref)

                if(istat_2dref .le. 0)goto 100

                l_radar_read = .true.

            endif

            if(istat_3dref .le. 0)then
                if(istat_2dref .eq. 1)then
                    write(6,*)
     1              ' Radar Xsect unavailable, try earlier time'
                endif
                goto 100
            endif

            goto1510

1500        continue

            call input_product_info(i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,3                       ! I
     1                             ,c_prodtype              ! O
     1                             ,ext                     ! O
     1                             ,directory               ! O
     1                             ,a9time                  ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                             ,istatus)                ! O

            var_2d = 'REF'

            if(c_prodtype .eq. 'A')then
                write(6,*)' Getting Radar data via get_laps_3dgrid'
                ext = 'lps'

                call get_laps_3dgrid(i4time_get,86400,i4time_radar
     1              ,NX_L,NY_L,NZ_L,ext,var_2d
     1              ,units_2d,comment_2d,grid_ra_ref,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' Could not read lps via get_laps_3dgrid'       
                    goto100
                endif

                call make_fnam_lp(i4time_radar,a9time,istatus)
                c33_label = 'LAPS  Reflectivity  Vert X-Sect  '

            elseif(c_prodtype .eq. 'F')then
                call get_lapsdata_3d(i4_initial,i4_valid,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,grid_ra_ref
     1                              ,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' Could not read forecast ref'       
                    goto100
                endif
                c33_label = 'LAPS  FUA Reflectivity '//fcst_hhmm
     1                    //'   dbz'

            else
                goto100

            endif

1510        continue

            call get_ref_base(ref_base,istatus)

            do i=1,NX_L
            do j=1,NY_L
            do k=1,NZ_L
                if(grid_ra_ref(i,j,k) .eq. r_missing_data)then
                    grid_ra_ref(i,j,k) = ref_base
                endif  
            enddo ! k
            enddo ! j
            enddo ! i

            if(c_field .ne. 'rs')then
                call interp_3d(grid_ra_ref,field_vert
     1                        ,xlow,xhigh,ylow,yhigh
     1                        ,NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)       

            else ! Get Spread Out Data
                call interp_3d_spread
     1          (grid_ra_ref,field_vert,xlow,xhigh,ylow,yhigh,
     1           NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            endif

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)

            write(6,*)' Generating Reflectivity Image'
            do i = 1,NX_C-1
            do k = ibottom,NZ_L-1
                x1 = i
                y1 = k

!               Perform a bi-linear interpolation to provide an image
!               This image consists of a set of line segments which will look
!               smoother than blocks the size of the grid.

                Z1=field_vert(i  , k  )
                Z2=field_vert(i+1, k  )
                Z3=field_vert(i+1, k+1)
                Z4=field_vert(i  , k+1)

                nii = 15
                njj = 60

                do ii = 0,nii
                  fraci = float(ii) / float(nii)
                  xx = x1 + fraci

                  do jj = 0,njj
                    fracj = float(jj) / float(njj)
                    yy = y1 + fracj
                    r_dbz =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                  - (Z2+Z4-Z3-Z1)*fraci*fracj

                    i_dbz = nint(r_dbz)
                    i_dbz = i_dbz/5 * 5 ! This reduces the resolution and saves on graphics

!                   Plot line segments of the same color
                    if(i_dbz .ne. i_dbz_ref
     1                        .and. fracj .gt. 0.)then

                        y_last = yy - 1./float(njj)

                        if(i_dbz_ref .ge. 0)then
                            icol = 180 + i_dbz_ref / 5

!                           write(6,432)i_dbz,i_dbz_ref,xx,y_ref,y_last
 432                        format(2i5,3f9.4)

                            call setusv_dum(2hIN,icol)

                            call line(xx,y_ref,xx,y_last)
                        endif

                        i_dbz_ref = i_dbz
                        y_ref = yy

                    elseif(fracj .eq. 0.0)then
                        i_dbz_ref = i_dbz
                        y_ref = yy

                    endif

                    if(fracj .eq. 1.0)then
                        if(i_dbz .ge. 0)then
                           icol = 180 + i_dbz_ref / 5

!                          write(6,432)i_dbz,i_dbz_ref,xx,y_ref,y_last        

                           call setusv_dum(2hIN,icol)

                           call line(xx,y_ref,xx,yy)

                        endif
                    endif

                  enddo
                enddo
            enddo ! k
            enddo ! i

            i_contour = 0

        elseif(c_field .eq. 'cf' )then ! Cloud Gridded Image
            i_image = 1

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)

            call setusv_dum(2hIN,2)

            write(6,*)' Plotting cloud gridded image'

            ext = 'lcp'
            var_2d = 'LCP'
            call get_laps_3dgrid(
     1                   i4time_ref,10000000,i4time_nearest,
     1                   NX_L,NY_L,NZ_L,ext,var_2d
     1                  ,units_2d,comment_2d,field_3d,istatus)
            if(istatus .ne. 1)then
                write(6,*)' No cloud grid available'
            endif

            call make_fnam_lp(i4time_nearest,a9time,istatus)

            call interp_3d(field_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            c33_label = 'LAPS Gridded Cloud Cover   X-Sect'

            call remap_field_2d(
     1                            NX_C,1,NX_C
     1                           ,NZ_C,ibottom,NZ_C
     1                           ,NX_P, ixl_remap, NX_P-ixh_remap
     1                           ,NX_P, iyl_remap, NX_P-iyh_remap
     1                           ,field_vert,field_vert3,r_missing_data)

            write(6,*)' calling solid fill cloud plot'
            call ccpfil(field_vert3,NX_P,NX_P,0.0,1.0,'linear',n_image)       

        elseif(c_field .eq. 'cg' )then ! Cloud Gridded Image
            i_image = 1

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)

            call setusv_dum(2hIN,2)

            write(6,*)' Plotting cloud gridded image'

            ext = 'lc3'
            call get_clouds_3dgrid(i4time_ref,i4time_nearest
     1               ,NX_L,NY_L,KCLOUD
     1               ,ext,clouds_3d,cld_hts,cld_pres,istatus)

            if(istatus .ne. 1)then
                write(6,*)' No cloud grid available'
            endif

            call make_fnam_lp(i4time_nearest,a9time,istatus)

            call interp_3dc(clouds_3d,clouds_vert,xlow,xhigh,ylow,yhigh,
     1               NX_L,NY_L,NX_C,r_missing_data)

            niii = 12 ! horizontal resolution of cloud plot

            do i = 1,NX_C

              i_eighths_ref = 0
              k_ref = ibottom+1

              do k = ibottom+1,KCLOUD-1

                if(clouds_vert(i,k) .ne. r_missing_data)then

                    if(l_atms)then ! Turn into a binary (bipolar) field
                        if(clouds_vert(i,k) .gt. 0.65)then
                            clouds_vert(i,k) = 1.0
                        else
                            clouds_vert(i,k) = 0.0
                        endif
                    endif

                    i_eighths = nint(clouds_vert(i,k)*8.)

                else
                    i_eighths = 0
                endif

c               if(clouds_vert(i,k) .gt. 0.01)
c       1       write(6,1100)i,k,nint(cld_hts(k)),clouds_vert(i,k)
c       1                                               ,i_eighths
1100            format(2i3,i6,f6.2,i3)

                if(i_eighths .ne. i_eighths_ref)then

!                 Remap clouds using standard atmosphere
                  chigh = (cld_hts(k) + cld_hts(k-1))/2.
                  clow  = (cld_hts(k_ref) + cld_hts(k_ref-1))/2.

!                 Remap clouds using ambient pressure in center of domain
                  phigh = (cld_pres(k) + cld_pres(k-1))/2.
                  plow  = (cld_pres(k_ref) + cld_pres(k_ref-1))/2.

c                 write(6,1101)i_eighths_ref,nint(clow),nint(chigh)
1101              format(1x,i2,2i6)

                  if(i_eighths_ref .ge. 1)then
                    do ii = -niii,+niii

                      fraci = float(ii) / float(2*niii)

                      if((fraci+0.5) .lt. i_eighths_ref/8.0)then

!                       Remap clouds using standard atmosphere
!                       r_low = height_to_zcoord(clow,istatus)
!                       r_high = height_to_zcoord(chigh,istatus)

!                       Remap clouds using ambient pressure in center of domain
                        r_low = zcoord_of_pressure(plow)
                        r_high = zcoord_of_pressure(phigh)

                        uu = i + fraci

                        if(r_low .lt. NZ_L)then ! Keep below the top of the domain
                          call line(uu, r_low, uu, r_high)
                        endif

                      endif ! This line is covered
                    enddo ! ii
                  endif ! Finite Cloud Cover in this grid segment

                  i_eighths_ref = i_eighths
                  k_ref = k

                endif ! Change in cloud cover

              enddo ! k

            enddo ! i

!           Generate Cloud Key
            if(.not. l_atms)then
              do i_eighths = 1,8
                i = 7 + 5 * i_eighths

                do ii = -niii,+niii
                      x = i-2
                      y = bottom - r_height * .055
                      write(c3_string,2013)i_eighths
                      call pwrity (x, y, c3_string, 3, 0, 0, 0)
2013                  format(i1,'/8')

                      fraci = float(ii) / float(2*niii)

                      if((fraci+0.5) .lt. i_eighths/8.0)then
                        r_low  = ibottom - 1. + 0.05
                        r_high = ibottom - 1. + 0.45
                        uu = i + fraci

                        if(r_low .lt. NZ_L)then ! Keep below the top of the domain
                          call line(uu, r_low, uu, r_high)
                        endif

                      endif ! This line is covered
                enddo ! ii
              enddo ! i

              c33_label = 'LAPS Gridded Cloud Cover   X-Sect'

            endif ! l_atms

        elseif(c_field .eq. 'pt')then
            iflag_temp = 0 ! Returns Potential Temperature
            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,NX_L,NY_L,NZ_L,temp_3d,istatus)
!           if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            clow = 200.
            chigh = +500.
            cint = 5.
            i_contour = 1
            c33_label = 'LAPS Potl Temp Vert X-Sect    K  '

        elseif(c_field .eq. 'pb')then
            iflag_temp = 3 ! Returns Balanced Potential Temperature
            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,NX_L,NY_L,NZ_L,temp_3d,istatus)
!           if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            clow = 200.
            chigh = +500.
            cint = 5.
            i_contour = 1
            c33_label = 'LAPS Potl Temp (Balanced)     K  '

        elseif(c_field .eq. 't ')then
            call input_product_info(i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,3                       ! I
     1                             ,c_prodtype              ! O
     1                             ,ext                     ! O
     1                             ,directory               ! O
     1                             ,a9time                  ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                             ,istatus)                ! O

            if(c_prodtype .eq. 'A')then
                iflag_temp = 1 ! Returns Ambient Temp (K)

                call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                          ,NX_L,NY_L,NZ_L,temp_3d,istatus)

                call make_fnam_lp(i4time_nearest,a9time,istatus)

            else
                write(6,*)' Sorry, not yet supported'
                goto100

            endif

            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            do k = 1,NZ_L
            do i = 1,NX_C
                field_vert(i,k) = field_vert(i,k) - 273.15 ! K to C
            enddo ! i
            enddo ! k

            clow = -100.
            chigh = +100.
            cint = 5.
            i_contour = 1
            c33_label = 'LAPS Temp      Vert X-Sect  Deg C'

        elseif(c_field .eq. 'tb')then
            var_2d = 'T3'
            call make_fnam_lp(i4time_ref,a9time,istatus)
            ext='lt1'

            call get_directory('balance',directory,lend)
            directory=directory(1:lend)//'lt1/'

            call get_3dgrid_dname(directory
     1           ,i4time_ref,laps_cycle_time*10000,i4time_nearest
     1           ,ext,var_2d,units_2d
     1           ,comment_2d,NX_L,NY_L,NZ_L,temp_3d,istatus)       

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            do k = 1,NZ_L
            do i = 1,NX_C
                field_vert(i,k) = field_vert(i,k) - 273.15 ! K to C
            enddo ! i
            enddo ! k

            clow = -100.
            chigh = +100.
            cint = 5.
            i_contour = 1
            c33_label = 'LAPS Temp (Balanced)        Deg C'

        elseif(c_field .eq. 'ts')then
            iflag_temp = 1 ! Returns Ambient Temp (K)
            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,NX_L,NY_L,NZ_L,temp_3d,istatus)
!           if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)


            call get_pres_3d(i4time_nearest,NX_L,NY_L,NZ_L,pres_3d
     1                                                    ,istatus)
            call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)


            do k = 1,NZ_L
            do i = 1,NX_C
                field_vert(i,k) =
     1           OS(field_vert(i,k)-273.15,pres_vert(i,k)/100.) + 273.15       
            enddo ! i
            enddo ! k

            clow = +200.
            chigh = +600.
            cint = 5.
            i_contour = 1
            c33_label = 'LAPS Theta(e) Sat   X-Sect  Deg K'

        elseif(c_field .eq. 'tw')then
            iflag_temp = 1 ! Returns Ambient Temp (K)
            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,NX_L,NY_L,NZ_L,temp_3d,istatus)
!           if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(temp_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            var_2d = 'RHL'
            ext = 'lh3'
            call get_laps_3dgrid
     1          (i4time_nearest,1000000,i4time_nearest,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,rh_3d,istatus)
            if(istatus .ne. 1)goto100

            call interp_3d(rh_3d,field_vert2,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            call get_pres_3d(i4time_nearest,NX_L,NY_L,NZ_L,pres_3d
     1                                                    ,istatus)
            call interp_3d(pres_3d,pres_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)


            do k = 1,NZ_L
            do i = 1,NX_C
                t_c         = field_vert(i,k) - 273.15
                td_c        = DWPT(t_c,field_vert2(i,k))
                pressure_mb = pres_vert(i,k)/100.

!               This function call here is fast but returns a t_wb_c
!               equal to t_c if pres < 500mb. This approximation should
!               not hurt the algorithm.

                t_wb_c = twet_fast(t_c,td_c,pressure_mb)
                field_vert(i,k) = t_wb_c
            enddo ! i
            enddo ! k

            clow = -100.
            chigh = +100.
            cint = 5.
            i_contour = 1
            c33_label = 'LAPS Wet Bulb       X-Sect  Deg K'

        elseif(c_field .eq. 'sh')then
            var_2d = 'SH '
            ext = 'lq3'
            call get_laps_3dgrid
     1      (i4time_ref,1000000,i4time_nearest,NX_L,NY_L,NZ_L
     1          ,ext,var_2d,units_2d,comment_2d
     1                                  ,q_3d,istatus)
            if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(q_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            do k = NZ_C,1,-1
            do i = 1,NX_C
                if(field_vert(i,k) .ne. r_missing_data)then
                    field_vert(i,k) = field_vert(i,k) * 1000.
                endif
            enddo ! i
            enddo ! k

            clow = 0.
            chigh = +40.
            cint = 0.4
            cint = -1.
            i_contour = 1
            c33_label = 'LAPS Specific Humidity     x1e3  '

        elseif(c_field .eq. 'rh' .or. c_field .eq. 'rl')then
            if(c_field .eq. 'rh')then
                var_2d = 'RH3'
            elseif(c_field .eq. 'rl')then
                var_2d = 'RHL'
            endif

            write(6,*)' Reading lh3 / ',var_2d

            ext = 'lh3'
            call get_laps_3dgrid
     1      (i4time_ref,1000000,i4time_nearest,NX_L,NY_L,NZ_L
     1          ,ext,var_2d,units_2d,comment_2d
     1                                  ,rh_3d,istatus)
            if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(rh_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            do k = NZ_C,1,-1
            do i = 1,NX_C
                if(field_vert(i,k) .ne. r_missing_data)then
                    field_vert(i,k) = field_vert(i,k)
                endif
            enddo ! i
            enddo ! k

            clow = 0.
            chigh = +100.
            cint = 10.
            i_contour = 1
            c33_label = 'LAPS Relative Humidity     %     '

            NULBLL = 1 ! for conrec (number of lines between labels)

        elseif(c_field .eq. 'cv')then
            var_2d = 'LCP'
            ext = 'lcp'
            call get_laps_3dgrid
     1  (i4time_ref,1000000,i4time_nearest,NX_L,NY_L,NZ_L
     1          ,ext,var_2d,units_2d,comment_2d
     1                                  ,rh_3d,istatus)
            if(istatus .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            call interp_3d(rh_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            clow = 0.
            chigh = +1.
            cint = 0.2
            i_contour = 1
            c33_label = 'LAPS Cloud Fraction              '

            NULBLL = 1 ! for conrec (number of lines between labels)

        elseif(c_field .eq. 'la' .or. c_field .eq. 'lj'
     1                         .or. c_field .eq. 'sj'
     1                         .or. c_field .eq. 'ls'
     1                         .or. c_field .eq. 'ss'
     1                         .or. c_field .eq. 'ci'
     1                         .or. c_field .eq. 'pc'
     1                                          )then

            if(c_field .eq. 'la')then
                iflag_slwc = 1
                c33_label = 'LAPS Adiabatic LWC      g/m^3    '
            elseif(c_field .eq. 'lj')then
                iflag_slwc = 2
                c33_label = 'LAPS Adjusted  LWC      g/m^3    '
            elseif(c_field .eq. 'sj')then
                iflag_slwc = 3
                c33_label = 'LAPS Adjusted  SLWC     g/m^3    '
            elseif(c_field .eq. 'ls')then
                iflag_slwc = 13
                c33_label = 'LAPS Smith-Feddes LWC   g/m^3    '
            elseif(c_field .eq. 'ci')then
                iflag_slwc = 13
                c33_label = 'LAPS Cloud Ice          g/m^3    '
            elseif(c_field .eq. 'ss')then
                iflag_slwc = 14
                c33_label = 'LAPS Smith-Feddes SLWC  g/m^3    '
            elseif(c_field .eq. 'pc')then
                c33_label = 'LAPS Precip Concen      g/m^3    '
            endif

            i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time

            l_pregen = lapsplot_pregen

            write(6,*)' Getting pregenerated LWC file'
            if(c_field .eq. 'ls')then
                var_2d = 'LWC'
            elseif(c_field .eq. 'ci')then
                var_2d = 'ICE'
            elseif(c_field .eq. 'pc')then
                var_2d = 'PCN'
            endif
            ext = 'lwc'
            call get_laps_3dgrid(i4time_lwc,86400,i4time_nearest,
     1          NX_L,NY_L,NZ_L,ext,var_2d
     1                  ,units_2d,comment_2d,slwc_3d,istatus)

            call interp_3d(slwc_3d,field_vert,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

            do k = 1,NZ_C
            do i = 1,NX_C
                field_vert(i,k) = field_vert(i,k) * 1e3
            enddo ! i
            enddo ! k

            clow = 0.
            chigh = 0.
            if(c_field .eq. 'pc')then
                cint = -0.01
            else
                cint = -0.01
            endif
            i_contour = 1
            call make_fnam_lp(i4time_nearest,a9time,istatus)

        elseif(c_field .eq. 'ic')then

            iflag_slwc = 13
            c33_label = '    LAPS Icing Index             '

            i4time_lrp = i4time_ref/laps_cycle_time * laps_cycle_time

            if(lapsplot_pregen)then
!           if(.false.)then
                write(6,*)' Getting pregenerated LRP file'
                var_2d = 'LRP'
                ext = 'lrp'
                call get_laps_3dgrid(i4time_lrp,86400,i4time_cloud,
     1          NX_L,NY_L,NZ_L,ext,var_2d
     1                  ,units_2d,comment_2d,slwc_3d,istatus)
                call interp_3dn(slwc_3d,field_vert,xlow,xhigh,ylow,yhigh
     1                         ,NX_L,NY_L,NZ_L,NX_C,NZ_C)

                do k = 1,NZ_C
                do i = 1,NX_C
                    iarg = nint(field_vert(i,k))
                    icing_index_2d(i,k) = iarg
                enddo ! i
                enddo ! k

            else ! Calculate on the Fly

            endif ! Read Pregenerated File

            clow = 0.
            chigh = 0.
            cint = -0.1
            i_contour = 4
            call make_fnam_lp(i4time_cloud,a9time,istatus)

        elseif(c_field .eq. 'mv')then
            iflag_slwc = 0
            c33_label = 'LAPS Mean Volume Diameter  m^-6  '

            i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time

            if(.false.)then ! Calculate on the Fly

            else
                write(6,*)' Reading pregenerated MVD'
                var_2d = 'LMD'
                ext =    'lmd'
                call get_laps_3dgrid(i4time_lwc,86400,i4time_cloud,
     1          NX_L,NY_L,NZ_L,ext,var_2d
     1          ,units_2d,comment_2d,pcp_type_3d,istatus)

                call interp_3dn(pcp_type_3d,field_vert,xlow,xhigh
     1                         ,ylow,yhigh,NX_L,NY_L,NZ_L,NX_C,NZ_C)

            endif ! (Read Pregenerated File)

            do k = 1,NZ_C
            do i = 1,NX_C
                field_vert(i,k) = field_vert(i,k) * 1e6 + .01
            enddo ! i
            enddo ! k

            clow = 10.
            chigh = 26.
            cint = 2.
            i_contour = 1
            call make_fnam_lp(i4time_cloud,a9time,istatus)

        elseif(c_field .eq. 'tc')then
            iflag_slwc = 0
            c33_label = '        LAPS Cloud Type          '

            i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time

            write(6,*)' Reading pregenerated cloud type'
            var_2d = 'CTY'
            ext =    'lty'
            call get_laps_3dgrid(i4time_lwc,86400,i4time_nearest,
     1          NX_L,NY_L,NZ_L,ext,var_2d
     1          ,units_2d,comment_2d,pcp_type_3d,istatus)

            call interp_3dn(pcp_type_3d,field_2d,xlow,xhigh
     1                     ,ylow,yhigh,NX_L,NY_L,NZ_L,NX_C,NZ_C)

            i_contour = 3
            call make_fnam_lp(i4time_nearest,a9time,istatus)

        elseif(c_field .eq. 'tp')then
            iflag_slwc = 0
            c33_label = 'LAPS Precip Type                  '

            i4time_pcp = i4time_ref

            if(.true.)then

                write(6,*)' Reading pregenerated precip type'
                var_2d = 'PTY'
                ext =    'lty'
                call get_laps_3dgrid(i4time_pcp,86400,i4time_radar,
     1          NX_L,NY_L,NZ_L,ext,var_2d
     1          ,units_2d,comment_2d,pcp_type_3d,istatus)

                call interp_3dn(pcp_type_3d,field_2d,xlow,xhigh
     1                         ,ylow,yhigh,NX_L,NY_L,NZ_L,NX_C,NZ_C)

            else ! Calculate Precip Type on the Fly


            endif ! pregen

            i_contour = 5
            call make_fnam_lp(i4time_radar,a9time,istatus)

        elseif(c_field .eq. 'ic')then
            iflag_slwc = 0
            c33_label = '     LAPS Icing Severity Index   '

            i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time

            if(lapsplot_pregen)then

            else

            endif ! Get Pregenerated file

        endif ! c_field

        write(6,1605)c33_label,a9time
1605    format(2x,a33,2x,a9)

        call lib$set_symbol('DATE_LAPSPLOT',a9time)
        istatus = lib$set_logical('DATE_LAPSPLOT',a9time)

        c_metacode = 'c '

        call i4time_fname_lp(a9time,i4time_dum,istatus)
        call cv_i4tim_asc_lp(i4time_dum,asc_tim_24,istatus)
        asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

        i_label_overlay = i_label_overlay + 1

        if(i_image .eq. 0)then
            i_graphics_overlay = i_graphics_overlay + 1
            call setusv_dum(2hIN,icolors(i_graphics_overlay))
        else
            call setusv_dum(2hIN,7) ! Yellow
        endif

        write(6,*)' Plotting, Overlay = ',i_graphics_overlay
     1                                   ,i_label_overlay

        call upcase(c33_label,c33_label)
        call set(0., 1., vymin2, vymax2, 0.,1.,0.,1.,1)

!       Write bottom label
!       if(i_label_overlay .le. 1)then
!           ity = 35
!           call pwrity(cpux(320),cpux(ity),c33_label,33,1,0,0)
!           call pwrity(cpux(800),cpux(ity),asc_tim_24(1:17),17,1,0,0)
!       endif

        call write_label_lplot(100,94,c33_label,a9time,i_label_overlay
     1                        ,'xsect')       

        if(i_map .eq. 0)then

            i_map = 1

            call setusv_dum(2hIN,7)

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)

!           This lets up plot outside the main box
            call set(.00, 1.0, vymin2 , vymax2
     1             , rleft - width/8., right + width/8.,
     1               bottom - r_height/8., top + r_height/8.,1)

!           Draw box enclosing graph
            xcoord(1) = rleft
            ycoord(1) = bottom
            xcoord(2) = right
            ycoord(2) = bottom
            xcoord(3) = right
            ycoord(3) = top
            xcoord(4) = rleft
            ycoord(4) = top
            xcoord(5) = xcoord(1)
            ycoord(5) = ycoord(1)
            npts = 5
            call curve (xcoord, ycoord, npts)

!           Label Left Axis
            if(.true.)then ! Label Height (kft msl) on Left Axis
                Do i = ibottom,NZ_C
                    y = i
                    call line (rleft, y, rleft + width * .015, y )

!                   Pressure
                    x = rleft - width * .030
                    ipres_mb = nint(zcoord_of_level(i)/100.)
                    write(c4_string,2014)ipres_mb
                    call pwrity (x, y, c4_string, 4, 0, 0, 0)
2014                format(i4)
                end do
                call pwrity (rleft - .070 * width,bottom + r_height*0.5,
     1          ' PRESSURE (HPA) ',16,1,90,0)
!           endif


!           if(.not. l_atms)then
!               Label Height (km - msl) on Right Axis
                Do i = 0,16
                    y = height_to_zcoord(i*1000.,istatus)

                    if(y .ge. bottom)then
                        call line (right, y, right - width * .015, y )

!                       Height
                        x = right + width * .015
                        iht_km = i
                        write(c4_string,2014)iht_km
                        call pwrity (x, y, c4_string, 4, 0, 0, 0)
                    endif

                end do
                call pwrity (right + .070 * width,bottom + r_height*0.5,
     1          'HEIGHT  (KM MSL)',16,1,270,0)
            endif


!           Put in lat/lon of endpoints
            x = rleft-1.0
            y = bottom - r_height * .03
            write(c7_string,2017)lat_1d(1)
            call pwrity (x, y, c7_string, 7, 0, 0, 0)
            y = bottom - r_height * .05
            write(c7_string,2017)lon_1d(1)
            call pwrity (x, y, c7_string, 7, 0, 0, 0)

            x = right+1.0
            y = bottom - r_height * .03
            write(c7_string,2017)lat_1d(NX_C)
            call pwrity (x, y, c7_string, 7, 0, 0, 0)
            y = bottom - r_height * .05
            write(c7_string,2017)lon_1d(NX_C)
            call pwrity (x, y, c7_string, 7, 0, 0, 0)

2017        format(f7.2)

        endif ! i_map .eq. 0


        if(i_contour .eq. 1)then
            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, bottom, top,1)

            mini = icolors(i_graphics_overlay)
            maxi = icolors(i_graphics_overlay)

            call setusv_dum(2hIN,icolors(i_graphics_overlay)) ! Is this effective?

            vmax = -1e30
            vmin = 1e30

            do k = ibottom,NZ_C
            do i = 1,NX_C
                vmax = max(vmax,field_vert(i,k))
                vmin = min(vmin,field_vert(i,k))
            enddo ! k
            enddo ! i

            write(6,*)' CLOW,HIGH,CINT ',clow,chigh,cint
            write(6,*)' Max/Min = ',vmax,vmin

            if(cint .ge. 0.)then
                if(.false.)then
!                   CALL CPSETI ('CLC - CONTOUR LINE COLOR INDEX'
!    1                           , icolors(i_graphics_overlay))
                    call conrec
     1              (field_vert(1,ibottom),NX_C,NX_C,(NZ_C-ibottom+1)
     1                             ,clow,chigh,cint,-1,0,-1848,0)

                else ! Can we make this work (anamorphically) for color plots?
!                   call get_border(ni,nj,x_1,x_2,y_1,y_2)
!                   call set(x_1,x_2,y_1,y_2,0.05,0.95,0.05,0.95,1)

                    call get_border(NX_C,NZ_C-ibottom+1,x_1,x_2,y_1,y_2)       
                    write(6,*)' Calling Set for conrec_line'
                    call set(0.10,0.90,0.05,0.95,0.10,0.90,0.05,0.95,1)

!                   call set(.20, .80, .20, .80
!    1                     , 1., 61., 1., 21. ,1)

                    call remap_field_2d(
     1                            NX_C,1,NX_C
     1                           ,NZ_C,ibottom,NZ_C
     1                           ,NX_P, ixl_remap, NX_P-ixh_remap
     1                           ,NX_P, iyl_remap, NX_P-iyh_remap
     1                           ,field_vert,field_vert3,r_missing_data)


                    call conrec_line
!    1              (field_vert3(1,ibottom),NX_P,NX_P,NX_P
     1              (field_vert3,NX_P,NX_P,NX_P
     1                             ,clow,chigh,cint,-1,0,-1848,0)


                endif

            else ! logarithmic plot
              if(.false.)then
                call conrec(field_vert(1,ibottom)
     1                     ,NX_C,NX_C,(NZ_C-ibottom+1)
     1                     ,0.,1e8,1e8,-1,0,-1848,0)

                do i = 1,N_CONTOURS
                    cvalue = factor(i)
                    if(cvalue .ge. abs(cint))then
                        call conrec(field_vert(1,ibottom)
     1                             ,NX_C,NX_C,(NZ_C-ibottom+1)
     1                             ,cvalue,cvalue,1e-6,-1,0,-1848,0)
                        call conrec(field_vert(1,ibottom)
     1                             ,NX_C,NX_C,(NZ_C-ibottom+1)
     1                             ,-cvalue,-cvalue,1e-6,-1,0,-1848,0)
                    endif
                enddo ! i

              else
                call remap_field_2d(
     1                            NX_C,1,NX_C
     1                           ,NZ_C,ibottom,NZ_C
     1                           ,NX_P, ixl_remap, NX_P-ixh_remap
     1                           ,NX_P, iyl_remap, NX_P-iyh_remap
     1                           ,field_vert,field_vert3,r_missing_data)


                call array_range(field_vert3,NX_P,NX_P,rmin,rmax
     1                          ,r_missing_data)

                cmax = max(abs(rmin),abs(rmax))

                call conrec_line
     1              (field_vert3(1,ibottom),NX_P,NX_P,NX_P
     1                             ,0.,0.,cint,-1,0,-1848,0)

                do i = 1,N_CONTOURS
                    cvalue = factor(i)
                    if(cvalue .ge. abs(cint) .and. cvalue .le. cmax)then       
                        cint_in = 2 * cvalue
                        call conrec_line
     1                      (field_vert3(1,ibottom),NX_P,NX_P,NX_P
     1                             ,-cvalue,cvalue,cint_in,-1,0,-1848,0)       
!                       call conrec_line
!    1                      (field_vert3(1,ibottom),NX_P,NX_P,NX_P
!    1                             ,-cvalue,-cvalue,cint,-1,0,-1848,0)
                    endif
                enddo ! i
 
             endif

           endif ! cint > 0
        endif ! i_contour = 1

        if(i_contour .eq. 2)then ! Plot Wind Barbs
            call setusv_dum(2hIN,icolors(i_graphics_overlay))

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, rleft, right,1)
            du=0.4
            rot = 0.
            do i = NX_C,1,-3
                rk_terrain = 
     1          max(height_to_zcoord(terrain_vert1d(i),istatus),1.0)
                do k = ibottom,NZ_C
                    if(u_vert(i,k) .ne. r_missing_data .and.
     1                 v_vert(i,k) .ne. r_missing_data .and.
     1                abs(u_vert(i,k)) .lt. 1e6      )then
                        x1 = i
                        y1 = (k-ibottom) * float(NX_C-1)
     1                                   / float(NZ_C-ibottom) + 1.
                        call uv_to_disp(u_vert(i,k),
     1                                  v_vert(i,k),
     1                                  dir,
     1                                  spd_ms)
                        if(dir .gt. -400. .and. 
     1                       k .ge. int(rk_terrain))then
                            call barbs(spd_ms/mspkt,dir,x1,y1,du,rot
     1                                  ,-1e10,+1e10,-1e10,+1e10)
                        endif
                    endif
                enddo ! k
            enddo ! i
        endif ! i_contour = 2

        if(i_contour .eq. 3)then ! Plot Cloud Type
            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, rleft, right,1)
            du=0.4
            rot = 0.
            do i = NX_C,1,-2
                do k = ibottom,NZ_C
                    i_cloud_type = field_2d(i,k)

                    if(i_cloud_type .ne. 0)then
                        x = i
                        y = (k-ibottom) * float(NX_C-1)
     1                                  / float(NZ_C-ibottom) + 1.

                        c2_cloud_type = c2_cloud_types(i_cloud_type)
!                       call upcase(c2_cloud_type,c2_cloud_type)
                        write(6,*)i_cloud_type,c2_cloud_type,x,y
!                       call pwrity (x, y, c2_cloud_type, 2, 0, 0, 0)
                        CALL PCMEQU(x, y, c2_cloud_type,.0063,0,0)       
                    endif

                enddo ! k
            enddo ! i
        endif ! i_contour = 3

        if(i_contour .eq. 4)then ! Plot Icing Index
            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, rleft, right,1)
            du=0.4
            rot = 0.
            do i = NX_C,1,-1
                do k = ibottom,NZ_C
                    if(icing_index_2d(i,k) .ne. 0)then
                        x = i
                        y = (k-ibottom) * float(NX_C-1)
     1                                  / float(NZ_C-ibottom) + 1.
!                       write(c1_string,2011)icing_index_2d(i,k)
                        c1_string = '*'
2011                    format(i1)

!                       Normal Overlay Colors
                        call setusv_dum(2hIN,icolors(i_graphics_over
     1lay))

!                       Multicolored ISI display
                        if(.false.)then
                        if(icing_index_2d(i,k) .gt. 3)then
                            i_color_value = icing_index_2d(i,k) - 3
                        else
                            i_color_value = icing_index_2d(i,k)
                        endif

                        if(i_color_value .eq. 1)call setusv_dum(2hIN
     1,7) ! Y
                        if(i_color_value .eq. 2)call setusv_dum(2hIN
     1,245) ! O
                        if(i_color_value .eq. 3)call setusv_dum(2hIN
     1,111) ! R
                        endif

                        call pwrity (x, y, c1_string, 1, 0, 0, 0)
                    endif
                enddo ! k
            enddo ! i
        endif ! i_contour = 4


        if(i_contour .eq. 5)then ! Plot Precip Type
            call setusv_dum(2hIN,icolors(i_graphics_overlay))

            call set(vxmin, vxmax, vymin, vymax
     1             , rleft, right, rleft, right,1)
            do i = NX_C,1,-1
                do k = ibottom+1,NZ_C
                    i_precip_type = field_2d(i,k)
                    if(i_precip_type .ne. 0)then
                        x = i
                        y = (k-ibottom) * float(NX_C-1)
     1                                  / float(NZ_C-ibottom) + 1.

                        c2_precip_type = c1_precip_types(i_precip_type)
                        call pwrity (x, y, c2_precip_type, 2, 0, 0, 0)
                    endif
                enddo ! k
            enddo ! i
        endif ! i_contour = 5


        goto100

9999    continue

!       Contour in the terrain surface
        call set(vxmin, vxmax, vymin, vymax
     1         , rleft, right, bottom, top,1)
        n_div = 20
        call setusv_dum(2hIN,3)

!       Read in sfc pressure
        i4time_tol = 43200
        var_2d = 'PS'
        ext = 'lsx'
        call get_laps_2dgrid(i4time_ref,i4time_tol,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                  ,pres_2d,0,istatus)
        IF(istatus .ne. 1)THEN
            write(6,*)' Error Reading Surface Pressure Analysis'
            write(6,*)
     1        ' Converting Terrain to Sfc Pressure with Std Atmosphere'       
            istat_sfc_pres = 0
        else
            call interp_2d(pres_2d,pres_1d,xlow,xhigh,ylow,yhigh,
     1                     NX_L,NY_L,NX_C,r_missing_data)
            istat_sfc_pres = 1
        endif

        do i = 1,NX_C
            xcoord(i) = i
            if(istat_sfc_pres .eq. 1)then
                ycoord(i) = max(zcoord_of_pressure(pres_1d(i)),1.0)
            else
                ycoord(i) = 
     1            max(height_to_zcoord(terrain_vert1d(i),istatus),1.0)
            endif

            if(i .gt. 1)then
                do j = 0,n_div
                    frac = float(j) / float(n_div)
                    ybottom = bottom
                    ytop = ycoord(i-1) * (1.-frac) + ycoord(i) * frac
                    xval = float(i-1) + frac
                    call line(xval,ybottom,xval,ytop)
                enddo ! j
            endif
        enddo ! i

        call setusv_dum(2hIN,7)

        if(l_sta)then ! Label location of station

!          This lets up plot outside the main box
           call set(.00, 1.0, vymin2 , vymax2
     1             , rleft  - width/8.   , right + width/8.
     1             , bottom - r_height/8., top   + r_height/8. ,1)

           x = pos_sta
           write(6,*)'     Labeling ',c3_sta,pos_sta

           call line(x,bottom,x,bottom - .004 * r_height)

           y = bottom - .010 * r_height
           call upcase(c3_type,c3_type)

           if(l_arrival_gate)then
               write(c9_string,2039)c3_sta,c5_arrival_gate(1:3)
           else
               write(c9_string,2039)c3_sta,c3_type
           endif

2039       format(1x,a3,'-',a3)

!          call pwrity ((x+0.15), y, c9_string, 9, 0, 0, 0)
           call pcmequ ((x+0.15), y, c9_string, .0050, 0, 0)

           if(.not. l_atms)then
               i4time_label = i4time_ref/laps_cycle_time*laps_cycle_time
     1                                          -laps_cycle_time

               call label_other_stations(i4time_label,standard_longitude       
     1                                  ,y,xsta,lat,lon,NX_L,NY_L
     1                                  ,xlow,xhigh,ylow,yhigh
     1                                  ,NX_C,bottom,r_height,maxstns)

           endif

        else ! l_sta = .false.

           if(.not. l_atms)then

!              This lets up plot outside the main box
               call set(.00, 1.0, vymin2 , vymax2
     1                , rleft - width/8., right + width/8.
     1                , bottom - r_height/8., top + r_height/8.,1)

               xsta = -10000.

               i4time_label = i4time_ref/laps_cycle_time*laps_cycle_time
     1                                          -laps_cycle_time

               y = bottom - .015 * r_height

               call label_other_stations(i4time_label,standard_longitude       
     1                                  ,y,xsta,lat,lon,NX_L,NY_L
     1                                  ,xlow,xhigh,ylow,yhigh,NX_C
     1                                  ,bottom,r_height,maxstns)

           endif

        endif ! l_sta = .true.

        call frame

        return
        end


        subroutine interp_3d(array_in,array_out,xlow,xhigh,ylow,yhigh,
     1                       NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, NZ_L, NX_C, NZ_C as
!                                   dummy arguments
!       97-Aug-14     Ken Dritz     Added r_missing_data as dummy argument
!       97-Aug-14     Ken Dritz     Removed (commented out) parameter
!                                   declarations for NX_C, NZ_C
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, NX_C, r_missing_data
!                                   to interp_2d

        integer*4 NX_C,NZ_C
!       parameter (NX_C = 61) ! NX_L
!       parameter (NZ_C = NZ_L)

        real*4 array_in(NX_L,NY_L,NZ_L)
        real*4 array_out(NX_C,NZ_C)

!       array_out = r_missing_data

        do k = 1,NZ_L
             call interp_2d(array_in(1,1,k),array_out(1,k)
     1                     ,xlow,xhigh,ylow,yhigh
     1                     ,NX_L,NY_L,NX_C,r_missing_data)
        enddo ! i

        return
        end

        subroutine interp_3dn(array_in,array_out,xlow,xhigh,ylow,yhigh,
     1                        NX_L,NY_L,NZ_L,NX_C,NZ_C)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, NZ_L, NX_C, NZ_C as
!                                   dummy arguments
!       97-Aug-14     Ken Dritz     Removed (commented out) parameter
!                                   declarations for NX_C, NZ_C
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, NX_C to interp_2dn

!       Nearest neighbor interpolation of a vertical X-sect

        integer*4 NX_C,NZ_C
!       parameter (NX_C = 61) ! NX_L
!       parameter (NZ_C = NZ_L)

        real*4 array_in(NX_L,NY_L,NZ_L)
        real*4 array_out(NX_C,NZ_C)

        do k = 1,NZ_L
             call interp_2dn(array_in(1,1,k),array_out(1,k)
     1                                 ,xlow,xhigh,ylow,yhigh,
     1                       NX_L,NY_L,NX_C)
        enddo ! i

        return
        end

        subroutine interp_3d_spread(array_in,array_out,xlow,xhigh,ylow,y
     1high,NX_L,NY_L,NZ_L,NX_C,NZ_C,r_missing_data)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, NZ_L, NX_C, NZ_C as
!                                   dummy arguments
!       97-Aug-14     Ken Dritz     Added r_missing_data as dummy argument
!       97-Aug-14     Ken Dritz     Removed (commented out) parameter
!                                   declarations for NX_C, NZ_C
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, NX_C, r_missing_data
!                                   to interp_2d_spread

        integer*4 NX_C,NZ_C
!       parameter (NX_C = 61) ! NX_L
!       parameter (NZ_C = NZ_L)

        real*4 array_in(NX_L,NY_L,NZ_L)
        real*4 array_out(NX_C,NZ_C)

        do k = 1,NZ_L
             call interp_2d_spread(array_in(1,1,k),array_out(1,k)
     1                                 ,xlow,xhigh,ylow,yhigh,
     1                             NX_L,NY_L,NX_C,r_missing_data)
        enddo ! i

        return
        end


        subroutine interp_2d(array_in,array_out,xlow,xhigh,ylow,yhigh,
     1                       NX_L,NY_L,NX_C,r_missing_data)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, NX_C, and r_missing_data
!                                   as dummy arguments
!       97-Aug-14     Ken Dritz     Removed (commented out) parameter
!                                   declarations for NX_C, NZ_C (the latter
!                                   was unused)
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for

        integer*4 NX_C,NZ_C
!       parameter (NX_C = 61) ! NX_L
!       parameter (NZ_C = NZ_L)

        real*4 array_in(NX_L,NY_L)
        real*4 array_out(NX_C)

        if(NX_C .gt. 1)then
            deltax = (xhigh - xlow) / (float(NX_C) - 1.)
            deltay = (yhigh - ylow) / (float(NX_C) - 1.)
        else
            deltax = 0.
            deltay = 0.
        endif

!       Bilinearly interpolate from 2d grid to 1d array
        do ii = 1,NX_C
            ri = xlow + (float(ii-1) * deltax)
            rj = ylow + (float(ii-1) * deltay)

            i = int(ri)
            if(i .eq. NX_L)i=i-1

            j = int(rj)
            if(j .eq. NY_L)j=j-1

            fraci = ri - i
            fracj = rj - j

            Z1=array_in(i  , j  )
            Z2=array_in(i+1, j  )
            Z3=array_in(i+1, j+1)
            Z4=array_in(i  , j+1)

            if(    z1 .ne. r_missing_data
     1     .and. z2 .ne. r_missing_data
     1     .and. z3 .ne. r_missing_data
     1     .and. z4 .ne. r_missing_data)then

                array_out(ii) =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                - (Z2+Z4-Z3-Z1)*fraci*fracj

            else
                array_out(ii) = r_missing_data

            endif

!            if(array_out(ii) .ne. r_missing_data
!       1       .and. abs(array_out(ii) .lt. 1e-20)then
!                array_out(ii) = r_missing_data
!            endif

        enddo ! ii

        return
        end


        subroutine interp_2dn(array_in,array_out,xlow,xhigh,ylow,yhigh,
     1                        NX_L,NY_L,NX_C)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, NX_C as dummy arguments
!       97-Aug-14     Ken Dritz     Removed (commented out) parameter
!                                   declarations for NX_C, NZ_C (the latter
!                                   was unused)
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for

!       Nearest neighbor interpolation of a 2d array to a line.

        integer*4 NX_C,NZ_C
!       parameter (NX_C = 61) ! NX_L
!       parameter (NZ_C = NZ_L)

        real*4 array_in(NX_L,NY_L)
        real*4 array_out(NX_C)

        deltax = (xhigh - xlow) / (float(NX_C) - 1.)
        deltay = (yhigh - ylow) / (float(NX_C) - 1.)

!       Interpolate from 2d grid to 1d array using the nearest neighbor
        do ii = 1,NX_C
            ri = xlow + (float(ii-1) * deltax)
            rj = ylow + (float(ii-1) * deltay)

            i = nint(ri)
            j = nint(rj)

            array_out(ii) = array_in(i,j)

        enddo ! ii

        return
        end


        subroutine interp_2d_spread(array_in,array_out,xlow,xhigh,ylow,y
     1high,NX_L,NY_L,NX_C,r_missing_data)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, NX_C, r_missing_data
!                                   as dummy arguments
!       97-Aug-14     Ken Dritz     Removed (commented out) parameter
!                                   declarations for NX_C and NZ_C (the
!                                   latter was unused)
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for

        integer*4 NX_C,NZ_C
!       parameter (NX_C = 61) ! NX_L
!       parameter (NZ_C = NZ_L)

        real*4 array_in(NX_L,NY_L)
        real*4 array_out(NX_C)

        deltax = (xhigh - xlow) / (float(NX_C) - 1.)
        deltay = (yhigh - ylow) / (float(NX_C) - 1.)

!       Bilinearly interpolate from 2d grid to 1d array
        do ii = 1,NX_C
            ri = xlow + (float(ii-1) * deltax)
            rj = ylow + (float(ii-1) * deltay)

            i = int(ri)
            if(i .eq. NX_L)i=i-1

            j = int(rj)
            if(j .eq. NY_L)j=j-1

            fraci = ri - i
            fracj = rj - j

            Z1=array_in(i  , j  )
            Z2=array_in(i+1, j  )
            Z3=array_in(i+1, j+1)
            Z4=array_in(i  , j+1)

            if(    z1 .ne. r_missing_data
     1     .and. z2 .ne. r_missing_data
     1     .and. z3 .ne. r_missing_data
     1     .and. z4 .ne. r_missing_data)then

                array_out(ii) =  Z1+(Z2-Z1)*fraci+(Z4-Z1)*fracj
     1                - (Z2+Z4-Z3-Z1)*fraci*fracj

            else
                array_out(ii) = r_missing_data

            endif

!            if(array_out(ii) .ne. r_missing_data
!       1       .and. abs(array_out(ii) .lt. 1e-20)then
!                array_out(ii) = r_missing_data
!            endif

        enddo ! ii

        return
        end

        subroutine interp_3dc(array_in,array_out,xlow,xhigh,ylow,yhigh,
     1                        NX_L,NY_L,NX_C,r_missing_data)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, NX_C, r_missing_data
!                                   as dummy arguments
!       97-Aug-14     Ken Dritz     Removed (commented out) parameter
!                                   declarations for NX_C, NZ_C (the
!                                   latter was unused)
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, NX_C, r_missing_data
!                                   to interp_2d

        include 'laps_cloud.inc'

        integer*4 NX_C,NZ_C
!       parameter (NX_C = 61) ! NX_L
!       parameter (NZ_C = kcloud)

        real*4 array_in(NX_L,NY_L,kcloud)
        real*4 array_out(NX_C,kcloud)

        do k = 1,kcloud
             call interp_2d(array_in(1,1,k),array_out(1,k)
     1                                 ,xlow,xhigh,ylow,yhigh,
     1                      NX_L,NY_L,NX_C,r_missing_data)
        enddo ! i

        return
        end


        subroutine xsect_endpoints(xsta,ysta,azi_xsect,
     1                  xlow,ylow,xhigh,yhigh,pos_sta,istatus,
     1                  NX_L,NY_L,NX_C)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, NX_C as dummy arguments
!       97-Aug-14     Ken Dritz     Removed (commented out) parameter
!                                   declarations for NX_C, NZ_C (the latter
!                                   was unused)
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for

        include 'trigd.inc'

        logical l_left, l_right, l_top, l_bottom

        integer*4 NX_C,NZ_C
!       parameter (NX_C = 61) ! NX_L
!       parameter (NZ_C = NZ_L)

        ANGDIF(X,Y)=MOD(X-Y+540.,360.)-180.
        COTAND(X) = TAND(90.-X)

!           Calculate endpoints of X-Sect from Waypoint and Azimuth
            azi_met = 90. - azi_xsect

!           Intersection with right edge
            l_right = .false.
            if(cotand(azi_met) .ne. 0.)then
                yright = ysta + (NX_L - xsta) * tand(azi_met)
                if(yright .le. float(NY_L) .and. yright .ge. 1.)then
                    l_right = .true.
                    write(6,311)float(NX_L),yright
311                 format('  Right Edge  ',2f7.2)
                endif
            endif

!           Intersection with left edge
            l_left = .false.
            if(cotand(azi_met) .ne. 0.)then
                yleft = ysta - (xsta - 1.) * tand(azi_met)
                if(yleft .le. float(NY_L) .and. yleft .ge. 1.)then
                    l_left = .true.
                    write(6,312)1.,yleft
312                 format('  Left Edge   ',2f7.2)
                endif
            endif

!           Intersection with top edge
            l_top = .false.
            if(tand(azi_met) .ne. 0.)then
                xtop = xsta + (NY_L - ysta) * cotand(azi_met)
                if(xtop .le. float(NX_L) .and. xtop .ge. 1.)then
                    l_top = .true.
                    write(6,313)xtop,float(NY_L)
313                 format('  Top Edge    ',2f7.2)
                endif
            endif

!           Intersection with bottom edge
            l_bottom = .false.
            if(tand(azi_met) .ne. 0.)then
                xbottom = xsta - (ysta - 1.) * cotand(azi_met)
                if(xbottom .le. float(NX_L) .and. xbottom .ge. 1.)then
                    l_bottom = .true.
                    write(6,314)xbottom,1.
314                 format('  Bottom Edge ',2f7.2)
                endif
            endif

!           Now we can get the endpoints
            if(abs(angdif(azi_xsect,45.)) .le. 45.)then
                if(l_left)then
                    xlow = 1.
                    ylow = yleft
                elseif(l_bottom)then
                    xlow = xbottom
                    ylow = 1.
                endif

                if(l_right)then
                    xhigh = NX_L
                    yhigh = yright
                elseif(l_top)then
                    xhigh = xtop
                    yhigh = NY_L
                endif

            elseif(abs(angdif(azi_xsect,135.)) .le. 45.)then
                if(l_left)then
                    xlow = 1.
                    ylow = yleft
                elseif(l_top)then
                    xlow = xtop
                    ylow = NY_L
                endif

                if(l_right)then
                    xhigh = NX_L
                    yhigh = yright
                elseif(l_bottom)then
                    xhigh = xbottom
                    yhigh = 1
                endif

            elseif(abs(angdif(azi_xsect,225.)) .le. 45.)then
                if(l_left)then
                    xhigh = 1.
                    yhigh = yleft
                elseif(l_bottom)then
                    xhigh = xbottom
                    yhigh = 1.
                endif

                if(l_right)then
                    xlow = NX_L
                    ylow = yright
                elseif(l_top)then
                    xlow = xtop
                    ylow = NY_L
                endif

            elseif(abs(angdif(azi_xsect,315.)) .le. 45.)then
                if(l_left)then
                    xhigh = 1.
                    yhigh = yleft
                elseif(l_top)then
                    xhigh = xtop
                    yhigh = NY_L
                endif

                if(l_right)then
                    xlow = NX_L
                    ylow = yright
                elseif(l_bottom)then
                    xlow = xbottom
                    ylow = 1
                endif

            endif

            if(xlow .ne. xhigh)then
                pos_sta = 1. + (NX_C-1.) * (xsta-xlow)/(xhigh-xlow)
            else
                pos_sta = 1. + (NX_C-1.) * (ysta-ylow)/(yhigh-ylow)
            endif

            write(6,91)xlow,ylow,xhigh,yhigh,pos_sta
91          format(1x,' End Points  ',2f7.2,2x,2f7.2,'  pos_sta',f7.2)


        return
        end

        subroutine label_other_stations(i4time,standard_longitude,y,xsta
     1          ,lat,lon,ni,nj
     1          ,xlow,xhigh,ylow,yhigh,nx_c,bottom,r_height,maxstns)

!       97-Aug-14     Ken Dritz     Added maxstns as dummy argument
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for
!       97-Aug-25     Steve Albers  Removed /read_sfc_cmn/.

!       This routine labels stations on the X-sect in a logical manner

        real*4 stapos_a(maxstns+1)

        real*4 lat(ni,nj),lon(ni,nj)

        real*4 lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real*4 cover_s(maxstns), hgt_ceil(maxstns), hgt_low(maxstns)
        real*4 t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real*4 dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns)
     1       , ffg_s(maxstns)
        real*4 vis_s(maxstns)
        character stations(maxstns)*3, wx_s(maxstns)*8    ! c5_stamus

!       Declarations for new read_surface routine
!       New arrays for reading in the SAO data from the LSO files
        real*4   pstn(maxstns),pmsl(maxstns),alt(maxstns)
     1          ,store_hgt(maxstns,5)
        real*4   ceil(maxstns),lowcld(maxstns),cover_a(maxstns)
     1          ,vis(maxstns)
     1                                          ,rad(maxstns)

        Integer*4   obstime(maxstns),kloud(maxstns),idp3(maxstns)

        Character   obstype(maxstns)*8
     1             ,store_emv(maxstns,5)*1,store_amt(maxstns,5)*4

!       common /read_sfc_cmn/ lat_s,lon_s,elev_s,cover_s,hgt_ceil,hgt_lo
!    1w
!    1                ,t_s,td_s,pr_s,sr_s,dd_s,ff_s,ddg_s,ffg_s,vis_s
c
        character atime*24, infile*255
        character directory*150,ext*31

        character*9 c9_string
        character*13 filename13

        character*2 icompass(8)
        data icompass/'N ','NE','E ','SE','S ','SW','W ','NW'/

        write(6,*)
     1  ' Reading Station locations from read_sfc for labelling '
        ext = 'lso'
        call get_directory(ext,directory,len_dir) ! Returns top level directory
        infile = directory(1:len_dir)//filename13(i4time,ext(1:3))

        call read_surface_old(infile,maxstns,atime,n_meso_g,n_meso_pos,
     &      n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,n_obs_pos_g,
     &      n_obs_b,n_obs_pos_b,stations,obstype,lat_s,lon_s,
     &      elev_s,wx_s,t_s,td_s,dd_s,ff_s,ddg_s,
     &      ffg_s,pstn,pmsl,alt,kloud,ceil,lowcld,cover_a,rad,idp3,
     &      store_emv,store_amt,store_hgt,vis,obstime,istatus)

100     write(6,*)'     n_obs_b',n_obs_b

        if(n_obs_b .gt. maxstns .or. istatus .ne. 1)then
            write(6,*)' Too many stations, or no file present'
            istatus = 0
            return
        endif

        stapos_a(1) = xsta

        sect_length = sqrt((xhigh-xlow)**2 + (yhigh-ylow)**2)

        write(6,*)'i,xsta,ysta,frac,rmin,stations(i),'
     1                            ,'stapos,dist_min,ran,azi'

        iplot = 1

        do isweep = 0,7
          do i = 1,n_obs_b
            call latlon_to_rlapsgrid(lat_s(i),lon_s(i),lat,lon
     1                          ,ni,nj,xsta,ysta,istatus)

!           Calculate distance from station to X-sect
!           Rmin is in terms of LAPS Grid Points
            call closest(xlow-xsta,ylow-ysta,xhigh-xlow,yhigh-ylow,frac,
     1rmin)

            stapos = 1. + frac * float(nx_c-1)

!           rmin = rmin * float(nx_c-1) / sect_length

            if(frac .ge. 0. .and. frac .le. 1.)then
              if(int(abs(rmin*2.)) .eq. isweep)then

!               Find distance from previously written out locations
                dist_min = min(abs(nx_c-stapos),abs(stapos-1.))
                do ii = 1,iplot
                    dist = abs(stapos-stapos_a(ii))
                    if(dist .lt. dist_min)then
                        dist_min = dist
                    endif
                enddo ! ii


                if(dist_min .gt. 5.0)then ! No other stations in the way
                    xclo = xlow + frac * (xhigh - xlow)
                    yclo = ylow + frac * (yhigh - ylow)

                    iplot = iplot + 1
                    stapos_a(iplot) = stapos

                    xdelt = xclo - xsta
                    ydelt = yclo - ysta

                    call xy_to_met_xm(xdelt,ydelt,ran,azi,istatus)

                    azi = azi + (lon_s(i) - standard_longitude)

                    i_cardinal_pt = mod(int((azi+22.5)/45.),8) + 1

                    if(icompass(i_cardinal_pt)(2:2) .eq. ' ')then
                        write(c9_string,2038)nint(ran*10.)
     1                    ,icompass(i_cardinal_pt),stations(i)(1:3)
2038                    format(i2,a2,a3)
!                       call pwrity ((stapos+0.6), y, c9_string
!    1                                                  , 9, 0, 0, 0)
                        call pcmequ ((stapos+0.8), y, c9_string
     1                               , .0060, 0, 0)
                    else
                        write(c9_string,2039)nint(ran*10.)
     1                    ,icompass(i_cardinal_pt),stations(i)(1:3)
2039                    format(i2,a2,' ',a3)

!                       call pwrity ((stapos+.15), y, c9_string
!    1                                                  , 9, 0, 0, 0)
                        call pcmequ ((stapos+.35), y, c9_string
     1                               , .0060, 0, 0)
                    endif


                    write(6,11)i,xsta,ysta,frac,rmin,stations(i)(1:3)
     1                  ,stapos,dist_min
     1                  ,ran,azi,c9_string,isweep
11                  format(i4,2f6.1,f6.3,f6.1,1x,a3,4f6.1,1x,a9,i2)

                    call line(stapos,bottom,stapos
     1                       ,bottom - .012 * r_height)
!                   write(6,*) ' Call line',stapos,bottom,r_height

                endif

              endif
            endif
          enddo ! stations
        enddo ! isweep

        return
        end

        subroutine closest(r1,r2,rd1,rd2,tclo,rmin)
        implicit real*4 (a-z)

c       write(6,*)

        rnum =  -(      R1 * RD1
     1  +       R2 * RD2)
        denom =  (      RD1 * RD1
     1  +       RD2 * RD2)

        tclo = rnum/denom ! + t0

        rdot = sqrt(RD1**2 + RD2**2)
        rnum = (R1*RD2) - (R2*RD1)
        rmin = rnum/rdot
c       write(6,*)' rnum,rdot = ',rnum,rdot

c       write(6,1)rmin,R1,R2,RD1,RD2,tclo
1       format(f12.10,5f12.8)
        return
        end

c

        subroutine remap_field_2d(nx_in,ixlow_in,ixhigh_in
     1                           ,ny_in,iylow_in,iyhigh_in
     1                           ,nx_out,ixlow_out,ixhigh_out
     1                           ,ny_out,iylow_out,iyhigh_out
     1                           ,field_in,field_out,r_missing_data)

        real*4 field_in(nx_in,ny_in)
        real*4 field_out(nx_out,ny_out)

        call constant(field_out,r_missing_data,nx_out,ny_out)

        rxlow_out  = ixlow_out
        rxhigh_out = ixhigh_out
        rylow_out  = iylow_out
        ryhigh_out = iyhigh_out

        rxlow_in   = ixlow_in
        rxhigh_in  = ixhigh_in 
        rylow_in   = iylow_in 
        ryhigh_in  = iyhigh_in

        do ixout = ixlow_out, ixhigh_out
        do iyout = iylow_out, iyhigh_out
            rxout = ixout
            ryout = iyout

            arg = rxout
            call stretch(rxlow_out,rxhigh_out,rxlow_in,rxhigh_in,arg)
            rxin = arg

            arg = ryout
            call stretch(rylow_out,ryhigh_out,rylow_in,ryhigh_in,arg)
            ryin = arg

            call bilinear_laps(rxin,ryin,nx_in,ny_in,field_in,result)       
            field_out(ixout,iyout) = result

        enddo
        enddo

        return
        end

        subroutine input_product_info(
     1                              i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,ndim                    ! I
     1                             ,c_prodtype              ! O
     1                             ,ext                     ! O
     1                             ,directory               ! O
     1                             ,a9time                  ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                             ,istatus)                ! O

        integer       maxbgmodels
        parameter     (maxbgmodels=10)

        integer       n_fdda_models
        integer       l,len_dir,lfdda
        integer       istatus
        character*9   c_fdda_mdl_src(maxbgmodels)
        character*(*) directory
        character*(*) ext
        character*20  c_model
        character*10  cmds
        character*1   cansw
        character*150 c_filenames(1000)

        character*1 c_prodtype
        character*4 fcst_hhmm
        character*9 a9time
        character*13 a13_time

        logical l_parse

        write(6,*)' Subroutine input_product_info...'

        write(6,1)
 1      format('  Product type: analysis [a], background [b]'
     1        ,', balance [n], forecast [f] ? ',$)        

        read(5,2)c_prodtype
 2      format(a)

        call upcase(c_prodtype,c_prodtype)

        if(c_prodtype .eq. 'A' .or. c_prodtype .eq. 'N')then
            return

        else
            if(c_prodtype .eq. 'B')then
                if(ndim .eq. 2)then
                    ext = 'lgb'
                else
                    ext = 'lga'
                endif

            elseif(c_prodtype .eq. 'F')then
                if(ndim .eq. 2)then
                    ext = 'fsf'
                else
                    ext = 'fua'
                endif

            endif

            call input_background_info(
     1                              ext                     ! I
     1                             ,directory               ! O
     1                             ,i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,a9time                  ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                             ,istatus)                ! O

        endif

        return
        end

