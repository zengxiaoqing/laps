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



        subroutine plot_station_locations(i4time,lat,lon,topo
     1                                   ,ni,nj,nk,iflag   
     1                                   ,maxsta,c_field,zoom
     1                                   ,namelist_parms,plot_parms
     1                                   ,MAX_SND_GRID,MAX_SND_LEVELS      ! I
     1                                   ,atime_s
     1                                   ,c_label,i_overlay)

        include 'lapsplot.inc'

!       Declarations for 'read_surface_snd', etc.
        include 'read_sfc.inc'

!       97-Aug-14     Ken Dritz     Added maxsta as dummy argument
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for
!       97-Aug-25     Steve Albers  Removed /read_sfc_cmn/.

!       This routine labels station locations on the H-sect

        real lat(ni,nj),lon(ni,nj),topo(ni,nj)

        character c_label*(*)
        character directory*150,ext*31,ext_lso*6
        character*255 c_filespec
        character*9 c9_string, asc_tim_9
        character*13 filename13
        character*4 c_field
        character*3 c3_presob
        character*5 c_staname

!       Declarations for 'read_surface_sa' call
!       New arrays for reading in the SAO data from the LSO files
        real   ceil(maxsta),lowcld(maxsta),cover_a(maxsta)
     1          ,vis(maxsta)

!       Integer   kloud(maxsta)

!       character atype(maxsta)*6

!       Declarations for 'read_sfc_precip' call
	character filetime*9, infile*256, btime*24
        character c20_stations*20
	character dum*132

        logical l_parse

        common /plot_field_cmn/ i_plotted_field
  
        i_plotted_field = 1

        ANGD = 0.
        CNTR = 0.

        call get_sfc_badflag(badflag,istatus)

        lun = 42

        if(c_field(1:2) .eq. 'mw')then ! Read Mesowx
            ext = 'mesowx'
            call get_directory(ext,directory,len_dir) ! Mesowx directory
            infile = 
     1      directory(1:len_dir)//'lfmpost_points.txt'

            open(lun,file=infile,status='unknown')

            i = 1
 4          read(lun,5,end=15)lat_s(i),lon_s(i)
 5          format(2x,1x,10x,1x,f8.0,f10.0)

            i = i+1

            go to 4            

 15         n_obs_b = i-1

        elseif(c_field(1:2) .eq. 'ov' .OR.
     1         c_field(1:2) .eq. 'or')then ! Read Cloud Obs
            write(6,*)' Calling read_cloud_obs...'
            i4time_lso = i4time
            call read_cloud_obs(i4time,maxsta,                           ! I
     &            n_obs_b,stations,reptype,                              ! O
     &            autostntype,                                           ! O
     &            lat_s,lon_s,elev_s,wx_s,t_s,td_s,vis_s,solar_s,        ! O
     &            kloud_s,store_amt,store_hgt,obstime,istatus)           ! O

            write(6,*)'      n_obs_b:',n_obs_b       
            call make_fnam_lp(i4time_lso,asc_tim_9,istatus)
            n_obs_g = n_obs_b

        else ! Read LSO or LSO_QC plus SND 

            call get_filespec('lso',2,c_filespec,istatus)
            call get_file_time(c_filespec,i4time,i4time_lso)

            if(i4time_lso .eq. 0)then
                write(6,*)
     1          ' No LSO files available for station plotting'       
                return
            endif

            call make_fnam_lp(i4time_lso,asc_tim_9,istatus)

            ext = 'lso'
            call get_directory(ext,directory,len_dir) ! Returns top level directory
            if(c_field(1:1) .eq. 'q')then ! LSO_QC file
                infile = 
     1          directory(1:len_dir)//filename13(i4time_lso,ext(1:3))
     1                              //'_qc'    

                ext_lso = 'lso_qc'

            else ! Regular LSO file
                infile = 
     1          directory(1:len_dir)//filename13(i4time_lso,ext(1:3))  

                ext_lso = 'lso'

            endif

            if(ext_lso .eq. 'lso')then ! phase in call to read_surface_data
              write(6,*)' Calling read_surface_data...',ext_lso
     1                                                 ,asc_tim_9      
              call read_surface_data(i4time,atime_s,n_obs_g,n_obs_b, !regular LSO
     &         obstime,wmoid,stations,provider,wx_s,reptype,autostntype,       
     &         lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &         alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &         pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &         td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &         sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,     
     &         maxsta,istatus)
            else
              write(6,*)' Calling read_surface_dataqc...',ext_lso
     1                                                 ,asc_tim_9      
              call read_surface_dataqc(i4time,atime_s,n_obs_g,n_obs_b, !regular LSO
     &         obstime,wmoid,stations,provider,wx_s,reptype,autostntype,       
     &         lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &         alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &         pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &         td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &         sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,
     &         maxsta,istatus)

            endif

            write(6,*)'     n_obs_g:',n_obs_g,'      n_obs_b:',n_obs_b       

            if(n_obs_b .gt. maxsta .or. istatus .ne. 1)then
                write(6,*)' Too many stations, or no file present'
                istatus = 0
                return
            endif

            write(6,*)' Calling read_sfc_snd...'
 	    call read_sfc_snd(i4time,atime_s,n_obs_g,n_obs_b, ! regular SND
     &        obstime,wmoid,stations,provider,wx_s,reptype,autostntype,       
     &        lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &        alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &        pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,       
     &        td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &        sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,
     &        maxsta,lat,lon,ni,nj,nk,                                    ! I
     &        MAX_SND_GRID,MAX_SND_LEVELS,                                ! I
     &        topo,                                                       ! I
     &        istatus)

            write(6,*)'     n_obs_g:',n_obs_g,'      n_obs_b:',n_obs_b       

            if(n_obs_b .gt. maxsta)then
                write(6,*)' Too many stations'
                istatus = 0
                return
            endif

            i_rh_convert = 0
	    do i=1,n_obs_b ! Preprocess the obs
!               Convert RH to dewpoint if dewpoint is missing 
                if(t_s(i) .ne. badflag .and. td_s(i) .eq. badflag 
     1                                 .and. rh_s(i) .ne. badflag )then       
                    t_c = f_to_c(t_s(i))
                    dwpt_c = dwpt(t_c,rh_s(i))
                    td_s(i) = c_to_f(dwpt_c)
                    i_rh_convert = i_rh_convert + 1
                endif
            enddo ! i

            if(i_rh_convert .gt. 0)then
                write(6,*)'# of dewpoints converted from RH = '
     1                   ,i_rh_convert       
            endif

        endif ! Mesowx, Cloud obs, or LSO/SND

        size = 0.5
        call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)
        du = float(nj) / 252.

        obs_size = plot_parms%contour_line_width

!       At zoom=1-zoom_max, make the obs plot larger if there are few obs
!       zoom_eff is inversely related to plot size
        zoom_max = 1.5
        if(zoom .lt. zoom_max)then
            if(n_obs_g .gt. 30)then
                zoom_eff = 1.0                 ! smaller obs 
            else
                zoom_eff = zoom / zoom_max     ! larger obs
            endif
        else
            zoom_eff = zoom / zoom_max         ! larger obs
        endif

        zoom_eff = zoom_eff / obs_size

        du2 = du / zoom_eff

!       call setusv_dum(2HIN,11)

        if(c_field(2:2) .eq. 'v')then ! Ceiling & Visibility
            iflag_cv = 1
        elseif(c_field(2:2) .eq. 'p')then     ! Precip
            write(6,*)' Reading precip obs'
	    call read_sfc_precip(i4time,btime,n_obs_g,n_obs_b,
     &        stations,provider,lat_s,lon_s,elev_s,
     &        pcp1,pcp3,pcp6,pcp24,snow,       
     &        maxsta,jstatus)
            iflag_cv = 2
        elseif(c_field(2:2) .eq. 'g')then ! Soil/Water Temp & Solar Radiation
            iflag_cv = 3
        elseif(c_field(2:2) .eq. 'r')then ! Solar Radiation
            iflag_cv = 4
        else
            iflag_cv = 0
        endif

        write(6,*)' plot_station_locations... ',iflag

        c3_presob = '   '

        if(iflag .ge. 1)then
            if(iflag_cv .eq. 0)then
                write(6,13)
13              format(' Select type of pressure ob [msl,alt,stn]'
     1                ,4x,'default=none      ? ',$)
                read(5,14)c3_presob
14              format(a)
            endif

            if(c_field(1:1) .eq. 'q')then ! LSO_QC file
                c_label = 'Sfc QC Obs   ('//c3_presob//' pres)'
            else
                c_label = 'Sfc Obs      ('//c3_presob//' pres)'
            endif

            if(iflag_cv .eq. 1)then
                c_label(14:50) = ' Cloud Base (hm) & Visibility (miles)'
            elseif(iflag_cv .eq. 2)then
                if(c_field(3:4) .eq. '24')then
                  if(namelist_parms%c_units_type .eq. 'english')then
                    c_label(13:34) = '24Hr Precip (in)'
                  else
                    c_label(13:39) = '24Hr Precip (in)'
                  endif
                elseif(c_field(3:3) .eq. '6')then
                  if(namelist_parms%c_units_type .eq. 'english')then
                    c_label(13:33) = '6Hr Precip (in)'
                  else
                    c_label(13:38) = '6Hr Precip (in)'
                  endif
                elseif(c_field(3:3) .eq. '3')then
                  if(namelist_parms%c_units_type .eq. 'english')then
                    c_label(13:33) = '3Hr Precip (in)'
                  else
                    c_label(13:38) = '3Hr Precip (in)'
                  endif
                elseif(c_field(3:3) .eq. '1')then
                  if(namelist_parms%c_units_type .eq. 'english')then
                    c_label(13:33) = '1Hr Precip (in)'
                  else
                    c_label(13:38) = '1Hr Precip (in)'
                  endif
                else
                  if(namelist_parms%c_units_type .eq. 'english')then
                    c_label(13:33) = '1Hr Pcp/Snw Dpth (in)'
                  else
                    c_label(13:38) = '1Hr Pcp (mm)/Snw Dpth (cm)'
                  endif
                endif
            elseif(iflag_cv .eq. 3)then
                c_label(14:51) =  
     1                     '   Sfc Temp, Soil Moisture & Solar Rad'    
            elseif(iflag_cv .eq. 4)then
                c_label(14:31) =  
     1                     '   Solar Radiation'    
            endif

            if(c_field(1:2) .eq. 'ow')then
                c_label = 'Sfc Obs      (wind)'
            endif

            call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
            call write_label_lplot(ni,nj,c_label,asc_tim_9
     1                            ,plot_parms,namelist_parms
     1                            ,i_overlay,'hsect')       

        endif

        call get_border(ni,nj,x_1,x_2,y_1,y_2)
        call set(x_1,x_2,y_1,y_2,1.,float(ni),1.,float(nj),1)

        call get_r_missing_data(r_missing_data,istatus)
        call get_sfc_badflag(badflag,istatus)

!       Plot Stations
        do i = 1,n_obs_b ! num_sfc
            call latlon_to_rlapsgrid(lat_s(i),lon_s(i),lat,lon
     1                          ,ni,nj,xsta,ysta,istatus)

            if(xsta .lt. 1. .or. xsta .gt. float(ni) .OR.
     1         ysta .lt. 1. .or. ysta .gt. float(nj)          )then       
                    goto80
            endif

            if(c_field(1:2) .eq. 'ow')then ! plot wind only
                if(dd_s(i) .eq. badflag .or. ff_s(i) .eq. badflag)then       
                    goto80
                else
                    t_s(i) = badflag
                    td_s(i) = badflag
                    write(6,*)' Wind only ',i,dd_s(i),ff_s(i),ffg_s(i)
                endif
            endif

!           call supcon(lat_s(i),lon_s(i),usta,vsta)

!           IFLAG = 0        --        Station locations only
!           IFLAG = 1        --        FSL Mesonet only (for WWW)
!           IFLAG = 2        --        All Sfc Obs

            if(iflag .ge. 1)then

                w1 = dd_s(i)
                w2 = ff_s(i)
                w3 = ffg_s(i)

                if(iflag .eq. 1)call setusv_dum(2HIN,14)

                c20_stations = stations(i)

                len_sta_plot = 5 ! maximum station name length allowed in plot

                call left_justify(c20_stations)
                call s_len(c20_stations,len_sta)

                if(len_sta .gt. len_sta_plot)then ! right justify
                    ifirst = len_sta - len_sta_plot + 1
                    c_staname = c20_stations(ifirst:len_sta)
                else
                    c_staname = c20_stations(1:len_sta)
                endif

                charsize = .0040 / zoom_eff

                if(iflag_cv .eq. 0 .and. autostntype(i) .ne. 'CUM')then       
!                   Plot station name & Wx String
                    CALL PCLOQU(xsta, ysta-du2*3.5, c_staname, 
     1                          charsize,ANGD,CNTR)
                    CALL PCLOQU(xsta+du2*1.1, ysta-du2*1.1, wx_s(i),        
     1                              charsize,ANGD,-1.0)
                endif

                relsize = 1.1

                if(iflag .eq. 1)call setusv_dum(2HIN,11)

                if(iflag_cv .eq. 1)then ! Ceiling & Visibility
                    if(vis_s(i) .ne. badflag)then
                        dewpoint = vis_s(i)
                    else
                        dewpoint = badflag
                    endif

                    nlyr = kloud_s(i)

!                   Plot station name
                    if(nlyr .ge. 1 .or. vis_s(i) .ne. badflag)then
                        CALL PCLOQU(xsta, ysta-du2*3.5, c_staname, 
     1                              charsize,ANGD,CNTR)
                    endif

                    if(nlyr .ge. 1)then
                        write(6,*)i,c20_stations(1:min(len_sta+1,20))       
     1                           ,c_staname,nint(xsta),nint(ysta)
     1                           ,reptype(i)(1:5),nlyr,vis_s(i)
     1                           ,store_amt(i,nlyr),store_hgt(i,1)
                    else
                        write(6,*)i,c20_stations(1:min(len_sta+1,20))       
     1                           ,c_staname,nint(xsta),nint(ysta)
     1                           ,reptype(i)(1:5),nlyr,vis_s(i)
                    endif

                    pressure = float(nlyr)        ! number of cloud layers
                    w1       = store_hgt(i,1)     ! height of 1st layer (m)
                    if(nlyr .ge. 2)w2 = store_hgt(i,2)      ! 2nd layer (m)
                    if(nlyr .ge. 3)w3 = store_hgt(i,3)      ! 3rd layer (m)

                    if(nlyr .ge. 1)then                    
                        if(l_parse(store_amt(i,nlyr),'CLR'))then
                            temp = 0.0
                        elseif(l_parse(store_amt(i,nlyr),'SKC'))then
                            temp = 0.0
                        elseif(l_parse(store_amt(i,nlyr),'FEW'))then
                            temp = 0.1
                        elseif(l_parse(store_amt(i,nlyr),'SCT'))then
                            temp = 0.3
                        elseif(l_parse(store_amt(i,nlyr),'BKN'))then
                            temp = 0.7
                        elseif(l_parse(store_amt(i,nlyr),'OVC'))then
                            temp = 1.0
                        else
                            write(6,*)' Unrecognized cloud fraction'
     1                               ,store_amt(i,nlyr)
                            temp = 0.0
                        endif

                        w1 = w1 / 100. ! convert meters MSL to Hecto-meters
                    else ! no cloud layers
                        w1 = r_missing_data

                    endif

                elseif(iflag_cv .eq. 2)then ! Precip
                    temp = badflag
                    if(c_field(3:4) .eq. '24')then
                        pcpval = pcp24(i)
                    elseif(c_field(3:3) .eq. '6')then
                        pcpval = pcp6(i)
                    elseif(c_field(3:3) .eq. '3')then
                        pcpval = pcp3(i)
                    elseif(c_field(3:3) .eq. '1')then
                        pcpval = pcp1(i)
                    else
                        pcpval = pcp1(i)
                    endif
                    if(pcpval .ne. badflag)then
                        if(namelist_parms%c_units_type 
     1                                             .eq. 'metric')then      
                            dewpoint = pcpval * 25.4
                        else
                            dewpoint = pcpval
                        endif
                        write(6,*)' Precip ob ',i,pcpval,dewpoint
                    else
                        dewpoint = badflag
                    endif

                    if(snow(i) .ne. badflag)then
                        if(namelist_parms%c_units_type 
     1                                             .eq. 'metric')then      
                            temp = snow(i) * 2.54
                        else
                            temp = snow(i)
                        endif
                        write(6,*)' Snow ob ',i,snow(i),temp
                    else
                        temp = badflag
                    endif

                    pressure = r_missing_data

                    call s_len(wx_s(i),lenwx)

                    call s_len(c_field,len_field)

!                   Test for valid report
                    if(pcpval .ne. badflag .or. 
     1                 snow(i) .ne. badflag .or. 
     1                 (wx_s(i) .ne. 'UNK' .and. lenwx .gt. 0) 
     1                                                         )then       

                        write(6,11)i,pcpval,snow(i),lenwx,c_staname
     1                                                   ,wx_s(i)
11                      format('  Plot Precip ob ',i4,2f8.2,i3,1x,a
     1                                                        ,1x,a)
!                       write(6,*)pcp1(i),pcp3(i),pcp6(i),pcp24(i)

!                       Plot Weather String
                        if(len_field .le. 2)then
                          if(lenwx .gt. 0 .and. wx_s(i) .ne. 'UNK')then
                            CALL PCLOQU(xsta-du2*0.9, ysta-du2*1.5
     1                        , wx_s(i)(1:lenwx), charsize,ANGD,+1.0)
                          endif
                        endif

!                       Plot name and Station Location
                        CALL PCLOQU(xsta, ysta-du2*3.5, c_staname, 
     1                              charsize,ANGD,CNTR)

                        call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)        
                        call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)

                    endif

                elseif(iflag_cv .eq. 3)then ! Soil/Water T (& solar radiation)
                    iplotsta = 0
                    temp = sfct(i)
                    dewpoint = sfcm(i)
                    pressure = badflag

                    if(temp .ne. badflag)then
                        write(6,*)' Sfc T  = ',i,temp,c_staname
                        iplotsta = 1
                    endif

                    if(dewpoint .ne. badflag)then
                        write(6,*)' Sfc RH = ',i,dewpoint,c_staname
                        iplotsta = 1
                    endif

                    if(solar_s(i) .ne. badflag .and.
     1                 solar_s(i) .ge. 0.            )then
                        write(6,*)' Solar Rad = ',i,solar_s(i),c_staname       
                        pressure = solar_s(i)
                        iplotsta = 1
                    endif
                  
                    if(iplotsta .eq. 1)then
!                       Plot name and Station Location
                        CALL PCLOQU(xsta, ysta-du2*3.5, c_staname, 
     1                              charsize,ANGD,CNTR)
                        call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
                        call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)
                    endif

                elseif(iflag_cv .eq. 4)then ! Solar Radiation                    
                    iplotsta = 0
                    temp = badflag 
                    dewpoint = badflag
                    pressure = badflag

                    if(solar_s(i) .ne. badflag .and.
     1                 solar_s(i) .ge. 0.            )then
                        write(6,*)' Solar Rad = ',i,solar_s(i),c_staname       
                        pressure = solar_s(i)
                        iplotsta = 1
                    endif
                  
                    if(iplotsta .eq. 1)then
!                       Plot name and Station Location
                        CALL PCLOQU(xsta, ysta-du2*2.5, c_staname(1:4), 
     1                              charsize,ANGD,CNTR)
                        call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
                        call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)
                    endif

                elseif(c_field(2:2) .ne. 'c')then ! Fahrenheit
                    temp = t_s(i)
                    dewpoint = td_s(i)

                else                          ! Celsius
                    if(t_s(i) .ne. badflag)then
                        temp = f_to_c(t_s(i))
                    else
                        temp = badflag
                    endif

                    if(td_s(i) .ne. badflag)then
                        dewpoint = f_to_c(td_s(i))
                    else
                        dewpoint = badflag
                    endif
                endif

                if(iflag_cv .eq. 0)then
                    if(c3_presob .eq. 'msl')then
                        pressure = pmsl_s(i)
                    elseif(c3_presob .eq. 'alt')then
                        pressure = alt_s(i)
                    elseif(c3_presob .eq. 'stn')then 
                        pressure = pstn_s(i)
                    else
                        pressure = r_missing_data
                    endif

                endif

                if(autostntype(i) .ne. 'CUM' .or. iflag_cv .eq. 2)then 
!                   exclude CWB precip stations for non-precip plots
                    call plot_mesoob(w1,w2,w3
     1                 ,temp,dewpoint
     1                 ,pressure,xsta,ysta
     1                 ,lat,lon,ni,nj,relsize,zoom,n_obs_g,11,du2
     1                 ,wx_s(i)
     1                 ,iflag,iflag_cv,namelist_parms,plot_parms)
                endif

                if(iflag .eq. 1)call setusv_dum(2HIN,33)

            else ! Write station location only
                if(c_field(1:2) .eq. 'mw')then ! Set Mesowx Red Color
                    call setusv_dum(2hIN,3)
                endif

                call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
                call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)

            endif

80      enddo ! i

        if(iflag .eq. 1)then ! special mesonet label 
            call setusv_dum(2hIN,2)
            call cv_i4tim_asc_lp(i4time,atime_s,istatus)
            atime_s = atime_s(1:14)//atime_s(16:17)//' '
            ix = 590
            iy = 270
            call pwrity(cpux(ix),cpux(iy),atime_s(1:17),17,-1,0,-1)
        endif

        call sflush

        return
        end

c
c
        subroutine plot_mesoob(dir,spd,gust,t,td,p,ri,rj
     1                        ,lat,lon,imax,jmax,relsize_in,zoom,nobs
     1                        ,icol_in,du2,wx,iflag,iflag_cv
     1                        ,namelist_parms,plot_parms)

        include 'lapsplot.inc'

        real lat(imax,jmax),lon(imax,jmax)
        character*3 t1,td1,p1
        character*4 c4_pcp
        character*25 wx

        call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)
!       write(6,1234) mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype
 1234   format(1x,4i5,4e12.4,i4)

        obs_size = plot_parms%contour_line_width

        ANGD = 0.
        CNTR = 0.

!       At zoom=1-zoom_max, make the obs plot larger if there are few obs
!       zoom_eff is inversely related to plot size
        zoom_max = 1.5
        if(zoom .lt. zoom_max)then
            if(nobs .gt. 30)then
                zoom_eff = 1.0                 ! smaller obs 
            else
                zoom_eff = zoom / zoom_max     ! larger obs
            endif
        else
            zoom_eff = zoom / zoom_max         ! larger obs
        endif

        zoom_eff = zoom_eff / obs_size

        relsize = relsize_in / zoom_eff

        du_b=(jmax)/240. * relsize

        jsize = nint(0.4 * relsize) - 1

!       write(6,*)' relsize,du_b,jsize,zoom,obs_size,zoom_eff = '
!    1             ,relsize,du_b,jsize,zoom,obs_size,zoom_eff

        call get_border(imax,jmax,x_1,x_2,y_1,y_2)
        call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax),1)

!       rot = (standard_longitude - lon(nint(ri),nint(rj))) / 57.295

        nri = min(max(nint(ri),1),imax)
        nrj = min(max(nint(rj),1),jmax)
        rot = projrot_latlon(lat(nri,nrj)
     1                      ,lon(nri,nrj),istatus) / 57.295

!       write(6,*)' lat/lon/rot ',lat(nri,nrj),lon(nri,nrj),rot

!       Convert ri and rj to x1 and y1 (U and V)
!       call supcon(alat,alon,x1,y1)
        x1 = umin + (umax - umin) * (ri-1.) / float(imax-1)
        y1 = vmin + (vmax - vmin) * (rj-1.) / float(jmax-1)

        xsta=ri
        ysta=rj

        u = ri
        v = rj

        if(iflag .eq. 3)then ! Plot on top of station location for 'tmg'
            du = 0.
        else
            du = du_b
        endif

        dv   = 1.2 * du
        du_t = 3.0 * du
        du_p = 3.0 * du

        call pcsetr('CE',0.)

        charsize = .0040 / zoom_eff

        if(iflag_cv .eq. 0)then ! Normal obs plot
            if(iflag .ne. 3)then ! Not a 'tmg' plot
                if(dir .ge. 0.  .and. spd .ge. 0. .and.
     1             dir .le. 360 .and. spd .le. 200.       )then
                    call barbs(spd,dir,ri,rj,du_b,rot
     1                        ,-1e10,+1e10,-1e10,+1e10,1.0,1.0)
                    if(spd .ge. 1.0)then
                        call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
                        call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)       
                    endif
                else
                    call line(xsta,ysta+du2*0.5,xsta,ysta-du2*0.5)
                    call line(xsta+du2*0.5,ysta,xsta-du2*0.5,ysta)
                endif
            endif

!           Plot Temperature       
            if(t.gt.-75. .and. t.lt.140.) then 
               write(t1,100,err=20) nint(t)
!              call pwrity(u-du_t,v+dv,t1,3,jsize,0,0)
               CALL PCLOQU(u-du_t,v+dv,t1,charsize,ANGD,CNTR)
            endif
 100        format(i3)
 20         continue

!           Plot Dew Point
            if(td.gt.-75. .and. td.lt.100.) then
               write(td1,100,err=30) nint(td)
               CALL PCLOQU(u-du_t,v-dv,td1,charsize,ANGD,CNTR)
            endif
 30         continue
 
!           Plot Pressure
            if(p .gt. 0. .and. p .lt. 10000.) then
               if(p .ge. 1000.) p = p - 1000.
               ip = ifix( p )
               write(p1,101,err=140) ip
 101           format(i3.3)
!              call pwrity(u+du_p,v+dv,p1,3,jsize,0,0)
               CALL PCLOQU(u+du_p,v+dv,p1,charsize,ANGD,CNTR)
            endif
 140        continue

!           Plot Gusts (FSL WWW)
            if(iflag .ge. 1 .and. iflag .ne. 3)then ! Not a 'tmg' plot
               if(gust .gt. 20)then
                   ig = int(gust)
                   write(p1,102,err=150) ig
!                  call setusv_dum(2HIN,4)
                   dg = 3.0 * du
!                  call pwrity(u,v+dg,p1,3,jsize,0,0)          ! On Top
!                  call pwrity(u+du_p,v-dv,p1,3,jsize,0,0)     ! Lower Right
                   CALL PCLOQU(u+du_p,v-dv,p1,charsize,ANGD,CNTR)
 102               format('G',i2)
!                  call setusv_dum(2HIN,icol_in)
               endif
            endif           
 150        continue

        elseif(iflag_cv .eq. 1)then ! C&V plot
!           Plot outer circle (use p for the number of layers?)
            if(p .ge. 1.0)then
                call plot_circle(u,v,du*0.8) 

!               Plot cloud cover (using t to hold the variable)
                call plot_circle_fill(u,v,du*0.8,t)

            endif ! number of layers

!           Plot Visibility (using td to hold the variable)
            if(td.gt.-75. .and. td.lt.100.) then
               write(td1,100,err=31) nint(td)
               call left_justify(td1)
               CALL PCLOQU(u+du_t,v-dv,td1,charsize,ANGD,CNTR)
            endif
 31         continue

!           Plot Cloud Height (dir/w1 variable)
            if(dir.gt.-75. .and. dir.lt.140.) then 
               write(t1,100,err=38) nint(dir) ! HectoMeters of lowest layer
               call left_justify(t1)
               CALL PCLOQU(u+du_t,v+dv,t1,charsize,ANGD,CNTR)
            endif
 38         continue

        elseif(iflag_cv .eq. 2)then ! Precip obs plot
!           Plot 1hr Precip Ob
            if(td.gt.-75. .and. td.lt.100.) then
               if(namelist_parms%c_units_type .eq. 'english')then
                   write(c4_pcp,103,err=32) td
 103               format(f4.2)
               else ! millimeters for metric
                   write(c4_pcp,104,err=32) td
 104               format(f4.1)
               endif
               call left_justify(c4_pcp)
               CALL PCLOQU(u+du_t,v-dv,c4_pcp,charsize,ANGD,CNTR)
               write(6,*)' Precip (td) plot = ',c4_pcp
            endif
 32         continue

        elseif(iflag_cv .eq. 3 .or. iflag_cv .eq. 4)then ! Sfc T, RH & Solar Radiation plot
!           Plot Temperature       
            if(t.gt.-75. .and. t.lt.140.) then 
               write(t1,100,err=40) nint(t)
               CALL PCLOQU(u-du_t,v+dv,t1,charsize,ANGD,CNTR)
            endif
 40         continue

!           Plot Sfc RH (td variable)
            if(td.ge.0. .and. td.le.100.) then
               write(td1,100,err=50) nint(td)
               CALL PCLOQU(u-du_t,v-dv,td1,charsize,ANGD,CNTR)
            endif
 50         continue

!           Plot Solar Radiation (pressure variable)
            if(p .gt. 0. .and. p .lt. 10000.) then
               if(p .ge. 1000.)then
                   p = p - 1000.
                   ip = ifix( p )
                   write(p1,201,err=60) ip
 201               format(i3.3)
               else
                   ip = ifix( p )
                   write(p1,202,err=60) ip
 202               format(i3)
               endif
               CALL PCLOQU(u+du_p,v+dv,p1,charsize,ANGD,CNTR)
            endif

 60         continue

        endif
c

        call sflush

        return
        end
