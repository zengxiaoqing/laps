
        subroutine plot_sounding(i4time_ref,lun,NX_L,NY_L,NZ_L
     1                          ,r_missing_data,laps_cycle_time,maxstns)       

        real*4 pres_3d(NX_L,NY_L,NZ_L)
        real*4 field_3d(NX_L,NY_L,NZ_L)
!       real*4 rh_3d(NX_L,NY_L,NZ_L)

        real*4 temp_vert(NZ_L)
        real*4 rh_vert(NZ_L)
        real*4 td_vert(NZ_L)

        real*4 pres_1d(NZ_L)
        real*4 logp_1d(NZ_L), logp_bottom, logp_top, logp
  
        real*4 k_to_c

        character*1 c_prodtype
        character*3 var_2d
        character*150  directory
        character*31  ext
        character*10  units_2d
        character*125 comment_2d
        character*9 a9time
        character*4 fcst_hhmm
        character*3 c3_string
        character*4 c4_string
        character*33 c33_label

        integer i_overlay
        save i_overlay
        data i_overlay /0/

        include 'icolors.inc'

        common /image/ n_image

        skewt(t_c,logp) = t_c - (logp - logp_bottom) * 32.

        n_image = 0

        write(6,*)' Input x grid point for sounding...'

        read(5,*)xsound

        write(6,*)' Input y grid point for sounding...'

        read(5,*)ysound

        tlow_c = -40.
        thigh_c = +40.

!       Read appropriate 3-D fields
50      call input_product_info(i4time_ref              ! I
     1                         ,laps_cycle_time         ! I
     1                         ,3                       ! I
     1                         ,c_prodtype              ! O
     1                         ,ext                     ! O
     1                         ,directory               ! O
     1                         ,a9time                  ! O
     1                         ,fcst_hhmm               ! O
     1                         ,i4_initial              ! O
     1                         ,i4_valid                ! O
     1                         ,istatus)                ! O

!       Read Temperature
        if(c_prodtype .eq. 'A')then
            iflag_temp = 1 ! Returns Ambient Temp (K)

            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                      ,NX_L,NY_L,NZ_L,field_3d,istatus)
            if(istatus .ne. 1)goto900

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            c33_label = 'LAPS Analysis Sounding'

        elseif(c_prodtype .eq. 'N')then
            call get_directory('balance',directory,len_dir)
            ext = 'lt1'
            directory = directory(1:len_dir)//ext(1:3)

            var_2d = 'T3'

            call get_3dgrid_dname(directory
     1                  ,i4time_ref,laps_cycle_time*10000,i4time_nearest       
     1                  ,ext,var_2d,units_2d
     1                  ,comment_2d,NX_L,NY_L,NZ_L,field_3d,istatus)       

            call make_fnam_lp(i4time_nearest,a9time,istatus)
            c33_label = 'LAPS Balanced Sounding'

        elseif(c_prodtype .eq. 'B' .or. c_prodtype .eq. 'F')then ! Bkg or Fcst
            var_2d = 'T3'
            call get_lapsdata_3d(i4_initial,i4_valid
     1                              ,NX_L,NY_L,NZ_L       
     1                              ,directory,var_2d
     1                              ,units_2d,comment_2d,field_3d
     1                              ,istatus)

            call make_fnam_lp(i4_valid,a9time,istatus)

            if(c_prodtype .eq. 'B')then
                c33_label = 'LAPS Background Sounding '//fcst_hhmm
            elseif(c_prodtype .eq. 'F')then
                c33_label = 'LAPS Forecast Sounding '//fcst_hhmm
            endif

        else
            write(6,*)' Unknown choice, will quit'
            go to 900

        endif

        i_overlay = i_overlay + 1

        call interp_3d(field_3d,temp_vert,xsound,xsound,ysound,ysound,       
     1                 NX_L,NY_L,NZ_L,1,NZ_L,r_missing_data)

!       Read RH   (a la xsect)
        if(c_prodtype .eq. 'A')then
            var_2d = 'RHL'
            ext = 'lh3'
            call get_laps_3dgrid
     1          (i4time_nearest,0,i4time_nearest,NX_L,NY_L,NZ_L       
     1          ,ext,var_2d,units_2d,comment_2d,field_3d,istat_rh)
            if(istat_rh .ne. 1)goto100

            call make_fnam_lp(i4time_nearest,a9time,istatus)

        else
            write(6,*)' Sorry, RH not yet supported for prodtype: '
     1               ,c_prodtype
            istat_rh = 0
            goto100

        endif

        call interp_3d(field_3d,rh_vert,xsound,xsound,ysound,ysound,       
     1                 NX_L,NY_L,NZ_L,1,NZ_L,r_missing_data)

        do iz = 1,NZ_L
            t_c         = k_to_c(temp_vert(iz)) 
            td_vert(iz) = DWPT(t_c,rh_vert(iz))
        enddo ! iz

!       Read Wind (a la xsect)

 100    continue

!       Get 3-D pressure field
        call get_pres_3d(i4_valid,NX_L,NY_L,NZ_L,pres_3d,istatus)
        if(istatus .ne. 1)go to 900

        ix = nint(xsound)
        iy = nint(ysound)

        do iz = 1,NZ_L
            pres_1d(iz) = pres_3d(ix,iy,iz)
            logp_1d(iz) = log(pres_1d(iz))
        enddo ! iz

        pbottom = 110000.
        ptop = 10000.

        logp_bottom = log(pbottom)
        logp_top = log(ptop)

!       Do plots on non-skew T, log P scale

        call set(0.1   , 0.9    , 0.1        , 0.9,
     1           tlow_c, thigh_c, logp_bottom, logp_top, 1)

        call setusv_dum(2hIN,icolors(i_overlay))

!       Plot temp sounding
        do iz = 2,NZ_L
            y1 = logp_1d(iz-1)
            y2 = logp_1d(iz)
            x1 = skewt(k_to_c(temp_vert(iz-1)),y1)
            x2 = skewt(k_to_c(temp_vert(iz)),y2)
            call line(x1,y1,x2,y2) 
        enddo ! iz

!       Plot dewpoint sounding
        if(istat_rh .eq. 1)then
            do iz = 2,NZ_L
                y1 = logp_1d(iz-1)
                y2 = logp_1d(iz)
                x1 = skewt(td_vert(iz-1),y1)
                x2 = skewt(td_vert(iz),y2)
                call line(x1,y1,x2,y2) 
            enddo ! iz
        endif

        call setusv_dum(2hIN,7) ! Yellow

!       Plot temp scale
        do it = -80, thigh_c, 10
            y1 = logp_bottom
            y2 = logp_top
            x1 = skewt(float(it),y1)
            x2 = skewt(float(it),y2)
            call line(x1,y1,x2,y2) 
        enddo ! it

!       Plot Horizontal lines at 100mb intervals
        do ip = 100,1100,100
            p_pa = float(ip) * 100.
            x1 = tlow_c
            x2 = thigh_c
            y1 = log(p_pa)
            y2 = log(p_pa)
            call line(x1,y1,x2,y2) 
        enddo ! ip

!       Plot Rest of Box
        x1 = tlow_c
        y1 = logp_bottom
        x2 = thigh_c
        y2 = logp_bottom
        call line(x1,y1,x2,y2) 

        x1 = tlow_c
        y1 = logp_top
        x2 = thigh_c
        y2 = logp_top
        call line(x1,y1,x2,y2) 

        x1 = tlow_c
        y1 = logp_bottom
        x2 = tlow_c
        y2 = logp_top
        call line(x1,y1,x2,y2) 

        x1 = thigh_c
        y1 = logp_bottom
        x2 = thigh_c
        y2 = logp_top
        call line(x1,y1,x2,y2) 


!       Plot Pressure labels at 100mb intervals
        width = thigh_c - tlow_c
        height = logp_top - logp_bottom
        call set(0.00   , 1.00    , 0.0        , 1.0,
     1           tlow_c-width/8., thigh_c+width/8., 
     1           logp_bottom-height/8., logp_top+height/8.   , 1)

        do ip = 100,1100,100
            p_pa = float(ip) * 100.
            y = log(p_pa) 

!           Pressure
            x = -43.
            write(c4_string,2014)ip
            call pwrity (x, y, c4_string, 4, 1, 0, 0)
2014        format(i4)
        enddo ! ip

        do it = -40, +40, 10
            t_c = float(it)
            write(c3_string,2013)it
2013        format(i3)
            x = t_c
            y = logp_bottom + .04
            call pwrity (x, y, c3_string, 3, 1, 0, 0)
        enddo

        call setusv_dum(2hIN,icolors(i_overlay))

        call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
        call write_label_lplot(100,94,c33_label,a9time
     1                                          ,i_overlay,'sound')


        write(6,*)' Sounding has been plotted...'

        go to 50

 900    continue

        call sflush

        return
        end
